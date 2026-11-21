use std::path::{Path, PathBuf};

const UNLIMITED_MIN_BYTES: u64 = 1 << 50;

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub(crate) enum MemorySource {
    NodeAvailable,
    CgroupLimit,
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub(crate) struct MemoryBudget {
    pub(crate) bytes: u64,
    pub(crate) source: MemorySource,
}

fn meminfo_available_bytes() -> Option<u64> {
    let meminfo = std::fs::read_to_string("/proc/meminfo").ok()?;
    parse_meminfo_available(&meminfo)
}

fn parse_meminfo_available(meminfo: &str) -> Option<u64> {
    for line in meminfo.lines() {
        if let Some(rest) = line.strip_prefix("MemAvailable:") {
            let kb: u64 = rest.split_whitespace().next()?.parse().ok()?;
            return Some(kb.saturating_mul(1024));
        }
    }
    None
}

pub(crate) fn cgroup_memory_limit() -> Option<u64> {
    #[cfg(target_os = "linux")]
    {
        let proc_self_cgroup = std::fs::read_to_string("/proc/self/cgroup").ok()?;
        cgroup_memory_limit_from(&proc_self_cgroup, |p| std::fs::read_to_string(p).ok())
    }
    #[cfg(not(target_os = "linux"))]
    {
        None
    }
}

#[cfg_attr(not(target_os = "linux"), allow(dead_code))]
pub(crate) fn cgroup_memory_limit_from(
    proc_self_cgroup: &str,
    read: impl Fn(&Path) -> Option<String>,
) -> Option<u64> {
    const V2_MOUNT: &str = "/sys/fs/cgroup";
    const V1_MEMORY_MOUNT: &str = "/sys/fs/cgroup/memory";
    const V2_FILES: &[&str] = &["memory.max", "memory.high"];
    const V1_FILES: &[&str] = &["memory.limit_in_bytes"];

    let mut best: Option<u64> = None;
    for line in proc_self_cgroup.lines() {
        let mut fields = line.splitn(3, ':');
        let (controllers, path) = match (fields.next(), fields.next(), fields.next()) {
            (Some(_), Some(controllers), Some(path)) => (controllers, path),
            _ => continue,
        };
        let found = if controllers.is_empty() {
            min_limit_over_chain(Path::new(V2_MOUNT), path, V2_FILES, &read)
        } else if controllers.split(',').any(|c| c == "memory") {
            min_limit_over_chain(Path::new(V1_MEMORY_MOUNT), path, V1_FILES, &read)
        } else {
            None
        };
        if let Some(v) = found {
            best = Some(best.map_or(v, |b| b.min(v)));
        }
    }
    best
}

fn min_limit_over_chain<F: Fn(&Path) -> Option<String>>(
    mount: &Path,
    cgroup_path: &str,
    files: &[&str],
    read: &F,
) -> Option<u64> {
    let mut dir = mount.join(cgroup_path.trim_start_matches('/'));
    if !dir.starts_with(mount) {
        return None;
    }
    let mut best: Option<u64> = None;
    loop {
        for file in files {
            if let Some(v) = read(&dir.join(file)).as_deref().and_then(parse_limit_value) {
                best = Some(best.map_or(v, |b| b.min(v)));
            }
        }
        if dir == mount {
            break;
        }
        match dir.parent() {
            Some(parent) if parent.starts_with(mount) => dir = parent.to_path_buf(),
            _ => break,
        }
    }
    best
}

fn parse_limit_value(raw: &str) -> Option<u64> {
    let text = raw.trim();
    if text.is_empty() || text == "max" {
        return None;
    }
    let bytes: u64 = text.parse().ok()?;
    if bytes == 0 || bytes >= UNLIMITED_MIN_BYTES {
        return None;
    }
    Some(bytes)
}

pub(crate) fn available_memory_budget() -> Option<MemoryBudget> {
    budget_from(meminfo_available_bytes(), cgroup_memory_limit())
}

fn budget_from(node_available: Option<u64>, cgroup_limit: Option<u64>) -> Option<MemoryBudget> {
    match (node_available, cgroup_limit) {
        (Some(node), Some(limit)) if limit < node => Some(MemoryBudget {
            bytes: limit,
            source: MemorySource::CgroupLimit,
        }),
        (Some(node), _) => Some(MemoryBudget {
            bytes: node,
            source: MemorySource::NodeAvailable,
        }),
        (None, Some(limit)) => Some(MemoryBudget {
            bytes: limit,
            source: MemorySource::CgroupLimit,
        }),
        (None, None) => None,
    }
}

fn process_resident_bytes() -> Option<u64> {
    let status = std::fs::read_to_string("/proc/self/status").ok()?;
    parse_status_resident(&status)
}

fn parse_status_resident(status: &str) -> Option<u64> {
    for line in status.lines() {
        if let Some(rest) = line.strip_prefix("VmRSS:") {
            let kb: u64 = rest.split_whitespace().next()?.parse().ok()?;
            return Some(kb.saturating_mul(1024));
        }
    }
    None
}

/// Share of the free memory left to the rest of the process when the rescoring budget is
/// sized automatically: allocator slack, buffers outside the budget and the estimate's error.
const RESCORE_MARGIN_SHARE: f64 = 0.1;

/// Smallest margin left when the rescoring budget is sized automatically.
const RESCORE_MARGIN_MIN_BYTES: u64 = 1 << 30;

/// Allowance assumed when neither the cgroup limit nor `MemAvailable` can be read.
const RESCORE_FALLBACK_BYTES: u64 = 16 << 30;

/// The native rescoring engine's memory limit. `gb` (`--rescore-memory-gb`) is a fixed
/// budget. Without it the engine sizes the budget to the PIN, never above
/// `DEFAULT_RESCORE_MEMORY_GB` nor above the memory free to this process now: the cgroup limit
/// less what the process holds (when there is a limit) and `MemAvailable`, whichever is lower,
/// less a margin.
pub(crate) fn rescore_memory(gb: Option<u64>, spill_dir: PathBuf) -> output::MemLimit {
    match gb {
        Some(gb) => output::MemLimit {
            bytes: gb.saturating_mul(1 << 30),
            spill_dir,
            sizing: output::Sizing::Fixed,
        },
        None => {
            let (bytes, why) = rescore_allowance(
                meminfo_available_bytes(),
                cgroup_memory_limit(),
                process_resident_bytes(),
            );
            output::MemLimit {
                bytes,
                spill_dir,
                sizing: output::Sizing::FitInput {
                    cap: crate::cli::DEFAULT_RESCORE_MEMORY_GB << 30,
                    why,
                },
            }
        }
    }
}

fn gb(bytes: u64) -> f64 {
    bytes as f64 / (1u64 << 30) as f64
}

fn rescore_allowance(
    node_available: Option<u64>,
    cgroup_limit: Option<u64>,
    resident: Option<u64>,
) -> (u64, String) {
    let in_cgroup = cgroup_limit.map(|l| {
        let held = resident.unwrap_or(0);
        (l.saturating_sub(held), l, held)
    });
    let (free, source) = match (node_available, in_cgroup) {
        (Some(node), Some((c, l, held))) if c < node => (
            c,
            format!(
                "the job's {:.1} GB limit less the {:.1} GB this process holds",
                gb(l),
                gb(held)
            ),
        ),
        (Some(node), _) => (node, format!("{:.1} GB available on the machine", gb(node))),
        (None, Some((c, l, held))) => (
            c,
            format!(
                "the job's {:.1} GB limit less the {:.1} GB this process holds",
                gb(l),
                gb(held)
            ),
        ),
        (None, None) => {
            return (
                RESCORE_FALLBACK_BYTES,
                format!(
                    "{:.1} GB assumed, the free memory could not be read",
                    gb(RESCORE_FALLBACK_BYTES)
                ),
            )
        }
    };
    let margin = ((free as f64 * RESCORE_MARGIN_SHARE) as u64).max(RESCORE_MARGIN_MIN_BYTES);
    let bytes = free.saturating_sub(margin);
    (
        bytes,
        format!(
            "{:.1} GB allowed: {source}, less a {:.1} GB margin",
            gb(bytes),
            gb(margin)
        ),
    )
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::collections::HashMap;

    const GIB: u64 = 1 << 30;

    fn fake_fs(files: &[(&str, &str)]) -> impl Fn(&Path) -> Option<String> {
        let map: HashMap<String, String> = files
            .iter()
            .map(|(k, v)| ((*k).to_string(), (*v).to_string()))
            .collect();
        move |p: &Path| map.get(&p.to_str()?.replace('\\', "/")).cloned()
    }

    #[test]
    fn cgroup_v2_numeric_memory_max_is_the_limit() {
        let limit = cgroup_memory_limit_from(
            "0::/slurm/uid_1000/job_42/step_0\n",
            fake_fs(&[(
                "/sys/fs/cgroup/slurm/uid_1000/job_42/step_0/memory.max",
                "257698037760\n",
            )]),
        );
        assert_eq!(limit, Some(257_698_037_760));
    }

    #[test]
    fn cgroup_v2_max_means_unlimited() {
        let limit = cgroup_memory_limit_from(
            "0::/user.slice\n",
            fake_fs(&[
                ("/sys/fs/cgroup/user.slice/memory.max", "max\n"),
                ("/sys/fs/cgroup/memory.max", "max\n"),
            ]),
        );
        assert_eq!(limit, None);
    }

    #[test]
    fn cgroup_v2_memory_high_below_memory_max_wins() {
        let limit = cgroup_memory_limit_from(
            "0::/pod\n",
            fake_fs(&[
                ("/sys/fs/cgroup/pod/memory.max", "8589934592\n"),
                ("/sys/fs/cgroup/pod/memory.high", "4294967296\n"),
            ]),
        );
        assert_eq!(limit, Some(4 * GIB));
    }

    #[test]
    fn cgroup_v1_numeric_limit_in_bytes_is_the_limit() {
        let limit = cgroup_memory_limit_from(
            "9:memory:/docker/abc123\n5:cpu,cpuacct:/docker/abc123\n",
            fake_fs(&[(
                "/sys/fs/cgroup/memory/docker/abc123/memory.limit_in_bytes",
                "2147483648\n",
            )]),
        );
        assert_eq!(limit, Some(2 * GIB));
    }

    #[test]
    fn cgroup_v1_huge_sentinel_means_unlimited() {
        let limit = cgroup_memory_limit_from(
            "9:memory:/\n",
            fake_fs(&[(
                "/sys/fs/cgroup/memory/memory.limit_in_bytes",
                "9223372036854771712\n",
            )]),
        );
        assert_eq!(limit, None);
    }

    #[test]
    fn missing_cgroup_files_yield_no_limit() {
        let limit = cgroup_memory_limit_from("0::/slurm/job_42\n", fake_fs(&[]));
        assert_eq!(limit, None);
    }

    #[test]
    fn malformed_cgroup_value_yields_no_limit() {
        let limit = cgroup_memory_limit_from(
            "0::/slurm/job_42\n",
            fake_fs(&[("/sys/fs/cgroup/slurm/job_42/memory.max", "not-a-number\n")]),
        );
        assert_eq!(limit, None);
    }

    #[test]
    fn limit_on_a_parent_cgroup_is_found_when_the_leaf_has_none() {
        let limit = cgroup_memory_limit_from(
            "0::/slurm/uid_1000/job_42/step_0\n",
            fake_fs(&[
                (
                    "/sys/fs/cgroup/slurm/uid_1000/job_42/step_0/memory.max",
                    "max\n",
                ),
                (
                    "/sys/fs/cgroup/slurm/uid_1000/job_42/memory.max",
                    "257698037760\n",
                ),
            ]),
        );
        assert_eq!(limit, Some(257_698_037_760));
    }

    #[test]
    fn budget_is_the_minimum_of_node_memory_and_the_cgroup_limit() {
        let budget = budget_from(Some(641 * GIB), Some(240 * GIB)).unwrap();
        assert_eq!(budget.bytes, 240 * GIB);
        assert_eq!(budget.source, MemorySource::CgroupLimit);

        let budget = budget_from(Some(8 * GIB), Some(64 * GIB)).unwrap();
        assert_eq!(budget.bytes, 8 * GIB);
        assert_eq!(budget.source, MemorySource::NodeAvailable);

        let budget = budget_from(Some(16 * GIB), None).unwrap();
        assert_eq!(budget.bytes, 16 * GIB);
        assert_eq!(budget.source, MemorySource::NodeAvailable);

        assert_eq!(budget_from(None, None), None);
    }

    #[test]
    fn rescore_allowance_is_the_tighter_of_job_and_machine_less_a_margin() {
        let (bytes, why) = rescore_allowance(Some(600 * GIB), Some(64 * GIB), Some(4 * GIB));
        assert_eq!(bytes, 60 * GIB - 6 * GIB);
        assert!(why.contains("limit") && why.contains("holds"), "{why}");

        let (bytes, why) = rescore_allowance(Some(20 * GIB), Some(64 * GIB), Some(4 * GIB));
        assert_eq!(bytes, 18 * GIB);
        assert!(why.contains("available on the machine"), "{why}");

        let (bytes, _) = rescore_allowance(Some(5 * GIB), None, None);
        assert_eq!(bytes, 4 * GIB);

        let (bytes, _) = rescore_allowance(None, Some(8 * GIB), Some(9 * GIB));
        assert_eq!(bytes, 0);

        let (bytes, why) = rescore_allowance(None, None, None);
        assert_eq!(bytes, RESCORE_FALLBACK_BYTES);
        assert!(why.contains("assumed"), "{why}");
    }

    #[test]
    fn a_fixed_rescore_budget_is_the_flag_in_gib() {
        let m = rescore_memory(Some(3), PathBuf::from("/s"));
        assert_eq!(m.bytes, 3 * GIB);
        assert_eq!(m.sizing, output::Sizing::Fixed);
        let m = rescore_memory(None, PathBuf::from("/s"));
        assert!(matches!(m.sizing, output::Sizing::FitInput { .. }));
    }

    #[test]
    fn status_resident_set_is_parsed_as_kilobytes() {
        let text = "Name:\tandes\nVmPeak:\t  900 kB\nVmRSS:\t  2048 kB\n";
        assert_eq!(parse_status_resident(text), Some(2048 * 1024));
        assert_eq!(parse_status_resident("Name: andes\n"), None);
    }

    #[test]
    fn meminfo_available_is_parsed_as_kilobytes() {
        let text = "MemTotal:       672000000 kB\nMemFree:          100 kB\n\
                    MemAvailable:   672103936 kB\n";
        assert_eq!(parse_meminfo_available(text), Some(672_103_936 * 1024));
        assert_eq!(parse_meminfo_available("MemTotal: 1 kB\n"), None);
    }
}
