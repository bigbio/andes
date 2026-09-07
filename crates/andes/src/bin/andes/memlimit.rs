//! How much memory this process may actually use.
//!
//! `MemAvailable` in `/proc/meminfo` is a NODE-wide figure. Under SLURM, Docker,
//! Kubernetes or a systemd slice the process sees the whole host while a cgroup
//! caps what it is allowed to touch, so planning against `MemAvailable` there
//! plans against memory that is not ours: a `--mem=240G` SLURM job on a 641 GiB
//! node budgeted 384 GiB for the in-RAM candidate index and was OOM-killed.
//!
//! The effective figure is therefore `min(MemAvailable, applicable cgroup limit)`,
//! and callers report which of the two they used.

use std::path::Path;

/// A cgroup limit at or above this is treated as "no limit": cgroup v1 writes a
/// page-aligned `LONG_MAX` sentinel (commonly 9223372036854771712) for unlimited,
/// and no real allowance is a petabyte.
const UNLIMITED_MIN_BYTES: u64 = 1 << 50; // 1 PiB

/// Where a memory budget came from, so a cluster log distinguishes the node's
/// free memory from this job's allowance.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub(crate) enum MemorySource {
    /// `MemAvailable` in `/proc/meminfo` — the whole node.
    NodeAvailable,
    /// A cgroup memory limit (v2 `memory.max`/`memory.high`, or v1
    /// `memory.limit_in_bytes`) that is lower than the node figure.
    CgroupLimit,
}

impl MemorySource {
    /// Short phrase for log lines, read as "… of 240 GiB from <this>".
    pub(crate) fn describe(self) -> &'static str {
        match self {
            MemorySource::NodeAvailable => "node MemAvailable",
            MemorySource::CgroupLimit => "cgroup memory limit",
        }
    }
}

/// The memory this process may plan against, and where the number came from.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub(crate) struct MemoryBudget {
    pub(crate) bytes: u64,
    pub(crate) source: MemorySource,
}

/// `MemAvailable` from Linux `/proc/meminfo`, in bytes.
///
/// `None` when it cannot be determined (non-Linux, restricted sandbox).
fn meminfo_available_bytes() -> Option<u64> {
    let meminfo = std::fs::read_to_string("/proc/meminfo").ok()?;
    parse_meminfo_available(&meminfo)
}

/// Pull `MemAvailable` (kB) out of `/proc/meminfo` text and return bytes.
fn parse_meminfo_available(meminfo: &str) -> Option<u64> {
    for line in meminfo.lines() {
        if let Some(rest) = line.strip_prefix("MemAvailable:") {
            // e.g. "MemAvailable:   12345678 kB"
            let kb: u64 = rest.split_whitespace().next()?.parse().ok()?;
            return Some(kb.saturating_mul(1024));
        }
    }
    None
}

/// The memory limit of the cgroup this process belongs to, in bytes.
///
/// `None` on non-Linux, when the cgroup files cannot be read, or when nothing in
/// the cgroup chain sets a limit.
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

/// Resolve the cgroup memory limit from the contents of `/proc/self/cgroup` and a
/// file reader, so it can be unit-tested without a real cgroup.
///
/// `proc_self_cgroup` lines are `hierarchy:controllers:path`. An empty controller
/// list marks the cgroup v2 unified line (`0::/some/path`), read under
/// `/sys/fs/cgroup`; a list containing `memory` marks the v1 memory controller
/// (`7:memory:/some/path`), read under `/sys/fs/cgroup/memory`.
///
/// A cgroup's effective limit is the tightest one along its chain of ancestors —
/// SLURM commonly sets the limit on a parent scope and leaves the leaf unlimited —
/// so every level from the leaf up to the mount root is read and the minimum of
/// everything found is returned.
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

/// Read `files` in `<mount>/<cgroup_path>` and in every ancestor directory up to
/// `mount`, returning the smallest limit found.
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

/// Parse one cgroup limit file. `None` for unlimited (`max`, the v1 sentinel, or
/// anything implausibly large), for zero, and for anything unparseable.
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

/// The memory this process may plan against: the node's `MemAvailable` capped by
/// any cgroup limit that applies to it.
///
/// `None` when neither figure is available (non-Linux, restricted sandbox);
/// callers should then keep their conservative default.
pub(crate) fn available_memory_budget() -> Option<MemoryBudget> {
    budget_from(meminfo_available_bytes(), cgroup_memory_limit())
}

/// Combine the node figure and the cgroup limit: the smaller wins, and the result
/// records which one it was.
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

#[cfg(test)]
mod tests {
    use super::*;
    use std::collections::HashMap;

    const GIB: u64 = 1 << 30;

    /// Build a reader over a fixed path → contents map, standing in for the
    /// cgroup filesystem.
    fn fake_fs(files: &[(&str, &str)]) -> impl Fn(&Path) -> Option<String> {
        let map: HashMap<String, String> = files
            .iter()
            .map(|(k, v)| ((*k).to_string(), (*v).to_string()))
            .collect();
        // `Path::join` uses the platform separator, so normalise to `/`
        // before the lookup; the keys above are written Linux-style.
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
        // SLURM sets --mem on the job scope; the step leaf is left unlimited.
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
    fn real_slurm_v2_layout_finds_the_job_limit() {
        // Captured from a live `--mem=240G` job on the SLURM cluster
        // (cgroup v2, unified mount). The LEAF is `max`: SLURM puts the limit on
        // the `user` scope and the job scope above it, so a reader that only
        // consulted the leaf would conclude "unlimited" and reproduce the OOM
        // this function exists to prevent. Two ancestors carry the same 240 GiB.
        let limit = cgroup_memory_limit_from(
            "0::/system.slice/slurmstepd.scope/job_<id>/step_batch/user/task_0\n",
            fake_fs(&[
                (
                    "/sys/fs/cgroup/system.slice/slurmstepd.scope/job_<id>/step_batch/user/task_0/memory.max",
                    "max\n",
                ),
                (
                    "/sys/fs/cgroup/system.slice/slurmstepd.scope/job_<id>/step_batch/user/task_0/memory.high",
                    "max\n",
                ),
                (
                    "/sys/fs/cgroup/system.slice/slurmstepd.scope/job_<id>/step_batch/user/memory.max",
                    "257698037760\n",
                ),
                (
                    "/sys/fs/cgroup/system.slice/slurmstepd.scope/job_<id>/step_batch/memory.max",
                    "max\n",
                ),
                (
                    "/sys/fs/cgroup/system.slice/slurmstepd.scope/job_<id>/memory.max",
                    "257698037760\n",
                ),
                (
                    "/sys/fs/cgroup/system.slice/slurmstepd.scope/memory.max",
                    "max\n",
                ),
                ("/sys/fs/cgroup/system.slice/memory.max", "max\n"),
            ]),
        );
        assert_eq!(limit, Some(257_698_037_760));

        // And the budget must take that allowance, not the node's 1.2 TiB.
        let node = 1_296_452_256u64 * 1024;
        let budget = budget_from(Some(node), limit).unwrap();
        assert_eq!(budget.bytes, 257_698_037_760);
        assert_eq!(budget.source, MemorySource::CgroupLimit);
    }

    #[test]
    fn budget_is_the_minimum_of_node_memory_and_the_cgroup_limit() {
        // The reported SLURM failure: 641 GiB node, 240 GiB allowance.
        let budget = budget_from(Some(641 * GIB), Some(240 * GIB)).unwrap();
        assert_eq!(budget.bytes, 240 * GIB);
        assert_eq!(budget.source, MemorySource::CgroupLimit);

        // A cgroup limit above the node figure never inflates the budget.
        let budget = budget_from(Some(8 * GIB), Some(64 * GIB)).unwrap();
        assert_eq!(budget.bytes, 8 * GIB);
        assert_eq!(budget.source, MemorySource::NodeAvailable);

        // No cgroup limit: unchanged behaviour.
        let budget = budget_from(Some(16 * GIB), None).unwrap();
        assert_eq!(budget.bytes, 16 * GIB);
        assert_eq!(budget.source, MemorySource::NodeAvailable);

        // Neither figure readable: no budget, callers keep their default.
        assert_eq!(budget_from(None, None), None);
    }

    #[test]
    fn meminfo_available_is_parsed_as_kilobytes() {
        let text = "MemTotal:       672000000 kB\nMemFree:          100 kB\n\
                    MemAvailable:   672103936 kB\n";
        assert_eq!(parse_meminfo_available(text), Some(672_103_936 * 1024));
        assert_eq!(parse_meminfo_available("MemTotal: 1 kB\n"), None);
    }
}
