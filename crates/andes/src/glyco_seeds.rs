//! Extract 1%-FDR Pass-1 glyco seeds (target + decoy) from a rescored PIN.
use andes_glyco::crossspectrum::Seed;

#[derive(Debug, Clone)]
pub struct SeedRow {
    pub scan: u32,
    pub is_decoy: bool,
    pub q_value: f64,
    pub score: f64,
}

/// Pull the integer scan from a glyco SpecId like
/// "controllerType=0 controllerNumber=1 scan=3000_glyco_3000_1". Public so the
/// driver (Task 8d) can join the native rescorer's `spec_id` strings back to a
/// scan number without duplicating this parsing logic.
pub fn extract_scan(id: &str) -> Option<u32> {
    let after = id.split("scan=").nth(1)?;
    let digits: String = after.chars().take_while(|c| c.is_ascii_digit()).collect();
    digits.parse().ok()
}

/// Keep rows at `q <= q_threshold` that map to a Pass-1 backbone, as `Seed`s.
pub fn seeds_at_fdr(
    rows: &[SeedRow],
    q_threshold: f64,
    lookup: impl Fn(u32) -> Option<(u32, f64, Option<f64>)>,
) -> Vec<Seed> {
    let mut seeds: Vec<Seed> = rows
        .iter()
        .filter(|r| r.q_value <= q_threshold)
        .filter_map(|r| {
            let (peptide_idx, backbone_mass, rt_seconds) = lookup(r.scan)?;
            Some(Seed {
                scan: r.scan,
                peptide_idx,
                backbone_mass,
                rt_seconds,
                seed_score: r.score,
                is_decoy: r.is_decoy,
            })
        })
        .collect();
    // Deterministic order (target/decoy seeds propagate identically downstream).
    seeds.sort_by(|a, b| a.scan.cmp(&b.scan).then(a.backbone_mass.total_cmp(&b.backbone_mass)));
    seeds
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn seeds_at_fdr_keeps_only_q_below_threshold_and_maps_backbone() {
        let rows = vec![
            SeedRow { scan: 10, is_decoy: false, q_value: 0.004, score: 3.1 },
            SeedRow { scan: 11, is_decoy: false, q_value: 0.05, score: 0.4 }, // fails FDR
            SeedRow { scan: 12, is_decoy: true, q_value: 0.008, score: 2.0 }, // decoy seed
        ];
        let lookup = |scan: u32| match scan {
            10 => Some((100u32, 1500.0f64, Some(900.0f64))),
            12 => Some((200u32, 1800.0f64, Some(905.0f64))),
            _ => None,
        };
        let seeds = seeds_at_fdr(&rows, 0.01, lookup);
        assert_eq!(seeds.len(), 2, "only q<=0.01 rows with a backbone: {seeds:?}");
        assert!(seeds.iter().any(|s| s.scan == 10 && !s.is_decoy && (s.backbone_mass - 1500.0).abs() < 1e-6));
        assert!(seeds.iter().any(|s| s.scan == 12 && s.is_decoy));
        assert!(!seeds.iter().any(|s| s.scan == 11));
    }
}
