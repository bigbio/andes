use std::sync::OnceLock;

use crate::param_model::Partition;
use crate::scoring::fragment_ions::{predict_cz_ions, IonKind};
use crate::scoring::rank_scorer::RankScorer;
use crate::scoring::scored_spectrum::ScoredSpectrum;
use crate::scoring::strong_score::DENSITY_HW;
use model::mass::nominal_from;
use model::peptide::Peptide;

pub fn edge_scoring_enabled(scorer: &RankScorer) -> bool {
    scorer.param().error_scaling_factor != 0 && !scorer.param().ion_existence_table.is_empty()
}

pub fn psm_edge_score(
    scored_spec: &ScoredSpectrum,
    peptide: &Peptide,
    scorer: &RankScorer,
    charge: u8,
) -> i32 {
    if charge == 0 {
        return 0;
    }
    let n = peptide.length();
    if n < 2 {
        return 0;
    }
    if !edge_scoring_enabled(scorer) {
        return 0;
    }

    let spectrum_parent_mass = scored_spec.parent_mass();
    let peptide_nominal = peptide.nominal_residue_mass();

    let mut prefix_mass_arr: Vec<f64> = Vec::with_capacity(n + 1);
    let mut prefix_nominal_arr: Vec<i32> = Vec::with_capacity(n + 1);
    prefix_mass_arr.push(0.0);
    prefix_nominal_arr.push(0);
    let mut prefix_mass_acc = 0.0_f64;
    for s in 1..=n {
        let aa = &peptide.residues[s - 1];
        let residue_mass = aa.mass + aa.mod_.as_ref().map_or(0.0, |m| m.mass_delta);
        prefix_mass_acc += residue_mass;
        if s < n {
            prefix_mass_arr.push(prefix_mass_acc);
            prefix_nominal_arr.push(nominal_from(prefix_mass_acc));
        } else {
            prefix_mass_arr.push(prefix_mass_acc);
            prefix_nominal_arr.push(peptide_nominal);
        }
    }

    let is_prefix_main = scored_spec.main_ion_direction();
    let mut edge_total: i32 = 0;
    if !is_prefix_main {
        let nominal_peptide_mass = prefix_nominal_arr[n];
        for i in (1..n).rev() {
            let cur_nominal = nominal_peptide_mass - prefix_nominal_arr[i];
            let prev_nominal = nominal_peptide_mass - prefix_nominal_arr[i + 1];
            let theo_mass = prefix_mass_arr[i + 1] - prefix_mass_arr[i];
            edge_total += scored_spec.edge_score(
                cur_nominal,
                prev_nominal,
                theo_mass,
                scorer,
                charge,
                spectrum_parent_mass,
            );
        }
    } else {
        for i in 1..n {
            let cur_nominal = prefix_nominal_arr[i];
            let prev_nominal = prefix_nominal_arr[i - 1];
            let theo_mass = prefix_mass_arr[i] - prefix_mass_arr[i - 1];
            edge_total += scored_spec.edge_score(
                cur_nominal,
                prev_nominal,
                theo_mass,
                scorer,
                charge,
                spectrum_parent_mass,
            );
        }
    }
    edge_total
}

pub fn psm_edge_existence_facts(
    scored_spec: &ScoredSpectrum,
    peptide: &Peptide,
    scorer: &RankScorer,
    charge: u8,
) -> Vec<(Partition, u32)> {
    if charge == 0 {
        return Vec::new();
    }
    let n = peptide.length();
    if n < 2 {
        return Vec::new();
    }

    let spectrum_parent_mass = scored_spec.parent_mass();
    let peptide_nominal = peptide.nominal_residue_mass();

    let mut prefix_mass_arr: Vec<f64> = Vec::with_capacity(n + 1);
    let mut prefix_nominal_arr: Vec<i32> = Vec::with_capacity(n + 1);
    prefix_mass_arr.push(0.0);
    prefix_nominal_arr.push(0);
    let mut prefix_mass_acc = 0.0_f64;
    for s in 1..=n {
        let aa = &peptide.residues[s - 1];
        let residue_mass = aa.mass + aa.mod_.as_ref().map_or(0.0, |m| m.mass_delta);
        prefix_mass_acc += residue_mass;
        if s < n {
            prefix_mass_arr.push(prefix_mass_acc);
            prefix_nominal_arr.push(nominal_from(prefix_mass_acc));
        } else {
            prefix_mass_arr.push(prefix_mass_acc);
            prefix_nominal_arr.push(peptide_nominal);
        }
    }

    let is_prefix_main = scored_spec.main_ion_direction();
    let mut facts: Vec<(Partition, u32)> = Vec::with_capacity(n - 1);
    if !is_prefix_main {
        let nominal_peptide_mass = prefix_nominal_arr[n];
        for i in (1..n).rev() {
            let cur_nominal = nominal_peptide_mass - prefix_nominal_arr[i];
            let prev_nominal = nominal_peptide_mass - prefix_nominal_arr[i + 1];
            let (part, idx, _, _) = scored_spec.edge_existence_facts(
                cur_nominal,
                prev_nominal,
                scorer,
                charge,
                spectrum_parent_mass,
            );
            facts.push((part, idx as u32));
        }
    } else {
        for i in 1..n {
            let cur_nominal = prefix_nominal_arr[i];
            let prev_nominal = prefix_nominal_arr[i - 1];
            let (part, idx, _, _) = scored_spec.edge_existence_facts(
                cur_nominal,
                prev_nominal,
                scorer,
                charge,
                spectrum_parent_mass,
            );
            facts.push((part, idx as u32));
        }
    }
    facts
}

pub fn score_psm(
    scored_spec: &ScoredSpectrum,
    peptide: &Peptide,
    scorer: &RankScorer,
    charge: u8,
    fragment_tolerance_da: f64,
) -> f32 {
    if charge == 0 {
        return 0.0;
    }
    let n = peptide.length();
    if n < 2 {
        return 0.0;
    }

    let spectrum_parent_mass = scored_spec.parent_mass();

    let peptide_nominal = peptide.nominal_residue_mass();

    let score_losses = scorer
        .param()
        .data_type
        .activation
        .predicts_neutral_losses()
        && scorer.has_loss_tables()
        && peptide.residues.iter().any(|aa| {
            aa.mod_
                .as_ref()
                .is_some_and(|m| m.loss_class != 0 && !m.neutral_losses.is_empty())
        });
    let mut total: i32 = 0;
    let mut prefix_mass_acc = 0.0_f64;
    for s in 1..n {
        let aa = &peptide.residues[s - 1];
        let residue_mass = aa.mass + aa.mod_.as_ref().map_or(0.0, |m| m.mass_delta);
        prefix_mass_acc += residue_mass;

        let prefix_nominal = nominal_from(prefix_mass_acc);
        let suffix_nominal = peptide_nominal - prefix_nominal;

        let contribution = scored_spec
            .cached_split_score(prefix_nominal, suffix_nominal)
            .unwrap_or_else(|| {
                scored_spec.node_score(
                    prefix_nominal as f64,
                    suffix_nominal as f64,
                    scorer,
                    charge,
                    spectrum_parent_mass,
                    fragment_tolerance_da,
                )
            });
        total += contribution;

        if score_losses {
            let prefix_losses = crate::scoring::scored_spectrum::span_losses(peptide, 0..s);
            let suffix_losses = crate::scoring::scored_spectrum::span_losses(peptide, s..n);
            total += scored_spec.loss_node_score(
                prefix_nominal as f64,
                suffix_nominal as f64,
                scorer,
                charge,
                spectrum_parent_mass,
                &prefix_losses,
                &suffix_losses,
            );
        }
    }
    total as f32
}

pub fn hyperscore_psm(scored_spec: &ScoredSpectrum, peptide: &Peptide, scorer: &RankScorer) -> f32 {
    hyperscore_psm_with_matches(scored_spec, peptide, scorer).0
}

pub fn hyperscore_psm_with_matches(
    scored_spec: &ScoredSpectrum,
    peptide: &Peptide,
    scorer: &RankScorer,
) -> (f32, u32) {
    if peptide.length() < 2 {
        return (0.0, 0);
    }
    let (mut n_b, mut n_y) = (0usize, 0usize);
    for f in scored_spec
        .ion_match_facts(peptide, scorer)
        .iter()
        .filter(|f| f.rank.is_some())
    {
        if f.ion_type.is_prefix() {
            n_b += 1;
        } else {
            n_y += 1;
        }
    }
    let ln_factorial = |n: usize| -> f64 { (2..=n).map(|i| (i as f64).ln()).sum() };
    (
        (ln_factorial(n_b) + ln_factorial(n_y)) as f32,
        (n_b + n_y) as u32,
    )
}

fn cz_fix_enabled() -> bool {
    static CELL: OnceLock<bool> = OnceLock::new();
    *CELL.get_or_init(|| true)
}

#[derive(Debug, Clone, Copy, PartialEq, Eq, Default)]
pub struct CzSettings {
    pub zmax_override: Option<u8>,
}

static CZ_SETTINGS: OnceLock<CzSettings> = OnceLock::new();

pub fn init_cz_settings(settings: CzSettings) {
    let _ = CZ_SETTINGS.set(settings);
}

#[inline]
fn cz_settings() -> CzSettings {
    CZ_SETTINGS.get().copied().unwrap_or_default()
}

fn cz_zmax_override() -> Option<u8> {
    use std::sync::OnceLock;
    static CELL: OnceLock<Option<u8>> = OnceLock::new();
    *CELL.get_or_init(|| cz_settings().zmax_override.filter(|&z| z >= 1))
}

pub fn cz_effective_zmax_for(default_zmax: u8, apply_deconvolution: bool) -> u8 {
    match cz_zmax_override() {
        Some(z) => default_zmax.min(z),
        None if apply_deconvolution => default_zmax.min(1),
        None => default_zmax,
    }
}

pub fn cz_hyperscore_psm(
    scored_spec: &ScoredSpectrum,
    peptide: &Peptide,
    glycan_mass: f64,
    glycosite: usize,
    max_frag_charge: u8,
    tol_ppm: f64,
) -> f32 {
    let n = peptide.length();
    if n < 2 || max_frag_charge == 0 {
        return 0.0;
    }
    let fix = cz_fix_enabled();
    let mut c_hit = vec![false; n];
    let mut z_hit = vec![false; n];
    for ion in predict_cz_ions(
        peptide,
        1..=cz_effective_zmax_for(max_frag_charge, scored_spec.is_deconvoluted()),
        glycan_mass,
        glycosite,
    ) {
        let tol_da = if fix {
            (ion.mz * tol_ppm * 1e-6).max(0.01)
        } else {
            0.5
        };
        if scored_spec.nearest_peak_full(ion.mz, tol_da).is_some() {
            let pos = ion.position as usize;
            match ion.kind {
                IonKind::C if pos < n => c_hit[pos] = true,
                IonKind::Z if pos < n => z_hit[pos] = true,
                _ => {}
            }
        }
    }
    let n_c = c_hit.iter().filter(|&&h| h).count();
    let n_z = z_hit.iter().filter(|&&h| h).count();
    let ln_factorial = |m: usize| -> f64 { (2..=m).map(|i| (i as f64).ln()).sum() };
    let raw = ln_factorial(n_c) + ln_factorial(n_z);
    if fix {
        let max_sites = n.saturating_sub(1).max(1);
        let norm = 2.0 * ln_factorial(max_sites);
        if norm > 0.0 {
            (10.0 * raw / norm) as f32
        } else {
            0.0
        }
    } else {
        raw as f32
    }
}

pub fn cz_matched_intensity_frac(
    scored_spec: &ScoredSpectrum,
    peptide: &Peptide,
    glycan_mass: f64,
    glycosite: usize,
    max_frag_charge: u8,
    tol_ppm: f64,
) -> f32 {
    let n = peptide.length();
    if n < 2 || max_frag_charge == 0 {
        return 0.0;
    }
    let (peaks, _) = scored_spec.active_peaks_and_ranks();
    let base = peaks
        .iter()
        .map(|&(_, i)| i)
        .fold(0.0f32, f32::max)
        .max(1e-9);
    let mut c_int = vec![0.0f32; n];
    let mut z_int = vec![0.0f32; n];
    for ion in predict_cz_ions(
        peptide,
        1..=cz_effective_zmax_for(max_frag_charge, scored_spec.is_deconvoluted()),
        glycan_mass,
        glycosite,
    ) {
        let tol_da = (ion.mz * tol_ppm * 1e-6).max(0.01);
        if let Some((_, intensity, _)) = scored_spec.nearest_peak_full(ion.mz, tol_da) {
            let pos = ion.position as usize;
            match ion.kind {
                IonKind::C if pos < n => c_int[pos] = c_int[pos].max(intensity),
                IonKind::Z if pos < n => z_int[pos] = z_int[pos].max(intensity),
                _ => {}
            }
        }
    }
    let summed: f32 = c_int.iter().chain(z_int.iter()).sum();
    summed / base
}

pub fn cz_structure_features(
    scored_spec: &ScoredSpectrum,
    peptide: &Peptide,
    glycan_mass: f64,
    glycosite: usize,
    max_frag_charge: u8,
    tol_ppm: f64,
) -> (f32, f32) {
    let n = peptide.length();
    if n < 2 || max_frag_charge == 0 {
        return (0.0, 0.0);
    }
    let zmax = cz_effective_zmax_for(max_frag_charge.min(3), scored_spec.is_deconvoluted());
    let mut sum_p = 0.0f64;
    let mut sum_matched_p = 0.0f64;
    let mut sum_chance_llr = 0.0f64;
    for ion in predict_cz_ions(peptide, 1..=zmax, glycan_mass, glycosite) {
        let pos = ion.position as usize;
        let spanning = match ion.kind {
            IonKind::C => glycosite < pos,
            IonKind::Z => glycosite >= n.saturating_sub(pos),
            _ => false,
        };
        let charge_w = if ion.charge <= 2 { 1.0 } else { 0.4 };
        let p = if spanning { 0.15 } else { 1.0 } * charge_w;
        sum_p += p;
        let tol_da = (ion.mz * tol_ppm * 1e-6).max(0.01);
        if scored_spec.nearest_peak_full(ion.mz, tol_da).is_some() {
            sum_matched_p += p;
            let rho = scored_spec.local_peak_density(ion.mz, DENSITY_HW);
            let p_chance = (rho * 2.0 * tol_da).clamp(1e-12, 1.0);
            let surprise = (-p_chance.ln()).max(0.0);
            sum_chance_llr += p * surprise;
        }
    }
    let explained = if sum_p > 0.0 {
        (sum_matched_p / sum_p) as f32
    } else {
        0.0
    };
    (explained, sum_chance_llr as f32)
}

pub fn score_psm_float(
    scored_spec: &ScoredSpectrum,
    peptide: &Peptide,
    scorer: &RankScorer,
    charge: u8,
    fragment_tolerance_da: f64,
) -> f32 {
    if charge == 0 {
        return 0.0;
    }
    let n = peptide.length();
    if n < 2 {
        return 0.0;
    }
    let spectrum_parent_mass = scored_spec.parent_mass();
    let peptide_nominal = peptide.nominal_residue_mass();
    let mut total = 0.0_f32;
    let mut prefix_mass_acc = 0.0_f64;
    for s in 1..n {
        let aa = &peptide.residues[s - 1];
        let residue_mass = aa.mass + aa.mod_.as_ref().map_or(0.0, |m| m.mass_delta);
        prefix_mass_acc += residue_mass;
        let prefix_nominal = nominal_from(prefix_mass_acc);
        let suffix_nominal = peptide_nominal - prefix_nominal;
        let contribution = scored_spec
            .cached_split_score_f32(prefix_nominal, suffix_nominal)
            .unwrap_or_else(|| {
                scored_spec.node_score_f32(
                    prefix_nominal as f64,
                    suffix_nominal as f64,
                    scorer,
                    charge,
                    spectrum_parent_mass,
                    fragment_tolerance_da,
                )
            });
        total += contribution;
    }
    total
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::param_model::{FragmentOffsetFrequency, IonType, Param, Partition, SpecDataType};
    use crate::scoring::rank_scorer::RankScorer;
    use crate::scoring::scored_spectrum::ScoredSpectrum;
    use crate::testutil::tiny_param;
    use model::amino_acid::AminoAcid;
    use model::peptide::Peptide;
    use model::spectrum::Spectrum;
    use rustc_hash::FxHashMap;

    fn pep(seq: &[u8]) -> Peptide {
        let residues: Vec<AminoAcid> = seq
            .iter()
            .map(|&r| AminoAcid::standard(r).unwrap())
            .collect();
        Peptide::new(residues, b'_', b'-')
    }

    fn empty_spectrum(title: &str) -> Spectrum {
        Spectrum {
            title: title.into(),
            precursor_mz: 0.0,
            precursor_intensity: None,
            precursor_charge: Some(2),
            rt_seconds: None,
            scan: None,
            peaks: vec![],
            activation_method: None,
            isolation_lower_offset: None,
            isolation_upper_offset: None,
            isolation_target_mz: None,
        }
    }

    fn any_mass_param() -> Param {
        use model::activation::ActivationMethod;
        use model::instrument::InstrumentType;
        use model::protocol::Protocol;

        let part = Partition {
            charge: 2,
            parent_mass: 0.0,
            seg_num: 0,
        };
        let prefix_ion = IonType::Prefix {
            charge: 1,
            offset_bits: 0.0_f32.to_bits(),
            loss_class: 0,
        };
        let noise_ion = IonType::Noise;

        let ion_freqs = vec![0.6_f32, 0.3, 0.05, 0.001];
        let noise_freqs = vec![0.1_f32, 0.2, 0.3, 0.4];

        let mut ion_table: FxHashMap<IonType, Vec<f32>> = FxHashMap::default();
        ion_table.insert(prefix_ion, ion_freqs);
        ion_table.insert(noise_ion, noise_freqs);

        let mut rank_dist_table: FxHashMap<Partition, FxHashMap<IonType, Vec<f32>>> =
            FxHashMap::default();
        rank_dist_table.insert(part, ion_table);

        let mut frag_off_table = FxHashMap::default();
        frag_off_table.insert(
            part,
            vec![FragmentOffsetFrequency {
                ion_type: prefix_ion,
                frequency: 0.7,
            }],
        );

        let mut p = Param {
            version: 10001,
            data_type: SpecDataType {
                activation: ActivationMethod::HCD,
                instrument: InstrumentType::QExactive,
                enzyme: None,
                protocol: Protocol::Automatic,
            },
            mme: model::tolerance::Tolerance::Da(0.2),
            apply_deconvolution: false,
            deconvolution_error_tolerance: 0.0,
            charge_hist: vec![(2, 100)],
            min_charge: 2,
            max_charge: 2,
            num_segments: 1,
            partitions: vec![part],
            num_precursor_off: 0,
            precursor_off_map: FxHashMap::default(),
            frag_off_table,
            max_rank: 3,
            rank_dist_table,
            error_scaling_factor: 0,
            ion_err_dist_table: FxHashMap::default(),
            noise_err_dist_table: FxHashMap::default(),
            ion_existence_table: FxHashMap::default(),
            partition_ion_types_cache: FxHashMap::default(),
            gbdt_peak_model: None,
            frag_intensity_model: None,
            rich_ion_model: None,
        };
        p.rebuild_cache();
        p
    }

    fn loss_table_param() -> Param {
        let part = Partition {
            charge: 2,
            parent_mass: 0.0,
            seg_num: 0,
        };
        let loss_ion = IonType::Prefix {
            charge: 1,
            offset_bits: 0.0_f32.to_bits(),
            loss_class: 1,
        };
        let mut p = any_mass_param();
        p.rank_dist_table
            .get_mut(&part)
            .unwrap()
            .insert(loss_ion, vec![0.6_f32, 0.3, 0.05, 0.001]);
        p
    }

    fn loss_peptide_twoheavy() -> Peptide {
        use model::modification::{ModLocation, Modification, ResidueSpec};
        let m = Modification {
            name: "L".into(),
            mass_delta: 0.0,
            residue: ResidueSpec::Specific(b'T'),
            location: ModLocation::Anywhere,
            fixed: false,
            accession: None,
            neutral_losses: vec![162.0528],
            loss_class: 1,
            annotated_site_mass: None,
            blocks_cleavage: false,
            parts: Vec::new(),
            diagnostic_ions: Vec::new(),
        };
        let t = AminoAcid::standard(b'T').unwrap().with_mod(m);
        let w1 = AminoAcid::standard(b'W').unwrap();
        let w2 = AminoAcid::standard(b'W').unwrap();
        Peptide::new(vec![t, w1, w2], b'_', b'-')
    }

    #[test]
    fn score_psm_adds_neutral_loss_contribution_and_is_inert_otherwise() {
        let _part = Partition {
            charge: 2,
            parent_mass: 0.0,
            seg_num: 0,
        };
        let loss_ion = IonType::Prefix {
            charge: 1,
            offset_bits: 0.0_f32.to_bits(),
            loss_class: 1,
        };

        let loss_pep = loss_peptide_twoheavy();
        let prefix_nominal = nominal_from(loss_pep.residues[0].mass + loss_pep.residues[1].mass);
        let loss_mz = loss_ion.mz(prefix_nominal as f64) - 162.0528;
        assert!(
            loss_mz > 0.0,
            "test setup: loss m/z must be positive, got {loss_mz}"
        );

        let s = Spectrum {
            title: "loss".into(),
            precursor_mz: 800.0,
            precursor_intensity: None,
            precursor_charge: Some(2),
            rt_seconds: None,
            scan: None,
            peaks: vec![(loss_mz, 100.0)],
            activation_method: None,
            isolation_lower_offset: None,
            isolation_upper_offset: None,
            isolation_target_mz: None,
        };

        let loss_scorer = RankScorer::new(&loss_table_param());
        let plain_scorer = RankScorer::new(&any_mass_param());
        let ss_loss = ScoredSpectrum::new(&s, &loss_scorer, 2);
        let ss_plain = ScoredSpectrum::new(&s, &plain_scorer, 2);

        let with_loss = score_psm(&ss_loss, &loss_pep, &loss_scorer, 2, 0.5);
        let without = score_psm(&ss_plain, &loss_pep, &plain_scorer, 2, 0.5);
        assert!(
            with_loss > without,
            "loss model must add positive loss-ion score: with={with_loss} without={without}"
        );

        let plain_pep = pep(b"AAA");
        let sp = Spectrum {
            title: "plain".into(),
            precursor_mz: 800.0,
            precursor_intensity: None,
            precursor_charge: Some(2),
            rt_seconds: None,
            scan: None,
            peaks: vec![(200.0, 50.0), (300.0, 80.0)],
            activation_method: None,
            isolation_lower_offset: None,
            isolation_upper_offset: None,
            isolation_target_mz: None,
        };
        let a = score_psm(
            &ScoredSpectrum::new(&sp, &loss_scorer, 2),
            &plain_pep,
            &loss_scorer,
            2,
            0.5,
        );
        let b = score_psm(
            &ScoredSpectrum::new(&sp, &plain_scorer, 2),
            &plain_pep,
            &plain_scorer,
            2,
            0.5,
        );
        assert_eq!(
            a, b,
            "non-loss peptide must be byte-identical across loss/no-loss models"
        );
    }

    #[test]
    fn empty_spectrum_returns_non_positive_score() {
        let peptide = pep(b"AGR");
        let spec = empty_spectrum("empty");
        let scored = ScoredSpectrum::new_without_filtering(&spec);
        let param = any_mass_param();
        let scorer = RankScorer::new(&param);
        let s = score_psm(&scored, &peptide, &scorer, 2, 0.2);
        assert!(s <= 0.0, "score should be ≤ 0 on empty spectrum, got {s}");
    }

    #[test]
    fn perfect_match_yields_positive_score() {
        let peptide = pep(b"AGR");
        let param = any_mass_param();

        let b_ion = IonType::Prefix {
            charge: 1,
            offset_bits: 0.0_f32.to_bits(),
            loss_class: 0,
        };
        let mut prefix_acc = 0.0_f64;
        let mut peaks = Vec::new();
        for s in 1..peptide.length() {
            let aa = &peptide.residues[s - 1];
            prefix_acc += aa.mass;
            let nom = model::mass::nominal_from(prefix_acc) as f64;
            let mz = b_ion.mz(nom);
            peaks.push((mz, 1000.0_f32 / s as f32));
        }
        peaks.sort_by(|a, b| a.0.partial_cmp(&b.0).unwrap());

        let spec = Spectrum {
            title: "match".into(),
            precursor_mz: 500.0,
            precursor_intensity: None,
            precursor_charge: Some(2),
            rt_seconds: None,
            scan: None,
            peaks,
            activation_method: None,
            isolation_lower_offset: None,
            isolation_upper_offset: None,
            isolation_target_mz: None,
        };
        let scored = ScoredSpectrum::new_without_filtering(&spec);
        let scorer = RankScorer::new(&param);
        let s = score_psm(&scored, &peptide, &scorer, 2, 0.2);
        assert!(
            s > 0.0,
            "score with matched b-ions should be positive, got {s}"
        );
    }

    #[test]
    fn perfect_match_outscores_empty_spectrum() {
        let peptide = pep(b"AGR");
        let param = any_mass_param();

        let b_ion = IonType::Prefix {
            charge: 1,
            offset_bits: 0.0_f32.to_bits(),
            loss_class: 0,
        };
        let mut prefix_acc = 0.0_f64;
        let mut match_peaks = Vec::new();
        for s in 1..peptide.length() {
            let aa = &peptide.residues[s - 1];
            prefix_acc += aa.mass;
            let nom = model::mass::nominal_from(prefix_acc) as f64;
            let mz = b_ion.mz(nom);
            match_peaks.push((mz, 1000.0_f32));
        }
        match_peaks.sort_by(|a, b| a.0.partial_cmp(&b.0).unwrap());

        let match_spec = Spectrum {
            title: "match".into(),
            precursor_mz: 500.0,
            precursor_intensity: None,
            precursor_charge: Some(2),
            rt_seconds: None,
            scan: None,
            peaks: match_peaks,
            activation_method: None,
            isolation_lower_offset: None,
            isolation_upper_offset: None,
            isolation_target_mz: None,
        };

        let scorer = RankScorer::new(&param);
        let scored_match = ScoredSpectrum::new_without_filtering(&match_spec);
        let empty_spec = empty_spectrum("empty");
        let scored_empty = ScoredSpectrum::new_without_filtering(&empty_spec);
        let s_match = score_psm(&scored_match, &peptide, &scorer, 2, 0.2);
        let s_empty = score_psm(&scored_empty, &peptide, &scorer, 2, 0.2);
        assert!(
            s_match > s_empty,
            "matched spectrum ({s_match}) should outscore empty ({s_empty})"
        );
    }

    #[test]
    fn score_psm_matches_sum_of_node_scores_across_splits() {
        use model::amino_acid::AminoAcid;
        use model::mass::nominal_from;

        let peptide = pep(b"AGR");
        let param = tiny_param();
        let scorer = RankScorer::new(&param);

        let empty_spec = empty_spectrum("empty");
        let scored = ScoredSpectrum::new_without_filtering(&empty_spec);

        let parent_mass = peptide.mass();
        let peptide_nominal = peptide.nominal_residue_mass();
        let charge = 2u8;
        let tolerance_da = 0.05;

        let mut manual_total: i32 = 0;
        let mut prefix_acc = 0.0_f64;
        for s in 1..peptide.length() {
            let aa: &AminoAcid = &peptide.residues[s - 1];
            prefix_acc += aa.mass + aa.mod_.as_ref().map_or(0.0, |m| m.mass_delta);
            let pref = nominal_from(prefix_acc);
            let suf = peptide_nominal - pref;
            manual_total += scored.node_score(
                pref as f64,
                suf as f64,
                &scorer,
                charge,
                parent_mass,
                tolerance_da,
            );
        }

        let computed = score_psm(&scored, &peptide, &scorer, charge, tolerance_da);
        assert_eq!(
            computed as i32, manual_total,
            "score_psm ({computed}) should equal manual split sum ({manual_total})"
        );
    }

    #[test]
    fn charge_zero_returns_zero() {
        let peptide = pep(b"AGR");
        let param = tiny_param();
        let scorer = RankScorer::new(&param);
        let spec = empty_spectrum("empty");
        let scored = ScoredSpectrum::new_without_filtering(&spec);
        assert_eq!(score_psm(&scored, &peptide, &scorer, 0, 0.1), 0.0);
    }

    #[test]
    fn single_residue_peptide_returns_zero() {
        let peptide = pep(b"A");
        let param = tiny_param();
        let scorer = RankScorer::new(&param);
        let spec = empty_spectrum("empty");
        let scored = ScoredSpectrum::new_without_filtering(&spec);
        assert_eq!(score_psm(&scored, &peptide, &scorer, 2, 0.1), 0.0);
    }
}
