//! `quant library`: join the runs' identifications, align retention times and write the
//! precursor library with its gates.

use std::collections::{BTreeMap, HashMap};
use std::path::{Path, PathBuf};

use model::isomers::bundled_unimod_mass;
use model::IsomerClass;
use output::library::{write_library, StoredLibrary};
use output::quantities::runs::{open_run_quant_table, read_run_quant_mods, RunQuantBatches};
use quant::align::AlignParams;
use quant::library::{is_donor_psm, Library, LibraryBuilder, PrecursorMeta, PsmIn};
use rayon::prelude::*;

use super::extract::{split_key, BestForm, DonorPsms};
use super::scores::{RunScores, ScoreLoader, ScoreSet, SpecIdKeys, MAX_JOIN_MISMATCH};
use super::{artifact, ENVELOPE_LEN, MS1_SUFFIX, TABLE_SUFFIX};

const META_COLUMNS: [&str; 7] = [
    "precursor_key",
    "is_decoy",
    "mz",
    "charge",
    "envelope",
    "proteins",
    "isotope_share",
];
const JOIN_COLUMNS: [&str; 7] = [
    "spec_id",
    "precursor_key",
    "peptidoform",
    "is_decoy",
    "rt_s",
    "faims_cv",
    "ion_mobility",
];

pub(super) struct KeyMeta {
    has_target: bool,
    pub(super) run: usize,
    pub(super) mz: f64,
    pub(super) charge: u8,
    envelope: Result<[f32; ENVELOPE_LEN], usize>,
    pub(super) proteins: String,
    isotope_share: f32,
}

fn envelope_of(v: &[f32]) -> Result<[f32; ENVELOPE_LEN], usize> {
    v.try_into().map_err(|_| v.len())
}

impl KeyMeta {
    fn of(row: &output::quantities::runs::RunPsmRow, run: usize) -> Self {
        Self {
            has_target: !row.is_decoy,
            run,
            mz: row.mz,
            charge: row.charge,
            envelope: envelope_of(&row.envelope),
            proteins: row.proteins.clone(),
            isotope_share: row.isotope_share,
        }
    }
}

const PARALLEL_RUNS: usize = 32;

#[derive(Default)]
pub(super) struct RunJoin {
    pub(super) rows: usize,
    joined: usize,
    mismatched: usize,
    unscored_t: usize,
    unscored_d: usize,
    shadowed: usize,
    pub(super) psms: Vec<(PsmIn, u32, Option<Box<str>>)>,
}

#[derive(Default)]
pub(super) struct RunKeys {
    pub(super) rows: usize,
    pub(super) keys: HashMap<String, KeyMeta>,
    spec: Vec<u128>,
}

fn align_params(dir: &Path, stems: &[String]) -> AlignParams {
    let mut ap = AlignParams::default();
    ap.run_min_spacing_s = stems
        .iter()
        .map(|stem| {
            let path = artifact(dir, stem, MS1_SUFFIX);
            let spacing = output::ms1::read_ms1_scan_spacing(&path).ok().flatten();
            let floor = spacing.map_or(ap.min_spacing_s, |s| s / 2.0);
            match spacing {
                Some(s) => log::debug!(
                    target: "mbr",
                    "run {stem}: MS1 scan spacing {s:.2} s, RT residual-width floor {floor:.2} s"
                ),
                None => log::debug!(
                    target: "mbr",
                    "run {stem}: MS1 scan spacing unknown, RT residual-width floor {floor:.2} s"
                ),
            }
            floor
        })
        .collect();
    ap
}

fn log_floor_hits(library: &Library, ap: &AlignParams, stems: &[String]) {
    for (r, knots) in library.alignment.residual_sd.iter().enumerate() {
        let Some(knots) = knots else { continue };
        let floor = ap.floor_s(r);
        let at_floor = knots.iter().filter(|k| k.1 <= floor * (1.0 + 1e-9)).count();
        log::debug!(
            target: "mbr",
            "run {}: {at_floor} of {} RT residual-width knots at the {floor:.2} s floor",
            stems.get(r).map_or("?", String::as_str),
            knots.len()
        );
    }
}

fn log_alignment(library: &Library, ap: &AlignParams, stems: &[String]) {
    let a = &library.alignment;
    let no_consensus = library
        .entries
        .iter()
        .filter(|e| !e.consensus_rt.is_finite())
        .count();
    if no_consensus > 0 {
        log::info!(
            target: "library",
            "{no_consensus} entries have no consensus RT (no observation in an aligned run) and are not transferred"
        );
    }
    if a.n_components() > 1 {
        let mut groups: Vec<Vec<&str>> = vec![Vec::new(); a.n_components()];
        for (r, stem) in stems.iter().enumerate() {
            if let Some(c) = a.component(r as u16) {
                groups[c as usize].push(stem);
            }
        }
        log::info!(
            target: "alignment",
            "the aligned runs form {} components (runs joined when they share at least {} precursors); transfers stay inside a component: {}",
            groups.len(),
            ap.min_shared,
            groups
                .iter()
                .enumerate()
                .map(|(c, g)| format!("component {c}: {}", g.join(", ")))
                .collect::<Vec<_>>()
                .join("; ")
        );
    }
    for (r, stem) in stems.iter().enumerate() {
        if a.maps[r].is_none() {
            let isolated = a.anchors[r] >= ap.min_anchors && a.coverage[r] >= ap.min_coverage;
            log::warn!(
                target: "mbr",
                "run {stem} is unaligned ({} anchors, {:.2} of the consensus RT range{}): it neither gives nor receives transfers and keeps its identified quantities",
                a.anchors[r],
                a.coverage[r],
                if isolated {
                    format!(
                        "; shares fewer than {} precursors with every other aligned run",
                        ap.min_shared
                    )
                } else {
                    String::new()
                }
            );
        }
    }
}

fn check_library_gates(max_q: f64, donor_q: f64, donor_pep: Option<f64>) -> Result<(), String> {
    if !(max_q > 0.0 && max_q <= 1.0) {
        return Err(format!("library entry gate {max_q} is outside (0, 1]"));
    }
    if !(donor_q > 0.0 && donor_q <= 1.0) {
        return Err(format!(
            "donor PSM q-value threshold {donor_q} is outside (0, 1]"
        ));
    }
    if let Some(c) = donor_pep.filter(|c| !(0.0..=1.0).contains(c)) {
        return Err(format!("donor PSM PEP cap {c} is outside [0, 1]"));
    }
    Ok(())
}

pub(crate) fn build_and_write_library(
    dir: &Path,
    stems: &[String],
    scores: ScoreLoader<'_>,
    max_q: f64,
    donor_q: f64,
    donor_pep: Option<f64>,
) -> Result<StoredLibrary, String> {
    check_library_gates(max_q, donor_q, donor_pep)?;
    let paths: Vec<PathBuf> = stems
        .iter()
        .map(|s| artifact(dir, s, TABLE_SUFFIX))
        .collect();
    let tables = RunTables {
        paths: &paths,
        stems,
    };
    let started = std::time::Instant::now();
    let merged = matches!(scores, ScoreLoader::Merged(_));
    let all_runs: Vec<usize> = (0..stems.len()).collect();
    let merge = read_keys(&tables, merged, &all_runs)?;
    let KeyMerge {
        mut keys,
        seen,
        shared,
        n_rows,
    } = merge;
    let n_spec_ids = seen.len();
    drop(seen);
    let (merged_scores, per_run) = match scores {
        ScoreLoader::Merged(load) => {
            let set = load(n_spec_ids)?;
            if let Some((a, b, row)) = shared {
                let id = tables
                    .open(b, &["spec_id"])?
                    .flat_map(|batch| batch.into_iter().flatten())
                    .nth(row)
                    .map_or_else(String::new, |r| r.spec_id);
                return Err(format!(
                    "SpecId {id:?} is in runs {} and {}, so one merged Percolator result cannot be joined; pass per-run results: --psms <stem>=<target.tsv> --decoy-psms <stem>=<decoy.tsv> for every run",
                    stems[a], stems[b]
                ));
            }
            (Some(set), None)
        }
        ScoreLoader::PerRun(load) => (None, Some(load)),
    };
    let mod_masses = quant_mod_masses(&paths)?;
    let resolve = |tag: &str| {
        mod_masses
            .get(tag)
            .copied()
            .or_else(|| bundled_unimod_mass(tag))
    };
    let site_keys: Vec<(String, KeyMeta)> = std::mem::take(&mut keys).into_iter().collect();
    let IsomerTable {
        mut meta,
        proteins,
        target_key,
        forms,
    } = isomer_classes(&paths, &site_keys, resolve)?;
    let id_of = site_ids(&site_keys, &forms);
    let mut builder = LibraryBuilder::new(meta.len(), stems.len(), donor_q, donor_pep);
    let mut donors = DonorPsms::new(meta.len());
    let mut best = vec![None::<BestForm>; meta.len()];
    let (mut unscored_t, mut unscored_d, mut shadowed, mut mismatched) = (0usize, 0, 0, 0);
    let (mut n_psms, mut donor_psms, mut n_joined) = (0usize, 0usize, 0usize);
    let inputs = JoinInputs {
        tables: &tables,
        merged_scores: &merged_scores,
        per_run: &per_run,
        id_of: &id_of,
        target_key: &target_key,
        resolve,
        donor_q,
        donor_pep,
    };
    for group in all_runs.chunks(PARALLEL_RUNS) {
        let joined: Vec<Result<RunJoin, String>> =
            group.par_iter().map(|&r| inputs.join(r)).collect();
        for run in joined {
            let run = run?;
            unscored_t += run.unscored_t;
            unscored_d += run.unscored_d;
            shadowed += run.shadowed;
            mismatched += run.mismatched;
            n_joined += run.joined;
            for (p, site, spec_id) in run.psms {
                n_psms += 1;
                let donor = spec_id.is_some();
                donor_psms += usize::from(donor && !p.is_decoy);
                builder.add(&p);
                let form = BestForm {
                    donor,
                    score: p.score,
                    rt: p.rt,
                    run: p.run,
                    site,
                };
                let slot = &mut best[p.precursor as usize];
                if slot.as_ref().is_none_or(|b| form.beats(b)) {
                    *slot = Some(form);
                }
                if let Some(spec_id) = spec_id {
                    donors.add(&p, &spec_id);
                }
            }
        }
    }
    let joined_s = started.elapsed().as_secs_f64();
    if n_psms == 0 && n_rows > 0 {
        return Err(format!(
            "none of the {n_rows} PSMs of the run tables has a rescoring score; \
             the PSM files must come from the PINs of the same searches"
        ));
    }
    drop(id_of);
    drop(merged_scores);
    let (multi, n_forms) = rekey_by_best_form(&mut meta, &best, &site_keys, &target_key, &forms);
    drop(forms);
    drop(site_keys);
    log::info!(
        target: "mbr",
        "positional isomers: {} target precursors hold {n_forms} site forms; {multi} precursors join several site forms and are keyed by the site form of their best PSM",
        target_key.iter().filter(|&&t| t).count()
    );
    let donor_gate = donor_gate_text(donor_q, donor_pep);
    let ap = align_params(dir, stems);
    let aligning = std::time::Instant::now();
    let library = builder.finish(meta, max_q, &ap);
    log_floor_hits(&library, &ap, stems);
    let align_s = aligning.elapsed().as_secs_f64();
    let spectra_started = std::time::Instant::now();
    let (spectra, donor_scans, sc) = donors.finish(dir, stems, &library)?;
    log::debug!(
        target: "library",
        "run tables and scores {joined_s:.1} s, alignment {align_s:.1} s, spectra {:.1} s",
        spectra_started.elapsed().as_secs_f64()
    );
    log::info!(
        target: "library",
        "join check: {mismatched} of {n_joined} PSMs joined on SpecId ({:.4}%) disagree with the rescoring result on the peptide's residues, modification masses or target/decoy label (a run fails above {}%)",
        100.0 * mismatched as f64 / n_joined.max(1) as f64,
        100.0 * MAX_JOIN_MISMATCH
    );
    log::info!(
        target: "library",
        "{} runs, {} PSMs scored, {donor_psms} target PSMs at {donor_gate} identify their precursor in their run; skipped: {} PSMs whose SpecId joins a rescoring PSM of another peptide, modification composition or target/decoy label, {} target and {} decoy PSMs without a rescoring score, {} decoy PSMs sharing a target's key; {} precursors, {} entries at the internal entry gate q<={max_q} (precursor-level, not an identification FDR) with a PSM at {donor_gate}, {} of them identity decoys (decoy precursors with decoy PSMs at {donor_gate}, measured as targets and never reported)",
        stems.len(),
        n_psms,
        mismatched,
        unscored_t,
        unscored_d,
        shadowed,
        library.meta.len(),
        library.entries.len(),
        library.entries.iter().filter(|e| e.identity_decoy).count()
    );
    log::info!(
        target: "library",
        "fragment spectra (best donor PSM, at most {} b/y fragments of {}+ residues matched on its scan) for {} of {} entries; none for {} whose scan is not in the MS2 table, {} without fragment masses, {} with fewer than {} matched fragments",
        quant::ms2::LIBRARY_FRAGMENTS,
        quant::ms2::MIN_ION_LENGTH,
        sc.with_spectrum,
        library.entries.len(),
        sc.no_scan,
        sc.no_masses,
        sc.too_few,
        quant::ms2::MIN_LIBRARY_FRAGMENTS
    );
    log_alignment(&library, &ap, stems);
    let stored = StoredLibrary {
        runs: stems.to_vec(),
        max_q,
        donor_q,
        donor_pep,
        library,
        proteins,
        spectra,
        donor_scans,
    };
    write_library(dir, &stored)?;
    log::info!(
        target: "mbr",
        "wrote {} and {}",
        dir.join(output::library::LIBRARY_FILE).display(),
        dir.join(output::library::ALIGNMENT_FILE).display()
    );
    Ok(stored)
}

/// The runs' quant tables, opened by run index and checked to be the table of that run.
struct RunTables<'a> {
    paths: &'a [PathBuf],
    stems: &'a [String],
}

impl RunTables<'_> {
    fn open(&self, r: usize, columns: &[&str]) -> Result<RunQuantBatches, String> {
        let batches = open_run_quant_table(&self.paths[r], columns)?;
        if batches.run() != self.stems[r] {
            return Err(format!(
                "{}: the table is of run {:?}, not {:?}",
                self.paths[r].display(),
                batches.run(),
                self.stems[r]
            ));
        }
        Ok(batches)
    }
}

/// Every run's precursor keys, and its SpecIds when the scores are one merged result, read
/// `PARALLEL_RUNS` runs at a time.
fn read_keys(tables: &RunTables<'_>, merged: bool, all_runs: &[usize]) -> Result<KeyMerge, String> {
    let mut merge = KeyMerge::default();
    let spec_keys = SpecIdKeys::default();
    let columns: Vec<&str> = META_COLUMNS
        .iter()
        .copied()
        .chain(merged.then_some("spec_id"))
        .collect();
    for group in all_runs.chunks(PARALLEL_RUNS) {
        let read: Vec<Result<RunKeys, String>> = group
            .par_iter()
            .map(|&r| {
                let mut run = RunKeys::default();
                for batch in tables.open(r, &columns)? {
                    for row in batch? {
                        run.rows += 1;
                        if merged {
                            run.spec.push(spec_keys.key(&row.spec_id));
                        }
                        let target = !row.is_decoy;
                        match run.keys.get_mut(row.precursor_key.as_str()) {
                            Some(m) => {
                                if target && !m.has_target {
                                    *m = KeyMeta::of(&row, r);
                                }
                            }
                            None => {
                                let meta = KeyMeta::of(&row, r);
                                run.keys.insert(row.precursor_key, meta);
                            }
                        }
                    }
                }
                Ok(run)
            })
            .collect();
        for (&r, run) in group.iter().zip(read) {
            merge.add(r, run?);
        }
    }
    Ok(merge)
}

/// What joining one run's PSMs to their rescoring scores reads.
struct JoinInputs<'a, R> {
    tables: &'a RunTables<'a>,
    merged_scores: &'a Option<ScoreSet>,
    per_run: &'a Option<RunScores<'a>>,
    id_of: &'a HashMap<&'a str, (u32, u32)>,
    target_key: &'a [bool],
    resolve: R,
    donor_q: f64,
    donor_pep: Option<f64>,
}

impl<R: Fn(&str) -> Option<f64> + Sync> JoinInputs<'_, R> {
    /// Run `r`'s PSMs joined on SpecId to their scores, with the counts of those skipped.
    fn join(&self, r: usize) -> Result<RunJoin, String> {
        let own;
        let set = match (self.merged_scores, self.per_run) {
            (Some(m), _) => m,
            (None, Some(load)) => {
                own = load(&self.tables.stems[r])?;
                &own
            }
            (None, None) => unreachable!("one score source"),
        };
        let mut class = IsomerClass::default();
        let mut run = RunJoin::default();
        for batch in self.tables.open(r, &JOIN_COLUMNS)? {
            for row in batch? {
                run.rows += 1;
                let (id, site) = *self.id_of.get(row.precursor_key.as_str()).ok_or_else(|| {
                    format!(
                        "{}: precursor {:?} changed between two reads of the table",
                        self.tables.paths[r].display(),
                        row.precursor_key
                    )
                })?;
                if row.is_decoy && self.target_key[id as usize] {
                    run.shadowed += 1;
                    continue;
                }
                let Some(&j) = set.get(&row.spec_id) else {
                    if row.is_decoy {
                        run.unscored_d += 1;
                    } else {
                        run.unscored_t += 1;
                    }
                    continue;
                };
                run.joined += 1;
                class
                    .parse_into(&row.peptidoform, &self.resolve)
                    .map_err(|e| format!("{}: {e}", self.tables.paths[r].display()))?;
                if !set.agrees(&j, row.is_decoy, &class) {
                    run.mismatched += 1;
                    continue;
                }
                let p = PsmIn {
                    q_gate: j.pooled_q,
                    run: r as u16,
                    precursor: id,
                    rt: row.rt_s,
                    score: j.score,
                    q: j.q,
                    pep: j.pep,
                    is_decoy: row.is_decoy,
                    faims_cv: row.faims_cv,
                    ion_mobility: row.ion_mobility,
                };
                let donor = is_donor_psm(&p, self.donor_q, self.donor_pep);
                run.psms
                    .push((p, site, donor.then(|| row.spec_id.as_str().into())));
            }
        }
        if self.per_run.is_some() && run.joined == 0 && run.rows > 0 && !set.psms.is_empty() {
            return Err(format!(
                "run {}: none of the {} rescored PSMs in {} has a SpecId of this run; are these the rescoring results of this run?",
                self.tables.stems[r],
                set.psms.len(),
                set.source
            ));
        }
        if run.mismatched as f64 > MAX_JOIN_MISMATCH * run.joined as f64 {
            return Err(format!(
                "run {}: {} of {} PSMs joined on SpecId disagree with {} on the peptide, its modifications or the target/decoy label; are these the rescoring results of this run?",
                self.tables.stems[r], run.mismatched, run.joined, set.source
            ));
        }
        Ok(run)
    }
}

/// Each site key's precursor id and site index, for joining PSMs to the library.
/// The modification masses the runs' quant tables declare, keyed by tag.
fn quant_mod_masses(paths: &[PathBuf]) -> Result<HashMap<String, f64>, String> {
    let mut mod_masses: HashMap<String, f64> = HashMap::new();
    for p in paths {
        mod_masses.extend(read_run_quant_mods(p)?);
    }
    Ok(mod_masses)
}

/// The donor-PSM gate as logged: q-value, and PEP when one is set.
fn donor_gate_text(donor_q: f64, donor_pep: Option<f64>) -> String {
    match donor_pep {
        Some(c) => format!("PSM q<={donor_q} and PEP<={c}"),
        None => format!("PSM q<={donor_q}"),
    }
}

fn site_ids<'k>(
    site_keys: &'k [(String, KeyMeta)],
    forms: &[Vec<u32>],
) -> HashMap<&'k str, (u32, u32)> {
    forms
        .iter()
        .enumerate()
        .flat_map(|(id, sites)| {
            sites
                .iter()
                .map(move |&i| (site_keys[i as usize].0.as_str(), (id as u32, i)))
        })
        .collect()
}

/// The precursor keys of all runs, merged: a key's metadata comes from its first target row, else
/// its first row; with merged scores, the spec ids seen and the first one shared by two runs.
#[derive(Default)]
struct KeyMerge {
    keys: BTreeMap<String, KeyMeta>,
    seen: HashMap<u128, u32, rustc_hash::FxBuildHasher>,
    shared: Option<(usize, usize, usize)>,
    n_rows: usize,
}

impl KeyMerge {
    /// Adds the keys and spec ids read from run `r`.
    fn add(&mut self, r: usize, run: RunKeys) {
        let KeyMerge {
            keys,
            seen,
            shared,
            n_rows,
        } = self;
        *n_rows += run.rows;
        for (i, h) in run.spec.iter().enumerate() {
            let other = *seen.entry(*h).or_insert(r as u32);
            if shared.is_none() && other as usize != r {
                *shared = Some((other as usize, r, i));
            }
        }
        for (key, m) in run.keys {
            match keys.get_mut(key.as_str()) {
                Some(g) => {
                    if m.has_target && !g.has_target {
                        *g = m;
                    }
                }
                None => {
                    keys.insert(key, m);
                }
            }
        }
    }
}

/// Keys each precursor by the site form of its best PSM; returns the precursors that join several
/// site forms and the site forms of the target precursors.
fn rekey_by_best_form(
    meta: &mut [PrecursorMeta],
    best: &[Option<BestForm>],
    site_keys: &[(String, KeyMeta)],
    target_key: &[bool],
    forms: &[Vec<u32>],
) -> (usize, usize) {
    let (mut multi, mut n_forms) = (0usize, 0usize);
    for (id, b) in best.iter().enumerate() {
        if let Some(b) = b {
            meta[id].key = site_keys[b.site as usize].0.clone();
        }
        if target_key[id] {
            n_forms += forms[id].len();
            multi += usize::from(forms[id].len() > 1);
        }
    }
    (multi, n_forms)
}

/// The library's precursor table: one precursor per isomer class and charge, its metadata,
/// proteins and whether a target form exists, and the site keys (forms) it groups.
struct IsomerTable {
    meta: Vec<PrecursorMeta>,
    proteins: Vec<String>,
    target_key: Vec<bool>,
    forms: Vec<Vec<u32>>,
}

/// Groups the site keys into isomer classes (per charge), resolving modification tags with the runs'
/// quant modifications and the bundled Unimod table, and checks each class's envelope and isotope share.
fn isomer_classes(
    paths: &[PathBuf],
    site_keys: &[(String, KeyMeta)],
    resolve: impl Fn(&str) -> Option<f64> + Copy,
) -> Result<IsomerTable, String> {
    let mut classes: BTreeMap<(IsomerClass, u8), Vec<u32>> = BTreeMap::new();
    for (i, (key, m)) in site_keys.iter().enumerate() {
        let class = IsomerClass::parse(split_key(key), resolve)
            .map_err(|e| format!("{}: precursor {key:?}: {e}", paths[m.run].display()))?;
        classes.entry((class, m.charge)).or_default().push(i as u32);
    }
    let mut meta = Vec::with_capacity(classes.len());
    let mut proteins = Vec::with_capacity(classes.len());
    let mut target_key = Vec::with_capacity(classes.len());
    let mut forms: Vec<Vec<u32>> = Vec::with_capacity(classes.len());
    for (_, sites) in std::mem::take(&mut classes) {
        let first = *sites
            .iter()
            .find(|&&i| site_keys[i as usize].1.has_target)
            .unwrap_or(&sites[0]);
        let (key, m) = &site_keys[first as usize];
        let envelope = m.envelope.map_err(|len| {
            format!(
                "{}: precursor {key:?} has an envelope of {len} isotopes, expected {ENVELOPE_LEN}",
                paths[m.run].display()
            )
        })?;
        if !(m.isotope_share > 0.0 && m.isotope_share <= 1.0) {
            return Err(format!(
                "{}: precursor {key:?} has an isotope share of {}, outside (0, 1]",
                paths[m.run].display(),
                m.isotope_share
            ));
        }
        meta.push(PrecursorMeta {
            key: key.clone(),
            mz: m.mz,
            charge: m.charge,
            envelope,
            isotope_share: m.isotope_share,
        });
        proteins.push(m.proteins.clone());
        target_key.push(m.has_target);
        forms.push(sites);
    }
    Ok(IsomerTable {
        meta,
        proteins,
        target_key,
        forms,
    })
}
