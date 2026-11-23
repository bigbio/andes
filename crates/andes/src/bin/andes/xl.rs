pub(crate) mod ambiguity;
pub(crate) mod fdr;

use std::io::Write;
use std::path::{Path, PathBuf};

use crate::cli::Cli;
use crate::xl::fdr::XlRun;
use model::Spectrum;
use search::crosslinks::{xl_search_run, PeptideHit, XlCandidates, XlConfig, XlSpectrumResult};
use search::PreparedSearch;
use xl::crosslinker::Crosslinker;
use xl::crosslinker::LinkForm;
use xl::decoy::PairClass;

fn sibling(report_base: &Path, suffix: &str) -> PathBuf {
    let mut s = if report_base
        .extension()
        .is_some_and(|e| e.eq_ignore_ascii_case("pin"))
    {
        report_base.with_extension("").into_os_string()
    } else {
        report_base.to_path_buf().into_os_string()
    };
    s.push(suffix);
    PathBuf::from(s)
}

fn peptide_cols(
    h: &PeptideHit,
    table: &XlCandidates,
    prepared: &PreparedSearch,
    linker: &Crosslinker,
    stub_cols: Option<usize>,
) -> Vec<String> {
    let p = &table.peptides()[h.pep as usize];
    let proteins: Vec<&str> = p
        .proteins
        .iter()
        .filter_map(|&i| prepared.idx.db.proteins.get(i as usize))
        .map(|pr| pr.accession.as_str())
        .collect();
    let fmt = |v: f32| {
        if v.is_finite() {
            format!("{v:.3}")
        } else {
            "NA".to_string()
        }
    };
    let site_kind = if h.site.n_term {
        "nterm".to_string()
    } else if h.site.c_term {
        "cterm".to_string()
    } else {
        p.peptide
            .residues
            .get(h.site.pos)
            .map_or("NA".to_string(), |aa| (aa.residue as char).to_string())
    };
    let mut cols = vec![
        p.peptide.to_proforma(),
        (h.site.pos + 1).to_string(),
        site_kind,
        proteins.join(";"),
        u8::from(p.is_decoy).to_string(),
        fmt(h.score),
        linker.form_name(h.form).to_string(),
    ];
    match stub_cols {
        None => cols.extend(linker.forms().map(|f| fmt(h.form_scores[f.index()]))),
        Some(n) => {
            let stubs: Vec<String> = linker
                .forms()
                .filter(|f| *f != LinkForm::Intact)
                .map(|f| fmt(h.form_scores[f.index()]))
                .collect();
            cols.extend((0..n).map(|k| stubs.get(k).cloned().unwrap_or_else(|| "NA".into())));
            cols.push(fmt(h.form_scores[LinkForm::Intact.index()]));
        }
    }
    cols.extend([
        h.n_sites.to_string(),
        h.gap_next.map_or("NA".to_string(), |g| format!("{g:.3}")),
        format!("{:.2}", h.mass_err_ppm),
        p.peptide.length().to_string(),
        h.ms3.map_or("NA".to_string(), |v| format!("{v:.3}")),
    ]);
    cols
}

fn pep_cols(linker: &Crosslinker, stub_cols: Option<usize>) -> Vec<String> {
    let mut v: Vec<String> = [
        "peptide",
        "site",
        "site_kind",
        "proteins",
        "decoy",
        "score",
        "form",
    ]
    .iter()
    .map(|s| s.to_string())
    .collect();
    match stub_cols {
        None => v.extend(
            linker
                .forms()
                .map(|f| format!("score_{}", linker.form_name(f))),
        ),
        Some(n) => {
            v.extend((1..=n).map(|k| format!("score_stub{k}")));
            v.push("score_intact".into());
        }
    }
    v.extend(
        ["n_sites", "gap_next", "mass_err_ppm", "length", "ms3"]
            .iter()
            .map(|s| s.to_string()),
    );
    v
}

/// The scan number in a native spectrum id (`... scan=<n>`).
fn native_scan(id: &str) -> Option<i32> {
    let rest = &id[id.rfind("scan=")? + "scan=".len()..];
    let digits: String = rest.chars().take_while(char::is_ascii_digit).collect();
    digits.parse().ok()
}

/// Files `ms3` under the MS2 scan each was taken from: `ms2` are one run's spectra, at
/// `offset` in the search's spectrum list, their titles the file's native ids after
/// `title_prefix`, and `out[offset + i]` receives the MS3 scans of `ms2[i]`. An MS3 scan is
/// linked by its parent's native id (`scan=`, `index=` or any other form the file uses),
/// else by the scan number in it. Returns how many were linked.
pub(crate) fn link_ms3(
    ms2: &[Spectrum],
    offset: usize,
    title_prefix: Option<&str>,
    ms3: Vec<model::scan::ProductScan>,
    out: &mut Vec<Vec<model::scan::ProductScan>>,
) -> usize {
    out.resize_with(offset + ms2.len(), Vec::new);
    let id_of = |s: &Spectrum| -> String {
        let t = s.title.as_str();
        title_prefix
            .and_then(|p| t.strip_prefix(p))
            .unwrap_or(t)
            .to_string()
    };
    let by_id: std::collections::HashMap<String, usize> =
        ms2.iter().enumerate().map(|(i, s)| (id_of(s), i)).collect();
    let by_scan: std::collections::HashMap<i32, usize> = ms2
        .iter()
        .enumerate()
        .filter_map(|(i, s)| s.scan.map(|n| (n, i)))
        .collect();
    let mut linked = 0;
    for p in ms3 {
        let Some(parent) = p.parent_id.as_deref() else {
            continue;
        };
        let hit = by_id
            .get(parent)
            .or_else(|| native_scan(parent).and_then(|n| by_scan.get(&n)));
        if let Some(&i) = hit {
            out[offset + i].push(p);
            linked += 1;
        }
    }
    linked
}

fn scan_of(spec: &Spectrum) -> String {
    spec.scan
        .map_or_else(|| spec.title.clone(), |s| s.to_string())
}

#[allow(clippy::too_many_arguments)]
/// Which link classes had their CSM scores replaced by this run's own Percolator scores
/// (`--xl-rescore-intra`, `--xl-rescore-inter`).
#[derive(Clone, Copy)]
struct Rescored {
    intra: bool,
    inter: bool,
}

impl Rescored {
    /// The `score_source` column: `pair` for the search's pair score, which is on one scale
    /// across runs, or `percolator` for a score learned on this run alone.
    fn source(self, intra: bool) -> &'static str {
        if (intra && self.intra) || (!intra && self.inter) {
            "percolator"
        } else {
            "pair"
        }
    }
}

#[allow(clippy::too_many_arguments)]
fn write_csms(
    path: &Path,
    spectra: &[Spectrum],
    results: &[XlSpectrumResult],
    score: &[Option<f64>],
    rescored: Rescored,
    q: &[Option<f64>],
    runs: &[XlRun],
    prepared: &PreparedSearch,
) -> std::io::Result<usize> {
    let linker = &runs[0].linker;
    let stub_cols = (runs.len() > 1).then(|| {
        runs.iter()
            .map(|r| r.linker.cleavage.as_ref().map_or(0, |c| c.stubs.len()))
            .max()
            .unwrap_or(0)
    });
    let mut w = std::io::BufWriter::new(std::fs::File::create(path)?);
    let mut head: Vec<String> = [
        "spectrum",
        "scan",
        "spec_id",
        "charge",
        "precursor_mass",
        "class",
        "link",
        "score",
        "score_source",
        "q_value",
        "pair_score",
        "second_pair_score",
        "pair_kind",
        "isotope",
        "precursor_err_ppm",
        "complement_stub_ions",
        "n_mass_pairs",
        "n_doublet_masses",
        "doublet_rel_intensity",
    ]
    .iter()
    .map(|s| s.to_string())
    .collect();
    for side in ["alpha", "beta"] {
        head.extend(
            pep_cols(linker, stub_cols)
                .iter()
                .map(|c| format!("{side}_{c}")),
        );
    }
    if stub_cols.is_some() {
        head.push("crosslinker".into());
    }
    head.extend(["residue_pair".to_string(), "protein_pair".to_string()]);
    writeln!(w, "{}", head.join("\t"))?;
    let mut n = 0;
    for ((r, s), q) in results.iter().zip(score).zip(q) {
        let Some(c) = &r.csm else { continue };
        let spec = &spectra[r.spec_idx];
        let mut row = vec![
            spec.title.clone(),
            scan_of(spec),
            crate::xl::fdr::spec_id(r),
            r.charge.to_string(),
            format!("{:.5}", r.precursor_mass),
            c.class.name().to_string(),
            xl::rollup::LinkClass::of(c.intra, c.homotypic)
                .name()
                .to_string(),
            s.map_or("NA".to_string(), |v| v.to_string()),
            rescored.source(c.intra).to_string(),
            q.map_or("NA".to_string(), |v| format!("{v:.5}")),
            format!("{:.3}", c.pair_score),
            c.second_pair_score
                .map_or("NA".to_string(), |v| format!("{v:.3}")),
            c.kind.name().to_string(),
            c.iso.to_string(),
            format!("{:.2}", c.precursor_err_ppm),
            c.b_stub_ions.to_string(),
            c.n_mass_pairs.to_string(),
            r.doublets.masses.len().to_string(),
            format!("{:.4}", r.doublets.top_rel_intensity),
        ];
        let run = &runs[r.linker];
        row.extend(peptide_cols(
            &c.alpha,
            &run.table,
            prepared,
            &run.linker,
            stub_cols,
        ));
        row.extend(peptide_cols(
            &c.beta,
            &run.table,
            prepared,
            &run.linker,
            stub_cols,
        ));
        if stub_cols.is_some() {
            row.push(run.linker.name.clone());
        }
        let (rp_key, ppi_key) = crate::xl::fdr::pair_keys(c, &run.table, prepared);
        row.extend([rp_key, ppi_key]);
        writeln!(w, "{}", row.join("\t"))?;
        n += 1;
    }
    w.flush()?;
    Ok(n)
}

fn write_doublets(
    path: &Path,
    spectra: &[Spectrum],
    results: &[XlSpectrumResult],
    runs: &[XlRun],
) -> std::io::Result<()> {
    let mut w = std::io::BufWriter::new(std::fs::File::create(path)?);
    writeln!(
        w,
        "spectrum\tscan\tcharge\tprecursor_mass\tn_peaks\tbase_peak\ttic\tn_doublets\t\
         n_doublets_iso\tn_doublet_masses\ttop_rel_intensity\tdoublet_masses\t\
         doublet_intensity_tic\tn_full_pairs\tn_half_pairs\toracle_pairs\toracle_score\t\
         oracle_class\toracle_alpha\toracle_beta\toracle_doublet_support\tcsm_score\t\
         oracle_csm_pair\tprobe"
    )?;
    for r in results {
        let spec = &spectra[r.spec_idx];
        let d = &r.doublets;
        let masses: Vec<String> = d
            .masses
            .iter()
            .map(|m| {
                format!(
                    "{:.4}:{}{}",
                    m.mass,
                    m.min_charge,
                    if m.isotope_confirmed { "i" } else { "" }
                )
            })
            .collect();
        let inten: Vec<String> = d
            .masses
            .iter()
            .map(|m| {
                format!(
                    "{:.4}",
                    if d.tic > 0.0 {
                        m.intensity / d.tic
                    } else {
                        0.0
                    }
                )
            })
            .collect();
        let (op, os, oc, oa, ob, osup) = match &r.oracle {
            Some(o) => {
                let table = &runs[r.linker].table;
                let seq = |i: u32| table.peptides()[i as usize].peptide.to_proforma();
                match o.best {
                    Some((i, j, s)) => (
                        o.n_pairs.to_string(),
                        format!("{s:.3}"),
                        o.class.map_or("NA", PairClass::name).to_string(),
                        seq(i),
                        seq(j),
                        format!(
                            "{}{}",
                            u8::from(o.doublet_support.0),
                            u8::from(o.doublet_support.1)
                        ),
                    ),
                    None => (
                        o.n_pairs.to_string(),
                        "NA".into(),
                        "NA".into(),
                        "NA".into(),
                        "NA".into(),
                        "NA".into(),
                    ),
                }
            }
            None => (
                "NA".into(),
                "NA".into(),
                "NA".into(),
                "NA".into(),
                "NA".into(),
                "NA".into(),
            ),
        };
        writeln!(
            w,
            "{}\t{}\t{}\t{:.5}\t{}\t{:.1}\t{:.1}\t{}\t{}\t{}\t{:.4}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}",
            spec.title,
            scan_of(spec),
            r.charge,
            r.precursor_mass,
            d.n_peaks,
            d.base_peak,
            d.tic,
            d.n_doublets,
            d.n_doublets_iso,
            d.masses.len(),
            d.top_rel_intensity,
            masses.join(";"),
            inten.join(";"),
            d.n_full_pairs,
            d.n_half_pairs,
            op,
            os,
            oc,
            oa,
            ob,
            osup,
            r.csm
                .as_ref()
                .map_or("NA".to_string(), |c| format!("{:.3}", c.pair_score)),
            r.oracle
                .as_ref()
                .and_then(|o| o.csm_pair)
                .map_or("NA".to_string(), |(found, s)| format!("{}:{s:.3}", u8::from(found))),
            d.probe
                .iter()
                .flat_map(|f| f.iter().flat_map(|i| i.iter()))
                .map(|v| format!("{v:.1}"))
                .collect::<Vec<_>>()
                .join(","),
        )?;
    }
    w.flush()
}

/// How many spectra had doublets explained by chance, and how many doublet masses were left
/// out of the pairs for it (`search::crosslinks::paired_masses`).
fn report_chance_doublets(results: &[XlSpectrumResult]) {
    let (mut spectra, mut dropped, mut total) = (0usize, 0usize, 0usize);
    for r in results {
        let all = r.doublets.masses.len();
        total += all;
        if r.doublets.n_paired_masses < all {
            spectra += 1;
            dropped += all - r.doublets.n_paired_masses;
        }
    }
    if spectra > 0 {
        log::info!(
            target: "xl",
            "{spectra} spectra have doublets explained by chance; {dropped} of {total} doublet \
             masses were not paired (the least intense)"
        );
    }
}

/// How many searched spectra had MS3 scans, and how many CSMs an MS3 scan supports on one
/// or both peptides.
fn report_ms3(results: &[XlSpectrumResult], ms3: &[Vec<model::scan::ProductScan>]) {
    let with = results
        .iter()
        .filter(|r| ms3.get(r.spec_idx).is_some_and(|v| !v.is_empty()))
        .count();
    if with == 0 {
        return;
    }
    let (mut one, mut both, mut csms) = (0usize, 0usize, 0usize);
    for c in results.iter().filter_map(|r| r.csm.as_ref()) {
        csms += 1;
        match (c.alpha.ms3.is_some(), c.beta.ms3.is_some()) {
            (true, true) => both += 1,
            (true, false) | (false, true) => one += 1,
            _ => {}
        }
    }
    log::info!(
        target: "xl",
        "{with} searched spectra have MS3 scans; of {csms} CSMs, an MS3 scan matches both \
         peptides in {both} and one peptide in {one}"
    );
}

fn report(results: &[XlSpectrumResult]) {
    let n = results.len();
    let mut hist = [0usize; 4];
    for r in results {
        hist[r.doublets.masses.len().min(3)] += 1;
    }
    log::info!(
        target: "xl",
        "doublets: {n} spectra (charge >= gate); doublet peptide masses per spectrum \
         0: {}, 1: {}, 2: {}, >=3: {}",
        hist[0], hist[1], hist[2], hist[3]
    );
    log::info!(
        target: "xl",
        "mass pairs: {} spectra with a full pair, {} with a half pair only",
        results
            .iter()
            .filter(|r| r.doublets.n_full_pairs > 0)
            .count(),
        results
            .iter()
            .filter(|r| r.doublets.n_full_pairs == 0 && r.doublets.n_half_pairs > 0)
            .count()
    );
}

fn report_linker(linker: &Crosslinker) {
    let stubs = linker.cleavage.as_ref().map_or(String::new(), |c| {
        let s: Vec<String> = c
            .stubs
            .iter()
            .map(|s| format!("{} {:+.6}", s.name, s.mass))
            .collect();
        let d: Vec<String> = (0..c.doublets.len())
            .map(|k| format!("{:.6}", c.doublet_delta(k)))
            .collect();
        format!(
            "; stubs {}; doublet spacing {} Da",
            s.join(", "),
            d.join(", ")
        )
    });
    log::info!(
        target: "xl",
        "linker {}{}: intact {:+.6} Da, sites {}{}",
        linker.name,
        linker
            .accession
            .as_deref()
            .map_or(String::new(), |a| format!(" ({a})")),
        linker.intact,
        linker.describe_sites(),
        stubs
    );
}

pub(crate) fn water_monolink_mods(linker: &Crosslinker) -> Vec<model::Modification> {
    use model::modification::{ModLocation, ResidueSpec};
    let Some(water) = linker
        .monolinks
        .iter()
        .find(|m| m.name.ends_with("water") || water_record(&m.name))
    else {
        return Vec::new();
    };
    let end = &linker.ends[0];
    let mut sites: Vec<(ResidueSpec, ModLocation)> = end
        .residues
        .iter()
        .map(|&r| (ResidueSpec::Specific(r), ModLocation::Anywhere))
        .collect();
    if end.protein_n_term {
        sites.push((ResidueSpec::Wildcard, ModLocation::ProtNTerm));
    }
    if end.peptide_n_term {
        sites.push((ResidueSpec::Wildcard, ModLocation::NTerm));
    }
    sites
        .into_iter()
        .map(|(residue, location)| model::Modification {
            accession: water.accession.clone(),
            ..model::Modification::new(water.name.clone(), water.mass, residue, location, false)
        })
        .collect()
}

fn water_record(name: &str) -> bool {
    model::xl_reagents::bundled().reagents.iter().any(|r| {
        r.monolinks
            .iter()
            .any(|m| m.name == "water" && m.unimod.as_deref() == Some(name))
    })
}

pub(crate) fn xl_default_variable_mods(
    explicit: Option<&crate::cli::ModList>,
    linkers: &[Crosslinker],
) -> Option<crate::cli::ModList> {
    if explicit.is_some() {
        return None;
    }
    let mono: Vec<model::Modification> = linkers.iter().flat_map(water_monolink_mods).collect();
    if mono.is_empty() {
        return None;
    }
    let mut mods = crate::models::default_variable_mods();
    mods.extend(mono);
    Some(crate::cli::ModList(mods))
}

#[allow(clippy::too_many_arguments)]
pub(crate) fn run_xl(
    cli: &Cli,
    spectra: &[Spectrum],
    ms3: &[Vec<model::scan::ProductScan>],
    prepared: &PreparedSearch,
    report_base: &Path,
    spectrum_paths: &[PathBuf],
    run_ends: &[usize],
    t_total: std::time::Instant,
) -> Result<(), Box<dyn std::error::Error>> {
    let linkers = cli.xl.clone().ok_or("--xl: no crosslinker")?.0;
    let several = linkers.len() > 1;
    let spectra_for_xl: &[Spectrum] = if cli.max_spectra > 0 {
        &spectra[..spectra.len().min(cli.max_spectra)]
    } else {
        spectra
    };
    let mut runs: Vec<XlRun> = Vec::new();
    let mut results: Vec<XlSpectrumResult> = Vec::new();
    for (li, linker) in linkers.into_iter().enumerate() {
        report_linker(&linker);
        let t = std::time::Instant::now();
        let table = XlCandidates::from_index(prepared, &linker);
        let n_decoy = table.peptides().iter().filter(|p| p.is_decoy).count();
        log::info!(
            target: "xl",
            "{} linkable peptidoforms ({} decoy){} [{:.2}s]",
            table.peptides().len(),
            n_decoy,
            if table.has_open_index() {
                "; open-mass fragment index built (non-cleavable linker)"
            } else {
                ""
            },
            t.elapsed().as_secs_f64()
        );
        let mut cfg = XlConfig::new(linker.clone());
        cfg.oracle_spectra = cli.xl_oracle;
        cfg.charge_by_pair_score = several;
        if let Some((lo, hi)) = cli.isotope_error {
            cfg.isotope_offsets = lo..=hi;
        }
        let t = std::time::Instant::now();
        let mut found = xl_search_run(spectra_for_xl, ms3, prepared, &table, &cfg);
        log::info!(
            target: "xl",
            "{}: searched {} spectra, {} CSMs [{:.2}s]",
            linker.name,
            found.len(),
            found.iter().filter(|r| r.csm.is_some()).count(),
            t.elapsed().as_secs_f64()
        );
        report(&found);
        report_chance_doublets(&found);
        report_ms3(&found, ms3);
        for r in found.iter_mut() {
            r.linker = li;
        }
        runs.push(XlRun { linker, table });
        if li == 0 {
            results = found;
        } else {
            results = merge_linkers(results, found);
        }
    }
    if runs.len() > 1 {
        settle_runner_ups(&mut results);
        let mut by: Vec<usize> = vec![0; runs.len()];
        for r in results.iter().filter(|r| r.csm.is_some()) {
            by[r.linker] += 1;
        }
        let parts: Vec<String> = runs
            .iter()
            .zip(&by)
            .map(|(r, n)| format!("{} {n}", r.linker.name))
            .collect();
        log::info!(
            target: "xl",
            "best CSM per spectrum across linkers: {}",
            parts.join(", ")
        );
    }
    let run_of = |i: usize| -> String {
        let k = run_ends.partition_point(|&end| end <= i);
        spectrum_paths
            .get(k)
            .and_then(|p| p.file_stem())
            .map_or_else(|| "run".to_string(), |s| s.to_string_lossy().into_owned())
    };
    let pin_path = sibling(report_base, ".xl.pin");
    let n_pin = crate::xl::fdr::write_xl_pin(
        &pin_path,
        spectra_for_xl,
        &results,
        &runs,
        prepared,
        &run_of,
        None,
    )?;
    log::info!(target: "xl", "wrote {} ({n_pin} CSM rows)", pin_path.display());
    let mut score: Vec<Option<f64>> = results
        .iter()
        .map(|r| r.csm.as_ref().map(|c| c.score() as f64))
        .collect();
    let mut rescore_class = |intra: bool, suffix: &str| -> Result<(), Box<dyn std::error::Error>> {
        let path = sibling(report_base, suffix);
        let n = crate::xl::fdr::write_xl_pin(
            &path,
            spectra_for_xl,
            &results,
            &runs,
            prepared,
            &run_of,
            Some(intra),
        )?;
        log::info!(target: "xl", "wrote {} ({n} CSM rows)", path.display());
        let s = crate::xl::fdr::percolator_scores(cli, &path, &results)?;
        for ((dst, src), r) in score.iter_mut().zip(s).zip(&results) {
            if r.csm.as_ref().is_some_and(|c| c.intra == intra) {
                *dst = src;
            }
        }
        Ok(())
    };
    if cli.xl_rescore_intra {
        rescore_class(true, ".xl.intra.pin")?;
    }
    if cli.xl_rescore_inter {
        rescore_class(false, ".xl.inter.pin")?;
    }
    let levels = crate::xl::fdr::roll_up(&results, &score, &runs, prepared);
    crate::xl::fdr::report_levels(&levels);
    let mut q = vec![None; results.len()];
    for (&row, v) in levels.csm_rows.iter().zip(&levels.csm_q) {
        q[row] = *v;
    }
    let csm_path = sibling(report_base, ".xl.tsv");
    let n = write_csms(
        &csm_path,
        spectra_for_xl,
        &results,
        &score,
        Rescored {
            intra: cli.xl_rescore_intra,
            inter: cli.xl_rescore_inter,
        },
        &q,
        &runs,
        prepared,
    )?;
    crate::xl::fdr::write_levels(
        &sibling(report_base, ".xl.residue_pairs.tsv"),
        &sibling(report_base, ".xl.ppi.tsv"),
        &levels,
        &|row| spectra_for_xl[results[row].spec_idx].title.clone(),
    )?;
    let dbl_path = sibling(report_base, ".xl.doublets.tsv");
    write_doublets(&dbl_path, spectra_for_xl, &results, &runs)?;
    log::info!(
        target: "xl",
        "wrote {} ({n} CSMs), the residue-pair, protein-pair and doublet tables; \
         total {:.2} s",
        csm_path.display(),
        t_total.elapsed().as_secs_f64()
    );
    Ok(())
}

/// One result per spectrum: the CSM with the best pair score over the linkers, with the
/// peptide pairs every linker found for it (see `settle_runner_ups`).
fn merge_linkers(
    into: Vec<XlSpectrumResult>,
    other: Vec<XlSpectrumResult>,
) -> Vec<XlSpectrumResult> {
    let mut by_spec: std::collections::HashMap<usize, XlSpectrumResult> =
        other.into_iter().map(|r| (r.spec_idx, r)).collect();
    let mut out: Vec<XlSpectrumResult> = into
        .into_iter()
        .map(|r| match by_spec.remove(&r.spec_idx) {
            Some(o) => merge_one(r, o),
            None => r,
        })
        .collect();
    out.extend(by_spec.into_values());
    out.sort_by_key(|r| r.spec_idx);
    out
}

fn merge_one(r: XlSpectrumResult, o: XlSpectrumResult) -> XlSpectrumResult {
    let (Some(a), Some(b)) = (&r.csm, &o.csm) else {
        return if r.csm.is_some() || o.csm.is_none() {
            r
        } else {
            o
        };
    };
    let (mut win, lose) = if b.pair_score > a.pair_score {
        (o, r)
    } else {
        (r, o)
    };
    win.pairs.extend(lose.pairs);
    win
}

/// After several linkers are merged: each CSM's runner-up becomes the best other peptide pair
/// any linker found for the spectrum at any charge searched. A single-linker search keeps the
/// runner-up of the chosen charge.
fn settle_runner_ups(results: &mut [XlSpectrumResult]) {
    for r in results.iter_mut() {
        if let Some(c) = r.csm.as_mut() {
            c.second_pair_score = search::crosslinks::runner_up(&r.pairs, &c.pair_key);
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn linker(spec: &str) -> Crosslinker {
        crate::cli::parse::parse_xl_one(spec).unwrap()
    }

    #[test]
    fn the_water_monolink_sits_on_the_link_sites_of_the_first_end() {
        let m = water_monolink_mods(&linker("DSSO"));
        assert_eq!(m.len(), 2);
        assert!(m.iter().all(|m| m.name == "Xlink:DSSO[176]" && !m.fixed));
        assert_eq!(m[0].accession.as_deref(), Some("UNIMOD:1878"));
        assert!((m[0].mass_delta - 176.01433).abs() < 1e-9);
        assert_eq!(water_monolink_mods(&linker("DSS (K)")).len(), 1);
        assert_eq!(water_monolink_mods(&linker("DSBU (K, N-term)")).len(), 2);
        let b = water_monolink_mods(&linker("DSBSO"));
        assert_eq!(b[0].name, "DSBSO:water");
        assert!(water_monolink_mods(&linker("CDI")).is_empty());
        assert!(water_monolink_mods(&linker("DTBP")).is_empty());
    }

    fn result(linker: usize, csm: Option<(&str, f32)>, pairs: &[(&str, f32)]) -> XlSpectrumResult {
        use search::crosslinks::{Csm, PeptideHit};
        use xl::crosslinker::LinkForm;
        use xl::sites::LinkSite;
        let key = |p: &str| [p.as_bytes().to_vec(), b"PEPTIDEK".to_vec()];
        let hit = PeptideHit {
            pep: 0,
            site: LinkSite {
                pos: 0,
                n_term: false,
                c_term: false,
            },
            form: LinkForm::Intact,
            score: 0.0,
            form_scores: Default::default(),
            n_sites: 1,
            gap_next: None,
            mass_err_ppm: 0.0,
            ms3: None,
        };
        let pairs: Vec<_> = pairs.iter().map(|(p, s)| (key(p), *s)).collect();
        XlSpectrumResult {
            spec_idx: 7,
            charge: 3,
            precursor_mass: 0.0,
            doublets: Default::default(),
            csm: csm.map(|(p, score)| Csm {
                alpha: hit,
                beta: hit,
                pair_key: key(p),
                pair_score: score,
                second_pair_score: search::crosslinks::runner_up(&pairs, &key(p)),
                kind: search::crosslinks::PairKind::Full,
                iso: 0,
                class: xl::decoy::PairClass::TT,
                intra: true,
                homotypic: false,
                precursor_err_ppm: 0.0,
                b_stub_ions: 0,
                n_mass_pairs: 1,
                alpha_doublet: None,
                beta_doublet: None,
                rescores: Default::default(),
                charge: 3,
            }),
            oracle: None,
            linker,
            pairs,
        }
    }

    fn merged(results: Vec<XlSpectrumResult>) -> (usize, f32, Option<f32>) {
        let mut it = results.into_iter();
        let mut acc = vec![it.next().unwrap()];
        for r in it {
            acc = merge_linkers(acc, vec![r]);
        }
        settle_runner_ups(&mut acc);
        let c = acc[0].csm.as_ref().unwrap();
        (acc[0].linker, c.pair_score, c.second_pair_score)
    }

    #[test]
    fn the_runner_up_covers_every_linker_and_charge() {
        // A different pair of the other linker at 99 is the runner-up of a 100 / 50 CSM.
        let a = result(0, Some(("AK", 100.0)), &[("AK", 100.0), ("BK", 50.0)]);
        let b = result(1, Some(("CK", 99.0)), &[("CK", 99.0), ("DK", 10.0)]);
        assert_eq!(merged(vec![a, b]), (0, 100.0, Some(99.0)));
        // The same peptide pair found by the other linker does not compete with itself.
        let a = result(0, Some(("AK", 100.0)), &[("AK", 100.0), ("BK", 50.0)]);
        let b = result(1, Some(("AK", 99.0)), &[("AK", 99.0), ("DK", 60.0)]);
        assert_eq!(merged(vec![a, b]), (0, 100.0, Some(60.0)));
        // Three linkers: B's pair at another charge (110) survives when C wins later.
        let a = result(0, Some(("AK", 105.0)), &[("AK", 105.0)]);
        let b = result(1, Some(("BK", 100.0)), &[("BK", 100.0), ("EK", 110.0)]);
        let c = result(2, Some(("CK", 120.0)), &[("CK", 120.0), ("DK", 80.0)]);
        assert_eq!(merged(vec![a, b, c]), (2, 120.0, Some(110.0)));
        // The winner's own pairs at other charges count; an alternative that is the
        // winner's own pair does not.
        let a = result(0, Some(("PK", 120.0)), &[("PK", 120.0), ("QK", 110.0)]);
        let b = result(1, Some(("RK", 100.0)), &[("RK", 100.0), ("PK", 115.0)]);
        assert_eq!(merged(vec![a, b]), (0, 120.0, Some(110.0)));
        // A spectrum with no CSM in one linker keeps the other's.
        let a = result(0, None, &[]);
        let b = result(1, Some(("AK", 90.0)), &[("AK", 90.0)]);
        assert_eq!(merged(vec![a, b]), (1, 90.0, None));
        // ... and its runner-up, 50 at its chosen charge, becomes 110 from its other charge.
        let a = result(0, None, &[]);
        let mut b = result(
            1,
            Some(("AK", 100.0)),
            &[("AK", 100.0), ("BK", 50.0), ("CK", 110.0)],
        );
        b.csm.as_mut().unwrap().second_pair_score = Some(50.0);
        assert_eq!(merged(vec![a, b]), (1, 100.0, Some(110.0)));
        // Without a precursor charge, linker A keeps 110 at charge 4 over 100 / 50 at charge
        // 3 (best pair score), and beats linker B's 105; its runner-up is A's own 109.
        let mut a = result(
            0,
            Some(("CK", 110.0)),
            &[("AK", 100.0), ("BK", 50.0), ("CK", 110.0), ("DK", 109.0)],
        );
        a.charge = 4;
        let b = result(1, Some(("EK", 105.0)), &[("EK", 105.0)]);
        let mut acc = merge_linkers(vec![a], vec![b]);
        settle_runner_ups(&mut acc);
        let c = acc[0].csm.as_ref().unwrap();
        assert_eq!(
            (
                acc[0].linker,
                acc[0].charge,
                c.pair_score,
                c.second_pair_score
            ),
            (0, 4, 110.0, Some(109.0))
        );
    }

    #[test]
    fn ms3_scans_are_filed_under_their_parent_ms2_scan() {
        use model::scan::ProductScan;
        let ms2 = |scan: i32| Spectrum {
            scan: Some(scan),
            ..Default::default()
        };
        let child = |id: &str, parent: Option<&str>| ProductScan {
            id: id.to_string(),
            parent_id: parent.map(str::to_string),
            ..Default::default()
        };
        let run = [ms2(10), ms2(14)];
        let mut out = vec![Vec::new(); 3];
        let linked = link_ms3(
            &run,
            3,
            None,
            vec![
                child("a", Some("controllerType=0 controllerNumber=1 scan=14")),
                child("b", Some("controllerType=0 controllerNumber=1 scan=10")),
                child("c", Some("controllerType=0 controllerNumber=1 scan=14")),
                child("d", Some("controllerType=0 controllerNumber=1 scan=99")),
                child("e", None),
            ],
            &mut out,
        );
        assert_eq!(linked, 3);
        assert_eq!(out.len(), 5);
        let ids = |v: &Vec<ProductScan>| v.iter().map(|p| p.id.clone()).collect::<Vec<_>>();
        assert_eq!(ids(&out[3]), ["b"]);
        assert_eq!(ids(&out[4]), ["a", "c"]);
        assert!(out[..3].iter().all(Vec::is_empty));
        assert_eq!(native_scan("scan=7"), Some(7));
        assert_eq!(native_scan("index=7"), None);

        let titled = |title: &str, scan: i32| Spectrum {
            title: title.to_string(),
            scan: Some(scan),
            ..Default::default()
        };
        let run = [titled("R1/index=100", 101), titled("R1/index=101", 102)];
        let mut out = Vec::new();
        let linked = link_ms3(
            &run,
            0,
            Some("R1/"),
            vec![child("m", Some("index=101")), child("n", Some("index=7"))],
            &mut out,
        );
        assert_eq!(linked, 1);
        assert_eq!(ids(&out[1]), ["m"]);
    }

    #[test]
    fn several_linkers_are_a_semicolon_list_and_share_the_monolink_default() {
        let set = crate::cli::parse_xl("DSSO; DSBU (K)").unwrap();
        let names: Vec<&str> = set.0.iter().map(|l| l.name.as_str()).collect();
        assert_eq!(names, ["DSSO", "DSBU"]);
        assert_eq!(set.0[1].ends[0].residues, b"K");
        assert!(crate::cli::parse_xl("DSSO; DSSO (K)")
            .unwrap_err()
            .contains("twice"));
        assert!(crate::cli::parse_xl(" ; ").is_err());
        let d = xl_default_variable_mods(None, &set.0).unwrap().0;
        assert!(d.iter().any(|m| m.name == "Xlink:DSSO[176]"));
        assert!(d.iter().any(|m| m.name == "Xlink:BuUrBu[214]"));
    }

    #[test]
    fn an_explicit_variable_mod_list_replaces_the_monolink_default() {
        let x = linker("DSS");
        let d = xl_default_variable_mods(None, std::slice::from_ref(&x))
            .unwrap()
            .0;
        assert_eq!(
            d.len(),
            crate::models::default_variable_mods().len() + 2
        );
        assert!(d.iter().any(|m| m.name == "Xlink:DSS[156]"));
        let explicit = crate::cli::parse_variable_mods("Oxidation (M)").unwrap();
        assert!(xl_default_variable_mods(Some(&explicit), std::slice::from_ref(&x)).is_none());
        let off = crate::cli::parse_variable_mods("").unwrap();
        assert!(off.0.is_empty());
        assert!(xl_default_variable_mods(Some(&off), std::slice::from_ref(&x)).is_none());
        assert!(xl_default_variable_mods(None, &[linker("CDI")]).is_none());
    }
}
