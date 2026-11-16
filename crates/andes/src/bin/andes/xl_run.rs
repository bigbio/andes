//! The `--xl` standalone driver: crosslink search, the pair-aware
//! `.xl.pin`, optional Percolator rescoring (`--rescore`), and the CSM,
//! residue-pair, protein-pair and doublet tables.

use std::io::Write;
use std::path::{Path, PathBuf};

use crate::cli::Cli;
use andes_xl::crosslinker::Crosslinker;
use andes_xl::decoy::PairClass;
use model::Spectrum;
use search::xl_search::{xl_search_run, PeptideHit, XlCandidates, XlConfig, XlSpectrumResult};
use search::PreparedSearch;

/// `<base>.<suffix>` next to the PIN path (a trailing `.pin` is replaced).
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
    cols.extend(linker.forms().map(|f| fmt(h.form_scores[f.index()])));
    cols.extend([
        h.n_sites.to_string(),
        h.gap_next.map_or("NA".to_string(), |g| format!("{g:.3}")),
        format!("{:.2}", h.mass_err_ppm),
        p.peptide.length().to_string(),
    ]);
    cols
}

/// Per-peptide columns: `score_<form>` for each form of the linker (its stubs,
/// then `intact`) sit between `form` and `n_sites`.
fn pep_cols(linker: &Crosslinker) -> Vec<String> {
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
    v.extend(
        linker
            .forms()
            .map(|f| format!("score_{}", linker.form_name(f))),
    );
    v.extend(
        ["n_sites", "gap_next", "mass_err_ppm", "length"]
            .iter()
            .map(|s| s.to_string()),
    );
    v
}

fn scan_of(spec: &Spectrum) -> String {
    spec.scan
        .map_or_else(|| spec.title.clone(), |s| s.to_string())
}

fn write_csms(
    path: &Path,
    spectra: &[Spectrum],
    results: &[XlSpectrumResult],
    score: &[Option<f64>],
    q: &[Option<f64>],
    table: &XlCandidates,
    prepared: &PreparedSearch,
    linker: &Crosslinker,
) -> std::io::Result<usize> {
    let mut w = std::io::BufWriter::new(std::fs::File::create(path)?);
    let mut head: Vec<String> = [
        "spectrum",
        "scan",
        "charge",
        "precursor_mass",
        "class",
        "link",
        "score",
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
        head.extend(pep_cols(linker).iter().map(|c| format!("{side}_{c}")));
    }
    writeln!(w, "{}", head.join("\t"))?;
    let mut n = 0;
    for ((r, s), q) in results.iter().zip(score).zip(q) {
        let Some(c) = &r.csm else { continue };
        let spec = &spectra[r.spec_idx];
        let mut row = vec![
            spec.title.clone(),
            scan_of(spec),
            r.charge.to_string(),
            format!("{:.5}", r.precursor_mass),
            c.class.name().to_string(),
            if c.intra { "intra" } else { "inter" }.to_string(),
            s.map_or("NA".to_string(), |v| format!("{v:.4}")),
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
        row.extend(peptide_cols(&c.alpha, table, prepared, linker));
        row.extend(peptide_cols(&c.beta, table, prepared, linker));
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
    table: &XlCandidates,
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

fn report(results: &[XlSpectrumResult]) {
    let n = results.len();
    let mut hist = [0usize; 4];
    for r in results {
        hist[r.doublets.masses.len().min(3)] += 1;
    }
    eprintln!(
        "[xl] doublets: {n} spectra (charge >= gate); doublet peptide masses per spectrum \
         0: {}, 1: {}, 2: {}, >=3: {}",
        hist[0], hist[1], hist[2], hist[3]
    );
    eprintln!(
        "[xl] mass pairs: {} spectra with a full pair, {} with a half pair only",
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

/// The linkers the search runs: MS-cleavable ones (candidates from reporter
/// doublets) with one site set for both ends.
fn check_searchable(linker: &Crosslinker) -> Result<(), String> {
    if linker.cleavage.is_none() {
        return Err(format!(
            "--xl {}: a non-cleavable linker; the crosslink search finds candidates from \
             the reporter doublets of an MS-cleavable linker (DSSO, DSBU, DSBSO, CDI or a \
             custom linker with stubs)",
            linker.name
        ));
    }
    if !linker.is_homobifunctional() {
        return Err(format!(
            "--xl {} ({}): a heterobifunctional linker; the crosslink search takes one site \
             set for both ends",
            linker.name,
            linker.describe_sites()
        ));
    }
    Ok(())
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
    eprintln!(
        "[xl] linker {}{}: intact {:+.6} Da, sites {}{}",
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

/// Crosslink mode: search every spectrum for crosslinks of `--xl`'s linker and
/// write `<base>.xl.pin`, `<base>.xl.tsv` (CSMs), `<base>.xl.residue_pairs.tsv`,
/// `<base>.xl.ppi.tsv` and `<base>.xl.doublets.tsv`. CSMs are scored by the
/// pair score; with `--rescore` the intra-links are scored by Percolator on
/// their own PIN (`<base>.xl.intra.pin`). The FDR is pooled over every input
/// run.
pub(crate) fn run_xl(
    cli: &Cli,
    spectra: &[Spectrum],
    prepared: &PreparedSearch,
    report_base: &Path,
    spectrum_paths: &[PathBuf],
    run_ends: &[usize],
    t_total: std::time::Instant,
) -> Result<(), Box<dyn std::error::Error>> {
    let linker: Crosslinker = cli.xl.clone().ok_or("--xl: no crosslinker")?;
    check_searchable(&linker)?;
    report_linker(&linker);
    let t = std::time::Instant::now();
    let table = XlCandidates::from_out_of_core(prepared, &linker);
    let n_decoy = table.peptides().iter().filter(|p| p.is_decoy).count();
    eprintln!(
        "[xl] {} linkable peptidoforms ({} decoy) [{:.2}s]",
        table.peptides().len(),
        n_decoy,
        t.elapsed().as_secs_f64()
    );
    let spectra_for_xl: &[Spectrum] = if cli.max_spectra > 0 {
        &spectra[..spectra.len().min(cli.max_spectra)]
    } else {
        spectra
    };
    let mut cfg = XlConfig::new(linker.clone());
    cfg.oracle_spectra = cli.xl_oracle;
    let t = std::time::Instant::now();
    let results = xl_search_run(spectra_for_xl, prepared, &table, &cfg);
    eprintln!(
        "[xl] searched {} spectra, {} CSMs [{:.2}s]",
        results.len(),
        results.iter().filter(|r| r.csm.is_some()).count(),
        t.elapsed().as_secs_f64()
    );
    report(&results);
    // The input run of a spectrum: the first run whose end lies past it.
    let run_of = |i: usize| -> String {
        let k = run_ends.partition_point(|&end| end <= i);
        spectrum_paths
            .get(k)
            .and_then(|p| p.file_stem())
            .map_or_else(|| "run".to_string(), |s| s.to_string_lossy().into_owned())
    };
    let pin_path = sibling(report_base, ".xl.pin");
    let n_pin = crate::xl_fdr::write_xl_pin(
        &pin_path,
        spectra_for_xl,
        &results,
        &table,
        prepared,
        &run_of,
        None,
        cfg.linker.intact,
    )?;
    eprintln!("[xl] wrote {} ({n_pin} CSM rows)", pin_path.display());
    // Inter-links keep the pair score: a Percolator model over all CSMs or
    // over the inter rows does not hold the inter FDR on the PXD042173 truth.
    // `--rescore` rescores the intra-links with a Percolator model trained on
    // the intra rows only; `--xl-rescore-inter` (a comparison) does the same
    // for the inter rows.
    let mut score: Vec<Option<f64>> = results
        .iter()
        .map(|r| r.csm.as_ref().map(|c| c.pair_score as f64))
        .collect();
    let mut rescore_class = |intra: bool, suffix: &str| -> Result<(), Box<dyn std::error::Error>> {
        let path = sibling(report_base, suffix);
        let n = crate::xl_fdr::write_xl_pin(
            &path,
            spectra_for_xl,
            &results,
            &table,
            prepared,
            &run_of,
            Some(intra),
            cfg.linker.intact,
        )?;
        eprintln!("[xl] wrote {} ({n} CSM rows)", path.display());
        let s = crate::xl_fdr::percolator_scores(cli, &path, &results)?;
        for ((dst, src), r) in score.iter_mut().zip(s).zip(&results) {
            if r.csm.as_ref().is_some_and(|c| c.intra == intra) {
                *dst = src;
            }
        }
        Ok(())
    };
    if cli.rescore {
        rescore_class(true, ".xl.intra.pin")?;
    }
    if cli.xl_rescore_inter {
        rescore_class(false, ".xl.inter.pin")?;
    }
    let levels = crate::xl_fdr::roll_up(&results, &score, &table, prepared);
    crate::xl_fdr::report_levels(&levels);
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
        &q,
        &table,
        prepared,
        &linker,
    )?;
    crate::xl_fdr::write_levels(
        &sibling(report_base, ".xl.residue_pairs.tsv"),
        &sibling(report_base, ".xl.ppi.tsv"),
        &levels,
        &results,
        spectra_for_xl,
    )?;
    let dbl_path = sibling(report_base, ".xl.doublets.tsv");
    write_doublets(&dbl_path, spectra_for_xl, &results, &table)?;
    eprintln!(
        "Wrote {} ({n} CSMs), the residue-pair, protein-pair and doublet tables \
         [PHASE TOTAL: {:.2}s]",
        csm_path.display(),
        t_total.elapsed().as_secs_f64()
    );
    Ok(())
}
