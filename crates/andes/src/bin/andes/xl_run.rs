use std::io::Write;
use std::path::{Path, PathBuf};

use crate::cli::Cli;
use andes_xl::crosslinker::Crosslinker;
use andes_xl::decoy::PairClass;
use model::Spectrum;
use search::xl_search::{xl_search_run, PeptideHit, XlCandidates, XlConfig, XlSpectrumResult};
use search::PreparedSearch;

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
    cols.extend(linker.forms().map(|f| fmt(h.form_scores[f.index()])));
    cols.extend([
        h.n_sites.to_string(),
        h.gap_next.map_or("NA".to_string(), |g| format!("{g:.3}")),
        format!("{:.2}", h.mass_err_ppm),
        p.peptide.length().to_string(),
    ]);
    cols
}

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

#[allow(clippy::too_many_arguments)]
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
            andes_xl::rollup::LinkClass::of(c.intra, c.homotypic)
                .name()
                .to_string(),
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
    linker: &Crosslinker,
) -> Option<crate::cli::ModList> {
    if explicit.is_some() {
        return None;
    }
    let mono = water_monolink_mods(linker);
    if mono.is_empty() {
        return None;
    }
    let mut mods = crate::model_select::default_variable_mods();
    mods.extend(mono);
    Some(crate::cli::ModList(mods))
}

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
    report_linker(&linker);
    let t = std::time::Instant::now();
    let mut cfg = XlConfig::new(linker.clone());
    let table = XlCandidates::from_out_of_core(prepared, &linker, cfg.open_cleavable);
    let n_decoy = table.peptides().iter().filter(|p| p.is_decoy).count();
    eprintln!(
        "[xl] {} linkable peptidoforms ({} decoy){} [{:.2}s]",
        table.peptides().len(),
        n_decoy,
        if table.has_open_index() {
            "; open-mass fragment index built (non-cleavable linker)"
        } else {
            ""
        },
        t.elapsed().as_secs_f64()
    );
    let spectra_for_xl: &[Spectrum] = if cli.max_spectra > 0 {
        &spectra[..spectra.len().min(cli.max_spectra)]
    } else {
        spectra
    };
    cfg.oracle_spectra = cli.xl_oracle;
    if let Some((lo, hi)) = cli.isotope_error {
        cfg.isotope_offsets = lo..=hi;
    }
    let t = std::time::Instant::now();
    let results = xl_search_run(spectra_for_xl, prepared, &table, &cfg);
    eprintln!(
        "[xl] searched {} spectra, {} CSMs [{:.2}s]",
        results.len(),
        results.iter().filter(|r| r.csm.is_some()).count(),
        t.elapsed().as_secs_f64()
    );
    report(&results);
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
    if cli.xl_rescore_intra {
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

#[cfg(test)]
mod tests {
    use super::*;

    fn linker(spec: &str) -> Crosslinker {
        crate::cli::parse_xl(spec).unwrap()
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

    #[test]
    fn an_explicit_variable_mod_list_replaces_the_monolink_default() {
        let x = linker("DSS");
        let d = xl_default_variable_mods(None, &x).unwrap().0;
        assert_eq!(
            d.len(),
            crate::model_select::default_variable_mods().len() + 2
        );
        assert!(d.iter().any(|m| m.name == "Xlink:DSS[156]"));
        let explicit = crate::cli::parse_variable_mods("Oxidation (M)").unwrap();
        assert!(xl_default_variable_mods(Some(&explicit), &x).is_none());
        let off = crate::cli::parse_variable_mods("").unwrap();
        assert!(off.0.is_empty());
        assert!(xl_default_variable_mods(Some(&off), &x).is_none());
        assert!(xl_default_variable_mods(None, &linker("CDI")).is_none());
    }
}
