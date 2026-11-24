//! The glycopeptide quantification tables: per-feature rows and the composition, site and
//! protein roll-ups.

use std::collections::{BTreeSet, HashMap};
use std::io::Write;
use std::path::{Path, PathBuf};

use glyco::glycan::db::GlycanComp;
use glyco::glycan::first::GlycanCore;
use search::candidates::Candidate;
use search::index::SearchIndex;

use super::{accepted, Evidence, FeatureRow, Settings};

pub(super) fn output_path(pin: &Path, kind: &str) -> PathBuf {
    let s = pin.to_string_lossy();
    match s.strip_suffix(".glyco.pin") {
        Some(base) => PathBuf::from(format!("{base}.glyco.{kind}.tsv")),
        None => PathBuf::from(format!("{s}.glyco.{kind}.tsv")),
    }
}

pub(super) fn write_feature_table(
    path: &Path,
    rows: &[FeatureRow],
    settings: &Settings,
) -> Result<(), String> {
    let mut out = std::io::BufWriter::new(
        std::fs::File::create(path).map_err(|e| format!("{}: {e}", path.display()))?,
    );
    writeln!(
        out,
        "Run\tSpecId\tScan\tPeptide\tGlycosite\tGlycan\tTopology\tProteins\tCharge\t\
         IsotopeOffset\tAnchorRT\tApexRT\tRTStart\tRTStop\tRTDelta_s\tMonoMz\tTheoMz\tPpm\t\
         Area\tMonoArea\tApexIntensity\tCosine\tNIsotopesFound\tNScans\tFeatureScore\t\
         FeatureQValue\tDecoyCosine\tAccepted\tEvidence\tMatchScore"
    )
    .map_err(|e| format!("{}: {e}", path.display()))?;
    for row in rows {
        let t = &row.target;
        let (apex, start, stop, area, mono_area, apex_int, cos, n_found, n_scans, delta) =
            match &row.feature {
                Some(f) => (
                    f.apex_rt,
                    f.rt_start,
                    f.rt_stop,
                    f.area,
                    f.isotope_areas.first().copied().unwrap_or(0.0),
                    f.apex_intensity as f64,
                    f.cosine as f64,
                    f.n_isotopes_found,
                    f.n_scans,
                    f.rt_delta_s,
                ),
                None => (
                    f64::NAN,
                    f64::NAN,
                    f64::NAN,
                    f64::NAN,
                    f64::NAN,
                    f64::NAN,
                    f64::NAN,
                    0,
                    0,
                    f64::NAN,
                ),
            };
        let decoy_cos = row.decoy.as_ref().map_or(f64::NAN, |f| f.cosine as f64);
        writeln!(
            out,
            "{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t\
             {}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}",
            t.run,
            t.spec_id,
            t.scan,
            t.peptide,
            t.site,
            t.glycan,
            t.topology,
            t.proteins,
            t.charge,
            t.isotope_offset,
            fmt(t.anchor_rt),
            fmt(apex),
            fmt(start),
            fmt(stop),
            fmt(delta),
            fmt(t.mono_mz),
            fmt(t.theo_mz),
            fmt(t.ppm),
            fmt(area),
            fmt(mono_area),
            fmt(apex_int),
            fmt(cos),
            n_found,
            n_scans,
            fmt(row.score.map_or(f64::NAN, f64::from)),
            fmt(row.q.unwrap_or(f64::NAN)),
            fmt(decoy_cos),
            u8::from(accepted(row, settings)),
            t.evidence.as_str(),
            fmt(row.match_score.map_or(f64::NAN, f64::from)),
        )
        .map_err(|e| format!("{}: {e}", path.display()))?;
    }
    out.flush().map_err(|e| format!("{}: {e}", path.display()))
}

pub(super) struct Kept {
    pub(super) glycopeptides: usize,
    pub(super) sites: usize,
}

pub(super) fn write_rollup_tables(
    pin: &Path,
    rows: &[FeatureRow],
    run_names: &[&str],
    settings: &Settings,
) -> Result<Kept, String> {
    #[derive(PartialEq, Eq, Hash)]
    struct Key(String, String, String, u8);
    let mut best: HashMap<Key, usize> = HashMap::new();
    for (i, row) in rows.iter().enumerate() {
        if !accepted(row, settings) || row.feature.is_none() {
            continue;
        }
        let key = Key(
            row.target.run.clone(),
            row.target.peptide.clone(),
            row.target.glycan.clone(),
            row.target.charge,
        );
        let better = match best.get(&key) {
            None => true,
            Some(&j) => {
                let sj = rows[j].score.unwrap_or(f32::NEG_INFINITY);
                let si = row.score.unwrap_or(f32::NEG_INFINITY);
                si > sj
                    || (si == sj
                        && row.feature.as_ref().map_or(0.0, |f| f.area)
                            > rows[j].feature.as_ref().map_or(0.0, |f| f.area))
            }
        };
        if better {
            best.insert(key, i);
        }
    }

    struct Glyco {
        area: f64,
        best_cosine: f32,
        best_q: f64,
        charges: BTreeSet<u8>,
        scans: BTreeSet<i32>,
        transferred: usize,
        proteins: String,
        site: String,
    }
    let mut glyco: HashMap<(String, String, String), Glyco> = HashMap::new();
    for &i in best.values() {
        let row = &rows[i];
        let Some(f) = row.feature.as_ref() else {
            continue;
        };
        let key = (
            row.target.run.clone(),
            row.target.peptide.clone(),
            row.target.glycan.clone(),
        );
        let e = glyco.entry(key).or_insert_with(|| Glyco {
            area: 0.0,
            best_cosine: 0.0,
            best_q: 1.0,
            charges: BTreeSet::new(),
            scans: BTreeSet::new(),
            transferred: 0,
            proteins: row.target.proteins.clone(),
            site: row.target.site.clone(),
        });
        e.area += f.area;
        e.best_cosine = e.best_cosine.max(f.cosine);
        e.best_q = e.best_q.min(row.q.unwrap_or(1.0));
        e.charges.insert(row.target.charge);
        match row.target.evidence {
            Evidence::Identified => {
                e.scans.insert(row.target.scan);
            }
            Evidence::Transferred => e.transferred += 1,
        }
    }
    let quant_path = output_path(pin, "quant");
    let mut out = std::io::BufWriter::new(
        std::fs::File::create(&quant_path).map_err(|e| format!("{}: {e}", quant_path.display()))?,
    );
    writeln!(
        out,
        "Run\tPeptide\tGlycosite\tGlycan\tProteins\tCharges\tNScans\tNTransferred\tBestCosine\t\
         BestQValue\tArea"
    )
    .map_err(|e| format!("{}: {e}", quant_path.display()))?;
    let mut glyco_rows: Vec<_> = glyco.iter().collect();
    glyco_rows.sort_by(|a, b| {
        a.0 .0
            .cmp(&b.0 .0)
            .then_with(|| a.0 .1.cmp(&b.0 .1))
            .then_with(|| a.0 .2.cmp(&b.0 .2))
    });
    for ((run, peptide, glycan), g) in &glyco_rows {
        writeln!(
            out,
            "{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}",
            run,
            peptide,
            g.site,
            glycan,
            g.proteins,
            g.charges
                .iter()
                .map(|c| c.to_string())
                .collect::<Vec<_>>()
                .join(";"),
            g.scans.len(),
            g.transferred,
            fmt(g.best_cosine as f64),
            fmt(g.best_q),
            fmt(g.area),
        )
        .map_err(|e| format!("{}: {e}", quant_path.display()))?;
    }
    out.flush()
        .map_err(|e| format!("{}: {e}", quant_path.display()))?;

    let mut site_total: HashMap<(String, String), f64> = HashMap::new();
    for ((run, peptide, _), g) in &glyco_rows {
        *site_total
            .entry((run.clone(), peptide.clone()))
            .or_insert(0.0) += g.area;
    }
    let site_path = output_path(pin, "site");
    let mut out = std::io::BufWriter::new(
        std::fs::File::create(&site_path).map_err(|e| format!("{}: {e}", site_path.display()))?,
    );
    writeln!(
        out,
        "Run\tPeptide\tGlycosite\tProteins\tGlycan\tArea\tSiteGlycoArea\tGlycoformShare"
    )
    .map_err(|e| format!("{}: {e}", site_path.display()))?;
    for ((run, peptide, glycan), g) in &glyco_rows {
        let total = site_total
            .get(&(run.clone(), peptide.clone()))
            .copied()
            .unwrap_or(0.0);
        let share = if total > 0.0 {
            g.area / total
        } else {
            f64::NAN
        };
        writeln!(
            out,
            "{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}",
            run,
            peptide,
            g.site,
            g.proteins,
            glycan,
            fmt(g.area),
            fmt(total),
            fmt(share),
        )
        .map_err(|e| format!("{}: {e}", site_path.display()))?;
    }
    out.flush()
        .map_err(|e| format!("{}: {e}", site_path.display()))?;

    struct MatrixRow {
        site: String,
        proteins: String,
        areas: HashMap<String, f64>,
    }
    let mut matrix: HashMap<(String, String), MatrixRow> = HashMap::new();
    for ((run, peptide, glycan), g) in &glyco_rows {
        let e = matrix
            .entry((peptide.clone(), glycan.clone()))
            .or_insert_with(|| MatrixRow {
                site: g.site.clone(),
                proteins: g.proteins.clone(),
                areas: HashMap::new(),
            });
        e.areas.insert(run.clone(), g.area);
    }
    let matrix_path = output_path(pin, "matrix");
    let mut out = std::io::BufWriter::new(
        std::fs::File::create(&matrix_path)
            .map_err(|e| format!("{}: {e}", matrix_path.display()))?,
    );
    write!(out, "Peptide\tGlycosite\tGlycan\tProteins\tNRuns")
        .map_err(|e| format!("{}: {e}", matrix_path.display()))?;
    for name in run_names {
        write!(out, "\t{name}").map_err(|e| format!("{}: {e}", matrix_path.display()))?;
    }
    writeln!(out).map_err(|e| format!("{}: {e}", matrix_path.display()))?;
    let mut matrix_rows: Vec<_> = matrix.iter().collect();
    matrix_rows.sort_by(|a, b| a.0 .0.cmp(&b.0 .0).then_with(|| a.0 .1.cmp(&b.0 .1)));
    for ((peptide, glycan), m) in matrix_rows {
        write!(
            out,
            "{}\t{}\t{}\t{}\t{}",
            peptide,
            m.site,
            glycan,
            m.proteins,
            m.areas.len()
        )
        .map_err(|e| format!("{}: {e}", matrix_path.display()))?;
        for name in run_names {
            match m.areas.get(*name) {
                Some(area) => write!(out, "\t{}", fmt(*area)),
                None => write!(out, "\t"),
            }
            .map_err(|e| format!("{}: {e}", matrix_path.display()))?;
        }
        writeln!(out).map_err(|e| format!("{}: {e}", matrix_path.display()))?;
    }
    out.flush()
        .map_err(|e| format!("{}: {e}", matrix_path.display()))?;

    let mut per_protein: HashMap<(String, String), Vec<f64>> = HashMap::new();
    for ((run, _, _), g) in &glyco_rows {
        for acc in g.proteins.split(';').filter(|s| !s.is_empty()) {
            per_protein
                .entry((run.clone(), acc.to_string()))
                .or_default()
                .push(g.area);
        }
    }
    let protein_path = output_path(pin, "protein");
    let mut out = std::io::BufWriter::new(
        std::fs::File::create(&protein_path)
            .map_err(|e| format!("{}: {e}", protein_path.display()))?,
    );
    writeln!(out, "Run\tProtein\tGlycopeptides\tSummedArea\tTop3")
        .map_err(|e| format!("{}: {e}", protein_path.display()))?;
    let mut protein_matrix: HashMap<String, HashMap<String, (f64, f64)>> = HashMap::new();
    for ((run, acc), areas) in &per_protein {
        let mut sorted = areas.clone();
        sorted.sort_by(|a, b| b.total_cmp(a));
        let top3: f64 = sorted.iter().take(3).sum();
        protein_matrix
            .entry(acc.clone())
            .or_default()
            .insert(run.clone(), (areas.iter().sum::<f64>(), top3));
    }
    let protein_matrix_path = output_path(pin, "protein_matrix");
    let mut pout = std::io::BufWriter::new(
        std::fs::File::create(&protein_matrix_path)
            .map_err(|e| format!("{}: {e}", protein_matrix_path.display()))?,
    );
    write!(pout, "Protein").map_err(|e| format!("{}: {e}", protein_matrix_path.display()))?;
    for name in run_names {
        write!(pout, "\t{name}_SummedArea\t{name}_Top3")
            .map_err(|e| format!("{}: {e}", protein_matrix_path.display()))?;
    }
    writeln!(pout).map_err(|e| format!("{}: {e}", protein_matrix_path.display()))?;
    let mut protein_matrix_rows: Vec<_> = protein_matrix.iter().collect();
    protein_matrix_rows.sort_by(|a, b| a.0.cmp(b.0));
    for (acc, per_run) in protein_matrix_rows {
        write!(pout, "{acc}").map_err(|e| format!("{}: {e}", protein_matrix_path.display()))?;
        for name in run_names {
            match per_run.get(*name) {
                Some((sum, top3)) => write!(pout, "\t{}\t{}", fmt(*sum), fmt(*top3)),
                None => write!(pout, "\t\t"),
            }
            .map_err(|e| format!("{}: {e}", protein_matrix_path.display()))?;
        }
        writeln!(pout).map_err(|e| format!("{}: {e}", protein_matrix_path.display()))?;
    }
    pout.flush()
        .map_err(|e| format!("{}: {e}", protein_matrix_path.display()))?;

    let mut protein_rows: Vec<_> = per_protein.into_iter().collect();
    protein_rows.sort_by(|a, b| a.0 .0.cmp(&b.0 .0).then_with(|| a.0 .1.cmp(&b.0 .1)));
    for ((run, acc), mut areas) in protein_rows {
        let total: f64 = areas.iter().sum();
        areas.sort_by(|a, b| b.total_cmp(a));
        let top3: f64 = areas.iter().take(3).sum();
        writeln!(
            out,
            "{}\t{}\t{}\t{}\t{}",
            run,
            acc,
            areas.len(),
            fmt(total),
            fmt(top3),
        )
        .map_err(|e| format!("{}: {e}", protein_path.display()))?;
    }
    out.flush()
        .map_err(|e| format!("{}: {e}", protein_path.display()))?;

    Ok(Kept {
        glycopeptides: glyco_rows.len(),
        sites: site_total.len(),
    })
}

pub(super) fn composition_label(comp: &GlycanComp) -> String {
    format!(
        "HexNAc{}Hex{}Fuc{}NeuAc{}NeuGc{}",
        comp.hexnac, comp.hex, comp.fuc, comp.neuac, comp.neugc
    )
}

pub(super) fn site_label(residues: &[u8], post: u8, core: GlycanCore) -> String {
    match core {
        GlycanCore::NGlycan => match glyco::sequon::glycosites(residues, post).as_slice() {
            [only] => format!("N{}", only + 1),
            _ => "N?".to_string(),
        },
        GlycanCore::OGlycan => match glyco::sequon::all_ser_thr_sites(residues).as_slice() {
            [only] => format!("{}{}", residues[*only] as char, only + 1),
            _ => "S/T?".to_string(),
        },
    }
}

pub(super) fn proteins_of(
    psm: &search::psm::PsmMatch,
    candidates: &[Candidate],
    index: &SearchIndex,
) -> String {
    psm.candidate_idxs
        .iter()
        .filter_map(|&c| candidates.get(c as usize))
        .map(|c| {
            index.protein_at(c.protein_index).map_or_else(
                || format!("PROT_{}", c.protein_index),
                |p| p.accession.clone(),
            )
        })
        .collect::<Vec<_>>()
        .join(";")
}

fn fmt(v: f64) -> String {
    if v.is_finite() {
        format!("{v:.6}")
    } else {
        String::new()
    }
}
