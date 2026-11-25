use std::io::Write;
use std::path::{Path, PathBuf};

use crate::cli::Cli;
use model::Spectrum;
use search::crosslinks::{Csm, PeptideHit, XlCandidates, XlPeptide, XlSpectrumResult};
use xl::crosslinker::Crosslinker;
use xl::crosslinker::{LinkForm, MAX_STUBS};
use xl::decoy::PairClass;
use xl::fdr::ClassCounts;
use xl::rollup::{as_items, group_by, split_q_values, Group, Item, LinkClass};

pub(crate) struct XlRun {
    pub linker: Crosslinker,
    pub table: XlCandidates,
}
use search::PreparedSearch;

const FEATURES: [&str; 46] = [
    "pair_score",
    "alpha_score",
    "beta_score",
    "score_sum",
    "score_diff",
    "pair_gap",
    "pair_gap_missing",
    "alpha_gap",
    "alpha_gap_missing",
    "beta_gap",
    "beta_gap_missing",
    "alpha_stub",
    "alpha_intact",
    "beta_stub",
    "beta_intact",
    "alpha_sites",
    "beta_sites",
    "alpha_len",
    "beta_len",
    "min_len",
    "abs_precursor_ppm",
    "alpha_abs_ppm",
    "beta_abs_ppm",
    "abs_ppm_diff",
    "isotope",
    "full_pair",
    "complement_stub_ions",
    "log_mass_pairs",
    "n_doublet_masses",
    "alpha_doublet_tic",
    "beta_doublet_tic",
    "doublet_tic_product",
    "alpha_doublet_iso",
    "beta_doublet_iso",
    "alpha_doublet_charges",
    "beta_doublet_charges",
    "charge3",
    "charge4",
    "charge5",
    "charge6plus",
    "alpha_no_stubs",
    "beta_no_stubs",
    "min_no_stubs",
    "alpha_decoupled",
    "beta_decoupled",
    "min_decoupled",
];

/// The PIN features that compare a CSM with the next-best candidate of its spectrum. With
/// them a rescorer separates target-decoy from false target-target inter-links, so the
/// inter-link model leaves them out.
pub(crate) const GAP_FEATURES: [&str; 6] = [
    "pair_gap",
    "pair_gap_missing",
    "alpha_gap",
    "alpha_gap_missing",
    "beta_gap",
    "beta_gap_missing",
];

/// The fragment-intensity features of a CSM the model has no row for: zero, marked missing.
const MS2_NONE: crate::xl::ms2::Ms2Row = [0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0];

/// The PIN `SpecId` of a spectrum's crosslink row; the CSM table repeats it in `spec_id`, so a
/// rescored PIN joins back to the table by it.
pub(crate) fn spec_id(r: &XlSpectrumResult) -> String {
    format!("xl{}", r.spec_idx)
}

fn best_stub(h: &PeptideHit) -> f64 {
    h.form_scores[..MAX_STUBS]
        .iter()
        .copied()
        .filter(|v| v.is_finite())
        .fold(f32::NEG_INFINITY, f32::max)
        .max(-1000.0) as f64
}

fn features(r: &XlSpectrumResult, c: &Csm, table: &XlCandidates) -> Vec<f64> {
    let pa = &table.peptides()[c.alpha.pep as usize];
    let pb = &table.peptides()[c.beta.pep as usize];
    let la = pa.peptide.length() as f64;
    let lb = pb.peptide.length() as f64;
    let gap = |g: Option<f32>| (g.map_or(0.0, |v| v as f64), u8::from(g.is_none()) as f64);
    let (pg, pgm) = gap(c.second_pair_score.map(|s| c.pair_score - s));
    let (ag, agm) = gap(c.alpha.gap_next);
    let (bg, bgm) = gap(c.beta.gap_next);
    let tic = r.doublets.tic.max(1.0);
    let dbl = |d: &Option<xl::doublet::DoubletMass>| match d {
        Some(d) => (
            d.intensity / tic,
            d.isotope_confirmed as u8 as f64,
            d.n_charges as f64,
        ),
        None => (0.0, 0.0, 0.0),
    };
    let (ad, ai, ac) = dbl(&c.alpha_doublet);
    let (bd, bi, bc) = dbl(&c.beta_doublet);
    let z = r.charge;
    let rs = &c.rescores;
    vec![
        c.pair_score as f64,
        c.alpha.score as f64,
        c.beta.score as f64,
        (c.alpha.score + c.beta.score) as f64,
        (c.alpha.score - c.beta.score) as f64,
        pg,
        pgm,
        ag,
        agm,
        bg,
        bgm,
        best_stub(&c.alpha),
        c.alpha.form_scores[LinkForm::Intact.index()].max(-1000.0) as f64,
        best_stub(&c.beta),
        c.beta.form_scores[LinkForm::Intact.index()].max(-1000.0) as f64,
        c.alpha.n_sites as f64,
        c.beta.n_sites as f64,
        la,
        lb,
        la.min(lb),
        c.precursor_err_ppm.abs(),
        c.alpha.mass_err_ppm.abs(),
        c.beta.mass_err_ppm.abs(),
        (c.alpha.mass_err_ppm - c.beta.mass_err_ppm).abs(),
        c.iso as f64,
        u8::from(c.kind == search::crosslinks::PairKind::Full) as f64,
        c.b_stub_ions as f64,
        (c.n_mass_pairs as f64).ln_1p(),
        r.doublets.masses.len() as f64,
        ad,
        bd,
        ad * bd,
        ai,
        bi,
        ac,
        bc,
        u8::from(z == 3) as f64,
        u8::from(z == 4) as f64,
        u8::from(z == 5) as f64,
        u8::from(z >= 6) as f64,
        rs.alpha_no_stubs as f64,
        rs.beta_no_stubs as f64,
        rs.alpha_no_stubs.min(rs.beta_no_stubs) as f64,
        rs.alpha_decoupled as f64,
        rs.beta_decoupled as f64,
        rs.alpha_decoupled.min(rs.beta_decoupled) as f64,
    ]
}

fn accessions<'a>(p: &XlPeptide, prepared: &'a PreparedSearch) -> Vec<(u32, &'a str)> {
    accessions_of(p, &p.occurrences, prepared)
}

fn accessions_of<'a>(
    p: &XlPeptide,
    occurrences: &[(u32, u32)],
    prepared: &'a PreparedSearch,
) -> Vec<(u32, &'a str)> {
    let idx = prepared.idx;
    occurrences
        .iter()
        .filter_map(|&(pi, start)| {
            let acc = idx.db.proteins.get(pi as usize)?.accession.as_str();
            (p.is_decoy || !idx.is_decoy_accession(acc)).then_some((start, acc))
        })
        .collect()
}

#[allow(clippy::too_many_arguments)]
pub(crate) fn write_xl_pin(
    path: &Path,
    spectra: &[Spectrum],
    results: &[XlSpectrumResult],
    runs: &[XlRun],
    prepared: &PreparedSearch,
    run_of: &dyn Fn(usize) -> String,
    intra: Option<bool>,
    ms2: Option<&[Option<crate::xl::ms2::Ms2Row>]>,
) -> std::io::Result<usize> {
    let mut w = std::io::BufWriter::new(std::fs::File::create(path)?);
    let mut head: Vec<&str> = vec![
        "SpecId", "Label", "ScanNr", "ExpMass", "CalcMass", "FileName",
    ];
    head.extend(FEATURES);
    if ms2.is_some() {
        head.extend(crate::xl::ms2::FEATURES);
    }
    head.extend(["Peptide", "Proteins"]);
    writeln!(w, "{}", head.join("\t"))?;
    let mut n = 0;
    for (ri, r) in results.iter().enumerate() {
        let Some(c) = &r.csm else { continue };
        if intra.is_some_and(|i| i != c.intra) {
            continue;
        }
        let spec = &spectra[r.spec_idx];
        let run = &runs[r.linker];
        let table = &run.table;
        let pa = &table.peptides()[c.alpha.pep as usize];
        let pb = &table.peptides()[c.beta.pep as usize];
        let label = if c.class == PairClass::TT { 1 } else { -1 };
        let calc = pa.peptide.mass() + pb.peptide.mass() + run.linker.intact;
        let f = features(r, c, table);
        debug_assert_eq!(f.len(), FEATURES.len());
        let bare = |p: &XlPeptide| String::from_utf8_lossy(&p.bare()).into_owned();
        let mut row: Vec<String> = vec![
            spec_id(r),
            label.to_string(),
            spec.scan.unwrap_or(r.spec_idx as i32).to_string(),
            format!("{:.5}", r.precursor_mass),
            format!("{calc:.5}"),
            run_of(r.spec_idx),
        ];
        row.extend(f.iter().map(|v| format!("{v:.5}")));
        if let Some(rows) = ms2 {
            let m = rows.get(ri).copied().flatten().unwrap_or(MS2_NONE);
            row.extend(m.iter().map(|v| format!("{v:.5}")));
        }
        row.push(format!(
            "-.{}{}-{}{}.-",
            bare(pa),
            c.alpha.site.pos + 1,
            bare(pb),
            c.beta.site.pos + 1
        ));
        let prots: Vec<&str> = accessions(pa, prepared)
            .into_iter()
            .chain(accessions(pb, prepared))
            .map(|(_, a)| a)
            .collect();
        row.push(prots.join("\t"));
        writeln!(w, "{}", row.join("\t"))?;
        n += 1;
    }
    w.flush()?;
    Ok(n)
}

pub(crate) fn percolator_scores(
    cli: &Cli,
    pin: &Path,
    results: &[XlSpectrumResult],
) -> Result<Vec<Option<f64>>, Box<dyn std::error::Error>> {
    let backend = output::resolve_backend(
        cli.percolator_bin.as_deref(),
        cli.percolator_docker,
        &cli.percolator_image,
    )
    .map_err(|e| format!("{e}"))?;
    let user_args: Vec<String> = cli
        .percolator_args
        .split_whitespace()
        .map(|s| s.to_string())
        .collect();
    let extra = output::percolator_args_with_threads(&user_args, cli.threads);
    let t = std::time::Instant::now();
    let targets =
        output::run_percolator_pooled(&backend, pin, &extra).map_err(|e| format!("{e}"))?;
    let decoy_text = std::fs::read_to_string(output::decoy_results_path(pin))?;
    let decoys = output::parse_psm_results(&decoy_text).map_err(|e| format!("{e}"))?;
    let mut missing = 0usize;
    let scores: Vec<Option<f64>> = results
        .iter()
        .map(|r| {
            r.csm.as_ref()?;
            let id = spec_id(r);
            let row = targets.get(&id).or_else(|| decoys.get(&id))?;
            if row.score.is_none() {
                missing += 1;
            }
            row.score
        })
        .collect();
    log::info!(
        target: "xl",
        "Percolator ({}): {} target and {} decoy rows scored, {missing} CSMs without a score \
         [{:.2}s]",
        backend.describe(),
        targets.len(),
        decoys.len(),
        t.elapsed().as_secs_f64()
    );
    Ok(scores)
}

fn rp_side(h: &PeptideHit, table: &XlCandidates, prepared: &PreparedSearch) -> String {
    let p = &table.peptides()[h.pep as usize];
    let mut v: Vec<String> = accessions_of(p, p.occurrences_at(&h.site), prepared)
        .into_iter()
        .map(|(start, acc)| format!("{acc}:{}", start as usize + h.site.pos + 1))
        .collect();
    v.sort();
    v.dedup();
    v.join(";")
}

fn ppi_side(h: &PeptideHit, table: &XlCandidates, prepared: &PreparedSearch) -> String {
    let p = &table.peptides()[h.pep as usize];
    let mut v: Vec<&str> = accessions_of(p, p.occurrences_at(&h.site), prepared)
        .into_iter()
        .map(|(_, a)| a)
        .collect();
    v.sort_unstable();
    v.dedup();
    v.join(";")
}

fn pair_key(a: String, b: String) -> String {
    if a <= b {
        format!("{a} -- {b}")
    } else {
        format!("{b} -- {a}")
    }
}

/// The residue-pair and protein-pair keys of a CSM: each side's accessions with the linked
/// residue's position (`P12345:120`), or the accessions alone, the two sides in sorted order.
pub(crate) fn pair_keys(
    c: &Csm,
    table: &XlCandidates,
    prepared: &PreparedSearch,
) -> (String, String) {
    (
        pair_key(rp_side(&c.alpha, table, prepared), rp_side(&c.beta, table, prepared)),
        pair_key(ppi_side(&c.alpha, table, prepared), ppi_side(&c.beta, table, prepared)),
    )
}

pub(crate) struct Levels {
    pub csm_items: Vec<Item>,
    pub csm_rows: Vec<usize>,
    pub csm_q: Vec<Option<f64>>,
    pub rp: Vec<Group>,
    pub rp_q: Vec<Option<f64>>,
    pub rp_q2: Vec<Option<f64>>,
    pub ppi: Vec<Group>,
    pub ppi_q: Vec<Option<f64>>,
}

pub(crate) fn roll_up(
    results: &[XlSpectrumResult],
    score: &[Option<f64>],
    runs: &[XlRun],
    prepared: &PreparedSearch,
) -> Levels {
    let mut csm_items = Vec::new();
    let mut csm_rows = Vec::new();
    let mut ppi_of_csm = Vec::new();
    for (i, (r, s)) in results.iter().zip(score).enumerate() {
        let (Some(c), Some(s)) = (&r.csm, s) else {
            continue;
        };
        let (rp_key, ppi_key) = pair_keys(c, &runs[r.linker].table, prepared);
        csm_items.push(Item {
            key: rp_key,
            score: *s,
            class: c.class,
            link: LinkClass::of(c.intra, c.homotypic),
            n_csms: 1,
        });
        ppi_of_csm.push(ppi_key);
        csm_rows.push(i);
    }
    levels_of(csm_items, ppi_of_csm, csm_rows)
}

/// The CSM, residue-pair and protein-pair levels with their q-values, from scored CSMs:
/// `csm_items` keyed by residue pair, `ppi_of_csm` the protein-pair key of each, and
/// `csm_rows` each CSM's row in the caller's own table.
pub(crate) fn levels_of(
    csm_items: Vec<Item>,
    ppi_of_csm: Vec<String>,
    csm_rows: Vec<usize>,
) -> Levels {
    let csm_groups: Vec<Group> = csm_items
        .iter()
        .enumerate()
        .map(|(i, it)| Group {
            key: it.key.clone(),
            score: it.score,
            class: it.class,
            link: it.link,
            n_members: 1,
            n_csms: 1,
            best: i,
        })
        .collect();
    let csm_q = split_q_values(&csm_groups, |_| true);
    let rp = group_by(&csm_items, |_, it| it.key.clone());
    let rp_q = split_q_values(&rp, |_| true);
    let rp_q2 = split_q_values(&rp, |g| g.n_csms >= 2);
    let rp_items = as_items(&rp, |g| ppi_of_csm[g.best].clone());
    let ppi = group_by(&rp_items, |_, it| it.key.clone());
    let ppi_q = split_q_values(&ppi, |g| g.link == LinkClass::Inter);
    Levels {
        csm_items,
        csm_rows,
        csm_q,
        rp,
        rp_q,
        rp_q2,
        ppi,
        ppi_q,
    }
}

fn counts_line(classes: impl Iterator<Item = (PairClass, LinkClass, Option<f64>)>) -> String {
    let rows: Vec<(PairClass, LinkClass, Option<f64>)> = classes.collect();
    let mut parts = Vec::new();
    for link in LinkClass::ALL {
        if !rows.iter().any(|r| r.1 == link) {
            continue;
        }
        let label = link.name();
        for (cut_label, cut) in [("all", f64::INFINITY), ("q<=5%", 0.05), ("q<=1%", 0.01)] {
            let mut c = ClassCounts::default();
            for (class, l, q) in &rows {
                if *l == link && (cut.is_infinite() || q.is_some_and(|v| v <= cut)) {
                    c.add(*class);
                }
            }
            parts.push(format!(
                "{label} {cut_label}: TT {} TD {} DD {}",
                c.tt, c.td, c.dd
            ));
        }
    }
    parts.join("; ")
}

pub(crate) fn report_levels(l: &Levels) {
    log::info!(
        target: "xl",
        "CSM: {}",
        counts_line(
            l.csm_items
                .iter()
                .zip(&l.csm_q)
                .map(|(it, q)| (it.class, it.link, *q))
        )
    );
    log::info!(
        target: "xl",
        "residue pair: {}",
        counts_line(l.rp.iter().zip(&l.rp_q).map(|(g, q)| (g.class, g.link, *q)))
    );
    log::info!(
        target: "xl",
        "residue pair (>=2 CSMs): {}",
        counts_line(
            l.rp.iter()
                .zip(&l.rp_q2)
                .filter(|(g, _)| g.n_csms >= 2)
                .map(|(g, q)| (g.class, g.link, *q))
        )
    );
    log::info!(
        target: "xl",
        "protein pair (inter only): {}",
        counts_line(
            l.ppi
                .iter()
                .zip(&l.ppi_q)
                .filter(|(g, _)| g.link == LinkClass::Inter)
                .map(|(g, q)| (g.class, g.link, *q))
        )
    );
}

fn fmt_q(q: Option<f64>) -> String {
    q.map_or("NA".to_string(), |v| format!("{v:.5}"))
}

/// Write the residue-pair and protein-pair tables; `title_of` names the spectrum of the CSM
/// at a `csm_rows` value.
pub(crate) fn write_levels(
    rp_path: &PathBuf,
    ppi_path: &PathBuf,
    l: &Levels,
    title_of: &dyn Fn(usize) -> String,
) -> std::io::Result<()> {
    let mut w = std::io::BufWriter::new(std::fs::File::create(rp_path)?);
    writeln!(
        w,
        "residue_pair\tclass\tlink\tscore\tq_value\tq_value_min2\tn_csms\tbest_spectrum"
    )?;
    for (g, (q, q2)) in l.rp.iter().zip(l.rp_q.iter().zip(&l.rp_q2)) {
        writeln!(
            w,
            "{}\t{}\t{}\t{:.4}\t{}\t{}\t{}\t{}",
            g.key,
            g.class.name(),
            g.link.name(),
            g.score,
            fmt_q(*q),
            fmt_q(*q2),
            g.n_csms,
            title_of(l.csm_rows[g.best])
        )?;
    }
    w.flush()?;
    let mut w = std::io::BufWriter::new(std::fs::File::create(ppi_path)?);
    writeln!(
        w,
        "protein_pair\tclass\tlink\tscore\tq_value\tn_residue_pairs\tn_csms\tbest_spectrum"
    )?;
    for (g, q) in l.ppi.iter().zip(&l.ppi_q) {
        writeln!(
            w,
            "{}\t{}\t{}\t{:.4}\t{}\t{}\t{}\t{}",
            g.key,
            g.class.name(),
            g.link.name(),
            g.score,
            fmt_q(*q),
            g.n_members,
            g.n_csms,
            title_of(l.csm_rows[l.rp[g.best].best])
        )?;
    }
    w.flush()
}
