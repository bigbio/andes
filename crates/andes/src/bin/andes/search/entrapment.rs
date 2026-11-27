//! Native entrapment diagnostics after rescoring.

use std::collections::HashMap;
use std::path::Path;

use output::fdr::entrapment::{self, Discovery, Estimate, Level};
use search::decoy::entrapment::{EntrapmentLookup, is_entrapment_accession};

use crate::cli::Cli;

mod shared;

/// Percentage points per unit fraction in the diagnostic summary.
const PERCENT_PER_FRACTION: f64 = 100.0;

struct OwnedLevel {
    name: &'static str,
    r: f64,
    rows: Vec<Discovery>,
    curve: Vec<Estimate>,
}

impl OwnedLevel {
    fn new(name: &'static str, r: f64, rows: Vec<Discovery>) -> Self {
        let mut rows = entrapment::best_per_key(rows);
        let curve = entrapment::curve(&mut rows, r);
        Self {
            name,
            r,
            rows,
            curve,
        }
    }
}

pub(super) struct Inputs<'a> {
    pub cli: &'a Cli,
    pub index: &'a search::SearchIndex,
    pub params: &'a search::SearchParams,
    pub rows: Option<&'a HashMap<String, output::PercolatorPsm>>,
    pub groups: Option<&'a output::Inference>,
    pub base: &'a Path,
    pub pin: &'a Path,
    pub true_charges: &'a HashMap<String, u8>,
    pub precursor_pin: bool,
}

pub(super) fn report(input: Inputs<'_>) -> Result<Option<String>, String> {
    let Inputs {
        cli,
        index,
        params,
        rows,
        groups,
        base,
        true_charges,
        precursor_pin,
        ..
    } = input;
    if !cli.add_entrapment {
        return Ok(None);
    }
    let Some(rows) = rows else {
        return Ok(None);
    };
    let digest = search::decoy::Proteases::of(params);
    let lookup = EntrapmentLookup::new(
        index
            .db
            .proteins
            .iter()
            .filter(|p| !index.is_decoy_accession(&p.accession)),
        &digest,
    )
    .ok_or("--add-entrapment: searched database contains no entrapment proteins")?;
    let shared = shared::peptides(index, rows.values())?;
    let mut psms = Vec::with_capacity(rows.len());
    let mut peptides = Vec::with_capacity(rows.len());
    let mut precursors = Vec::with_capacity(rows.len());
    for p in rows.values() {
        let score = p.score.filter(|s| s.is_finite()).ok_or_else(|| {
            format!(
                "entrapment estimate: rescored row {} has no finite score",
                p.psm_id
            )
        })?;
        let peptide = output::offsets::stripped_sequence(&p.peptide);
        let proteins: Vec<&str> = p.proteins.split(';').filter(|s| !s.is_empty()).collect();
        if proteins.is_empty() || peptide.is_empty() || !p.q_value.is_finite() {
            return Err(format!(
                "entrapment estimate: invalid peptide, proteins or q-value for {}",
                p.psm_id
            ));
        }
        let entrapment = proteins.iter().all(|a| is_entrapment_accession(a))
            && !shared.contains(&peptide.replace('I', "L"));
        let pair = peptide_pair(&lookup, &proteins, &peptide, entrapment);
        let discovery = Discovery {
            key: peptide.replace('I', "L"),
            pair,
            entrapment,
            score,
            q_value: p.q_value,
        };
        peptides.push(discovery.clone());
        let charge = true_charges
            .get(&p.psm_id)
            .ok_or_else(|| format!("no true charge for rescored row {}", p.psm_id))?;
        precursors.push(Discovery {
            key: format!("{}|{charge}", peptidoform(&p.peptide)),
            pair: None,
            ..discovery.clone()
        });
        psms.push(Discovery {
            key: p.psm_id.clone(),
            pair: None,
            ..discovery
        });
    }
    let r = lookup.ratio();
    let mut levels = vec![
        OwnedLevel::new(
            if precursor_pin {
                "dia_precursor_rows"
            } else {
                "psm"
            },
            r,
            psms,
        ),
        OwnedLevel::new("precursor_at_psm_q", r, precursors),
        OwnedLevel::new("peptide_at_psm_q", r, peptides),
    ];
    if let Some(groups) = groups {
        let proteins = index
            .db
            .proteins
            .iter()
            .filter(|p| !index.is_decoy_accession(&p.accession));
        let (mut target_count, mut entrapment_count) = (0usize, 0usize);
        for p in proteins {
            if is_entrapment_accession(&p.accession) {
                entrapment_count += 1;
            } else {
                target_count += 1;
            }
        }
        let rows = groups
            .groups
            .iter()
            .map(|g| {
                let entrapment = !g.accessions.is_empty()
                    && g.accessions.iter().all(|a| is_entrapment_accession(a));
                Discovery {
                    key: g.lead.clone(),
                    pair: None,
                    entrapment,
                    score: g.score,
                    q_value: g.q_value,
                }
            })
            .collect();
        levels.push(OwnedLevel::new(
            "protein_group",
            entrapment_count as f64 / target_count.max(1) as f64,
            rows,
        ));
    }
    let path = base.with_extension("entrapment.parquet");
    let borrowed: Vec<_> = levels
        .iter()
        .map(|l| Level {
            name: l.name,
            r: l.r,
            rows: &l.rows,
            curve: &l.curve,
        })
        .collect();
    entrapment::write_curve(&path, &borrowed, cli.entrapment_seed)?;
    log::info!(target: "output", "Wrote entrapment diagnostics: {}", path.display());
    let q = cli.psm_fdr();
    let mut summary = format!("Entrapment  experimental q-only diagnostics at q <= {q}\n");
    if !shared.is_empty() {
        summary.push_str(&format!(
            "  {} distinct peptide sequences rescued by target-substring matches (I = L)\n",
            output::statistics::thousands(shared.len()),
        ));
    }
    for level in &levels {
        let e = entrapment::at_q(&level.rows, &level.curve, q);
        summary.push_str(&format!(
            "  {}: targets {}  entrapment {}  r {:.3}\n\
             \x20 FDP lower {:.3}%  combined {:.3}%  paired {:.3}%\n",
            level.name,
            output::statistics::thousands(e.targets),
            output::statistics::thousands(e.entrapment),
            level.r,
            e.lower_bound * PERCENT_PER_FRACTION,
            e.combined * PERCENT_PER_FRACTION,
            e.paired * PERCENT_PER_FRACTION,
        ));
    }
    summary.push_str(
        "  Peptide/precursor counts use PSM q, not independent level FDR.\n\
         \x20 Counts do not apply the additional --pep or --protein-fdr report filters.\n\
                        \x20 Combined FDP uses the searched digest ratio as a background proxy.\n\
                        \x20 NaN paired FDP means pairing is unavailable, not zero error.\n",
    );
    Ok(Some(summary))
}

fn peptide_pair(
    lookup: &EntrapmentLookup<'_>,
    proteins: &[&str],
    peptide: &str,
    entrapment: bool,
) -> Option<String> {
    if !entrapment {
        return Some(peptide.replace('I', "L"));
    }
    let mut pairs: Vec<String> = proteins
        .iter()
        .map(|a| {
            lookup
                .paired_target(a, peptide.as_bytes())
                .map(|p| String::from_utf8_lossy(p).replace('I', "L"))
        })
        .collect::<Option<Vec<_>>>()?;
    pairs.sort();
    pairs.dedup();
    (pairs.len() == 1).then(|| pairs.remove(0))
}

fn peptidoform(peptide: &str) -> String {
    let core = match (peptide.find('.'), peptide.rfind('.')) {
        (Some(a), Some(b)) if a <= 1 && b > a + 1 && b + 2 >= peptide.len() => &peptide[a + 1..b],
        _ => peptide,
    };
    core.replace('I', "L")
}

#[cfg(test)]
#[path = "tests/entrapment.rs"]
mod tests;
