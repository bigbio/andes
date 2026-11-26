//! Native entrapment diagnostics after rescoring.

use std::collections::HashMap;
use std::io::{BufRead, BufReader};
use std::path::Path;

use output::fdr::entrapment::{self, Discovery, Estimate, Level};
use search::decoy::entrapment::{is_entrapment_accession, EntrapmentLookup};

use crate::cli::Cli;

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
        pin,
        precursor_pin,
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
    let charges = charges(pin, rows)?;
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
            && !lookup.is_target_peptide(peptide.as_bytes());
        let pair = peptide_pair(&lookup, &proteins, &peptide, entrapment);
        let discovery = Discovery {
            key: peptide.replace('I', "L"),
            pair,
            entrapment,
            score,
            q_value: p.q_value,
        };
        peptides.push(discovery.clone());
        let charge = charges
            .get(&p.psm_id)
            .ok_or_else(|| format!("{}: no charge for rescored row {}", pin.display(), p.psm_id))?;
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
    let mut summary = format!("Entrapment  experimental diagnostics at q <= {q}\n");
    for level in &levels {
        let e = entrapment::at_q(&level.rows, &level.curve, q);
        summary.push_str(&format!(
            "  {}: targets {}  entrapment {}  r {:.3}\n\
             \x20 FDP lower {:.3}%  combined {:.3}%  paired {:.3}%\n",
            level.name,
            output::statistics::thousands(e.targets),
            output::statistics::thousands(e.entrapment),
            level.r,
            e.lower_bound * 100.0,
            e.combined * 100.0,
            e.paired * 100.0,
        ));
    }
    summary.push_str(
        "  Peptide/precursor counts use PSM q, not independent level FDR.\n\
                        \x20 PSM/precursor combined FDP uses the digest peptide ratio as a proxy.\n\
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

fn charges(
    pin: &Path,
    rows: &HashMap<String, output::PercolatorPsm>,
) -> Result<HashMap<String, u8>, String> {
    let err = |e: &dyn std::fmt::Display| format!("{}: {e}", pin.display());
    let mut reader = BufReader::new(std::fs::File::open(pin).map_err(|e| err(&e))?);
    let mut line = String::new();
    reader.read_line(&mut line).map_err(|e| err(&e))?;
    let header: Vec<&str> = line.trim_end().split('\t').collect();
    let id = header
        .iter()
        .position(|s| *s == "SpecId")
        .ok_or_else(|| err(&"no SpecId column"))?;
    let charge_cols: Vec<(usize, u8)> = header
        .iter()
        .enumerate()
        .filter_map(|(i, c)| Some((i, c.strip_prefix("charge")?.parse().ok()?)))
        .collect();
    let mut out = HashMap::with_capacity(rows.len());
    loop {
        line.clear();
        if reader.read_line(&mut line).map_err(|e| err(&e))? == 0 {
            break;
        }
        let fields: Vec<&str> = line.trim_end().split('\t').collect();
        let Some(key) = fields.get(id).filter(|s| rows.contains_key(**s)) else {
            continue;
        };
        let active: Vec<u8> = charge_cols
            .iter()
            .filter_map(|(c, z)| (fields.get(*c) == Some(&"1")).then_some(*z))
            .collect();
        if active.len() != 1 || active[0] == 0 {
            return Err(err(&format!("row {key} has no unique positive charge")));
        }
        if out.insert((*key).to_string(), active[0]).is_some() {
            return Err(err(&format!("duplicate SpecId {key}")));
        }
    }
    Ok(out)
}

#[cfg(test)]
#[path = "tests/entrapment.rs"]
mod tests;
