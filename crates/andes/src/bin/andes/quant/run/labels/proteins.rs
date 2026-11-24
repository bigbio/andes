//! Protein ratios from the channel ratios of their peptides, and their TSV output.

use quant::label::median;

pub(crate) const MIN_RATIO_COUNT: usize = 2;

#[derive(Debug, Clone, PartialEq)]
pub(crate) struct ProteinRatio {
    pub protein: String,
    pub lead: String,
    pub run: String,
    pub channel: u8,
    pub channel_name: String,
    pub log2_ratio: f64,
    pub log2_ratio_normalized: Option<f64>,
    pub ratio_count: usize,
    pub inferred_count: usize,
}

pub(crate) fn protein_ratios(
    rows: &[output::LabelRow],
    group_of: impl Fn(&str) -> Option<(String, String)>,
    min_count: usize,
) -> Vec<ProteinRatio> {
    use std::collections::BTreeMap;
    type Key = (String, String, String, u8, String);
    type Evidence = (Vec<f64>, Vec<f64>, bool);
    let mut by: BTreeMap<Key, BTreeMap<String, Evidence>> = BTreeMap::new();
    for r in rows {
        let Some(v) = r.log2_ratio else { continue };
        let (protein, lead) = group_of(&r.sequence).unwrap_or_else(|| {
            let p = r.proteins.join(";");
            (p.clone(), r.proteins.first().cloned().unwrap_or(p))
        });
        let e = by
            .entry((
                r.run.clone(),
                protein,
                lead,
                r.channel,
                r.channel_name.clone(),
            ))
            .or_default()
            .entry(r.peptide.clone())
            .or_default();
        e.0.push(v);
        if let Some(n) = r.log2_ratio_normalized {
            e.1.push(n);
        }
        e.2 |= r.status == "identified";
    }
    let mut out = Vec::new();
    for ((run, protein, lead, channel, channel_name), peptides) in by {
        let ratios: Vec<f64> = peptides.values().map(|p| median(p.0.clone())).collect();
        if ratios.len() < min_count {
            continue;
        }
        let norm: Vec<f64> = peptides
            .values()
            .filter(|p| !p.1.is_empty())
            .map(|p| median(p.1.clone()))
            .collect();
        out.push(ProteinRatio {
            protein,
            lead,
            run,
            channel,
            channel_name,
            log2_ratio: median(ratios.clone()),
            log2_ratio_normalized: (norm.len() >= min_count).then(|| median(norm)),
            ratio_count: ratios.len(),
            inferred_count: peptides.values().filter(|p| !p.2).count(),
        });
    }
    out
}

pub(crate) fn protein_ratio_tsv(ratios: &[ProteinRatio]) -> String {
    let mut s = String::from(
        "protein_group\tlead\tfilename\tchannel\tchannel_name\tlog2_ratio\t\
         log2_ratio_normalized\tratio_count\tinferred_ratio_count\n",
    );
    for r in ratios {
        s.push_str(&format!(
            "{}\t{}\t{}\t{}\t{}\t{:.6}\t{}\t{}\t{}\n",
            r.protein,
            r.lead,
            r.run,
            r.channel,
            r.channel_name,
            r.log2_ratio,
            r.log2_ratio_normalized
                .map_or_else(|| "NA".to_string(), |v| format!("{v:.6}")),
            r.ratio_count,
            r.inferred_count
        ));
    }
    s
}
