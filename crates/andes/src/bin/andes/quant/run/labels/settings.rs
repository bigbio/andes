//! Label settings of a run: the label spec, the modifications the search adds for it, and the
//! number of label sites a peptide can carry.

use quant::label::{LabelSpec, Site};

#[derive(Debug, Clone)]
pub(crate) struct LabelSettings {
    pub spec: LabelSpec,
    pub fdr: f64,
    pub tol_ppm: f64,
    pub mbr_fdr: Option<f64>,
    pub arg_pro: bool,
}

impl LabelSettings {
    pub(crate) fn describe(&self) -> String {
        format!(
            "labelled pairs {} (tol {} ppm, match-between-channels q<={}{})",
            self.spec.describe(),
            self.tol_ppm,
            self.fdr,
            match self.mbr_fdr {
                Some(q) => format!(", match between runs at transfer q<={q}"),
                None => String::new(),
            }
        ) + if self.arg_pro {
            ", arginine-to-proline conversion fitted"
        } else {
            ""
        }
    }
}

pub(crate) fn search_label_mods(spec: &LabelSpec) -> (Vec<model::Modification>, Vec<(String, u8)>) {
    use model::modification::{ModLocation, ResidueSpec};
    let mut mods = Vec::new();
    let mut channels = Vec::new();
    for (c, l) in spec.search_labels() {
        let (residue, location) = match l.site {
            Site::Residue(r) => (ResidueSpec::Specific(r), ModLocation::Anywhere),
            Site::NTerm => (ResidueSpec::Wildcard, ModLocation::NTerm),
        };
        mods.push(model::Modification {
            name: l.name.clone(),
            mass_delta: l.delta,
            residue,
            location,
            fixed: false,
            accession: l.unimod.map(|u| format!("UNIMOD:{u}")),
            neutral_losses: Vec::new(),
            loss_class: 0,
            annotated_site_mass: None,
            blocks_cleavage: false,
            diagnostic_ions: Vec::new(),
            parts: Vec::new(),
        });
        if !channels.iter().any(|(n, _): &(String, u8)| *n == l.name) {
            channels.push((l.name.clone(), c as u8 + 1));
        }
    }
    (mods, channels)
}

pub(crate) fn max_label_sites(
    spec: &LabelSpec,
    enzymes: &[model::enzyme::Enzyme],
    num_tolerable_termini: u8,
    missed_cleavages: u32,
    db: &model::ProteinDb,
    max_length: u32,
) -> u32 {
    use model::enzyme::Enzyme;
    let sites = spec.sites();
    let specific = num_tolerable_termini > 0
        && !enzymes.is_empty()
        && enzymes
            .iter()
            .all(|e| !matches!(e, Enzyme::NonSpecific | Enzyme::NoCleavage));
    let mut cleavage = false;
    let mut free: Vec<u8> = Vec::new();
    for s in &sites {
        if let Site::Residue(r) = s {
            if specific && enzymes.iter().all(|e| e.is_cleavable_after(*r)) {
                cleavage = true;
            } else {
                free.push(*r);
            }
        }
    }
    let nterm = sites.contains(&Site::NTerm);
    let bound = if cleavage { missed_cleavages + 1 } else { 0 };
    let free = if free.is_empty() {
        0
    } else {
        most_in_window(db, &free, max_length as usize)
    };
    bound + free + u32::from(nterm)
}

fn most_in_window(db: &model::ProteinDb, residues: &[u8], len: usize) -> u32 {
    let len = len.max(1);
    let hit = |r: &u8| residues.contains(&r.to_ascii_uppercase());
    let mut best = 0u32;
    for p in &db.proteins {
        let s = &p.sequence;
        let mut n = 0u32;
        for (i, r) in s.iter().enumerate() {
            n += u32::from(hit(r));
            if i >= len {
                n -= u32::from(hit(&s[i - len]));
            }
            best = best.max(n);
        }
    }
    best
}
