//! The modification set of the `--refine` pass and the modification count and class of a
//! pass-2 PSM: the base set plus the refinement tier's modifications that do not clash with a
//! fixed modification or repeat one already present.

use std::collections::HashSet;

use model::aa_set::{AminoAcidSet, AminoAcidSetBuilder};
use model::modification::{ModLocation, Modification, ResidueSpec};
use model::peptide::Peptide;

use crate::refine_config::RefineConfig;

fn parse_location(loc: &str) -> ModLocation {
    match loc.trim().to_ascii_lowercase().as_str() {
        "n_term" | "n-term" | "nterm" => ModLocation::NTerm,
        "c_term" | "c-term" | "cterm" => ModLocation::CTerm,
        "protein_n_term" | "prot-n-term" | "prot_n_term" => ModLocation::ProtNTerm,
        "protein_c_term" | "prot-c-term" | "prot_c_term" => ModLocation::ProtCTerm,
        _ => ModLocation::Anywhere,
    }
}

fn parse_residue_spec(res: &str) -> ResidueSpec {
    let b = res.as_bytes();
    if b.len() == 1 && b[0].is_ascii_uppercase() {
        ResidueSpec::Specific(b[0])
    } else {
        ResidueSpec::Wildcard
    }
}

pub fn refinement_aa_set(
    base: &AminoAcidSet,
    cfg: &RefineConfig,
    high_res: bool,
) -> Result<AminoAcidSet, String> {
    let mut builder = AminoAcidSetBuilder::new_standard();

    let mod_key = |residue: ResidueSpec, location: ModLocation, delta: f64| {
        (residue, location, (delta * 1e4).round() as i64)
    };
    let mut present: HashSet<(ResidueSpec, ModLocation, i64)> = HashSet::new();
    let mut fixed_slots: Vec<(ResidueSpec, ModLocation)> = Vec::new();

    for m in base.distinct_mods() {
        if m.annotated_site_mass.is_some() {
            continue;
        }
        present.insert(mod_key(m.residue, m.location, m.mass_delta));
        if m.fixed {
            fixed_slots.push((m.residue, m.location));
            builder = builder.add_fixed_mod((*m).clone());
        } else {
            builder = builder.add_variable_mod((*m).clone());
        }
    }

    let residue_overlaps = |a: ResidueSpec, b: ResidueSpec| {
        a == ResidueSpec::Wildcard || b == ResidueSpec::Wildcard || a == b
    };

    for rm in &cfg.mods {
        if !high_res && rm.class == "deamidation" {
            continue;
        }
        let location = parse_location(&rm.location);
        for res in &rm.residues {
            let residue = parse_residue_spec(res);
            if fixed_slots
                .iter()
                .any(|&(fr, fl)| fl == location && residue_overlaps(fr, residue))
            {
                continue;
            }
            if !present.insert(mod_key(residue, location, rm.delta)) {
                continue;
            }
            let m = Modification {
                name: rm.name.clone(),
                mass_delta: rm.delta,
                residue,
                location,
                fixed: false,
                accession: None,
                neutral_losses: Vec::new(),
                loss_class: 0,
                annotated_site_mass: None,
                blocks_cleavage: false,
                parts: Vec::new(),
                diagnostic_ions: Vec::new(),
            };
            builder = builder.add_variable_mod(m);
        }
    }

    builder
        .with_label_channels(base.label_channels().to_vec())
        .build()
        .map_err(|e| format!("refinement_aa_set: {e}"))
}

pub fn mod_count_and_class(
    peptide: &Peptide,
    cfg: &RefineConfig,
    aa_set: &AminoAcidSet,
) -> (u32, u32) {
    let mut num_mods = 0u32;
    let mut class_id = 0u32;

    for aa in &peptide.residues {
        let Some(m) = aa.mod_ else { continue };
        let single = [m];
        let parts: &[&Modification] = if m.parts.is_empty() {
            aa_set.stack_parts(m).unwrap_or(&single)
        } else {
            &m.parts
        };
        for part in parts {
            if part.fixed {
                continue;
            }
            let matched = cfg.mods.iter().find(|rm| {
                (rm.delta - part.mass_delta).abs() < 1e-4
                    && parse_location(&rm.location) == part.location
                    && rm
                        .residues
                        .iter()
                        .any(|r| r == "*" || (r.len() == 1 && r.as_bytes()[0] == aa.residue))
            });
            let Some(rm) = matched else { continue };
            num_mods += 1;
            let this_class = RefineConfig::class_id(&rm.class);
            let pending = class_id == 0 || class_id == 99;
            if pending && (class_id == 0 || this_class != 99) {
                class_id = this_class;
            }
        }
    }

    (num_mods, class_id)
}
