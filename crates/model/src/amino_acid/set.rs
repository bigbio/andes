use rustc_hash::FxHashMap;

use crate::amino_acid::AminoAcid;
use crate::enzyme::Enzyme;
use crate::modification::{ModLocation, Modification, ResidueSpec};

const STANDARD_RESIDUES: &[u8] = b"ACDEFGHIKLMNPQRSTVWY";
const IMPLAUSIBLE_MASS_THRESHOLD: f64 = 5000.0;

#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
pub enum Terminus {
    No,
    Peptide,
    Protein,
}

impl Terminus {
    pub const ALL: [Terminus; 3] = [Terminus::No, Terminus::Peptide, Terminus::Protein];

    #[inline]
    pub fn rank(self) -> usize {
        match self {
            Terminus::No => 0,
            Terminus::Peptide => 1,
            Terminus::Protein => 2,
        }
    }

    #[inline]
    pub fn at(is_peptide_terminus: bool, is_protein_terminus: bool) -> Self {
        match (is_peptide_terminus, is_protein_terminus) {
            (false, _) => Terminus::No,
            (true, false) => Terminus::Peptide,
            (true, true) => Terminus::Protein,
        }
    }
}

#[derive(Debug, Clone, Default)]
struct FormList {
    variants: Vec<AminoAcid>,
    slots: Vec<u8>,
    label_slots: Vec<u8>,
    labels: Vec<u8>,
}

pub mod label_code {
    pub const LIGHT: u8 = 0;
    pub const MAX_CHANNEL: u8 = 250;
    pub const NEUTRAL: u8 = 254;
    pub const MIXED: u8 = 255;

    #[derive(Debug, Clone, Copy, PartialEq, Eq, Default)]
    pub struct LabelState(u8);

    impl LabelState {
        #[inline]
        pub fn step(self, code: u8) -> Option<Self> {
            match code {
                NEUTRAL => Some(self),
                MIXED => None,
                c => {
                    let target = c + 1;
                    if self.0 == 0 || self.0 == target {
                        Some(LabelState(target))
                    } else {
                        None
                    }
                }
            }
        }

        pub fn channel(self) -> Option<u8> {
            self.0.checked_sub(1)
        }
    }
}

#[inline]
pub fn form_index(residue: u8, n: Terminus, c: Terminus) -> usize {
    residue as usize * 9 + n.rank() * 3 + c.rank()
}

pub const FORM_TABLE_LEN: usize = 256 * 9;

#[derive(Debug, Clone)]
pub struct AminoAcidSet {
    forms: Vec<FormList>,
    declared: Vec<&'static Modification>,
    stack_parts: FxHashMap<usize, Vec<&'static Modification>>,
    aa_lists_cache: FxHashMap<ModLocation, Vec<AminoAcid>>,
    has_cterm_mods: bool,
    min_aa_mass: f64,
    max_aa_mass: f64,
    max_residue_mod_mass: f64,
    max_fixed_term_mod_mass: f64,
    peptide_cleavage_credit: i32,
    peptide_cleavage_penalty: i32,
    neighboring_aa_cleavage_credit: i32,
    neighboring_aa_cleavage_penalty: i32,
    has_label_channels: bool,
    label_channels: Vec<(String, u8)>,
}

impl AminoAcidSet {
    pub fn variants_for(&self, residue: u8, location: ModLocation) -> &[AminoAcid] {
        let (n, c) = match location {
            ModLocation::Anywhere => (Terminus::No, Terminus::No),
            ModLocation::NTerm => (Terminus::Peptide, Terminus::No),
            ModLocation::ProtNTerm => (Terminus::Protein, Terminus::No),
            ModLocation::CTerm => (Terminus::No, Terminus::Peptide),
            ModLocation::ProtCTerm => (Terminus::No, Terminus::Protein),
        };
        self.forms(residue, n, c)
    }

    #[inline]
    pub fn forms(&self, residue: u8, n: Terminus, c: Terminus) -> &[AminoAcid] {
        &self.forms[form_index(residue, n, c)].variants
    }

    #[inline]
    pub fn form_slots(&self, residue: u8, n: Terminus, c: Terminus) -> &[u8] {
        &self.forms[form_index(residue, n, c)].slots
    }

    #[inline]
    pub fn form_label_slots(&self, residue: u8, n: Terminus, c: Terminus) -> &[u8] {
        &self.forms[form_index(residue, n, c)].label_slots
    }

    #[inline]
    pub fn form_labels(&self, residue: u8, n: Terminus, c: Terminus) -> Option<&[u8]> {
        let l = &self.forms[form_index(residue, n, c)].labels;
        (!l.is_empty()).then_some(l.as_slice())
    }

    pub fn label_channels(&self) -> &[(String, u8)] {
        &self.label_channels
    }

    pub fn has_label_channels(&self) -> bool {
        self.has_label_channels
    }

    pub fn form_lists(&self) -> impl Iterator<Item = (&[AminoAcid], &[u8])> {
        self.forms
            .iter()
            .filter(|f| !f.variants.is_empty())
            .map(|f| (f.variants.as_slice(), f.slots.as_slice()))
    }

    pub fn form_lists_with_labels(
        &self,
    ) -> impl Iterator<Item = (&[AminoAcid], &[u8], &[u8], &[u8])> {
        self.forms
            .iter()
            .filter(|f| !f.variants.is_empty())
            .map(|f| {
                (
                    f.variants.as_slice(),
                    f.slots.as_slice(),
                    f.label_slots.as_slice(),
                    f.labels.as_slice(),
                )
            })
    }

    pub fn standard(&self, residue: u8) -> Option<&AminoAcid> {
        self.variants_for(residue, ModLocation::Anywhere)
            .iter()
            .find(|aa| !aa.is_modified())
    }

    pub fn contains_cterm_mods(&self) -> bool {
        self.has_cterm_mods
    }
    pub fn min_aa_mass(&self) -> f64 {
        self.min_aa_mass
    }
    pub fn max_aa_mass(&self) -> f64 {
        self.max_aa_mass
    }
    pub fn max_residue_mod_mass(&self) -> f64 {
        self.max_residue_mod_mass
    }
    pub fn max_fixed_term_mod_mass(&self) -> f64 {
        self.max_fixed_term_mod_mass
    }

    pub fn iter_variants(&self) -> impl Iterator<Item = &AminoAcid> {
        self.forms.iter().flat_map(|f| f.variants.iter())
    }

    pub fn fixed_terminal_deltas(&self) -> Vec<(ModLocation, f64)> {
        let mut out: Vec<(ModLocation, f64)> = Vec::new();
        for m in self.declared.iter().filter(|m| m.fixed && m.location != ModLocation::Anywhere) {
            if !out.iter().any(|(l, d)| *l == m.location && d.to_bits() == m.mass_delta.to_bits()) {
                out.push((m.location, m.mass_delta));
            }
        }
        out
    }

    pub fn fixed_anywhere_deltas(&self) -> Vec<(u8, f64)> {
        let mut out: Vec<(u8, f64)> = Vec::new();
        for m in self.declared.iter().filter(|m| m.fixed && m.location == ModLocation::Anywhere) {
            for &r in STANDARD_RESIDUES {
                if m.applies_to(r, m.location) && !out.iter().any(|(x, d)| *x == r && d.to_bits() == m.mass_delta.to_bits()) {
                    out.push((r, m.mass_delta));
                }
            }
        }
        out
    }

    pub fn fixed_mod_deltas(&self) -> Vec<(u8, f64)> {
        let mut seen: std::collections::HashSet<(u8, u64)> = std::collections::HashSet::new();
        let mut out: Vec<(u8, f64)> = Vec::new();
        for m in self.declared.iter().filter(|m| m.fixed) {
            for &r in STANDARD_RESIDUES {
                if m.applies_to(r, m.location) && seen.insert((r, m.mass_delta.to_bits())) {
                    out.push((r, m.mass_delta));
                }
            }
        }
        out
    }

    pub fn distinct_mods(&self) -> Vec<&'static Modification> {
        self.declared.clone()
    }

    pub fn stack_parts(&self, m: &Modification) -> Option<&[&'static Modification]> {
        self.stack_parts
            .get(&(std::ptr::from_ref(m) as usize))
            .map(Vec::as_slice)
    }

    pub fn fixed_mod_fingerprint(&self) -> Vec<(u8, u8, u64)> {
        fn loc_rank(l: ModLocation) -> u8 {
            match l {
                ModLocation::Anywhere => 0,
                ModLocation::NTerm => 1,
                ModLocation::CTerm => 2,
                ModLocation::ProtNTerm => 3,
                ModLocation::ProtCTerm => 4,
            }
        }
        let mut out: Vec<(u8, u8, u64)> = Vec::new();
        for m in self.declared.iter().filter(|m| m.fixed) {
            for &r in STANDARD_RESIDUES {
                if m.applies_to(r, m.location) {
                    out.push((r, loc_rank(m.location), m.mass_delta.to_bits()));
                }
            }
        }
        out.sort_unstable();
        out.dedup();
        out
    }

    pub fn aa_list_for(&self, location: ModLocation) -> Vec<&AminoAcid> {
        self.aa_lists_cache
            .get(&location)
            .map(|v| v.iter().collect())
            .unwrap_or_default()
    }

    pub fn cached_aa_list(&self, location: ModLocation) -> &[AminoAcid] {
        self.aa_lists_cache
            .get(&location)
            .map(|v| v.as_slice())
            .unwrap_or(&[])
    }

    pub fn peptide_cleavage_credit(&self) -> i32 {
        self.peptide_cleavage_credit
    }

    pub fn peptide_cleavage_penalty(&self) -> i32 {
        self.peptide_cleavage_penalty
    }

    pub fn neighboring_aa_cleavage_credit(&self) -> i32 {
        self.neighboring_aa_cleavage_credit
    }

    pub fn neighboring_aa_cleavage_penalty(&self) -> i32 {
        self.neighboring_aa_cleavage_penalty
    }

    pub fn prob_cleavage_sites(&self, enzyme: Enzyme) -> f32 {
        let residues = enzyme.residues();
        if residues.is_empty() {
            return 0.0;
        }
        let prob_per_aa = 1.0_f32 / STANDARD_RESIDUES.len() as f32;
        residues
            .iter()
            .filter(|&&r| STANDARD_RESIDUES.contains(&r))
            .count() as f32
            * prob_per_aa
    }

    pub fn register_enzyme(
        &mut self,
        enzyme: Enzyme,
        peptide_efficiency: f32,
        neighboring_efficiency: f32,
    ) {
        let prob = self.prob_cleavage_sites(enzyme);
        if prob <= 0.0 || prob >= 1.0 || peptide_efficiency == 0.0 {
            return;
        }
        let credit = |eff: f32| -> i32 { ((eff as f64 / prob as f64).ln()).round() as i32 };
        let penalty =
            |eff: f32| -> i32 { (((1.0 - eff) as f64 / (1.0 - prob) as f64).ln()).round() as i32 };
        self.peptide_cleavage_credit = credit(peptide_efficiency);
        self.peptide_cleavage_penalty = penalty(peptide_efficiency);
        self.neighboring_aa_cleavage_credit = credit(neighboring_efficiency);
        self.neighboring_aa_cleavage_penalty = penalty(neighboring_efficiency);
    }
}

#[derive(Debug, Clone)]
pub struct AminoAcidSetBuilder {
    fixed_mods: Vec<Modification>,
    variable_mods: Vec<Modification>,
    label_channels: Vec<(String, u8)>,
}

impl AminoAcidSetBuilder {
    pub fn new_standard() -> Self {
        Self {
            fixed_mods: vec![],
            variable_mods: vec![],
            label_channels: vec![],
        }
    }

    pub fn new_standard_with_carbamidomethyl_c() -> Self {
        let cam = Modification {
            name: "Carbamidomethyl".to_string(),
            mass_delta: 57.02146,
            residue: ResidueSpec::Specific(b'C'),
            location: ModLocation::Anywhere,
            fixed: true,
            accession: Some("UNIMOD:4".to_string()),
            neutral_losses: Vec::new(),
            loss_class: 0,
            annotated_site_mass: None,
            blocks_cleavage: false,
            parts: Vec::new(),
            diagnostic_ions: Vec::new(),
        };
        Self {
            fixed_mods: vec![cam],
            variable_mods: vec![],
            label_channels: vec![],
        }
    }

    pub fn with_label_channels(mut self, labels: Vec<(String, u8)>) -> Self {
        self.label_channels = labels;
        self
    }

    pub fn add_fixed_mod(mut self, m: Modification) -> Self {
        self.fixed_mods.push(m);
        self
    }

    pub fn add_variable_mod(mut self, m: Modification) -> Self {
        self.variable_mods.push(m);
        self
    }

    pub fn add_mods(mut self, mods: impl IntoIterator<Item = Modification>) -> Self {
        for m in mods {
            if m.fixed {
                self.fixed_mods.push(m);
            } else {
                self.variable_mods.push(m);
            }
        }
        self
    }

    pub fn build(self) -> Result<AminoAcidSet, AaSetError> {
        for m in self.fixed_mods.iter().chain(self.variable_mods.iter()) {
            if m.mass_delta.abs() > IMPLAUSIBLE_MASS_THRESHOLD {
                return Err(AaSetError::ImplausibleMassDelta {
                    name: m.name.clone(),
                    delta: m.mass_delta,
                });
            }
        }

        for fm in &self.fixed_mods {
            for vm in &self.variable_mods {
                if mods_target_same_slot(fm, vm) {
                    let res_char = match fm.residue {
                        ResidueSpec::Specific(r) => r as char,
                        ResidueSpec::Wildcard => '*',
                    };
                    return Err(AaSetError::ConflictingMods {
                        residue: res_char,
                        location: fm.location,
                    });
                }
            }
        }

        let fixed_mods: Vec<&'static Modification> = self
            .fixed_mods
            .iter()
            .cloned()
            .map(crate::modification::leak_mod)
            .collect();
        let variable_mods: Vec<&'static Modification> = self
            .variable_mods
            .iter()
            .cloned()
            .map(crate::modification::leak_mod)
            .collect();

        let label_channels = self.label_channels.clone();
        let label_of = |m: &Modification| -> Option<u8> {
            label_channels
                .iter()
                .find(|(name, _)| *name == m.name)
                .map(|(_, c)| *c)
        };
        let mut stacked: FxHashMap<(u8, Vec<usize>), &'static Modification> = FxHashMap::default();
        let mut stack_parts: FxHashMap<usize, Vec<&'static Modification>> = FxHashMap::default();
        let mut forms = vec![FormList::default(); FORM_TABLE_LEN];
        for &r in STANDARD_RESIDUES {
            let std_aa = AminoAcid::standard(r).expect("STANDARD_RESIDUES has only valid residues");
            let side = site_choices(r, &fixed_mods, &variable_mods, &[ModLocation::Anywhere]);
            for n in Terminus::ALL {
                let n_choices = terminal_choices(
                    r,
                    n,
                    ModLocation::NTerm,
                    ModLocation::ProtNTerm,
                    &fixed_mods,
                    &variable_mods,
                );
                for c in Terminus::ALL {
                    let c_choices = terminal_choices(
                        r,
                        c,
                        ModLocation::CTerm,
                        ModLocation::ProtCTerm,
                        &fixed_mods,
                        &variable_mods,
                    );
                    let list = &mut forms[form_index(r, n, c)];
                    let labelable = |choices: &[(Vec<&'static Modification>, u8)]| {
                        choices
                            .iter()
                            .any(|(mods, _)| mods.iter().any(|m| label_of(m).is_some()))
                    };
                    let lab_sites = [
                        labelable(&side),
                        labelable(&n_choices),
                        labelable(&c_choices),
                    ];
                    let position_labelable = lab_sites.iter().any(|&b| b);
                    for (n_mods, n_slots) in &n_choices {
                        for (c_mods, c_slots) in &c_choices {
                            for (s_mods, s_slots) in &side {
                                let code = form_label_code(
                                    [
                                        (s_mods.as_slice(), *s_slots, lab_sites[0]),
                                        (n_mods.as_slice(), *n_slots, lab_sites[1]),
                                        (c_mods.as_slice(), *c_slots, lab_sites[2]),
                                    ],
                                    &label_of,
                                );
                                let parts: Vec<&'static Modification> =
                                    s_mods.iter().chain(n_mods).chain(c_mods).copied().collect();
                                if parts
                                    .iter()
                                    .filter(|m| m.annotated_site_mass.is_some())
                                    .count()
                                    > 1
                                {
                                    continue;
                                }
                                let slots = s_slots + n_slots + c_slots;
                                let label_slots =
                                    parts.iter().filter(|m| label_of(m).is_some()).count() as u8;
                                let aa = match parts.as_slice() {
                                    [] => std_aa.clone(),
                                    [m] => std_aa.clone().with_mod(*m),
                                    _ => {
                                        let key = (
                                            r,
                                            parts
                                                .iter()
                                                .map(|m| std::ptr::from_ref(*m) as usize)
                                                .collect(),
                                        );
                                        let m = *stacked.entry(key).or_insert_with(|| {
                                            crate::modification::leak_mod(stack_mods(
                                                r, &parts, slots,
                                            ))
                                        });
                                        stack_parts
                                            .entry(std::ptr::from_ref(m) as usize)
                                            .or_insert_with(|| parts.clone());
                                        std_aa.clone().with_mod(m)
                                    }
                                };
                                let restricted = |a: &AminoAcid| {
                                    a.mod_.is_some_and(|m| m.annotated_site_mass.is_some())
                                };
                                match list.variants.iter().position(|u| *u == aa) {
                                    None => {
                                        list.variants.push(aa);
                                        list.slots.push(slots);
                                        list.label_slots.push(label_slots);
                                        if position_labelable {
                                            list.labels.push(code);
                                        }
                                    }
                                    Some(i)
                                        if restricted(&list.variants[i]) && !restricted(&aa) =>
                                    {
                                        list.variants[i] = aa;
                                        list.slots[i] = slots;
                                        list.label_slots[i] = label_slots;
                                        if position_labelable {
                                            list.labels[i] = code;
                                        }
                                    }
                                    Some(_) => {}
                                }
                            }
                        }
                    }
                }
            }
        }

        let standard_masses: Vec<f64> = STANDARD_RESIDUES
            .iter()
            .filter_map(|&r| AminoAcid::standard(r).map(|aa| aa.mass))
            .collect();
        let min_aa_mass = standard_masses
            .iter()
            .copied()
            .fold(f64::INFINITY, f64::min);
        let max_aa_mass = standard_masses
            .iter()
            .copied()
            .fold(f64::NEG_INFINITY, f64::max);

        let mut max_mod_delta = 0.0_f64;
        for m in self.fixed_mods.iter().chain(self.variable_mods.iter()) {
            if m.mass_delta > max_mod_delta {
                max_mod_delta = m.mass_delta;
            }
        }
        let max_residue_mod_mass = max_aa_mass + max_mod_delta;

        let max_fixed_term_mod_mass = self
            .fixed_mods
            .iter()
            .filter(|m| {
                matches!(
                    m.location,
                    ModLocation::NTerm
                        | ModLocation::CTerm
                        | ModLocation::ProtNTerm
                        | ModLocation::ProtCTerm
                )
            })
            .map(|m| m.mass_delta)
            .fold(0.0_f64, f64::max);

        let has_cterm_mods = self
            .fixed_mods
            .iter()
            .chain(self.variable_mods.iter())
            .any(|m| matches!(m.location, ModLocation::CTerm | ModLocation::ProtCTerm));

        let mut aa_lists_cache: FxHashMap<ModLocation, Vec<AminoAcid>> = FxHashMap::default();
        for (loc, n, c) in [
            (ModLocation::Anywhere, Terminus::No, Terminus::No),
            (ModLocation::NTerm, Terminus::Peptide, Terminus::No),
            (ModLocation::CTerm, Terminus::No, Terminus::Peptide),
            (ModLocation::ProtNTerm, Terminus::Protein, Terminus::No),
            (ModLocation::ProtCTerm, Terminus::No, Terminus::Protein),
        ] {
            let list: Vec<AminoAcid> = STANDARD_RESIDUES
                .iter()
                .flat_map(|&r| forms[form_index(r, n, c)].variants.iter().cloned())
                .collect();
            aa_lists_cache.insert(loc, list);
        }

        Ok(AminoAcidSet {
            forms,
            declared: fixed_mods.into_iter().chain(variable_mods).collect(),
            stack_parts,
            aa_lists_cache,
            has_cterm_mods,
            min_aa_mass,
            max_aa_mass,
            max_residue_mod_mass,
            max_fixed_term_mod_mass,
            peptide_cleavage_credit: 0,
            peptide_cleavage_penalty: 0,
            neighboring_aa_cleavage_credit: 0,
            neighboring_aa_cleavage_penalty: 0,
            has_label_channels: !self.label_channels.is_empty(),
            label_channels: self.label_channels,
        })
    }
}

fn form_label_code(
    sites: [(&[&'static Modification], u8, bool); 3],
    label_of: &impl Fn(&Modification) -> Option<u8>,
) -> u8 {
    let mut channel: Option<u8> = None;
    let mut light = false;
    for (mods, slots, labelable) in sites {
        if !labelable {
            continue;
        }
        match mods.iter().find_map(|m| label_of(m)) {
            Some(c) => match channel {
                Some(prev) if prev != c => return label_code::MIXED,
                _ => channel = Some(c),
            },
            None if slots > 0 || !mods.is_empty() => {}
            None => light = true,
        }
    }
    match (channel, light) {
        (Some(_), true) => label_code::MIXED,
        (Some(c), false) => c.min(label_code::MAX_CHANNEL),
        (None, true) => label_code::LIGHT,
        (None, false) => label_code::NEUTRAL,
    }
}

fn site_choices(
    r: u8,
    fixed_mods: &[&'static Modification],
    variable_mods: &[&'static Modification],
    locations: &[ModLocation],
) -> Vec<(Vec<&'static Modification>, u8)> {
    let at =
        |m: &&&'static Modification| locations.contains(&m.location) && m.applies_to(r, m.location);
    let mut fixed: Vec<&'static Modification> = fixed_mods.iter().filter(at).copied().collect();
    fixed.sort_by(|a, b| {
        a.name
            .cmp(&b.name)
            .then(a.mass_delta.to_bits().cmp(&b.mass_delta.to_bits()))
    });
    let mut out = vec![(fixed, 0u8)];
    out.extend(variable_mods.iter().filter(at).map(|m| (vec![*m], 1u8)));
    out
}

fn terminal_choices(
    r: u8,
    t: Terminus,
    peptide_loc: ModLocation,
    protein_loc: ModLocation,
    fixed_mods: &[&'static Modification],
    variable_mods: &[&'static Modification],
) -> Vec<(Vec<&'static Modification>, u8)> {
    match t {
        Terminus::No => vec![(Vec::new(), 0)],
        Terminus::Peptide => site_choices(r, fixed_mods, variable_mods, &[peptide_loc]),
        Terminus::Protein => {
            let protein_fixed = fixed_mods
                .iter()
                .any(|m| m.location == protein_loc && m.applies_to(r, protein_loc));
            let fixed_loc = if protein_fixed {
                protein_loc
            } else {
                peptide_loc
            };
            let mut out = site_choices(r, fixed_mods, &[], &[fixed_loc]);
            out.extend(
                site_choices(r, &[], variable_mods, &[peptide_loc, protein_loc])
                    .into_iter()
                    .skip(1),
            );
            out
        }
    }
}

fn stack_mods(r: u8, parts: &[&'static Modification], slots: u8) -> Modification {
    let mut neutral_losses: Vec<f64> = Vec::new();
    for loss in parts.iter().flat_map(|m| m.neutral_losses.iter()) {
        if !neutral_losses.iter().any(|l| l.to_bits() == loss.to_bits()) {
            neutral_losses.push(*loss);
        }
    }
    Modification {
        name: parts
            .iter()
            .map(|m| m.name.as_str())
            .collect::<Vec<_>>()
            .join("+"),
        mass_delta: parts.iter().map(|m| m.mass_delta).sum(),
        residue: ResidueSpec::Specific(r),
        location: parts
            .iter()
            .map(|m| m.location)
            .find(|l| *l != ModLocation::Anywhere)
            .unwrap_or(ModLocation::Anywhere),
        fixed: slots == 0,
        accession: None,
        neutral_losses,
        loss_class: parts
            .iter()
            .map(|m| m.loss_class)
            .find(|c| *c != 0)
            .unwrap_or(0),
        annotated_site_mass: parts.iter().find_map(|m| m.annotated_site_mass),
        blocks_cleavage: parts.iter().any(|m| m.blocks_cleavage),
        parts: parts.to_vec(),
        diagnostic_ions: Vec::new(),
    }
}

fn mods_target_same_slot(a: &Modification, b: &Modification) -> bool {
    a.residue == b.residue && a.location == b.location
}

#[derive(thiserror::Error, Debug)]
pub enum AaSetError {
    #[error("conflicting fixed and variable mod for residue {residue:?} at {location:?}")]
    ConflictingMods {
        residue: char,
        location: ModLocation,
    },
    #[error("mod {name:?} mass delta {delta} is implausible (>5000 Da)")]
    ImplausibleMassDelta { name: String, delta: f64 },
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::amino_acid::AminoAcid;
    use crate::enzyme::Enzyme;
    use crate::modification::{ModLocation, Modification, ResidueSpec};

    fn spec(list: &str, fixed: bool) -> Vec<Modification> {
        crate::modification::spec::parse_mod_list(list, fixed).unwrap()
    }

    fn annotated(item: &str) -> Modification {
        let mut m = spec(item, false).remove(0);
        m.annotated_site_mass = Some(m.mass_delta);
        m
    }

    fn carbamidomethyl_c() -> Modification {
        Modification {
            name: "Carbamidomethyl".to_string(),
            mass_delta: 57.02146,
            residue: ResidueSpec::Specific(b'C'),
            location: ModLocation::Anywhere,
            fixed: true,
            accession: None,
            neutral_losses: Vec::new(),
            loss_class: 0,
            annotated_site_mass: None,
            blocks_cleavage: false,
            parts: Vec::new(),
            diagnostic_ions: Vec::new(),
        }
    }

    fn oxidation_m() -> Modification {
        Modification {
            name: "Oxidation".to_string(),
            mass_delta: 15.99491,
            residue: ResidueSpec::Specific(b'M'),
            location: ModLocation::Anywhere,
            fixed: false,
            accession: None,
            neutral_losses: Vec::new(),
            loss_class: 0,
            annotated_site_mass: None,
            blocks_cleavage: false,
            parts: Vec::new(),
            diagnostic_ions: Vec::new(),
        }
    }

    #[test]
    fn standard_set_has_20_residues() {
        let set = AminoAcidSetBuilder::new_standard().build().unwrap();
        let mut seen = std::collections::HashSet::new();
        for aa in set.iter_variants() {
            seen.insert(aa.residue);
        }
        assert_eq!(seen.len(), 20);
    }

    #[test]
    fn standard_set_no_mods() {
        let set = AminoAcidSetBuilder::new_standard().build().unwrap();
        for aa in set.iter_variants() {
            assert!(!aa.is_modified());
        }
    }

    #[test]
    fn fixed_mod_replaces_residue() {
        let set = AminoAcidSetBuilder::new_standard()
            .add_fixed_mod(carbamidomethyl_c())
            .build()
            .unwrap();
        let c_variants = set.variants_for(b'C', ModLocation::Anywhere);
        assert_eq!(c_variants.len(), 1);
        assert!(c_variants[0].is_modified());
    }

    #[test]
    fn variable_mod_adds_residue_variant() {
        let set = AminoAcidSetBuilder::new_standard()
            .add_variable_mod(oxidation_m())
            .build()
            .unwrap();
        let m_variants = set.variants_for(b'M', ModLocation::Anywhere);
        assert_eq!(m_variants.len(), 2);
        assert!(m_variants.iter().any(|aa| !aa.is_modified()));
        assert!(m_variants.iter().any(|aa| aa.is_modified()));
    }

    #[test]
    fn conflicting_fixed_and_variable_errors() {
        let cam_fixed = carbamidomethyl_c();
        let mut cam_variable = carbamidomethyl_c();
        cam_variable.fixed = false;

        let err = AminoAcidSetBuilder::new_standard()
            .add_fixed_mod(cam_fixed)
            .add_variable_mod(cam_variable)
            .build()
            .unwrap_err();
        assert!(matches!(
            err,
            AaSetError::ConflictingMods {
                residue: 'C',
                location: ModLocation::Anywhere
            }
        ));
    }

    #[test]
    fn implausible_mass_errors() {
        let bad = Modification {
            name: "Bad".to_string(),
            mass_delta: 8000.0,
            residue: ResidueSpec::Specific(b'C'),
            location: ModLocation::Anywhere,
            fixed: true,
            accession: None,
            neutral_losses: Vec::new(),
            loss_class: 0,
            annotated_site_mass: None,
            blocks_cleavage: false,
            parts: Vec::new(),
            diagnostic_ions: Vec::new(),
        };
        let err = AminoAcidSetBuilder::new_standard()
            .add_fixed_mod(bad)
            .build()
            .unwrap_err();
        assert!(matches!(err, AaSetError::ImplausibleMassDelta { .. }));
    }

    #[test]
    fn standard_lookup() {
        let set = AminoAcidSetBuilder::new_standard().build().unwrap();
        let g = set.standard(b'G').unwrap();
        assert_eq!(g.residue, b'G');
        assert!(set.standard(b'!').is_none());
    }

    #[test]
    fn min_max_aa_mass() {
        let set = AminoAcidSetBuilder::new_standard().build().unwrap();
        let g = AminoAcid::standard(b'G').unwrap().mass;
        let w = AminoAcid::standard(b'W').unwrap().mass;
        assert_eq!(set.min_aa_mass(), g);
        assert_eq!(set.max_aa_mass(), w);
    }

    #[test]
    fn max_residue_mod_mass_includes_mods() {
        let set = AminoAcidSetBuilder::new_standard()
            .add_variable_mod(oxidation_m())
            .build()
            .unwrap();
        let w = AminoAcid::standard(b'W').unwrap().mass;
        let expected = w + 15.99491;
        assert!((set.max_residue_mod_mass() - expected).abs() < 1e-9);
    }

    #[test]
    fn contains_cterm_mods_default_false() {
        let set = AminoAcidSetBuilder::new_standard().build().unwrap();
        assert!(!set.contains_cterm_mods());
    }

    #[test]
    fn contains_cterm_mods_when_added() {
        let cterm_mod = Modification {
            name: "Amide".to_string(),
            mass_delta: -0.984016,
            residue: ResidueSpec::Wildcard,
            location: ModLocation::CTerm,
            fixed: false,
            accession: None,
            neutral_losses: Vec::new(),
            loss_class: 0,
            annotated_site_mass: None,
            blocks_cleavage: false,
            parts: Vec::new(),
            diagnostic_ions: Vec::new(),
        };
        let set = AminoAcidSetBuilder::new_standard()
            .add_variable_mod(cterm_mod)
            .build()
            .unwrap();
        assert!(set.contains_cterm_mods());
    }

    #[test]
    fn standard_with_carbamidomethyl_c_convenience() {
        let set = AminoAcidSetBuilder::new_standard_with_carbamidomethyl_c()
            .build()
            .unwrap();
        let c_variants = set.variants_for(b'C', ModLocation::Anywhere);
        assert_eq!(c_variants.len(), 1);
        assert!(c_variants[0].is_modified());
    }

    #[test]
    fn fixed_mod_deltas_empty_for_standard_set() {
        let set = AminoAcidSetBuilder::new_standard().build().unwrap();
        assert!(set.fixed_mod_deltas().is_empty());
    }

    #[test]
    fn fixed_mod_deltas_reports_cam_c_once() {
        let set = AminoAcidSetBuilder::new_standard_with_carbamidomethyl_c()
            .build()
            .unwrap();
        let deltas = set.fixed_mod_deltas();
        assert_eq!(
            deltas.len(),
            1,
            "exactly one distinct fixed mod, got {deltas:?}"
        );
        let (res, delta) = deltas[0];
        assert_eq!(res, b'C');
        assert!(
            (delta - 57.02146).abs() < 1e-6,
            "expected CAM-C delta, got {delta}"
        );
    }

    #[test]
    fn fixed_mod_deltas_reports_each_distinct_fixed_mod() {
        let tmt_k = Modification {
            name: "TMT6plex".to_string(),
            mass_delta: 229.162932,
            residue: ResidueSpec::Specific(b'K'),
            location: ModLocation::Anywhere,
            fixed: true,
            accession: None,
            neutral_losses: Vec::new(),
            loss_class: 0,
            annotated_site_mass: None,
            blocks_cleavage: false,
            parts: Vec::new(),
            diagnostic_ions: Vec::new(),
        };
        let set = AminoAcidSetBuilder::new_standard()
            .add_fixed_mod(carbamidomethyl_c())
            .add_fixed_mod(tmt_k)
            .add_variable_mod(oxidation_m())
            .build()
            .unwrap();
        let mut deltas = set.fixed_mod_deltas();
        deltas.sort_by(|a, b| a.0.cmp(&b.0));
        assert_eq!(deltas.len(), 2, "two distinct fixed mods, got {deltas:?}");
        assert_eq!(deltas[0].0, b'C');
        assert!((deltas[0].1 - 57.02146).abs() < 1e-6);
        assert_eq!(deltas[1].0, b'K');
        assert!((deltas[1].1 - 229.162932).abs() < 1e-6);
    }

    #[test]
    fn fixed_mod_fingerprint_distinguishes_location_where_deltas_collide() {
        let mk = |loc: ModLocation| Modification {
            name: "TMT6plex".to_string(),
            mass_delta: 229.162932,
            residue: ResidueSpec::Specific(b'K'),
            location: loc,
            fixed: true,
            accession: None,
            neutral_losses: Vec::new(),
            loss_class: 0,
            annotated_site_mass: None,
            blocks_cleavage: false,
            parts: Vec::new(),
            diagnostic_ions: Vec::new(),
        };
        let anywhere = AminoAcidSetBuilder::new_standard()
            .add_fixed_mod(mk(ModLocation::Anywhere))
            .build()
            .unwrap();
        let nterm = AminoAcidSetBuilder::new_standard()
            .add_fixed_mod(mk(ModLocation::NTerm))
            .build()
            .unwrap();
        assert_eq!(
            anywhere.fixed_mod_deltas(),
            nterm.fixed_mod_deltas(),
            "fixed_mod_deltas is expected to collide on location"
        );
        assert_ne!(
            anywhere.fixed_mod_fingerprint(),
            nterm.fixed_mod_fingerprint(),
            "fingerprint must distinguish Anywhere vs NTerm fixed mod"
        );
    }

    #[test]
    fn fixed_mods_split_into_anywhere_and_terminal() {
        let mk = |r: ResidueSpec, loc: ModLocation, d: f64| Modification {
            name: format!("m{d}"),
            mass_delta: d,
            residue: r,
            location: loc,
            fixed: true,
            accession: Some(format!("UNIMOD:{}", d as u32)),
            neutral_losses: Vec::new(),
            loss_class: 0,
            annotated_site_mass: None,
            blocks_cleavage: false,
            parts: Vec::new(),
            diagnostic_ions: Vec::new(),
        };
        let set = AminoAcidSetBuilder::new_standard()
            .add_fixed_mod(mk(ResidueSpec::Specific(b'K'), ModLocation::Anywhere, 229.0))
            .add_fixed_mod(mk(ResidueSpec::Wildcard, ModLocation::NTerm, 229.0))
            .build()
            .unwrap();
        assert_eq!(set.fixed_anywhere_deltas(), vec![(b'K', 229.0)]);
        assert_eq!(set.fixed_terminal_deltas(), vec![(ModLocation::NTerm, 229.0)]);
        let stacked = set
            .variants_for(b'K', ModLocation::NTerm)
            .iter()
            .filter_map(|aa| aa.mod_)
            .find(|m| m.parts.len() == 2)
            .expect("a stacked N-terminal K");
        assert_eq!(stacked.accession, None);
        assert!(stacked.components().iter().all(|m| m.accession.is_some()));
    }

    #[test]
    fn add_mods_sorts_fixed_and_variable() {
        let set = AminoAcidSetBuilder::new_standard()
            .add_mods(spec("Carbamidomethyl (C)", true))
            .add_mods(spec("Oxidation (M)", false))
            .build()
            .unwrap();

        assert_eq!(set.variants_for(b'C', ModLocation::Anywhere).len(), 1);
        assert!(set.variants_for(b'C', ModLocation::Anywhere)[0].is_modified());
        assert_eq!(set.variants_for(b'M', ModLocation::Anywhere).len(), 2);
    }

    #[test]
    fn peptide_cleavage_credit_default_is_zero() {
        let set = AminoAcidSetBuilder::new_standard().build().unwrap();
        assert_eq!(set.peptide_cleavage_credit(), 0);
    }

    #[test]
    fn prob_cleavage_sites_for_trypsin_is_approximately_0_1() {
        let set = AminoAcidSetBuilder::new_standard().build().unwrap();
        let prob = set.prob_cleavage_sites(Enzyme::Trypsin);
        assert!((prob - 0.1_f32).abs() < 1e-5, "expected ~0.1, got {prob}");
    }

    #[test]
    fn aa_list_for_anywhere_returns_20_residues() {
        let set = AminoAcidSetBuilder::new_standard().build().unwrap();
        let list = set.aa_list_for(ModLocation::Anywhere);
        assert_eq!(
            list.len(),
            20,
            "standard set should have exactly 20 standard residues"
        );
        for &r in STANDARD_RESIDUES {
            assert!(
                list.iter().any(|aa| aa.residue == r),
                "residue {} missing from aa_list_for(Anywhere)",
                r as char
            );
        }
    }

    #[test]
    fn aa_list_for_nterm_returns_at_least_20_residues() {
        let set = AminoAcidSetBuilder::new_standard().build().unwrap();
        let list = set.aa_list_for(ModLocation::NTerm);
        assert_eq!(list.len(), 20);
    }

    #[test]
    fn aa_list_for_nterm_includes_terminal_mods() {
        let nterm_mod = Modification {
            name: "TMT6plex".to_string(),
            mass_delta: 229.16293,
            residue: ResidueSpec::Wildcard,
            location: ModLocation::NTerm,
            fixed: false,
            accession: None,
            neutral_losses: Vec::new(),
            loss_class: 0,
            annotated_site_mass: None,
            blocks_cleavage: false,
            parts: Vec::new(),
            diagnostic_ions: Vec::new(),
        };
        let set = AminoAcidSetBuilder::new_standard()
            .add_variable_mod(nterm_mod)
            .build()
            .unwrap();
        let list = set.aa_list_for(ModLocation::NTerm);
        assert_eq!(
            list.len(),
            40,
            "expected 20 standard + 20 NTerm-mod variants, got {}",
            list.len()
        );
    }

    #[test]
    fn prob_cleavage_sites_for_lysc_is_0_05() {
        let set = AminoAcidSetBuilder::new_standard().build().unwrap();
        let prob = set.prob_cleavage_sites(Enzyme::LysC);
        assert!((prob - 0.05_f32).abs() < 1e-5, "expected ~0.05, got {prob}");
    }

    #[test]
    fn prob_cleavage_sites_for_nocleavage_is_zero() {
        let set = AminoAcidSetBuilder::new_standard().build().unwrap();
        let prob = set.prob_cleavage_sites(Enzyme::NoCleavage);
        assert_eq!(prob, 0.0);
    }

    #[test]
    fn register_enzyme_sets_cleavage_scores() {
        let mut set = AminoAcidSetBuilder::new_standard().build().unwrap();
        set.register_enzyme(Enzyme::Trypsin, 0.99999, 0.99999);
        assert_eq!(set.peptide_cleavage_credit(), 2);
        assert_eq!(set.peptide_cleavage_penalty(), -11);
    }

    #[test]
    fn tmt_style_mod_lists_build() {
        let set = AminoAcidSetBuilder::new_standard()
            .add_mods(spec(
                "TMT6plex (K), TMT6plex (N-term), Carbamidomethyl (C)",
                true,
            ))
            .add_mods(spec("Oxidation (M)", false))
            .build()
            .unwrap();

        let k_variants = set.variants_for(b'K', ModLocation::Anywhere);
        assert_eq!(
            k_variants.len(),
            1,
            "K should have exactly one variant (TMT-modified)"
        );
        assert!(
            k_variants[0].is_modified(),
            "K's Anywhere variant must carry the TMT mod"
        );

        let a_nterm = set.variants_for(b'A', ModLocation::NTerm);
        assert!(
            a_nterm.iter().any(|aa| aa.is_modified()),
            "A at N-term should have a TMT variant"
        );

        let c_variants = set.variants_for(b'C', ModLocation::Anywhere);
        assert_eq!(c_variants.len(), 1);
        assert!(c_variants[0].is_modified());

        let m_variants = set.variants_for(b'M', ModLocation::Anywhere);
        assert_eq!(m_variants.len(), 2);
        assert!(m_variants.iter().any(|aa| aa.is_modified()));
        assert!(m_variants.iter().any(|aa| !aa.is_modified()));
    }

    #[test]
    fn acetyl_prot_n_term_appears_in_source_aas_for_scoring() {
        let acetyl = Modification {
            name: "Acetyl".to_string(),
            mass_delta: 42.010565,
            residue: ResidueSpec::Wildcard,
            location: ModLocation::ProtNTerm,
            fixed: false,
            accession: None,
            neutral_losses: Vec::new(),
            loss_class: 0,
            annotated_site_mass: None,
            blocks_cleavage: false,
            parts: Vec::new(),
            diagnostic_ions: Vec::new(),
        };
        let set = AminoAcidSetBuilder::new_standard()
            .add_fixed_mod(carbamidomethyl_c())
            .add_variable_mod(oxidation_m())
            .add_variable_mod(acetyl)
            .build()
            .unwrap();

        let anywhere = set.cached_aa_list(ModLocation::Anywhere);
        let prot_n = set.cached_aa_list(ModLocation::ProtNTerm);

        let n_any_modified = anywhere.iter().filter(|aa| aa.is_modified()).count();
        let n_any_acetyl = anywhere
            .iter()
            .filter(|aa| aa.mod_.as_ref().is_some_and(|m| m.name == "Acetyl"))
            .count();
        assert_eq!(
            n_any_acetyl, 0,
            "Acetyl Prot-N-term must NOT appear in Anywhere AA list"
        );

        let n_pn_acetyl = prot_n
            .iter()
            .filter(|aa| aa.mod_.as_ref().is_some_and(|m| m.name.contains("Acetyl")))
            .count();
        assert_eq!(
            n_pn_acetyl, 21,
            "Prot-N-term AA list must include 21 acetyl-bearing forms (one per Anywhere form)"
        );

        let c_pn = prot_n
            .iter()
            .find(|aa| {
                aa.residue == b'C' && aa.mod_.as_ref().is_some_and(|m| m.name.contains("Acetyl"))
            })
            .expect("Cys must have an acetyl-bearing Prot-N-term variant");
        let c_delta = c_pn.mod_.as_ref().unwrap().mass_delta;
        assert!(
            (c_delta - (57.02146 + 42.010565)).abs() < 1e-4,
            "Prot-N-term Cys acetyl variant must be CAM+Acetyl (+99.032), got +{c_delta} \
             (Acetyl-alone +42.01 would drop the mandatory fixed CAM)"
        );

        assert_eq!(
            prot_n.len(),
            anywhere.len() + 21,
            "Prot-N-term list = Anywhere list + 21 acetyl-bearing forms; \
             actual Anywhere len = {}, Prot-N-term len = {}, Anywhere modified = {}",
            anywhere.len(),
            prot_n.len(),
            n_any_modified
        );
    }

    fn form_deltas(set: &AminoAcidSet, r: u8, n: Terminus, c: Terminus) -> Vec<(f64, u8)> {
        let mut v: Vec<(f64, u8)> = set
            .forms(r, n, c)
            .iter()
            .zip(set.form_slots(r, n, c))
            .map(|(aa, &s)| {
                let d = aa.mod_.as_ref().map_or(0.0, |m| m.mass_delta);
                ((d * 1e4).round() / 1e4, s)
            })
            .collect();
        v.sort_by(|a, b| a.partial_cmp(b).unwrap());
        v
    }

    #[test]
    fn terminal_forms_stack_side_chain_and_terminal_mods() {
        let set = AminoAcidSetBuilder::new_standard()
            .add_mods(spec(
                "TMT6plex (K), TMT6plex (N-term), Carbamidomethyl (C)",
                true,
            ))
            .add_mods(spec("Oxidation (M), Acetyl (Protein N-term)", false))
            .build()
            .unwrap();
        let (no, pep, prot) = (Terminus::No, Terminus::Peptide, Terminus::Protein);
        let (tmt, ac) = (229.1629, 42.0106);
        let r4 = |x: f64| (x * 1e4).round() / 1e4;
        assert_eq!(form_deltas(&set, b'K', no, no), vec![(tmt, 0)]);
        assert_eq!(
            form_deltas(&set, b'K', pep, no),
            vec![(r4(2.0 * 229.162932), 0)]
        );
        assert_eq!(
            form_deltas(&set, b'C', pep, no),
            vec![(r4(57.021464 + 229.162932), 0)]
        );
        assert_eq!(
            form_deltas(&set, b'M', pep, no),
            vec![(tmt, 0), (r4(229.162932 + 15.994915), 1)]
        );
        assert_eq!(
            form_deltas(&set, b'M', prot, no),
            vec![
                (ac, 1),
                (r4(15.994915 + 42.010565), 2),
                (tmt, 0),
                (r4(229.162932 + 15.994915), 1),
            ]
        );
        assert_eq!(
            form_deltas(&set, b'K', prot, no),
            vec![(r4(229.162932 + 42.010565), 1), (r4(2.0 * 229.162932), 0)]
        );
        assert_eq!(
            form_deltas(&set, b'K', pep, pep),
            form_deltas(&set, b'K', pep, no)
        );
        for (_, slots) in set.form_lists() {
            assert_eq!(slots[0], 0);
            assert!(slots[1..].iter().all(|&s| s > 0));
        }
    }

    #[test]
    fn stacks_carry_at_most_one_restricted_part_and_prefer_unrestricted_forms() {
        let (prot, no) = (Terminus::Protein, Terminus::No);
        let mod_named = |set: &AminoAcidSet, r: u8, name: &str| {
            set.forms(r, prot, no)
                .iter()
                .filter_map(|aa| aa.mod_)
                .filter(|m| m.name == name)
                .collect::<Vec<_>>()
        };
        let set = AminoAcidSetBuilder::new_standard()
            .add_variable_mod(annotated("Phospho (S)"))
            .add_mods(spec("Oxidation (M)", false))
            .add_variable_mod(annotated("Acetyl (Protein N-term)"))
            .build()
            .unwrap();
        assert!(mod_named(&set, b'S', "Phospho+Acetyl").is_empty());
        let ox_ac = mod_named(&set, b'M', "Oxidation+Acetyl");
        assert_eq!(ox_ac.len(), 1);
        assert_eq!(ox_ac[0].annotated_site_mass, Some(42.010565));

        let set = AminoAcidSetBuilder::new_standard()
            .add_variable_mod(annotated("Acetyl (Protein N-term)"))
            .add_mods(spec("Acetyl (N-term)", false))
            .build()
            .unwrap();
        let ac = mod_named(&set, b'A', "Acetyl");
        assert_eq!(ac.len(), 1);
        assert_eq!(ac[0].annotated_site_mass, None);
    }

    #[test]
    fn stacked_terminal_mod_keeps_the_restricted_mods_own_site_mass() {
        let acetyl = annotated("Acetyl (Protein N-term)");
        let set = AminoAcidSetBuilder::new_standard()
            .add_fixed_mod(carbamidomethyl_c())
            .add_variable_mod(acetyl)
            .build()
            .unwrap();
        let c = set
            .variants_for(b'C', ModLocation::ProtNTerm)
            .iter()
            .find_map(|aa| aa.mod_.filter(|m| m.name.contains("Acetyl")))
            .expect("Cys has a stacked acetyl variant");
        assert!((c.mass_delta - (57.02146 + 42.010565)).abs() < 1e-6);
        assert_eq!(c.annotated_site_mass, Some(42.010565));
        let a = set
            .variants_for(b'A', ModLocation::ProtNTerm)
            .iter()
            .find_map(|aa| aa.mod_)
            .expect("Ala has an acetyl variant");
        assert_eq!(a.annotated_site_mass, Some(42.010565));
    }
}
