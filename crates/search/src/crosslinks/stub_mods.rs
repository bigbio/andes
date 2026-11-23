//! The stub-form modifications a crosslink search puts on link sites.

use std::collections::HashMap;
use std::sync::Mutex;

use model::modification::{leak_mod, ModLocation, Modification, ResidueSpec};
use xl::crosslinker::{Crosslinker, LinkForm, N_FORMS};

use super::XlCandidates;

/// The modifications that put a stub form's mass on a link site, by form and the mass of the
/// modification already on that residue. `known` holds every pair the candidate table can ask
/// for, built before the search so the threads read it without a lock; `cache` takes any other
/// pair.
pub(super) struct StubMods {
    names: Vec<String>,
    known: rustc_hash::FxHashMap<(LinkForm, u64), &'static Modification>,
    cache: Mutex<HashMap<(LinkForm, u64), &'static Modification>>,
}

impl StubMods {
    pub(super) fn new(linker: &Crosslinker) -> Self {
        let mut names = vec![String::new(); N_FORMS];
        for form in linker.forms() {
            names[form.index()] = format!("XL:{}", linker.form_name(form));
        }
        StubMods {
            names,
            known: rustc_hash::FxHashMap::default(),
            cache: Mutex::new(HashMap::new()),
        }
    }

    /// The stub modifications of `linker` for every link site of `table`'s peptides.
    pub(super) fn for_table(linker: &Crosslinker, table: &XlCandidates) -> Self {
        let mut mods = StubMods::new(linker);
        let mut existing: Vec<u64> = vec![0.0_f64.to_bits()];
        for p in &table.peptides {
            let ends = std::iter::once(&p.sites).chain(p.sites_end2.as_ref());
            for site in ends.flatten() {
                if let Some(m) = p.peptide.residues.get(site.pos).and_then(|aa| aa.mod_) {
                    existing.push(m.mass_delta.to_bits());
                }
            }
        }
        existing.sort_unstable();
        existing.dedup();
        for form in linker.forms() {
            let Some(stub) = form.delta(linker, 0.0) else {
                continue;
            };
            if form == LinkForm::Intact {
                continue;
            }
            for &bits in &existing {
                let m = mods.make(form, stub, f64::from_bits(bits));
                mods.known.insert((form, bits), m);
            }
        }
        mods
    }

    pub(super) fn make(&self, form: LinkForm, stub: f64, existing: f64) -> &'static Modification {
        let mut m = Modification::new(
            self.names[form.index()].clone(),
            existing + stub,
            ResidueSpec::Wildcard,
            ModLocation::Anywhere,
            false,
        );
        m.blocks_cleavage = true;
        leak_mod(m)
    }

    pub(super) fn get(&self, form: LinkForm, stub: f64, existing: f64) -> &'static Modification {
        let key = (form, existing.to_bits());
        if let Some(m) = self.known.get(&key) {
            return m;
        }
        let mut c = self
            .cache
            .lock()
            .unwrap_or_else(std::sync::PoisonError::into_inner);
        c.entry(key)
            .or_insert_with(|| self.make(form, stub, existing))
    }
}
