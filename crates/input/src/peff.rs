//! PSI Extended FASTA Format (PEFF) databases: FASTA whose headers carry
//! `\Key=value` annotations. The variant, modification and processing keys
//! become the same claims an annotation file holds (see
//! [`crate::annotations_file`]), checked and built the same way.
//!
//! | key | value | becomes |
//! |---|---|---|
//! | `\VariantSimple` | `(pos\|aa[\|tag])…` | a one-residue substitution |
//! | `\VariantComplex` | `(start\|end\|seq[\|tag])…` | an edit (empty `seq`: a deletion) |
//! | `\ModResUnimod` | `(pos\|UNIMOD:n\|name)…` | a site, mass from the bundled UniMod |
//! | `\ModRes` | `(pos\|\|name)…` | a site, mass from UniProt's PTM vocabulary |
//! | `\ModResPsi` | `(pos\|MOD:n\|name)…` | a site when its name is in UniProt's PTM vocabulary |
//! | `\Processed` | `(start\|end[\|accession]\|name)…` | a mature chain; signal, transit and propeptides are left out |
//!
//! Positions are 1-based and inclusive. `\PName` is the description; other
//! keys are ignored. Lines starting with `#` are the file header.

use std::io::BufRead;

use model::{PriorTier, Protein, ProteinDb, SiteAnchor};

use crate::annotation::{change_label, KnockoffOptions};
use crate::annotations_file::{
    add_claims, sequence_md5, AnnotationLoadStats, AnnotationsError, Row,
};
use crate::uniprot_ptm::{clean_description, ptm_mass};

/// Counts from a PEFF load, for the run log.
#[derive(Debug, Clone, Default, PartialEq)]
pub struct PeffLoadStats {
    /// Entries read (= reference proteins).
    pub entries: usize,
    /// Sites whose modification has no known mass (left out).
    pub sites_without_mass: usize,
    /// `\Processed` spans that are not mature chains (left out).
    pub processing_not_chain: usize,
    /// The claims, as for an annotation file.
    pub claims: AnnotationLoadStats,
}

#[derive(thiserror::Error, Debug)]
pub enum PeffError {
    #[error("reading PEFF: {0}")]
    Io(#[from] std::io::Error),
    #[error("PEFF header at line {line} has no accession")]
    EmptyAccession { line: usize },
    #[error("PEFF entry {accession}: malformed {key} value {value:?}")]
    Malformed {
        accession: String,
        key: String,
        value: String,
    },
    #[error(transparent)]
    Claims(#[from] AnnotationsError),
}

/// Load a PEFF database: one reference protein per entry, then the entries its
/// header keys make (with `knockoffs`, knockoffs of single-residue variants).
pub fn load_peff<R: BufRead>(
    reader: R,
    knockoffs: Option<KnockoffOptions>,
) -> Result<(ProteinDb, PeffLoadStats), PeffError> {
    let mut stats = PeffLoadStats::default();
    let mut db = ProteinDb::default();
    let mut rows = Vec::new();
    let mut header: Option<(String, Vec<(String, String)>)> = None;
    let mut sequence: Vec<u8> = Vec::new();
    let finish = |header: Option<(String, Vec<(String, String)>)>,
                  sequence: &mut Vec<u8>,
                  db: &mut ProteinDb,
                  rows: &mut Vec<Row>,
                  stats: &mut PeffLoadStats|
     -> Result<(), PeffError> {
        let Some((accession, keys)) = header else {
            return Ok(());
        };
        let seq = std::mem::take(sequence);
        let description = keys
            .iter()
            .find(|(k, _)| k == "PName")
            .map_or_else(String::new, |(_, v)| v.clone());
        rows.extend(entry_rows(&accession, &seq, &keys, stats)?);
        db.proteins.push(Protein {
            accession,
            description,
            sequence: seq,
            annotation: None,
        });
        stats.entries += 1;
        Ok(())
    };
    for (n, line) in reader.lines().enumerate() {
        let line = line?;
        let line = line.trim_end();
        if line.starts_with('#') || line.is_empty() {
            continue;
        }
        if let Some(h) = line.strip_prefix('>') {
            finish(header.take(), &mut sequence, &mut db, &mut rows, &mut stats)?;
            let (accession, keys) = parse_header(h);
            if accession.is_empty() {
                return Err(PeffError::EmptyAccession { line: n + 1 });
            }
            header = Some((accession, keys));
        } else {
            sequence.extend(
                line.bytes()
                    .filter(|b| !b.is_ascii_whitespace())
                    .map(|b| b.to_ascii_uppercase()),
            );
        }
    }
    finish(header.take(), &mut sequence, &mut db, &mut rows, &mut stats)?;
    add_claims(rows, &mut db, knockoffs, &mut stats.claims)?;
    Ok((db, stats))
}

/// The accession (first token) and the `\Key=value` pairs of a header line
/// without its `>`.
fn parse_header(h: &str) -> (String, Vec<(String, String)>) {
    let h = h.trim();
    let (accession, rest) = h.split_once(char::is_whitespace).unwrap_or((h, ""));
    // A header that starts with a key has no accession.
    let (accession, rest) = if accession.starts_with('\\') {
        ("", h)
    } else {
        (accession, rest)
    };
    let keys = rest
        .split(" \\")
        .map(|part| part.trim().trim_start_matches('\\'))
        .filter_map(|part| part.split_once('='))
        .map(|(k, v)| (k.trim().to_string(), v.trim().to_string()))
        .collect();
    (accession.to_string(), keys)
}

/// The `(…)` groups of a value, each split on `|`.
fn groups(value: &str) -> Option<Vec<Vec<&str>>> {
    let mut out = Vec::new();
    let mut rest = value.trim();
    while !rest.is_empty() {
        let inner = rest.strip_prefix('(')?;
        let end = inner.find(')')?;
        out.push(inner[..end].split('|').map(str::trim).collect());
        rest = inner[end + 1..].trim_start();
    }
    Some(out)
}

/// The claims of one entry's header keys.
fn entry_rows(
    accession: &str,
    seq: &[u8],
    keys: &[(String, String)],
    stats: &mut PeffLoadStats,
) -> Result<Vec<Row>, PeffError> {
    let md5 = sequence_md5(seq);
    let malformed = |key: &str, value: &str| PeffError::Malformed {
        accession: accession.to_string(),
        key: key.to_string(),
        value: value.to_string(),
    };
    let pos = |s: &str, key: &str, value: &str| -> Result<u32, PeffError> {
        s.parse::<u32>()
            .ok()
            .filter(|&p| p >= 1)
            .ok_or_else(|| malformed(key, value))
    };
    // The residues of a 1-based inclusive span, empty when it does not fit.
    let residues = |start: u32, end: u32| -> String {
        let (s, e) = (start as usize - 1, end as usize);
        if s <= e && e <= seq.len() {
            String::from_utf8_lossy(&seq[s..e]).into_owned()
        } else {
            String::new()
        }
    };
    let row = |id: String| Row {
        id,
        reference: accession.to_string(),
        md5: md5.clone(),
        edits: Vec::new(),
        mods: Vec::new(),
        processing: None,
        source: "PEFF".to_string(),
        tier: PriorTier::Catalogue,
    };
    let mut rows = Vec::new();
    for (key, value) in keys {
        let parsed = || groups(value).ok_or_else(|| malformed(key, value));
        match key.as_str() {
            "VariantSimple" => {
                for g in parsed()? {
                    let [p, alt, rest @ ..] = g.as_slice() else {
                        return Err(malformed(key, value));
                    };
                    let p = pos(p, key, value)?;
                    let original = residues(p, p);
                    let id = rest.first().filter(|t| !t.is_empty()).map_or_else(
                        || {
                            format!(
                                "VAR_{}",
                                change_label(original.as_bytes(), p, p, alt.as_bytes())
                            )
                        },
                        |t| t.to_string(),
                    );
                    let mut r = row(id);
                    r.edits = vec![(p - 1, p, original, alt.to_string())];
                    rows.push(r);
                }
            }
            "VariantComplex" => {
                for g in parsed()? {
                    let [s, e, alt, rest @ ..] = g.as_slice() else {
                        return Err(malformed(key, value));
                    };
                    let (s, e) = (pos(s, key, value)?, pos(e, key, value)?);
                    let original = residues(s, e);
                    let id = rest.first().filter(|t| !t.is_empty()).map_or_else(
                        || {
                            format!(
                                "VAR_{}",
                                change_label(original.as_bytes(), s, e, alt.as_bytes())
                            )
                        },
                        |t| t.to_string(),
                    );
                    let mut r = row(id);
                    r.edits = vec![(s - 1, e, original, alt.to_string())];
                    rows.push(r);
                }
            }
            "ModResUnimod" | "ModRes" | "ModResPsi" => {
                for g in parsed()? {
                    let [p, acc, name, ..] = g.as_slice() else {
                        return Err(malformed(key, value));
                    };
                    let p = pos(p, key, value)?;
                    let site = if key == "ModResUnimod" {
                        model::unimod::bundled()
                            .by_accession(acc)
                            .map(|m| (m.mono_mass, SiteAnchor::Residue))
                    } else {
                        ptm_mass(clean_description(name)).map(|(_, mass, anchor)| (mass, anchor))
                    };
                    let Some((mass, anchor)) = site else {
                        stats.sites_without_mass += 1;
                        continue;
                    };
                    let anchor = match anchor {
                        SiteAnchor::Residue => "residue",
                        SiteAnchor::NTerminus => "protein_n_term",
                        SiteAnchor::CTerminus => "protein_c_term",
                    };
                    let label = if acc.is_empty() { name } else { acc };
                    let mut r = row(format!("SITE_{p}_{label}"));
                    r.mods = vec![(p - 1, residues(p, p), mass, anchor.to_string())];
                    rows.push(r);
                }
            }
            "Processed" => {
                for g in parsed()? {
                    let [s, e, .., name] = g.as_slice() else {
                        return Err(malformed(key, value));
                    };
                    let name = name.to_ascii_lowercase();
                    if ["signal", "transit", "propeptide", "initiator"]
                        .iter()
                        .any(|w| name.contains(w))
                    {
                        stats.processing_not_chain += 1;
                        continue;
                    }
                    let (s, e) = (pos(s, key, value)?, pos(e, key, value)?);
                    let mut r = row(format!("CHAIN_{s}-{e}"));
                    r.processing = Some((s - 1, e));
                    rows.push(r);
                }
            }
            _ => {}
        }
    }
    Ok(rows)
}

#[cfg(test)]
mod tests {
    use super::*;
    use model::{EntryKind, Enzyme};

    const KO: KnockoffOptions = KnockoffOptions {
        per_variant: 4,
        enzyme: Enzyme::Trypsin,
    };

    fn peff() -> String {
        let seq = String::from_utf8(b"MGGGDLGGGK".repeat(8)).unwrap();
        format!(
            "# PEFF 1.0\n# //\n# DbName=test\n# //\n\
             >sp:P00009 \\DbUniqueId=P00009 \\PName=Knockoff test \\Length=80 \
             \\VariantSimple=(5|N|VAR_000009)(15|E) \\VariantComplex=(25|26|) \
             \\ModResUnimod=(10|UNIMOD:121|GG) \\ModResPsi=(20|MOD:00000|not a UniProt name) \
             \\Processed=(1|10|PEFF:0001021|signal peptide)(11|80|PEFF:0001022|mature protein)\n\
             {}\n{}\n",
            &seq[..40],
            &seq[40..]
        )
    }

    #[test]
    fn header_keys_become_entries_and_sites() {
        let (db, stats) = load_peff(peff().as_bytes(), Some(KO)).unwrap();
        assert_eq!(stats.entries, 1);
        assert_eq!(db.proteins[0].accession, "sp:P00009");
        assert_eq!(db.proteins[0].description, "Knockoff test");
        let c = &stats.claims;
        assert_eq!((c.variant_entries, c.chain_entries, c.sites), (3, 1, 1));
        assert!(c.knockoff_entries > 0);
        assert_eq!(
            (stats.sites_without_mass, stats.processing_not_chain),
            (1, 1)
        );
        let accs: Vec<&str> = db.proteins.iter().map(|p| p.accession.as_str()).collect();
        assert!(accs.contains(&"sp:P00009_VAR_000009"));
        assert!(accs.contains(&"sp:P00009_VAR_D15E"));
        assert!(accs.contains(&"sp:P00009_VAR_DL25-26del"));
        let gg = db.proteins[0].mod_sites();
        assert_eq!((gg.len(), gg[0].position), (1, 9));
        let chain = db
            .proteins
            .iter()
            .find(|p| p.entry_kind() == EntryKind::Chain)
            .unwrap();
        assert_eq!(chain.sequence.len(), 70);
        let v = db
            .proteins
            .iter()
            .find(|p| p.accession == "sp:P00009_VAR_000009")
            .unwrap();
        assert_eq!(&*v.origin().unwrap().family, "P00009");
    }

    #[test]
    fn malformed_values_are_errors() {
        let text = ">sp:P1 \\VariantSimple=(x|N)\nMGGGDLGGGK\n";
        assert!(matches!(
            load_peff(text.as_bytes(), None),
            Err(PeffError::Malformed { .. })
        ));
        let text = "> \\PName=x\nMGGGDLGGGK\n";
        assert!(matches!(
            load_peff(text.as_bytes(), None),
            Err(PeffError::EmptyAccession { line: 1 })
        ));
    }

    #[test]
    fn a_variant_that_does_not_fit_is_counted() {
        // Position 5 holds D, so a variant at 81 is outside the sequence.
        let text = ">sp:P1 \\VariantSimple=(81|N)\nMGGGDLGGGK\n";
        let (db, stats) = load_peff(text.as_bytes(), None).unwrap();
        assert_eq!(db.proteins.len(), 1);
        assert_eq!(stats.claims.dropped.values().sum::<usize>(), 1);
    }
}
