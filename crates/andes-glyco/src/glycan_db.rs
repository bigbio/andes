// Clean-room pGlyco `.gdb` glycan database loader.
//
// Loads a pGlyco-style canonical-string glycan database and preserves the tree
// structure just far enough to recover the one structural fact the glycan-first
// index needs: whether a fucose is **core-fucose** (a direct child of the
// reducing-end monosaccharide, retained on every core Y ion) or antenna fucose.
// Branching topology is otherwise collapsed to residue counts, because the
// downstream search indexes glycans by composition.

use crate::glycan_mass::{FUC, HEX, HEXNAC, NEUAC, NEUGC};

/// A single glycan (residue counts + the core-fucose flag + monoisotopic mass).
#[derive(Debug, Clone, PartialEq)]
pub struct GlycanComp {
    pub hexnac: u8,
    pub hex: u8,
    pub fuc: u8,
    pub neuac: u8,
    pub neugc: u8,
    /// 1 when one Fuc is core-fucose (a direct child of the reducing-end
    /// monosaccharide), 0 otherwise. Only meaningful for N-glycans; the search
    /// core gates it behind `GlycanCore::core_fucose()`.
    pub core_fuc: u8,
    pub mass: f64,
}

/// Error loading a `.gdb` glycan database.
#[derive(Debug, Clone, PartialEq)]
pub enum GdbLoadError {
    /// The file has no header line.
    Empty,
    /// The file had a valid header line but no glycan trees (e.g. header-only).
    NoGlycans,
    /// A glycan line contained a monosaccharide symbol outside the N-glycan set.
    UnknownSymbol { line: usize, symbol: String },
    /// A glycan line was not a well-formed parenthesized tree.
    MalformedTree { line: usize, reason: &'static str },
    /// `--glyco-species` named a species we do not bundle.
    UnknownSpecies { species: String },
}

impl std::fmt::Display for GdbLoadError {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self {
            GdbLoadError::Empty => write!(f, "glycan .gdb is empty (missing header line)"),
            GdbLoadError::NoGlycans => {
                write!(f, "glycan .gdb contains no glycan trees (header-only?)")
            }
            GdbLoadError::UnknownSymbol { line, symbol } => {
                write!(f, "glycan .gdb line {line}: unknown monosaccharide symbol {symbol:?}")
            }
            GdbLoadError::MalformedTree { line, reason } => {
                write!(f, "glycan .gdb line {line}: malformed tree ({reason})")
            }
            GdbLoadError::UnknownSpecies { species } => {
                write!(f, "unknown --glyco-species {species:?}")
            }
        }
    }
}

impl std::error::Error for GdbLoadError {}

/// A node of a parsed glycan tree: a residue symbol plus its child subtrees.
struct Node {
    sym: u8,
    children: Vec<Node>,
}

/// Map a residue symbol to its index in the count tally `[H, N, F, A, G]`.
#[inline]
fn symbol_index(sym: u8) -> Option<usize> {
    match sym {
        b'H' => Some(0),
        b'N' => Some(1),
        b'F' => Some(2),
        b'A' => Some(3),
        b'G' => Some(4),
        _ => None,
    }
}

/// Maximum nesting depth of a glycan tree. Real N-glycans nest a handful of
/// levels; a corrupted or malicious `.gdb` with thousands of open parens would
/// otherwise overflow the recursive-descent parser's call stack.
const MAX_GLYCAN_DEPTH: usize = 64;

/// Recursive-descent parse of one parenthesized node `( sym child* )`.
fn parse_node(
    bytes: &[u8],
    pos: &mut usize,
    line: usize,
    depth: usize,
) -> Result<Node, GdbLoadError> {
    if depth > MAX_GLYCAN_DEPTH {
        return Err(GdbLoadError::MalformedTree {
            line,
            reason: "glycan tree exceeds maximum depth",
        });
    }
    if *pos >= bytes.len() || bytes[*pos] != b'(' {
        return Err(GdbLoadError::MalformedTree {
            line,
            reason: "expected '('",
        });
    }
    *pos += 1;
    if *pos >= bytes.len() {
        return Err(GdbLoadError::MalformedTree {
            line,
            reason: "missing residue symbol",
        });
    }
    let sym = bytes[*pos];
    if symbol_index(sym).is_none() {
        return Err(GdbLoadError::UnknownSymbol {
            line,
            symbol: (sym as char).to_string(),
        });
    }
    *pos += 1;
    let mut children = Vec::new();
    loop {
        if *pos >= bytes.len() {
            return Err(GdbLoadError::MalformedTree {
                line,
                reason: "missing ')'",
            });
        }
        match bytes[*pos] {
            b'(' => children.push(parse_node(bytes, pos, line, depth + 1)?),
            b')' => {
                *pos += 1;
                return Ok(Node { sym, children });
            }
            _ => {
                return Err(GdbLoadError::MalformedTree {
                    line,
                    reason: "child residue must be parenthesized",
                })
            }
        }
    }
}

/// Sum every residue of a subtree into `counts`.
fn tally(node: &Node, counts: &mut [u8; 5]) {
    if let Some(i) = symbol_index(node.sym) {
        counts[i] += 1;
    }
    for c in &node.children {
        tally(c, counts);
    }
}

/// Parse one glycan tree starting at `*pos`, recovering `core_fuc` from the tree.
/// Advances `*pos` past the tree. A `.gdb` line may carry several concatenated
/// trees (the pGlyco `*-multi` databases do), so callers loop until the line is
/// consumed rather than assuming one tree per line.
fn parse_glycan_at(
    bytes: &[u8],
    pos: &mut usize,
    line_no: usize,
) -> Result<GlycanComp, GdbLoadError> {
    let root = parse_node(bytes, pos, line_no, 0)?;
    let mut counts = [0u8; 5];
    tally(&root, &mut counts);
    let hex = counts[0];
    let hexnac = counts[1];
    let fuc = counts[2];
    let neuac = counts[3];
    let neugc = counts[4];
    // Core-fucose = a Fuc that is a direct child of the reducing-end residue
    // (the outermost symbol). Antenna/Lewis fucose sits on an inner node.
    let core_fuc = root.children.iter().any(|c| c.sym == b'F');
    let mass = hexnac as f64 * HEXNAC
        + hex as f64 * HEX
        + fuc as f64 * FUC
        + neuac as f64 * NEUAC
        + neugc as f64 * NEUGC;
    Ok(GlycanComp {
        hexnac,
        hex,
        fuc,
        neuac,
        neugc,
        core_fuc: core_fuc as u8,
        mass,
    })
}

/// Load a pGlyco-style `.gdb` glycan database, preserving the tree structure.
///
/// Format:
/// ```text
/// H,N,A,G,F          <- header: monosaccharide symbols (order irrelevant)
/// (N(F)(N(H(H)(H)))) <- glycan 1 (nested tree of symbols; reducing end outermost)
/// ```
/// Symbols: `N`=HexNAc, `H`=Hex, `F`=Fuc, `A`=NeuAc, `G`=NeuGc. Each glycan is a
/// nested parenthesized tree; a line may hold several concatenated trees (the
/// pGlyco `*-multi` databases concatenate structures), each parsed separately. We
/// walk each tree once to count residues and recover `core_fuc` (a `F` that is a
/// direct child of the outermost symbol). Returns
/// glycans sorted by mass ascending (tiebroken by composition, then `core_fuc`)
/// and deduplicated by `(counts, core_fuc)` — so a core-fucosylated glycan and an
/// antenna-fucosylated isomer of the same composition stay distinct.
pub fn load_glycan_gdb(content: &str) -> Result<Vec<GlycanComp>, GdbLoadError> {
    let mut lines = content.lines();
    let _header = lines.next().ok_or(GdbLoadError::Empty)?;
    let mut out: Vec<GlycanComp> = Vec::new();
    let mut seen: std::collections::HashSet<(u8, u8, u8, u8, u8, u8)> =
        std::collections::HashSet::new();
    for (i, raw) in lines.enumerate() {
        let line = raw.trim();
        if line.is_empty() {
            continue;
        }
        let line_no = i + 2; // header is line 1
        let bytes = line.as_bytes();
        let mut pos = 0usize;
        while pos < bytes.len() {
            // Skip stray whitespace between concatenated trees (none in practice,
            // but harmless to tolerate).
            while pos < bytes.len() && (bytes[pos] as char).is_whitespace() {
                pos += 1;
            }
            if pos >= bytes.len() {
                break;
            }
            let comp = parse_glycan_at(bytes, &mut pos, line_no)?;
            if seen.insert((
                comp.hexnac,
                comp.hex,
                comp.fuc,
                comp.neuac,
                comp.neugc,
                comp.core_fuc,
            )) {
                out.push(comp);
            }
        }
    }
    out.sort_by(|a, b| {
        a.mass
            .to_bits()
            .cmp(&b.mass.to_bits())
            .then(a.hexnac.cmp(&b.hexnac))
            .then(a.hex.cmp(&b.hex))
            .then(a.fuc.cmp(&b.fuc))
            .then(a.neuac.cmp(&b.neuac))
            .then(a.neugc.cmp(&b.neugc))
            .then(a.core_fuc.cmp(&b.core_fuc))
    });
    if out.is_empty() {
        return Err(GdbLoadError::NoGlycans);
    }
    Ok(out)
}

/// Load one of the bundled species-specific N-glycan databases by name.
///
/// The databases are compiled in via `include_str!` (see `glycan-db/`); `species`
/// is the kebab-case name exposed by the `--glyco-species` flag (`human`,
/// `human-multi`, `mouse`, `mouse-large`, `plant`, `plant-multi`,
/// `high-mannose`). This is the `--glyco-species` path of the glycan-first
/// search; an explicit `--glyco-glycan-gdb` file bypasses it entirely.
pub fn load_species_glycan_db(species: &str) -> Result<Vec<GlycanComp>, GdbLoadError> {
    let content = match species {
        "human" => include_str!("../glycan-db/pGlyco-N-Human.gdb"),
        "human-multi" => include_str!("../glycan-db/pGlyco-N-Human-multi.gdb"),
        "mouse" => include_str!("../glycan-db/pGlyco-N-Mouse.gdb"),
        "mouse-large" => include_str!("../glycan-db/pGlyco-N-Mouse-large.gdb"),
        "high-mannose" => include_str!("../glycan-db/pGlyco-N-HighMannose.gdb"),
        _ => {
            return Err(GdbLoadError::UnknownSpecies {
                species: species.to_string(),
            })
        }
    };
    load_glycan_gdb(content)
}

#[cfg(test)]
mod tests {
    use super::*;

    fn load_one(s: &str) -> GlycanComp {
        let content = format!("H,N,A,G,F\n{s}\n");
        let list = load_glycan_gdb(&content).unwrap();
        assert_eq!(list.len(), 1, "expected one glycan, got {list:?}");
        list[0].clone()
    }

    #[test]
    fn trimannosyl_core_has_no_core_fucose() {
        let g = load_one("(N(N(H(H)(H))))");
        assert_eq!((g.hexnac, g.hex, g.fuc, g.core_fuc), (2, 3, 0, 0));
    }

    #[test]
    fn core_fucose_is_direct_child_of_root() {
        // (N(F)(N(H(H)(H)))) — F on the innermost GlcNAc.
        let g = load_one("(N(F)(N(H(H)(H))))");
        assert_eq!((g.hexnac, g.hex, g.fuc, g.core_fuc), (2, 3, 1, 1));
    }

    #[test]
    fn antenna_fucose_is_not_core_fucose() {
        // (N(N(H(H)(H(N(F)))))) — F on an antenna GlcNAc (Lewis), not the root.
        let g = load_one("(N(N(H(H)(H(N(F))))))");
        assert_eq!((g.hexnac, g.hex, g.fuc, g.core_fuc), (3, 3, 1, 0));
    }

    #[test]
    fn core_and_antenna_fucose_coexist() {
        // A core-Fuc AND an antenna-Fuc on the same glycan.
        let g = load_one("(N(F)(N(H(H)(H(N(F))))))");
        assert_eq!((g.hexnac, g.hex, g.fuc, g.core_fuc), (3, 3, 2, 1));
    }

    #[test]
    fn sialylated_biantennary_core_fucosylated() {
        // HexNAc4 Hex5 Fuc1 NeuAc2: two sialylated antennae + core-Fuc.
        let g = load_one("(N(F)(N(H(H(N(H(A))))(H(N(H(A)))))))");
        assert_eq!(
            (g.hexnac, g.hex, g.fuc, g.neuac, g.core_fuc),
            (4, 5, 1, 2, 1)
        );
    }

    #[test]
    fn core_fuc_isomers_are_kept_distinct() {
        // Same composition (N2H3F1) but core-Fuc vs man-Fuc: two distinct glycans.
        let content = "H,N,A,G,F\n(N(F)(N(H(H)(H))))\n(N(N(H(H)(H)(F))))\n";
        let list = load_glycan_gdb(content).unwrap();
        assert_eq!(list.len(), 2);
        assert_eq!(list.iter().filter(|g| g.core_fuc == 1).count(), 1);
        assert_eq!(list.iter().filter(|g| g.core_fuc == 0).count(), 1);
    }

    #[test]
    fn unbalanced_tree_is_rejected() {
        let err = load_glycan_gdb("H,N,A,G,F\n(N(N(H)\n").unwrap_err();
        assert!(matches!(err, GdbLoadError::MalformedTree { .. }));
    }

    #[test]
    fn trailing_characters_are_rejected() {
        let err = load_glycan_gdb("H,N,A,G,F\n(N))\n").unwrap_err();
        assert!(matches!(err, GdbLoadError::MalformedTree { .. }));
    }

    #[test]
    fn concatenated_trees_on_one_line_are_split() {
        // pGlyco `*-multi` databases concatenate two trees with no separator.
        let content = "H,N,A,G,F\n(N(N(H(H)(H))))(N(N(H(H(H))(H(H)))))\n";
        let list = load_glycan_gdb(content).unwrap();
        assert_eq!(list.len(), 2);
    }

    #[test]
    fn deeply_nested_tree_is_rejected() {
        // A pathological tree nested far past any real glycan; must not overflow
        // the parser's call stack.
        let mut s = String::from("H,N,A,G,F\n");
        for _ in 0..200 {
            s.push_str("(N");
        }
        for _ in 0..200 {
            s.push(')');
        }
        s.push('\n');
        let err = load_glycan_gdb(&s).unwrap_err();
        assert!(matches!(err, GdbLoadError::MalformedTree { .. }));
    }

    #[test]
    fn unknown_symbol_is_rejected() {
        let err = load_glycan_gdb("H,N,A,G,F\n(X)\n").unwrap_err();
        assert!(matches!(err, GdbLoadError::UnknownSymbol { symbol, .. } if symbol == "X"));
    }

    #[test]
    fn empty_file_is_rejected() {
        assert_eq!(load_glycan_gdb("").unwrap_err(), GdbLoadError::Empty);
    }

    #[test]
    fn header_only_file_is_rejected() {
        let err = load_glycan_gdb("H,N,A,G,F\n").unwrap_err();
        assert!(matches!(err, GdbLoadError::NoGlycans));
    }

    #[test]
    fn bundled_species_databases_load() {
        for sp in [
            "human",
            "human-multi",
            "mouse",
            "mouse-large",
            "high-mannose",
        ] {
            let list = load_species_glycan_db(sp).unwrap();
            assert!(!list.is_empty(), "{sp}: empty glycan list");
        }
        // The mouse and human databases should deduplicate to well over a thousand
        // compositions; a truncated bundle would collapse far below this.
        assert!(load_species_glycan_db("mouse").unwrap().len() >= 1000);
        assert!(load_species_glycan_db("human").unwrap().len() >= 1000);
        assert!(load_species_glycan_db("high-mannose").unwrap().len() >= 10);
    }

    #[test]
    fn unknown_species_is_rejected() {
        let err = load_species_glycan_db("zebra").unwrap_err();
        assert!(matches!(err, GdbLoadError::UnknownSpecies { .. }));
    }
}
