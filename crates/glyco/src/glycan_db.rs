// Clean-room glycan database loader: canonical-string `.gdb` and WURCS 2.0.
//
// Loads a canonical-string `.gdb` glycan database (or a WURCS 2.0 CSV)
// and retains the FULL tree topology, normalised at parse time (children sorted by
// canonical form). The topology is the glycan's identity: the downstream search
// dedups by tree and generates Y-ions from the rooted tree, so same-composition
// isomers (bisecting vs antenna GlcNAc, different antenna branching, core vs
// Lewis fucose) stay distinct candidates and score differently.

use std::sync::Arc;

use crate::glycan_mass::{FUC, HEX, HEXNAC, NEUAC, NEUGC};

/// Reserved `topology_id` for the de-novo (unannotated) candidate branch. Real
/// ids are FNV-1a hashes of canonical strings; the loader asserts none of them
/// equal this sentinel.
pub const DENOVO_TOPOLOGY_ID: u64 = u64::MAX;

/// A single glycan: residue counts + core-fucose flag + monoisotopic mass + the
/// FULL tree topology. Equality/hash are topology identity (the count/mass fields
/// are functions of the tree, so they are stored for cheap access but do not
/// participate in equality).
#[derive(Debug, Clone)]
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
    /// The full tree topology (reducing end outermost). Shared (`Arc`) so
    /// `GlycanComp::clone` — done per candidate in `db_branch` — is O(1).
    pub tree: Arc<GlycanTree>,
    /// Stable, owned, `Copy` topology identity for the winner-collapse key.
    pub topology_id: u64,
}

impl PartialEq for GlycanComp {
    fn eq(&self, other: &Self) -> bool {
        self.tree == other.tree
    }
}

impl Eq for GlycanComp {}

impl std::hash::Hash for GlycanComp {
    fn hash<H: std::hash::Hasher>(&self, state: &mut H) {
        self.tree.hash(state);
    }
}

impl GlycanComp {
    /// Derive every count/mass/identity field from a (normalised) tree.
    pub fn from_tree(tree: GlycanTree) -> GlycanComp {
        let counts = tree.counts();
        let hexnac = counts[1];
        let hex = counts[0];
        let fuc = counts[2];
        let neuac = counts[3];
        let neugc = counts[4];
        let core_fuc = tree.children.iter().any(|c| c.sym == b'F') as u8;
        let mass = hexnac as f64 * HEXNAC
            + hex as f64 * HEX
            + fuc as f64 * FUC
            + neuac as f64 * NEUAC
            + neugc as f64 * NEUGC;
        let topology_id = tree.topology_id();
        GlycanComp {
            hexnac,
            hex,
            fuc,
            neuac,
            neugc,
            core_fuc,
            mass,
            tree: Arc::new(tree),
            topology_id,
        }
    }

    /// Parse a canonical S-expression and build the glycan in one step.
    pub fn from_tree_str(s: &str) -> Result<GlycanComp, GdbLoadError> {
        Ok(GlycanComp::from_tree(GlycanTree::parse(s)?))
    }

    /// Build a composition from raw residue counts when the topology is unknown
    /// (composition-only consumers: RT-offset prediction, tests). Synthesises a
    /// deterministic canonical tree — the reducing-end HexNAc, its optional core
    /// Fuc, then a single LINEAR chain in the biosynthetic order chitobiose →
    /// trimannose → remaining residues — so `canonical()`/`topology_id()`/
    /// `y_ion_retained_masses()` stay well-defined, stable for a given count
    /// vector, and match the legacy cumulative Y-ladder. `fuc` is the TOTAL
    /// fucose count and `core_fuc` (0/1) says whether one of them sits on the
    /// reducing-end HexNAc.
    pub fn from_counts(
        hexnac: u8,
        hex: u8,
        fuc: u8,
        neuac: u8,
        neugc: u8,
        core_fuc: u8,
    ) -> GlycanComp {
        let mass = hexnac as f64 * HEXNAC
            + hex as f64 * HEX
            + fuc as f64 * FUC
            + neuac as f64 * NEUAC
            + neugc as f64 * NEUGC;
        let core_fuc = core_fuc.min(1);
        let core_hn = hexnac.min(2);
        let core_h = hex.min(3);
        // Chain residues after the reducing-end HexNAc, in biosynthetic order.
        let mut chain: Vec<u8> = Vec::with_capacity((hexnac + hex + fuc + neuac + neugc) as usize);
        chain.extend(std::iter::repeat_n(
            b'N',
            core_hn.saturating_sub(1) as usize,
        )); // second chitobiose GlcNAc
        chain.extend(std::iter::repeat_n(b'H', core_h as usize)); // trimannosyl core
        chain.extend(std::iter::repeat_n(b'N', (hexnac - core_hn) as usize)); // remaining HexNAc (antenna)
        chain.extend(std::iter::repeat_n(b'H', (hex - core_h) as usize)); // remaining Hex (antenna)
        chain.extend(std::iter::repeat_n(b'F', (fuc - core_fuc) as usize)); // antenna Fuc (core Fuc is already a root child)
        chain.extend(std::iter::repeat_n(b'A', neuac as usize));
        chain.extend(std::iter::repeat_n(b'G', neugc as usize));

        let mut s = String::from("(N");
        if core_fuc > 0 {
            s.push_str("(F)");
        }
        // Nest the chain innermost-first: `(sym0 (sym1 (sym2 …)))`.
        let mut chain_str = String::new();
        for &sym in chain.iter().rev() {
            chain_str = format!("({}{})", sym as char, chain_str);
        }
        s.push_str(&chain_str);
        s.push(')');

        let tree = GlycanTree::parse(&s).expect("synthetic counts tree is always well-formed");
        let topology_id = tree.topology_id();
        GlycanComp {
            hexnac,
            hex,
            fuc,
            neuac,
            neugc,
            core_fuc,
            mass,
            tree: Arc::new(tree),
            topology_id,
        }
    }

    /// Y-ion retained masses (single-glycosidic-bond-cleavage fragments), ascending.
    /// See [`GlycanTree::y_ion_retained_masses`].
    pub fn y_ion_retained_masses(&self) -> Vec<f64> {
        self.tree.y_ion_retained_masses()
    }
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
    /// A glycan line held more than one complete tree. Upstream `*-multi`
    /// databases list combinations of glycans for peptides with several
    /// glycosites; andes places one glycan per peptide, so such a line has no
    /// meaning here.
    MultiSite { line: usize, trees: usize },
    /// A WURCS 2.0 structure string could not be parsed into a tree.
    MalformedWurcs { line: usize, reason: &'static str },
    /// A WURCS 2.0 residue code is outside the supported N-glycan set
    /// (`H`/`N`/`F`/`A`/`G`), e.g. a sulfate, uronic acid or pentose.
    UnknownWurcsResidue { line: usize, code: String },
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
                write!(
                    f,
                    "glycan .gdb line {line}: unknown monosaccharide symbol {symbol:?}"
                )
            }
            GdbLoadError::MalformedTree { line, reason } => {
                write!(f, "glycan .gdb line {line}: malformed tree ({reason})")
            }
            GdbLoadError::MultiSite { line, trees } => write!(
                f,
                "glycan .gdb line {line}: {trees} glycan trees on one line (a multi-site \
                 combination); andes places one glycan per peptide, so use a single-site \
                 database"
            ),
            GdbLoadError::MalformedWurcs { line, reason } => {
                write!(
                    f,
                    "glycan WURCS line {line}: malformed structure ({reason})"
                )
            }
            GdbLoadError::UnknownWurcsResidue { line, code } => {
                write!(
                    f,
                    "glycan WURCS line {line}: unsupported residue code {code:?}"
                )
            }
            GdbLoadError::UnknownSpecies { species } => {
                write!(f, "unknown --glyco-species {species:?}")
            }
        }
    }
}

impl std::error::Error for GdbLoadError {}

/// A glycan tree node: a residue symbol plus its child subtrees. Trees are
/// NORMALISED at parse time (children sorted by canonical serialization), so
/// `PartialEq`/`Hash`/`canonical()`/`topology_id()` are all order-independent.
#[derive(Debug, Clone, PartialEq, Eq, Hash)]
pub struct GlycanTree {
    pub sym: u8,
    pub children: Vec<GlycanTree>,
}

impl GlycanTree {
    /// Parse one canonical S-expression `(N(F)(N(H(H)(H))))` (reducing end outermost).
    pub fn parse(s: &str) -> Result<GlycanTree, GdbLoadError> {
        let bytes = s.as_bytes();
        let mut pos = 0usize;
        let tree = parse_node(bytes, &mut pos, 1, 1)?.check_counts(1)?;
        if pos < bytes.len() {
            return Err(GdbLoadError::MalformedTree {
                line: 1,
                reason: "trailing characters after the glycan tree",
            });
        }
        Ok(tree)
    }

    /// Parse one WURCS 2.0 structure string (a `StructureInformation`
    /// value) into the same canonical tree `parse` builds from an S-expression.
    ///
    /// Only the N-glycan backbone is retained: each residue is reduced to its
    /// monosaccharide class (`H`=Hex, `N`=HexNAc, `F`=Fuc, `A`=NeuAc, `G`=NeuGc),
    /// and anomericity, ring form, substituents and linkage positions are dropped.
    /// The linkage block is resolved into parent→child edges (the child links
    /// through its anomeric carbon — C1 for aldoses, C2 for sialic acid) and the
    /// connected graph is rooted at the single atom with children and no parent
    /// (the reducing end). Isolated atoms — WURCS `+` repeating-unit placeholders —
    /// are ignored.
    pub fn from_wurcs(wurcs: &str) -> Result<GlycanTree, GdbLoadError> {
        wurcs_to_tree(wurcs, 1)
    }

    /// Residue mass of a symbol (`H`/`N`/`F`/`A`/`G`). The parser rejects any other
    /// symbol upstream, so a non-matching byte is a bug and panics.
    fn residue_mass(sym: u8) -> f64 {
        match sym {
            b'H' => HEX,
            b'N' => HEXNAC,
            b'F' => FUC,
            b'A' => NEUAC,
            b'G' => NEUGC,
            _ => unreachable!("unknown glycan symbol (parser rejects it upstream)"),
        }
    }

    /// Monoisotopic mass of this subtree (sum of residue masses).
    pub fn mass(&self) -> f64 {
        let mut m = GlycanTree::residue_mass(self.sym);
        for c in &self.children {
            m += c.mass();
        }
        m
    }

    /// Tally `[H, N, F, A, G]` counts over this subtree, each saturating at
    /// 255. The loaders reject trees with more (see [`GlycanTree::checked_counts`]).
    pub fn counts(&self) -> [u8; 5] {
        let mut counts = [0u8; 5];
        tally(self, &mut counts);
        counts
    }

    /// The `[H, N, F, A, G]` counts, or `None` when a residue occurs more than
    /// 255 times in the subtree.
    pub fn checked_counts(&self) -> Option<[u8; 5]> {
        let mut counts = [0u8; 5];
        tally(self, &mut counts).then_some(counts)
    }

    /// `tree` itself, or [`GdbLoadError::MalformedTree`] on `line` when a
    /// residue occurs more than 255 times.
    fn check_counts(self, line: usize) -> Result<GlycanTree, GdbLoadError> {
        match self.checked_counts() {
            Some(_) => Ok(self),
            None => Err(GdbLoadError::MalformedTree {
                line,
                reason: "more than 255 residues of one kind",
            }),
        }
    }

    /// Canonical serialization `(sym child*)` with children sorted by their
    /// canonical forms (deterministic for any child order).
    pub fn canonical(&self) -> String {
        let mut out = String::new();
        self.write_canonical(&mut out);
        out
    }

    fn write_canonical(&self, out: &mut String) {
        out.push('(');
        out.push(self.sym as char);
        let mut kids: Vec<&GlycanTree> = self.children.iter().collect();
        kids.sort_by_key(|a| a.canonical());
        for k in kids {
            k.write_canonical(out);
        }
        out.push(')');
    }

    /// Deterministic 64-bit identity of the topology: FNV-1a over the canonical
    /// string. NOT `DefaultHasher` (whose output is unspecified and version-unstable).
    pub fn topology_id(&self) -> u64 {
        let mut h: u64 = 0xcbf2_9ce4_8422_2325;
        for b in self.canonical().as_bytes() {
            h ^= *b as u64;
            h = h.wrapping_mul(0x0000_0100_0000_01b3);
        }
        h
    }

    /// Y-ion retained masses: the reducing-end fragment from each SINGLE
    /// glycosidic-bond cleavage, ascending and deduplicated by tolerance.
    ///
    /// This follows the standard B/Y convention (Domon &
    /// Costello): a Y ion keeps the reducing end, so cleaving the bond above a
    /// child `c` retains everything except `c`'s subtree — mass `total −
    /// subtree(c)`. There is one such ion per glycosidic bond, so the set is
    /// O(n) for n residues — NOT the combinatorial set of multi-cleavage
    /// "internal" fragments, which would over-generate physically-unlikely ions
    /// and blow up the scoring ladder's rung count for branched glycans.
    ///
    /// Y0 (the bare peptide, retained mass 0) is NOT included — callers add it.
    /// The intact glycan (retained == full mass) IS included as the top rung.
    pub fn y_ion_retained_masses(&self) -> Vec<f64> {
        let total = self.mass();
        let mut v: Vec<f64> = vec![total]; // the intact glycan (top rung)
        single_cleavage_y_masses(self, total, &mut v);
        v.sort_by(|a, b| a.total_cmp(b));
        let mut out: Vec<f64> = Vec::with_capacity(v.len());
        for m in v {
            match out.last() {
                // Merge only true duplicates (cutting two symmetric antenna
                // bonds yields the same reducing-end mass; float addition-order
                // noise is ~1e-11).
                Some(&last) if (m - last) <= 1e-6 => {}
                _ => out.push(m),
            }
        }
        out
    }
}

/// Append the reducing-end masses from every SINGLE glycosidic-bond cleavage of
/// `node`'s subtree to `out`. For each edge (parent → child), cleaving it retains
/// `total − subtree(child)` on the reducing-end side.
fn single_cleavage_y_masses(node: &GlycanTree, total: f64, out: &mut Vec<f64>) {
    for child in &node.children {
        out.push(total - child.mass());
        single_cleavage_y_masses(child, total, out);
    }
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
) -> Result<GlycanTree, GdbLoadError> {
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
                // Normalise: sort children by canonical form so equal topologies
                // compare/hash equal regardless of child order in the source.
                children.sort_by_cached_key(|c: &GlycanTree| c.canonical());
                return Ok(GlycanTree { sym, children });
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

/// Sum every residue of a subtree into `counts`; `false` when a count reached
/// its limit (it stays at 255).
fn tally(node: &GlycanTree, counts: &mut [u8; 5]) -> bool {
    let mut fits = true;
    if let Some(i) = symbol_index(node.sym) {
        match counts[i].checked_add(1) {
            Some(n) => counts[i] = n,
            None => fits = false,
        }
    }
    for c in &node.children {
        fits &= tally(c, counts);
    }
    fits
}

/// Parse one glycan tree starting at `*pos`. `GlycanComp::from_tree` recovers the
/// counts, `core_fuc`, mass and topology identity from the tree. Advances `*pos`
/// past the tree.
fn parse_glycan_at(
    bytes: &[u8],
    pos: &mut usize,
    line_no: usize,
) -> Result<GlycanComp, GdbLoadError> {
    let root = parse_node(bytes, pos, line_no, 0)?.check_counts(line_no)?;
    Ok(GlycanComp::from_tree(root))
}

/// Load a `.gdb` glycan database, retaining the full tree topology.
///
/// Format:
/// ```text
/// H,N,A,G,F          <- header: monosaccharide symbols (order irrelevant)
/// (N(F)(N(H(H)(H)))) <- glycan 1 (nested tree of symbols; reducing end outermost)
/// ```
/// Symbols: `N`=HexNAc, `H`=Hex, `F`=Fuc, `A`=NeuAc, `G`=NeuGc. Each line is one
/// nested parenthesized tree. A line with several concatenated trees is an upstream
/// `*-multi` multi-site combination (one glycan per glycosite of the same
/// peptide); the search places one glycan per peptide, so such a line is
/// rejected with [`GdbLoadError::MultiSite`] rather than read as separate
/// glycans. Returns glycans sorted by mass ascending (tiebroken by composition,
/// then `core_fuc`, then `topology_id`) and deduplicated by **topology** — so
/// same-composition isomers (bisecting vs antenna, core vs Lewis fucose, …)
/// stay distinct.
pub fn load_glycan_gdb(content: &str) -> Result<Vec<GlycanComp>, GdbLoadError> {
    let mut lines = content.lines();
    let _header = lines.next().ok_or(GdbLoadError::Empty)?;
    let mut out: Vec<GlycanComp> = Vec::new();
    let mut seen: std::collections::HashSet<Arc<GlycanTree>> = std::collections::HashSet::new();
    for (i, raw) in lines.enumerate() {
        let line = raw.trim();
        if line.is_empty() {
            continue;
        }
        let line_no = i + 2; // header is line 1
        let bytes = line.as_bytes();
        let mut pos = 0usize;
        let comp = parse_glycan_at(bytes, &mut pos, line_no)?;
        if pos < bytes.len() {
            if bytes[pos] != b'(' {
                return Err(GdbLoadError::MalformedTree {
                    line: line_no,
                    reason: "trailing characters after the glycan tree",
                });
            }
            let mut trees = 1usize;
            while pos < bytes.len() {
                parse_glycan_at(bytes, &mut pos, line_no)?;
                trees += 1;
            }
            return Err(GdbLoadError::MultiSite {
                line: line_no,
                trees,
            });
        }
        // Dedup by topology (the tree), not by composition — same-composition
        // isomers are distinct glycans.
        if seen.insert(comp.tree.clone()) {
            out.push(comp);
        }
    }
    finish_glycan_list(out)
}

/// Sort a loaded glycan list by mass ascending (tiebroken by composition, then
/// `core_fuc`, then `topology_id`) and run the topology-id integrity guards.
/// Shared by the `.gdb` and WURCS loaders so both produce the same order and the
/// same loud failure on a topology-id collision.
fn finish_glycan_list(mut out: Vec<GlycanComp>) -> Result<Vec<GlycanComp>, GdbLoadError> {
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
            .then(a.topology_id.cmp(&b.topology_id))
    });
    if out.is_empty() {
        return Err(GdbLoadError::NoGlycans);
    }
    // A topology-id collision would silently merge two isomers at winner
    // collapse; fail loudly instead. (FNV-1a over ~10⁴ distinct trees has a
    // ~10⁻¹¹ collision probability, so this is a guard, not a real risk.)
    debug_assert!(
        out.iter().all(|g| g.topology_id != DENOVO_TOPOLOGY_ID),
        "glycan topology_id collides with the de-novo sentinel"
    );
    let mut ids: std::collections::HashSet<u64> = std::collections::HashSet::new();
    debug_assert!(
        out.iter().all(|g| ids.insert(g.topology_id)),
        "two distinct glycan topologies share a topology_id (FNV-1a collision)"
    );
    Ok(out)
}

/// Reduce a WURCS 2.0 residue code to an andes monosaccharide symbol
/// (`H`=Hex, `N`=HexNAc, `F`=Fuc, `A`=NeuAc, `G`=NeuGc), or `None` when the
/// residue is outside the supported N-glycan set (e.g. a sulfate, uronic acid
/// or pentose).
///
/// The discriminators are the ones that distinguish the five N-glycan residue
/// classes in a `glycanDatabase.csv` WURCS database:
/// - an uppercase `A` prefix marks a nonulosonic acid (sialic); the `*NCCO`
///   N-glycolyl marker separates NeuGc (`G`) from NeuAc (`A`);
/// - `*NCC` marks an N-acetyl hexosamine (`N`, GlcNAc / GalNAc);
/// - the `m` (6-deoxy) marker in the sugar-type token marks fucose (`F`);
/// - the `h` (hexose) marker marks hexose (`H`);
/// - any *other* `*` substituent (sulfate, phosphate, …) is unsupported.
///
/// The `A`-prefix check comes first because the sialic codes also carry `*NCC` /
/// `*NCCO`; the `*NCC` check comes before the sugar-type fall-through because a
/// HexNAc code is `a…h…_2*NCC`, not a plain hexose.
fn classify_wurcs_residue(code: &str) -> Option<u8> {
    if code.starts_with('A') {
        return Some(if code.contains("*NCCO") { b'G' } else { b'A' });
    }
    if code.contains("*NCC") {
        return Some(b'N');
    }
    // Any other substituent marker (sulfate `*OSO`, phosphate, …) is outside the
    // five-symbol N-glycan set and must not be silently dropped onto a hexose.
    if code.contains('*') {
        return None;
    }
    let sugar = code.split('-').next().unwrap_or(code);
    if sugar.contains('m') {
        return Some(b'F');
    }
    if sugar.contains('h') {
        return Some(b'H');
    }
    None
}

/// Parse one WURCS linkage `Xp-Yq[|Yr…][}]` into its FIRST alternative: `(atom_x,
/// pos_x, atom_y, pos_y)` with each `pos` a raw `?`/digit substring. The `|`- and
/// `}`-separated alternatives are dropped because only the first contributes to
/// the tree topology (they encode ambiguous linkage positions, not extra bonds).
fn parse_wurcs_link(link: &str) -> Option<(usize, &str, usize, &str)> {
    let b = link.as_bytes();
    let mut i = 0usize;
    let xi = wurcs_atom_index(*b.get(i)?)?;
    i += 1;
    let (xp, ni) = parse_wurcs_pos(link, i)?;
    i = ni;
    if *b.get(i)? != b'-' {
        return None;
    }
    i += 1;
    let yi = wurcs_atom_index(*b.get(i)?)?;
    i += 1;
    let (yq, _ni) = parse_wurcs_pos(link, i)?;
    Some((xi, xp, yi, yq))
}

/// Decode a WURCS atom letter: `a`–`z` are backbone atoms 0–25 and `A`–`Z`
/// continue at 26–51 (structures with more than 26 atoms).
fn wurcs_atom_index(c: u8) -> Option<usize> {
    match c {
        b'a'..=b'z' => Some((c - b'a') as usize),
        b'A'..=b'Z' => Some((c - b'A') as usize + 26),
        _ => None,
    }
}

/// Parse a linkage carbon position: `?` (unknown) or a run of digits. Stops at any
/// other byte (`-`, `|`, `}`, …). Returns `(position, new_index)`.
fn parse_wurcs_pos(s: &str, mut i: usize) -> Option<(&str, usize)> {
    let b = s.as_bytes();
    if b.get(i) == Some(&b'?') {
        return Some(("?", i + 1));
    }
    let start = i;
    while let Some(&c) = b.get(i) {
        if c.is_ascii_digit() {
            i += 1;
        } else {
            break;
        }
    }
    if i == start {
        return None;
    }
    Some((&s[start..i], i))
}

/// Does a WURCS linkage position (`n`, or `?`) name the anomeric carbon `anomeric`?
fn wurcs_pos_is_anomeric(pos: &str, anomeric: u32) -> bool {
    pos.parse::<u32>().ok() == Some(anomeric)
}

/// Parse one WURCS 2.0 structure string into a rooted glycan tree. See
/// [`GlycanTree::from_wurcs`] for the reduction rules; `line` is the 1-based
/// source line used in error messages.
fn wurcs_to_tree(wurcs: &str, line: usize) -> Result<GlycanTree, GdbLoadError> {
    let rest = wurcs
        .trim()
        .strip_prefix("WURCS=2.0/")
        .ok_or(GdbLoadError::MalformedWurcs {
            line,
            reason: "missing WURCS=2.0/ prefix",
        })?;

    // WURCS=2.0/<nres>,<natom>,<nlink>[+]/[<residue codes>]/<backbone>/<linkages>.
    // The counts field is not needed: the structure is recovered from the residue
    // block, the backbone and the linkages. Residue codes may contain '/' (their
    // substituent block, e.g. `…*NCC/3=O`), so the residue block is delimited by
    // brackets, not by the next '/'.
    let (_, after_counts) = rest.split_once('/').ok_or(GdbLoadError::MalformedWurcs {
        line,
        reason: "missing '/' after the counts field",
    })?;
    let open = after_counts.find('[').ok_or(GdbLoadError::MalformedWurcs {
        line,
        reason: "missing '[' before the residue block",
    })?;
    // The residue block is `[code][code]…[code]`: its closing bracket is the LAST
    // `]` before the backbone (residue codes themselves contain no brackets).
    let close = after_counts[open..]
        .rfind(']')
        .ok_or(GdbLoadError::MalformedWurcs {
            line,
            reason: "missing ']' after the residue block",
        })?;
    let close = open + close;
    let residue_block = &after_counts[open + 1..close];
    let after_res =
        after_counts[close + 1..]
            .strip_prefix('/')
            .ok_or(GdbLoadError::MalformedWurcs {
                line,
                reason: "missing '/' after the residue block",
            })?;
    let (backbone, linkages) = after_res
        .split_once('/')
        .ok_or(GdbLoadError::MalformedWurcs {
            line,
            reason: "missing '/' before the linkage block",
        })?;

    // Residue codes -> symbols.
    let residues: Vec<u8> = residue_block
        .split("][")
        .map(|code| {
            classify_wurcs_residue(code).ok_or(GdbLoadError::UnknownWurcsResidue {
                line,
                code: code.to_string(),
            })
        })
        .collect::<Result<_, _>>()?;

    // Backbone: one 1-based residue-type index per atom (a, b, c, …).
    let mut atom_res: Vec<usize> = Vec::new();
    for tok in backbone.split('-') {
        let idx = tok
            .parse::<usize>()
            .map_err(|_| GdbLoadError::MalformedWurcs {
                line,
                reason: "non-numeric residue index in the backbone",
            })?;
        let idx = idx.checked_sub(1).ok_or(GdbLoadError::MalformedWurcs {
            line,
            reason: "zero residue index in the backbone",
        })?;
        if idx >= residues.len() {
            return Err(GdbLoadError::MalformedWurcs {
                line,
                reason: "backbone residue index out of range",
            });
        }
        atom_res.push(idx);
    }
    let n = atom_res.len();

    // Classify a linkage carbon for its atom: `1` = anomeric (the non-reducing
    // end / donor, C1 for aldoses, C2 for sialic acid), `0` = a known non-anomeric
    // carbon (the reducing end), `-1` = `?` unknown. A glycosidic bond joins the
    // child's anomeric carbon to the parent's non-anomeric carbon, so a lone
    // anomeric carbon decides the child, and a lone non-anomeric carbon (with the
    // other side `?`) decides the parent.
    let sym_of = |ai: usize| residues[atom_res[ai]];
    let role = |ai: usize, pos: &str| -> i8 {
        if pos == "?" {
            return -1;
        }
        let anomeric = if sym_of(ai) == b'A' || sym_of(ai) == b'G' {
            2
        } else {
            1
        };
        if wurcs_pos_is_anomeric(pos, anomeric) {
            1
        } else {
            0
        }
    };

    // Parent→child edges, where the child is the non-reducing end.
    let mut edges: Vec<(usize, usize)> = Vec::new();
    for link in linkages.split('_') {
        if link.is_empty() {
            continue;
        }
        let (xi, xp, yi, yq) = parse_wurcs_link(link).ok_or(GdbLoadError::MalformedWurcs {
            line,
            reason: "malformed linkage",
        })?;
        if xi >= n || yi >= n {
            return Err(GdbLoadError::MalformedWurcs {
                line,
                reason: "linkage references an atom beyond the backbone",
            });
        }
        let rx = role(xi, xp);
        let ry = role(yi, yq);
        let sx = sym_of(xi) == b'A' || sym_of(xi) == b'G';
        let sy = sym_of(yi) == b'A' || sym_of(yi) == b'G';
        if rx == 1 && ry != 1 {
            edges.push((yi, xi)); // x is the donor (child), y the parent
        } else if ry == 1 && rx != 1 {
            edges.push((xi, yi)); // y is the donor (child), x the parent
        } else if sx && !sy {
            edges.push((yi, xi)); // x is sialic → terminal (child) in an N-glycan
        } else if sy && !sx {
            edges.push((xi, yi)); // y is sialic → terminal (child)
        } else if rx == 0 && ry == -1 {
            edges.push((xi, yi)); // x at a non-anomeric carbon → x is the parent
        } else if ry == 0 && rx == -1 {
            edges.push((yi, xi)); // y at a non-anomeric carbon → y is the parent
        } else if xi < yi {
            // Both `?` (or both on the same kind of carbon): fall back to backbone
            // order — the lower-index atom is closer to the reducing end.
            edges.push((xi, yi));
        } else if yi < xi {
            edges.push((yi, xi));
        }
    }

    // Build the parent→children adjacency and an undirected adjacency for
    // component-finding. The real glycan is the connected component containing
    // atom 0 (the first backbone atom, conventionally the reducing end); WURCS
    // `+` repeating-unit structures add extra placeholder atoms — either isolated
    // or forming a separate component — which are ignored by rooting there.
    let mut children: Vec<Vec<usize>> = vec![Vec::new(); n];
    let mut has_parent = vec![false; n];
    let mut adj: Vec<Vec<usize>> = vec![Vec::new(); n];
    for &(p, c) in &edges {
        if p != c {
            children[p].push(c);
            has_parent[c] = true;
            adj[p].push(c);
            adj[c].push(p);
        }
    }
    if n == 1 {
        // A single-residue glycan (e.g. bare HexNAc, `nlink=0`): no linkages, so
        // the lone atom is its own (leaf) root.
        return Ok(GlycanTree {
            sym: residues[atom_res[0]],
            children: Vec::new(),
        });
    }
    // Connected component of atom 0.
    let mut in_main = vec![false; n];
    let mut stack = vec![0usize];
    in_main[0] = true;
    while let Some(u) = stack.pop() {
        for &v in &adj[u] {
            if !in_main[v] {
                in_main[v] = true;
                stack.push(v);
            }
        }
    }
    // Root the tree at the single atom of the main component with children but no
    // parent (the reducing end).
    let mut roots = (0..n).filter(|&i| in_main[i] && !children[i].is_empty() && !has_parent[i]);
    let root = roots.next().ok_or(GdbLoadError::MalformedWurcs {
        line,
        reason: "no root residue (no connected acyclic backbone)",
    })?;
    if roots.next().is_some() {
        return Err(GdbLoadError::MalformedWurcs {
            line,
            reason: "multiple roots (cyclic or disconnected backbone)",
        });
    }

    build_wurcs_node(root, &children, &residues, &atom_res).check_counts(line)
}

/// Recursively build a `GlycanTree` from a WURCS parent→children adjacency,
/// normalising child order so topology identity is order-independent.
fn build_wurcs_node(
    i: usize,
    children: &[Vec<usize>],
    residues: &[u8],
    atom_res: &[usize],
) -> GlycanTree {
    let sym = residues[atom_res[i]];
    let mut kids: Vec<GlycanTree> = children[i]
        .iter()
        .map(|&c| build_wurcs_node(c, children, residues, atom_res))
        .collect();
    kids.sort_by_cached_key(|k| k.canonical());
    GlycanTree {
        sym,
        children: kids,
    }
}

/// Split one CSV line into fields, honouring double-quoted fields (RFC 4180
/// quotes, with `""` as an escaped quote). Embedded newlines are not handled
/// because callers feed one line at a time and the glycan database has none.
fn split_csv_line(line: &str) -> Vec<String> {
    let mut fields = Vec::new();
    let mut cur = String::new();
    let mut in_quotes = false;
    let mut chars = line.chars().peekable();
    while let Some(c) = chars.next() {
        if in_quotes {
            if c == '"' {
                if chars.peek() == Some(&'"') {
                    chars.next();
                    cur.push('"');
                } else {
                    in_quotes = false;
                }
            } else {
                cur.push(c);
            }
        } else if c == '"' {
            in_quotes = true;
        } else if c == ',' {
            fields.push(std::mem::take(&mut cur));
        } else {
            cur.push(c);
        }
    }
    fields.push(cur);
    fields
}

/// Load a WURCS 2.0 glycan database CSV (a `glycanDatabase.csv`) into glycans,
/// retaining the full tree topology.
///
/// Columns: `GlycanID,Mass,GlycanComposition,StructureInformation`.
/// `StructureInformation` is a quoted WURCS 2.0 string parsed into the same
/// canonical tree as a `.gdb` line; the `Mass` and `GlycanComposition` columns
/// are ignored (both are derived from the tree, so a rounded or stale mass in the
/// CSV cannot desynchronise the search space). Deduplicated by topology and
/// sorted identically to [`load_glycan_gdb`].
pub fn load_glycan_wurcs_csv(content: &str) -> Result<Vec<GlycanComp>, GdbLoadError> {
    let mut lines = content.lines();
    let header = lines.next().ok_or(GdbLoadError::Empty)?;
    let header_fields = split_csv_line(header);
    let wurcs_col = header_fields
        .iter()
        .position(|h| h.trim() == "StructureInformation")
        .unwrap_or(header_fields.len().saturating_sub(1));
    let mut out: Vec<GlycanComp> = Vec::new();
    let mut seen: std::collections::HashSet<Arc<GlycanTree>> = std::collections::HashSet::new();
    for (i, raw) in lines.enumerate() {
        let line = raw.trim();
        if line.is_empty() {
            continue;
        }
        let line_no = i + 2; // header is line 1
        let fields = split_csv_line(line);
        let wurcs = fields.get(wurcs_col).ok_or(GdbLoadError::MalformedWurcs {
            line: line_no,
            reason: "missing StructureInformation column",
        })?;
        let tree = wurcs_to_tree(wurcs, line_no)?;
        let comp = GlycanComp::from_tree(tree);
        if seen.insert(comp.tree.clone()) {
            out.push(comp);
        }
    }
    finish_glycan_list(out)
}

/// Load one of the bundled species-specific N-glycan databases by name.
///
/// The databases are compiled in via `include_str!` (see `glycan-db/`); `species`
/// is the kebab-case name exposed by the `--glyco-species` flag (`human`,
/// `mouse`, `mouse-large`, `high-mannose`). This is the `--glyco-species` path of the glycan-first
/// search; an explicit `--glyco-glycan-gdb` file bypasses it entirely.
pub fn load_species_glycan_db(species: &str) -> Result<Vec<GlycanComp>, GdbLoadError> {
    let content = match species {
        "human" => include_str!("../glycan-db/N-glycan-human.gdb"),
        "mouse" => include_str!("../glycan-db/N-glycan-mouse.gdb"),
        "mouse-large" => include_str!("../glycan-db/N-glycan-mouse-large.gdb"),
        "high-mannose" => include_str!("../glycan-db/N-glycan-high-mannose.gdb"),
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

    #[test]
    fn a_residue_count_past_255_is_rejected() {
        let wide = |n: usize| format!("(N{})", "(H)".repeat(n));
        let ok = GlycanTree::parse(&wide(255)).unwrap();
        assert_eq!(ok.counts()[0], 255);
        assert_eq!(ok.checked_counts().map(|c| c[0]), Some(255));
        for n in [256, 300] {
            let err = GlycanTree::parse(&wide(n)).unwrap_err();
            assert!(
                matches!(err, GdbLoadError::MalformedTree { line: 1, .. }),
                "{err:?}"
            );
            let gdb = format!("H,N,A,G,F\n(N(F))\n{}\n", wide(n));
            let err = load_glycan_gdb(&gdb).unwrap_err();
            assert!(
                matches!(err, GdbLoadError::MalformedTree { line: 3, .. }),
                "{err:?}"
            );
        }
        // The count also holds across branches.
        let split = format!("(N(N{})(N{}))", "(H)".repeat(128), "(H)".repeat(128));
        assert!(GlycanTree::parse(&split).is_err());
    }

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
    fn multi_site_line_is_rejected_not_split() {
        // Upstream `*-multi` lines concatenate one complete tree per glycosite (each
        // with its own reducing-end GlcNAc). Splitting them would add each half as
        // a single glycan; andes places one glycan per peptide, so refuse.
        let content = "H,N,A,G,F\n(N(N(H(H)(H))))\n(N(N(H(H)(H))))(N(N(H(H(H))(H(H)))))\n";
        let err = load_glycan_gdb(content).unwrap_err();
        assert_eq!(err, GdbLoadError::MultiSite { line: 3, trees: 2 });
        assert!(err.to_string().contains("multi-site"));
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
        for sp in ["human", "mouse", "mouse-large", "high-mannose"] {
            let list = load_species_glycan_db(sp).unwrap();
            assert!(!list.is_empty(), "{sp}: empty glycan list");
        }
        // The mouse and human databases must deduplicate by TOPOLOGY, keeping
        // thousands of distinct trees (human ~2920 trees vs ~1004 compositions;
        // mouse ~6660 vs ~1833). A regression to composition dedup would collapse
        // far below these lower bounds.
        assert!(load_species_glycan_db("mouse-large").unwrap().len() >= 6500);
        assert!(load_species_glycan_db("mouse").unwrap().len() >= 5500);
        assert!(load_species_glycan_db("human").unwrap().len() >= 2500);
        assert!(load_species_glycan_db("high-mannose").unwrap().len() >= 10);
    }

    #[test]
    fn unknown_species_is_rejected() {
        let err = load_species_glycan_db("zebra").unwrap_err();
        assert!(matches!(err, GdbLoadError::UnknownSpecies { .. }));
        let err = load_species_glycan_db("human-multi").unwrap_err();
        assert!(matches!(err, GdbLoadError::UnknownSpecies { .. }));
    }

    #[test]
    fn topology_isomers_are_kept_distinct() {
        // Same composition (N2 H3 F1) but core-fucose vs antenna-fucose: two
        // distinct topologies that must not collapse under tree dedup.
        let content = "H,N,A,G,F\n(N(F)(N(H(H)(H))))\n(N(N(H(H)(H)(F))))\n";
        let list = load_glycan_gdb(content).unwrap();
        assert_eq!(list.len(), 2);
        assert_ne!(list[0].topology_id, list[1].topology_id);
        assert_ne!(list[0].tree.canonical(), list[1].tree.canonical());
        assert_eq!(
            (list[0].hexnac, list[0].hex, list[0].fuc),
            (list[1].hexnac, list[1].hex, list[1].fuc)
        );
    }

    #[test]
    fn y_ion_retained_masses_is_single_cleavage_set() {
        // The bare trimannosyl core yields one ion per glycosidic bond (plus the
        // intact glycan): cutting root→GlcNAc keeps GlcNAc (Y1), cutting
        // GlcNAc→Man keeps chitobiose (Y2), cutting either Man→Man keeps
        // chitobiose+2Man. The chitobiose+1Man rung (2N+H) is a TWO-cleavage
        // fragment, so it is NOT in the single-cleavage set.
        let g = load_one("(N(N(H(H)(H))))");
        let y = g.y_ion_retained_masses();
        assert_eq!(y.len(), 4);
        let expect = [
            HEXNAC,
            2.0 * HEXNAC,
            2.0 * HEXNAC + 2.0 * HEX,
            2.0 * HEXNAC + 3.0 * HEX,
        ];
        for (a, b) in y.iter().zip(expect) {
            assert!((a - b).abs() < 1e-6, "{a} != {b}");
        }
    }

    #[test]
    fn y_ion_retained_masses_bisecting_vs_antenna_differ() {
        // Bisecting GlcNAc on the beta-Man vs GlcNAc on an antenna alpha-Man:
        // same composition (N3 H3), different topology, different Y-ion masses.
        let bisecting = GlycanComp::from_tree_str("(N(N(H(N)(H)(H))))").unwrap();
        let antenna = GlycanComp::from_tree_str("(N(N(H(H(N))(H))))").unwrap();
        assert_eq!(
            (bisecting.hexnac, bisecting.hex),
            (antenna.hexnac, antenna.hex)
        );
        assert_ne!(bisecting.topology_id, antenna.topology_id);
        assert_ne!(
            bisecting.y_ion_retained_masses(),
            antenna.y_ion_retained_masses()
        );
    }

    #[test]
    fn canonical_is_order_independent() {
        // Children are normalised at parse, so a hand-shuffled equivalent tree
        // still compares/hashes/canonicalises identically.
        let a = GlycanTree::parse("(N(H(H))(H))").unwrap();
        let b = GlycanTree::parse("(N(H)(H(H)))").unwrap();
        assert_eq!(a, b);
        assert_eq!(a.canonical(), b.canonical());
        assert_eq!(a.topology_id(), b.topology_id());
    }

    #[test]
    fn topology_id_is_deterministic_and_distinct() {
        let a = GlycanTree::parse("(N(N(H(H)(H))))").unwrap();
        let b = GlycanTree::parse("(N(F)(N(H(H)(H))))").unwrap();
        assert_eq!(a.topology_id(), a.topology_id());
        assert_ne!(a.topology_id(), b.topology_id());
        assert_ne!(a.topology_id(), DENOVO_TOPOLOGY_ID);
    }

    fn wurcs(w: &str) -> GlycanTree {
        GlycanTree::from_wurcs(w).unwrap()
    }

    #[test]
    fn wurcs_trimannosyl_core_matches_gdb() {
        // The WURCS trimannosyl core (HexNAc2 Hex3) must reduce to the exact
        // tree the `.gdb` loader builds: reducing-end GlcNAc → chitobiose GlcNAc →
        // the three mannose residues.
        let t = wurcs(
            "WURCS=2.0/3,5,4/[a2122h-1b_1-5_2*NCC/3=O][a1122h-1b_1-5][a1122h-1a_1-5]\
             /1-1-2-3-3/a4-b1_b4-c1_c3-d1_c6-e1",
        );
        assert_eq!(t.canonical(), "(N(N(H(H)(H))))");
        let g = GlycanComp::from_tree(t);
        assert_eq!((g.hexnac, g.hex, g.fuc, g.core_fuc), (2, 3, 0, 0));
    }

    #[test]
    fn wurcs_core_fucose_is_direct_child_of_root() {
        // HexNAc + Fuc, linked through the fucose anomeric carbon: Fuc sits
        // directly on the reducing-end HexNAc, so it is a core fucose.
        let t = wurcs("WURCS=2.0/2,2,1/[a2122h-1x_1-5_2*NCC/3=O][a1221m-1x_1-5]/1-2/a?-b1");
        assert_eq!(t.canonical(), "(N(F))");
        let g = GlycanComp::from_tree(t);
        assert_eq!((g.hexnac, g.hex, g.fuc, g.core_fuc), (1, 0, 1, 1));
    }

    #[test]
    fn wurcs_sialic_acid_is_classified() {
        // NeuAc and NeuGc differ only in the N-glycolyl marker `*NCCO`; the A/G
        // split must not collapse them.
        let neuac = wurcs(
            "WURCS=2.0/5,8,7/[a2122h-1b_1-5_2*NCC/3=O][a1122h-1b_1-5][a1122h-1a_1-5]\
             [a2112h-1b_1-5][Aad21122h-2a_2-6_5*NCC/3=O]/1-1-2-3-1-4-5-3\
             /a4-b1_b4-c1_c3-d1_c6-h1_d2-e1_e4-f1_f3-g2",
        );
        let g = GlycanComp::from_tree(neuac);
        assert_eq!((g.neuac, g.neugc), (1, 0));
        let neugc = wurcs(
            "WURCS=2.0/5,8,7/[a2122h-1b_1-5_2*NCC/3=O][a1122h-1b_1-5][a1122h-1a_1-5]\
             [a2112h-1b_1-5][Aad21122h-2a_2-6_5*NCCO/3=O]/1-1-2-3-1-4-5-3\
             /a4-b1_b4-c1_c3-d1_c6-h1_d2-e1_e4-f1_f3-g2",
        );
        let g = GlycanComp::from_tree(neugc);
        assert_eq!((g.neuac, g.neugc), (0, 1));
    }

    #[test]
    fn wurcs_repeat_placeholder_atoms_are_ignored() {
        // Hex(6)HexNAc(2) carries a `+` repeat marker: two isolated backbone atoms
        // (the repeating-unit placeholders) must be ignored, leaving exactly
        // Hex6 HexNAc2.
        let t = wurcs(
            "WURCS=2.0/3,10,7+/[a2122h-1b_1-5_2*NCC/3=O][a1122h-1b_1-5][a1122h-1a_1-5]\
             /1-1-2-3-3-3-3-3-3-3/a4-b1_b4-c1_c3-d1_c6-f1_d2-e1_f3-g1_f6-h1",
        );
        let g = GlycanComp::from_tree(t);
        assert_eq!((g.hexnac, g.hex), (2, 6));
    }

    #[test]
    fn wurcs_alternative_linkage_uses_first_position() {
        // `b1-a3|a4` / `c1-b3|b4` encode ambiguous (alternative) linkage positions;
        // the first alternative decides the topology, yielding the same
        // trimannosyl core as the unambiguous encoding.
        let t = wurcs(
            "WURCS=2.0/4,5,4/[u2122h_2*NCC/3=O][a2122h-1b_1-5_2*NCC/3=O][a1122h-1b_1-5]\
             [a1122h-1a_1-5]/1-2-3-4-4/c3-d1_c6-e1_b1-a3|a4_c1-b3|b4",
        );
        assert_eq!(t.canonical(), "(N(N(H(H)(H))))");
    }

    #[test]
    fn wurcs_repeat_disconnected_component_is_ignored() {
        // The `+` repeat marker adds a second, disconnected component (`f4-g1`) to
        // the backbone; only the component containing the reducing end (atoms a–e)
        // is the real glycan, so the placeholder must be dropped, leaving the
        // trimannosyl core.
        let t = wurcs(
            "WURCS=2.0/4,7,5+/[a2122h-1b_1-5_2*NCC/3=O][a1122h-1b_1-5][a1122h-1a_1-5]\
             [a2112h-1b_1-5]/1-1-2-3-3-1-4/a4-b1_b4-c1_c3-d1_c6-e1_f4-g1",
        );
        assert_eq!(t.canonical(), "(N(N(H(H)(H))))");
        let g = GlycanComp::from_tree(t);
        assert_eq!((g.hexnac, g.hex), (2, 3));
    }

    #[test]
    fn wurcs_fuzzy_linkage_keeps_the_residue() {
        // `f?-a?|b?|…}` has an unknown linkage position on BOTH sides; the extra
        // Hex `f` must still be kept (attached to the lower-index reducing end)
        // rather than silently dropped as an isolated atom.
        let t = wurcs(
            "WURCS=2.0/4,6,5/[a2122h-1b_1-5_2*NCC/3=O][a1122h-1b_1-5][a1122h-1a_1-5]\
             [a1122h-1x_1-5]/1-1-2-3-3-4/a4-b1_b4-c1_c3-d1_c6-e1_f?-a?|b?|c?|d?|e?}",
        );
        let g = GlycanComp::from_tree(t);
        assert_eq!((g.hexnac, g.hex), (2, 4));
    }

    #[test]
    fn wurcs_sialic_at_position_one_is_terminal() {
        // A sialic acid fuzzily linked at position 1 (`k1-f?|i?}`) is still the
        // terminal child: sialic caps an antenna, so it must not be made a parent
        // just because position 1 is non-anomeric for it.
        let t = wurcs(
            "WURCS=2.0/7,11,10/[a2122h-1x_1-5_2*NCC/3=O][a2122h-1b_1-5_2*NCC/3=O]\
             [a1122h-1b_1-5][a1122h-1a_1-5][a2112h-1b_1-5][a1221m-1a_1-5]\
             [Aad21122h-2a_2-6_5*NCC/3=O]\
             /1-2-3-4-2-5-4-2-5-6-7/a4-b1_a6-j1_b4-c1_c3-d1_c6-g1_d2-e1_e4-f1_g2-h1_h4-i1_k1-f?|i?}",
        );
        let g = GlycanComp::from_tree(t);
        assert_eq!((g.hexnac, g.hex, g.fuc, g.neuac), (4, 5, 1, 1));
    }

    #[test]
    fn wurcs_atom_letters_continue_past_z() {
        // >26 backbone atoms use uppercase letters (`A` = atom 26, `D` = 29), so
        // linkages like `a6-D1`, `z4-A1`, … must parse and yield the full
        // Hex(15) HexNAc(14) Fuc(1) composition.
        let t = wurcs(
            "WURCS=2.0/5,30,29/[a2122h-1b_1-5_2*NCC/3=O][a1122h-1b_1-5][a1122h-1a_1-5]\
             [a2112h-1b_1-5][a1221m-1a_1-5]\
             /1-1-2-3-1-4-1-4-1-4-1-4-3-1-4-1-4-1-4-1-4-1-4-1-4-1-4-1-4-5\
             /a4-b1_a6-D1_b4-c1_c3-d1_c6-m1_d2-e1_e4-f1_f3-g1_g4-h1_h3-i1_i4-j1_j3-k1_k4-l1\
             _m2-n1_m6-v1_n4-o1_o3-p1_p4-q1_q3-r1_r4-s1_s3-t1_t4-u1_v4-w1_w3-x1_x4-y1_y3-z1\
             _z4-A1_A3-B1_B4-C1",
        );
        let g = GlycanComp::from_tree(t);
        assert_eq!((g.hexnac, g.hex, g.fuc), (14, 15, 1));
    }

    #[test]
    fn wurcs_unknown_residue_is_rejected() {
        // A sulfate residue (no A prefix, no *NCC, no hexose/deoxy marker) must be
        // rejected, not silently mapped to a hexose.
        let err = GlycanTree::from_wurcs(
            "WURCS=2.0/1,1,0/[a2122h-1x_1-5_2*NCC/3=O][s2122h-1x_1-5_6*OSO/3=O]/1-2/a?-b1",
        )
        .unwrap_err();
        assert!(matches!(err, GdbLoadError::UnknownWurcsResidue { .. }));
    }

    #[test]
    fn wurcs_malformed_is_rejected() {
        assert!(GlycanTree::from_wurcs("not a wurcs").is_err());
        assert!(GlycanTree::from_wurcs("WURCS=2.0/1,1,0/[a2122h-1x_1-5_2*NCC/3=O]/1/").is_ok());
    }

    #[test]
    fn wurcs_csv_loads_and_derives_from_tree() {
        // The Mass / GlycanComposition columns are ignored: the loader must derive
        // both from the WURCS tree, and it must handle the quoted WURCS field
        // (which itself contains commas in the counts).
        let csv = "\
GlycanID,Mass,GlycanComposition,StructureInformation
G64581RP,999.0,HexNAc(1),\"WURCS=2.0/1,1,0/[a2122h-1x_1-5_2*NCC/3=O]/1/\"
G39631KX,999.0,HexNAc(1)Fuc(1),\"WURCS=2.0/2,2,1/[a2122h-1x_1-5_2*NCC/3=O][a1221m-1x_1-5]/1-2/a?-b1\"
";
        let list = load_glycan_wurcs_csv(csv).unwrap();
        assert_eq!(list.len(), 2);
        // Sorted by mass ascending: HexNAc before HexNAc+Fuc.
        assert_eq!((list[0].hexnac, list[0].fuc), (1, 0));
        assert_eq!((list[1].hexnac, list[1].fuc), (1, 1));
        assert!((list[0].mass - HEXNAC).abs() < 1e-6);
        assert!((list[1].mass - (HEXNAC + FUC)).abs() < 1e-6);
    }

    #[test]
    fn wurcs_csv_duplicate_topologies_are_deduplicated() {
        // Two rows with different GlycanIDs but the same tree reduce to one glycan.
        let csv = "\
GlycanID,Mass,GlycanComposition,StructureInformation
GAAA0001,203.0794,HexNAc(1),\"WURCS=2.0/1,1,0/[a2122h-1x_1-5_2*NCC/3=O]/1/\"
GAAA0002,203.0794,HexNAc(1),\"WURCS=2.0/1,1,0/[a2122h-1x_1-5_2*NCC/3=O]/1/\"
";
        assert_eq!(load_glycan_wurcs_csv(csv).unwrap().len(), 1);
    }
}
