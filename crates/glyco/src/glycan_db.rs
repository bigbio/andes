use std::sync::Arc;

use model::hash::fnv1a64;

use crate::glycan_mass::{FUC, HEX, HEXNAC, NEUAC, NEUGC};

pub const DENOVO_TOPOLOGY_ID: u64 = u64::MAX;

#[derive(Debug, Clone)]
pub struct GlycanComp {
    pub hexnac: u8,
    pub hex: u8,
    pub fuc: u8,
    pub neuac: u8,
    pub neugc: u8,
    pub core_fuc: u8,
    pub mass: f64,
    pub tree: Arc<GlycanTree>,
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

    pub fn from_tree_str(s: &str) -> Result<GlycanComp, GdbLoadError> {
        Ok(GlycanComp::from_tree(GlycanTree::parse(s)?))
    }

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
        let mut chain: Vec<u8> = Vec::with_capacity((hexnac + hex + fuc + neuac + neugc) as usize);
        chain.extend(std::iter::repeat_n(
            b'N',
            core_hn.saturating_sub(1) as usize,
        ));
        chain.extend(std::iter::repeat_n(b'H', core_h as usize));
        chain.extend(std::iter::repeat_n(b'N', (hexnac - core_hn) as usize));
        chain.extend(std::iter::repeat_n(b'H', (hex - core_h) as usize));
        chain.extend(std::iter::repeat_n(b'F', (fuc - core_fuc) as usize));
        chain.extend(std::iter::repeat_n(b'A', neuac as usize));
        chain.extend(std::iter::repeat_n(b'G', neugc as usize));

        let mut s = String::from("(N");
        if core_fuc > 0 {
            s.push_str("(F)");
        }
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

    pub fn y_ion_retained_masses(&self) -> Vec<f64> {
        self.tree.y_ion_retained_masses()
    }
}

#[derive(Debug, Clone, PartialEq)]
pub enum GdbLoadError {
    Empty,
    NoGlycans,
    UnknownSymbol { line: usize, symbol: String },
    MalformedTree { line: usize, reason: &'static str },
    MultiSite { line: usize, trees: usize },
    MalformedWurcs { line: usize, reason: &'static str },
    UnknownWurcsResidue { line: usize, code: String },
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

#[derive(Debug, Clone, PartialEq, Eq, Hash)]
pub struct GlycanTree {
    pub sym: u8,
    pub children: Vec<GlycanTree>,
}

impl GlycanTree {
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

    pub fn from_wurcs(wurcs: &str) -> Result<GlycanTree, GdbLoadError> {
        wurcs_to_tree(wurcs, 1)
    }

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

    pub fn mass(&self) -> f64 {
        let mut m = GlycanTree::residue_mass(self.sym);
        for c in &self.children {
            m += c.mass();
        }
        m
    }

    pub fn counts(&self) -> [u8; 5] {
        let mut counts = [0u8; 5];
        tally(self, &mut counts);
        counts
    }

    pub fn checked_counts(&self) -> Option<[u8; 5]> {
        let mut counts = [0u8; 5];
        tally(self, &mut counts).then_some(counts)
    }

    fn check_counts(self, line: usize) -> Result<GlycanTree, GdbLoadError> {
        match self.checked_counts() {
            Some(_) => Ok(self),
            None => Err(GdbLoadError::MalformedTree {
                line,
                reason: "more than 255 residues of one kind",
            }),
        }
    }

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

    pub fn topology_id(&self) -> u64 {
        fnv1a64(self.canonical().as_bytes())
    }

    pub fn y_ion_retained_masses(&self) -> Vec<f64> {
        let total = self.mass();
        let mut v: Vec<f64> = vec![total];
        single_cleavage_y_masses(self, total, &mut v);
        v.sort_by(|a, b| a.total_cmp(b));
        let mut out: Vec<f64> = Vec::with_capacity(v.len());
        for m in v {
            match out.last() {
                Some(&last) if (m - last) <= 1e-6 => {}
                _ => out.push(m),
            }
        }
        out
    }
}

fn single_cleavage_y_masses(node: &GlycanTree, total: f64, out: &mut Vec<f64>) {
    for child in &node.children {
        out.push(total - child.mass());
        single_cleavage_y_masses(child, total, out);
    }
}

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

const MAX_GLYCAN_DEPTH: usize = 64;

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

fn parse_glycan_at(
    bytes: &[u8],
    pos: &mut usize,
    line_no: usize,
) -> Result<GlycanComp, GdbLoadError> {
    let root = parse_node(bytes, pos, line_no, 0)?.check_counts(line_no)?;
    Ok(GlycanComp::from_tree(root))
}

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
        let line_no = i + 2;
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
        if seen.insert(comp.tree.clone()) {
            out.push(comp);
        }
    }
    finish_glycan_list(out)
}

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

fn classify_wurcs_residue(code: &str) -> Option<u8> {
    if code.starts_with('A') {
        return Some(if code.contains("*NCCO") { b'G' } else { b'A' });
    }
    if code.contains("*NCC") {
        return Some(b'N');
    }
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

fn wurcs_atom_index(c: u8) -> Option<usize> {
    match c {
        b'a'..=b'z' => Some((c - b'a') as usize),
        b'A'..=b'Z' => Some((c - b'A') as usize + 26),
        _ => None,
    }
}

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

fn wurcs_pos_is_anomeric(pos: &str, anomeric: u32) -> bool {
    pos.parse::<u32>().ok() == Some(anomeric)
}

fn wurcs_to_tree(wurcs: &str, line: usize) -> Result<GlycanTree, GdbLoadError> {
    let rest = wurcs
        .trim()
        .strip_prefix("WURCS=2.0/")
        .ok_or(GdbLoadError::MalformedWurcs {
            line,
            reason: "missing WURCS=2.0/ prefix",
        })?;

    let (_, after_counts) = rest.split_once('/').ok_or(GdbLoadError::MalformedWurcs {
        line,
        reason: "missing '/' after the counts field",
    })?;
    let open = after_counts.find('[').ok_or(GdbLoadError::MalformedWurcs {
        line,
        reason: "missing '[' before the residue block",
    })?;
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

    let residues: Vec<u8> = residue_block
        .split("][")
        .map(|code| {
            classify_wurcs_residue(code).ok_or(GdbLoadError::UnknownWurcsResidue {
                line,
                code: code.to_string(),
            })
        })
        .collect::<Result<_, _>>()?;

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
            edges.push((yi, xi));
        } else if ry == 1 && rx != 1 {
            edges.push((xi, yi));
        } else if sx && !sy {
            edges.push((yi, xi));
        } else if (sy && !sx) || (rx == 0 && ry == -1) {
            edges.push((xi, yi));
        } else if ry == 0 && rx == -1 {
            edges.push((yi, xi));
        } else if xi < yi {
            edges.push((xi, yi));
        } else if yi < xi {
            edges.push((yi, xi));
        }
    }

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
        return Ok(GlycanTree {
            sym: residues[atom_res[0]],
            children: Vec::new(),
        });
    }
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
        let line_no = i + 2;
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
        let g = load_one("(N(F)(N(H(H)(H))))");
        assert_eq!((g.hexnac, g.hex, g.fuc, g.core_fuc), (2, 3, 1, 1));
    }

    #[test]
    fn antenna_fucose_is_not_core_fucose() {
        let g = load_one("(N(N(H(H)(H(N(F))))))");
        assert_eq!((g.hexnac, g.hex, g.fuc, g.core_fuc), (3, 3, 1, 0));
    }

    #[test]
    fn core_and_antenna_fucose_coexist() {
        let g = load_one("(N(F)(N(H(H)(H(N(F))))))");
        assert_eq!((g.hexnac, g.hex, g.fuc, g.core_fuc), (3, 3, 2, 1));
    }

    #[test]
    fn sialylated_biantennary_core_fucosylated() {
        let g = load_one("(N(F)(N(H(H(N(H(A))))(H(N(H(A)))))))");
        assert_eq!(
            (g.hexnac, g.hex, g.fuc, g.neuac, g.core_fuc),
            (4, 5, 1, 2, 1)
        );
    }

    #[test]
    fn core_fuc_isomers_are_kept_distinct() {
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
        let content = "H,N,A,G,F\n(N(N(H(H)(H))))\n(N(N(H(H)(H))))(N(N(H(H(H))(H(H)))))\n";
        let err = load_glycan_gdb(content).unwrap_err();
        assert_eq!(err, GdbLoadError::MultiSite { line: 3, trees: 2 });
        assert!(err.to_string().contains("multi-site"));
    }

    #[test]
    fn deeply_nested_tree_is_rejected() {
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
        let t = wurcs("WURCS=2.0/2,2,1/[a2122h-1x_1-5_2*NCC/3=O][a1221m-1x_1-5]/1-2/a?-b1");
        assert_eq!(t.canonical(), "(N(F))");
        let g = GlycanComp::from_tree(t);
        assert_eq!((g.hexnac, g.hex, g.fuc, g.core_fuc), (1, 0, 1, 1));
    }

    #[test]
    fn wurcs_sialic_acid_is_classified() {
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
        let t = wurcs(
            "WURCS=2.0/3,10,7+/[a2122h-1b_1-5_2*NCC/3=O][a1122h-1b_1-5][a1122h-1a_1-5]\
             /1-1-2-3-3-3-3-3-3-3/a4-b1_b4-c1_c3-d1_c6-f1_d2-e1_f3-g1_f6-h1",
        );
        let g = GlycanComp::from_tree(t);
        assert_eq!((g.hexnac, g.hex), (2, 6));
    }

    #[test]
    fn wurcs_alternative_linkage_uses_first_position() {
        let t = wurcs(
            "WURCS=2.0/4,5,4/[u2122h_2*NCC/3=O][a2122h-1b_1-5_2*NCC/3=O][a1122h-1b_1-5]\
             [a1122h-1a_1-5]/1-2-3-4-4/c3-d1_c6-e1_b1-a3|a4_c1-b3|b4",
        );
        assert_eq!(t.canonical(), "(N(N(H(H)(H))))");
    }

    #[test]
    fn wurcs_repeat_disconnected_component_is_ignored() {
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
        let t = wurcs(
            "WURCS=2.0/4,6,5/[a2122h-1b_1-5_2*NCC/3=O][a1122h-1b_1-5][a1122h-1a_1-5]\
             [a1122h-1x_1-5]/1-1-2-3-3-4/a4-b1_b4-c1_c3-d1_c6-e1_f?-a?|b?|c?|d?|e?}",
        );
        let g = GlycanComp::from_tree(t);
        assert_eq!((g.hexnac, g.hex), (2, 4));
    }

    #[test]
    fn wurcs_sialic_at_position_one_is_terminal() {
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
        let csv = "\
GlycanID,Mass,GlycanComposition,StructureInformation
G64581RP,999.0,HexNAc(1),\"WURCS=2.0/1,1,0/[a2122h-1x_1-5_2*NCC/3=O]/1/\"
G39631KX,999.0,HexNAc(1)Fuc(1),\"WURCS=2.0/2,2,1/[a2122h-1x_1-5_2*NCC/3=O][a1221m-1x_1-5]/1-2/a?-b1\"
";
        let list = load_glycan_wurcs_csv(csv).unwrap();
        assert_eq!(list.len(), 2);
        assert_eq!((list[0].hexnac, list[0].fuc), (1, 0));
        assert_eq!((list[1].hexnac, list[1].fuc), (1, 1));
        assert!((list[0].mass - HEXNAC).abs() < 1e-6);
        assert!((list[1].mass - (HEXNAC + FUC)).abs() < 1e-6);
    }

    #[test]
    fn wurcs_csv_duplicate_topologies_are_deduplicated() {
        let csv = "\
GlycanID,Mass,GlycanComposition,StructureInformation
GAAA0001,203.0794,HexNAc(1),\"WURCS=2.0/1,1,0/[a2122h-1x_1-5_2*NCC/3=O]/1/\"
GAAA0002,203.0794,HexNAc(1),\"WURCS=2.0/1,1,0/[a2122h-1x_1-5_2*NCC/3=O]/1/\"
";
        assert_eq!(load_glycan_wurcs_csv(csv).unwrap().len(), 1);
    }
}
