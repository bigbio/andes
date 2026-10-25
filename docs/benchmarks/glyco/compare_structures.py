#!/usr/bin/env python3
"""Structure-level (topology) comparison: andes vs engine N, MouseLiver-Z-T-1.

Where `compare_engines.py` collapses every engine to a glycan COMPOSITION
`(HexNAc, Hex, Fuc, NeuAc, NeuGc)`, this one keeps the glycan TREE and asks the
narrower question: for spectra that andes and engine N both identify as the same
peptide, do they agree on the *structure* (branching topology), not just the
residue counts?

Both sides are normalised to the same five-symbol alphabet (H=Hex, N=HexNAc,
F=Fuc, A=NeuAc, G=NeuGc) and the same canonical nested-paren form with the
reducing-end residue as root and children sorted, e.g. `(N(N(H(H)(H))))` for the
trimannosyl core. An andes tree is read straight from its Peptide tag; an engine N
tree is reconstructed from its `GlycanID` via the WURCS 2.0 `StructureInformation`
in the glycan database CSV.

Inputs (paths on the command line):
  * engine N GPSM text export   (PeptideFDR <= 0.01)
  * glycan database CSV         (GlycanID, Mass, GlycanComposition, StructureInformation)
  * andes <stem>.glyco.pin      (top-1-per-scan production PIN)
  * andes native-rescore .psms  (q <= 0.01)
"""
import argparse, csv, re, collections

# ---------- CLI ----------
_ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
_ap.add_argument('--engine-n', required=True, help='engine N GPSM text export')
_ap.add_argument('--glycan-db', required=True, help='glycan database CSV (WURCS)')
_ap.add_argument('--andes-pin', required=True)
_ap.add_argument('--andes-psms', required=True)
_ap.add_argument('--q', type=float, default=0.01)
_args = _ap.parse_args()


# ---------- WURCS 2.0 -> tree ----------
def classify_residue(code):
    if code.startswith('A'):                 # nonulosonic acid (sialic)
        return 'G' if '*NCCO' in code else 'A'
    if '*NCC' in code:                       # N-acetyl hexosamine
        return 'N'
    if 'm' in code.split('-')[0]:            # 6-deoxy (fucose)
        return 'F'
    return 'H'                               # hexose


def wurcs_to_tree(wurcs):
    """WURCS=2.0 backbone form -> canonical nested-paren tree string."""
    m = re.match(r'WURCS=2\.0/(\d+),(\d+),(\d+)\+?/\[(.+)\]/' + r'([^/]+)/(.+)', wurcs)
    if not m:
        return None
    resblock, backbone, linkages = m.group(4), m.group(5), m.group(6)
    residues = [classify_residue(r) for r in resblock.split('][')]
    atom_res = [residues[int(x) - 1] for x in backbone.split('-')]
    n = len(atom_res)

    edges = []                              # (parent_idx, child_idx)
    for link in linkages.split('_'):
        lm = re.match(r'([a-z])(\d+|\?)-([a-z])(\d+|\?)', link)
        if not lm:
            continue
        x, xp, y, yq = lm.group(1), lm.group(2), lm.group(3), lm.group(4)
        xi, yi = ord(x) - 97, ord(y) - 97

        def is_donor(ai, pos):
            # the non-reducing end links through its anomeric carbon: C1 for
            # aldoses, C2 for sialic acid.
            if pos == '?':
                return False
            return int(pos) == (2 if atom_res[ai] in 'AG' else 1)

        dx, dy = is_donor(xi, xp), is_donor(yi, yq)
        if dx and not dy:
            edges.append((yi, xi))
        elif dy and not dx:
            edges.append((xi, yi))
        else:
            if xp == '1' and yq != '1':
                edges.append((yi, xi))
            elif yq == '1' and xp != '1':
                edges.append((xi, yi))

    children = collections.defaultdict(list)
    parent_of = {}
    for p, c in edges:
        children[p].append(c)
        parent_of[c] = p
    roots = [i for i in range(n) if i in children and i not in parent_of]
    if len(roots) != 1:
        return None

    def node_str(i):
        return f"({atom_res[i]}" + ''.join(sorted(node_str(c) for c in children[i])) + ")"
    return node_str(roots[0])


def load_glycan_db(path):
    db = {}
    with open(path) as f:
        r = csv.reader(f)
        next(r)
        for row in r:
            if len(row) < 4:
                continue
            db[row[0]] = (row[2], wurcs_to_tree(row[3]))  # (composition, tree)
    return db


def comp_of_tree(t):
    c = collections.Counter(ch for ch in t if ch in 'HNAGF')
    return (c['H'], c['N'], c['F'], c['A'], c['G'])


def bare(pep):
    """PIN peptides look like `K.SEQ[+glycan].K[tag]`; strip flanks, mods, tag."""
    m = re.match(r'^[A-Za-z_]\.(.*)\.[A-Za-z_](\[[^\]]*\])?$', pep, re.DOTALL)
    core = m.group(1) if m else pep
    core = re.sub(r'\[[^\]]*\]', '', core)
    core = re.sub(r'[+-]\d+\.\d+', '', core)
    return re.sub(r'[^A-Z]', '', core.upper()).replace('I', 'L')


def tree_of_andes(pep):
    m = re.search(r'\[([^\]]*)\]$', pep)
    if not m:
        return None
    tag = m.group(1)
    i = tag.find('(')
    return tag[i:] if i >= 0 else None


# ---------- engine N (GPSM) ----------
def load_engine_n(path, db):
    rows = {}                               # scan -> (peptide, tree, comp)
    with open(path) as f:
        r = csv.reader(f, delimiter='\t')
        next(r)
        for line in r:
            if len(line) < 19:
                continue
            try:
                if float(line[11]) > 0.01:  # PeptideFDR
                    continue
            except ValueError:
                continue
            scan = int(line[1])
            pep = bare(line[5])
            gid = line[16]
            if gid not in db:
                continue
            _, tree = db[gid]
            if tree is None:
                continue
            rows[scan] = (pep, tree, comp_of_tree(tree))
    return rows


# ---------- andes ----------
def load_andes(pin_path, psms_path, qmax):
    q = {}
    with open(psms_path) as f:
        r = csv.reader(f, delimiter='\t')
        qi = next(r).index('q-value')  # a `filename` column may precede it
        for line in r:
            if len(line) <= qi:
                continue
            try:
                if float(line[qi]) <= qmax:
                    q[line[0]] = float(line[qi])
            except ValueError:
                continue
    rows = {}
    with open(pin_path) as f:
        r = csv.reader(f, delimiter='\t')
        hdr = next(r)
        ci = {h: i for i, h in enumerate(hdr)}
        si, li, pi = ci['SpecId'], ci['Label'], ci['Peptide']
        for line in r:
            if len(line) <= max(si, li, pi):
                continue
            if line[li] != '1':
                continue
            spec = line[si]
            if spec not in q:
                continue
            m = re.search(r'scan=(\d+)', spec)
            if not m:
                continue
            scan = int(m.group(1))
            pep = line[pi]
            tree = tree_of_andes(pep)
            if tree is None:
                continue
            rows[scan] = (bare(pep), tree, comp_of_tree(tree))
    return rows


db = load_glycan_db(_args.glycan_db)
ref = load_engine_n(_args.engine_n, db)
andes = load_andes(_args.andes_pin, _args.andes_psms, _args.q)

co = set(ref) & set(andes)
same_pep = [s for s in co if ref[s][0] == andes[s][0]]
same_comp = [s for s in same_pep if ref[s][2] == andes[s][2]]
same_tree = [s for s in same_pep if ref[s][1] == andes[s][1]]
same_comp_diff_tree = [s for s in same_comp if ref[s][1] != andes[s][1]]

print("=" * 72)
print("MouseLiver-Z-T-1  |  glycan STRUCTURE (topology) agreement, andes vs engine N")
print("=" * 72)
print(f"  engine N accepted PSMs (PeptideFDR<=0.01): {len(ref):5d} distinct scans")
print(f"  andes accepted PSMs (q<={_args.q}):            {len(andes):5d} distinct scans")
print(f"  co-identified scans (both accepted):          {len(co):5d}")
print()
print(f"  same backbone peptide (I/L-folded):           {len(same_pep):5d}  ({100*len(same_pep)/max(len(co),1):4.1f}% of co-identified)")
print(f"    -- same glycan COMPOSITION:                 {len(same_comp):5d}  ({100*len(same_comp)/max(len(same_pep),1):4.1f}% of same-peptide)")
print(f"    -- same glycan STRUCTURE (tree):            {len(same_tree):5d}  ({100*len(same_tree)/max(len(same_pep),1):4.1f}% of same-peptide)")
print()
print(f"  structure agreement (headline)                {100*len(same_tree)/max(len(same_pep),1):4.1f}%")
print(f"  same composition but DIFFERENT structure:     {len(same_comp_diff_tree):5d}  ({100*len(same_comp_diff_tree)/max(len(same_comp),1):4.1f}% of same-composition)")

# set-level: distinct (peptide, tree) and (peptide, comp)
ref_pep_tree = set((ref[s][0], ref[s][1]) for s in ref)
and_pep_tree = set((andes[s][0], andes[s][1]) for s in andes)
ref_pep_comp = set((ref[s][0], ref[s][2]) for s in ref)
and_pep_comp = set((andes[s][0], andes[s][2]) for s in andes)
print()
print("  distinct (peptide, structure):  engine N %5d   andes %5d   shared %5d"
      % (len(ref_pep_tree), len(and_pep_tree), len(ref_pep_tree & and_pep_tree)))
print("  distinct (peptide, composition): engine N %5d   andes %5d   shared %5d"
      % (len(ref_pep_comp), len(and_pep_comp), len(ref_pep_comp & and_pep_comp)))

# a few example disagreements (same peptide, same comp, different tree)
print()
print("  example same-composition / different-structure disagreements:")
for s in sorted(same_comp_diff_tree)[:8]:
    print(f"    scan {s}:  {ref[s][0]}  comp {ref[s][2]}")
    print(f"        engine N {ref[s][1]}")
    print(f"        andes    {andes[s][1]}")
