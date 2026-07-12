#!/usr/bin/env python3
"""Entrapment-FDR validation for the andes --glyco path.

Implements the standard published entrapment-FDP estimators (lower-bound,
combined, and paired). Entrapment estimates the false-discovery PROPORTION
directly from known-false discoveries, independent of the decoy count, so it
answers the objection that "two decoys at 1% is too few for a stable
calibration": averaging over random entrapment draws yields a stable
empirical-FDR curve.

Two subcommands:

  build-fasta   Generate a 1:1 PAIRED shuffled-entrapment FASTA (r=1) from a
                target-only FASTA. Each target protein is paired with one
                entrapment protein whose sequence is a per-protein seeded shuffle
                with the C-terminal residue fixed (tryptic specificity preserved)
                and N-X-S/T sequon density preserved (so entrapment backbones are
                searchable glyco candidates -- otherwise they can never be hit and
                the entrapment axis is dead). Entrapment accessions are prefixed
                'ENT_'; a pairing TSV (target_acc <TAB> entrapment_acc) is written.
                Concatenate this with the target FASTA (and your decoys) to search.

  fdp           Given an andes glyco identification table (Percolator .psms, or the
                glyco PIN) whose protein column marks entrapment hits with 'ENT_',
                compute the lower-bound, combined, and PAIRED FDP at each q cutoff,
                per the formulas below, and compare to the nominal q at 1%.

Estimators (per score cutoff s; n_tau = target discoveries, n_eps = entrapment
discoveries, r = entrapment:target ratio):
    lower_bound = n_eps / (n_eps + n_tau)                       # proves FAILURE only
    combined    = n_eps * (1 + 1/r) / (n_eps + n_tau)           # valid upper bound, conservative
    paired      = (n_eps + 2*N_pts + N_pst) / (n_eps + n_tau)   # valid upper bound, TIGHTER (r=1)
        N_pst = # entrapment discoveries whose PAIRED TARGET was NOT discovered (< s)
        N_pts = # entrapment discoveries whose PAIRED TARGET was discovered but scored LOWER
Ordering always: lower_bound <= paired <= combined. Verdict at the q=0.01 line:
    combined <= 0.01         -> evidence andes controls FDR (conservative)
    lower_bound > 0.01       -> andes FAILS to control FDR
NOTE: the exact paired estimator needs an explicit peptide<->entrapment map
(--pairs, not yet wired); until then the tool reports the COMBINED upper bound as
a conservative proxy in the paired column, so a <=1% reading still evidences
control -- it is just looser than the true paired estimator would be.
"""
import argparse
import re
import sys

SEQUON = re.compile(r"N[^P][ST]")


# ---------------------------------------------------------------- build-fasta ---
def _read_fasta(path):
    acc, seq, cur = None, [], []
    out = []
    for line in open(path):
        if line.startswith(">"):
            if acc is not None:
                out.append((acc, "".join(cur)))
            acc = line[1:].strip().split()[0]
            cur = []
        else:
            cur.append(line.strip())
    if acc is not None:
        out.append((acc, "".join(cur)))
    return out


def _xorshift(state):
    state ^= (state << 13) & 0xFFFFFFFFFFFFFFFF
    state ^= state >> 7
    state ^= (state << 17) & 0xFFFFFFFFFFFFFFFF
    return state & 0xFFFFFFFFFFFFFFFF


def _seeded_shuffle_fix_cterm(seq, seed):
    """Fisher-Yates on seq[:-1] with a per-protein xorshift RNG; C-term fixed."""
    if len(seq) < 3:
        return seq
    body = list(seq[:-1])
    state = seed or 0xD1B54A32D192ED03
    for j in range(len(body) - 1, 0, -1):
        state = _xorshift(state)
        k = state % (j + 1)
        body[j], body[k] = body[k], body[j]
    return "".join(body) + seq[-1]


def _count_sequon_starts(s):
    return sum(
        1 for i in range(len(s) - 2)
        if s[i] == "N" and s[i + 1] != "P" and s[i + 2] in "ST"
    )


def _restore_sequon_density(shuf, target_count):
    """Greedily raise the shuffle's sequon-start count toward target_count via
    composition-preserving swaps (same rule as andes' SequonReverse decoy). Each
    swap can also break/create a sequon at the donor window, so re-count the actual
    total and keep the swap only if it strictly increases the count without
    overshooting (a naive +1 is unreliable and would leave entrapment proteins with
    FEWER sequons than targets, biasing the FDP optimistically -- adversarial review)."""
    r = list(shuf)
    n = len(r)
    cur = _count_sequon_starts(shuf)
    i = 0
    while cur < target_count and i + 2 < n:
        if r[i] == "N" and r[i + 1] != "P" and r[i + 2] not in "ST":
            j = i + 3
            while j < n:
                if r[j] in "ST" and not (j >= 2 and r[j - 2] == "N" and r[j - 1] != "P"):
                    r[i + 2], r[j] = r[j], r[i + 2]
                    new_count = _count_sequon_starts("".join(r))
                    if new_count > cur and new_count <= target_count:
                        cur = new_count
                    else:
                        r[i + 2], r[j] = r[j], r[i + 2]  # revert
                    break
                j += 1
        i += 1
    return "".join(r)


def build_fasta(args):
    prots = _read_fasta(args.target)
    seen = set(s for _, s in prots)
    out = open(args.out_fasta, "w")
    pair = open(args.out_pairs, "w")
    pair.write("target_acc\tentrapment_acc\n")
    n_made = 0
    for idx, (acc, seq) in enumerate(prots):
        target_sequons = _count_sequon_starts(seq)
        ent = None
        base = (args.seed ^ (idx * 0x9E3779B97F4A7C15)) & 0xFFFFFFFFFFFFFFFF
        for attempt in range(args.retries):
            cand = _seeded_shuffle_fix_cterm(seq, (base + attempt) & 0xFFFFFFFFFFFFFFFF)
            cand = _restore_sequon_density(cand, target_sequons)
            if cand not in seen and cand != seq:
                ent = cand
                break
        if ent is None:
            continue  # could not build a distinct entrapment; drop (rare)
        seen.add(ent)
        ent_acc = f"ENT_{acc}"
        out.write(f">{ent_acc}\n{ent}\n")
        pair.write(f"{acc}\t{ent_acc}\n")
        n_made += 1
    out.close()
    pair.close()
    print(f"built {n_made}/{len(prots)} paired entrapment proteins -> {args.out_fasta}")
    print(f"pairing -> {args.out_pairs}. Concatenate {args.target} + {args.out_fasta} "
          f"(+ decoys) as the search DB.")


# ------------------------------------------------------------------------ fdp ---
def _bare_peptide(pep):
    s = re.sub(r"\[[^\]]*\]", "", pep)
    s = re.sub(r"[+-]\d+\.\d+", "", s)
    parts = s.split(".")
    core = parts[1] if len(parts) >= 3 else parts[0]
    return re.sub(r"[^A-Z]", "", core.upper())


def fdp(args):
    # Read discoveries: (score_for_sort, q, is_entrapment, peptide_key).
    # score: we sort by q ascending (best first). is_entrapment: protein all ENT_.
    rows = []
    with open(args.ids) as f:
        hdr = next(f).rstrip("\n").split("\t")
        # tolerate either Percolator .psms (PSMId/q-value/proteinIds) or a PIN.
        qi = hdr.index("q-value") if "q-value" in hdr else None
        pi = hdr.index("proteinIds") if "proteinIds" in hdr else None
        pepi = hdr.index("peptide") if "peptide" in hdr else (
            hdr.index("Peptide") if "Peptide" in hdr else None)
        if qi is None or pi is None:
            sys.exit("need a Percolator .psms with q-value + proteinIds columns")
        dec = args.decoy_prefix
        n_decoy = 0
        for ln in f:
            p = ln.rstrip("\n").split("\t")
            if len(p) <= pi:
                continue
            try:
                q = float(p[qi])
            except ValueError:
                continue
            accs = [a for a in p[pi:] if a]
            # Classify each discovery. DECOY rows (all accessions carry the decoy
            # prefix) are EXCLUDED entirely — counting them in n_tau would inflate the
            # denominator and understate the FDP (adversarial review). A Percolator
            # --results-psms table already excludes decoys, but be robust if a raw
            # combined table is passed.
            if bool(accs) and all(a.startswith(dec) for a in accs):
                n_decoy += 1
                continue
            is_ent = bool(accs) and all(a.startswith("ENT_") for a in accs)
            pepkey = _bare_peptide(p[pepi]) if pepi is not None else str(len(rows))
            rows.append((q, is_ent, pepkey))
    rows.sort(key=lambda r: r[0])  # best (lowest q) first
    if n_decoy:
        print(f"# excluded {n_decoy} decoy rows (prefix '{dec}') from the target count")

    # Optional peptide<->entrapment pairing (enables the EXACT paired estimator).
    # File: two columns, target_peptide <TAB> entrapment_peptide (bare sequences).
    ent_to_tgt = {}
    if args.pairs:
        for ln in open(args.pairs):
            if not ln.strip() or ln.startswith("#") or ln.lower().startswith("target"):
                continue
            a, b = ln.rstrip("\n").split("\t")[:2]
            ent_to_tgt[re.sub(r"[^A-Z]", "", b.upper())] = re.sub(r"[^A-Z]", "", a.upper())
    # Best (lowest) q at which each TARGET peptide was discovered — for the paired terms.
    tgt_best_q = {}
    for q, e, pep in rows:
        if not e:
            if pep not in tgt_best_q or q < tgt_best_q[pep]:
                tgt_best_q[pep] = q

    r_ratio = args.ratio
    have_paired = bool(ent_to_tgt)
    col = "paired" if have_paired else "comb(=paired proxy)"
    print(f"{'q<=':>7} {'n_tau':>7} {'n_eps':>6} {'lower':>8} {'combined':>9} {col:>20}")
    for thr in (0.005, 0.01, 0.02, 0.05):
        ntau = sum(1 for q, e, _ in rows if q <= thr and not e)
        neps = sum(1 for q, e, _ in rows if q <= thr and e)
        denom = neps + ntau
        lower = neps / denom if denom else 0.0
        combined = neps * (1 + 1.0 / r_ratio) / denom if denom else 0.0
        if have_paired and denom:
            # EXACT paired (r=1): for each entrapment discovery at q<=thr, inspect its
            # PAIRED target. N_pst = paired target NOT discovered at thr; N_pts = paired
            # target discovered but scored WORSE (higher q) than the entrapment.
            n_pst = n_pts = 0
            for q, e, pep in rows:
                if not e or q > thr:
                    continue
                t = ent_to_tgt.get(pep)
                tq = tgt_best_q.get(t) if t is not None else None
                if tq is None or tq > thr:
                    n_pst += 1
                elif tq > q:
                    n_pts += 1
            paired = (neps + 2 * n_pts + n_pst) / denom
        else:
            paired = combined  # conservative proxy when no pairing supplied
        flag = "  <-- 1% line" if abs(thr - 0.01) < 1e-9 else ""
        print(f"{thr:>7.3f} {ntau:>7d} {neps:>6d} {100*lower:>7.2f}% "
              f"{100*combined:>8.2f}% {100*paired:>17.2f}%{flag}")
    v = "paired" if have_paired else "combined"
    print(f"\nVerdict: {v}<=1% at the 1% line -> evidence of FDR control; "
          "lower>1% -> FDR NOT controlled." + ("" if have_paired else
          " (Combined is the conservative proxy for the tighter paired estimator, "
          "which needs --pairs from a peptide-level entrapment build-pep.)") +
          " Average over >=20 random draws for a stable empirical-FDR curve.")


def _tryptic(seq):
    peps, start = [], 0
    for i, a in enumerate(seq):
        if a in "KR" and (i + 1 >= len(seq) or seq[i + 1] != "P"):
            peps.append(seq[start:i + 1])
            start = i + 1
    if start < len(seq):
        peps.append(seq[start:])
    return peps


def build_pep(args):
    """PEPTIDE-level 1:1 entrapment: for each unique tryptic SEQUON-bearing target
    peptide, emit one shuffled entrapment peptide (C-term fixed, sequon density
    preserved), each as a single-peptide protein `>ENT_pep_<i>`. Writes the search
    FASTA (concatenate with the target) + a target_peptide<TAB>entrapment_peptide
    pairing for the EXACT paired FDP estimator (`fdp --pairs`)."""
    prots = _read_fasta(args.target)
    tgt_peps = {}
    for _, seq in prots:
        for p in _tryptic(seq):
            if 6 <= len(p) <= 50 and SEQUON.search(p) and p not in tgt_peps:
                tgt_peps[p] = None
    out = open(args.out_fasta, "w")
    pair = open(args.out_pairs, "w")
    pair.write("target_peptide\tentrapment_peptide\n")
    seen = set(tgt_peps)
    n = 0
    for idx, tp in enumerate(tgt_peps):
        tc = _count_sequon_starts(tp)
        base = (args.seed ^ (idx * 0x9E3779B97F4A7C15)) & 0xFFFFFFFFFFFFFFFF
        ep = None
        for attempt in range(args.retries):
            cand = _restore_sequon_density(
                _seeded_shuffle_fix_cterm(tp, (base + attempt) & 0xFFFFFFFFFFFFFFFF), tc)
            if cand not in seen:
                ep = cand
                break
        if ep is None:
            continue
        seen.add(ep)
        out.write(f">ENT_pep_{idx}\n{ep}\n")
        pair.write(f"{tp}\t{ep}\n")
        n += 1
    out.close()
    pair.close()
    print(f"built {n}/{len(tgt_peps)} peptide-level entrapment peptides -> {args.out_fasta}")
    print(f"pairing -> {args.out_pairs}. Search {args.target} + {args.out_fasta} (+ decoys), "
          f"then: fdp --ids <psms> --pairs {args.out_pairs}")


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    sub = ap.add_subparsers(dest="cmd", required=True)
    b = sub.add_parser("build-fasta", help="build a paired shuffled entrapment FASTA")
    b.add_argument("--target", required=True, help="target-only FASTA")
    b.add_argument("--out-fasta", required=True)
    b.add_argument("--out-pairs", required=True)
    b.add_argument("--seed", type=int, default=42)
    b.add_argument("--retries", type=int, default=20)
    b.set_defaults(func=build_fasta)
    bp = sub.add_parser("build-pep", help="build a PEPTIDE-level 1:1 entrapment FASTA + pairing")
    bp.add_argument("--target", required=True, help="target-only FASTA")
    bp.add_argument("--out-fasta", required=True)
    bp.add_argument("--out-pairs", required=True)
    bp.add_argument("--seed", type=int, default=42)
    bp.add_argument("--retries", type=int, default=20)
    bp.set_defaults(func=build_pep)
    f = sub.add_parser("fdp", help="compute entrapment FDP from an ID table")
    f.add_argument("--ids", required=True, help="Percolator .psms (q-value + proteinIds)")
    f.add_argument("--ratio", type=float, default=1.0, help="entrapment:target ratio r")
    f.add_argument("--decoy-prefix", default="XXX_",
                   help="accession prefix marking decoy rows to exclude from n_tau")
    f.add_argument("--pairs", default=None,
                   help="target_peptide<TAB>entrapment_peptide map (from build-pep) "
                        "→ enables the EXACT paired estimator instead of the combined proxy")
    f.set_defaults(func=fdp)
    args = ap.parse_args()
    args.func(args)


if __name__ == "__main__":
    main()
