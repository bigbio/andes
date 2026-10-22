#!/usr/bin/env python3
"""Build a peptide-level shuffled entrapment database.

    build_shuffled_entrap.py <targets.fasta> <out.fasta> [--seed N] [--proline-rule]

The output holds every target protein verbatim, then one entrapment protein per target,
tagged `ENTRAP_` in the accession and the entry name (`>sp|ENTRAP_Q99JY4|ENTRAP_TRABD_MOUSE`).
The script ends by printing the database's T/E ratio, measured in peptide space by
`entrap_common.py te`. Use that ratio in FDP = (entrapment hits / accepted) x (1 + T/E);
it is close to 1 but never exactly 1.

CONSTRUCTION (the peptide-level entrapment of Wen et al., Nat. Methods 2025)
  1. Each target protein is cut into fully cleaved tryptic peptides.
  2. Each peptide is shuffled with its C-terminal residue fixed, so every entrapment peptide
     has the same length, composition and mass as its target, and the cleavage sites stay
     where they were. Missed-cleavage peptides therefore pair up as well.
  3. A peptide that occurs more than once is shuffled once and replaced the same way
     everywhere.
  4. A shuffle that reproduces a target peptide (I and L folded, as they are isobaric) is
     drawn again, up to 100 times. Peptides that cannot avoid it (too short, or all residues
     alike) stay as they are; `entrap_common.py` counts hits on them as target.
  5. The shuffled peptides are joined back in their original order.
  The seed makes the output deterministic: the same input and seed give the same file.

ENZYME
  The default cuts after every K and R, including before P. That is andes's `trypsin` (it has
  no proline rule) and what the T/E digest in `entrap_common.py` assumes. `--proline-rule`
  does not cut before P, for searches with a proline rule. It then also keeps every K, R
  and P of a peptide in place, so a shuffle neither creates nor removes a cleavage site.

WHY NOT SHUFFLE WHOLE PROTEINS
  Shuffling a whole protein moves its K and R, so its peptides have no composition- or
  mass-matched target: on a human/yeast/E. coli database 92% of such entrapment peptides
  have no target of the same composition. They score unlike real false matches and the
  FDP comes out low (by 16-36% against a foreign proteome or a peptide-level shuffle on the
  same results). The size ratio is not 1:1 either: whole-protein shuffles of that database
  give T/E about 0.90.

FOREIGN PROTEOMES
  `glyco/build_entrap.py` appends an unrelated proteome instead. Both designs are valid; a
  foreign proteome has its own size ratio, and for glyco its own sequon density, so measure
  T/E (and for glyco the sequon ratio, `glyco/eval_entrap.py`) on the database you built.
"""
import argparse
import os
import random
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import entrap_common as ec  # noqa: E402

MAX_TRIES = 100


def cleave(seq, proline_rule=False):
    """Fully cleaved tryptic peptides of `seq`, in order; they concatenate back to `seq`."""
    out, start = [], 0
    for i, c in enumerate(seq):
        if c in "KR" and i + 1 < len(seq) and not (proline_rule and seq[i + 1] == "P"):
            out.append(seq[start:i + 1])
            start = i + 1
    if start < len(seq):
        out.append(seq[start:])
    return out


def shuffle_peptide(pep, rng, proline_rule=False):
    """`pep` with its residues shuffled, except the C-terminal one (and, with the proline
    rule, every K, R and P)."""
    fixed = "KRP" if proline_rule else ""
    free = [i for i in range(len(pep) - 1) if pep[i] not in fixed]
    chars = [pep[i] for i in free]
    rng.shuffle(chars)
    out = list(pep)
    for i, c in zip(free, chars):
        out[i] = c
    return "".join(out)


def shuffleable(pep, proline_rule=False):
    """Can a shuffle change `pep` at all (two different residues among the free ones)?"""
    fixed = "KRP" if proline_rule else ""
    return len({c for c in pep[:-1] if c not in fixed}) > 1


def build(entries, seed=1, proline_rule=False):
    """Entrapment sequences for `entries` [(header, sequence)], in order, and statistics."""
    rng = random.Random(seed)
    peptides = [cleave(seq, proline_rule) for _, seq in entries]
    target = {ec.fold_il(p) for peps in peptides for p in peps}
    mapping = {}
    stats = {"peptides": 0, "distinct": 0, "unchanged": 0, "collisions": 0, "collisions_6": 0}
    out = []
    for peps in peptides:
        parts = []
        for p in peps:
            stats["peptides"] += 1
            s = mapping.get(p)
            if s is None:
                stats["distinct"] += 1
                s = p
                if not shuffleable(p, proline_rule):
                    stats["unchanged"] += 1
                else:
                    for _ in range(MAX_TRIES):
                        s = shuffle_peptide(p, rng, proline_rule)
                        if ec.fold_il(s) not in target:
                            break
                    else:
                        stats["collisions"] += 1
                        stats["collisions_6"] += len(p) >= 6
                mapping[p] = s
            parts.append(s)
        out.append("".join(parts))
    return out, stats


def tag(hdr):
    """`>sp|Q99JY4|TRABD_MOUSE ...` -> `>sp|ENTRAP_Q99JY4|ENTRAP_TRABD_MOUSE ...`"""
    parts = hdr.split("|")
    if len(parts) >= 3:
        parts[1] = "ENTRAP_" + parts[1]
        parts[2] = "ENTRAP_" + parts[2]
        return "|".join(parts)
    return ">ENTRAP_" + hdr.lstrip(">")


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__.split("\n\n")[0],
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("targets")
    ap.add_argument("out")
    ap.add_argument("--seed", type=int, default=1)
    ap.add_argument("--proline-rule", action="store_true",
                    help="do not cut before P, and keep K, R and P in place")
    ap.add_argument("--no-te", action="store_true", help="skip the T/E measurement")
    a = ap.parse_args(argv)

    entries = list(ec.read_fasta(a.targets))
    if any(ec.protein_class(h) != "sample" for h, _ in entries):
        sys.exit(f"{a.targets}: holds decoy or ENTRAP_ entries; give the target proteins only")
    seqs, st = build(entries, a.seed, a.proline_rule)
    with open(a.out, "w") as fh:
        for hdr, seq in entries:
            fh.write(f"{hdr}\n{seq}\n")
        for (hdr, _), seq in zip(entries, seqs):
            fh.write(f"{tag(hdr)}\n{seq}\n")
    enzyme = "no cleavage before P" if a.proline_rule else "cleavage before P (andes trypsin)"
    print(f"wrote {a.out}: {len(entries)} targets + {len(seqs)} entrapment proteins "
          f"(peptide-level shuffle, C-terminus fixed, {enzyme}, seed={a.seed})")
    print(f"    peptides {st['peptides']}, distinct {st['distinct']}, unchanged "
          f"{st['unchanged']} (nothing to shuffle), still equal to a target peptide after "
          f"{MAX_TRIES} draws {st['collisions']} ({st['collisions_6']} of 6 or more residues)")
    if not a.no_te:
        t, e, moved = ec.te_ratio(a.out)
        r = t / max(e, 1)
        print(f"    target peptides {t}  entrapment-only peptides {e}  ({moved} found in a "
              f"sample protein, counted as target)")
        print(f"    T/E = {r:.3f}  ->  FDP = (hits/total) x {1 + r:.3f}")
    print("run andes on the output as on any database; with --glyco use "
          "--decoy-strategy sequon-reverse")
    return 0


if __name__ == "__main__":
    sys.exit(main())
