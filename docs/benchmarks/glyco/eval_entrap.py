#!/usr/bin/env python3
"""Measure the entrapment false-discovery proportion of a glyco search.

`build_entrap.py` appends an unrelated proteome to the search database as TARGETS.
Those proteins are not in the sample, and prokaryotes/yeast lack the mammalian
N-glycosylation machinery, so any accepted glycopeptide on an `ENTRAP_` protein is
false by construction. Counting them gives an absolute error estimate that does not
depend on a reference identification set, and that no amount of search-space expansion
can flatter.

This is the check that must accompany `eval_yield.py`. Yield alone rises whenever the
candidate space grows, real or not: the full 4034-composition glycan list looked like a
gain of 59 compositions on yield and turned out to inflate this number 5.4x.

RAW `entrap / (entrap + sample)` UNDERSTATES the true FDP, and badly. False IDs land in
the sample proteome and the entrapment proteome roughly in proportion to how much of the
searchable space each contributes, so only the entrapment SHARE of the errors is visible:

    FDP_true ~= (n_entrap / total) * (1 + S/E)

where S and E size the two halves. For a glyco search the right denominator is NOT the
protein count but the count of N-X-S/T SEQUONS, because a glycopeptide can only be
reported at a sequon. Measured on the plasma benchmark database: human 20,411 proteins /
58,955 sequons against E. coli 4,531 proteins / 6,032 sequons -- so the protein-count
ratio is 4.50 but the sequon ratio is 9.77, and the correction factors are 5.50x and
10.77x respectively. Reporting the raw fraction called a run "CONSERVATIVE (0.39%)" whose
sequon-corrected FDP was ~4.3%.

Pass the search FASTA to get the corrected number. Without it the raw fraction is printed
with a loud warning and NO verdict, because a verdict on the raw fraction is misleading.

The search-space denominator is core-dependent: N-glycan reports only at N-X-S/T sequons,
O-glycan at any Ser/Thr residue. Pass `--core o` for an O-glycan search (`--glyco-core o`).

Isoleucine and leucine are isobaric, so an entrapment-only PSM whose bare peptide, with I
folded to L, occurs anywhere in a folded sample protein is indistinguishable from a sample
peptide and is counted as target, not entrapment (the rule lives in ../entrap_common.py and
is shared with reproduce/run.sh). This rescue needs the search FASTA (argv[4]) -- the same
file that sizes the S/E correction.

Usage:
    eval_entrap.py [--core n|o] <pooled.pin> <percolator.psms> [q-threshold] [search.fasta]
"""
import os
import re
import sys
from collections import defaultdict

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from entrap_common import SampleSequences, classify  # noqa: E402


def _chi2_ppf(p, k):
    """Chi-square quantile via bisection on a series gamma CDF. k is a positive int.

    Avoids a scipy dependency: the benchmark VM has only the stdlib.
    """
    import math

    def cdf(x):
        if x <= 0:
            return 0.0
        a = k / 2.0
        # lower regularised incomplete gamma, series expansion
        s, term = 1.0 / a, 1.0 / a
        for n in range(1, 500):
            term *= (x / 2.0) / (a + n)
            s += term
            if term < s * 1e-14:
                break
        return s * math.exp(-x / 2.0 + a * math.log(x / 2.0) - math.lgamma(a))

    lo, hi = 0.0, 1.0
    while cdf(hi) < p:
        hi *= 2.0
        if hi > 1e6:
            return hi
    for _ in range(200):
        mid = 0.5 * (lo + hi)
        if cdf(mid) < p:
            lo = mid
        else:
            hi = mid
    return 0.5 * (lo + hi)


if len(sys.argv) < 3:
    sys.exit(__doc__)
core = "n"
args = sys.argv[1:]
if args and args[0] == "--core":
    core = args[1]
    args = args[2:]
if len(args) < 2:
    sys.exit(__doc__)
pin, psms = args[0], args[1]
qmax = float(args[2]) if len(args) > 2 else 0.01
fasta = args[3] if len(args) > 3 else None


def space_ratio(path, core="n"):
    """(sample, entrap) site counts, and the protein counts as a fallback.

    The site count is the correct denominator for a glyco search: a glycopeptide can only
    be reported at an N-X-S/T sequon (N-glycan) or a Ser/Thr residue (O-glycan), so a
    proteome with few sites contributes little to the space in which a FALSE glyco ID can
    be made, regardless of its protein count.

    The counts are unaffected by I/L folding: neither site pattern tests I or L.
    """
    # N-glycan: N-X-S/T sequons. O-glycan: every Ser/Thr is a potential glycosite.
    pat = re.compile(r"[ST]" if core == "o" else r"N[^P][ST]")
    s_seq = e_seq = s_prot = e_prot = 0
    hdr, buf = None, []

    def flush(hdr, seq):
        nonlocal s_seq, e_seq, s_prot, e_prot
        if hdr is None:
            return
        n = len(pat.findall(seq))
        # SUBSTRING, not prefix: the accession that carries the tag may be inside a
        # UniProt-style header (`>sp|ENTRAP_Q9XXXX|...`), which is what the cluster's
        # entrapment databases emit. The hit-counting pass below already tests
        # `"ENTRAP_" in accession`, so a prefix test here made the two halves of this
        # script disagree: entrapment hits were counted while the space that sizes
        # their correction factor came out empty, and every run silently reported
        # `entrapment FDP: UNKNOWN` / `NO VERDICT` on a database that plainly had
        # entrapment proteins in it.
        # ORDER MATTERS: exclude decoys FIRST. A decoy built from an entrapment
        # protein keeps the tag inside the accession (">XXX_sp|ENTRAP_Q9XXXX|...");
        # testing "ENTRAP_" before the decoy prefixes counted those decoys as
        # entrapment SPACE, which inflates the denominator that sizes the
        # correction factor and therefore reports every FDP too LOW -- a run would
        # silently pass a verdict it should fail.
        if hdr.startswith(">XXX_") or hdr.startswith(">DECOY"):
            return
        if "ENTRAP_" in hdr:
            e_prot += 1
            e_seq += n
        else:
            s_prot += 1
            s_seq += n

    with open(path) as fh:
        for line in fh:
            if line.startswith(">"):
                flush(hdr, "".join(buf))
                hdr, buf = line.strip(), []
            else:
                buf.append(line.strip())
    flush(hdr, "".join(buf))
    return s_seq, e_seq, s_prot, e_prot


# SpecId -> accessions. Proteins occupy every column after "Peptide". The Peptide cell is
# kept per SpecId for the I/L-equivalence rescue in the counting pass.
prot_of = {}
pep_of = {}
header = None
with open(pin) as fh:
    for line in fh:
        f = line.rstrip("\n").split("\t")
        if header is None:
            header = f
            i_pep = header.index("Peptide")
            continue
        prot_of[f[0]] = f[i_pep + 1:]
        pep_of[f[0]] = f[i_pep]

sample = SampleSequences(fasta) if fasta is not None else None

n_entrap = n_sample = n_il = 0
per_charge = defaultdict(lambda: [0, 0])
QI = open(psms).readline().rstrip("\n").split("\t").index("q-value")  # Percolator column order varies (filename column)
with open(psms) as fh:
    for line in fh:
        f = line.rstrip("\n").split("\t")
        if f[0] == "PSMId" or len(f) < 3:
            continue
        try:
            if float(f[QI]) > qmax:
                continue
        except ValueError:
            continue
        accs = prot_of.get(f[0])
        if accs is None:
            continue
        # Entrapment only if EVERY accession is an entrapment protein and the peptide is not
        # an I/L variant of a sequence in a sample protein (see entrap_common.classify).
        c = classify(accs, pep_of[f[0]], sample)
        if c == "entrap":
            n_entrap += 1
        else:
            n_sample += 1
            if c == "rescued":
                n_il += 1

total = n_entrap + n_sample
raw = (n_entrap / total) if total else 0.0
print(f"q threshold      : {qmax:.4f}")
print(f"sample glycoPSMs : {n_sample}")
print(f"entrapment hits  : {n_entrap}")
if n_il:
    print(f"I/L-equivalent   : {n_il} entrapment PSMs rescued as target (I/L-folded peptide occurs in a sample protein)")
print(f"raw entrap frac  : {raw * 100:.2f}%   (entrap / accepted -- NOT the FDP)")

if total == 0:
    print("verdict          : NO ACCEPTED IDs — nothing to judge")
elif fasta is None:
    print("entrapment FDP   : UNKNOWN — pass the search FASTA as argv[4]")
    print("verdict          : NO VERDICT — the raw fraction understates FDP by (1 + S/E);")
    print("                   on the plasma benchmark that factor is ~10x by sequon count.")
else:
    s_seq, e_seq, s_prot, e_prot = space_ratio(fasta, core)
    if e_seq == 0 or e_prot == 0:
        print(f"entrapment FDP   : UNKNOWN — no ENTRAP_ proteins found in {fasta}")
        print("verdict          : NO VERDICT")
    else:
        f_seq = 1.0 + s_seq / e_seq
        f_prot = 1.0 + s_prot / e_prot
        fdp = raw * f_seq
        print(f"space (sites)    : sample {s_seq} / entrap {e_seq}  -> correction {f_seq:.2f}x")
        print(f"space (proteins) : sample {s_prot} / entrap {e_prot}  -> correction {f_prot:.2f}x")
        print(f"entrapment FDP   : {fdp * 100:.2f}%   (site-corrected)")
        # Poisson 95% CI on the observed count: with a handful of hits the point
        # estimate is far less certain than its two decimal places suggest.
        if n_entrap == 0:
            lo, hi = 0.0, 3.689
        else:
            import math
            lo = 0.5 * _chi2_ppf(0.025, 2 * n_entrap)
            hi = 0.5 * _chi2_ppf(0.975, 2 * n_entrap + 2)
        print(f"  95% CI on FDP  : {100 * lo / total * f_seq:.2f}% - {100 * hi / total * f_seq:.2f}%"
              f"   ({n_entrap} hit{'s' if n_entrap != 1 else ''})")
        if fdp <= qmax:
            print(f"verdict          : CONSERVATIVE (FDP <= nominal {qmax * 100:.1f}%)")
        else:
            print("verdict          : OPTIMISTIC — reported FDR understates the real error")
