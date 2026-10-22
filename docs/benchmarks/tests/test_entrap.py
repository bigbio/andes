#!/usr/bin/env python3
"""Tests for the shared entrapment rules (../entrap_common.py) and their three users.

Fixture (fixtures/il_entrap.*): one sample protein containing ...TLEADKPWQR... at a
non-tryptic position, one ENTRAP_ protein containing ...TIEADKPWQR..., a glyco-style PIN and
the matching Percolator target PSM file. Accepted rows at q <= 0.01:

    s1  sample peptide                                  -> sample
    s2  ENTRAP_ only, TIEADKPWQR (I/L variant)          -> rescued
    s3  ENTRAP_ only, FFFFNYTWWR                        -> entrapment hit
    s4  ENTRAP_ and sample protein (shared)             -> sample (ALL rule)
    s5  as s2, with a bracketed mod inside the sequence -> rescued
    s6  entrapment at q 0.05                            -> not accepted

Standard library only, no network:  python3 -I docs/benchmarks/tests/test_entrap.py
"""
import ast
import csv
import os
import re
import stat
import subprocess
import sys
import tempfile
import unittest

HERE = os.path.dirname(os.path.abspath(__file__))
BENCH = os.path.dirname(HERE)
FIX = os.path.join(HERE, "fixtures")
FASTA = os.path.join(FIX, "il_entrap.fasta")
PIN = os.path.join(FIX, "il_entrap.pin")
PSMS = os.path.join(FIX, "il_entrap.t.psms")

sys.path.insert(0, BENCH)
import entrap_common as ec  # noqa: E402


def run(args, **kw):
    return subprocess.run(args, capture_output=True, text=True, check=True, **kw).stdout


class BarePeptide(unittest.TestCase):
    def test_formats(self):
        cases = {
            "K.TAANGTR.K[HexNAc2Hex8Fuc0NeuAc0NeuGc0@N4]": "TAANGTR",
            "K.PEPS+79.96633TIDE.K": "PEPSTIDE",
            "R.A+229.16293VFVK+229.16293.D": "AVFVK",
            "-.Q-17.02655PEPTIDEK.-": "QPEPTIDEK",
            "K.TIEADK[UNIMOD:1]PWQR.F[HexNAc2Hex5Fuc0NeuAc0NeuGc0@N?]": "TIEADKPWQR",
            "K.M[UNIMOD:35]PEPK.R": "MPEPK",
            "K.NM+15.99491ST.R[HexNAc1Hex0Fuc0NeuAc0NeuGc0@S/T?]": "NMST",
            "PEPTIDE": "PEPTIDE",
        }
        for cell, want in cases.items():
            self.assertEqual(ec.bare_peptide(cell), want, cell)


class SharedRules(unittest.TestCase):
    def setUp(self):
        self.sample = ec.SampleSequences(FASTA)

    def test_substring_rescue_is_enzyme_independent(self):
        # TLEADKPWQR starts at a non-tryptic position of the sample protein.
        self.assertTrue(self.sample.contains("TIEADKPWQR"))
        self.assertTrue(self.sample.contains("TLEADKPWQR"))
        self.assertFalse(self.sample.contains("FFFFNYTWWR"))
        # Entrapment sequences are not part of the sample text.
        self.assertFalse(self.sample.contains("MPPKTLEADK"))

    def test_no_match_across_proteins(self):
        two = os.path.join(tempfile.mkdtemp(), "two.fasta")
        with open(two, "w") as fh:
            fh.write(">sp|A|A\nPEPTAAK\n>sp|B|B\nRGGGEEK\n")
        self.assertFalse(ec.SampleSequences(two).contains("AAKRGG"))

    def test_all_rule(self):
        e = "ENTRAP_sp|P99999|ENT1_ECOLI"
        s = "sp|P00001|SAMP1_YEAST"
        self.assertEqual(ec.classify([e, s], "R.FFFFNYTWWR.-", self.sample), "sample")
        self.assertEqual(ec.classify([e, ""], "R.FFFFNYTWWR.-", self.sample), "entrap")
        self.assertEqual(ec.classify([e], "K.TIEADKPWQR.F", self.sample), "rescued")
        self.assertEqual(ec.classify([e], "K.TIEADKPWQR.F", None), "entrap")
        self.assertEqual(ec.classify([], "K.TIEADKPWQR.F", None), "sample")

    def test_count_percolator(self):
        self.assertEqual(ec.count_percolator_psms(PSMS, 0.01, self.sample), (5, 1, 2))
        self.assertEqual(ec.count_percolator_psms(PSMS, 0.01, None), (5, 3, 0))

    def test_te_moves_sample_substrings_to_target(self):
        # T: 6 sample peptides. E: 5 entrapment peptides, of which TLEADK and TLEADKPWQR
        # occur in the folded sample protein and move to T.
        self.assertEqual(ec.te_ratio(FASTA), (8, 3, 2))
        out = run([sys.executable, "-I", os.path.join(BENCH, "entrap_common.py"), "te", FASTA])
        self.assertIn("T/E = 2.667", out)


class ShuffledEntrap(unittest.TestCase):
    """build_shuffled_entrap.py: the peptide-level shuffle."""

    script = os.path.join(BENCH, "build_shuffled_entrap.py")
    targets = [
        (">sp|P00001|A_HUMAN", "MSDEQWLAGHKTPEFNMLCRGGYSTNDEPWVKAHMSEQK"),
        (">sp|P00002|B_HUMAN", "MPEDNVWQYTRSAHGFLEKTPEFNMLCRDWQHNEPSGAMYLK"),
        (">sp|P00003|C_YEAST", "MKRPAGWDYENSTLFKPMNWEHQRPDSAGLTVK"),
    ]

    @classmethod
    def setUpClass(cls):
        sys.path.insert(0, BENCH)
        import build_shuffled_entrap as bse
        cls.bse = bse

    def write(self, entries):
        path = os.path.join(tempfile.mkdtemp(), "t.fasta")
        with open(path, "w") as fh:
            fh.write("".join(f"{h}\n{s}\n" for h, s in entries))
        return path

    def build(self, entries, *flags):
        src = self.write(entries)
        out = src + ".out.fasta"
        log = run([sys.executable, "-I", self.script, src, out, *flags])
        with open(out) as fh:
            text = fh.read()
        db = list(ec.read_fasta(out))
        return db[:len(entries)], db[len(entries):], text, log

    def test_cterm_composition_and_cleavage_sites(self):
        for flags, rule in (((), False), (("--proline-rule",), True)):
            tgt, ent, _, _ = self.build(self.targets, "--seed", "7", *flags)
            self.assertEqual(tgt, self.targets)
            self.assertEqual(len(ent), len(tgt))
            for (th, ts), (eh, es) in zip(tgt, ent):
                self.assertIn("ENTRAP_", eh)
                self.assertNotIn("ENTRAP_", th)
                tp, ep = self.bse.cleave(ts, rule), self.bse.cleave(es, rule)
                # Same cleavage sites, so the same peptide boundaries.
                self.assertEqual([len(p) for p in tp], [len(p) for p in ep])
                for a, b in zip(tp, ep):
                    self.assertEqual(a[-1], b[-1], (a, b))
                    self.assertEqual(sorted(a), sorted(b), (a, b))
                    if rule:
                        for i, c in enumerate(a):
                            if c in "KRP":
                                self.assertEqual(b[i], c, (a, b))
                self.assertNotEqual(ts, es)

    def test_no_target_collisions_il_folded(self):
        # LAVK's free residues have six orders; five are targets, so the only allowed
        # shuffle of all five is VALK.
        five = ["LAVK", "LVAK", "ALVK", "AVLK", "VLAK"]
        entries = [(f">sp|Q{i}|P{i}_HUMAN", p) for i, p in enumerate(five)]
        _, ent, _, log = self.build(entries)
        self.assertEqual({s for _, s in ent}, {"VALK"})
        self.assertIn("after 100 draws 0", log)
        # VAIK folds to VALK: now no shuffle avoids a target, and that is reported.
        _, _, _, log = self.build(entries + [(">sp|Q9|P9_HUMAN", "VAIK")])
        self.assertIn("after 100 draws 6", log)

    def test_identical_peptides_shuffled_alike(self):
        tgt, ent, _, _ = self.build(self.targets)
        # TPEFNMLCR is peptide 1 of A and peptide 2 of B.
        a = self.bse.cleave(ent[0][1])[1]
        b = self.bse.cleave(ent[1][1])[2]
        self.assertEqual(self.bse.cleave(tgt[0][1])[1], "TPEFNMLCR")
        self.assertEqual(self.bse.cleave(tgt[1][1])[2], "TPEFNMLCR")
        self.assertEqual(a, b)
        self.assertNotEqual(a, "TPEFNMLCR")

    def test_deterministic(self):
        _, _, one, _ = self.build(self.targets, "--seed", "3")
        _, _, two, _ = self.build(self.targets, "--seed", "3")
        _, _, other, _ = self.build(self.targets, "--seed", "4")
        self.assertEqual(one, two)
        self.assertNotEqual(one, other)

    def test_reports_measured_te(self):
        _, _, _, log = self.build(self.targets)
        self.assertRegex(log, r"T/E = \d+\.\d{3}")
        self.assertNotIn("2.0 (it is exactly 1:1)", log)

    def test_rejects_entrapment_input(self):
        src = self.write(self.targets + [(">sp|ENTRAP_X|ENTRAP_X", "AAAK")])
        r = subprocess.run([sys.executable, "-I", self.script, src, src + ".o"],
                           capture_output=True, text=True)
        self.assertNotEqual(r.returncode, 0)


class EvalEntrap(unittest.TestCase):
    script = os.path.join(BENCH, "glyco", "eval_entrap.py")

    def field(self, out, name):
        m = re.search(rf"^{re.escape(name)}\s*:\s*(\d+)", out, re.M)
        return int(m.group(1)) if m else 0

    def test_rescue_with_fasta(self):
        out = run([sys.executable, "-I", self.script, PIN, PSMS, "0.01", FASTA])
        self.assertEqual(self.field(out, "sample glycoPSMs"), 4)
        self.assertEqual(self.field(out, "entrapment hits"), 1)
        self.assertEqual(self.field(out, "I/L-equivalent"), 2)

    def test_no_rescue_without_fasta(self):
        out = run([sys.executable, "-I", self.script, PIN, PSMS, "0.01"])
        self.assertEqual(self.field(out, "sample glycoPSMs"), 2)
        self.assertEqual(self.field(out, "entrapment hits"), 3)


class CompareEngines(unittest.TestCase):
    """compare_engines.py runs end to end at import, so its andes loader is lifted out of
    the source and executed on its own."""

    @classmethod
    def setUpClass(cls):
        with open(os.path.join(BENCH, "glyco", "compare_engines.py")) as fh:
            src = fh.read()
        wanted = {"acc", "parse_comp_andes", "load_andes"}
        tree = ast.parse(src)
        body = [n for n in tree.body if isinstance(n, ast.FunctionDef) and n.name in wanted]
        ns = {"csv": csv, "re": re, "bare_peptide": ec.bare_peptide, "classify": ec.classify}
        exec(compile(ast.Module(body=body, type_ignores=[]), "compare_engines", "exec"), ns)
        cls.load_andes = staticmethod(ns["load_andes"])

    def test_rescue_and_all_rule(self):
        rows, n_entrap = self.load_andes(PIN, PSMS, 0.01, ec.SampleSequences(FASTA))
        self.assertEqual(len(rows), 5)
        self.assertEqual(n_entrap, 1)
        flags = {p: e for p, _, _, e in rows}
        self.assertFalse(flags["TIEADKPWQR"])
        self.assertTrue(flags["FFFFNYTWWR"])
        # Shared peptide whose FIRST protein is ENTRAP_: not entrapment under the ALL rule.
        self.assertFalse(flags["PWQRSAMPLEK"])
        self.assertEqual(rows[0][1], (2, 5, 0, 0, 0))

    def test_without_fasta(self):
        _, n_entrap = self.load_andes(PIN, PSMS, 0.01, None)
        self.assertEqual(n_entrap, 3)


class RunSh(unittest.TestCase):
    """reproduce/run.sh end to end, with stand-ins for the search engine and Percolator."""

    def test_ups1_counts(self):
        tmp = tempfile.mkdtemp()
        data, binp = os.path.join(tmp, "data"), os.path.join(tmp, "bin")
        for d in ("ups1", "databases", "results"):
            os.makedirs(os.path.join(data, d))
        os.makedirs(binp)
        with open(os.path.join(data, "ups1", "x.mzML"), "w") as fh:
            fh.write("stub\n")
        with open(FASTA) as src, open(os.path.join(data, "databases", "yeast_entrap.fasta"), "w") as dst:
            dst.write(src.read())
        stubs = {
            "andes": "#!/bin/sh\necho 'andes stub'\n",
            # Percolator stand-in: write the fixture where run.sh expects the target PSMs.
            "docker": '#!/bin/sh\ncp "$FIXTURE_PSMS" "$DATA_DIR/results/ups1.t.psms"\n',
        }
        for name, body in stubs.items():
            p = os.path.join(binp, name)
            with open(p, "w") as fh:
                fh.write(body)
            os.chmod(p, os.stat(p).st_mode | stat.S_IEXEC)
        env = dict(os.environ, PATH=binp + os.pathsep + os.environ["PATH"],
                   FIXTURE_PSMS=PSMS, DATA_DIR=data)
        out = run(["bash", os.path.join(BENCH, "reproduce", "run.sh"), data, "ups1"], env=env)
        row = [ln.split() for ln in out.splitlines() if ln.startswith("ups1")]
        self.assertEqual(len(row), 1, out)
        self.assertEqual(row[0][2:], ["5", "1", "2"])


if __name__ == "__main__":
    unittest.main(verbosity=2)
