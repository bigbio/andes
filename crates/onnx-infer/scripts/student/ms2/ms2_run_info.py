#!/usr/bin/env python3
"""Instrument, activation and NCE of an mzML run, in andes-ms2 vocabulary.

Prints the instrument model and its family, and the activation / NCE counts over
the first --n MS2 scans, i.e. the ``--instrument``, ``--activation`` and ``--nce``
to pass to ``rescore_pin --model andes-ms2`` for that run.

Usage::

  python ms2_run_info.py run.mzML [--n 2000]
"""
from __future__ import annotations

import argparse
import collections
import os
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from ms2_common import INSTRUMENTS, instrument_class, parse_filter  # noqa: E402
from mzml_reader import instrument_name, read_scans  # noqa: E402


def main():
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("mzml")
    p.add_argument("--n", type=int, default=2000)
    a = p.parse_args()
    name = instrument_name(a.mzml)
    print(f"instrument {name!r} -> family {INSTRUMENTS[instrument_class(name)]}")
    c = collections.Counter()
    for k, s in enumerate(read_scans(a.mzml)):
        if k >= a.n:
            break
        c[parse_filter(s["filter"])] += 1
    for (act, nce), n in c.most_common():
        print(f"activation {act} NCE {nce:g}: {n} scans")


if __name__ == "__main__":
    main()
