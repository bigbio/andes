#!/usr/bin/env python3
"""Fill a MaxQuant 2.0.3 parameter template for a SILAC run.

    make_mqpar.py TEMPLATE OUT FASTA LABELS FILE...

TEMPLATE comes from `maxquant --create mqpar_template.xml`. LABELS lists the
channels separated by '|', with ';' between a channel's labels, e.g.
'|Arg6;Lys4|Arg10;Lys8' (triplex) or '|Lys8' (duplex). FILE are the mzXML inputs
(see mzml_to_mzxml.py). Sets 4 threads, Mono instead of .NET Core, an
identifier rule that takes the accession up to the first space, and at most 3
labelled residues.
"""
import sys, re
from xml.sax.saxutils import escape
tpl, out, fasta, labels, *files = sys.argv[1:]
s = open(tpl).read()
s = s.replace("<fastaFilePath>example.fasta</fastaFilePath>", f"<fastaFilePath>{escape(fasta)}</fastaFilePath>")
s = s.replace("<identifierParseRule>>(.*)</identifierParseRule>", "<identifierParseRule>>([^\\s]*)</identifierParseRule>")
s = s.replace("<useDotNetCore>True</useDotNetCore>", "<useDotNetCore>False</useDotNetCore>")
s = s.replace("<numThreads>1</numThreads>", "<numThreads>4</numThreads>")
def block(tag, items, kind):
    return f"<{tag}>\n" + "".join(f"      <{kind}>{escape(str(i))}</{kind}>\n" for i in items) + f"   </{tag}>"
s = re.sub(r"<filePaths>.*?</filePaths>", lambda _m: block("filePaths", files, "string"), s, flags=re.S)
s = re.sub(r"<experiments>.*?</experiments>", lambda _m: block("experiments", [f.rsplit('/',1)[-1].rsplit('.',1)[0] for f in files], "string"), s, flags=re.S)
s = re.sub(r"<fractions>.*?</fractions>", lambda _m: block("fractions", ["32767"] * len(files), "short"), s, flags=re.S)
s = re.sub(r"<ptms>.*?</ptms>", lambda _m: block("ptms", ["False"] * len(files), "boolean"), s, flags=re.S)
s = re.sub(r"<paramGroupIndices>.*?</paramGroupIndices>", lambda _m: block("paramGroupIndices", ["0"] * len(files), "int"), s, flags=re.S)
s = re.sub(r"<referenceChannel>.*?</referenceChannel>", lambda _m: block("referenceChannel", [""] * len(files), "string"), s, flags=re.S)
chans = labels.split("|")
s = re.sub(r"<labelMods>.*?</labelMods>", lambda _m: "<labelMods>\n" + "".join(f"            <string>{escape(c)}</string>\n" for c in chans) + "         </labelMods>", s, count=1, flags=re.S)
s = s.replace("<multiplicity>1</multiplicity>", f"<multiplicity>{len(chans)}</multiplicity>")
s = s.replace("<maxLabeledAa>0</maxLabeledAa>", "<maxLabeledAa>3</maxLabeledAa>")
open(out, "w").write(s)
