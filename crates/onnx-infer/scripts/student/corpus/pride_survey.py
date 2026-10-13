#!/usr/bin/env python3
"""Survey PRIDE Archive for training-corpus candidates (metadata only, no data download).

Subcommands
-----------
``fetch``     page through ``/v3/projects`` (100 per page; larger pageSize is silently
              clamped) into a compact JSONL. Submitter/PI contact details are dropped.
              Resumable; finishes by checking the unique-accession count against
              ``/v3/projects/count`` and refuses to report success on a short listing.
``classify``  derive the corpus strata of every project (analyser, fragmentation,
              labelling, enrichment, sample type, enzyme, organism, year, licence) from
              its CV terms + free text, write ``projects_classified.tsv`` and the
              per-stratum counts ``strata_counts.tsv``.
``files``     for a list of accessions, page through ``/v3/projects/<acc>/files`` (checked
              against the ``total_records`` header) and summarise which identification
              formats are present (mzIdentML, mzTab, PRIDE XML, pepXML, MaxQuant txt,
              FragPipe/Sage/DIA-NN outputs, peak lists, raw), i.e. whether the submission
              is usable without re-searching.

Every value derived from text is a *heuristic* and is flagged as such in the output
(``*_src`` columns say ``cv`` or ``text``). Shortlisted datasets must be verified
against the paper before use.

Usage::

  python pride_survey.py fetch --out projects.jsonl
  python pride_survey.py classify --projects projects.jsonl --out-dir survey/
  python pride_survey.py files --accessions acc.txt --out files_summary.tsv
"""
from __future__ import annotations

import argparse
import collections
import csv
import json
import os
import re
import sys
import time
import urllib.error
import urllib.request

API = "https://www.ebi.ac.uk/pride/ws/archive/v3"
UA = "andes-corpus-survey/1.0 (python-urllib)"
KEEP = [
    "accession", "title", "projectDescription", "sampleProcessingProtocol", "dataProcessingProtocol",
    "projectTags", "keywords", "submissionType", "license", "submissionDate", "publicationDate",
    "instruments", "softwares", "experimentTypes", "quantificationMethods", "sampleAttributes",
    "organisms", "organismParts", "diseases", "references", "identifiedPTMStrings", "totalFileDownloads",
]


def get(url: str, tries: int = 6):
    for k in range(tries):
        try:
            with urllib.request.urlopen(url, timeout=120) as r:
                return json.load(r), dict(r.headers)
        except (urllib.error.URLError, TimeoutError, json.JSONDecodeError) as e:  # transient
            if k == tries - 1:
                raise
            time.sleep(2 ** k)
            print(f"retry {k + 1} {url}: {e}", file=sys.stderr)


# ----------------------------------------------------------------------------- fetch
def cmd_fetch(a):
    total, _ = get(f"{API}/projects/count")
    seen = set()
    page = 0
    if os.path.exists(a.out):
        with open(a.out) as fh:
            for line in fh:
                seen.add(json.loads(line)["accession"])
        page = a.start_page if a.start_page is not None else len(seen) // 100  # duplicates are skipped
    with open(a.out, "a") as out:
        while True:
            rows, _ = get(f"{API}/projects?pageSize=100&page={page}")
            if not rows:
                break
            for r in rows:
                if r["accession"] in seen:
                    continue
                seen.add(r["accession"])
                out.write(json.dumps({k: r.get(k) for k in KEEP}) + "\n")
            page += 1
            if page % 20 == 0:
                print(f"page {page}: {len(seen)}/{total}", file=sys.stderr, flush=True)
            if len(rows) < 100:
                break
    print(f"{len(seen)} unique projects fetched; API count {total}")
    if len(seen) < total and not a.allow_short:
        sys.exit(f"SHORT LISTING: {total - len(seen)} projects missing; rerun (resumes) before using")


# ----------------------------------------------------------------------------- classify
ANALYSER_BY_INSTRUMENT = [
    # (regex on the instrument name, analyser of the MS2 scan most often used, family)
    (r"astral", "Astral", "Astral"),
    (r"timstof|tims tof", "TOF(timsTOF)", "timsTOF"),
    (r"tripletof|triple tof|qstar|5600|6600|7600|zenotof|x500|synapt|xevo|impact|maxis|compact|q-tof|qtof|6550|6545|6560|tof", "TOF", "QTOF"),
    (r"exploris", "Orbitrap", "Exploris"),
    (r"q exactive|q-exactive|qexactive|exactive", "Orbitrap", "QExactive"),
    (r"fusion|lumos|eclipse|ascend|tribrid", "Orbitrap/IT (tribrid)", "Tribrid"),
    (r"orbitrap velos|orbitrap elite|ltq orbitrap|orbitrap xl|orbitrap classic|orbitrap", "Orbitrap/IT (hybrid)", "LTQ-Orbitrap"),
    (r"ltq ft|ft-icr|fticr|solarix|apex", "FT-ICR/IT", "FT-ICR"),
    (r"ltq|velos pro|lcq|ion trap|amazon|hct|4000 qtrap|qtrap|lxq", "IonTrap", "IonTrap"),
]

TEXT_FLAGS = {
    # fragmentation
    "frag_etd": r"\b(etd|electron[- ]transfer dissociation)\b",
    "frag_ethcd": r"\bet[-]?hcd\b|\betcid\b",
    "frag_ecd": r"\becd\b|electron[- ]capture dissociation",
    "frag_uvpd": r"\buvpd\b|ultraviolet photodissociation",
    "frag_cid": r"\bcid\b|collision[- ]induced dissociation|\bcad\b",
    "frag_hcd": r"\bhcd\b|higher[- ]energy collision|beam[- ]type",
    "ms2_ion_trap": r"(ion trap|linear trap|\blit\b|\bitms\b|rapid scan|ion-trap|turbo scan).{0,40}(ms2|ms/ms|ms\/ms|fragment|detect|acquir)|(ms2|ms/ms|fragment).{0,60}(in the|linear) ion trap",
    "ms2_orbitrap": r"(ms2|ms/ms|fragment).{0,80}orbitrap|orbitrap.{0,40}(ms2|ms/ms)|resolution of 1[57],?[05]00|17,?500|15,?000 (resolution|at m/z)",
    # labelling
    "lab_tmt": r"\btmt|tandem mass tag|tmtpro|tmt10|tmt11|tmt16|tmt18",
    "lab_tmtpro": r"tmt ?pro|16-?plex tmt|tmt ?16|tmt ?18|18-?plex",
    "lab_itraq": r"\bitraq",
    "lab_silac": r"\bsilac\b|stable isotope labell?ing (by|with) amino",
    "lab_dimethyl": r"dimethyl labell?ing|reductive dimethyl",
    # enrichment / PTM focus
    "enr_phospho": r"phosphoproteom|phosphopeptide|phospho-?enrich|\btio2\b|\bimac\b|fe-nta|ti-imac|ti4\+|phosphosite|phospho-?tyrosine (enrich|immuno)|py(100|1000)",
    "enr_glyco": r"glycoproteom|glycopeptide|glyco-?site|n-glyco(syl|pep|form)|o-glyco(syl|pep|form)|intact glyco|glycoform|sialoglyco|\bhilic\b.{0,40}enrich|lectin.{0,40}enrich",
    # note: no bare "kegg"-like pattern - "k-?e-?gg" matched every KEGG pathway analysis
    "enr_ubiquitin": r"ubiquitinom|ubiquitylom|ubiquitin remnant|diglycine|di-glycine|gly-gly|k-[eε]-gg|\bgg[- ]remnant|ubiscan|ptmscan ubiquitin|ubiquitin.{0,30}enrich",
    "enr_acetyl": r"acetyl[- ]?(ome|lysine|ation enrich)|acetyl-lysine|lysine acetylation",
    "enr_methyl": r"methyl[- ]?(ome|arginine|lysine)|arginine methylation",
    "enr_succinyl_acyl": r"succinyl|crotonyl|malonyl|lactyl|glutaryl|hydroxybutyryl|propionyl",
    "enr_citrullin": r"citrullin",
    "enr_sumo": r"sumo",
    "enr_crosslink": r"cross-?link|xl-ms|\bdss\b|\bdsso\b|\bbs3\b",
    "enr_redox": r"cysteine oxidation|redox|s-nitrosyl|persulfid|sulfenyl",
    "enr_adpr": r"adp-ribosyl",
    "enr_histone": r"histone",
    "hla_class1": r"hla[- ]?(class[- ]?)?(i\b|a\b|b\b|c\b|-a|-b|-c)|mhc[- ]?(class[- ]?)?i\b|w6/32|hla class i\b|class i hla",
    "hla_class2": r"hla[- ]?(class[- ]?)?ii\b|hla-dr|hla-dq|hla-dp|mhc[- ]?(class[- ]?)?ii\b|l243|class ii",
    # strong evidence only: "neoantigen" / "antigen presentation" alone flagged proteogenomics papers
    "immunopeptidome": r"immunopeptidom|hla[- ]?(class[- ]?i{1,2}[- ])?peptidom|mhc[- ]?(class[- ]?i{1,2}[- ])?peptidom|hla ligand|mhc ligand|ligandom|hla-bound|mhc-bound|hla-associated peptid|mhc-associated peptid|eluted (hla|mhc)|(hla|mhc).{0,40}(immunoprecipitat|immunoaffinity|elut)",
    # enzyme
    "enz_trypsin": r"tryps",
    "enz_lysc": r"lys-?c\b|lysc",
    "enz_gluc": r"\bglu-?c\b|endoproteinase glu|v8 protease",
    "enz_aspn": r"asp-?n\b|aspn",
    "enz_chymo": r"chymotryp",
    "enz_argc": r"arg-?c\b|argc",
    "enz_lysn": r"lys-?n\b|lysn",
    "enz_pepsin_ntryp": r"pepsin|elastase|proteinase k|alpha-lytic|non-?tryptic|unspecific|no enzyme|nonspecific|semi-?specific",
    "dia": r"\bdia\b|data[- ]independent|swath|diapasef",
    "synthetic": r"synthetic peptide|proteometools|peptide library",
    "plasma": r"plasma|serum",
    "single_cell": r"single[- ]cell",
    "metaproteomics": r"metaproteom|microbiome",
    "crosslink_ms": r"cross-?linking mass spec",
}
_TEXT_RE = {k: re.compile(v, re.I) for k, v in TEXT_FLAGS.items()}

ID_SOFTWARE = {
    "maxquant": "MaxQuant", "andromeda": "MaxQuant", "mascot": "Mascot", "sequest": "SEQUEST",
    "proteome discoverer": "ProteomeDiscoverer", "msfragger": "MSFragger", "fragpipe": "MSFragger",
    "ms-gf": "MS-GF+", "msgf": "MS-GF+", "comet": "Comet", "x!tandem": "X!Tandem", "xtandem": "X!Tandem",
    "peaks": "PEAKS", "spectrum mill": "SpectrumMill", "byonic": "Byonic", "pfind": "pFind",
    "sage": "Sage", "dia-nn": "DIA-NN", "diann": "DIA-NN", "spectronaut": "Spectronaut",
    "scaffold": "Scaffold", "percolator": "Percolator", "skyline": "Skyline", "omssa": "OMSSA",
    "chimerys": "CHIMERYS", "metamorpheus": "MetaMorpheus", "msamanda": "MSAmanda", "ms amanda": "MSAmanda",
    "netmhc": "NetMHC", "pglyco": "pGlyco", "protein pilot": "ProteinPilot", "proteinpilot": "ProteinPilot",
    "philosopher": "MSFragger", "quantms": "quantms", "openms": "OpenMS",
}


def names(lst):
    return [x.get("name", "") for x in (lst or []) if isinstance(x, dict)]


def classify_one(r: dict) -> dict:
    instr = names(r.get("instruments"))
    text = " ".join(str(r.get(k) or "") for k in (
        "title", "projectDescription", "sampleProcessingProtocol", "dataProcessingProtocol"))
    text += " " + " ".join(r.get("keywords") or []) + " " + " ".join(r.get("projectTags") or [])
    t = text.lower()
    # instrument names ("Orbitrap Fusion ETD", "Orbitrap Elite ETD") are not evidence of ETD spectra
    t = re.sub(r"(fusion|lumos|tribrid|elite|velos|orbitrap|ltq)\s*(lumos\s*)?(tribrid\s*)?etd", " ", t)
    flags = {k: bool(rx.search(t)) for k, rx in _TEXT_RE.items()}

    analysers, families = set(), set()
    for name in instr:
        n = name.lower()
        for rx, an, fam in ANALYSER_BY_INSTRUMENT:
            if re.search(rx, n):
                analysers.add(an)
                families.add(fam)
                break
        else:
            families.add("other")
    ptms = [x.get("accession", "") + ":" + x.get("name", "") for x in (r.get("identifiedPTMStrings") or [])]
    ptm_names = " ".join(ptms).lower()
    quant = " ".join(names(r.get("quantificationMethods"))).lower()
    sw = sorted({v for k, v in ID_SOFTWARE.items() if k in (" ".join(names(r.get("softwares"))) + " " +
                                                         str(r.get("dataProcessingProtocol") or "")).lower()})
    orgs = names(r.get("organisms"))
    # Labelling: CV PTM / quantification first (cv), then text.
    lab = []
    if "tmt" in ptm_names or "tmt" in quant:
        lab.append("TMT")
    if "itraq" in ptm_names or "itraq" in quant:
        lab.append("iTRAQ")
    if "silac" in quant or "label:13c(6)" in ptm_names:
        lab.append("SILAC")
    lab_src = "cv" if lab else "text"
    if "TMT" in lab and (flags["lab_tmtpro"] or "unimod:2016" in ptm_names):
        lab.append("TMTpro")
    if not lab:
        lab = [n for n, f in (("TMT", "lab_tmt"), ("iTRAQ", "lab_itraq"), ("SILAC", "lab_silac"),
                              ("dimethyl", "lab_dimethyl")) if flags[f]]
        if "TMT" in lab and flags["lab_tmtpro"]:
            lab.append("TMTpro")
    enr = []
    for name, f, cv in (("phospho", "enr_phospho", "unimod:21:"), ("glyco", "enr_glyco", "hexnac"),
                        ("ubiquitin", "enr_ubiquitin", "unimod:121:"), ("acetyl", "enr_acetyl", None),
                        ("methyl", "enr_methyl", None), ("acyl", "enr_succinyl_acyl", None),
                        ("citrullination", "enr_citrullin", None), ("crosslink", "enr_crosslink", None),
                        ("ADP-ribosyl", "enr_adpr", None), ("histone", "enr_histone", None)):
        # a CV PTM alone is enough only for phospho / GG (they are rarely searched without enrichment)
        if flags[f] or (cv and cv in ptm_names and name in ("phospho", "ubiquitin") and flags.get("enr_" + name, True)) \
                or (cv and cv in ptm_names and name == "glyco"):
            enr.append(name)
    hla = []
    if flags["immunopeptidome"]:
        hla.append("HLA")
        if flags["hla_class1"]:
            hla.append("HLA-I")
        if flags["hla_class2"]:
            hla.append("HLA-II")
    frag = [n for n, f in (("ETD", "frag_etd"), ("EThcD", "frag_ethcd"), ("ECD", "frag_ecd"), ("UVPD", "frag_uvpd"),
                           ("CID", "frag_cid"), ("HCD", "frag_hcd")) if flags[f]]
    ms2 = []
    if flags["ms2_ion_trap"]:
        ms2.append("IT")
    if flags["ms2_orbitrap"]:
        ms2.append("OT")
    enz = [n for n, f in (("trypsin", "enz_trypsin"), ("LysC", "enz_lysc"), ("GluC", "enz_gluc"), ("AspN", "enz_aspn"),
                          ("chymotrypsin", "enz_chymo"), ("ArgC", "enz_argc"), ("LysN", "enz_lysn"),
                          ("nonspecific", "enz_pepsin_ntryp")) if flags[f]]
    year = (r.get("publicationDate") or r.get("submissionDate") or "")[:4]
    refs = r.get("references") or []
    return {
        "accession": r["accession"],
        "year": year,
        "submission_type": r.get("submissionType") or "",
        "licence": r.get("license") or "",
        "organisms": "|".join(orgs),
        "human": int(any("9606" in (x.get("accession") or "") for x in (r.get("organisms") or []))),
        "instruments": "|".join(instr),
        "instrument_family": "|".join(sorted(families)),
        "analyser": "|".join(sorted(analysers)),
        "fragmentation_text": "|".join(frag),
        "ms2_analyser_text": "|".join(ms2),
        "labelling": "|".join(lab) or "LFQ?",
        "labelling_src": lab_src,
        "enrichment": "|".join(enr),
        "hla": "|".join(hla),
        "enzymes_text": "|".join(enz),
        "dia": int(flags["dia"]),
        "synthetic": int(flags["synthetic"]),
        "plasma": int(flags["plasma"]),
        "single_cell": int(flags["single_cell"]),
        "metaproteomics": int(flags["metaproteomics"]),
        "software": "|".join(sw),
        "ptms_cv": "|".join(ptms),
        "n_ptms_cv": len(ptms),
        "experiment_types": "|".join(names(r.get("experimentTypes"))),
        "pubmed": "|".join(str(x.get("pubmedID") or "") for x in refs if x.get("pubmedID")),
        "doi": "|".join(str(x.get("doi") or "") for x in refs if x.get("doi")),
        "downloads": r.get("totalFileDownloads") or 0,
        "title": (r.get("title") or "").replace("\t", " ").replace("\n", " ")[:200],
    }


def cmd_classify(a):
    os.makedirs(a.out_dir, exist_ok=True)
    rows = []
    with open(a.projects) as fh:
        for line in fh:
            rows.append(classify_one(json.loads(line)))
    cols = list(rows[0])
    with open(os.path.join(a.out_dir, "projects_classified.tsv"), "w", newline="") as fh:
        w = csv.DictWriter(fh, cols, delimiter="\t")
        w.writeheader()
        w.writerows(rows)
    # strata counts: overall and restricted to COMPLETE submissions (identifications inside)
    def explode(v):
        return [x for x in v.split("|") if x] or ["(none)"]

    counts = collections.Counter()
    for r in rows:
        dda = not r["dia"]
        for subset, ok in (("all", True), ("COMPLETE", r["submission_type"] == "COMPLETE"),
                           ("COMPLETE_DDA", r["submission_type"] == "COMPLETE" and dda)):
            if not ok:
                continue
            counts[(subset, "total", "all")] += 1
            counts[(subset, "submission_type", r["submission_type"])] += 1
            counts[(subset, "year", r["year"])] += 1
            counts[(subset, "human", str(r["human"]))] += 1
            for k in ("analyser", "instrument_family", "fragmentation_text", "ms2_analyser_text", "labelling",
                      "enrichment", "hla", "enzymes_text", "licence"):
                for v in explode(r[k]):
                    counts[(subset, k, v)] += 1
            for k in ("dia", "synthetic", "plasma", "single_cell", "metaproteomics"):
                counts[(subset, k, str(r[k]))] += 1
            for p in explode(r["ptms_cv"]):
                counts[(subset, "ptm_cv", p)] += 1
    with open(os.path.join(a.out_dir, "strata_counts.tsv"), "w") as fh:
        fh.write("subset\tfield\tvalue\tn_projects\n")
        for (s, f, v), n in sorted(counts.items(), key=lambda kv: (kv[0][0], kv[0][1], -kv[1])):
            fh.write(f"{s}\t{f}\t{v}\t{n}\n")
    print(f"{len(rows)} projects classified -> {a.out_dir}")


# ----------------------------------------------------------------------------- files
FORMAT_RULES = [
    ("mzIdentML", r"\.mzid(\.gz|\.zip)?$|mzidentml[^/]*\.xml(\.gz)?$"),
    ("mzTab", r"\.mztab(\.gz|\.zip)?$"),
    # PRIDE XML files are PRIDE_Exp_*.xml; a bare ".xml" is usually mqpar.xml or a Mascot export
    ("PRIDE_XML", r"pride_exp[^/]*\.xml(\.gz)?$|\.pride\.xml(\.gz)?$"),
    ("Mascot_xml_dat", r"mascot[^/]*\.xml$|\.dat$|f\d{6}\.xml$"),
    ("pepXML", r"\.pep\.xml(\.gz)?$|\.pepxml(\.gz)?$"),
    ("MaxQuant", r"msms\.txt|evidence\.txt|txt\.zip$|combined.*\.zip$|maxquant|_txt\.zip$|/txt/"),
    ("FragPipe", r"psm\.tsv$|fragpipe|\.pin$"),
    ("Sage", r"results\.sage\.tsv"),
    ("PD_msf", r"\.msf$|\.pdresult$"),
    ("mzML", r"\.mzml(\.gz|\.zip)?$"),
    ("MGF", r"\.mgf(\.gz|\.zip)?$"),
    ("MS2/apl", r"\.ms2$|\.apl$|\.mzxml(\.gz)?$|\.dta"),
    ("thermo_raw", r"\.raw(\.zip)?$"),
    ("bruker_d", r"\.d\.zip$|\.d\.tar|\.tdf$|\.baf$|\.d$"),
    ("sciex_wiff", r"\.wiff(\.scan)?$|\.wiff2$"),
    ("dia_lib", r"\.speclib$"),
    # generic search-result archives ("search.zip", "*_txt.zip", ProteomeTools per-pool zips) are
    # usually MaxQuant txt folders; verify with peek-zip before relying on it
    ("search_zip?", r"search[^/]*\.zip$|txt[^/]*\.zip$|[-_](un)?specific\.zip$|[-_](semi)?tryptic\.zip$"),
    ("sdrf", r"sdrf"),
]
_FMT = [(k, re.compile(v, re.I)) for k, v in FORMAT_RULES]


def file_formats(name: str):
    return [k for k, rx in _FMT if rx.search(name)]


def list_files(acc: str):
    out, page, total = [], 0, None
    while True:
        rows, hdr = get(f"{API}/projects/{acc}/files?pageSize=100&page={page}")
        if total is None:
            total = int({k.lower(): v for k, v in hdr.items()}.get("total_records", -1))
        out += rows or []
        if not rows or len(rows) < 100:
            break
        page += 1
    return out, total


FILE_COLS = ["accession", "n_files", "total_records", "listing_complete", "categories", "formats",
             "n_raw", "raw_gb", "n_peak", "n_result", "n_search", "usable_without_research", "example_result"]
IDENT = {"mzIdentML", "mzTab", "PRIDE_XML", "pepXML", "MaxQuant", "FragPipe", "Sage", "PD_msf", "Mascot_xml_dat",
         "search_zip?"}


def summarise_files(acc, files, total):
    """files: [(name, category, size)] -> one FILE_COLS row."""
    cats = collections.Counter(c or "?" for _, c, _ in files)
    fmts = collections.Counter(k for n, _, _ in files for k in file_formats(n or ""))
    raw = [f for f in files if f[1] == "RAW"]
    peaks = fmts["mzML"] + fmts["MGF"] + fmts["MS2/apl"]
    has_id = sorted(k for k in IDENT if fmts[k])
    if fmts["mzIdentML"] or fmts["mzTab"] or fmts["PRIDE_XML"]:
        usable = "standard_ids" + ("+peaks" if peaks else "+raw_only" if raw else "")
    elif has_id:
        usable = "engine_ids(" + ",".join(has_id) + ")" + ("+peaks" if peaks else "+raw" if raw else "")
    else:
        usable = "research_needed"
    ex = next((n for n, _, _ in files if set(file_formats(n or "")) & IDENT), "")
    return [acc, len(files), total, int(total == len(files)),
            ",".join(f"{k}:{v}" for k, v in cats.most_common()), ",".join(f"{k}:{v}" for k, v in fmts.most_common()),
            len(raw), round(sum(sz or 0 for _, _, sz in raw) / 1e9, 1),
            cats.get("PEAK", 0), cats.get("RESULT", 0), cats.get("SEARCH", 0), usable, ex]


def cmd_files(a):
    if a.from_listing:  # re-summarise saved listings (after a format-rule change) without the API
        rows = {}
        for line in open(a.from_listing):
            d = json.loads(line)
            rows[d["accession"]] = summarise_files(d["accession"], d["files"], len(d["files"]))
        with open(a.out, "w") as fh:
            fh.write("\t".join(FILE_COLS) + "\n")
            for r in rows.values():
                fh.write("\t".join(map(str, r)) + "\n")
        print(f"{len(rows)} projects re-summarised -> {a.out}")
        return
    accs = [x.strip() for x in open(a.accessions) if x.strip() and not x.startswith("#")]
    done = set()
    if os.path.exists(a.out):
        done = {r["accession"] for r in csv.DictReader(open(a.out), delimiter="\t")}
    new = not os.path.exists(a.out)
    with open(a.out, "a") as fh:
        if new:
            fh.write("\t".join(FILE_COLS) + "\n")
        for acc in accs:
            if acc in done:
                continue
            try:
                files, total = list_files(acc)
            except Exception as e:  # noqa: BLE001 - keep going; record the failure
                fh.write(f"{acc}\t-1\t-1\t0\tERROR {e}\t\t0\t0\t0\t0\t0\tunknown\t\n")
                continue
            slim = [[f.get("fileName"), (f.get("fileCategory") or {}).get("value"), f.get("fileSizeBytes")] for f in files]
            with open(a.out + ".listing.jsonl", "a") as lf:
                lf.write(json.dumps({"accession": acc, "files": slim}) + "\n")
            fh.write("\t".join(map(str, summarise_files(acc, slim, total))) + "\n")
            fh.flush()
            time.sleep(a.sleep)
    print(f"files summary -> {a.out}")


# ----------------------------------------------------------------------------- rank
# Gap value relative to the corpus we already hold (qpx: human, tryptic, Orbitrap HCD/CID,
# LFQ + TMT10, 4 phospho sets; MsNet flats: PXD009875 HeLa FTMS-HCD + ITMS). See CORPUS.md.
GAP_WEIGHTS = {
    "hla": 6.0, "hla2_extra": 1.0, "nontryptic": 3.0, "etd_family": 3.0, "timsTOF": 4.0, "Astral": 4.0,
    "TOF": 2.0, "IonTrap": 2.0, "ms2_IT": 2.0, "iTRAQ": 2.0, "TMTpro": 2.0, "TMT": 0.5, "SILAC": 0.5,
    "ptm_rare": 3.0, "histone": 2.0, "glyco": 1.0, "phospho": 1.0, "synthetic": 3.0, "nonhuman": 0.5,
    "many_ptms": 1.0, "complete": 3.0, "parseable_engine": 1.0, "open_licence": 0.5, "recent": 0.5,
    "dia": -6.0, "crosslink": -3.0, "metaproteomics": -1.0,
}
PARSEABLE = {"MaxQuant", "MSFragger", "Sage", "Comet", "MS-GF+", "Mascot", "PEAKS", "ProteomeDiscoverer", "pFind",
             "MetaMorpheus", "Byonic", "X!Tandem", "quantms", "OpenMS"}


def gap_score(r: dict) -> tuple[float, list]:
    g, why = 0.0, []

    def add(k, cond=True):
        nonlocal g
        if cond:
            g += GAP_WEIGHTS[k]
            why.append(k)
    enr = r["enrichment"].split("|")
    add("hla", bool(r["hla"]))
    add("hla2_extra", "HLA-II" in r["hla"])
    enz = set(r["enzymes_text"].split("|")) - {""}
    add("nontryptic", bool(enz - {"trypsin", "LysC"}) or bool(r["hla"]))
    add("etd_family", bool(set(r["fragmentation_text"].split("|")) & {"ETD", "EThcD", "ECD", "UVPD"}))
    fam = r["instrument_family"].split("|")
    add("timsTOF", "timsTOF" in fam)
    add("Astral", "Astral" in fam)
    add("TOF", "QTOF" in fam)
    add("IonTrap", "IonTrap" in fam)
    add("ms2_IT", "IT" in r["ms2_analyser_text"])
    add("iTRAQ", "iTRAQ" in r["labelling"])
    add("TMTpro", "TMTpro" in r["labelling"])
    add("TMT", "TMT" in r["labelling"])
    add("SILAC", "SILAC" in r["labelling"])
    add("ptm_rare", bool(set(enr) & {"ubiquitin", "acetyl", "methyl", "acyl", "citrullination", "ADP-ribosyl"}))
    add("histone", "histone" in enr)
    add("glyco", "glyco" in enr)
    add("phospho", "phospho" in enr)
    add("synthetic", r["synthetic"] == "1")
    add("nonhuman", r["human"] == "0")
    add("many_ptms", int(r["n_ptms_cv"] or 0) >= 6)
    add("complete", r["submission_type"] == "COMPLETE")
    add("parseable_engine", r["submission_type"] != "COMPLETE" and bool(set(r["software"].split("|")) & PARSEABLE))
    add("open_licence", "Creative Commons" in r["licence"])
    add("recent", (r["year"] or "0") >= "2018")
    add("dia", r["dia"] == "1")
    add("crosslink", "crosslink" in enr)
    add("metaproteomics", r["metaproteomics"] == "1")
    return g, why


def cmd_rank(a):
    rows = list(csv.DictReader(open(a.classified), delimiter="\t"))
    excluded = {"PXD070049", "PXD001819", "PXD007683", "PXD007653", "PXD005553", "PXD002047"} | set(a.exclude)
    out = []
    for r in rows:
        if r["accession"] in excluded:
            continue
        g, why = gap_score(r)
        out.append((g, r["accession"], ",".join(why)))
    out.sort(reverse=True)
    with open(a.out, "w") as fh:
        fh.write("accession\tgap_score\treasons\n")
        for g, acc, why in out[: a.top]:
            fh.write(f"{acc}\t{g:.1f}\t{why}\n")
    print(f"top {a.top} of {len(out)} -> {a.out}")


# ----------------------------------------------------------------------------- summary
def cmd_summary(a):
    """Headline stratum counts (projects) for all / COMPLETE / COMPLETE-and-not-DIA submissions."""
    import pandas as pd

    d = pd.read_csv(a.classified, sep="\t", dtype=str).fillna("")

    def has(col, v):
        return lambda x: x[col].str.split("|").apply(lambda l: v in l)

    def con(col, pat):
        return lambda x: x[col].str.contains(pat, regex=True)
    strata = {
        "total": lambda x: pd.Series(True, index=x.index), "human": lambda x: x.human == "1",
        "DIA (text)": lambda x: x.dia == "1",
        "instr: Q Exactive family": has("instrument_family", "QExactive"),
        "instr: Exploris": has("instrument_family", "Exploris"),
        "instr: Tribrid (Fusion/Lumos/Eclipse/Ascend)": has("instrument_family", "Tribrid"),
        "instr: LTQ-Orbitrap hybrid (Velos/Elite/XL)": has("instrument_family", "LTQ-Orbitrap"),
        "instr: Astral": has("instrument_family", "Astral"), "instr: timsTOF": has("instrument_family", "timsTOF"),
        "instr: other Q-TOF (SCIEX/Waters/Agilent/Bruker)": has("instrument_family", "QTOF"),
        "instr: ion trap only": has("instrument_family", "IonTrap"), "instr: FT-ICR": has("instrument_family", "FT-ICR"),
        "MS2 in ion trap (text)": con("ms2_analyser_text", "IT"),
        "frag: HCD (text)": has("fragmentation_text", "HCD"), "frag: CID (text)": has("fragmentation_text", "CID"),
        "frag: ETD (text)": has("fragmentation_text", "ETD"), "frag: EThcD (text)": has("fragmentation_text", "EThcD"),
        "frag: ECD/UVPD (text)": con("fragmentation_text", "ECD|UVPD"),
        "lab: TMT (any)": has("labelling", "TMT"), "lab: TMTpro": has("labelling", "TMTpro"),
        "lab: iTRAQ": has("labelling", "iTRAQ"), "lab: SILAC": has("labelling", "SILAC"),
        "enr: phospho": has("enrichment", "phospho"), "enr: glyco": has("enrichment", "glyco"),
        "enr: ubiquitin/GG": has("enrichment", "ubiquitin"), "enr: acetyl": has("enrichment", "acetyl"),
        "enr: methyl": has("enrichment", "methyl"), "enr: acyl (succ/cro/mal/lac/...)": has("enrichment", "acyl"),
        "enr: citrullination": has("enrichment", "citrullination"), "enr: ADP-ribosyl": has("enrichment", "ADP-ribosyl"),
        "enr: histone": has("enrichment", "histone"), "enr: crosslink": has("enrichment", "crosslink"),
        "HLA / immunopeptidome": lambda x: x.hla != "", "HLA-I": con("hla", r"HLA-I(?:\||$)"),
        "HLA-II": con("hla", "HLA-II"),
        "enz: GluC": has("enzymes_text", "GluC"), "enz: AspN": has("enzymes_text", "AspN"),
        "enz: chymotrypsin": has("enzymes_text", "chymotrypsin"), "enz: ArgC": has("enzymes_text", "ArgC"),
        "enz: LysN": has("enzymes_text", "LysN"), "enz: nonspecific (text)": has("enzymes_text", "nonspecific"),
        "synthetic peptides": lambda x: x.synthetic == "1", "plasma/serum": lambda x: x.plasma == "1",
        "licence CC0": con("licence", "CC0"),
    }
    sub = {"all": d, "COMPLETE": d[d.submission_type == "COMPLETE"],
           "COMPLETE_DDA": d[(d.submission_type == "COMPLETE") & (d.dia != "1")]}
    rows = [[k] + [int(f(v).sum()) for v in sub.values()] for k, f in strata.items()]
    t = pd.DataFrame(rows, columns=["stratum"] + list(sub))
    t.to_csv(a.out, sep="\t", index=False)
    print(t.to_string(index=False))
    print(pd.crosstab(d.year, d.submission_type).to_string())


# ----------------------------------------------------------------------------- peek-zip
def _range(url, start, end=None):
    h = {"Range": f"bytes={start}-{'' if end is None else end}", "User-Agent": UA}
    with urllib.request.urlopen(urllib.request.Request(url, headers=h), timeout=120) as r:
        return r.read()


def zip_members(url):
    """Member names + sizes of a remote zip, read from its central directory with two HTTP
    Range requests (no download). Handles zip64."""
    import struct

    with urllib.request.urlopen(urllib.request.Request(url, method="HEAD", headers={"User-Agent": UA}), timeout=120) as r:
        size = int(r.headers["Content-Length"])
    tail = _range(url, max(0, size - 65536))
    i = tail.rfind(b"PK\x05\x06")
    cd_size, cd_off = struct.unpack("<II", tail[i + 12:i + 20])
    if cd_off == 0xFFFFFFFF:
        j = tail.rfind(b"PK\x06\x06")
        cd_size, cd_off = struct.unpack("<QQ", tail[j + 40:j + 56])
    cd = _range(url, cd_off, cd_off + cd_size - 1)
    out, p = [], 0
    while cd[p:p + 4] == b"PK\x01\x02":
        method = struct.unpack("<H", cd[p + 10:p + 12])[0]
        csize, usize = struct.unpack("<II", cd[p + 20:p + 28])
        n, e, c = struct.unpack("<HHH", cd[p + 28:p + 34])
        off = struct.unpack("<I", cd[p + 42:p + 46])[0]
        extra = cd[p + 46 + n:p + 46 + n + e]
        q = 0
        while q + 4 <= len(extra):  # zip64 extended information
            tag, ln = struct.unpack("<HH", extra[q:q + 4])
            if tag == 1:
                vals, k = extra[q + 4:q + 4 + ln], 0
                if usize == 0xFFFFFFFF:
                    usize = struct.unpack("<Q", vals[k:k + 8])[0]
                    k += 8
                if csize == 0xFFFFFFFF:
                    csize = struct.unpack("<Q", vals[k:k + 8])[0]
                    k += 8
                if off == 0xFFFFFFFF:
                    off = struct.unpack("<Q", vals[k:k + 8])[0]
            q += 4 + ln
        out.append((cd[p + 46:p + 46 + n].decode(errors="replace"), usize, csize, off, method))
        p += 46 + n + e + c
    return size, out


def zip_extract(url, member, dest):
    """Fetch one member of a remote zip by HTTP Range (deflate or stored) into dest."""
    import struct
    import zlib

    _, members = zip_members(url)
    hit = [m for m in members if m[0] == member]
    if not hit:
        raise SystemExit(f"{member} not in zip")
    _, usize, csize, off, method = hit[0]
    head = _range(url, off, off + 29)
    n, e = struct.unpack("<HH", head[26:30])
    start = off + 30 + n + e
    d = zlib.decompressobj(-15) if method == 8 else None
    with open(dest, "wb") as fh:
        pos, step = start, 64 << 20
        while pos < start + csize:
            chunk = _range(url, pos, min(pos + step, start + csize) - 1)
            fh.write(d.decompress(chunk) if d else chunk)
            pos += len(chunk)
        if d:
            fh.write(d.flush())
    return usize


def cmd_peek_zip(a):
    """Verify what a submission's zip holds (e.g. MaxQuant txt/msms.txt) without downloading it."""
    page = 0
    while True:
        rows, _ = get(f"{API}/projects/{a.accession}/files?pageSize=100&page={page}")
        hit = [x for x in rows if x["fileName"] == a.file]
        if hit or len(rows) < 100:
            break
        page += 1
    if not hit:
        sys.exit(f"{a.file} not in {a.accession}")
    loc = [x["value"] for x in hit[0]["publicFileLocations"] if x["value"].startswith("ftp://")][0]
    url = loc.replace("ftp://ftp.pride.ebi.ac.uk", "https://ftp.pride.ebi.ac.uk")
    size, members = zip_members(url)
    print(f"{url}\t{size}")
    for m in members:
        print(f"{m[1]}\t{m[0]}")
    if a.extract:
        n = zip_extract(url, a.extract, a.out or os.path.basename(a.extract))
        print(f"extracted {a.extract} ({n} bytes) -> {a.out or os.path.basename(a.extract)}")


def main():
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    sub = p.add_subparsers(dest="cmd", required=True)
    f = sub.add_parser("fetch")
    f.add_argument("--out", required=True)
    f.add_argument("--start-page", type=int, default=None,
                   help="page to start from (0 = full re-scan that only appends missing projects; "
                        "the listing order is not stable across a long scan, so a second pass is normal)")
    f.add_argument("--allow-short", action="store_true",
                   help="accept a listing shorter than /projects/count after a full re-scan found nothing new "
                        "(2026-10-13: two passes gave 41,465 of 42,027; sampled unlisted PXDs return 404)")
    c = sub.add_parser("classify")
    c.add_argument("--projects", required=True)
    c.add_argument("--out-dir", required=True)
    fl = sub.add_parser("files")
    fl.add_argument("--accessions", help="one accession per line")
    fl.add_argument("--from-listing", help="re-summarise a saved <out>.listing.jsonl instead of calling the API")
    fl.add_argument("--out", required=True)
    fl.add_argument("--sleep", type=float, default=0.2)
    rk = sub.add_parser("rank")
    rk.add_argument("--classified", required=True)
    rk.add_argument("--out", required=True)
    rk.add_argument("--top", type=int, default=800)
    rk.add_argument("--exclude", nargs="*", default=[])
    sm = sub.add_parser("summary")
    sm.add_argument("--classified", required=True)
    sm.add_argument("--out", required=True)
    pz = sub.add_parser("peek-zip")
    pz.add_argument("--accession", required=True)
    pz.add_argument("--file", required=True)
    pz.add_argument("--extract", help="member to fetch by HTTP Range (e.g. txt/msms.txt)")
    pz.add_argument("--out", help="destination of --extract")
    a = p.parse_args()
    {"fetch": cmd_fetch, "classify": cmd_classify, "files": cmd_files, "rank": cmd_rank,
     "summary": cmd_summary, "peek-zip": cmd_peek_zip}[a.cmd](a)


if __name__ == "__main__":
    main()
