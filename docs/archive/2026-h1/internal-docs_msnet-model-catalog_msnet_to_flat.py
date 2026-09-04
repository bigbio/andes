#!/usr/bin/env python3
"""MSNet PSM parquet -> flat Andes training parquet (quality-gated).
Schema-adaptive; filters to confident, well-fragmented, non-decoy HCD tryptic PSMs;
resolves UNIMOD mods -> mass deltas; isotope-corrects the precursor; emits the flat
schema the `andes train-from-msnet` reader expects.
Args: <in_parquet> <out_parquet> [pep_max=0.001] [sample] [min_consensus=0] [min_matched=6] [min_explained=0.50]
"""
import os, sys, duckdb, pyarrow as pa, pyarrow.parquet as pq

UNIMOD = {
    "Carbamidomethyl": 57.021464, "Oxidation": 15.994915, "Acetyl": 42.010565,
    "Phospho": 79.966331, "Deamidated": 0.984016, "Methyl": 14.015650,
    "Dimethyl": 28.031300, "Trimethyl": 42.046950, "Gln->pyro-Glu": -17.026549,
    "Glu->pyro-Glu": -18.010565, "Ammonia-loss": -17.026549, "Carbamyl": 43.005814,
    "Formyl": 27.994915, "Dehydrated": -18.010565, "Nitro": 44.985078,
    "TMT6plex": 229.162932, "TMTpro": 304.207146, "iTRAQ4plex": 144.102063, "iTRAQ8plex": 304.205360,
}
_UNIMOD_LC = {k.lower(): v for k, v in UNIMOD.items()}
NTERM = {"Acetyl", "Carbamyl", "Formyl", "TMT6plex", "TMTpro", "iTRAQ4plex", "iTRAQ8plex",
         "Gln->pyro-Glu", "Glu->pyro-Glu", "Ammonia-loss"}


# --- monoisotopic residue masses + b/y fragment-offset recalibration ---
AA = {"G":57.021464,"A":71.037114,"S":87.032028,"P":97.052764,"V":99.068414,
 "T":101.047678,"C":103.009185,"L":113.084064,"I":113.084064,"N":114.042927,
 "D":115.026943,"Q":128.058578,"K":128.094963,"E":129.042593,"M":131.040485,
 "H":137.058912,"F":147.068414,"R":156.101111,"Y":163.063329,"W":186.079313}
PROTON = 1.007276; H2O = 18.010565
import re as _re
def spectrum_offset_ppm(seq, res_pos, res_delta, nterm, cterm, mz, iontypes, charges):
    """Median ppm offset (observed - theoretical) over annotated b/y ions; None if too few."""
    n=len(seq); rmod=[0.0]*(n+1)
    for p,d in zip(res_pos,res_delta):
        if 1<=p<=n: rmod[p]+=d
    cum=[0.0]*(n+1)
    for i in range(1,n+1):
        m=AA.get(seq[i-1]); 
        if m is None: return None
        cum[i]=cum[i-1]+m+rmod[i]
    diffs=[]
    for t,obs,z in zip(iontypes or [], mz or [], charges or []):
        if not t: continue
        mm=_re.match(r'^([by])(\d+)', t)
        if not mm: continue
        kind=mm.group(1); idx=int(mm.group(2))
        zz=int(z) if z else 1
        if idx<1 or idx>=n: continue
        if kind=="b": neutral=cum[idx]+nterm
        else: neutral=(cum[n]-cum[n-idx])+H2O+cterm
        # neutral-loss adjustment if annotated (e.g. y3-NH3)
        if "-NH3" in t: neutral-=17.026549
        elif "-H2O" in t: neutral-=18.010565
        theo=(neutral+zz*PROTON)/zz
        if theo>0:
            diffs.append((float(obs)-theo)/theo*1e6)
    if len(diffs)<4: return None
    diffs.sort(); return diffs[len(diffs)//2]

def base(name):
    return name.split(" (")[0].split("[")[0].strip()

def build_query(url, pep_max, lmin, lmax, min_consensus=0.0, frag_filter="hcd"):
    cols = {r[0] for r in duckdb.sql(f"DESCRIBE SELECT * FROM read_parquet('{url}')").fetchall()}
    charge = "precursor_charge" if "precursor_charge" in cols else "charge"
    pmz = "exp_mass_to_charge" if "exp_mass_to_charge" in cols else "observed_mz"
    if "cv_params" in cols:
        shape = duckdb.sql(f"SELECT typeof(cv_params) FROM read_parquet('{url}') LIMIT 1").fetchone()[0]
        method = (frag_filter or "hcd").lower()
        if method == "all":
            frag = "TRUE"
        elif "cv_name" in shape:
            frag = (
                f"list_contains(list_transform(cv_params, x -> lower(x.cv_value)), '{method}')"
            )
        else:
            frag = f"lower(cv_params.\"Fragmentation Method\") LIKE '%{method}%'"
        cvsel = "cv_params AS cvp"
    else:
        frag = "TRUE"; cvsel = "NULL AS cvp"
    consensus = (f"consensus_support >= {min_consensus}"
                 if (min_consensus and "consensus_support" in cols) else "TRUE")
    if "posterior_error_probability" in cols:
        pep = (
            f"(posterior_error_probability IS NULL "
            f"OR posterior_error_probability <= {pep_max})"
        )
    else:
        pep = "TRUE"
    return f"""
      SELECT sequence AS seq, CAST({charge} AS INTEGER) AS charge,
             CAST({pmz} AS DOUBLE) AS prec_mz, modifications AS mods,
             mz_array AS mz, intensity_array AS intensity, ion_type_array AS iontypes,
             charge_array AS charge_arr, {cvsel}
      FROM read_parquet('{url}')
      WHERE is_decoy = FALSE
        AND {pep}
        AND {charge} BETWEEN 2 AND 4
        AND length(sequence) BETWEEN {lmin} AND {lmax}
        AND mz_array IS NOT NULL
        AND {consensus}
        AND ({frag})
    """

def get_iso(cvp):
    if cvp is None: return 0
    try:
        if isinstance(cvp, dict):
            return int(cvp.get("isotope_error") or 0)
        if isinstance(cvp, list):
            for kv in cvp:
                if "isotope" in (kv.get("cv_name") or "").lower():
                    return int(float(kv.get("cv_value") or 0))
    except Exception:
        return 0
    return 0

def resolve_mods(seq, mods):
    res_pos, res_delta, nterm, cterm = [], [], 0.0, 0.0
    n = len(seq)
    for m in (mods or []):
        bn = base(m["name"]); d = UNIMOD.get(bn) or UNIMOD.get(bn.capitalize()) or _UNIMOD_LC.get(bn.lower())
        if d is None:
            return None
        for p in (m["positions"] or []):
            pos = int(p["position"])
            if bn in NTERM and pos <= 1: nterm += d
            elif pos == n: res_pos.append(pos); res_delta.append(d)
            elif pos > n: cterm += d
            else: res_pos.append(pos); res_delta.append(d)
    return res_pos, res_delta, nterm, cterm

def main():
    url, out = sys.argv[1], sys.argv[2]
    pep_max = float(sys.argv[3]) if len(sys.argv) > 3 else 0.001
    limit = sys.argv[4] if len(sys.argv) > 4 else None
    min_consensus = float(sys.argv[5]) if len(sys.argv) > 5 else 0.0
    min_matched = int(sys.argv[6]) if len(sys.argv) > 6 else 6
    min_explained = float(sys.argv[7]) if len(sys.argv) > 7 else 0.50
    recal = (len(sys.argv) > 8 and sys.argv[8] == "recal")
    frag_filter = (sys.argv[9] if len(sys.argv) > 9 else os.environ.get("FRAG_FILTER", "hcd"))
    inner = build_query(url, pep_max, 7, 40, min_consensus, frag_filter=frag_filter)
    q = f"SELECT * FROM ({inner})" + (f" USING SAMPLE {limit} ROWS" if limit else "")
    con = duckdb.connect(); con.execute("INSTALL httpfs; LOAD httpfs;")
    con.execute("SET memory_limit='6GB'; SET threads=2; SET preserve_insertion_order=false;")
    print(f"filters: PEP<={pep_max} frag={frag_filter} consensus>={min_consensus} matched_by>={min_matched} explained>={min_explained} isotope-corrected non-decoy")
    rel = con.execute(q); kept = skipped = nrecal = 0
    schema = pa.schema([("seq", pa.string()), ("charge", pa.int32()), ("prec_mz", pa.float64()),
        ("res_mod_pos", pa.list_(pa.int32())), ("res_mod_delta", pa.list_(pa.float64())),
        ("nterm_delta", pa.float64()), ("cterm_delta", pa.float64()),
        ("mz", pa.list_(pa.float32())), ("intensity", pa.list_(pa.float32()))])
    writer = pq.ParquetWriter(out, schema)
    while True:
        batch = rel.fetchmany(20000)
        if not batch: break
        c = {k: [] for k in schema.names}
        for seq, charge, prec_mz, mods, mz, inten, iontypes, charge_arr, cvp in batch:
            if mz is None or len(mz) < 10: skipped += 1; continue
            its = iontypes or []; nby = 0; expl = 0.0; tot = 0.0
            for t, ii in zip(its, (inten or [])):
                fi = float(ii); tot += fi
                if t and (t[0] in ("b", "y")): nby += 1; expl += fi
            if nby < min_matched: skipped += 1; continue
            if tot > 0 and (expl / tot) < min_explained: skipped += 1; continue
            r = resolve_mods(seq, mods)
            if r is None: skipped += 1; continue
            rp, rd, nt, ct = r
            iso = get_iso(cvp)
            pm = float(prec_mz) - (iso * 1.0033548 / int(charge)) if iso else float(prec_mz)
            mzc = mz
            if recal:
                off = spectrum_offset_ppm(seq, rp, rd, nt, ct, mz, iontypes, charge_arr)
                if off is not None:
                    f = 1.0 / (1.0 + off * 1e-6)
                    mzc = [float(x) * f for x in mz]
                    nrecal += 1
            c["seq"].append(seq); c["charge"].append(int(charge)); c["prec_mz"].append(pm)
            c["res_mod_pos"].append(rp); c["res_mod_delta"].append(rd)
            c["nterm_delta"].append(nt); c["cterm_delta"].append(ct)
            c["mz"].append([float(x) for x in mzc]); c["intensity"].append([float(x) for x in inten])
            kept += 1
        if c["seq"]: writer.write_table(pa.table(c, schema=schema))
    writer.close()
    print(f"kept={kept:,} skipped={skipped:,} recalibrated={nrecal:,} -> {out}")

if __name__ == "__main__": main()
