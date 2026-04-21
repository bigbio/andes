# Isobaric-labelling searches: TMT, TMTpro, iTRAQ

[MS-GF+ Documentation home](readme.md) · [ChangeLog](changelog.md)

MS-GF+ supports isobaric-labelled samples (TMT-6/10/11/16, TMTpro-18, iTRAQ-4/8) natively — **no code change is needed for new plex levels**. The label is expressed as a fixed (static) modification on the peptide N-terminus and on lysine, plus the matching `-protocol` flag so the scoring model is tuned for enriched spectra.

This page collects the mod-file recipes users most often ask for. The numbers come from Unimod.

---

## Command line

```bash
java -Xmx8G -jar MSGFPlus.jar \
  -s spectra.mzML -d db.fasta \
  -mod Mods_TMTpro.txt \
  -protocol 4 \
  -tda 1 -t 10ppm -ti -1,2 -ntt 2 -inst 1 -m 3
```

The `-protocol` values relevant to labelled samples:

| `-protocol` | Use for |
|-------------|---------|
| 2           | iTRAQ-labelled samples |
| 3           | iTRAQ + phosphopeptide enrichment |
| 4           | TMT / TMTpro labelled samples |
| 0           | Automatic (inferred from modification names in the mod file; see below) |

With `-protocol 0`, MS-GF+ detects the label from case-insensitive name prefixes (`itraq`, `tmt`, `phospho`) in the mod file. Naming the fixed mod `TMT6plex` or `iTRAQ8plex` is enough; you don't need to also pass `-protocol 4`.

---

## Modification recipes

Drop the appropriate block into your Mods.txt (see [examples/Mods.txt](examples/Mods.txt) for the full file format). Mass values are monoisotopic.

### TMT-6plex / TMT-10plex / TMT-11plex

Same reagent mass (229.162932) — the channels differ in isotope distribution but not monoisotopic mass.

```text
229.162932, K,    fix, any,    TMT6plex   # Fixed TMT on lysine
229.162932, *,    fix, N-term, TMT6plex   # Fixed TMT on peptide N-term
```

### TMT-16plex / TMTpro-16 and TMTpro-18

Both plex levels use the same reagent mass (304.207146).

```text
304.207146, K,    fix, any,    TMTpro     # Fixed TMTpro on lysine
304.207146, *,    fix, N-term, TMTpro     # Fixed TMTpro on peptide N-term
```

Use `-protocol 4` (same as TMT-6/10/11).

### iTRAQ-4plex

```text
144.102063, K,    fix, any,    iTRAQ4plex
144.102063, *,    fix, N-term, iTRAQ4plex
```

### iTRAQ-8plex

```text
304.205360, K,    fix, any,    iTRAQ8plex
304.205360, *,    fix, N-term, iTRAQ8plex
```

⚠ **TMTpro-16/18 and iTRAQ-8 have the same nominal mass (≈304.2) but different monoisotopic masses** (304.207146 vs 304.205360). Use the correct Unimod value for your reagent.

---

## Typical full Mods.txt for TMTpro + Carbamidomethyl + Ox-M

```text
NumMods=3

# Isobaric label (fixed)
304.207146, K,    fix, any,    TMTpro
304.207146, *,    fix, N-term, TMTpro

# Alkylation (fixed)
C2H3N1O1,   C,    fix, any,    Carbamidomethyl

# Common variable mods
O1,         M,    opt, any,    Oxidation
42.010565,  *,    opt, N-term, Acetyl
```

---

## Tips

- For phospho-enriched, TMT-labelled samples use `-protocol 3` or include both `TMT*` and `Phospho` mods with `-protocol 0`.
- If you see unexpectedly few PSMs with a labelled dataset, the single most common cause is forgetting to add the N-terminal fixed mod (only lysine is set). MS-GF+ needs both.
- For correct mzIdentML output, keep the Unimod **PSI-MS name** in the last column (e.g. `TMTpro`, `iTRAQ8plex`), not a free-text description.

Related upstream issue asking how to configure TMT-16: [#82](https://github.com/MSGFPlus/msgfplus/issues/82).
