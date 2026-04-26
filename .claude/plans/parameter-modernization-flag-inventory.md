# MS-GF+ flag inventory (Phase 1 input)

Snapshot of every flag registered by `ParamManager.addMSGFPlusParams()`
plus the parsing semantics each one currently relies on. This is the
foundation document for the Phase 1 picocli rewrite described in
`parameter-modernization.md`. Total: 34 flags (27 visible + 7 hidden).
Required: `-s`, `-d`.

## Visible flags

| Short | Canonical name | Type | Default | Bounds | Notes |
|---|---|---|---|---|---|
| `-conf` | `ConfigurationFile` | file | — | exists | Config file; CLI overrides config |
| `-s` | `SpectrumFile` | file/dir | — | exists | **Required.** mzML/mzXML/mgf/ms2/pkl/_dta.txt or directory |
| `-d` | `DatabaseFile` | file | — | exists | **Required.** *.fasta / *.fa / *.faa |
| `-decoy` | `DecoyPrefix` | string | `DECOY_` | — | Decoy protein prefix |
| `-o` | `OutputFile` | file | `<spec>.pin` | — | *.pin (default) or *.tsv |
| `-t` | `PrecursorMassTolerance` | tolerance | `20ppm` | ≥0 | Symmetric (`20ppm`) or asymmetric (`0.5Da,2.5Da`); units must match |
| `-ti` | `IsotopeErrorRange` | int range | `0,1` | ≥0, max-incl | Isotope-error window, both ends inclusive |
| `-m` | `FragmentationMethodID` | dyn-enum | `ASWRITTEN` | — | 0=as-written, 1=CID, 2=ETD, 3=HCD |
| `-inst` | `InstrumentID` | dyn-enum | `LOW_RES_LTQ` | registry | `InstrumentType` registry-driven |
| `-e` | `EnzymeID` | dyn-enum | `TRYPSIN` | registry | `Enzyme` registry-driven |
| `-protocol` | `ProtocolID` | dyn-enum | `AUTOMATIC` | registry | `Protocol` registry-driven |
| `-ntt` | `NTT` | enum | `2` | 0..2 | Number of tolerable termini |
| `-mod` | `ModificationFile` | file | built-in (C+57) | exists | Mod file; config-file path also accepts `StaticMod=`/`DynamicMod=`/`CustomAA=` |
| `-minLength` | `MinPepLength` | int | `6` | ≥1 | |
| `-maxLength` | `MaxPepLength` | int | `40` | ≥1 | |
| `-minCharge` | `MinCharge` | int | `2` | ≥1 | |
| `-maxCharge` | `MaxCharge` | int | `3` | ≥1 | |
| `-n` | `NumMatchesPerSpec` | int | `1` | ≥1 | |
| `-thread` | `NumThreads` | int | `Runtime.availableProcessors()` | ≥1 | |
| `-tasks` | `NumTasks` | int | `0` (auto) | ≥-10 | 0=auto, >0=fixed, <0=N×threads |
| `-minSpectraPerThread` | `MinSpectraPerThread` | int | `250` | ≥1 | |
| `-verbose` | `Verbose` | enum | `0` | 0..1 | 0=total, 1=per-thread |
| `-tda` | `TDA` | enum | `0` | 0..1 | 0=no decoy, 1=concat decoy search |
| `-addFeatures` | `AddFeatures` | enum | `0` | 0..1 | Percolator extra features |
| `-outputFormat` | `OutputFormat` | enum | `pin` | pin/tsv | mzIdentML removed |
| `-precursorCal` | `PrecursorCal` | string | `auto` | auto/on/off | Case-insensitive |
| `-ccm` | `ChargeCarrierMass` | double | `1.00727649` | >0.1 | Proton mass default |
| `-maxMissedCleavages` | `MaxMissedCleavages` | int | `-1` | ≥-1 | -1 = unlimited |
| `-numMods` | `NumMods` | int | `3` | ≥0 | Max dynamic mods per peptide |
| `-allowDenseCentroidedPeaks` | `AllowDenseCentroidedPeaks` | enum | `0` | 0..1 | |
| `-msLevel` | `MSLevel` | int range | `2,2` | ≥1, max-incl | `min,max` or single |
| `-u` | `PrecursorMassToleranceUnits` | enum | `2` | 0..2 | **Hidden** — legacy; 0=Da, 1=ppm, 2=as-written |

## Hidden flags

| Short | Canonical name | Type | Default | Notes |
|---|---|---|---|---|
| `-dd` | `DBIndexDir` | dir | — | Database index dir |
| `-index` | `SpecIndex` | int range | `1,INT_MAX-1` | Spectrum index range, both inclusive |
| `-edgeScore` | `EdgeScore` | enum | `0` | 0=use, 1=skip |
| `-minNumPeaks` | `MinNumPeaks` | int | `Constants.MIN_NUM_PEAKS_PER_SPECTRUM` | |
| `-iso` | `NumIsoforms` | int | `Constants.NUM_VARIANTS_PER_PEPTIDE` | |
| `-ignoreMetCleavage` | `IgnoreMetCleavage` | enum | `0` | 0=consider, 1=ignore |
| `-minDeNovoScore` | `MinDeNovoScore` | int | `Constants.MIN_DE_NOVO_SCORE` | |

## Sharp edges the picocli rewrite must preserve

1. **Asymmetric tolerance.** `-t 0.5Da,2.5Da` → left tolerance (observed < theoretical) ≠ right tolerance. Both sides must use the same unit. Numeric-only value (e.g. `20`) defaults to Da. Trailing unit suffix is case-insensitive (`Da`/`ppm`/`Th`).
2. **Range inclusivity is per-flag.** `IntRangeParameter` defaults to `min` inclusive / `max` exclusive, but `-ti`, `-index`, `-msLevel` flip max to inclusive via `.setMaxInclusive()`.
3. **Dynamic enums.** `-inst`, `-e`, `-protocol`, `-m` are registry-driven (`InstrumentType`, `Enzyme`, `Protocol`, `ActivationMethod`). Numeric indices depend on registry load order; help text is generated at startup. Picocli converters must read from the same registries, not hardcode indices.
4. **`OutputFormat` legacy mapping is gone.** Old `0=mzIdentML`, `2=both` are no longer accepted; only `pin` (0) and `tsv` (1) remain. Numeric indices are deprecated but still parse internally.
5. **`-precursorCal` is a string, not an enum class.** Values: `auto` / `on` / `off` (case-insensitive, `.trim()`-ed). `auto` means "run pre-pass, apply only if ≥200 confident PSMs collected".
6. **Trailing `!` on numbers.** `IntParameter` and `DoubleParameter` strip trailing `!` (legacy DMS config-file integration). Decide if Phase 1 keeps this quirk.
7. **`-tasks` semantics.** `0` = auto, `>0` = fixed, `<0` = `N × threads`. Range allows down to `-10`.
8. **Config-file-only entries.** `StaticMod=`, `DynamicMod=`, `CustomAA=` are not CLI flags. They're parsed from `-mod` file and `-conf` config file only. Repeated entries are *expected* (each line is a separate mod). Config parser preserves order.
9. **Config-file aliases (canonical-name normalization in `ParamNameEnum.getParamNameFromLine()`).** Auto-renames at least 13 deprecated keys:
   - `IsotopeError` → `IsotopeErrorRange`
   - `TargetDecoyAnalysis` → `TDA`
   - `FragmentationMethod` → `FragmentationMethodID`
   - `Instrument` → `InstrumentID`
   - `Enzyme` → `EnzymeID`
   - `Protocol` → `ProtocolID`
   - `NumTolerableTermini` → `NTT`
   - `MinNumPeaks` → `MinNumPeaksPerSpectrum`
   - `MaxNumMods` / `MaxNumModsPerPeptide` → `NumMods`
   - `minLength` / `MinPeptideLength` → `MinPepLength`
   - `maxLength` / `MaxPeptideLength` → `MaxPepLength`
   - `PMTolerance` / `ParentMassTolerance` → `PrecursorMassTolerance`
10. **File-format validation chain.** Order: directory-vs-file → format-suffix match → existence → no-reuse. Suffix matching is case-insensitive for `.pin`/`.tsv`/`.fasta`. Spec parameter auto-allows directories.
11. **Defaults that depend on runtime.** `-thread` defaults to `Runtime.getRuntime().availableProcessors()` (includes hyperthreading; per CLAUDE.md, physical cores often give better wall-time).
12. **Help-text drift.** Existing tests likely compare exact `--help` output. picocli's formatter is different. Decide: snapshot-update vs. custom renderer that mimics current format.

## Out-of-scope reminders for Phase 1

- `MSGFDB`, `MSGF`, `MSGFLib` entry points share `ParamManager`. Phase 1 only modernizes `MSGFPlus`; the other three keep using `ParamManager.parseParams()` until Phase 4.
- Config-file parsing is Phase 2. Phase 1 covers CLI only.
- The `Parameter` / `IntParameter` / `IntRangeParameter` / `ToleranceParameter` / etc. hierarchy is **not** removed in Phase 1. Removal is Phase 3.
- `ParamManager` itself stays. Phase 1 adds an adapter that produces a populated `ParamManager` from the typed `MSGFPlusOptions`, so `SearchParams.parse(ParamManager)` is unchanged.
