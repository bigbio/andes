# Other N-glycopeptide tools — clean-room extraction for andes

Companion to a glyco search engine / a cross-spectrum glyco engine / O-Pair notes. Focus: what maps onto
andes' failure mode (identification = ranking + FDR, not generation). Papers
cited; code NOT read/copied.

## 1. GlycReSoft (Klein & Zaia; *Nat Commun* 2024, PMC11263600; DOI 10.1038/s41467-024-50338-5) — **Apache-2.0**, Python/Cython ([repo](https://github.com/mobiusklein/glycresoft))

The richest published, permissively-licensed reference. **Glycan-first**: builds a
**complement glycan-ion index** per precursor mass-shift, scans it, scores the
**top k = 70** (glycan-mass, peptide-mass) combinations/spectrum. Glycans =
**compositions** (from a glyco search engine DB or biosynthetic sim), one site/peptide. Y-ladder
is *semi-structured*: emit core fragments (± dHex/pentose side-chains) then every
combination of remaining non-labile monosaccharides.

Scoring is separated then mixed:
- Peptide `ScoreP = Σ log(intᵢ)·(1−ppmᵢ/tol)⁴ · coverageP^γ`, γ=1.
- Glycan `ScoreG = Σ log(intᵢ)·(1−ppmᵢ/tol)⁴ · covG^α · covG,core^β`, α=0.5, β=0.4.
- Combined `ScoreGP = w·ScoreP + (1−w)·ScoreG + SignIon() + MassAcc(...)`, **w=0.65**.

The **peak-relation fragmentation model**: a **Bayesian/multinomial log-linear**
model (UniNovo-style) yielding per-peak reliability φᵢ and intensity-correlation
Cor(θ); folded in as `ScoreModelP = Σ log₁₀(intᵢ)(1−ppmᵢ/tol)⁴(φᵢ+1) + φᵢ +
CorP(θ)·covP^γ`. They **tried gradient-boosted trees, saw overfitting, and chose
logistic regression** — a direct warning for andes' learned model. **FDR**: reverse
proteins keeping sequons at disrupted sites; **glycan decoy = give peptide+Y
fragments beyond Y1 a uniform random mass shift ∈ [1.0, 30.0] Da**; peptide and
glycan FDR estimated on **separate axes** (2D), joint 1%. Two post-hoc re-rankers:
**glycome network smoothing** (biosynthetic graph, propagate confidence to
unobserved glycoforms at a glycosite) and an **RT-revision model** (monosaccharide
RT contributions resolve dHex/NeuAc quasi-isobars fragmentation can't).

## 2. GlycoPAT (Liu et al.; *MCP* 2017, PMC5672007) — MATLAB, open-source ([GitHub](https://github.com/kaichengub/GlycoPAT))

**Peptide-first** (MS1 within 10 ppm → candidates). **Ensemble Score**
`ES = Σ sᵢwᵢ` over 4 terms: SEQUEST-style Xcorr, %-ion-match, top-10 peaks, and a
**Poisson match probability** `P(K₁,p)=C(N₁,K₁)·p^K₁·K₁!·e^(−N₁p)`, p=K/N. **Weights
switch by activation**: HCD w=(⅓,⅓,0,⅓); ETD w=(0.2,0,0.1,0.7) — activation-aware
weighting, relevant to andes' stepped-HCD sparsity. **FDR**: decoy = scramble
peptide **plus** perturb per-monosaccharide mass ±50 Da at fixed total glycan mass;
`FDR = #decoy(ES>c)/#target(ES>c)`. Notation: **SmallGlyPep** inline string.

## 3. GPSeeker (Xiao & Tian; *J. Proteome Res.* 2019, PMID 31117584; DOI 10.1021/acs.jproteome.9b00191) — structure-specific

**Structure** (not composition) engine. Key idea = **G score = count of
structure-diagnostic Y/Z fragment ions that uniquely discriminate isomers of the
same composition** (e.g. core- vs branch-fucose). Precursor/fragment matching via
**isotopic m/z + envelope fingerprinting (iMEF)** — strict per-isotope m/z *and*
abundance deviation. **Spectrum-level target-decoy FDR ≤ 1%**. Borrowable: a cheap
integer "how many peaks in this spectrum could *only* come from this candidate"
discriminator.

## 4. a glyco search engineQuant (Kong et al.; *Nat Commun* 2022, PMC9729625) — **CC-BY-4.0** ([GitHub](https://github.com/Power-Quant/a glyco search engineQuant))

A **quantifier/MBR layer over a glyco search engine/a commercial glyco engine/the reference engine outputs**, not an ID/FDR
engine — least relevant to andes' ranking gap. Uses **ResNet18** to embed a
glycopeptide chromatographic evidence vector and a **glycan-subnet** structure
(glycans grouped by backbone excluding sialics) for match-in-run. Its
target-decoy Gaussian-mixture **FQR** `= π₀P(X>t|T=0) / Σ πₖP(X>t|T=k)` is a clean
EM template for a Percolator-post-process mixture cutoff.

## 5. SugarPy (Schulze et al.; *Bioinformatics* 2020, PMID 33325487) — **GPL-3.0** ([GitHub](https://github.com/SugarPy/SugarPy))

Glycan-**database-independent**, discovery-driven: from an identified peptide it
builds and matches **peptide+Y glycan ladders** in MS1/ISCID, scoring by matched-Y
coverage. Confirms andes' Y-first premise but adds nothing new on ranking/FDR.
**GPL-3.0 → do NOT read its code** (copyleft, incompatible with andes/Apache).

## 6. Single most valuable borrow + what to AVOID

**BORROW (clean-room, Apache-safe): GlycReSoft's separated-axis scoring +
`Y-beyond-Y1 random-mass-shift [1,30] Da` glycan decoy.** This is exactly andes'
missing target/decoy-discriminating glyco feature and its 2D-FDR: keep peptide b/y
(Percolator target axis) and glycan-Y as an **independent** axis, generate glycan
decoys by mass-shifting the Y-ladder beyond Y1, and let **Percolator do peptide-axis
FDR while a thin post-process applies the glycan-axis cut** — directly matching
andes' hard constraint (Percolator-only, 2D as post-process) and fixing the G3
"unified pile crashed recovery 29→4%" result. Pair it with GlycoPAT's
**activation-conditioned weighting** for stepped-HCD sparsity and GPSeeker's
**diagnostic-Y-ion-count** as a cheap discriminator feature.

**AVOID:** (a) GlycReSoft's own empirical finding — **gradient-boosted trees
overfit** the fragmentation model; prefer a regularized/log-linear or
regime-matched learned model for SP-B. (b) A **single unified score axis** (GlycoPAT
ES / one Percolator Label) — collapses target vs decoy at one backbone mass (andes
already proved this). (c) **Structure-level enumeration** (GPSeeker) — combinatorial
blow-up andes doesn't need; stay composition-level. (d) Any a glyco search engineQuant/SugarPy
code (CC-BY narrative fine to cite; GPL-3.0 code is off-limits).
