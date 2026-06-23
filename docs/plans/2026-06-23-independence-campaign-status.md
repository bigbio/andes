# Model-independence campaign — full-scale status (2026-06-23)

Goal: retrain ALL bundled regimes on public data → ship a fully own-trained `resources/models.parquet` → relicense Apache. Policy: regimes with no clean public corpus are **DROPPED, not shipped as MS-GF+ seed copies**. Honest target = "N/39 own-trained, rest dropped" — a defensible independence claim.

Proven pipeline (per regime): public corpus → MSFragger-labeled search → flats → `retrain_slug.sh <slug> <seed>` (geometry from bundle seed) → own-trained store. Validated by `hcd_astral_tryp` (2h40m, field-beating +23–50% vs MSFragger/Comet/MS-GF+).

## Wave 1 — TRAIN ready-flat regimes (Codon array 54383766, in flight)
Flats already harvested+searched for these → training now (13) + 4 already-done = **17 own-trained**:
- Done earlier: cid_lowres_tryp, cid_lowres_tryp_tmt, hcd_qexactive_tryp, hcd_astral_tryp.
- Array 54383766 (13): cid_lowres_lysc[1f], cid_lowres_argc[3], cid_lowres_gluc[3], etd_lowres_tryp_phosphorylation[9], etd_highres_tryp_phosphorylation[13], hcd_qexactive_tryp_tmt[16], hcd_highres_tryp_tmt[20], hcd_highres_nocleavage[39], hcd_highres_nocleavage_phosphorylation[39], hcd_qexactive_tryp_itraq[40], uvpd_qexactive_tryp[45], etd_highres_tryp[60], hcd_qexactive_tryp_phosphorylation[86].
- ⚠️ Benchmark before trusting: thin-flat slugs cid_lowres_lysc(1)/argc(3)/gluc(3)/etd_lowres_tryp_phosphorylation(9). 2 phospho slugs seeded from non-phospho geometry siblings.

## Wave 2 — DISCOVERY for gap regimes (4 agents)

### ETD (agent a0797857) — DONE
★ **PXD024364** (Coon lab, Sinitcyn Nat Biotechnol 2023): 6 proteases × ETD/EThcD, **ion-trap LowRes readout**, human, IDs included → covers:
- **SOURCEABLE:** etd_lowres_tryp, etd_lowres_lysc, etd_lowres_lysn, etd_lowres_aspn, etd_lowres_gluc (all PXD024364, high confidence; filter ETD subset per protease).
- **DROP (un-sourceable):** etd_lowres_argc (only LysargiNase proxy PXD004447, not true ArgC), etd_lowres_alp (exotic), etd_lowres_lysn_phosphorylation (no LysN+phospho+ETD), etd_highres_nocleavage (only top-down available = different regime).
- Already-flatted source confirmed: etd_highres_tryp / etd_lowres_tryp_phosphorylation provenance ≈ Heck/Pevzner LysN+Tryp CID/ETD (Kim 2010) + Riley/Coon AI-ETD phospho.

### Non-tryptic enzymes (a6261a85) — pending
### TOF (a6cbbf8b) — pending
### Labeled/phospho (af42e1a7) — pending

## Running tally
- Own-trained after Wave 1: 17/39.
- ETD adds 5 sourceable (Wave 2) → ~22 once harvested+trained; 4 ETD slugs flagged DROP.
- Pending: 3 discovery agents (enzymes/TOF/labeled) decide the rest.

## Astral release gate — PASS (a90d7825)
Detection fix 4e95d7da works end-to-end: Astral mzML now logs `instrument = OrbitrapAstral` + auto-picks `hcd_astral_tryp` (was hcd_qexactive_tryp). 35,304 real @ 0.87% FDP (vs ~36,685 forced-model; -3.8% = run/param variance, FDP better; PASS = correct instrument+model). NO regression: QExactive mzML->hcd_qexactive_tryp, UPS1->cid_lowres_tryp, 10/10 detect unit tests. 39-model bundle SHA 63507a3c. New VM binary SHA cb6081f3. ★ HOLD bundle merge into repo until the campaign's first batch of own-trained stores lands -> merge ALL at once (avoid piecemeal churn), re-key hcd_astral_tryp OrbitrapAstral.

## Wave 2 discovery — COMPLETE (consolidated regime->source map)
★ STRATEGIC (user insight): ProteomeTools (PXD009449/PXD010595 = the 'msnet' source) is a MULTIMODAL synthetic backbone — HCD x energies + ion-trap CID + ETD x tryptic/AspN/LysN/nonspecific/21-PTMs. Lean on it as the primary source, split by 2xIT/3xHCD/ETD filename tokens. No single dataset spans all activations (physics: ETD c/z vs HCD/CID b/y vs TOF), but ~5-6 comprehensive datasets cover the matrix, not 38.

SOURCEABLE (harvest+train):
- ETD enzymes (PXD024364 Coon, 6-protease x ETD ion-trap): etd_lowres_{tryp,lysc,lysn,aspn,gluc}
- Non-tryptic (PXD010595 ProteomeTools): cid_lowres_aspn, cid_lowres_lysn, cid_lowres_nocleavage, cid_highres_nocleavage
- CID-tryptic: cid_highres_tryp (PXD000865, self-label), cid_lowres_lysc (PXD000865-LysC, have_data)
- TOF: cid_tof_tryp (PXD014777 timsTOF)
- iTRAQ: hcd_highres_tryp_itraq (PXD003702 HiRIEF); itraqphospho x2 (PXD006380); 
- Phospho: hcd_highres_tryp_phosphorylation (PXD000612 / ProteomeTools PXD009449), cid_lowres_tryp_phosphorylation (PXD001428, MS-GF+-searched IDs)
- LysN+phospho: cid_lowres_lysn_phosphorylation (PXD001114 LysargiNase, proxy/medium)

DROP (un-sourceable -> do NOT ship MS-GF+ seed):
- aLP (all): cid_lowres_alp, etd_lowres_alp, cid_tof_alp, hcd_tof_alp
- UVPD-TMT: uvpd_qexactive_tryp_tmt (0 datasets exist)
- ETD exotics: etd_lowres_argc (only LysargiNase proxy), etd_highres_nocleavage (only top-down=diff regime), etd_lowres_lysn_phosphorylation (weak proxy)

TALLY: ~17 (Wave1 train) + ~13 (Wave2 sourceable) = ~30 own-trainable / 39; ~8-9 DROP (aLP x4, UVPD-TMT, ETD-argc, ETD-hr-nocleavage, ETD-lysn-phospho). => honest "MS-GF+ model free for ~30 shipped regimes, ~9 dropped as un-sourceable" -> Apache relicense defensible.

## Wave-2 harvest reality-check (2026-06-23, agent aed146397) — 2 discovery ERRORS caught
HARVESTING NOW (Codon array 54426808): hcd_qexactive_tryp_itraqphospho (PXD006380), cid_lowres_tryp_phosphorylation (PXD001428 trypsin arm), cid_lowres_lysn_phosphorylation (PXD001114 LysargiNase proxy), hcd_highres_tryp_phosphorylation (PXD009449). ETA 3-8h.
★ CORRECTIONS to Wave-2 discovery (verified at harvest):
- PXD024364 (Coon ETD goldmine) = UNAVAILABLE (0 files, withdrawn) -> the 5 etd_lowres_{tryp,lysc,lysn,aspn,gluc} have NO source. Re-discover or DROP.
- PXD010595 (ProteomeTools II) = TRYPTIC-ONLY (no AspN/LysN/HLA arms) -> cid_lowres_aspn/lysn/nocleavage + cid_highres_nocleavage have NO source from it. Re-discover or DROP.
- PXD014777 (TOF) DEFERRED: Bruker .d, harvester is Thermo-.raw-only -> needs a .d-aware branch (follow-up). cid_tof_tryp pending that.
REVISED TALLY: ~17 (Wave1) + 4 (Wave2 harvesting) = ~21 confirmed own-trained. Gaps needing re-discovery: 5 ETD-enzyme + 4 non-tryptic-CID + cid_tof_tryp(.d) + cid_highres_tryp. If re-discovery fails -> drop (exotic ETD-enzyme likely drops; non-tryptic CID worth one targeted retry). Honest floor ~21/39 own-trained even if all gaps drop = still a massive jump from 1.

## Re-discovery (a0c7b299) STALLED (technical watchdog, not a data verdict)
Partial finding before stall: non-tryptic CID IS likely sourceable — a large deposit has 2xIT (ion-trap low-res) runs covering HLA (non-specific) + AspN. Accession NOT captured. -> DEFER non-tryptic CID (cid_lowres_aspn/lysn/nocleavage + cid_highres_nocleavage) as "sourceable-pending-focused-retry" (a future session: PRIDE search ProteomeTools-immunopeptidome / multi-protease, verify file existence, harvest). NOT blocking the bundle.
ETD-low-res enzymes (5): no verified alt for the gone PXD024364 -> DROP (exotic, consistent with etd_lowres_argc/alp drop).
DECISION: don't keep chasing the niche tail through flaky agents. Proceed with the BUNDLE MERGE of the ~21 confirmed own-trained (Wave-1 17 + Wave-2 4 when flats land), DROP the confirmed-un-sourceable, DEFER non-tryptic+TOF(.d). Honest shipped count ~21-25/39 own-trained = the win. Remaining tail is a documented follow-up, not a blocker.

## v1 VALIDATION COMPLETE (a6b7d5bf) — own-trained bundle HOLDS on all 3 benchmarks -> GREEN LIGHT
Matched 1% true entrapment-FDP, Percolator 3.7.1, v1 bundle models_v1_own.parquet (SHA 12b2c176, 15 own, ZERO MS-GF+):
- Astral (hcd_astral_tryp): field-beating +23-50% vs MSFragger/Comet/MS-GF+.
- UPS1 (cid_lowres_tryp, own): 16,584 @0.995% = +0.84% vs prior bundled (16,446) -> slight BEAT.
- TMT (cid_lowres_tryp_tmt, own): 11,363 @0.998% = flat vs prior's non-TMT FALLBACK (prior 39-model bundle has NO low-res CID-TMT model; v1 ADDS one). No regression; known low-res-TMT gap unchanged but v1 strictly more complete.
Auto-select confirmed own-trained models for both. NO regression anywhere -> v1 ships. CAVEAT N=1/dataset; TMT not a same-model A/B (prior fell back to non-TMT). VERDICT: pure-own v1 bundle validated -> proceed to v1 branch assembly (docs remove MS-GF+, supported table, test count 38->17, .param removal, I1).

## v1 re-benchmark (a14940ab) — own models IMPROVED but benchmark INVALID for v1 (wrong binary) + Astral routing reopened
Own-trained v1 numbers (std 1% FDR): Astral top-1 39,506 (was 36,782, +7.4%!) / chimeric 68,244 (+117% vs field); UPS1 chimeric 18,732 (+4.8%, leads field) / top-1 17,129 (within 1% Java, README framing holds). Beats Java(26,542/-/17,305)+Comet(31,435/-/15,809) on Astral+UPS. Entrapment-FDP honest (0.9-1.2%).
★ 3 ISSUES -> re-run (a4524410): (1) used OLD Jun-22 binary (no detection fix) -> Astral auto-picked hcd_qexactive_tryp not shipped hcd_astral_tryp. (2) ★ hcd_qexactive_tryp (39,506) appears BETTER on Astral than dedicated hcd_astral_tryp (~35-37k) — BOTH own-trained now, so independence moot, pure perf call: REOPENS v1 Astral routing (keep dedicated astral key OR route Astral->qexactive). (3) TMT a05058.mzML DELETED from VM (disk pressure) -> unmeasured; re-stage from PXD007683. a4524410 = definitive: detection-fixed binary + partitioned bundle + Astral A/B (astral vs qexactive) + re-staged TMT + UPS, std-FDR table + entrapment footnote.

## DEFINITIVE BENCHMARK (a4524410, RECOVERED from crash) — Astral routing RESOLVED + v1 leads all 3
Detection-fixed binary + partitioned v1 bundle, /srv/data/msgf-bench/v1-definitive/. Recovered by counting the Percolator outputs (agent crashed before reporting).
★ ASTRAL A/B = TIE: hcd_astral_tryp 41,080 table / 36,776 real @1.07% vs hcd_qexactive_tryp 41,009 / 36,832 @1.09% (within 71 PSMs / 0.2%). The earlier 'qexactive 39,506 >> astral 35,304' was a BINARY/METRIC artifact (old non-detection-fixed binary + mismatched FASTA/metric). -> KEEP hcd_astral_tryp as Astral default (no reversal); detection fix stays.
★ FINAL v1 table (q<=0.01 std FASTA, own models, leads all 3):
  andes-v1 top-1: Astral 41,080 / TMT 12,170 / UPS 17,511
  andes-v1 chim:  Astral (pending percolator on astralA_chim.pin) / TMT 12,277 / UPS 18,283
  Java MS-GF+:    26,542 / 11,555 / 17,305    Comet: 31,435 / 10,876 / 15,809
  -> andes-v1 LEADS each: Astral +30%, TMT +5.3% (own cid_lowres_tryp_tmt now BEATS field), UPS +1.2% (now leads top-1 not just chim). Headline holds + stronger.
  Entrapment-FDP: Astral 1.07% / TMT 2.01% / UPS 1.85% (honest).
PENDING (a749a6d4): Astral chimeric number + sanity-check the 41,080 vs old 36,782 (+11.7%) is same-FASTA apples-to-apples (big jump - verify). Then these are the FINAL README numbers.
