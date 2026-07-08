use crate::glycan_mass::{CORE_Y_STEPS, PROTON};

/// Water molecule mass (monoisotopic).
pub(crate) const H2O: f64 = 18.010565;

/// Constant for b/y complement-pair sum: for a backbone of neutral mass `bb`,
/// a singly-charged b-ion and its complementary singly-charged y-ion satisfy:
///   b_mz + y_mz = bb + H2O + 2*PROTON
/// so we check |p1_mz + p2_mz - (bb + BY_COMPLEMENT_OFFSET)| <= tol.
const BY_COMPLEMENT_OFFSET: f64 = H2O + 2.0 * PROTON; // 18.010565 + 2*1.0072765 = 20.025118

#[derive(Debug, Clone)]
pub struct BackboneCandidate {
    pub backbone_mass: f64,
    pub core_y_hits: u8,
    pub votes: u32,
    /// Sum of sqrt-compressed, base-peak-normalised intensities of all peaks
    /// that voted for this backbone cluster. Primary ranking after core_y_hits.
    pub intensity_score: f64,
    /// Complement-pair (b/y) confirmation score: sum of min(i1, i2) for peak
    /// pairs (p1_mz, p2_mz) where |p1_mz + p2_mz - (bb + BY_COMPLEMENT_OFFSET)|
    /// ≤ max(|bb + BY_COMPLEMENT_OFFSET| * 20 ppm, 0.02). Spurious clusters have
    /// few/no complement pairs; real backbones with b/y ladders score higher.
    pub complement_score: f64,
}

/// Count complement (b/y) pairs for a candidate backbone neutral mass `bb`.
///
/// Returns the sum of min(sqrt_norm_i1, sqrt_norm_i2) over unique singly-charged
/// (b, y) peak pairs where:
///   - b_mz ∈ [50, bb/2]  (singly-charged b-ion, below midpoint — avoids glycan-heavy region)
///   - y_mz = target - b_mz  (complementary y-ion)
///   - |b_mz + y_mz - target| ≤ tol
///
/// Restricting b_mz < bb/2 ensures each pair is counted only once and avoids the
/// large-mz glycan-ion region (glycan Y-ions dominate the upper half of the spectrum
/// and create incidental pair noise for every candidate backbone).
///
/// `peaks` must be sorted by m/z ascending.
/// `norm_sqrt` is the sqrt-compressed base-peak-normalised intensity array
/// (parallel to `peaks`).
fn complement_score(peaks: &[(f64, f32)], norm_sqrt: &[f64], bb: f64) -> f64 {
    let target = bb + BY_COMPLEMENT_OFFSET;
    // tolerance: ±20 ppm of the target sum, floor 0.02 Da
    let tol = (target * 20e-6).max(0.02);

    let n = peaks.len();
    if n < 2 {
        return 0.0;
    }

    // b-ion window: singly-charged peptide backbone b-ions are in [50, bb/2].
    // This cuts the glycan-heavy region (y-ion series > bb/2) which creates
    // noise pairs for every candidate regardless of whether it is the true backbone.
    let b_lo = 50.0_f64;
    let b_hi = bb / 2.0; // midpoint — each pair (b, y) has b < y when b < bb/2

    let mut score = 0.0_f64;

    // For each candidate b-ion (mz ∈ [b_lo, b_hi]), look up its complement y at
    // y_mz = target - b_mz using a binary search into the sorted peak list.
    for (i, &(b_mz, _)) in peaks.iter().enumerate() {
        if b_mz < b_lo {
            continue;
        }
        if b_mz > b_hi {
            break; // peaks are sorted; nothing further qualifies as b-ion
        }
        let y_mz = target - b_mz;
        // y_mz must be > b_mz (ensured by b_mz < bb/2 since target ≈ bb+20), > 0
        if y_mz <= b_mz || y_mz <= 0.0 {
            continue;
        }
        // Binary search for y_mz ± tol in sorted peaks (skip self)
        let y_lo = y_mz - tol;
        let y_hi = y_mz + tol;
        let start = peaks.partition_point(|&(m, _)| m < y_lo);
        let end = peaks.partition_point(|&(m, _)| m <= y_hi);
        for j in start..end {
            if j == i {
                continue; // self-pair guard (shouldn't happen given b<y constraint)
            }
            let w = norm_sqrt[i].min(norm_sqrt[j]);
            score += w;
        }
    }
    score
}

/// Y-ion-first backbone solver at the default core-Y quorum (≥2 distinct core-Y
/// rungs). See [`solve_backbone_min`] for the quorum-parameterised variant used
/// by the confidence-cascade rescue.
pub fn solve_backbone(
    peaks: &[(f64, f32)],
    precursor_neutral: f64,
    precursor_z: u8,
    tol_ppm: f64,
    top_k: usize,
) -> Vec<BackboneCandidate> {
    solve_backbone_min(peaks, precursor_neutral, precursor_z, tol_ppm, top_k, 2)
}

/// As [`solve_backbone`], but with a configurable minimum core-Y rung quorum
/// (`min_core_y`). The primary path uses `2` (conservative, few candidates); the
/// rescue path relaxes to `1` (Y0/Y1-only, weak-ladder spectra) — still
/// evidence-driven and bounded to `top_k`, unlike a blind precursor−glycan
/// enumeration. `min_core_y` is clamped to at least 1 (a candidate always has
/// its own Y0/anchor rung, so `0` would admit every mass bin).
pub fn solve_backbone_min(
    peaks: &[(f64, f32)],
    precursor_neutral: f64,
    precursor_z: u8,
    tol_ppm: f64,
    top_k: usize,
    min_core_y: u8,
) -> Vec<BackboneCandidate> {
    use std::collections::HashMap;

    let min_core_y = min_core_y.max(1);

    if peaks.is_empty() {
        return Vec::new();
    }

    // --- intensity normalisation: sqrt-compress then normalise to base peak ---
    // sqrt compression prevents a single huge peak from dominating (e.g. oxonium
    // at m/z 204 that happens to vote for a spurious low-mass backbone), while
    // still weighting true Y-ions (which tend to be higher-intensity) over noise.
    let base_peak_intensity = peaks
        .iter()
        .map(|&(_, i)| i)
        .fold(0.0_f32, f32::max)
        .max(1.0_f32); // avoid divide-by-zero on all-zero spectra

    // Precompute normalised sqrt intensities once.
    let norm_sqrt: Vec<f64> = peaks
        .iter()
        .map(|&(_, i)| (i / base_peak_intensity).max(0.0).sqrt() as f64)
        .collect();

    // Minimum plausible backbone mass: real N-glycopeptides are >=500 Da
    // (smallest tryptic glycopeptide with NXS/T sequon). This cuts the very common
    // spurious cluster at ~203 Da (HexNAc oxonium m/z 204 at z=1 votes as Y0).
    const MIN_BB: f64 = 500.0;
    // Minimum implied glycan mass for N-glycopeptides: the core requires at least
    // 2×HexNAc (~406 Da). This cuts spurious high-mass backbone candidates that
    // leave no room for a real glycan.
    const MIN_GLYCAN: f64 = 406.0; // 2×HexNAc = minimum N-glycan core

    let prec_z_max = precursor_z.max(1);
    let neutral = |mz: f64, z: f64| (mz - PROTON) * z;

    // Per bin key, accumulate:
    //   .0  raw vote count (u32)             — for centroid/dedup
    //   .1  per-rung best-intensity seen      — [f64; 6], one slot per Y-ion rung
    //   .2  vote-weighted mass sum (f64)      — for centroid
    //   .3  sum of normalised sqrt intensities — intensity ranking score
    //
    // KEY CHANGE: each rung slot records the MAXIMUM intensity that voted for it
    // (de-duplicated across charges). This means a single peak seen at z=2 and z=3
    // does NOT double-count; only its best weight contributes to the intensity score.
    // This collapses the spurious inflation of low-mass clusters.
    //
    // Per-key per-rung accumulation: we first collect per (key, rung) the best
    // weight from any charge, then aggregate into the bucket.
    //
    // Data layout: outer HashMap key = backbone bin key; inner = rung index (0–5),
    // value = (best_w_for_rung, mass_accum, vote_count). We compute the final
    // intensity_score as the sum of per-rung best weights.

    // Two-pass: collect best-weight per (backbone_bin, rung) then aggregate.
    // We use a flat HashMap keyed by (backbone_bin, rung_idx) for efficiency.
    let mut rung_best: HashMap<(i64, u8), (f64 /*best_w*/, f64 /*mass for centroid*/, u32 /*vote count*/)> =
        HashMap::new();

    // Per-rung weights: Y0 (bare backbone) and Y1 (first HexNAc) are the most
    // diagnostic evidence for a specific backbone mass; they are unique to the
    // backbone and not part of a repeated glycan motif. Y3-Y5 (Hex additions)
    // are less specific because hexose residues are common in many glycan
    // antennae structures and can arise coincidentally.
    //   ri: 0=Y0, 1=Y1, 2=Y2, 3=Y3, 4=Y4, 5=Y5
    const RUNG_WEIGHT: [f64; 6] = [2.0, 2.0, 1.5, 1.0, 1.0, 1.0];

    for (peak_idx, &(pmz, _pi)) in peaks.iter().enumerate() {
        let w = norm_sqrt[peak_idx];
        for z in 1..=prec_z_max {
            let pn = neutral(pmz, z as f64);
            if pn <= 0.0 || pn > precursor_neutral {
                continue;
            }
            for (ri, step) in std::iter::once(0.0)
                .chain(CORE_Y_STEPS.iter().copied())
                .enumerate()
            {
                let bb = pn - step;
                // Plausibility gates:
                //  1. backbone must be within [MIN_BB, precursor_neutral)
                //  2. implied glycan (precursor − backbone) must be ≥ MIN_GLYCAN
                //     to avoid spurious candidates where the backbone consumes nearly
                //     the entire precursor mass, leaving no room for a real glycan.
                // `bb` is the Y0-derived peptide NEUTRAL mass (water INCLUDED), but
                // `precursor_neutral` is the RESIDUE convention (water excluded). The
                // implied-glycan / upper-bound gates MUST compare in the residue
                // convention (`bb - H2O`), else they undercount the implied glycan by
                // H2O and silently drop minimal-core (2×HexNAc, 406.16 Da)
                // glycopeptides (the off-by-H2O gate audit finding).
                let bb_res = bb - H2O;
                if bb < MIN_BB || bb_res >= precursor_neutral {
                    continue;
                }
                if precursor_neutral - bb_res < MIN_GLYCAN {
                    continue;
                }
                let key = (bb * 100.0).round() as i64;
                let tol_key = ((bb * tol_ppm / 1e6).max(0.01) * 100.0).round() as i64;
                let ri_u8 = ri.min(5) as u8;
                // Rung-weighted intensity: Y0 and Y1 are more diagnostic than Y3-Y5.
                let rw = RUNG_WEIGHT[ri_u8 as usize];
                let weighted_w = w * rw;
                for k in (key - tol_key)..=(key + tol_key) {
                    let e = rung_best
                        .entry((k, ri_u8))
                        .or_insert((0.0, 0.0, 0));
                    // keep the best rung-weighted intensity seen (charge-dedup)
                    if weighted_w > e.0 {
                        e.0 = weighted_w;
                    }
                    e.1 += bb; // for vote-weighted centroid
                    e.2 += 1;
                }
            }
        }
    }

    // Aggregate per (bin, rung) into per-bin stats.
    // Per-bin stats: (raw_votes: u32, rung_hit_mask: [bool;6], mass_sum: f64, intensity_score: f64)
    let mut bins: HashMap<i64, (u32, [bool; 6], f64, f64)> = HashMap::new();

    for ((k, ri), (best_w, mass_sum, vote_count)) in rung_best {
        let e = bins.entry(k).or_insert((0, [false; 6], 0.0, 0.0));
        e.0 += vote_count;
        if (ri as usize) < 6 {
            e.1[ri as usize] = true;
        }
        e.2 += mass_sum;
        // intensity score = sum of per-rung best weights; each rung contributes
        // at most once (the best weight seen across all charges + all peaks
        // that could have generated this rung).
        e.3 += best_w;
    }

    // Sort peaks by m/z for the two-pointer complement-score sweep.
    // Note: peaks passed in may not be sorted; we need a sorted copy.
    // Avoid cloning the whole peaks vec if already sorted (most mzML parsers sort).
    let sorted_peaks: Vec<(f64, f32)>;
    let sorted_norm: Vec<f64>;
    let (sp, sn): (&[(f64, f32)], &[f64]) = {
        let already_sorted = peaks
            .windows(2)
            .all(|w| w[0].0 <= w[1].0);
        if already_sorted {
            (peaks, &norm_sqrt)
        } else {
            let mut idx: Vec<usize> = (0..peaks.len()).collect();
            idx.sort_by(|&a, &b| peaks[a].0.partial_cmp(&peaks[b].0).unwrap_or(std::cmp::Ordering::Equal));
            sorted_peaks = idx.iter().map(|&i| peaks[i]).collect();
            sorted_norm = idx.iter().map(|&i| norm_sqrt[i]).collect();
            (&sorted_peaks, &sorted_norm)
        }
    };

    // Build candidates from bins. Core-Y quorum: the primary path requires ≥2
    // distinct core-Y rungs; the relaxed rescue (min_core_y=1) additionally
    // requires b/y COMPLEMENT support for a single-rung candidate (Codex #3).
    //
    // PERF (profiled to 45% of the glyco phase): compute the O(peaks) complement
    // score LAZILY — only for bins that already clear the cheap core-Y count.
    // The old code computed it for EVERY bin (thousands of spurious single-peak
    // bins) via .map() before the .filter(), which dominated runtime. This
    // filter_map is result-IDENTICAL (same survivors, same scores).
    let mut cands: Vec<BackboneCandidate> = bins
        .into_iter()
        .filter_map(|(_k, (v, hits, mass_sum, int_score))| {
            let core_y_hits = hits.iter().filter(|&&h| h).count() as u8;
            // `min_core_y` is already clamped to ≥1. A bin below it can never pass
            // the quorum, so reject BEFORE the expensive complement_score.
            if core_y_hits < min_core_y {
                return None;
            }
            let bb = mass_sum / v as f64;
            let cscore = complement_score(sp, sn, bb);
            if core_y_hits >= 2 || cscore > 0.0 {
                Some(BackboneCandidate {
                    backbone_mass: bb,
                    core_y_hits,
                    votes: v,
                    intensity_score: int_score,
                    complement_score: cscore,
                })
            } else {
                None
            }
        })
        .collect();

    // Sort: PRIMARY = core_y_hits (more distinct rungs = stronger evidence),
    // SECONDARY = combined score: intensity_score + complement_score * COMPLEMENT_WEIGHT
    //   (complement pairs from real b/y ladders break ties; spurious clusters have none),
    // TERTIARY = backbone_mass ascending (deterministic tiebreak).
    //
    // COMPLEMENT_WEIGHT tuned so complement evidence has meaningful influence but does
    // not override a large core_y_hits gap. A weight of 0.3 means a candidate needs
    // ~3 extra complement pairs per unit of intensity score deficit to win the tie.
    const COMPLEMENT_WEIGHT: f64 = 0.3;
    cands.sort_by(|a, b| {
        b.core_y_hits
            .cmp(&a.core_y_hits)
            .then({
                let sa = a.intensity_score + a.complement_score * COMPLEMENT_WEIGHT;
                let sb = b.intensity_score + b.complement_score * COMPLEMENT_WEIGHT;
                sb.partial_cmp(&sa)
                    .unwrap_or(std::cmp::Ordering::Equal)
            })
            .then(
                a.backbone_mass
                    .partial_cmp(&b.backbone_mass)
                    .unwrap_or(std::cmp::Ordering::Equal),
            )
    });

    // Merge near-mass clusters into a single vote-weighted centroid. We keep the
    // higher-ranked representative (kept) and fold the removed neighbour (cur) into it.
    // dedup_by calls the closure as (cur, kept): `cur` is the element under
    // inspection and `kept` is the surviving representative. Returning true
    // removes `cur`; we first fold its votes/mass into `kept` so the survivor
    // reports the cluster CENTROID rather than its lowest-mass edge.
    cands.dedup_by(|cur, kept| {
        if (kept.backbone_mass - cur.backbone_mass).abs() < 0.05 {
            let total = kept.votes as f64 + cur.votes as f64;
            kept.backbone_mass = (kept.backbone_mass * kept.votes as f64
                + cur.backbone_mass * cur.votes as f64)
                / total;
            kept.votes += cur.votes;
            kept.intensity_score += cur.intensity_score;
            kept.complement_score = kept.complement_score.max(cur.complement_score);
            kept.core_y_hits = kept.core_y_hits.max(cur.core_y_hits);
            true
        } else {
            false
        }
    });
    cands.truncate(top_k);
    cands
}

/// Count core-Y ladder support ions for a given backbone mass in a spectrum.
///
/// For a candidate backbone mass `bb`, the core-Y ion ladder produces peaks at:
///   Y0 = bb + PROTON
///   Y1 = bb + PROTON + CORE_Y_STEPS[0]
///   ...
///   Y5 = bb + PROTON + CORE_Y_STEPS[4]
///
/// Returns the number of those 6 ions (Y0 through Y5) that have a matching peak
/// in `peaks` within the match tolerance: max(mz * tol_ppm / 1e6, 0.01) Da.
///
/// `peaks` need not be sorted. This function performs binary-search if sorted,
/// and falls back to a linear scan (but for typical spectrum sizes, linear is fine).
/// Intensity-weighted core-Y ladder match. Like [`count_core_y_hits`] but sums
/// the base-peak-normalised intensity of each matched core-Y ion (Y0..Y5)
/// instead of counting them. A real glycopeptide backbone matches the ladder
/// with high intensity; a wrong backbone matches little or nothing, so this
/// discriminates true from false better than the bare count.
///
/// `bb` is the NEUTRAL peptide backbone mass (Y0 = bb + PROTON), matching
/// `solve_backbone`'s output convention (NOT the residue mass).
/// Best raw (un-normalised) intensity of a fragment of NEUTRAL mass `neutral`,
/// searched at charges 1..=`max_charge.clamp(1,3)`. Multiply-charged Y-ions are
/// common in stepped-HCD glyco spectra — the intact glycopeptide is 2–6+, so its
/// glycan-retaining Y-ions are frequently 2+/3+; matching only +1 leaves that
/// signal unmatched. `sorted` = peaks are m/z-ascending (enables binary search).
/// Returns 0.0 for no match at any charge.
#[inline]
fn best_frag_intensity(
    peaks: &[(f64, f32)],
    sorted: bool,
    neutral: f64,
    tol_ppm: f64,
    max_charge: u8,
) -> f32 {
    let zmax = max_charge.clamp(1, 3);
    let mut best = 0.0f32;
    for z in 1..=zmax {
        let ion_mz = (neutral + z as f64 * PROTON) / z as f64;
        let tol = (ion_mz * tol_ppm / 1e6).max(0.01);
        let (lo, hi) = (ion_mz - tol, ion_mz + tol);
        let found = if sorted {
            let start = peaks.partition_point(|&(m, _)| m < lo);
            peaks[start..]
                .iter()
                .take_while(|&&(m, _)| m <= hi)
                .map(|&(_, i)| i)
                .fold(0.0f32, f32::max)
        } else {
            peaks
                .iter()
                .filter(|&&(m, _)| m >= lo && m <= hi)
                .map(|&(_, i)| i)
                .fold(0.0f32, f32::max)
        };
        best = best.max(found);
    }
    best
}

pub fn core_y_intensity(peaks: &[(f64, f32)], bb: f64, tol_ppm: f64, max_charge: u8) -> f64 {
    let base = peaks
        .iter()
        .map(|&(_, i)| i)
        .fold(0.0f32, f32::max)
        .max(1e-9) as f64;
    // Core-Y ion NEUTRAL masses (Y0..Y5): Y0 = bb, rungs = bb + CORE_Y_STEPS[k].
    let neutrals: [f64; 6] = [
        bb,
        bb + CORE_Y_STEPS[0],
        bb + CORE_Y_STEPS[1],
        bb + CORE_Y_STEPS[2],
        bb + CORE_Y_STEPS[3],
        bb + CORE_Y_STEPS[4],
    ];
    let sorted = peaks.windows(2).all(|w| w[0].0 <= w[1].0);
    let mut score = 0.0f64;
    for &neutral in &neutrals {
        score += best_frag_intensity(peaks, sorted, neutral, tol_ppm, max_charge) as f64 / base;
    }
    score
}

/// Composition-specific glycan Y-ladder intensity match.
///
/// Unlike [`core_y_intensity`] (the glycan-INDEPENDENT trimannosyl core, shared
/// by every N-glycan), this walks the ASSIGNED composition's Y-ion ladder — Y0
/// then each cumulative partial-glycan mass in biosynthetic-ish order (chitobiose
/// core → trimannose → remaining antennae → fucose → sialic acids) — and sums the
/// base-peak-normalised intensity of the matched ions. Because the intermediate
/// Y-ions depend on the composition, a wrong (decoy) glycan of the SAME total
/// mass produces a DIFFERENT ladder and scores lower. That difference is what
/// lets a glycan-axis decoy discriminate on the glycan axis (2D FDR).
///
/// `bb_neutral` is the NEUTRAL peptide backbone (Y0 = bb_neutral + PROTON).
pub fn glycan_y_intensity(
    peaks: &[(f64, f32)],
    bb_neutral: f64,
    comp: &crate::glycan_db::GlycanComp,
    tol_ppm: f64,
    max_charge: u8,
) -> f64 {
    use crate::glycan_mass::{FUC, HEX, HEXNAC, NEUAC, NEUGC};
    let base = peaks
        .iter()
        .map(|&(_, i)| i)
        .fold(0.0f32, f32::max)
        .max(1e-9) as f64;

    // Cumulative monosaccharide additions in a canonical order (composition-specific).
    let core_hexnac = comp.hexnac.min(2);
    let core_hex = comp.hex.min(3);
    let mut adds: Vec<f64> = Vec::new();
    for _ in 0..core_hexnac {
        adds.push(HEXNAC);
    }
    for _ in 0..core_hex {
        adds.push(HEX);
    }
    for _ in 0..(comp.hexnac - core_hexnac) {
        adds.push(HEXNAC);
    }
    for _ in 0..(comp.hex - core_hex) {
        adds.push(HEX);
    }
    for _ in 0..comp.fuc {
        adds.push(FUC);
    }
    for _ in 0..comp.neuac {
        adds.push(NEUAC);
    }
    for _ in 0..comp.neugc {
        adds.push(NEUGC);
    }

    let sorted = peaks.windows(2).all(|w| w[0].0 <= w[1].0);
    // Match at charges 1..=max_charge (multiply-charged Y-ions), by NEUTRAL mass.
    let match_int = |neutral: f64| -> f64 {
        best_frag_intensity(peaks, sorted, neutral, tol_ppm, max_charge) as f64 / base
    };

    let mut score = match_int(bb_neutral); // Y0 (bare backbone), neutral = bb_neutral
    let mut cum = 0.0;
    for m in adds {
        cum += m;
        score += match_int(bb_neutral + cum);
    }
    score
}

/// Canonical cumulative monosaccharide add-order shared by the ladder scorers:
/// core HexNAc(≤2), core Hex(≤3), then the remaining HexNAc/Hex/Fuc/NeuAc/NeuGc.
/// Factored out so the full and intermediate ladders stay bit-for-bit in lock-step.
fn glycan_cumulative_adds(comp: &crate::glycan_db::GlycanComp) -> Vec<f64> {
    use crate::glycan_mass::{FUC, HEX, HEXNAC, NEUAC, NEUGC};
    let core_hexnac = comp.hexnac.min(2);
    let core_hex = comp.hex.min(3);
    let mut adds: Vec<f64> = Vec::new();
    for _ in 0..core_hexnac {
        adds.push(HEXNAC);
    }
    for _ in 0..core_hex {
        adds.push(HEX);
    }
    for _ in 0..(comp.hexnac - core_hexnac) {
        adds.push(HEXNAC);
    }
    for _ in 0..(comp.hex - core_hex) {
        adds.push(HEX);
    }
    for _ in 0..comp.fuc {
        adds.push(FUC);
    }
    for _ in 0..comp.neuac {
        adds.push(NEUAC);
    }
    for _ in 0..comp.neugc {
        adds.push(NEUGC);
    }
    adds
}

/// INTERMEDIATE glycan-Y ladder (Y2 and beyond) — the composition-DISCRIMINATING part
/// of the ladder, for the glycan FDR axis.
///
/// [`glycan_y_intensity`] sums Y0 + every cumulative rung, but Y0 (bare backbone) and
/// Y1 (peptide + first core HexNAc) are shared by EVERY N-glycan AND preserved by the
/// glycan-decoy, and they carry most of the intensity — so they DOMINATE the summed
/// score and wash out the target−decoy difference (measured target/decoy separation of
/// the full ladder is only Cohen d≈0.12). Summing ONLY Y2..Yn isolates exactly the
/// rungs where a true glycan and a same-mass decoy diverge; removing the shared Y0/Y1
/// additive term can only raise the target−decoy separation the 2D-FDR glycan axis
/// depends on. Same neutral-mass multi-charge matcher as [`glycan_y_intensity`].
pub fn glycan_intermediate_y_intensity(
    peaks: &[(f64, f32)],
    bb_neutral: f64,
    comp: &crate::glycan_db::GlycanComp,
    tol_ppm: f64,
    max_charge: u8,
) -> f64 {
    let base = peaks
        .iter()
        .map(|&(_, i)| i)
        .fold(0.0f32, f32::max)
        .max(1e-9) as f64;
    let adds = glycan_cumulative_adds(comp);
    let sorted = peaks.windows(2).all(|w| w[0].0 <= w[1].0);
    let match_int =
        |neutral: f64| -> f64 { best_frag_intensity(peaks, sorted, neutral, tol_ppm, max_charge) as f64 / base };
    // Skip Y0 (bb_neutral) and Y1 (bb_neutral + adds[0]); sum Y2..Yn.
    let mut score = 0.0;
    let mut cum = 0.0;
    for (ai, m) in adds.iter().enumerate() {
        cum += m;
        if ai >= 1 {
            score += match_int(bb_neutral + cum);
        }
    }
    score
}

/// Glycan-axis DECOY of [`glycan_intermediate_y_intensity`]: sums the SAME intermediate
/// rungs (Y2..Yn) but each shifted by the deterministic [`glycan_y_intensity_decoy`]
/// offset, so on a real-glycan spectrum the shifted rungs miss the true Y-ions and the
/// decoy sums far less than the target — the per-row target−decoy gap the glycan FDR
/// axis turns into 2D FDR. `seed` (e.g. the glycan index) makes it reproducible.
pub fn glycan_intermediate_y_intensity_decoy(
    peaks: &[(f64, f32)],
    bb_neutral: f64,
    comp: &crate::glycan_db::GlycanComp,
    tol_ppm: f64,
    seed: u64,
) -> f64 {
    let base = peaks
        .iter()
        .map(|&(_, i)| i)
        .fold(0.0f32, f32::max)
        .max(1e-9) as f64;
    let adds = glycan_cumulative_adds(comp);
    let sorted = peaks.windows(2).all(|w| w[0].0 <= w[1].0);
    let match_int = |ion_mz: f64| -> f64 {
        let tol = (ion_mz * tol_ppm / 1e6).max(0.01);
        let (lo, hi) = (ion_mz - tol, ion_mz + tol);
        let best = if sorted {
            let start = peaks.partition_point(|&(m, _)| m < lo);
            peaks[start..]
                .iter()
                .take_while(|&&(m, _)| m <= hi)
                .map(|&(_, i)| i)
                .fold(0.0f32, f32::max)
        } else {
            peaks
                .iter()
                .filter(|&&(m, _)| m >= lo && m <= hi)
                .map(|&(_, i)| i)
                .fold(0.0f32, f32::max)
        };
        (best as f64) / base
    };
    // Skip Y0 and Y1 (kept, shared anchor); sum the shifted Y2..Yn at 1+.
    let mut score = 0.0;
    let mut cum = 0.0;
    for (ai, m) in adds.iter().enumerate() {
        cum += m;
        if ai >= 1 {
            let h = splitmix64(seed ^ (ai as u64).wrapping_mul(0x9E37_79B9_7F4A_7C15));
            let shift = 1.0 + (h % 29_000) as f64 / 1000.0;
            score += match_int(bb_neutral + cum + shift + PROTON);
        }
    }
    score
}

/// Deterministic mixing (splitmix64) → a reproducible per-(seed,rung) offset.
#[inline]
fn splitmix64(mut x: u64) -> u64 {
    x = x.wrapping_add(0x9E37_79B9_7F4A_7C15);
    let mut z = x;
    z = (z ^ (z >> 30)).wrapping_mul(0xBF58_476D_1CE4_E5B9);
    z = (z ^ (z >> 27)).wrapping_mul(0x94D0_49BB_1331_11EB);
    z ^ (z >> 31)
}

/// Glycan-axis DECOY of [`glycan_y_intensity`]: same composition (same TOTAL mass
/// → still matches the precursor), same Y0 and Y1 (the bare-backbone anchor + the
/// first HexNAc, shared by every N-glycan so preserving them keeps the decoy a
/// plausible glycopeptide), but every INTERMEDIATE Y-rung (Y2..Y_{n-1}) is shifted
/// by a deterministic pseudo-random offset in [1, 30] Da (a glyco search engine/an open-search PTM tool
/// decoy-glycan recipe). On a real-glycan spectrum the shifted rungs miss the true
/// Y-ions, so the decoy sums far less matched intensity than the target ladder —
/// that target−decoy gap is the glycan-axis signal Percolator turns into 2D FDR.
///
/// `seed` (e.g. the glycan index) makes the decoy reproducible. `bb_neutral` is
/// the NEUTRAL peptide backbone (Y0 = bb_neutral + PROTON), as in the target.
pub fn glycan_y_intensity_decoy(
    peaks: &[(f64, f32)],
    bb_neutral: f64,
    comp: &crate::glycan_db::GlycanComp,
    tol_ppm: f64,
    seed: u64,
) -> f64 {
    use crate::glycan_mass::{FUC, HEX, HEXNAC, NEUAC, NEUGC};
    let base = peaks
        .iter()
        .map(|&(_, i)| i)
        .fold(0.0f32, f32::max)
        .max(1e-9) as f64;

    // Same canonical cumulative-adds order as glycan_y_intensity.
    let core_hexnac = comp.hexnac.min(2);
    let core_hex = comp.hex.min(3);
    let mut adds: Vec<f64> = Vec::new();
    for _ in 0..core_hexnac {
        adds.push(HEXNAC);
    }
    for _ in 0..core_hex {
        adds.push(HEX);
    }
    for _ in 0..(comp.hexnac - core_hexnac) {
        adds.push(HEXNAC);
    }
    for _ in 0..(comp.hex - core_hex) {
        adds.push(HEX);
    }
    for _ in 0..comp.fuc {
        adds.push(FUC);
    }
    for _ in 0..comp.neuac {
        adds.push(NEUAC);
    }
    for _ in 0..comp.neugc {
        adds.push(NEUGC);
    }

    let sorted = peaks.windows(2).all(|w| w[0].0 <= w[1].0);
    let match_int = |ion_mz: f64| -> f64 {
        let tol = (ion_mz * tol_ppm / 1e6).max(0.01);
        let (lo, hi) = (ion_mz - tol, ion_mz + tol);
        let best = if sorted {
            let start = peaks.partition_point(|&(m, _)| m < lo);
            peaks[start..]
                .iter()
                .take_while(|&&(m, _)| m <= hi)
                .map(|&(_, i)| i)
                .fold(0.0f32, f32::max)
        } else {
            peaks
                .iter()
                .filter(|&&(m, _)| m >= lo && m <= hi)
                .map(|&(_, i)| i)
                .fold(0.0f32, f32::max)
        };
        (best as f64) / base
    };

    // Y0 (bare backbone) is always kept.
    let mut score = match_int(bb_neutral + PROTON);
    let mut cum = 0.0;
    // `ai` is the add index: ai==0 → Y1 (kept, core anchor); ai>=1 → shifted.
    for (ai, m) in adds.iter().enumerate() {
        cum += m;
        // The LAST cumulative add reconstructs the full glycan (= precursor − bb),
        // which is the precursor, not a fragment; skip it as glycan_y_intensity's
        // ladder does implicitly (it walks all adds, but the top rung coincides
        // with the precursor and matches for target and decoy alike — so to make
        // the decoy discriminate we only need the INTERMEDIATE rungs shifted).
        let shift = if ai == 0 {
            0.0 // keep Y1
        } else {
            // deterministic offset in [1, 30] Da at mDa resolution
            let h = splitmix64(seed ^ (ai as u64).wrapping_mul(0x9E37_79B9_7F4A_7C15));
            1.0 + (h % 29_000) as f64 / 1000.0
        };
        score += match_int(bb_neutral + cum + shift + PROTON);
    }
    score
}

/// G2 Y0/Y1 peptide-mass ANCHOR score.
///
/// Y0 (the bare peptide backbone) and Y1 (peptide + one innermost HexNAc) are the
/// ONE peptide-mass-conditioned signal in a stepped-HCD N-glyco spectrum: they
/// are high-intensity even when the interior b/y ladder is dead, and — unlike the
/// backbone/spectrum-level oxonium / core-Y / glycan-mass features (which are
/// identical for two peptides competing at the same backbone window) — their m/z
/// is a function of the PEPTIDE mass, so they discriminate competing peptides.
/// This is the SP-B lever (see docs/plans/glyco/20-theory/why-andes-fails-and-succeed.md §2).
///
/// Returns the summed base-peak-normalised intensity of Y0 and Y1 matched at
/// charges `1..=min(prec_z, 2)`. `bb_neutral` is the NEUTRAL peptide backbone
/// mass (Y0 neutral = `bb_neutral`; observed 1+ = `bb_neutral + PROTON`).
///
/// It is an ADDITIVE PIN feature only — never fold it into the ranking score.
pub fn y0y1_anchor_intensity(
    peaks: &[(f64, f32)],
    bb_neutral: f64,
    prec_z: u8,
    tol_ppm: f64,
) -> f64 {
    use crate::glycan_mass::HEXNAC;
    let base = peaks
        .iter()
        .map(|&(_, i)| i)
        .fold(0.0f32, f32::max)
        .max(1e-9) as f64;
    let sorted = peaks.windows(2).all(|w| w[0].0 <= w[1].0);
    let match_int = |ion_mz: f64| -> f64 {
        let tol = (ion_mz * tol_ppm / 1e6).max(0.01);
        let (lo, hi) = (ion_mz - tol, ion_mz + tol);
        let best = if sorted {
            let start = peaks.partition_point(|&(m, _)| m < lo);
            peaks[start..]
                .iter()
                .take_while(|&&(m, _)| m <= hi)
                .map(|&(_, i)| i)
                .fold(0.0f32, f32::max)
        } else {
            peaks
                .iter()
                .filter(|&&(m, _)| m >= lo && m <= hi)
                .map(|&(_, i)| i)
                .fold(0.0f32, f32::max)
        };
        (best as f64) / base
    };

    let zmax = prec_z.max(1).min(2);
    let y0_neutral = bb_neutral;
    let y1_neutral = bb_neutral + HEXNAC;
    let mut score = 0.0;
    for z in 1..=zmax {
        let zf = z as f64;
        score += match_int((y0_neutral + zf * PROTON) / zf);
        score += match_int((y1_neutral + zf * PROTON) / zf);
    }
    score
}

pub fn count_core_y_hits(peaks: &[(f64, f32)], bb: f64, tol_ppm: f64, max_charge: u8) -> u8 {
    // Core-Y NEUTRAL masses (Y0..Y5): Y0 = bb, rungs = bb + CORE_Y_STEPS[k].
    let neutrals: [f64; 6] = [
        bb,
        bb + CORE_Y_STEPS[0],
        bb + CORE_Y_STEPS[1],
        bb + CORE_Y_STEPS[2],
        bb + CORE_Y_STEPS[3],
        bb + CORE_Y_STEPS[4],
    ];
    let sorted = peaks.windows(2).all(|w| w[0].0 <= w[1].0);
    let mut hits: u8 = 0;
    for &neutral in &neutrals {
        // A rung hits if a matching peak exists at ANY charge 1..=max_charge.
        if best_frag_intensity(peaks, sorted, neutral, tol_ppm, max_charge) > 0.0 {
            hits += 1;
        }
    }
    hits
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn core_y_matches_multiply_charged_ions() {
        // Stepped-HCD glyco Y-ions are frequently 2+/3+; matching only +1 misses them.
        let bb = 2000.0_f64; // water-included neutral backbone
        let mz2 = |neutral: f64| (neutral + 2.0 * PROTON) / 2.0; // 2+ observed m/z
        let mut peaks = vec![
            (mz2(bb), 100.0f32),                    // Y0 as 2+
            (mz2(bb + CORE_Y_STEPS[0]), 100.0f32),  // Y1 as 2+
            (mz2(bb + CORE_Y_STEPS[1]), 100.0f32),  // Y2 as 2+
        ];
        peaks.sort_by(|a, b| a.0.total_cmp(&b.0));
        // +1-only matching misses all three 2+ ions.
        assert_eq!(count_core_y_hits(&peaks, bb, 20.0, 1), 0, "z=1-only must miss 2+ Y-ions");
        // Matching up to +2 recovers them.
        assert_eq!(count_core_y_hits(&peaks, bb, 20.0, 2), 3, "z<=2 must match the 2+ Y-ions");
        assert!(core_y_intensity(&peaks, bb, 20.0, 2) > 0.0);
    }

    #[test]
    fn solve_finds_backbone_with_minimal_2hexnac_core_glycan() {
        // Off-by-H2O gate regression (audit finding): a de-novo backbone whose true
        // glycan is the minimal N-glycan core (2×HexNAc = 406.16 Da) must be found.
        // The gate compared the water-INCLUDED bb against residue-convention
        // precursor_neutral, undercounting the implied glycan by H2O → 406 rejected.
        use crate::glycan_mass::HEXNAC;
        let bb_neutral = 1500.0_f64; // Y0-derived, water-included
        let bb_res = bb_neutral - H2O; // residue backbone
        let glycan = 2.0 * HEXNAC; // 406.159, the minimal N-glycan core
        let precursor_neutral = bb_res + glycan; // residue convention
        let mut peaks = vec![
            (bb_neutral + PROTON, 1000.0f32),               // Y0
            (bb_neutral + HEXNAC + PROTON, 900.0f32),       // Y1
            (bb_neutral + 2.0 * HEXNAC + PROTON, 500.0f32), // Y2 (= intact-ish)
        ];
        peaks.sort_by(|a, b| a.0.total_cmp(&b.0));
        let cands = solve_backbone_min(&peaks, precursor_neutral, 2, 20.0, 50, 2);
        assert!(
            cands
                .iter()
                .any(|c| (c.backbone_mass - bb_neutral).abs() < 0.05),
            "minimal 2xHexNAc-core backbone dropped by the off-by-H2O gate; got {:?}",
            cands.iter().map(|c| c.backbone_mass).collect::<Vec<_>>()
        );
    }

    #[test]
    fn solve_backbone_recovers_known_mass_from_core_ladder() {
        let bb = 1500.0;
        let mut peaks: Vec<(f64, f32)> = crate::glycan_mass::CORE_Y_STEPS
            .iter()
            .map(|&s| (bb + s + crate::glycan_mass::PROTON, 50.0))
            .collect();
        peaks.push((bb + crate::glycan_mass::PROTON, 40.0)); // Y0 = bare backbone+H
        peaks.push((999.9, 5.0)); // noise
        let precursor = bb + 1444.53; // backbone + a HexNAc2Hex5 glycan (~1444.53)
        let out = solve_backbone(&peaks, precursor, 2, 20.0, 5);
        assert!(!out.is_empty());
        // dedup_by merges candidates within 0.05 Da; the surviving rep may be up to 0.05 Da off
        assert!((out[0].backbone_mass - bb).abs() < 0.05, "got {}", out[0].backbone_mass);
        assert!(out[0].core_y_hits >= 2);
    }

    #[test]
    fn solve_backbone_empty_without_core_quorum() {
        let peaks = vec![(700.0, 50.0_f32), (1234.5, 50.0)]; // no core-Y ladder
        assert!(solve_backbone(&peaks, 2500.0, 2, 20.0, 5).is_empty());
    }

    /// Codex re-review #3: the relaxed quorum-1 rescue must NOT accept a lone
    /// unrelated peak as a spurious single-rung backbone. Without b/y complement
    /// support, a single core-Y hit is not enough — otherwise noise peaks become
    /// candidates that enter peptide scoring.
    /// Composition-specific glycan Y-ladder: the intensity match must be high
    /// for the TRUE composition (whose stepwise Y-ions line up with the spectrum)
    /// and lower for a DECOY composition of the same total mass (different
    /// intermediate Y-ions). This is what lets a glycan-axis decoy discriminate.
    #[test]
    fn glycan_y_intensity_prefers_true_composition_over_same_mass_decoy() {
        use crate::glycan_db::GlycanComp;
        use crate::glycan_mass::{HEX, HEXNAC};
        let proton = PROTON;
        let bb = 1500.0_f64; // neutral backbone

        // True composition HexNAc2Hex3 (trimannosyl core).
        let truth = GlycanComp { hexnac: 2, hex: 3, fuc: 0, neuac: 0, neugc: 0,
                                 mass: 2.0 * HEXNAC + 3.0 * HEX };
        // Build the spectrum from the TRUE composition's Y-ladder.
        let mut cum = 0.0;
        let mut peaks: Vec<(f64, f32)> = vec![(bb + proton, 1000.0)]; // Y0
        for m in [HEXNAC, HEXNAC, HEX, HEX, HEX] {
            cum += m;
            peaks.push((bb + cum + proton, 500.0));
        }
        peaks.push((277.0, 20.0)); // noise
        peaks.sort_by(|a, b| a.0.partial_cmp(&b.0).unwrap());

        // Decoy of the SAME total mass but different composition order:
        // HexNAc3Hex2Fuc? keep same mass by construction is hard; instead use a
        // composition whose intermediate Y-ions differ (Hex-first is not built,
        // but HexNAc5 has different steps). Use HexNAc1Hex4 (~same count, diff steps).
        let decoy = GlycanComp { hexnac: 1, hex: 4, fuc: 0, neuac: 0, neugc: 0,
                                 mass: 1.0 * HEXNAC + 4.0 * HEX };

        let s_true = glycan_y_intensity(&peaks, bb, &truth, 20.0, 3);
        let s_decoy = glycan_y_intensity(&peaks, bb, &decoy, 20.0, 3);
        assert!(s_true > s_decoy,
            "true composition Y-ladder ({s_true}) must beat the decoy ({s_decoy})");
        assert!(s_true > 0.0, "true composition must have positive Y-ladder intensity");
    }

    /// Intermediate (Y2+) ladder: (a) EXCLUDES Y0/Y1 — a spectrum carrying only Y0/Y1
    /// scores exactly 0; (b) fires on the Y2.. rungs; (c) its shifted decoy scores
    /// below the target on a true-glycan spectrum (the isolated glycan-axis signal).
    #[test]
    fn glycan_intermediate_y_excludes_anchor_and_decoy_scores_below_target() {
        use crate::glycan_db::GlycanComp;
        use crate::glycan_mass::{HEX, HEXNAC};
        let bb = 1500.0_f64;
        let truth = GlycanComp { hexnac: 2, hex: 3, fuc: 0, neuac: 0, neugc: 0,
                                 mass: 2.0 * HEXNAC + 3.0 * HEX };
        // Spectrum with ONLY Y0 and Y1 (bare backbone + first HexNAc) → no Y2+ evidence.
        let anchor_only: Vec<(f64, f32)> =
            vec![(bb + PROTON, 1000.0), (bb + HEXNAC + PROTON, 800.0), (277.0, 20.0)];
        let mut ao = anchor_only.clone();
        ao.sort_by(|a, b| a.0.partial_cmp(&b.0).unwrap());
        assert_eq!(
            glycan_intermediate_y_intensity(&ao, bb, &truth, 20.0, 3), 0.0,
            "intermediate ladder must exclude the shared Y0/Y1 anchor"
        );
        // Full ladder present (Y0..Y5) → intermediate fires and its decoy scores below.
        let mut cum = 0.0;
        let mut peaks: Vec<(f64, f32)> = vec![(bb + PROTON, 1000.0)];
        for m in [HEXNAC, HEXNAC, HEX, HEX, HEX] {
            cum += m;
            peaks.push((bb + cum + PROTON, 500.0));
        }
        peaks.sort_by(|a, b| a.0.partial_cmp(&b.0).unwrap());
        let s_t = glycan_intermediate_y_intensity(&peaks, bb, &truth, 20.0, 3);
        let s_d = glycan_intermediate_y_intensity_decoy(&peaks, bb, &truth, 20.0, 42);
        assert!(s_t > 0.0, "intermediate ladder must fire on Y2+ evidence, got {s_t}");
        assert!(s_t > s_d, "intermediate target ({s_t}) must beat its shifted decoy ({s_d})");
    }

    /// G2 Y0/Y1 peptide-mass ANCHOR: the score must be (a) high when the Y0 (bare
    /// peptide) and Y1 (peptide+HexNAc) ions for THIS backbone are present, (b)
    /// zero when they are absent, and — the key property SPA2 features lack — (c)
    /// PEPTIDE-MASS-DISCRIMINATING: for a spectrum whose Y0/Y1 belong to backbone
    /// A, `anchor(A) > anchor(B)` for a different-mass backbone B. This is what
    /// lets it separate competing peptides at one backbone window.
    #[test]
    fn y0y1_anchor_scores_and_discriminates() {
        use crate::glycan_mass::HEXNAC;
        let bb_a = 1500.0_f64; // neutral peptide mass A
        let bb_b = 1650.0_f64; // a different peptide mass B
        // Spectrum carries Y0/Y1 (1+ and 2+) for backbone A only, plus noise.
        let peaks: Vec<(f64, f32)> = vec![
            (bb_a + PROTON, 1000.0),                       // Y0 1+
            (bb_a + HEXNAC + PROTON, 800.0),               // Y1 1+
            ((bb_a + 2.0 * PROTON) / 2.0, 500.0),          // Y0 2+
            (204.087, 300.0),                              // oxonium noise
        ];
        let a = y0y1_anchor_intensity(&peaks, bb_a, 2, 20.0);
        let b = y0y1_anchor_intensity(&peaks, bb_b, 2, 20.0);
        assert!(a > 0.0, "anchor must fire for the true backbone, got {a}");
        assert!(a > b, "anchor must discriminate: A ({a}) > B ({b})");

        // No Y0/Y1 anywhere → zero.
        let empty = vec![(204.087, 100.0), (366.14, 80.0)];
        assert_eq!(y0y1_anchor_intensity(&empty, bb_a, 2, 20.0), 0.0);
    }

    /// G3 glycan-axis decoy: a DECOY glycan ladder keeps Y0/Y1 (core anchor +
    /// precursor match) but shifts the intermediate Y-rungs, so on a TRUE-glycan
    /// spectrum it matches far less intensity than the target ladder. That gap is
    /// the glycan-axis discrimination Percolator needs for 2D FDR. It must also be
    /// deterministic (same seed → same decoy) for reproducibility.
    #[test]
    fn glycan_y_intensity_decoy_scores_below_target_on_true_spectrum() {
        use crate::glycan_db::GlycanComp;
        use crate::glycan_mass::{HEX, HEXNAC};
        let proton = PROTON;
        let bb = 1500.0_f64;
        let truth = GlycanComp { hexnac: 2, hex: 3, fuc: 0, neuac: 0, neugc: 0,
                                 mass: 2.0 * HEXNAC + 3.0 * HEX };
        // Spectrum built from the TRUE composition's full Y-ladder.
        let mut cum = 0.0;
        let mut peaks: Vec<(f64, f32)> = vec![(bb + proton, 1000.0)]; // Y0
        for m in [HEXNAC, HEXNAC, HEX, HEX, HEX] {
            cum += m;
            peaks.push((bb + cum + proton, 500.0));
        }
        peaks.sort_by(|a, b| a.0.partial_cmp(&b.0).unwrap());

        let s_true = glycan_y_intensity(&peaks, bb, &truth, 20.0, 3);
        let s_decoy = glycan_y_intensity_decoy(&peaks, bb, &truth, 20.0, 12345);
        assert!(s_decoy < s_true,
            "decoy ladder ({s_decoy}) must score below the target ({s_true})");
        // Decoy still retains Y0 + Y1 (bb + first HexNAc), so it is > 0 but small.
        assert!(s_decoy > 0.0, "decoy keeps Y0/Y1 anchor, should be > 0");
        // Determinism: same seed → identical score.
        let s_decoy2 = glycan_y_intensity_decoy(&peaks, bb, &truth, 20.0, 12345);
        assert!((s_decoy - s_decoy2).abs() < 1e-12, "decoy must be deterministic");
    }

    /// Convention + intensity regression: the core-Y ladder ions live at the
    /// NEUTRAL peptide mass (Y0 = neutral + PROTON). `count_core_y_hits` and the
    /// new `core_y_intensity` must both be given the NEUTRAL backbone; feeding
    /// the residue mass (H2O too low) misses the ladder. `core_y_intensity` sums
    /// matched base-peak-normalised intensity instead of counting.
    #[test]
    fn core_y_intensity_and_count_use_neutral_convention() {
        let proton = PROTON;
        let steps = CORE_Y_STEPS;
        let neutral = 1500.0_f64; // neutral peptide backbone (Y0 = neutral + proton)
        let mut peaks: Vec<(f64, f32)> = vec![(neutral + proton, 1000.0)]; // Y0, base peak
        for &s in steps.iter() {
            peaks.push((neutral + s + proton, 500.0));
        }
        peaks.push((321.0, 10.0)); // noise
        peaks.sort_by(|a, b| a.0.partial_cmp(&b.0).unwrap());

        // NEUTRAL input → full ladder (6 rungs), positive intensity.
        assert_eq!(count_core_y_hits(&peaks, neutral, 20.0, 3), 6, "neutral must match all 6 rungs");
        let inten = core_y_intensity(&peaks, neutral, 20.0, 3);
        assert!(inten > 0.0, "core_y_intensity must be positive for a real ladder, got {inten}");
        // Y0 (1.0 normalised) + 5×0.5 = 3.5 expected.
        assert!((inten - 3.5).abs() < 1e-6, "expected 3.5 normalised intensity, got {inten}");

        // RESIDUE input (neutral − H2O) → the ladder is missed (wrong convention).
        let residue = neutral - H2O;
        assert!(
            count_core_y_hits(&peaks, residue, 20.0, 3) < 6,
            "residue mass (H2O too low) must NOT match the neutral-anchored ladder"
        );
    }

    #[test]
    fn quorum1_rejects_single_rung_without_complement() {
        // Scattered noise peaks: some will coincidentally hit a single Y-rung
        // bin at various masses, but none form a b/y complement pair.
        let peaks = vec![
            (700.0, 50.0_f32),
            (911.3, 40.0),
            (1301.7, 45.0),
            (1502.9, 30.0),
        ];
        let out = solve_backbone_min(&peaks, 3000.0, 2, 20.0, 50, 1);
        assert!(
            out.is_empty(),
            "quorum-1 must reject single-rung candidates lacking b/y complement, got {out:?}"
        );
    }

    /// The reported backbone mass must be the vote-weighted CENTROID of a merged
    /// near-mass cluster, not its lowest-mass edge. We build a Y-ladder anchored
    /// at a backbone whose mass is offset by a non-integer amount from the 0.01-Da
    /// bin grid; the tolerance window spreads votes symmetrically across several
    /// bins. With a centroid, the recovered mass must land within ~1 mDa of the
    /// true center; a lowest-edge tiebreak would bias it low by the tol window.
    #[test]
    fn solve_backbone_reports_cluster_centroid_not_low_edge() {
        let proton = crate::glycan_mass::PROTON;
        let steps = crate::glycan_mass::CORE_Y_STEPS;
        // True backbone chosen to sit between 0.01-Da bins (…x.xx5 region).
        let bb = 1500.005_f64;
        let mut peaks: Vec<(f64, f32)> = vec![(bb + proton, 50.0)]; // Y0
        for &s in steps.iter() {
            peaks.push((bb + s + proton, 50.0));
        }
        let precursor = bb + 1444.53;
        let out = solve_backbone(&peaks, precursor, 2, 20.0, 5);
        assert!(!out.is_empty());
        // Centroid must be within 2 mDa of the true backbone (unbiased).
        // The old lowest-edge behaviour skewed low by up to the 0.03-Da window.
        assert!(
            (out[0].backbone_mass - bb).abs() < 0.002,
            "centroid off by {:.4} (got {})",
            (out[0].backbone_mass - bb).abs(),
            out[0].backbone_mass
        );
    }

    /// Regression test: solve_backbone must return the same ordered result on repeated calls.
    /// We build two synthetic peptide+core-Y ladders with backbone masses 1200.0 and 1800.0 Da,
    /// each yielding the same number of core-Y hits and votes, producing a tie on
    /// (core_y_hits, intensity_score). Without the backbone_mass tiebreaker the HashMap
    /// iteration order could permute the two candidates across runs.
    #[test]
    fn solve_backbone_is_deterministic_under_ties() {
        let proton = crate::glycan_mass::PROTON;
        let steps = crate::glycan_mass::CORE_Y_STEPS;

        // Build identical core-Y ladders anchored at two different backbone masses.
        // Both backbones get exactly the same number of Y-ladder peaks → tied hits & scores.
        let bb1 = 1200.0_f64;
        let bb2 = 1800.0_f64;

        let mut peaks: Vec<(f64, f32)> = Vec::new();
        for &bb in &[bb1, bb2] {
            // Y0 (bare backbone)
            peaks.push((bb + proton, 50.0));
            // Y_r for each core step
            for &s in steps.iter() {
                peaks.push((bb + s + proton, 50.0));
            }
        }

        // Precursor must be larger than any peak neutral mass.
        let precursor_neutral = bb2 + steps.last().copied().unwrap_or(0.0) + 500.0;

        // Call solve_backbone twice with identical args.
        let run1 = solve_backbone(&peaks, precursor_neutral, 2, 20.0, 5);
        let run2 = solve_backbone(&peaks, precursor_neutral, 2, 20.0, 5);

        assert!(!run1.is_empty(), "expected candidates in run1");
        assert_eq!(run1.len(), run2.len(), "result length differs between runs");
        for (c1, c2) in run1.iter().zip(run2.iter()) {
            assert_eq!(
                c1.core_y_hits, c2.core_y_hits,
                "core_y_hits differ between runs"
            );
            assert_eq!(c1.votes, c2.votes, "votes differ between runs");
            assert!(
                (c1.backbone_mass - c2.backbone_mass).abs() < 1e-9,
                "backbone_mass differs between runs: {} vs {}",
                c1.backbone_mass,
                c2.backbone_mass
            );
        }
    }

    /// Complement-pair (b/y) confirmation: a candidate WITH complement-pair support
    /// must outrank an otherwise-equal candidate WITHOUT. We build two synthetic
    /// backbones at the same precursor:
    ///  - `true_bb`: full core-Y ladder (same rung count and intensity as spurious)
    ///    PLUS synthetic b/y fragment ion pairs that satisfy b+y ≈ bb+20.025118.
    ///  - `spurious_bb`: same core-Y ladder (identical rung count and intensity)
    ///    but NO complement b/y pairs.
    /// After solve_backbone, `true_bb` must be first (complement_score tips the tie).
    #[test]
    fn solve_backbone_complement_pairs_promote_true_backbone() {
        let proton = crate::glycan_mass::PROTON;
        let steps = crate::glycan_mass::CORE_Y_STEPS;

        // Two backbones at different masses; same number of core-Y rungs (6 each)
        // and same intensity — so core_y_hits and intensity_score are tied.
        // Use intensity 500.0 for both so sqrt-norm is equal.
        let true_bb = 1500.0_f64;
        let spurious_bb = 1300.0_f64;

        let intensity = 500.0_f32;
        let mut peaks: Vec<(f64, f32)> = Vec::new();

        // Core-Y ladder for true_bb
        peaks.push((true_bb + proton, intensity)); // Y0
        for &s in steps.iter() {
            peaks.push((true_bb + s + proton, intensity));
        }

        // Core-Y ladder for spurious_bb (identical count and intensity)
        peaks.push((spurious_bb + proton, intensity)); // Y0
        for &s in steps.iter() {
            peaks.push((spurious_bb + s + proton, intensity));
        }

        // Add complement pairs only for true_bb.
        // A valid b-ion at b_mz and its complement y at y_mz = (true_bb + BY_COMPLEMENT_OFFSET) - b_mz.
        // Choose b_mz values far from the core-Y peaks to avoid confusion.
        let comp_target = true_bb + BY_COMPLEMENT_OFFSET; // 1500 + 20.025118 = 1520.025118
        // Three b/y pairs: b at 300, 450, 600 Da → y = comp_target - b
        for &b_mz in &[300.0_f64, 450.0, 600.0] {
            let y_mz = comp_target - b_mz;
            // Only add if y_mz > 0 and distinct from existing peaks
            if y_mz > 0.0 {
                peaks.push((b_mz, intensity));
                peaks.push((y_mz, intensity));
            }
        }

        // Precursor must be large enough: true_bb + glycan (1444.53 ≥ MIN_GLYCAN=406)
        let precursor_neutral = true_bb + 1444.53;

        // Both spurious_bb (1300) and true_bb (1500) are within [MIN_BB, precursor_neutral)
        // and leave room for MIN_GLYCAN:
        //   precursor_neutral - spurious_bb = 1500+1444.53 - 1300 = 1644.53 ≥ 406 ✓
        let out = solve_backbone(&peaks, precursor_neutral, 2, 20.0, 5);

        assert!(!out.is_empty(), "expected candidates");
        // true_bb must be ranked first because its complement_score > 0 while spurious_bb = 0.
        assert!(
            (out[0].backbone_mass - true_bb).abs() < 0.05,
            "expected true_bb ({}) first, got backbone_mass={:.4} complement_score={:.4}",
            true_bb,
            out[0].backbone_mass,
            out[0].complement_score,
        );
        assert!(
            out[0].complement_score > 0.0,
            "true_bb must have positive complement_score, got {}",
            out[0].complement_score
        );
    }

    /// High-intensity true Y-ions should rank above equal-rung-count low-intensity noise.
    /// Build a known backbone ladder with high intensity, and a spurious same-rung-count
    /// cluster built from very low-intensity peaks. The true backbone must rank first.
    #[test]
    fn solve_backbone_intensity_weighted_prefers_bright_ladder() {
        let proton = crate::glycan_mass::PROTON;
        let steps = crate::glycan_mass::CORE_Y_STEPS;

        // True backbone at 1500 Da — bright peaks
        let true_bb = 1500.0_f64;
        let mut peaks: Vec<(f64, f32)> = vec![(true_bb + proton, 10000.0)]; // Y0 bright
        for &s in steps.iter() {
            peaks.push((true_bb + s + proton, 8000.0));
        }

        // Spurious backbone at 800 Da — dim peaks that by flat-vote coincidence
        // would accumulate the same rung count but much lower intensity.
        let spurious_bb = 800.0_f64;
        peaks.push((spurious_bb + proton, 1.0)); // Y0 dim
        for &s in steps.iter() {
            peaks.push((spurious_bb + s + proton, 1.0));
        }

        let precursor_neutral = true_bb + 1444.53;
        let out = solve_backbone(&peaks, precursor_neutral, 2, 20.0, 5);
        assert!(!out.is_empty());
        // The bright true-backbone cluster must be ranked first.
        assert!(
            (out[0].backbone_mass - true_bb).abs() < 0.05,
            "expected true backbone first, got {}",
            out[0].backbone_mass
        );
    }
}
