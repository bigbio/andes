//! Calibrations of a search: DIA retention times and MS1 masses, the fragment mass error per file,
//! and the end of the monoisotopic correction.

use model::Spectrum;
use search::fragments::index::{index_chunk_len, index_order_mass};
use search::precursor::monoisotopic::{MonoCorrection, MonoParams};
use search::{PreparedSearch, SearchParams, TopNQueue};

use super::GlycoIsotopeOverride;
use crate::cli::Cli;
use crate::monoisotopic::{MonoDump, MonoStats};

/// Share of a DIA run's spectra searched first, with the plain vote rank, to fit the RT prior.
const DIA_RT_CAL_FRACTION: f64 = 0.02;
/// Fewest calibration spectra (all of them in a smaller run).
const DIA_RT_CAL_MIN_SPECTRA: usize = 2_000;
/// Target-decoy q-value of the calibration winners that become RT anchors.
const DIA_RT_CAL_Q: f64 = 0.01;
const SECS_PER_MIN: f64 = 60.0;

/// Calibration sample sizes in search order: the first sample (`DIA_RT_CAL_FRACTION` of the
/// run, at least `DIA_RT_CAL_MIN_SPECTRA`), then doubling up to every spectrum of the run.
pub(super) fn dia_rt_cal_sizes(n_spectra: usize) -> Vec<usize> {
    let mut n = ((n_spectra as f64 * DIA_RT_CAL_FRACTION).round() as usize)
        .max(DIA_RT_CAL_MIN_SPECTRA)
        .min(n_spectra);
    let mut sizes = vec![n];
    while n < n_spectra {
        n = (2 * n).min(n_spectra);
        sizes.push(n);
    }
    sizes
}

/// Fits the run's RT prior for the DIA vote cut and installs it on `prepared`: a deterministic
/// subsample of the run is searched with the plain vote rank, the top match of each sample
/// spectrum that is a target at q <= `DIA_RT_CAL_Q` becomes an anchor (andes-rt raw prediction
/// of its peptide, scan minutes), and [`search::dia::DiaRtPrior::fit`] maps raw
/// predictions to minutes. While the anchors are fewer than the map needs for its full
/// resolution ([`search::dia::FULL_MAP_ANCHORS`]) the sample grows (sizes of
/// `dia_rt_cal_sizes`, each adding the next spectra of the same hash order), so a short or
/// sparse run gets the anchors of a long one when its spectra hold them. With too few anchors
/// no prior is installed and the cut keeps the vote rank.
pub(super) fn calibrate_dia_rt(prepared: &PreparedSearch, spectra: &[Spectrum], span_da: f64) {
    let Some(model) = prepared.form_rt_model() else {
        return;
    };
    if spectra.is_empty() {
        return;
    }
    let t = std::time::Instant::now();
    let mut order: Vec<(u64, usize)> = (0..spectra.len())
        .map(|i| (splitmix64(i as u64), i))
        .collect();
    order.sort_unstable();
    let mut tops: Vec<(u32, f64)> = Vec::new();
    let mut labels: Vec<search::tdc::ScoredLabel> = Vec::new();
    let mut anchors: Vec<(f64, f64)> = Vec::new();
    let mut searched = 0;
    for n in dia_rt_cal_sizes(spectra.len()) {
        let mut part: Vec<&Spectrum> = order[searched..n]
            .iter()
            .map(|&(_, i)| &spectra[i])
            .collect();
        part.sort_by(|a, b| index_order_mass(a).total_cmp(&index_order_mass(b)));
        let part: Vec<Spectrum> = part.into_iter().cloned().collect();
        let mut start = 0;
        while start < part.len() {
            let take = index_chunk_len(&part[start..], index_order_mass, span_da);
            let chunk = &part[start..start + take];
            for (s, q) in chunk
                .iter()
                .zip(prepared.run_chunk(chunk, searched + start))
            {
                if let Some(p) = q.peek_top() {
                    let idx = p.primary_candidate_idx();
                    let rt = s
                        .rt_seconds
                        .filter(|t| t.is_finite())
                        .map_or(f64::NAN, |t| t / SECS_PER_MIN);
                    tops.push((idx, rt));
                    labels.push(search::tdc::ScoredLabel {
                        score: p.score,
                        is_decoy: prepared.candidate_is_decoy(idx),
                    });
                }
            }
            start += take;
        }
        searched = n;
        let q = search::tdc::qvalues(&labels);
        anchors = tops
            .iter()
            .zip(&labels)
            .zip(&q)
            .filter(|((top, l), q)| !l.is_decoy && **q <= DIA_RT_CAL_Q && top.1.is_finite())
            .map(|((top, _), _)| (f64::from(model(&prepared.candidate_peptide(top.0))), top.1))
            .collect();
        if anchors.len() >= search::dia::FULL_MAP_ANCHORS {
            break;
        }
    }
    match search::dia::DiaRtPrior::fit(&anchors) {
        Some(prior) => {
            log::info!(
                target: "dia-rt-rank",
                "RT prior from {searched} of {} spectra: {} anchors at q <= {DIA_RT_CAL_Q} ({} \
                 within the trimming band), {} knots, sigma {:.3} min; the vote cut ranks by \
                 votes and RT [{:.1}s]",
                spectra.len(),
                prior.anchors(),
                prior.kept(),
                prior.knots(),
                prior.sigma_min(),
                t.elapsed().as_secs_f64()
            );
            prepared.set_dia_rt_prior(Some(prior));
        }
        None => log::info!(
            target: "dia-rt-rank",
            "no RT prior: {} anchors at q <= {DIA_RT_CAL_Q} from {searched} sample spectra (at \
             least {} needed); the vote cut ranks by votes [{:.1}s]",
            anchors.len(),
            search::dia::MIN_ANCHORS,
            t.elapsed().as_secs_f64()
        ),
    }
}

fn splitmix64(x: u64) -> u64 {
    let mut z = x.wrapping_add(0x9E37_79B9_7F4A_7C15);
    z = (z ^ (z >> 30)).wrapping_mul(0xBF58_476D_1CE4_E5B9);
    z = (z ^ (z >> 27)).wrapping_mul(0x94D0_49BB_1331_11EB);
    z ^ (z >> 31)
}

pub(super) fn recalibrate_dia_ms1(
    prepared: &mut PreparedSearch,
    pending: &[Spectrum],
    span_da: f64,
) {
    const SAMPLE: usize = 20_000;
    const WINDOW_PPM: f64 = 25.0;
    let t = std::time::Instant::now();
    let mut refs: Vec<&Spectrum> = pending.iter().filter(|s| s.peaks.len() >= 10).collect();
    refs.sort_by(|a, b| index_order_mass(a).total_cmp(&index_order_mass(b)));
    let n = SAMPLE.min(refs.len());
    let lo = (refs.len() / 2).saturating_sub(n / 2);
    let sample: Vec<Spectrum> = refs[lo..lo + n].iter().map(|s| (*s).clone()).collect();
    drop(refs);
    let mut wide = prepared.params.clone();
    wide.dia_ms1_tol_ppm = WINDOW_PPM;
    let mut queues: Vec<TopNQueue> = Vec::with_capacity(sample.len());
    let mut start = 0;
    while start < sample.len() {
        let take = index_chunk_len(&sample[start..], index_order_mass, span_da);
        queues.extend(prepared.run_chunk_with_params(&sample[start..start + take], start, &wide));
        start += take;
    }
    match search::fragments::error::dia_ms1_offset_ppm(prepared, &sample, &queues, WINDOW_PPM) {
        Some((ppm, n)) => {
            prepared.shift_dia_ms1_ppm(ppm);
            log::info!(
                target: "dia",
                "MS1 m/z offset {ppm:+.2} ppm from {n} confident PSMs; survey peaks recentred ({:.1}s)",
                t.elapsed().as_secs_f64()
            );
        }
        None => log::info!(
            target: "dia",
            "MS1 m/z offset not measured (too few confident PSMs) ({:.1}s)",
            t.elapsed().as_secs_f64()
        ),
    }
}

pub(super) struct FragmentErrorMeasurer {
    pub(super) enabled: bool,
    pub(super) done: bool,
    pub(super) per_file: Vec<Option<scoring_crate::psm::FragmentErrorModel>>,
}

impl FragmentErrorMeasurer {
    pub(super) fn start_file(&mut self, prepared: &PreparedSearch) {
        self.done = false;
        self.per_file.push(None);
        if self.enabled {
            prepared.set_fragment_error(None);
        }
    }

    pub(super) fn measure<'s>(
        &mut self,
        prepared: &PreparedSearch,
        spectra: impl Iterator<Item = &'s Spectrum>,
        index_span_da: Option<f64>,
    ) {
        if !self.enabled || self.done {
            return;
        }
        self.done = true;
        let t = std::time::Instant::now();
        let mut refs: Vec<&Spectrum> = spectra.filter(|s| s.peaks.len() >= 10).collect();
        refs.sort_by(|a, b| index_order_mass(a).total_cmp(&index_order_mass(b)));
        let n = search::fragments::error::SAMPLE_SPECTRA.min(refs.len());
        let lo = (refs.len() / 2).saturating_sub(n / 2);
        let sample: Vec<Spectrum> = refs[lo..lo + n].iter().map(|s| (*s).clone()).collect();
        drop(refs);
        let mut queues: Vec<TopNQueue> = Vec::with_capacity(sample.len());
        match index_span_da {
            Some(span_da) => {
                let mut start = 0;
                while start < sample.len() {
                    let take = index_chunk_len(&sample[start..], index_order_mass, span_da);
                    queues.extend(prepared.run_chunk(&sample[start..start + take], start));
                    start += take;
                }
            }
            None => queues.extend(prepared.run_chunk(&sample, 0)),
        }
        match search::fragments::error::measure_on_queues(prepared, &sample, &queues) {
            Ok(m) => {
                prepared.set_fragment_error(Some(m));
                if let Some(last) = self.per_file.last_mut() {
                    *last = Some(m);
                }
                log::info!(
                    target: "fragment-error",
                    "offset {:+.2} ppm, spread {:.2} ppm from {} ions of {} PSMs \
                     at 1% q in {} sampled spectra -> fragment peaks recentred by {:+.2} ppm, \
                     PIN feature tolerance {:.1} ppm, ppm-Gaussian width {:.2} ppm, \
                     node-score window {:?} ({:.1}s)",
                    m.offset_ppm,
                    m.sigma_ppm,
                    m.ions,
                    m.psms,
                    sample.len(),
                    prepared.scorer.fragment_recentre_ppm(),
                    m.feature_tolerance_ppm(),
                    m.gaussian_sigma_ppm(),
                    prepared.scorer.node_score_window(),
                    t.elapsed().as_secs_f64()
                );
            }
            Err(why) => log::info!(
                target: "fragment-error",
                "not measured ({why}); the PIN features keep {:?} ({:.1}s)",
                prepared.scorer.feature_match_tolerance(),
                t.elapsed().as_secs_f64()
            ),
        }
    }
}

/// Ends the monoisotopic correction of the search: logs what it fitted (or turns it off when no MS2
/// could be fitted), flushes `--precursor-mono-dump`, and narrows the glyco isotope-error window for
/// the spectra with a verified monoisotope. Returns whether the correction stays on and that window.
pub(super) fn finish_mono_correction(
    cli: &Cli,
    params: &SearchParams,
    mono_params: &MonoParams,
    mono_active: bool,
    mono_stats: &MonoStats,
    mono_dump: &mut Option<MonoDump>,
    mono_table: &mut Vec<Option<MonoCorrection>>,
) -> Result<(bool, GlycoIsotopeOverride), Box<dyn std::error::Error>> {
    let mono_active = mono_active && {
        if mono_stats.fitted == 0 {
            log::warn!(
                "--precursor-mono auto: no MS2 could be fitted ({} MS2 seen, none with \
                 a linked MS1 and a known charge); precursors are left as recorded and the \
                 glyco PIN Mono columns are 0.",
                mono_stats.seen
            );
            if let Some(d) = mono_dump.take() {
                d.finish()
                    .map_err(|e| format!("flush --precursor-mono-dump: {e}"))?;
            }
            mono_table.clear();
        }
        mono_stats.fitted > 0
    };
    if mono_active {
        let hist: Vec<String> = mono_stats
            .by_shift
            .iter()
            .enumerate()
            .skip(1)
            .filter(|(_, &n)| n > 0)
            .map(|(k, n)| format!("{k}:{n}"))
            .collect();
        log::info!(
            target: "precursor-mono",
            "{} of {} MS2 fitted (linked MS1 + known charge), {} shifted \
             down [shift:count {}]",
            mono_stats.fitted,
            mono_stats.seen,
            mono_stats.shifted(),
            if hist.is_empty() {
                "-".to_string()
            } else {
                hist.join(" ")
            }
        );
        if let Some(d) = mono_dump.take() {
            d.finish()
                .map_err(|e| format!("flush --precursor-mono-dump: {e}"))?;
        }
    }
    let (window, narrowed) = crate::monoisotopic::coupled_isotope_window(
        mono_active,
        cli.isotope_error,
        params.isotope_error_range.clone(),
    );
    let glyco_isotope_override = if narrowed {
        let mask: Vec<bool> = mono_table
            .iter()
            .map(|c| {
                c.as_ref()
                    .is_some_and(|c| c.verified_within(mono_params.min_fit, *window.end()))
            })
            .collect();
        log::info!(
            target: "precursor-mono",
            "isotope-error window {}..={} -> {}..={} for the {} of {} MS2 \
             with a verified monoisotope (the others keep {}..={}; a verified monoisotope \
             needs no +2 step, which only admits the Hex+Fuc/NeuGc composition degeneracy — \
             pass --isotope-error to override)",
            params.isotope_error_range.start(),
            params.isotope_error_range.end(),
            window.start(),
            window.end(),
            mask.iter().filter(|&&v| v).count(),
            mono_stats.seen,
            params.isotope_error_range.start(),
            params.isotope_error_range.end()
        );
        Some((window, mask))
    } else {
        None
    };
    Ok((mono_active, glyco_isotope_override))
}
