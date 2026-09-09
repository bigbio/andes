//! `--precursor-mono`: apply the MS1 isotope-envelope precursor correction to
//! a streamed chunk of MS2 spectra, and keep the per-scan record the glyco PIN
//! and the diagnostic dump read. The decision itself lives in
//! `search::precursor_mono`; this is the driver-side plumbing.

use std::fs::File;
use std::io::{self, BufWriter, Write};
use std::ops::RangeInclusive;
use std::path::Path;

use input::Ms1Link;
use model::Spectrum;
use search::precursor_mono::{
    correct_precursor, median_nonzero_intensity, MonoCorrection, MonoParams,
};

/// The isotope-error window to search once `--precursor-mono auto` has corrected
/// precursors. With corrected precursors the `+2` step of the sweep only carries
/// the (k, X) ≡ (k−1, X + Hex + Fuc − NeuGc) composition degeneracy: on pGlyco2
/// liver T-1 the `0..1` window confirmed 59 of the 62 firmware-mispicked +4
/// reference scans against 55 with `0..2`, at 96.9% vs 96.8% peptidoform
/// agreement and 0.98% vs 0.91% entrapment FDP, and the `+2` tier went from 402
/// accepted PSMs to 21; the same held on all five liver fractions, heart and
/// lung (bigbio/andes#64, arm F). So `auto` narrows the window to `0..=1` —
/// but only when the correction actually fitted something (`mono_fitted`), so
/// MGF and MS1-less mzML stay byte-identical to `off`, and never over an
/// explicit `--isotope-error` or a non-default `--glyco-isotope-error`.
///
/// The window is applied per spectrum by the glyco driver: only spectra whose
/// envelope was fitted (a `Some` in the correction table) take it, the rest keep
/// the configured window (`search::glyco_search::isotope_window_for`).
///
/// Returns the window to search and whether it was narrowed.
pub(crate) fn coupled_isotope_window(
    mono_fitted: bool,
    explicit: Option<(i8, i8)>,
    glyco_flag_is_default: bool,
    current: RangeInclusive<i8>,
) -> (RangeInclusive<i8>, bool) {
    if mono_fitted && explicit.is_none() && glyco_flag_is_default && current != (0..=1) {
        (0..=1, true)
    } else {
        (current, false)
    }
}

/// Running tally of what the correction did, for the end-of-stream log line.
#[derive(Debug, Default)]
pub(crate) struct MonoStats {
    /// MS2 spectra seen while the correction was active.
    pub seen: usize,
    /// MS2 with a linked MS1 and a known charge (fit attempted).
    pub fitted: usize,
    /// Applied shifts, indexed by shift (0 = left as recorded).
    pub by_shift: Vec<usize>,
}

impl MonoStats {
    pub fn shifted(&self) -> usize {
        self.by_shift.iter().skip(1).sum()
    }
}

/// Diagnostic TSV writer (`--precursor-mono-dump`): one row per MS2.
pub(crate) struct MonoDump {
    w: BufWriter<File>,
}

impl MonoDump {
    pub fn create(path: &Path, max_shift: u8) -> io::Result<Self> {
        let mut w = BufWriter::new(File::create(path)?);
        let fits: Vec<String> = (0..=max_shift).map(|k| format!("fit{k}")).collect();
        writeln!(
            w,
            "scan\tcharge\trt_s\trecorded_mz\tcorrected_mz\tshift\tbest_shift\tfit_recorded\t\
             fit_best\tfit_up\tsnr\tn_iso\t{}",
            fits.join("\t")
        )?;
        Ok(Self { w })
    }

    fn row(&mut self, spec: &Spectrum, c: Option<&MonoCorrection>) -> io::Result<()> {
        let scan = spec.scan.unwrap_or(-1);
        let z = spec.precursor_charge.unwrap_or(0);
        let rt = spec.rt_seconds.unwrap_or(0.0);
        match c {
            Some(c) => {
                write!(
                    self.w,
                    "{scan}\t{z}\t{rt:.3}\t{:.6}\t{:.6}\t{}\t{}\t{:.4}\t{:.4}\t{:.4}\t{:.3}\t{}",
                    c.recorded_mz,
                    c.corrected_mz,
                    c.shift,
                    c.best_shift,
                    c.fit_recorded,
                    c.fit_best,
                    c.fit_up,
                    c.snr,
                    c.n_iso
                )?;
                for f in &c.fits {
                    write!(self.w, "\t{f:.4}")?;
                }
                writeln!(self.w)
            }
            None => writeln!(
                self.w,
                "{scan}\t{z}\t{rt:.3}\t{:.6}\t{:.6}\tNA\tNA\tNA\tNA\tNA\tNA\tNA",
                spec.precursor_mz, spec.precursor_mz
            ),
        }
    }

    pub fn finish(mut self) -> io::Result<()> {
        self.w.flush()
    }
}

/// Correct the precursors of one streamed chunk in place, appending one
/// `Option<MonoCorrection>` per spectrum to `table` (so `table` stays aligned
/// with the driver's `all_spectra`). A spectrum with no linked MS1 or no
/// charge gets `None` and is left untouched.
pub(crate) fn correct_chunk(
    chunk: &mut [Spectrum],
    link: &Ms1Link,
    params: &MonoParams,
    table: &mut Vec<Option<MonoCorrection>>,
    stats: &mut MonoStats,
    dump: Option<&mut MonoDump>,
) -> io::Result<()> {
    if stats.by_shift.len() <= params.max_shift as usize {
        stats.by_shift.resize(params.max_shift as usize + 1, 0);
    }
    // The MS1 median noise is a property of the MS1 scan, shared by every MS2
    // it precedes: compute it once per MS1, lazily.
    let mut medians: Vec<Option<Option<f32>>> = vec![None; link.ms1_peaks.len()];
    let mut dump = dump;
    for (i, spec) in chunk.iter_mut().enumerate() {
        stats.seen += 1;
        let ms1_idx = link.ms2_to_ms1.get(i).copied().flatten();
        let charge = spec
            .precursor_charge
            .filter(|&z| z > 0 && z <= u8::MAX as i32)
            .map(|z| z as u8);
        let correction = match (ms1_idx, charge) {
            (Some(mi), Some(z)) => {
                let ms1 = &link.ms1_peaks[mi];
                let med = *medians[mi].get_or_insert_with(|| median_nonzero_intensity(ms1));
                correct_precursor(ms1, spec.precursor_mz, z, med, params)
            }
            _ => None,
        };
        if let Some(c) = &correction {
            stats.fitted += 1;
            stats.by_shift[c.shift as usize] += 1;
            if c.applied() {
                spec.precursor_mz = c.corrected_mz;
            }
        }
        if let Some(d) = dump.as_deref_mut() {
            d.row(spec, correction.as_ref())?;
        }
        table.push(correction);
    }
    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn coupled_window_narrows_only_when_fitted_and_unset() {
        // Fitted, nothing explicit: the glyco default 0..=2 becomes 0..=1.
        assert_eq!(
            coupled_isotope_window(true, None, true, 0..=2),
            (0..=1, true)
        );
        // Nothing fitted (MGF, MS1-less mzML): untouched, byte-identical to `off`.
        assert_eq!(
            coupled_isotope_window(false, None, true, 0..=2),
            (0..=2, false)
        );
        // An explicit --isotope-error is honoured verbatim, whatever it is.
        assert_eq!(
            coupled_isotope_window(true, Some((0, 2)), true, 0..=2),
            (0..=2, false)
        );
        assert_eq!(
            coupled_isotope_window(true, Some((-1, 3)), true, -1..=3),
            (-1..=3, false)
        );
        // --glyco-isotope-error negative / wide are explicit choices too.
        assert_eq!(
            coupled_isotope_window(true, None, false, -1..=2),
            (-1..=2, false)
        );
        assert_eq!(
            coupled_isotope_window(true, None, false, 0..=5),
            (0..=5, false)
        );
        // Already 0..=1: nothing to report.
        assert_eq!(
            coupled_isotope_window(true, None, true, 0..=1),
            (0..=1, false)
        );
    }
    use model::isotope::glycopeptide_isotope_envelope;
    use model::mass::{ISOTOPE, PROTON};

    fn spec(scan: i32, mz: f64, z: Option<i32>) -> Spectrum {
        Spectrum {
            title: format!("scan={scan}"),
            precursor_mz: mz,
            precursor_charge: z,
            scan: Some(scan),
            ..Default::default()
        }
    }

    /// One MS1 carrying a 3 kDa z=2 glycopeptide envelope with its monoisotope
    /// at `mono`.
    fn ms1(mono: f64) -> Vec<(f64, f32)> {
        let env = glycopeptide_isotope_envelope((mono - PROTON) * 2.0, 8);
        let mut p: Vec<(f64, f32)> = (0..30)
            .map(|i| (mono - 50.0 + i as f64 * 1.7, 50.0))
            .collect();
        p.extend(
            env.iter()
                .enumerate()
                .map(|(k, &e)| (mono + k as f64 * ISOTOPE / 2.0, 1e6 * e as f32)),
        );
        p.sort_by(|a, b| a.0.total_cmp(&b.0));
        p
    }

    #[test]
    fn correct_chunk_shifts_only_the_linked_mispicked_scan_and_keeps_the_table_aligned() {
        let mono = 1500.25;
        let link = Ms1Link {
            ms1_peaks: vec![ms1(mono)],
            // scan 1: linked, recorded on M+4; scan 2: linked, on the mono;
            // scan 3: no MS1 link; scan 4: linked but no charge.
            ms2_to_ms1: vec![Some(0), Some(0), None, Some(0)],
        };
        let m4 = mono + 4.0 * ISOTOPE / 2.0;
        let mut chunk = vec![
            spec(1, m4, Some(2)),
            spec(2, mono, Some(2)),
            spec(3, m4, Some(2)),
            spec(4, m4, None),
        ];
        let mut table = Vec::new();
        let mut stats = MonoStats::default();
        correct_chunk(
            &mut chunk,
            &link,
            &MonoParams::default(),
            &mut table,
            &mut stats,
            None,
        )
        .unwrap();
        assert_eq!(table.len(), 4, "one table entry per spectrum");
        // Default back-off: M+4 is searched from M+1 (shift 3).
        let c1 = table[0].as_ref().expect("linked scan is fitted");
        assert_eq!(c1.shift, 3);
        assert!((chunk[0].precursor_mz - (mono + ISOTOPE / 2.0)).abs() < 1e-9);
        assert_eq!(table[1].as_ref().unwrap().shift, 0);
        assert_eq!(
            chunk[1].precursor_mz, mono,
            "a correct precursor is untouched"
        );
        assert!(table[2].is_none(), "no MS1 link → no fit");
        assert_eq!(chunk[2].precursor_mz, m4, "unlinked scan is untouched");
        assert!(table[3].is_none(), "no charge → no fit");
        assert_eq!(chunk[3].precursor_mz, m4);
        assert_eq!(stats.seen, 4);
        assert_eq!(stats.fitted, 2);
        assert_eq!(stats.shifted(), 1);
        assert_eq!(stats.by_shift[3], 1);
    }

    #[test]
    fn dump_writes_one_row_per_spectrum_with_na_for_unfitted() {
        let dir = tempfile::tempdir().unwrap();
        let path = dir.path().join("dump.tsv");
        let mono = 1500.25;
        let link = Ms1Link {
            ms1_peaks: vec![ms1(mono)],
            ms2_to_ms1: vec![Some(0), None],
        };
        let mut chunk = vec![
            spec(7, mono + 2.0 * ISOTOPE, Some(2)),
            spec(8, mono, Some(2)),
        ];
        let mut dump = MonoDump::create(&path, 6).unwrap();
        let mut table = Vec::new();
        let mut stats = MonoStats::default();
        correct_chunk(
            &mut chunk,
            &link,
            &MonoParams::default(),
            &mut table,
            &mut stats,
            Some(&mut dump),
        )
        .unwrap();
        dump.finish().unwrap();
        let text = std::fs::read_to_string(&path).unwrap();
        let lines: Vec<&str> = text.lines().collect();
        assert_eq!(lines.len(), 3, "header + 2 rows:\n{text}");
        assert!(lines[0].starts_with("scan\tcharge\trt_s\trecorded_mz\tcorrected_mz\tshift"));
        assert!(lines[0].ends_with("fit6"));
        let r7: Vec<&str> = lines[1].split('\t').collect();
        assert_eq!(r7[0], "7");
        assert_eq!(r7[5], "3", "M+4 → shift 3 with back-off: {}", lines[1]);
        assert_eq!(r7.len(), 12 + 7, "12 fixed columns + fit0..fit6");
        assert!(lines[2].starts_with("8\t2\t"));
        assert!(
            lines[2].ends_with("\tNA"),
            "unlinked scan writes NA: {}",
            lines[2]
        );
    }
}
