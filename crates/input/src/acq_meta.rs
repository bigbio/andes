use std::io::BufRead;

use quick_xml::{events::Event, Reader};

const CV_COLLISION_ENERGY: &str = "MS:1000045";
const CV_NORMALIZED_COLLISION_ENERGY: &str = "MS:1000138";
const CV_MS_LEVEL: &str = "MS:1000511";
const CV_ISOLATION_LOWER_OFFSET: &str = "MS:1000828";
const CV_ISOLATION_UPPER_OFFSET: &str = "MS:1000829";
const CV_FAIMS_CV: &str = "MS:1001581";

#[derive(Debug, Clone, Default, PartialEq)]
pub struct Ms2MetaSample {
    pub n_ms2: usize,
    pub collision_energies: Vec<f64>,
    pub isolation_widths: Vec<f64>,
    pub faims_cvs: Vec<f32>,
}

fn median(values: &[f64]) -> Option<f64> {
    let mut v: Vec<f64> = values.iter().copied().filter(|x| x.is_finite()).collect();
    if v.is_empty() {
        return None;
    }
    v.sort_by(f64::total_cmp);
    let n = v.len();
    Some(if n % 2 == 1 {
        v[n / 2]
    } else {
        0.5 * (v[n / 2 - 1] + v[n / 2])
    })
}

impl Ms2MetaSample {
    pub fn median_collision_energy(&self) -> Option<f64> {
        median(&self.collision_energies)
    }

    pub fn collision_energy_range(&self) -> Option<(f64, f64)> {
        let lo = self
            .collision_energies
            .iter()
            .copied()
            .filter(|x| x.is_finite())
            .reduce(f64::min)?;
        let hi = self
            .collision_energies
            .iter()
            .copied()
            .filter(|x| x.is_finite())
            .reduce(f64::max)?;
        Some((lo, hi))
    }

    pub fn median_isolation_width(&self) -> Option<f64> {
        median(&self.isolation_widths)
    }

    pub fn distinct_faims_cvs(&self) -> Vec<f32> {
        let mut v: Vec<f32> = self
            .faims_cvs
            .iter()
            .copied()
            .filter(|x| x.is_finite())
            .collect();
        v.sort_by(f32::total_cmp);
        v.dedup();
        v
    }
}

pub fn sample_mzml_ms2_meta<R: BufRead>(reader: R, n: usize) -> Ms2MetaSample {
    let mut xml = Reader::from_reader(reader);
    xml.trim_text(true);
    let mut buf: Vec<u8> = Vec::with_capacity(4096);
    let mut out = Ms2MetaSample::default();

    let mut in_spectrum = false;
    let mut in_activation = false;
    let mut in_product = false;
    let mut ms_level: Option<u32> = None;
    let mut ce: Option<f64> = None;
    let mut lower: Option<f64> = None;
    let mut upper: Option<f64> = None;
    let mut faims: Option<f32> = None;

    while out.n_ms2 < n {
        buf.clear();
        let event = match xml.read_event_into(&mut buf) {
            Ok(e) => e,
            Err(_) => break,
        };
        match event {
            Event::Eof => break,
            Event::Start(ref e) | Event::Empty(ref e) => {
                let is_start = matches!(event, Event::Start(_));
                match e.local_name().as_ref() {
                    b"spectrum" if is_start => {
                        in_spectrum = true;
                        in_activation = false;
                        in_product = false;
                        ms_level = None;
                        ce = None;
                        lower = None;
                        upper = None;
                        faims = None;
                    }
                    b"activation" if is_start && in_spectrum => in_activation = true,
                    b"product" if is_start && in_spectrum => in_product = true,
                    b"cvParam" | b"userParam" if in_spectrum => {
                        let Some(acc) = attr(e, b"accession") else {
                            continue;
                        };
                        let value = attr(e, b"value").and_then(|v| v.trim().parse::<f64>().ok());
                        match acc.as_str() {
                            CV_MS_LEVEL => ms_level = value.map(|v| v as u32),
                            CV_COLLISION_ENERGY | CV_NORMALIZED_COLLISION_ENERGY
                                if in_activation && ce.is_none() =>
                            {
                                ce = value.filter(|v| v.is_finite() && *v > 0.0);
                            }
                            CV_ISOLATION_LOWER_OFFSET if !in_product && lower.is_none() => {
                                lower = value;
                            }
                            CV_ISOLATION_UPPER_OFFSET if !in_product && upper.is_none() => {
                                upper = value;
                            }
                            CV_FAIMS_CV if faims.is_none() => faims = value.map(|v| v as f32),
                            _ => {}
                        }
                    }
                    _ => {}
                }
            }
            Event::End(ref e) => match e.local_name().as_ref() {
                b"activation" => in_activation = false,
                b"product" => in_product = false,
                b"spectrum" if in_spectrum => {
                    in_spectrum = false;
                    if ms_level == Some(2) {
                        out.n_ms2 += 1;
                        if let Some(v) = ce {
                            out.collision_energies.push(v);
                        }
                        if let (Some(lo), Some(up)) = (lower, upper) {
                            let w = lo.abs() + up.abs();
                            if w > 0.0 {
                                out.isolation_widths.push(w);
                            }
                        }
                        if let Some(cv) = faims {
                            out.faims_cvs.push(cv);
                        }
                    }
                }
                _ => {}
            },
            _ => {}
        }
    }
    out
}

const ISOTOPE_SPACING_DA: f64 = model::mass::ISOTOPE;
const PAIR_TOL_PPM: f64 = 10.0;
const CONTROL_OFFSETS_DA: [f64; 6] = [-0.21, -0.14, -0.07, 0.07, 0.14, 0.21];
const PAIR_PEAKS_PER_SPECTRUM: usize = 200;
pub const HIGH_RES_MIN_ENRICHMENT: f64 = 4.0;
pub const LOW_RES_MAX_ENRICHMENT: f64 = 1.8;
const MIN_PAIRS_HIGH: usize = 100;
const MIN_PEAKS: usize = 5_000;

#[derive(Debug, Clone, PartialEq)]
pub struct ResolutionEvidence {
    pub spectra: usize,
    pub peaks: usize,
    pub isotope_pairs: usize,
    pub chance_pairs: f64,
    pub enrichment: f64,
    pub high_res: Option<bool>,
}

impl ResolutionEvidence {
    pub fn describe(&self) -> String {
        format!(
            "{} spectra, {} peaks: {} isotope pairs at 1.00335/z within {PAIR_TOL_PPM} ppm vs {:.1} \
             by chance (enrichment {:.1}x; high-res at >= {HIGH_RES_MIN_ENRICHMENT}x, low-res below \
             {LOW_RES_MAX_ENRICHMENT}x)",
            self.spectra, self.peaks, self.isotope_pairs, self.chance_pairs, self.enrichment
        )
    }
}

fn count_partners(mzs: &[f64], delta: f64) -> usize {
    let mut n = 0;
    for &mz in mzs {
        let target = mz + delta;
        let tol = target * PAIR_TOL_PPM * 1e-6;
        let i = mzs.partition_point(|&x| x < target - tol);
        if i < mzs.len() && mzs[i] <= target + tol {
            n += 1;
        }
    }
    n
}

pub fn infer_fragment_resolution(spectra: &[model::Spectrum]) -> ResolutionEvidence {
    let mut peaks = 0usize;
    let mut pairs = 0usize;
    let mut control = 0usize;
    let mut used = 0usize;
    let mut mzs: Vec<f64> = Vec::new();
    for s in spectra {
        if s.peaks.len() < 2 {
            continue;
        }
        used += 1;
        let mut by_int: Vec<&(f64, f32)> = s.peaks.iter().collect();
        if by_int.len() > PAIR_PEAKS_PER_SPECTRUM {
            by_int.select_nth_unstable_by(PAIR_PEAKS_PER_SPECTRUM, |a, b| b.1.total_cmp(&a.1));
            by_int.truncate(PAIR_PEAKS_PER_SPECTRUM);
        }
        mzs.clear();
        mzs.extend(
            by_int
                .iter()
                .map(|p| p.0)
                .filter(|m| m.is_finite() && *m > 0.0),
        );
        mzs.sort_by(f64::total_cmp);
        peaks += mzs.len();
        for z in [1.0, 2.0] {
            let delta = ISOTOPE_SPACING_DA / z;
            pairs += count_partners(&mzs, delta);
            for off in CONTROL_OFFSETS_DA {
                control += count_partners(&mzs, delta + off);
            }
        }
    }
    let chance = control as f64 / CONTROL_OFFSETS_DA.len() as f64;
    let enrichment = pairs as f64 / chance.max(1.0);
    let high_res = if peaks < MIN_PEAKS {
        None
    } else if enrichment >= HIGH_RES_MIN_ENRICHMENT && pairs >= MIN_PAIRS_HIGH {
        Some(true)
    } else if enrichment < LOW_RES_MAX_ENRICHMENT {
        Some(false)
    } else {
        None
    };
    ResolutionEvidence {
        spectra: used,
        peaks,
        isotope_pairs: pairs,
        chance_pairs: chance,
        enrichment,
        high_res,
    }
}

fn attr(e: &quick_xml::events::BytesStart<'_>, name: &[u8]) -> Option<String> {
    e.attributes()
        .filter_map(|a| a.ok())
        .find(|a| a.key.local_name().as_ref() == name)
        .and_then(|a| {
            std::str::from_utf8(a.value.as_ref())
                .ok()
                .map(str::to_owned)
        })
}

#[cfg(test)]
mod tests {
    use super::*;

    fn spectrum(level: u32, ce: Option<f64>, lo: f64, up: f64, faims: Option<f64>) -> String {
        let ce_line = ce
            .map(|v| {
                format!(r#"<cvParam cvRef="MS" accession="MS:1000045" name="collision energy" value="{v}"/>"#)
            })
            .unwrap_or_default();
        let faims_line = faims
            .map(|v| {
                format!(r#"<cvParam cvRef="MS" accession="MS:1001581" name="FAIMS compensation voltage" value="{v}"/>"#)
            })
            .unwrap_or_default();
        format!(
            r#"<spectrum index="0" id="scan=1" defaultArrayLength="0">
<cvParam cvRef="MS" accession="MS:1000511" name="ms level" value="{level}"/>
<scanList count="1"><scan>{faims_line}</scan></scanList>
<precursorList count="1"><precursor>
<isolationWindow>
<cvParam cvRef="MS" accession="MS:1000827" name="isolation window target m/z" value="500.0"/>
<cvParam cvRef="MS" accession="MS:1000828" name="isolation window lower offset" value="{lo}"/>
<cvParam cvRef="MS" accession="MS:1000829" name="isolation window upper offset" value="{up}"/>
</isolationWindow>
<activation>
<cvParam cvRef="MS" accession="MS:1000422" name="beam-type collision-induced dissociation" value=""/>
{ce_line}
</activation>
</precursor></precursorList>
<productList count="1"><product><isolationWindow>
<cvParam cvRef="MS" accession="MS:1000828" name="isolation window lower offset" value="9.0"/>
<cvParam cvRef="MS" accession="MS:1000829" name="isolation window upper offset" value="9.0"/>
</isolationWindow></product></productList>
</spectrum>"#
        )
    }

    fn mzml(spectra: &[String]) -> String {
        format!(
            "<?xml version=\"1.0\"?><mzML><run><spectrumList count=\"{}\">{}</spectrumList></run></mzML>",
            spectra.len(),
            spectra.join("\n")
        )
    }

    #[test]
    fn reads_collision_energy_isolation_and_faims_of_ms2_only() {
        let doc = mzml(&[
            spectrum(1, None, 0.0, 0.0, Some(-45.0)),
            spectrum(2, Some(28.0), 0.7, 0.7, Some(-45.0)),
            spectrum(2, Some(30.0), 0.7, 0.7, Some(-65.0)),
            spectrum(2, None, 0.8, 0.8, None),
        ]);
        let s = sample_mzml_ms2_meta(doc.as_bytes(), 100);
        assert_eq!(s.n_ms2, 3);
        assert_eq!(s.collision_energies, vec![28.0, 30.0]);
        assert_eq!(s.median_collision_energy(), Some(29.0));
        assert_eq!(s.collision_energy_range(), Some((28.0, 30.0)));
        let w = s.median_isolation_width().unwrap();
        assert!((w - 1.4).abs() < 1e-9, "{w}");
        assert_eq!(s.distinct_faims_cvs(), vec![-65.0, -45.0]);
    }

    #[test]
    fn stops_after_n_ms2_spectra() {
        let doc = mzml(&[
            spectrum(2, Some(27.0), 1.0, 1.0, None),
            spectrum(2, Some(35.0), 1.0, 1.0, None),
        ]);
        let s = sample_mzml_ms2_meta(doc.as_bytes(), 1);
        assert_eq!(s.n_ms2, 1);
        assert_eq!(s.collision_energies, vec![27.0]);
    }

    #[test]
    fn no_metadata_gives_an_empty_sample() {
        let s = sample_mzml_ms2_meta(mzml(&[]).as_bytes(), 10);
        assert_eq!(s, Ms2MetaSample::default());
        assert_eq!(s.median_collision_energy(), None);
        assert_eq!(s.median_isolation_width(), None);
        assert!(s.distinct_faims_cvs().is_empty());
    }

    struct Lcg(u64);
    impl Lcg {
        fn next(&mut self) -> f64 {
            self.0 = self
                .0
                .wrapping_mul(6_364_136_223_846_793_005)
                .wrapping_add(1_442_695_040_888_963_407);
            (self.0 >> 11) as f64 / (1u64 << 53) as f64
        }
        fn normal(&mut self) -> f64 {
            let (u1, u2) = (self.next().max(1e-12), self.next());
            (-2.0 * u1.ln()).sqrt() * (2.0 * std::f64::consts::PI * u2).cos()
        }
    }

    fn synthetic(n_spectra: usize, mz_sd_da: f64, seed: u64) -> Vec<model::Spectrum> {
        let mut r = Lcg(seed);
        (0..n_spectra)
            .map(|_| {
                let mut peaks: Vec<(f64, f32)> = Vec::new();
                for _ in 0..30 {
                    let mono = 150.0 + 1200.0 * r.next();
                    let z = if r.next() < 0.8 { 1.0 } else { 2.0 };
                    let i = (1000.0 * r.next()) as f32 + 10.0;
                    peaks.push((mono + mz_sd_da * r.normal(), i));
                    peaks.push((
                        mono + ISOTOPE_SPACING_DA / z + mz_sd_da * r.normal(),
                        0.6 * i,
                    ));
                }
                for _ in 0..60 {
                    peaks.push((150.0 + 1200.0 * r.next(), (300.0 * r.next()) as f32 + 1.0));
                }
                peaks.sort_by(|a, b| a.0.total_cmp(&b.0));
                model::Spectrum {
                    peaks,
                    ..Default::default()
                }
            })
            .collect()
    }

    #[test]
    fn orbitrap_like_peaks_are_called_high_resolution() {
        let ev = infer_fragment_resolution(&synthetic(300, 0.0015, 7));
        assert_eq!(ev.high_res, Some(true), "{}", ev.describe());
        assert!(ev.enrichment > 10.0, "{}", ev.describe());
    }

    #[test]
    fn ion_trap_like_peaks_are_called_low_resolution() {
        let ev = infer_fragment_resolution(&synthetic(300, 0.1, 11));
        assert_eq!(ev.high_res, Some(false), "{}", ev.describe());
    }

    #[test]
    fn a_small_sample_makes_no_call() {
        let ev = infer_fragment_resolution(&synthetic(10, 0.0015, 3));
        assert_eq!(ev.high_res, None, "{}", ev.describe());
        let ev = infer_fragment_resolution(&[]);
        assert_eq!(ev.high_res, None);
        assert_eq!(ev.peaks, 0);
    }
}
