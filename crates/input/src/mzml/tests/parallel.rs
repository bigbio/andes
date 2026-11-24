//! Tests of the parallel piece readers against the sequential and chunked readers.

use std::io::Cursor;

use model::Spectrum;

use super::*;
use crate::mzml::parallel::{read_ms1_linked_in_pieces, read_spectra_in_pieces};
use crate::mzml::{Ms1Link, MzMLReader};

fn spectrum_fingerprint(s: &Spectrum) -> String {
    format!(
        "{}|{:?}|{:?}|{:?}|{:?}|{:?}|{:?}|{:?}|{:?}|{:?}|{:?}",
        s.title,
        s.precursor_mz.to_bits(),
        s.precursor_charge,
        s.rt_seconds.map(f64::to_bits),
        s.scan,
        s.activation_method,
        s.isolation_lower_offset.map(f64::to_bits),
        s.isolation_upper_offset.map(f64::to_bits),
        s.isolation_target_mz.map(f64::to_bits),
        s.precursor_intensity.map(f32::to_bits),
        s.peaks
            .iter()
            .map(|p| (p.0.to_bits(), p.1.to_bits()))
            .collect::<Vec<_>>()
    )
}

#[test]
fn parallel_read_matches_the_sequential_reader() {
    use std::io::Read;
    let gz = std::path::Path::new(env!("CARGO_MANIFEST_DIR"))
        .join("../../test-fixtures/orbitrap_lumos_120.mzML.gz");
    let mut xml = Vec::new();
    flate2::read::GzDecoder::new(std::fs::File::open(&gz).expect("open fixture"))
        .read_to_end(&mut xml)
        .expect("decompress fixture");
    let path =
        std::env::temp_dir().join(format!("andes-parallel-mzml-{}.mzML", std::process::id()));
    std::fs::write(&path, &xml).expect("write fixture copy");
    for (min, max) in [(2, 2), (1, 2), (1, 1)] {
        let sequential: Vec<String> = MzMLReader::new(xml.as_slice())
            .with_ms_level_range(min, max)
            .map(|r| spectrum_fingerprint(&r.expect("parse")))
            .collect();
        if (min, max) == (2, 2) {
            assert!(!sequential.is_empty(), "the fixture carries MS2 spectra");
        }
        for piece in [1_000u64, 7_919, 65_536, 1 << 30] {
            for threads in [1, 3] {
                let (spectra, errors) =
                    read_spectra_in_pieces(&path, min, max, threads, piece).expect("read");
                assert!(errors.is_empty(), "{errors:?}");
                let parallel: Vec<String> = spectra.iter().map(spectrum_fingerprint).collect();
                assert_eq!(
                    parallel, sequential,
                    "levels {min}-{max}, piece {piece}, {threads} threads"
                );
            }
        }
    }
    let _ = std::fs::remove_file(&path);
}

#[test]
fn an_ms2_scan_links_to_the_latest_ms1_at_its_faims_voltage() {
    let spectrum = |n: u32, level: u32, cv: Option<f64>, mz: f64| {
        let cv = cv.map_or(String::new(), |v| {
            format!(
                r#"<cvParam accession="MS:1001581" name="FAIMS compensation voltage" value="{v}"/>"#
            )
        });
        let precursor = if level == 2 {
            r#"<precursorList count="1"><precursor><selectedIonList count="1"><selectedIon>
                    <cvParam accession="MS:1000744" name="selected ion m/z" value="500.5"/>
                    <cvParam accession="MS:1000041" name="charge state" value="2"/>
                  </selectedIon></selectedIonList></precursor></precursorList>"#
        } else {
            ""
        };
        format!(
            r#"<spectrum index="{n}" id="scan={n}" defaultArrayLength="1">
              <cvParam accession="MS:1000511" name="ms level" value="{level}"/>
              <scanList count="1"><scan>
                <cvParam accession="MS:1000016" name="scan start time" value="1.0"
                         unitAccession="UO:0000031" unitName="minute"/>
                {cv}
              </scan></scanList>
              {precursor}
              <binaryDataArrayList count="2">{mzb}{intb}</binaryDataArrayList>
            </spectrum>"#,
            mzb = bda_plain("MS:1000514", &encode_f64_b64(&[mz])),
            intb = bda_plain("MS:1000515", &encode_f64_b64(&[7.0])),
        )
    };
    // MS1 at -50 (m/z 100), MS1 at -60 (200), MS2 -50, MS2 -60, MS1 at -50 (300),
    // MS2 -60, MS2 without a voltage.
    let xml = wrap_spectra(
        &[
            spectrum(1, 1, Some(-50.0), 100.0),
            spectrum(2, 1, Some(-60.0), 200.0),
            spectrum(3, 2, Some(-50.0), 400.0),
            spectrum(4, 2, Some(-60.0), 400.0),
            spectrum(5, 1, Some(-50.0), 300.0),
            spectrum(6, 2, Some(-60.0), 400.0),
            spectrum(7, 2, None, 400.0),
        ]
        .concat(),
    );
    let want = vec![Some(100.0), Some(200.0), Some(200.0), Some(300.0)];
    let first_mz = |link: &Ms1Link| -> Vec<Option<f64>> {
        link.ms2_to_ms1
            .iter()
            .map(|l| l.map(|m| link.ms1_peaks[m][0].0))
            .collect()
    };
    let (_, link) = MzMLReader::new(Cursor::new(xml.clone()))
        .with_ms1_capture(true)
        .read_with_ms1()
        .expect("read");
    assert_eq!(first_mz(&link), want);
    for chunk in [1, 2, 3] {
        let mut got = Vec::new();
        MzMLReader::new(Cursor::new(xml.clone())).read_with_ms1_chunked(
            chunk,
            usize::MAX,
            |_, link| got.extend(first_mz(&link)),
        );
        assert_eq!(got, want, "chunk {chunk}");
    }
    let path = std::env::temp_dir().join(format!("andes-faims-link-{}.mzML", std::process::id()));
    std::fs::write(&path, &xml).expect("write");
    for piece in [300u64, 1 << 30] {
        let read = read_ms1_linked_in_pieces(&path, 2, false, false, piece)
            .expect("read")
            .expect("no parse error");
        let got: Vec<Option<f64>> = read
            .ms2_to_ms1
            .iter()
            .map(|l| l.map(|m| read.run_ms1[m].peaks[0].0))
            .collect();
        assert_eq!(got, want, "piece {piece}");
    }
    let _ = std::fs::remove_file(&path);
}

#[test]
fn ms1_linked_parallel_read_matches_the_chunked_reader() {
    use std::io::Read;
    let gz = std::path::Path::new(env!("CARGO_MANIFEST_DIR"))
        .join("../../test-fixtures/mbr_run1.mzML.gz");
    let mut xml = Vec::new();
    flate2::read::GzDecoder::new(std::fs::File::open(&gz).expect("open fixture"))
        .read_to_end(&mut xml)
        .expect("decompress fixture");
    let path = std::env::temp_dir().join(format!(
        "andes-parallel-ms1-mzml-{}.mzML",
        std::process::id()
    ));
    std::fs::write(&path, &xml).expect("write fixture copy");
    for (context, noise) in [(true, false), (false, true)] {
        let mut reader = MzMLReader::new(xml.as_slice())
            .with_ms_level_range(2, 2)
            .with_ms1_capture(true)
            .with_run_ms1_capture(true)
            .with_scan_context_capture(context)
            .with_noise_capture(noise);
        let mut spectra = Vec::new();
        let mut linked = Vec::new();
        let (errors, _) = reader.read_with_ms1_chunked(7, usize::MAX, |chunk, link| {
            for (s, l) in chunk.into_iter().zip(&link.ms2_to_ms1) {
                spectra.push(format!("{s:?}"));
                linked.push(l.map(|m| link.ms1_peaks[m].clone()));
            }
        });
        assert_eq!(errors, 0);
        let run_ms1 = reader.take_run_ms1();
        let contexts = reader.take_scan_contexts();
        assert!(
            run_ms1.len() > 1 && !spectra.is_empty(),
            "MS1 and MS2 scans"
        );
        for piece in [1_000u64, 7_919, 65_536, 1 << 30] {
            for threads in [1, 3] {
                let read = read_ms1_linked_in_pieces(&path, threads, context, noise, piece)
                    .expect("read")
                    .expect("no parse error");
                let what = format!("piece {piece}, {threads} threads");
                let got: Vec<String> = read.spectra.iter().map(|s| format!("{s:?}")).collect();
                assert_eq!(got, spectra, "{what}");
                assert_eq!(read.run_ms1, run_ms1, "{what}");
                assert_eq!(read.contexts, contexts, "{what}");
                let got_linked: Vec<_> = read
                    .ms2_to_ms1
                    .iter()
                    .map(|l| l.map(|m| read.run_ms1[m].peaks.clone()))
                    .collect();
                assert_eq!(got_linked, linked, "{what}");
            }
        }
    }
    let _ = std::fs::remove_file(&path);
}
