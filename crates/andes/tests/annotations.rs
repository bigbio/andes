//! End-to-end: an annotation file next to a FASTA database searches like the
//! same claims in a UniProt XML database.

use std::path::{Path, PathBuf};
use std::process::{Command, Output};
use std::sync::Arc;

use arrow::array::{Array, ArrayRef, Int64Array, ListArray, RecordBatch, StringArray, StructArray};
use arrow::buffer::OffsetBuffer;
use arrow::datatypes::{DataType, Field, Schema};
use parquet::arrow::ArrowWriter;
use parquet::file::properties::WriterProperties;
use parquet::format::KeyValue;

fn fixture(rel: &str) -> PathBuf {
    PathBuf::from(env!("CARGO_MANIFEST_DIR"))
        .join("../..")
        .join(rel)
        .canonicalize()
        .unwrap_or_else(|e| panic!("canonicalize {rel}: {e}"))
}

fn search(database: &str, pin: &Path, extra: &[&str]) -> Output {
    let out = Command::new(env!("CARGO_BIN_EXE_andes"))
        .arg("--spectrum")
        .arg(fixture("test-fixtures/test.mgf.gz"))
        .arg("--database")
        .arg(fixture(database))
        .arg("--output-pin")
        .arg(pin)
        .args(["--decoy-prefix", "XXX_"])
        .args(extra)
        .output()
        .expect("run andes");
    assert!(
        out.status.success(),
        "andes failed on {database}: {}",
        String::from_utf8_lossy(&out.stderr)
    );
    out
}

/// One substitution row: (id, reference, md5, start, end, ref, alt).
type Substitution<'a> = (&'a str, &'a str, &'a str, i64, i64, &'a str, &'a str);

/// An annotations file with one substitution per row.
fn write_substitutions(path: &Path, rows: &[Substitution]) {
    let col = |f: &dyn Fn(&Substitution) -> String| -> ArrayRef {
        Arc::new(StringArray::from(rows.iter().map(f).collect::<Vec<_>>()))
    };
    let st = StructArray::from(vec![
        (
            Arc::new(Field::new("start", DataType::Int64, false)),
            Arc::new(Int64Array::from(
                rows.iter().map(|r| r.3).collect::<Vec<_>>(),
            )) as ArrayRef,
        ),
        (
            Arc::new(Field::new("end", DataType::Int64, false)),
            Arc::new(Int64Array::from(
                rows.iter().map(|r| r.4).collect::<Vec<_>>(),
            )),
        ),
        (
            Arc::new(Field::new("ref", DataType::Utf8, false)),
            Arc::new(StringArray::from(
                rows.iter().map(|r| r.5).collect::<Vec<_>>(),
            )),
        ),
        (
            Arc::new(Field::new("alt", DataType::Utf8, false)),
            Arc::new(StringArray::from(
                rows.iter().map(|r| r.6).collect::<Vec<_>>(),
            )),
        ),
    ]);
    let item = Arc::new(Field::new("item", st.data_type().clone(), true));
    let edits: ArrayRef = Arc::new(ListArray::new(
        item,
        OffsetBuffer::from_lengths(vec![1; rows.len()]),
        Arc::new(st),
        None,
    ));
    let columns: Vec<(&str, ArrayRef)> = vec![
        ("annotation_id", col(&|r| r.0.to_string())),
        ("reference_accession", col(&|r| r.1.to_string())),
        ("reference_md5", col(&|r| r.2.to_string())),
        ("source_db", col(&|_| "custom".to_string())),
        ("edits", edits),
    ];
    let schema = Arc::new(Schema::new(
        columns
            .iter()
            .map(|(n, a)| Field::new(*n, a.data_type().clone(), true))
            .collect::<Vec<_>>(),
    ));
    let batch = RecordBatch::try_new(
        schema.clone(),
        columns.into_iter().map(|(_, a)| a).collect(),
    )
    .unwrap();
    let props = WriterProperties::builder()
        .set_key_value_metadata(Some(vec![KeyValue::new(
            "andes.annotations.schema".to_string(),
            "1".to_string(),
        )]))
        .build();
    let mut w =
        ArrowWriter::try_new(std::fs::File::create(path).unwrap(), schema, Some(props)).unwrap();
    w.write(&batch).unwrap();
    w.close().unwrap();
}

#[test]
fn an_annotated_fasta_searches_like_the_uniprot_xml() {
    let dir = tempfile::tempdir().unwrap();
    // BSA's one curated variant, A214T, written as an annotation of the FASTA.
    let fasta = input::FastaReader::load_all(std::io::BufReader::new(
        std::fs::File::open(fixture("test-fixtures/BSA.fasta")).unwrap(),
    ))
    .unwrap();
    let bsa = &fasta.proteins[0];
    assert_eq!(bsa.sequence[213], b'A');
    let md5 = input::sequence_md5(&bsa.sequence);
    let ann = dir.path().join("bsa.annotations.parquet");
    write_substitutions(
        &ann,
        &[("VAR_A214T", &bsa.accession, &md5, 213, 214, "A", "T")],
    );

    let xml_pin = dir.path().join("xml.pin");
    search(
        "test-fixtures/BSA.xml.gz",
        &xml_pin,
        &["--uniprot-variants", "--uniprot-evidence", "curated"],
    );
    let fasta_pin = dir.path().join("fasta.pin");
    let out = search(
        "test-fixtures/BSA.fasta",
        &fasta_pin,
        &["--annotations", ann.to_str().unwrap()],
    );
    let stderr = String::from_utf8_lossy(&out.stderr);
    assert!(
        stderr.contains("Annotations: 1 rows from 1 files: 1 variant, "),
        "{stderr}"
    );
    let (xml, fasta) = (
        std::fs::read_to_string(&xml_pin).unwrap(),
        std::fs::read_to_string(&fasta_pin).unwrap(),
    );
    assert!(xml.contains("P02769_VAR_A214T"), "the variant is searched");
    assert_eq!(xml, fasta, "same entries, same PIN");
}

#[test]
fn rows_that_do_not_fit_are_reported() {
    let dir = tempfile::tempdir().unwrap();
    let ann = dir.path().join("bad.annotations.parquet");
    write_substitutions(
        &ann,
        &[(
            "X1",
            "sp|P02769|ALBU_BOVIN",
            &"0".repeat(32),
            213,
            214,
            "A",
            "T",
        )],
    );
    let out = search(
        "test-fixtures/BSA.fasta",
        &dir.path().join("x.pin"),
        &["--annotations", ann.to_str().unwrap()],
    );
    let stderr = String::from_utf8_lossy(&out.stderr);
    assert!(
        stderr.contains("left out: reference sequence differs (md5) (1)"),
        "{stderr}"
    );
}

/// A PIN's rows with the UniProt accession spelling `sp|P02769<x>|ALBU_BOVIN`
/// turned into the PEFF one `sp:P02769<x>`, and each row's proteins sorted.
fn normalised_rows(pin: &str) -> Vec<String> {
    let header: Vec<&str> = pin.lines().next().unwrap().split('\t').collect();
    let proteins = header.iter().position(|h| *h == "Proteins").unwrap();
    let mut rows: Vec<String> = pin
        .lines()
        .skip(1)
        .map(|l| {
            let cols: Vec<&str> = l.split('\t').collect();
            let mut prots: Vec<String> = cols[proteins..]
                .iter()
                .map(|p| p.replace("sp|", "sp:").replace("|ALBU_BOVIN", ""))
                .collect();
            prots.sort();
            format!("{}\t{}", cols[..proteins].join("\t"), prots.join("\t"))
        })
        .collect();
    rows.sort();
    rows
}

#[test]
fn a_peff_database_searches_like_the_uniprot_xml() {
    let dir = tempfile::tempdir().unwrap();
    let fasta = input::FastaReader::load_all(std::io::BufReader::new(
        std::fs::File::open(fixture("test-fixtures/BSA.fasta")).unwrap(),
    ))
    .unwrap();
    let seq = String::from_utf8(fasta.proteins[0].sequence.clone()).unwrap();
    let peff = dir.path().join("bsa.peff");
    std::fs::write(
        &peff,
        format!(
            "# PEFF 1.0\n# //\n>sp:P02769 \\PName=Serum albumin \\VariantSimple=(214|T)\n{seq}\n"
        ),
    )
    .unwrap();
    let xml_pin = dir.path().join("xml.pin");
    search(
        "test-fixtures/BSA.xml.gz",
        &xml_pin,
        &["--uniprot-variants", "--uniprot-evidence", "curated"],
    );
    let peff_pin = dir.path().join("peff.pin");
    let out = Command::new(env!("CARGO_BIN_EXE_andes"))
        .arg("--spectrum")
        .arg(fixture("test-fixtures/test.mgf.gz"))
        .arg("--database")
        .arg(&peff)
        .arg("--output-pin")
        .arg(&peff_pin)
        .args(["--decoy-prefix", "XXX_"])
        .output()
        .expect("run andes");
    let stderr = String::from_utf8_lossy(&out.stderr);
    assert!(out.status.success(), "{stderr}");
    assert!(stderr.contains("PEFF: 1 rows: 1 variant, "), "{stderr}");
    let (xml, peff) = (
        std::fs::read_to_string(&xml_pin).unwrap(),
        std::fs::read_to_string(&peff_pin).unwrap(),
    );
    assert!(
        peff.contains("sp:P02769_VAR_A214T"),
        "the variant is searched"
    );
    assert_eq!(normalised_rows(&xml), normalised_rows(&peff));
}
