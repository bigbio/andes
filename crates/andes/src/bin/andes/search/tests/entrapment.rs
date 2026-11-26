use super::*;
use clap::Parser;

#[derive(Parser)]
struct Args {
    #[command(flatten)]
    search: Cli,
}

fn psm(id: &str, peptide: &str, proteins: &str) -> output::PercolatorPsm {
    output::PercolatorPsm {
        psm_id: id.into(),
        peptide: peptide.into(),
        proteins: proteins.into(),
        q_value: 0.001,
        pooled_q_value: None,
        pep: 0.001,
        score: Some(5.0),
    }
}

#[test]
fn precursor_keys_ignore_flanks_but_keep_modifications() {
    assert_eq!(peptidoform("K.PEPTIM+15.99492K.G"), "PEPTLM+15.99492K");
    assert_eq!(peptidoform("R.PEPTIM+15.99492K.A"), "PEPTLM+15.99492K");
    assert_ne!(
        peptidoform("K.PEPTIMK.A"),
        peptidoform("K.PEPTIM+15.99492K.A")
    );
}

#[test]
fn charge_reader_rejects_multiple_active_charge_columns() {
    let dir = tempfile::tempdir().unwrap();
    let pin = dir.path().join("bad.pin");
    std::fs::write(&pin, "SpecId\tcharge2\tcharge3\ns1\t1\t1\n").unwrap();
    let rows = HashMap::from([("s1".into(), psm("s1", "PEPTIDEK", "P1"))]);
    assert!(charges(&pin, &rows)
        .unwrap_err()
        .contains("unique positive charge"));
}

#[test]
fn dda_and_dia_reporting_write_versioned_entrapment_curves() {
    for precursor_pin in [false, true] {
        let dir = tempfile::tempdir().unwrap();
        let pin = dir.path().join("run.pin");
        std::fs::write(&pin, "SpecId\tcharge2\tcharge3\nt\t1\t0\ne\t0\t1\n").unwrap();
        let proteins = [("P1", "MAGLEVTDPLSK"), ("ENTRAP_P1", "MPLGAVDTELSK")]
            .into_iter()
            .map(|(accession, sequence)| model::protein::Protein {
                accession: accession.into(),
                sequence: sequence.as_bytes().to_vec(),
                description: String::new(),
                annotation: None,
            })
            .collect();
        let index =
            search::SearchIndex::from_prebuilt(model::protein::ProteinDb { proteins }, "XXX", None);
        let cli = Args::try_parse_from(["andes", "--add-entrapment", "--rescore"])
            .unwrap()
            .search;
        let rows = HashMap::from([
            ("t".into(), psm("t", "K.MAGLEVTDPLSK.A", "P1")),
            ("e".into(), psm("e", "K.MPLGAVDTELSK.A", "ENTRAP_P1")),
        ]);
        let base = dir.path().join("run.pin");
        let text = report(Inputs {
            cli: &cli,
            index: &index,
            params: &search::SearchParams::default_tryptic(
                model::amino_acid::set::AminoAcidSetBuilder::new_standard().build().unwrap(),
            ),
            rows: Some(&rows),
            groups: None,
            base: &base,
            pin: &pin,
            precursor_pin,
        })
        .unwrap()
        .unwrap();
        assert!(text.contains("targets 1  entrapment 1"));
        assert!(text.contains(if precursor_pin {
            "dia_precursor_rows"
        } else {
            "psm:"
        }));
        let file = std::fs::File::open(base.with_extension("entrapment.parquet")).unwrap();
        let reader =
            parquet::arrow::arrow_reader::ParquetRecordBatchReaderBuilder::try_new(file).unwrap();
        assert_eq!(reader.metadata().file_metadata().num_rows(), 6);
        assert!(reader
            .metadata()
            .file_metadata()
            .key_value_metadata()
            .unwrap()
            .iter()
            .any(|m| m.key == "andes.entrapment.schema"));
    }
}
