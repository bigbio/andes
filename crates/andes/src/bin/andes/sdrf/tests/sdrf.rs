use super::*;

fn sdrf(header: &[&str], rows: &[&[&str]]) -> String {
    let mut s = header.join("\t");
    for r in rows {
        s.push('\n');
        s.push_str(&r.join("\t"));
    }
    s
}

const H: &[&str] = &[
    "source name",
    "Comment[data file]",
    "comment[label]",
    "comment[cleavage agent details]",
    "comment[modification parameters]",
    "comment[modification parameters]",
    "comment[precursor mass tolerance]",
    "comment[fragment mass tolerance]",
    "comment[dissociation method]",
];

#[test]
fn label_free_trypsin_phospho() {
    let t = sdrf(
        H,
        &[
            &[
                "s1",
                "a.raw",
                "label free sample",
                "NT=Trypsin;AC=MS:1001251",
                "NT=Carbamidomethyl;AC=UNIMOD:4;TA=C;MT=Fixed;PP=Anywhere",
                "NT=Phospho;TA=S,T,Y;MT=Variable;PP=Anywhere",
                "10 ppm",
                "0.02 Da",
                "NT=HCD;AC=PRIDE:0000590",
            ],
            &[
                "s2",
                "b.raw",
                "label free sample",
                "NT=Trypsin;AC=MS:1001251",
                "NT=Carbamidomethyl;AC=UNIMOD:4;TA=C;MT=Fixed;PP=Anywhere",
                "NT=Phospho;TA=S,T,Y;MT=Variable;PP=Anywhere",
                "10 ppm",
                "0.02 Da",
                "NT=HCD;AC=PRIDE:0000590",
            ],
        ],
    );
    let s = parse(&t, Path::new("/d")).unwrap();
    assert_eq!(
        s.data_files,
        vec![PathBuf::from("/d/a.raw"), PathBuf::from("/d/b.raw")]
    );
    assert_eq!(s.enzyme.as_deref(), Some("trypsin"));
    assert_eq!(s.precursor_tol.as_deref(), Some("10ppm"));
    assert_eq!(s.fixed_mods, Some(vec!["Carbamidomethyl (C)".to_string()]));
    assert_eq!(s.variable_mods, Some(vec!["Phospho (STY)".to_string()]));
    assert_eq!(s.quant, Some(Quant::Lfq));
    assert_eq!(s.dissociation, Some(ActivationMethod::HCD));
    assert_eq!(s.fragment_tol.as_deref(), Some("0.02 Da"));
    model::mod_spec::parse_mod_list(&s.fixed_mods.unwrap().join(", "), true).unwrap();
    model::mod_spec::parse_mod_list(&s.variable_mods.unwrap().join(", "), false).unwrap();
}

#[test]
fn tmt_plex_is_the_smallest_holding_the_labels() {
    let row = |f: &str, l: &str| -> Vec<String> {
        vec![
            "s".into(),
            f.into(),
            l.into(),
            "NT=Trypsin".into(),
            "NT=TMT6plex;TA=K;MT=Fixed;PP=Anywhere".into(),
            "NT=TMT6plex;PP=Any N-term;MT=Fixed".into(),
            "20 ppm".into(),
            "not available".into(),
            "HCD".into(),
        ]
    };
    let labels10 = [
        "TMT126", "TMT127N", "TMT127C", "TMT128N", "TMT128C", "TMT129N", "TMT129C", "TMT130N",
        "TMT130C", "TMT131",
    ];
    let rows: Vec<Vec<String>> = labels10.iter().map(|l| row("a.raw", l)).collect();
    let refs: Vec<Vec<&str>> = rows
        .iter()
        .map(|r| r.iter().map(String::as_str).collect())
        .collect();
    let slices: Vec<&[&str]> = refs.iter().map(|r| r.as_slice()).collect();
    let s = parse(&sdrf(H, &slices), Path::new(".")).unwrap();
    assert_eq!(s.quant, Some(Quant::Isobaric("tmt10".into())));
    assert_eq!(s.data_files.len(), 1);
    assert_eq!(
        s.fixed_mods,
        Some(vec![
            "TMT6plex (K)".to_string(),
            "TMT6plex (N-term)".to_string()
        ])
    );
    let mut rows11 = rows.clone();
    rows11.push(row("a.raw", "TMT131C"));
    let refs: Vec<Vec<&str>> = rows11
        .iter()
        .map(|r| r.iter().map(String::as_str).collect())
        .collect();
    let slices: Vec<&[&str]> = refs.iter().map(|r| r.as_slice()).collect();
    let s = parse(&sdrf(H, &slices), Path::new(".")).unwrap();
    assert_eq!(s.quant, Some(Quant::Isobaric("tmt11".into())));
}

#[test]
fn silac_channels_come_from_the_label_mods_and_leave_the_search_lists() {
    let silac = |mods: [&str; 2]| {
        let rows = ["SILAC light", "SILAC heavy"]
            .map(|l| ["s", "a.raw", l, "NT=Trypsin", mods[0], mods[1], "10 ppm", "", ""]);
        let slices: Vec<&[&str]> = rows.iter().map(|r| r.as_slice()).collect();
        parse(&sdrf(H, &slices), Path::new(".")).unwrap()
    };
    let s = silac([
        "NT=Label:13C(6)15N(2);TA=K;MT=Variable;PP=Anywhere",
        "NT=Label:13C(6)15N(4);TA=R;MT=Variable;PP=Anywhere",
    ]);
    assert_eq!(s.quant, Some(Quant::Label("[][Lys8,Arg10]".into())));
    assert_eq!(s.variable_mods, Some(vec![]));
    let s = silac(["NT=Label:13C(6);TA=K,R;MT=Variable;PP=Anywhere", ""]);
    assert_eq!(s.quant, Some(Quant::Label("[][Lys6,Arg6]".into())));
    assert_eq!(s.variable_mods, Some(vec![]));
}

#[test]
fn rows_that_disagree_are_an_error_naming_the_column_and_files() {
    let a: &[&str] = &[
        "s",
        "a.raw",
        "label free sample",
        "Trypsin",
        "",
        "",
        "10 ppm",
        "",
        "",
    ];
    let b: &[&str] = &[
        "s",
        "b.raw",
        "label free sample",
        "Lys-C",
        "",
        "",
        "10 ppm",
        "",
        "",
    ];
    let e = parse(&sdrf(H, &[a, b]), Path::new(".")).unwrap_err();
    assert!(
        e.contains("cleavage agent details") && e.contains("a.raw") && e.contains("b.raw"),
        "{e}"
    );
}

#[test]
fn the_name_or_the_accession_whichever_resolves() {
    let item = |v: &str| ModParam::parse(v).unwrap().item();
    assert_eq!(
        item("NT=Acetyl;AC=UNIMOD:67;PP=Protein N-term;MT=Variable").unwrap(),
        "Acetyl (Protein N-term)"
    );
    assert_eq!(
        item("NT=Phosphorylation;AC=UNIMOD:21;TA=S,T;MT=Variable").unwrap(),
        "UNIMOD:21 (ST)"
    );
    assert!(item("NT=Nothing;TA=S;MT=Variable").is_err());
}

#[test]
fn a_converted_data_file_is_found_by_its_stem() {
    let dir = tempfile::tempdir().unwrap();
    std::fs::write(dir.path().join("run1.mzML"), "").unwrap();
    assert_eq!(
        resolve_data_file(dir.path(), "run1.raw"),
        dir.path().join("run1.mzML")
    );
    assert_eq!(
        resolve_data_file(dir.path(), "run2.raw"),
        dir.path().join("run2.raw")
    );
}

#[test]
fn terminal_mods_and_reserved_words() {
    assert_eq!(
        ModParam::parse("NT=Acetyl;AC=UNIMOD:1;PP=Protein N-term;MT=Variable")
            .unwrap()
            .item()
            .unwrap(),
        "Acetyl (Protein N-term)"
    );
    assert_eq!(
        ModParam::parse("NT=Gln->pyro-Glu;TA=Q;PP=Any N-term;MT=Variable")
            .unwrap()
            .item()
            .unwrap(),
        "Gln->pyro-Glu (N-term Q)"
    );
    assert!(
        ModParam::parse("NT=Oxidation;TA=M").is_err(),
        "MT is required"
    );
    assert_eq!(value(" not applicable "), None);
    assert_eq!(tolerance("20 mmu").unwrap(), "0.02da");
    assert_eq!(norm_header("Comment [ Data File ]"), "comment[data file]");
    assert!(enzyme_of("NT=Pepsin").is_err());
    let two = enzyme_of("NT=Trypsin/P;AC=MS:1001313;NT=Lys-C;AC=MS:1001309").unwrap_err();
    assert!(
        two.contains("NT= appears twice") && two.contains(CLEAVAGE),
        "{two}"
    );
    assert_eq!(enzyme_of("NT=Trypsin/P;AC=MS:1001313").unwrap(), "trypsin");
    for both in ["NT=Trypsin;Lys-C", "Lys-C;NT=Trypsin"] {
        assert!(enzyme_of(both).is_err(), "{both}");
    }
    assert_eq!(enzyme_of("Trypsin").unwrap(), "trypsin");
}

const HX: &[&str] = &[
    "source name",
    "comment[data file]",
    "comment[label]",
    "comment[cross-linker]",
];

#[test]
fn the_crosslinkers_are_read_checked_and_listed() {
    let row = |f: &'static str, x: &'static str| -> Vec<&'static str> {
        vec!["s", f, "label free sample", x]
    };
    let dsso = "NT=DSSO;AC=XLMOD:02126;CL=yes;TA=K,S,T,Y,nterm;MH=54.01;ML=85.98";
    let t = sdrf(HX, &[&row("a.raw", dsso), &row("b.raw", dsso)]);
    let s = parse(&t, Path::new(".")).unwrap();
    assert_eq!(s.xl, vec![dsso.to_string()]);
    let t = sdrf(
        HX,
        &[&row("a.raw", dsso), &row("b.raw", "NT=DSS;AC=XLMOD:02001")],
    );
    let s = parse(&t, Path::new(".")).unwrap();
    assert_eq!(
        s.xl,
        vec![dsso.to_string(), "NT=DSS;AC=XLMOD:02001".to_string()]
    );
    let t = sdrf(HX, &[&row("a.raw", dsso), &row("b.raw", "not applicable")]);
    let e = parse(&t, Path::new(".")).unwrap_err();
    assert!(
        e.contains("cross-linker") && e.contains("some files"),
        "{e}"
    );
    let t = sdrf(HX, &[&row("a.raw", "NT=DSS;CL=yes")]);
    let e = parse(&t, Path::new(".")).unwrap_err();
    assert!(e.contains("CL=yes but DSS is not MS-cleavable"), "{e}");
    let t = sdrf(HX, &[&row("a.raw", "not applicable")]);
    assert!(parse(&t, Path::new(".")).unwrap().xl.is_empty());
}

#[test]
fn the_sdrf_crosslinker_applies_below_the_command_line_and_the_config() {
    use clap::{CommandFactory, FromArgMatches};
    let s = SdrfSettings {
        xl: vec!["NT=DSBU;AC=XLMOD:02120;TA=K,nterm".into()],
        ..Default::default()
    };
    let run = |extra: &[&str], yaml: Option<&str>| -> crate::SearchArgs {
        let base = [
            "andes",
            "--spectrum",
            "x",
            "--database",
            "y",
            "--output-pin",
            "o",
        ];
        let m = crate::TopCli::command()
            .try_get_matches_from(base.iter().chain(extra))
            .unwrap();
        let mut a = crate::SearchArgs::from_arg_matches(&m).unwrap();
        let mut src = ParamSources::default();
        if let Some(y) = yaml {
            let c = crate::config::RunConfig::from_yaml_str(y).unwrap();
            crate::config::apply(c, &mut a, &m, &mut src).unwrap();
        }
        apply(&s, &mut a, &m, &mut src).unwrap();
        a
    };
    let x = run(&[], None).xl.unwrap().0.remove(0);
    assert_eq!(x.name, "DSBU");
    assert_eq!(x.describe_sites(), "K, Protein N-term");
    assert_eq!(
        run(&["--xl", "DSS"], None).xl.unwrap().0.remove(0).name,
        "DSS"
    );
    let cfg = Some("search:\n  xl: DSSO\n");
    assert_eq!(run(&[], cfg).xl.unwrap().0.remove(0).name, "DSSO");
    assert_eq!(
        run(&["--xl", "DSS"], cfg).xl.unwrap().0.remove(0).name,
        "DSS"
    );
}

#[test]
fn equal_sdrf_crosslinkers_written_differently_are_one_linker() {
    use clap::{CommandFactory, FromArgMatches};
    let apply_xl = |cells: &[&str]| -> Result<Vec<String>, String> {
        let s = SdrfSettings {
            xl: cells.iter().map(|c| c.to_string()).collect(),
            ..Default::default()
        };
        let base = [
            "andes",
            "--spectrum",
            "x",
            "--database",
            "y",
            "--output-pin",
            "o",
        ];
        let m = crate::TopCli::command().try_get_matches_from(base).unwrap();
        let mut a = crate::SearchArgs::from_arg_matches(&m).unwrap();
        apply(&s, &mut a, &m, &mut ParamSources::default())?;
        Ok(a.xl.unwrap().0.into_iter().map(|x| x.name).collect())
    };
    assert_eq!(
        apply_xl(&["NT=DSSO;AC=XLMOD:02126", "AC=XLMOD:02126;NT=DSSO"]).unwrap(),
        vec!["DSSO"]
    );
    assert_eq!(
        apply_xl(&["NT=DSSO;TA=K,S", "NT=DSSO;TA=S,K"]).unwrap(),
        vec!["DSSO"]
    );
    let e = apply_xl(&["NT=DSSO;TA=K", "NT=DSSO;TA=K,S"]).unwrap_err();
    assert!(e.contains("different settings"), "{e}");
}

#[test]
fn a_spectrum_missing_from_the_sdrf_is_a_line_of_its_own_marked_as_a_warning() {
    use clap::{CommandFactory, FromArgMatches};
    let s = SdrfSettings {
        data_files: vec![PathBuf::from("run1.mzML")],
        ..Default::default()
    };
    let base = [
        "andes",
        "--spectrum",
        "other.mzML",
        "--database",
        "y",
        "--output-pin",
        "o",
    ];
    let m = crate::TopCli::command().try_get_matches_from(base).unwrap();
    let mut a = crate::SearchArgs::from_arg_matches(&m).unwrap();
    let lines = apply(&s, &mut a, &m, &mut ParamSources::default()).unwrap();
    assert!(!lines[0].contains("WARN"), "{lines:?}");
    assert!(
        lines
            .iter()
            .any(|l| l.starts_with("WARN: --spectrum files not in the SDRF")
                && l.ends_with("other.mzML")),
        "{lines:?}"
    );
}

const DESIGN_H: &[&str] = &[
    "source name",
    "characteristics[organism part]",
    "characteristics[biological replicate]",
    "comment[instrument]",
    "comment[label]",
    "comment[fraction identifier]",
    "comment[technical replicate]",
    "comment[data file]",
];

fn design_row<'a>(
    source: &'a str,
    part: &'a str,
    bio: &'a str,
    label: &'a str,
    fraction: &'a str,
    file: &'a str,
) -> [&'a str; 8] {
    [
        source,
        part,
        bio,
        "NT=Orbitrap Fusion;AC=MS:1002416",
        label,
        fraction,
        "1",
        file,
    ]
}

#[test]
fn reference_and_blank_channels_get_their_roles_and_fractions_one_plex() {
    let rows = [
        design_row("S1", "Lung", "1", "TMT126", "1", "p1_Fr01.raw"),
        design_row("S2", "Liver", "1", "TMT127N", "1", "p1_Fr01.raw"),
        design_row(
            "PXD-GIS",
            "Global Internal Standard",
            "pooled",
            "TMT131",
            "1",
            "p1_Fr01.raw",
        ),
        design_row("PXD-blank", "blank", "1", "TMT127C", "1", "p1_Fr01.raw"),
        design_row("S1", "Lung", "1", "TMT126", "2", "p1_Fr02.raw"),
        design_row("S2", "Liver", "1", "TMT127N", "2", "p1_Fr02.raw"),
        design_row(
            "PXD-GIS",
            "Global Internal Standard",
            "pooled",
            "TMT131",
            "2",
            "p1_Fr02.raw",
        ),
        design_row("PXD-blank", "blank", "1", "TMT127C", "2", "p1_Fr02.raw"),
        design_row("S3", "Lung", "1", "TMT126", "1", "p2_Fr01.mzML.gz"),
    ];
    let rows: Vec<&[&str]> = rows.iter().map(|r| &r[..]).collect();
    let d = parse_isobaric_design(&sdrf(DESIGN_H, &rows)).unwrap();
    assert_eq!(d.plexes, 2);
    let r1 = &d.runs["p1_Fr01"];
    let r2 = &d.runs["p1_Fr02"];
    assert_eq!((r1.plex, r2.plex, d.runs["p2_Fr01"].plex), (1, 1, 2));
    assert_eq!(r2.fraction.as_deref(), Some("2"));
    assert_eq!(r1.instrument.as_deref(), Some("Orbitrap Fusion"));
    assert_eq!(r1.channels["126"].role, ChannelRole::Sample);
    assert_eq!(r1.channels["131"].role, ChannelRole::Reference);
    assert_eq!(r1.channels["131"].sample, "PXD-GIS");
    assert_eq!(r1.channels["127C"].role, ChannelRole::Blank);
}

#[test]
fn the_pooled_sample_column_and_an_empty_source_name_word_set_the_roles() {
    let h = [
        "source name",
        "characteristics[pooled sample]",
        "comment[label]",
        "comment[data file]",
    ];
    let rows: [&[&str]; 3] = [
        &["ref", "SN=a,b", "TMT126", "x.raw"],
        &["s", "not pooled", "TMT127N", "x.raw"],
        &["chan_empty", "not pooled", "TMT127C", "x.raw"],
    ];
    let d = parse_isobaric_design(&sdrf(&h, &rows)).unwrap();
    let c = &d.runs["x"].channels;
    assert_eq!(c["126"].role, ChannelRole::Reference);
    assert_eq!(c["127N"].role, ChannelRole::Sample);
    assert_eq!(c["127C"].role, ChannelRole::Blank);
    assert_eq!(d.runs["x"].fraction, None);
}

#[test]
fn two_files_of_one_channel_map_and_fraction_are_refused() {
    let rows = [
        design_row("S1", "Lung", "1", "TMT126", "1", "a.raw"),
        design_row("S1", "Lung", "1", "TMT126", "1", "b.raw"),
    ];
    let rows: Vec<&[&str]> = rows.iter().map(|r| &r[..]).collect();
    let e = parse_isobaric_design(&sdrf(DESIGN_H, &rows)).unwrap_err();
    assert!(e.contains("cannot be told apart"), "{e}");
    let rows = [design_row(
        "S1",
        "Lung",
        "1",
        "label free sample",
        "1",
        "a.raw",
    )];
    let rows: Vec<&[&str]> = rows.iter().map(|r| &r[..]).collect();
    let e = parse_isobaric_design(&sdrf(DESIGN_H, &rows)).unwrap_err();
    assert!(e.contains("not an isobaric"), "{e}");
}

#[test]
fn a_technical_replicate_per_channel_keeps_the_fractions_of_a_plex_together() {
    let h = [
        "source name",
        "comment[label]",
        "comment[technical replicate]",
        "comment[fraction identifier]",
        "comment[data file]",
    ];
    let rows: [&[&str]; 8] = [
        &["S1", "TMT126", "2", "1", "p_Fr01.raw"],
        &["S2", "TMT127N", "3", "1", "p_Fr01.raw"],
        &["S1", "TMT126", "2", "2", "p_Fr02.raw"],
        &["S2", "TMT127N", "3", "2", "p_Fr02.raw"],
        &["S1", "TMT126", "4", "1", "q_Fr01.raw"],
        &["S2", "TMT127N", "4", "1", "q_Fr01.raw"],
        &["S1", "TMT126", "2", "3", "p_Fr03.raw"],
        &["S2", "TMT127N", "3", "3", "p_Fr03.raw"],
    ];
    let d = parse_isobaric_design(&sdrf(&h, &rows)).unwrap();
    let plex = |r: &str| d.runs[r].plex;
    assert_eq!((plex("p_Fr01"), plex("p_Fr02"), plex("p_Fr03")), (1, 1, 1));
    assert_eq!(
        plex("q_Fr01"),
        2,
        "a repeat injection of the plex is its own plex"
    );
}
