use std::path::{Path, PathBuf};
use std::process::{Command, Output};

const SPECTRUM: &str = "test-fixtures/test.mgf.gz";
const DATABASE: &str = "test-fixtures/BSA.fasta";

fn fixture(rel: &str) -> PathBuf {
    PathBuf::from(env!("CARGO_MANIFEST_DIR"))
        .join("../..")
        .join(rel)
        .canonicalize()
        .unwrap_or_else(|e| panic!("canonicalize {rel}: {e}"))
}

fn andes(args: &[&std::ffi::OsStr]) -> Output {
    Command::new(env!("CARGO_BIN_EXE_andes"))
        .args(args)
        .output()
        .expect("run andes")
}

fn os(s: &str) -> &std::ffi::OsStr {
    std::ffi::OsStr::new(s)
}

fn printed(args: &[&std::ffi::OsStr]) -> String {
    let mut all: Vec<&std::ffi::OsStr> = args.to_vec();
    all.push(os("--print-config"));
    all.push(os("-"));
    let out = andes(&all);
    assert!(
        out.status.success(),
        "print-config failed: {}",
        String::from_utf8_lossy(&out.stderr)
    );
    String::from_utf8(out.stdout).expect("utf-8 config")
}

fn parse(text: &str) -> serde_yaml::Value {
    serde_yaml::from_str(text).unwrap_or_else(|e| panic!("printed config is not YAML: {e}\n{text}"))
}

fn key_line<'a>(text: &'a str, key: &str) -> Option<&'a str> {
    let prefix = format!("{key}:");
    text.lines().find(|l| l.trim_start().starts_with(&prefix))
}

const NON_DEFAULT: [&str; 6] = [
    "--precursor-tol",
    "10ppm",
    "--max-missed-cleavages",
    "1",
    "--fixed-mods",
    "Carbamidomethyl (C)",
];

#[test]
fn printed_config_reproduces_the_run() {
    let dir = tempfile::tempdir().expect("tempdir");
    let pin_direct = dir.path().join("direct.pin");
    let pin_replay = dir.path().join("replay.pin");
    let printed_path = dir.path().join("printed.yaml");
    let spectrum = fixture(SPECTRUM);
    let database = fixture(DATABASE);

    let common = |pin: &Path| -> Vec<std::ffi::OsString> {
        let mut v: Vec<std::ffi::OsString> = vec![
            "--spectrum".into(),
            spectrum.clone().into(),
            "--database".into(),
            database.clone().into(),
            "--output-pin".into(),
            pin.into(),
        ];
        v.extend(NON_DEFAULT.iter().map(std::ffi::OsString::from));
        v
    };

    let direct = common(&pin_direct);
    let direct_refs: Vec<&std::ffi::OsStr> = direct.iter().map(|s| s.as_os_str()).collect();
    let out = andes(&direct_refs);
    assert!(out.status.success(), "{}", String::from_utf8_lossy(&out.stderr));

    let mut print = common(&pin_replay);
    print.push("--print-config".into());
    print.push(printed_path.clone().into());
    let print_refs: Vec<&std::ffi::OsStr> = print.iter().map(|s| s.as_os_str()).collect();
    let out = andes(&print_refs);
    assert!(out.status.success(), "{}", String::from_utf8_lossy(&out.stderr));
    assert!(!pin_replay.exists(), "--print-config must not search");
    assert!(printed_path.exists(), "config not written");

    let out = andes(&[os("--config"), printed_path.as_os_str()]);
    assert!(out.status.success(), "{}", String::from_utf8_lossy(&out.stderr));

    let a = std::fs::read(&pin_direct).expect("direct pin");
    let b = std::fs::read(&pin_replay).expect("replayed pin");
    assert!(a.len() > 100, "PIN looks empty");
    assert!(a == b, "the replayed PIN differs from the direct run's PIN");
}

#[test]
fn print_config_round_trips_lists_and_mods() {
    let dir = tempfile::tempdir().expect("tempdir");
    let first = fixture(SPECTRUM);
    let second = dir.path().join("second.mgf.gz");
    std::fs::copy(&first, &second).expect("copy spectrum");
    let database = fixture(DATABASE);
    let pin = dir.path().join("out.pin");
    let args = [
        os("--spectrum"),
        first.as_os_str(),
        os("--spectrum"),
        second.as_os_str(),
        os("--database"),
        database.as_os_str(),
        os("--output-pin"),
        pin.as_os_str(),
        os("--fixed-mods"),
        os("Carbamidomethyl (C), TMT6plex (K)"),
        os("--variable-mods"),
        os("Oxidation (M), Phospho (STY)"),
    ];
    let text = printed(&args);
    let doc = parse(&text);

    let spectra = doc["io"]["spectrum"].as_sequence().expect("spectrum list");
    let got: Vec<&str> = spectra.iter().filter_map(|v| v.as_str()).collect();
    assert_eq!(
        got,
        vec![first.to_str().expect("utf-8"), second.to_str().expect("utf-8")]
    );
    for (key, name) in [("fixed_mods", "Carbamidomethyl"), ("variable_mods", "Oxidation")] {
        let items = doc["search"][key].as_sequence().expect(key);
        assert_eq!(items.len(), 2, "{key}: {text}");
        let first_item = serde_yaml::to_string(&items[0]).expect("yaml");
        assert!(first_item.contains(name), "{key}[0] = {first_item}");
    }

    let cfg = dir.path().join("printed.yaml");
    std::fs::write(&cfg, &text).expect("write config");
    let again = parse(&printed(&[os("--config"), cfg.as_os_str()]));
    assert_eq!(doc["io"], again["io"]);
    assert_eq!(doc["search"], again["search"]);
}

#[test]
fn print_config_quotes_paths() {
    let dir = tempfile::tempdir().expect("tempdir");
    let odd = dir.path().join("my db: v1.fasta");
    std::fs::copy(fixture(DATABASE), &odd).expect("copy database");
    let spectrum = fixture(SPECTRUM);
    let pin = dir.path().join("out.pin");
    let text = printed(&[
        os("--spectrum"),
        spectrum.as_os_str(),
        os("--database"),
        odd.as_os_str(),
        os("--output-pin"),
        pin.as_os_str(),
    ]);
    let doc = parse(&text);
    assert_eq!(doc["io"]["database"].as_str(), odd.to_str());

    let cfg = dir.path().join("printed.yaml");
    std::fs::write(&cfg, &text).expect("write config");
    let again = parse(&printed(&[os("--config"), cfg.as_os_str()]));
    assert_eq!(again["io"]["database"].as_str(), odd.to_str());
}

#[test]
fn print_config_omits_default_internal_flags() {
    let dir = tempfile::tempdir().expect("tempdir");
    let spectrum = fixture(SPECTRUM);
    let database = fixture(DATABASE);
    let pin = dir.path().join("out.pin");
    let base = [
        os("--spectrum"),
        spectrum.as_os_str(),
        os("--database"),
        database.as_os_str(),
        os("--output-pin"),
        pin.as_os_str(),
    ];
    let text = printed(&base);
    for key in ["optimize_sample", "max_spectra", "seed"] {
        assert!(key_line(&text, key).is_none(), "{key} printed:\n{text}");
    }
    assert!(!text.contains("print_config") && !text.contains("config:"));

    let mut set = base.to_vec();
    set.push(os("--max-spectra"));
    set.push(os("100"));
    let text = printed(&set);
    let line = key_line(&text, "max_spectra").expect("max_spectra printed once set");
    assert!(!line.trim_start().starts_with('#'), "{line}");
}

#[test]
fn print_config_names_sources() {
    let dir = tempfile::tempdir().expect("tempdir");
    let cfg = dir.path().join("in.yaml");
    std::fs::write(&cfg, "search:\n  max_missed_cleavages: 2\n").expect("write config");
    let spectrum = fixture(SPECTRUM);
    let database = fixture(DATABASE);
    let pin = dir.path().join("out.pin");
    let text = printed(&[
        os("--config"),
        cfg.as_os_str(),
        os("--spectrum"),
        spectrum.as_os_str(),
        os("--database"),
        database.as_os_str(),
        os("--output-pin"),
        pin.as_os_str(),
        os("--precursor-tol"),
        os("10ppm"),
    ]);
    let line = key_line(&text, "max_missed_cleavages").expect("config key");
    assert!(!line.trim_start().starts_with('#') && line.trim_end().ends_with("# config"), "{line}");
    let line = key_line(&text, "precursor_tol").expect("cli key");
    assert!(
        !line.trim_start().starts_with('#') && line.trim_end().ends_with("# command line"),
        "{line}"
    );
    assert!(!text.contains("print_config"));
}
