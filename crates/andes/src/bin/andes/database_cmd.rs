//! `andes database`: the target + decoy database of a study or cohort, built once so that every
//! search of the cohort uses the same decoys, and the search side that reads it back.
//!
//! The decoys depend on the target proteins, the decoy strategy, seed, prefix and suffix, and on
//! the digests the peptide-level check covered (enzymes, specificity, peptide lengths, missed
//! cleavages): that is the decoy key (`search::decoy::decoy_key`). Modifications, tolerances,
//! charges and every other search setting are not part of it. The file layout is described in
//! `output::decoy_database`.

use std::io::Write;
use std::path::Path;

use model::protein::ProteinDb;
use output::decoy_database::{self, DatabaseSettings, Digest, OverlapRow};
use search::decoy::{DecoyStrategy, Proteases};

use crate::cli::{CheckDigest, DatabaseCmdArgs, EnzymeSpecificity};

/// Residues per sequence line of the written FASTA.
const FASTA_LINE: usize = 60;

fn termini(s: EnzymeSpecificity) -> u8 {
    match s {
        EnzymeSpecificity::Fully => 2,
        EnzymeSpecificity::Semi => 1,
        EnzymeSpecificity::NonSpecific => 0,
    }
}

fn proteases_of(d: &CheckDigest) -> Result<Proteases, Box<dyn std::error::Error>> {
    let (primary, extra) = crate::model_select::parse_enzymes(&d.enzyme)?;
    Ok(Proteases {
        primary,
        extra,
        num_tolerable_termini: termini(d.specificity),
        min_length: d.min_length,
        max_length: d.max_length,
        max_missed_cleavages: d.max_missed_cleavages,
    })
}

/// The digest as the settings table writes it.
pub(crate) fn digest_row(p: &Proteases) -> Digest {
    Digest {
        enzymes: p.enzymes_name(),
        specificity: p.specificity_name().to_string(),
        min_length: p.min_length,
        max_length: p.max_length,
        max_missed_cleavages: p.max_missed_cleavages,
    }
}

fn write_fasta(path: &Path, targets: &ProteinDb, decoys: &ProteinDb) -> std::io::Result<()> {
    output::write_atomic(path, |tmp| {
        let mut w = std::io::BufWriter::new(std::fs::File::create(tmp)?);
        for p in targets.proteins.iter().chain(&decoys.proteins) {
            if p.description.is_empty() {
                writeln!(w, ">{}", p.accession)?;
            } else {
                writeln!(w, ">{} {}", p.accession, p.description)?;
            }
            for line in p.sequence.chunks(FASTA_LINE) {
                w.write_all(line)?;
                w.write_all(b"\n")?;
            }
        }
        w.flush()
    })
}

/// Runs `andes database`.
pub(crate) fn run(args: DatabaseCmdArgs) -> Result<(), Box<dyn std::error::Error>> {
    let _ = rayon::ThreadPoolBuilder::new()
        .num_threads(args.threads.max(1))
        .build_global();
    let t0 = std::time::Instant::now();
    let strategy = args.decoy_strategy;
    if strategy == DecoyStrategy::None {
        return Err("--decoy-strategy none builds no decoys; choose a strategy".into());
    }
    let build = proteases_of(&CheckDigest {
        enzyme: args.enzyme.clone(),
        specificity: args.enzyme_specificity,
        min_length: args.min_length,
        max_length: args.max_length,
        max_missed_cleavages: args.max_missed_cleavages,
    })?;
    let mut digests = vec![build.clone()];
    for c in &args.check_digest {
        let d = proteases_of(c)?;
        if !digests.contains(&d) {
            digests.push(d);
        }
    }
    if strategy == DecoyStrategy::Mutate
        && digests
            .iter()
            .any(|d| d.is_non_specific() || d.num_tolerable_termini < 2)
    {
        return Err("--decoy-strategy mutate needs fully specific digests: a semi- or \
                    non-specific search reads sub-sequences of each decoy segment, and the ones \
                    without the changed residue are also target peptides."
            .into());
    }
    if input::DatabaseFormat::of_path(&args.fasta) != input::DatabaseFormat::Fasta {
        return Err(format!(
            "andes database reads a FASTA file; {} is not one",
            args.fasta.display()
        )
        .into());
    }
    let targets = input::load_protein_database(&args.fasta, &input::UniProtOptions::default())?.db;
    let suffix = args.decoy_suffix.as_deref().filter(|s| !s.is_empty());
    let marked = targets
        .proteins
        .iter()
        .filter(|p| search::decoy::is_decoy_accession_affix(&p.accession, &args.decoy_prefix, suffix))
        .count();
    if marked > 0 {
        return Err(format!(
            "{marked} of the {} entries of {} already carry the decoy prefix or suffix; \
             andes database takes target sequences only",
            targets.proteins.len(),
            args.fasta.display()
        )
        .into());
    }
    if let Some(warning) = search::decoy::ForeignDecoys::scan(&targets).verdict()? {
        log::warn!(target: "database", "{warning}");
    }
    log::info!(target: "database",
        "{} target proteins from {}; {} decoys for {}",
        targets.proteins.len(),
        args.fasta.display(),
        strategy.name(),
        describe(&digests)
    );
    let built = search::decoy::checked_decoy_db(
        &targets,
        &args.decoy_prefix,
        strategy,
        args.decoy_seed,
        &build,
        &digests[1..],
    );
    write_fasta(&args.output, &targets, &built.decoys)
        .map_err(|e| format!("{}: {e}", args.output.display()))?;
    reads_back(&args.output, &targets, &built.decoys)?;

    let targets_sha256 = search::decoy::targets_sha256(&targets, &args.decoy_prefix, suffix);
    let decoy_key = search::decoy::decoy_key(
        &targets_sha256,
        strategy,
        args.decoy_seed,
        &args.decoy_prefix,
        suffix,
        &digests,
    );
    let settings = DatabaseSettings {
        digests: digests.iter().map(digest_row).collect(),
        decoy_strategy: strategy.name().to_string(),
        decoy_seed: args.decoy_seed,
        decoy_prefix: search::normalize_decoy_prefix(&args.decoy_prefix),
        decoy_suffix: suffix.unwrap_or("").to_string(),
        targets_sha256,
        fasta_sha256: search::decoy::file_sha256(&args.output)?,
        decoy_key: decoy_key.clone(),
        targets: targets.proteins.len() as u64,
        decoys: built.decoys.proteins.len() as u64,
    };
    let settings_path = decoy_database::settings_path(&args.output);
    decoy_database::write_settings(&settings_path, &settings)?;

    let mut rows = Vec::new();
    for d in &digests {
        let by_length = search::decoy::peptide_overlap(&targets, &built.decoys, d);
        let (n, shared) = by_length
            .iter()
            .fold((0, 0), |(n, s), r| (n + r.decoy_peptides, s + r.shared_peptides));
        log::info!(target: "database",
            "{}: {shared} of {n} distinct decoy peptides equal a target peptide (I = L, {:.4} %)",
            describe(std::slice::from_ref(d)),
            100.0 * shared as f64 / n.max(1) as f64
        );
        rows.extend(by_length.into_iter().map(|r| OverlapRow {
            digest: digest_row(d),
            length: r.length,
            target_peptides: r.target_peptides,
            decoy_peptides: r.decoy_peptides,
            shared_peptides: r.shared_peptides,
        }));
    }
    let overlap_path = decoy_database::overlap_path(&args.output);
    decoy_database::write_overlap(&overlap_path, &rows, built.redrawn, built.left)?;
    log::info!(target: "database",
        "wrote {} ({} targets, {} decoys), {} and {}; decoy key {decoy_key}; {:.1} s",
        args.output.display(),
        settings.targets,
        settings.decoys,
        settings_path.display(),
        overlap_path.display(),
        t0.elapsed().as_secs_f64()
    );
    Ok(())
}

/// Checks that the written FASTA reads back as exactly the targets and decoys, in order, so a
/// search of it holds the proteins the search would generate.
fn reads_back(path: &Path, targets: &ProteinDb, decoys: &ProteinDb) -> Result<(), String> {
    let read = input::load_protein_database(path, &input::UniProtOptions::default())
        .map_err(|e| format!("{}: {e}", path.display()))?
        .db;
    let want = targets.proteins.iter().chain(&decoys.proteins);
    if read.proteins.len() != targets.proteins.len() + decoys.proteins.len() {
        return Err(format!(
            "{}: reads back {} entries, not {}",
            path.display(),
            read.proteins.len(),
            targets.proteins.len() + decoys.proteins.len()
        ));
    }
    for (got, want) in read.proteins.iter().zip(want) {
        if (&got.accession, &got.description, &got.sequence)
            != (&want.accession, &want.description, &want.sequence)
        {
            return Err(format!(
                "{}: entry {} does not read back as written (a residue or header the FASTA \
                 reader changes)",
                path.display(),
                want.accession
            ));
        }
    }
    Ok(())
}

fn describe(digests: &[Proteases]) -> String {
    digests
        .iter()
        .map(|d| {
            format!(
                "{} {} {}-{} aa, {} missed",
                d.enzymes_name(),
                d.specificity_name(),
                d.min_length,
                d.max_length,
                d.max_missed_cleavages
            )
        })
        .collect::<Vec<_>>()
        .join("; ")
}

/// The search's decoy settings, checked against a prebuilt database's settings table.
pub(crate) struct SearchDecoys<'a> {
    pub(crate) strategy: DecoyStrategy,
    pub(crate) seed: u64,
    pub(crate) prefix: &'a str,
    pub(crate) suffix: Option<&'a str>,
    pub(crate) digest: &'a Proteases,
}

/// The settings table beside `database` when there is one: the database was written by `andes
/// database`, and the search must search it with the same decoy key. The error names each
/// setting that differs.
pub(crate) fn prebuilt_settings(
    database: &Path,
    search: &SearchDecoys,
) -> Result<Option<DatabaseSettings>, String> {
    let path = decoy_database::settings_path(database);
    if !path.is_file() {
        return Ok(None);
    }
    let s = decoy_database::read_settings(&path)?;
    let fasta = search::decoy::file_sha256(database)
        .map_err(|e| format!("{}: {e}", database.display()))?;
    if fasta != s.fasta_sha256 {
        return Err(format!(
            "{} changed after andes database wrote it (its SHA-256 is not the one in {}); \
             rebuild the database",
            database.display(),
            path.display()
        ));
    }
    let mut differ = Vec::new();
    let mut field = |name: &str, db: String, run: String| {
        if db != run {
            differ.push(format!("{name}: database {db}, search {run}"));
        }
    };
    field(
        "--decoy-strategy",
        s.decoy_strategy.clone(),
        search.strategy.name().to_string(),
    );
    field("--decoy-seed", s.decoy_seed.to_string(), search.seed.to_string());
    field(
        "--decoy-prefix",
        s.decoy_prefix.clone(),
        search::normalize_decoy_prefix(search.prefix),
    );
    field(
        "--decoy-suffix",
        s.decoy_suffix.clone(),
        search.suffix.unwrap_or("").to_string(),
    );
    let run = digest_row(search.digest);
    if !s.digests.contains(&run) {
        let b = &s.digests[0];
        field("--enzyme", b.enzymes.clone(), run.enzymes.clone());
        field(
            "--enzyme-specificity",
            b.specificity.clone(),
            run.specificity.clone(),
        );
        field("--min-length", b.min_length.to_string(), run.min_length.to_string());
        field("--max-length", b.max_length.to_string(), run.max_length.to_string());
        field(
            "--max-missed-cleavages",
            b.max_missed_cleavages.to_string(),
            run.max_missed_cleavages.to_string(),
        );
        if s.digests.len() > 1 {
            differ.push(format!(
                "(the database was also checked for {} further digest(s), none of them this search's)",
                s.digests.len() - 1
            ));
        }
    }
    if differ.is_empty() {
        Ok(Some(s))
    } else {
        Err(format!(
            "{} is a decoy database written by andes database for other settings ({}); search \
             with its settings or rebuild it: {}",
            database.display(),
            path.display(),
            differ.join("; ")
        ))
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use model::enzyme::Enzyme;

    fn settings(fasta_sha256: String) -> DatabaseSettings {
        let lysc = Proteases {
            primary: Enzyme::LysC,
            ..Proteases::default()
        };
        DatabaseSettings {
            digests: vec![digest_row(&Proteases::default()), digest_row(&lysc)],
            decoy_strategy: "reverse".into(),
            decoy_seed: 42,
            decoy_prefix: "XXX".into(),
            decoy_suffix: String::new(),
            targets_sha256: "t".into(),
            fasta_sha256,
            decoy_key: "k".into(),
            targets: 1,
            decoys: 1,
        }
    }

    fn written(dir: &Path) -> std::path::PathBuf {
        let fasta = dir.join("td.fasta");
        std::fs::write(&fasta, ">P1\nPEPTIDEK\n>XXX_P1\nEDITPEPK\n").unwrap();
        let sha = search::decoy::file_sha256(&fasta).unwrap();
        decoy_database::write_settings(&decoy_database::settings_path(&fasta), &settings(sha))
            .unwrap();
        fasta
    }

    fn check(fasta: &Path, strategy: DecoyStrategy, digest: &Proteases) -> Result<(), String> {
        prebuilt_settings(
            fasta,
            &SearchDecoys {
                strategy,
                seed: 42,
                prefix: "XXX_",
                suffix: None,
                digest,
            },
        )
        .map(|_| ())
    }

    #[test]
    fn a_fasta_without_settings_is_not_prebuilt() {
        let dir = tempfile::tempdir().unwrap();
        let fasta = dir.path().join("plain.fasta");
        std::fs::write(&fasta, ">P1\nPEPTIDEK\n").unwrap();
        let got = prebuilt_settings(
            &fasta,
            &SearchDecoys {
                strategy: DecoyStrategy::Reverse,
                seed: 42,
                prefix: "XXX_",
                suffix: None,
                digest: &Proteases::default(),
            },
        );
        assert!(matches!(got, Ok(None)));
    }

    #[test]
    fn the_build_digest_and_a_checked_digest_are_accepted() {
        let dir = tempfile::tempdir().unwrap();
        let fasta = written(dir.path());
        check(&fasta, DecoyStrategy::Reverse, &Proteases::default()).unwrap();
        let lysc = Proteases {
            primary: Enzyme::LysC,
            ..Proteases::default()
        };
        check(&fasta, DecoyStrategy::Reverse, &lysc).unwrap();
    }

    #[test]
    fn each_differing_setting_is_named() {
        let dir = tempfile::tempdir().unwrap();
        let fasta = written(dir.path());
        let e = check(&fasta, DecoyStrategy::Shuffle, &Proteases::default()).unwrap_err();
        assert!(e.contains("--decoy-strategy: database reverse, search shuffle"), "{e}");
        let glyco = Proteases {
            max_missed_cleavages: 3,
            ..Proteases::default()
        };
        let e = check(&fasta, DecoyStrategy::Reverse, &glyco).unwrap_err();
        assert!(e.contains("--max-missed-cleavages: database 2, search 3"), "{e}");
        assert!(!e.contains("--enzyme:"), "{e}");
        assert!(e.contains("1 further digest"), "{e}");
    }

    #[test]
    fn a_changed_fasta_is_refused() {
        let dir = tempfile::tempdir().unwrap();
        let fasta = written(dir.path());
        std::fs::write(&fasta, ">P1\nPEPTIDEK\n>XXX_P1\nEDITPEKP\n").unwrap();
        let e = check(&fasta, DecoyStrategy::Reverse, &Proteases::default()).unwrap_err();
        assert!(e.contains("changed after andes database wrote it"), "{e}");
    }
}
