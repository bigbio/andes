//! Loading the search database (FASTA, PEFF or UniProt XML, with decoys and claims) and logging
//! what it holds.

use std::path::Path;

use search::SearchIndex;

use super::super::{load_claim_database, ClaimDatabase};
use super::checks::warn_search;
use super::SearchDatabase;
use crate::cli::{Cli, EnzymeSpecificity};
use crate::models::{parse_enzymes, warn_if_universal_protease_combo};
use crate::EXPLICIT_MISSED_CLEAVAGES;

/// Missed cleavages per peptide: --max-missed-cleavages, raised to 3 under --glyco unless it was
/// given explicitly.
pub(super) fn max_missed_cleavages(cli: &Cli) -> u32 {
    let explicit = *EXPLICIT_MISSED_CLEAVAGES.get().unwrap_or(&false);
    if cli.glyco && !explicit {
        cli.max_missed_cleavages.max(3)
    } else {
        cli.max_missed_cleavages
    }
}

/// Parses the enzymes, loads the protein database and builds the search index.
pub(crate) fn load_search_database(
    cli: &Cli,
    database_path: &Path,
) -> Result<SearchDatabase, Box<dyn std::error::Error>> {
    let (search_enzyme, extra_enzymes) = parse_enzymes(&cli.db.enzyme)?;
    warn_if_universal_protease_combo(search_enzyme, &extra_enzymes);
    let t_phase = std::time::Instant::now();
    let uniprot_opts = input::UniProtOptions {
        variants: cli.db.uniprot_variants,
        isoforms: cli.db.uniprot_isoforms,
        chains: cli.db.uniprot_chains,
        evidence: cli.db.uniprot_evidence,
        knockoffs: Some(input::KnockoffOptions {
            per_variant: cli.db.knockoffs,
            enzyme: search_enzyme,
        }),
    };
    let ClaimDatabase {
        loaded,
        annotations,
        site_knockoffs,
    } = load_claim_database(database_path, &cli.db.annotations, &uniprot_opts)?;
    let target_db = loaded.db;
    if let Some(peff) = &loaded.peff {
        log_claims("PEFF", &peff.claims);
        if peff.sites_without_mass + peff.processing_not_chain > 0 {
            log::info!(target: "search",
                "PEFF: {} entries; left out {} sites with no known mass and {} processed spans \
                 that are not mature chains",
                peff.entries,
                peff.sites_without_mass,
                peff.processing_not_chain
            );
        }
    }
    if let Some(stats) = &annotations {
        log_claims("Annotations", stats);
    }
    if site_knockoffs.knockoffs > 0 {
        log::info!(
            target: "site-claims",
            "{} annotated sites, {} knockoff sites ({} sites with no position for \
             one)",
            site_knockoffs.sites, site_knockoffs.knockoffs, site_knockoffs.untested
        );
    }
    log::info!(
        target: "database",
        "Loaded {} target proteins from {} in {:.2} s",
        target_db.proteins.len(),
        database_path.display(),
        t_phase.elapsed().as_secs_f64()
    );
    match &loaded.uniprot {
        Some(stats) => log_uniprot_load(stats, &uniprot_opts),
        None => {
            for (set, flag) in [
                (cli.db.uniprot_variants, "--uniprot-variants"),
                (cli.db.uniprot_isoforms, "--uniprot-isoforms"),
                (cli.db.uniprot_chains, "--uniprot-chains"),
            ] {
                if set {
                    log::warn!(target: "search",
                        "{flag} needs a UniProt XML database (*.xml[.gz]); {} is FASTA, \
                         so it does nothing.",
                        database_path.display()
                    );
                }
            }
        }
    }

    let decoy_strategy = cli.decoy_strategy;
    let proteases = search::decoy::Proteases {
        primary: search_enzyme,
        extra: extra_enzymes.clone(),
        num_tolerable_termini: match cli.enzyme_specificity {
            EnzymeSpecificity::Fully => 2,
            EnzymeSpecificity::Semi => 1,
            EnzymeSpecificity::NonSpecific => 0,
        },
        min_length: cli.min_length,
        max_length: cli.max_length,
        max_missed_cleavages: max_missed_cleavages(cli),
    };
    let prebuilt = crate::database::prebuilt_settings(
        database_path,
        &crate::database::SearchDecoys {
            strategy: decoy_strategy,
            seed: cli.decoy_seed,
            prefix: &cli.db.decoy_prefix,
            suffix: cli.db.decoy_suffix.as_deref(),
            digest: &proteases,
        },
    )?;
    if let Some(settings) = prebuilt {
        if cli.add_entrapment {
            return Err(format!(
                "--add-entrapment cannot be used with {}, a database written by andes database: \
                 its decoys were built without the entrapment proteins. Search the target FASTA \
                 instead.",
                database_path.display()
            )
            .into());
        }
        if !cli.db.annotations.is_empty() {
            return Err(format!(
                "--annotations cannot be used with {}, a database written by andes database:                  annotated sites change the decoys. Search the target FASTA instead.",
                database_path.display()
            )
            .into());
        }
        let t_phase = std::time::Instant::now();
        let idx = SearchIndex::from_prebuilt(
            target_db,
            &cli.db.decoy_prefix,
            cli.db.decoy_suffix.as_deref(),
        );
        let target_db = idx.target_db();
        let decoys = idx.db.proteins.len() - target_db.proteins.len();
        if (target_db.proteins.len() as u64, decoys as u64) != (settings.targets, settings.decoys) {
            return Err(format!(
                "{}: {} targets and {decoys} decoys, but andes database wrote {} and {}",
                database_path.display(),
                target_db.proteins.len(),
                settings.targets,
                settings.decoys
            )
            .into());
        }
        log::info!(target: "database",
            "{} targets and {decoys} decoys ({}) written by andes database, searched as written; \
             decoy key {}",
            target_db.proteins.len(),
            settings.decoy_strategy,
            settings.decoy_key
        );
        log::info!(
            target: "search_index_build",
            "done in {:.2} s",
            t_phase.elapsed().as_secs_f64()
        );
        return Ok(SearchDatabase {
            target_db,
            idx,
            search_enzyme,
            extra_enzymes,
            proteases,
            decoy_strategy,
            decoy_key: Some(settings.decoy_key),
        });
    }
    if decoy_strategy == search::decoy::DecoyStrategy::Mutate
        && (proteases.is_non_specific() || proteases.num_tolerable_termini < 2)
    {
        return Err(
            "--decoy-strategy mutate needs fully specific digestion: a semi- or non-specific \
             search reads sub-sequences of each decoy segment, and the ones without the changed \
             residue are also target peptides. Use --decoy-strategy reverse or shuffle."
                .into(),
        );
    }
    if decoy_strategy == search::decoy::DecoyStrategy::Shuffle
        && proteases.is_non_specific()
        && target_db.proteins.iter().any(|p| {
            matches!(
                p.entry_kind(),
                model::EntryKind::Isoform | model::EntryKind::Chain
            )
        })
    {
        return Err(
            "--decoy-strategy shuffle cannot search isoform or mature-chain entries \
                    (--uniprot-isoforms, --uniprot-chains) non-specifically: shuffled \
                    whole-sequence decoys do not keep the peptides these entries share with \
                    their canonical entry. Use --decoy-strategy reverse."
                .into(),
        );
    }
    let mut decoy_key = None;
    let target_db = with_entrapment(cli, target_db, &proteases)?;
    if decoy_strategy != search::decoy::DecoyStrategy::None {
        if let Some(warning) = search::decoy::ForeignDecoys::scan(&target_db).verdict()? {
            warn_search(&warning);
        }
        let suffix = cli.db.decoy_suffix.as_deref().filter(|s| !s.is_empty());
        let key = search::decoy::decoy_key(
            &search::decoy::targets_sha256(&target_db, &cli.db.decoy_prefix, suffix),
            decoy_strategy,
            cli.decoy_seed,
            &cli.db.decoy_prefix,
            suffix,
            std::slice::from_ref(&proteases),
        );
        log::info!(target: "database",
            "decoy key {key} (the key of these decoys; andes database writes the same decoys \
             for these settings)"
        );
        decoy_key = Some(key);
    }
    let t_phase = std::time::Instant::now();
    let idx = SearchIndex::from_target_db_for(
        &target_db,
        &cli.db.decoy_prefix,
        cli.db.decoy_suffix.as_deref(),
        decoy_strategy,
        cli.decoy_seed,
        &proteases,
    );
    log::info!(
        target: "search_index_build",
        "done in {:.2} s",
        t_phase.elapsed().as_secs_f64()
    );
    Ok(SearchDatabase {
        target_db,
        idx,
        search_enzyme,
        extra_enzymes,
        proteases,
        decoy_strategy,
        decoy_key,
    })
}

/// The target database with one shuffled entrapment protein per target appended under
/// `--add-entrapment` (`search::decoy::entrapment`), logged with its size and redraws; the
/// database itself otherwise.
fn with_entrapment(
    cli: &Cli,
    mut target_db: model::ProteinDb,
    proteases: &search::decoy::Proteases,
) -> Result<model::ProteinDb, String> {
    use search::decoy::entrapment;
    if !cli.add_entrapment {
        return Ok(target_db);
    }
    if let Some(p) = target_db.proteins.iter().find(|p| {
        search::decoy::is_decoy_accession_affix(
            &p.accession,
            &cli.db.decoy_prefix,
            cli.db.decoy_suffix.as_deref(),
        )
    }) {
        return Err(format!(
            "--add-entrapment needs a target-only database: {} is a decoy entry. \
             Remove prebuilt decoys so target and entrapment decoys are generated together",
            p.accession
        ));
    }
    if let Some(p) = target_db
        .proteins
        .iter()
        .find(|p| entrapment::is_entrapment_accession(&p.accession))
    {
        return Err(format!(
            "--add-entrapment: the database already holds entrapment entries ({}); give the \
             target proteins only",
            p.accession
        ));
    }
    let t0 = std::time::Instant::now();
    let Some(built) = entrapment::entrapment_db(&target_db, cli.entrapment_seed, proteases) else {
        return Err(
            "--add-entrapment needs a specific digest: a non-specific search has no cleavage \
             sites to keep, so its shuffles have no paired target peptide"
                .into(),
        );
    };
    let s = built.stats;
    log::info!(target: "database",
        "{} entrapment proteins (seed {}; {} variant, isoform, chain or knockoff entries have \
         none); {} entrapment peptides against {} target peptides (r = {:.3}); {} equalled a \
         target peptide and were redrawn, {} still do; {:.2} s",
        output::statistics::thousands(s.proteins),
        cli.entrapment_seed,
        output::statistics::thousands(s.skipped),
        output::statistics::thousands(s.entrapment_peptides),
        output::statistics::thousands(s.target_peptides),
        s.peptide_ratio(),
        output::statistics::thousands(s.redrawn),
        output::statistics::thousands(s.left),
        t0.elapsed().as_secs_f64()
    );
    if let Some(p) = built.proteins.proteins.iter().find(|p| {
        search::decoy::is_decoy_accession_affix(
            &p.accession,
            &cli.db.decoy_prefix,
            cli.db.decoy_suffix.as_deref(),
        )
    }) {
        return Err(format!(
            "--add-entrapment generated accession {} that matches the configured decoy prefix or \
             suffix; choose a different --decoy-prefix or --decoy-suffix",
            p.accession
        ));
    }
    target_db.proteins.extend(built.proteins.proteins);
    Ok(target_db)
}

fn log_claims(source: &str, stats: &input::AnnotationLoadStats) {
    let dropped: Vec<String> = stats
        .dropped
        .iter()
        .map(|(reason, n)| format!("{reason} ({n})"))
        .collect();
    let files = if stats.files > 0 {
        format!(" from {} files", stats.files)
    } else {
        String::new()
    };
    log::info!(
        target: source,
        "{} rows{files}: {} variant, {} knockoff and {} chain entries, {} sites on \
         reference proteins, {} rows merged into an equal claim, {} into an equal database \
         entry{}{}",
        stats.rows,
        stats.variant_entries,
        stats.knockoff_entries,
        stats.chain_entries,
        stats.sites,
        stats.merged,
        stats.merged_into_database,
        if dropped.is_empty() {
            ""
        } else {
            "; left out: "
        },
        dropped.join(", "),
    );
    if stats.knockoff_entries > 0 {
        log::info!(
            target: source,
            "{} single-residue variants with no position for a knockoff",
            stats.variants_without_knockoff
        );
    }
    if stats.knockoffs_dropped > 0 {
        log::info!(
            target: source,
            "{} earlier knockoff entries dropped, within the knockoff distance of a \
             variant these files list",
            stats.knockoffs_dropped
        );
    }
}

fn log_uniprot_load(stats: &input::UniProtLoadStats, opts: &input::UniProtOptions) {
    log::info!(
        target: "uniprot",
        "{} entries, {} annotated modification sites, {} variant entries \
         (evidence >= {}; {} features below it, {} not fitting their sequence, \
         {} on isoforms not searched)",
        stats.entries,
        stats.mod_sites,
        stats.variant_proteins,
        opts.evidence.as_str(),
        stats.below_evidence,
        stats.unusable_features,
        stats.isoform_features,
    );
    if opts.variants && opts.knockoffs.is_some_and(|k| k.per_variant > 0) {
        log::info!(
            target: "uniprot",
            "{} knockoff variant entries ({} single-residue variants with no \
             position for a knockoff)",
            stats.knockoff_proteins, stats.variants_without_knockoff,
        );
    }
    if opts.isoforms || opts.chains {
        log::info!(
            target: "uniprot",
            "{} isoform entries ({} isoforms not rebuildable), {} mature-chain entries",
            stats.isoform_proteins, stats.unbuilt_isoforms, stats.chain_proteins,
        );
    }
    if !stats.unknown_modifications.is_empty() {
        let mut names: Vec<(&String, &usize)> = stats.unknown_modifications.iter().collect();
        names.sort_by(|a, b| b.1.cmp(a.1).then(a.0.cmp(b.0)));
        let total: usize = names.iter().map(|(_, n)| **n).sum();
        let examples: Vec<String> = names
            .iter()
            .take(5)
            .map(|(name, n)| format!("{name} ({n})"))
            .collect();
        log::info!(
            target: "uniprot",
            "{total} modification features have no single known mass and are not \
             used as sites, e.g. {}",
            examples.join(", ")
        );
    }
}

pub(super) fn report_annotated_site_mods(aa: &model::AminoAcidSet, target_db: &model::ProteinDb) {
    let mut mods: Vec<(&str, f64)> = aa
        .distinct_mods()
        .iter()
        .filter_map(|m| m.annotated_site_mass.map(|mass| (m.name.as_str(), mass)))
        .collect();
    mods.sort_by(|a, b| a.0.cmp(b.0).then(a.1.total_cmp(&b.1)));
    mods.dedup_by(|b, a| a.0 == b.0 && a.1.to_bits() == b.1.to_bits());
    for (name, mass) in mods {
        let (mut sites, mut proteins) = (0usize, 0usize);
        let same = |m: f64| (m - mass).abs() <= model::protein::SITE_MASS_TOLERANCE;
        for p in &target_db.proteins {
            let knockoffs: Vec<u32> = p
                .annotation
                .as_deref()
                .map(|a| {
                    a.site_claims
                        .iter()
                        .filter(|c| c.knockoff_of.is_some() && same(c.mass_delta))
                        .map(|c| c.position)
                        .collect()
                })
                .unwrap_or_default();
            let n = p
                .mod_sites()
                .iter()
                .filter(|s| same(s.mass_delta) && !knockoffs.contains(&s.position))
                .count();
            if n > 0 {
                sites += n;
                proteins += 1;
            }
        }
        if sites == 0 {
            log::warn!(target: "search",
                "{name} ({mass:+.4} Da) is restricted to annotated sites (restrict_to: annotated), \
                 but no target protein annotates a site of that mass, so it is never placed. \
                 Annotated sites come from a UniProt XML database (--database *.xml[.gz]); \
                 --uniprot-evidence sets which are used."
            );
        } else {
            log::info!(target: "search",
                "Annotated-site mod {name} ({mass:+.4} Da): {sites} sites on {proteins} target \
                 proteins"
            );
        }
    }
}
