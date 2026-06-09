//! Diagnostic trace binary: scores a single scan against the same FASTA + param
//! used by the production search, prints candidate-window bounds, top-K PSMs,
//! and a per-split node_score breakdown for both Rust's top-1 and a
//! user-supplied "Java top-1" peptide. Use to localize Java/Rust scoring
//! divergences without rebuilding the full PXD001819 run.

use std::fs::File;
use std::io::BufReader;
use std::path::PathBuf;
use std::process::ExitCode;

// ─── Per-PSM JSON trace output (additive; no new deps) ─────────────────────
//
// Small hand-written JSON via `write!`. The diff harness parses on the
// Python side where stdlib `json` is sufficient.

struct TraceJson<W: std::io::Write> {
    out: W,
    first_psm: bool,
}

impl<W: std::io::Write> TraceJson<W> {
    fn new(mut out: W) -> std::io::Result<Self> {
        out.write_all(b"[\n")?;
        Ok(Self { out, first_psm: true })
    }

    fn begin_psm(
        &mut self,
        scan: i32,
        peptide: &str,
        charge: u8,
        rust_rank_score: i32,
    ) -> std::io::Result<()> {
        if !self.first_psm {
            self.out.write_all(b",\n")?;
        }
        self.first_psm = false;
        write!(
            self.out,
            "  {{\n    \"scan\": {},\n    \"peptide\": \"{}\",\n    \"charge\": {},\n    \"rust_rank_score\": {},\n    \"ions\": [",
            scan, escape_json(peptide), charge, rust_rank_score
        )
    }

    fn end_psm(&mut self) -> std::io::Result<()> {
        self.out.write_all(b"\n    ]\n  }")
    }

    #[allow(clippy::too_many_arguments)]
    fn ion(
        &mut self,
        first_ion: bool,
        ion_type: &str,
        theo_mz: f64,
        rank_assigned: Option<u32>,
        max_rank: u32,
        log_prob: f32,
        contribution: f32,
    ) -> std::io::Result<()> {
        if !first_ion {
            self.out.write_all(b",")?;
        }
        let rank_str = rank_assigned
            .map(|r| r.to_string())
            .unwrap_or_else(|| "null".to_string());
        write!(
            self.out,
            "\n      {{\"ion_type\": \"{}\", \"theo_mz\": {:.6}, \"rank\": {}, \"max_rank\": {}, \"log_prob\": {:.6}, \"contribution\": {:.6}}}",
            escape_json(ion_type), theo_mz, rank_str, max_rank, log_prob, contribution
        )
    }

    fn finish(mut self) -> std::io::Result<()> {
        self.out.write_all(b"\n]\n")
    }
}

fn escape_json(s: &str) -> String {
    s.replace('\\', "\\\\")
        .replace('"', "\\\"")
        .replace('\n', "\\n")
        .replace('\t', "\\t")
}

use clap::Parser;
use input::{FastaReader, MgfReader, MzMLReader};
use model::{
    AminoAcid, AminoAcidSetBuilder, ModLocation, Modification, PrecursorTolerance,
    ResidueSpec, Tolerance,
};
use model::mass::{nominal_from, H2O, PROTON};
use model::peptide::Peptide;
use scoring_crate::{Param, RankScorer};
use scoring_crate::scoring::{score_psm, ScoredSpectrum};
use scoring_crate::scoring::fragment_ions::ions_for_node;
use search::{enumerate_candidates, match_spectra, SearchIndex, SearchParams};

#[derive(Parser, Debug)]
#[command(name = "andes-trace", about = "Single-scan parity diagnostic for andes")]
struct Cli {
    /// Spectrum file (MGF or mzML — format auto-detected by extension).
    #[arg(long)]
    spectrum: PathBuf,
    /// Target FASTA database.
    #[arg(long)]
    database: PathBuf,
    /// Param file.
    #[arg(long)]
    param: PathBuf,
    /// Scan number to trace.
    #[arg(long)]
    scan: i32,
    /// Java top-1 peptide in `K.PEPTIDE.D` form (with flanking residues).
    /// Optional — when omitted, only Rust's top-1 is shown.
    #[arg(long)]
    java_top1: Option<String>,
    /// Decoy prefix.
    #[arg(long, default_value = "XXX")]
    decoy_prefix: String,
    /// Top-N PSMs per spectrum.
    #[arg(long, default_value = "10")]
    top_n: u32,
    /// Precursor tolerance (ppm).
    #[arg(long, default_value = "5.0")]
    precursor_tol_ppm: f64,
    /// Min isotope error.
    #[arg(long, default_value = "0")]
    isotope_error_min: i8,
    /// Max isotope error.
    #[arg(long, default_value = "1")]
    isotope_error_max: i8,
    /// Charge range min.
    #[arg(long, default_value = "2")]
    charge_min: u8,
    /// Charge range max.
    #[arg(long, default_value = "4")]
    charge_max: u8,
    /// Number of tolerable termini.
    #[arg(long, default_value = "2")]
    ntt: u8,
    /// Max missed cleavages.
    #[arg(long, default_value = "2")]
    max_missed_cleavages: u32,
    /// Min peaks.
    #[arg(long, default_value = "10")]
    min_peaks: u32,
    /// Min peptide length.
    #[arg(long, default_value = "6")]
    min_length: u32,
    /// Max peptide length.
    #[arg(long, default_value = "40")]
    max_length: u32,
    /// Output structured per-PSM per-ion JSON to this path. Additive: the
    /// existing human-readable stderr trace is unaffected.
    #[arg(long)]
    trace_json: Option<PathBuf>,
    /// Dump the post-filter, post-deconvolution active peak list (sorted by
    /// rank ascending) for this scan/charge as `rank<TAB>mz<TAB>intensity`
    /// lines, preceded by a `DUMP_PEAKS` header. Read-only diagnostic.
    #[arg(long)]
    dump_peaks: bool,
}

fn main() -> ExitCode {
    let cli = Cli::parse();
    match run(cli) {
        Ok(()) => ExitCode::SUCCESS,
        Err(e) => {
            eprintln!("andes-trace: {e}");
            ExitCode::from(1)
        }
    }
}

fn run(cli: Cli) -> Result<(), Box<dyn std::error::Error>> {
    // Load target db, build target+decoy SearchIndex.
    let target_db = FastaReader::load_all(BufReader::new(File::open(&cli.database)?))?;
    let idx = SearchIndex::from_target_db(&target_db, &cli.decoy_prefix);
    println!(
        "DB: {} target proteins, {} total (target+decoy)",
        target_db.proteins.len(),
        idx.db.proteins.len()
    );

    // Build aa_set with standard mods (CAM fixed C, Oxidation variable M).
    let cam = Modification {
        name: "Carbamidomethyl".into(),
        mass_delta: 57.02146,
        residue: ResidueSpec::Specific(b'C'),
        location: ModLocation::Anywhere,
        fixed: true,
        accession: None,
    };
    let ox = Modification {
        name: "Oxidation".into(),
        mass_delta: 15.99491,
        residue: ResidueSpec::Specific(b'M'),
        location: ModLocation::Anywhere,
        fixed: false,
        accession: None,
    };
    let aa = AminoAcidSetBuilder::new_standard()
        .add_fixed_mod(cam)
        .add_variable_mod(ox)
        .build()?;

    // Param + scorer.
    let param = Param::load_from_file(&cli.param)?;
    let scorer = RankScorer::new(&param);
    println!(
        "Param: activation={:?} instrument={:?} mme={:?} num_segments={} num_partitions={} error_scaling_factor={} max_rank={}",
        param.data_type.activation,
        param.data_type.instrument,
        param.mme,
        param.num_segments,
        param.partitions.len(),
        param.error_scaling_factor,
        param.max_rank
    );
    // Dump rank_dist values for the FIRST partition's first non-noise ion +
    // Noise frequencies, so we can compare against expected Java output.
    if let Some((part, ion_table)) = param.rank_dist_table.iter().next() {
        println!("\n  --- Sample rank_dist (partition {:?}) ---", part);
        let noise = ion_table.get(&scoring_crate::param_model::IonType::Noise);
        if let Some(noise) = noise {
            println!("    Noise freqs (first 5 ranks): {:?}", &noise[..5.min(noise.len())]);
            println!("    Noise freq at max_rank ({}): {}", param.max_rank, noise[param.max_rank as usize]);
        }
        for (ion, freqs) in ion_table.iter().take(3) {
            if matches!(ion, scoring_crate::param_model::IonType::Noise) { continue; }
            println!("    Ion {:?}: first 5 freqs = {:?}", ion, &freqs[..5.min(freqs.len())]);
            println!("                missing slot ({}): {}", param.max_rank, freqs[param.max_rank as usize]);
        }
        // Sanity: dump scorer.node_score for a known (partition, ion, rank).
        if let Some((ion, _)) = ion_table.iter().find(|(i, _)| !matches!(i, scoring_crate::param_model::IonType::Noise)) {
            for rank in [1, 5, 20, 100, 150] {
                let s = scorer.node_score(*part, *ion, rank);
                println!("    scorer.node_score({:?}, rank={}) = {:.4}", ion, rank, s);
            }
            let miss = scorer.missing_ion_score(*part, *ion);
            println!("    scorer.missing_ion_score = {:.4}", miss);
        }
    }
    // Diagnostic: ion type counts per (segment, all-partitions-union) vs per-partition-only.
    // Rust's `ions_for_node` iterates the union; Java's NewScoredSpectrum iterates per-partition.
    for seg in 0..param.num_segments as usize {
        let union_ions = param.ion_types_for_segment(seg);
        let prefix_n = union_ions.iter().filter(|i| matches!(i, scoring_crate::param_model::IonType::Prefix { .. })).count();
        let suffix_n = union_ions.iter().filter(|i| matches!(i, scoring_crate::param_model::IonType::Suffix { .. })).count();
        println!(
            "  seg={}: ion_types_for_segment(union) = {} ion types (prefix={}, suffix={})",
            seg, union_ions.len(), prefix_n, suffix_n
        );
    }
    // Count partitions per (charge, seg) so we know how much the union differs from a single partition.
    let mut partition_counts: std::collections::BTreeMap<(i32, i32), usize> = std::collections::BTreeMap::new();
    for p in &param.partitions {
        *partition_counts.entry((p.charge, p.seg_num)).or_insert(0) += 1;
    }
    println!("  Partition counts per (charge, seg):");
    for ((c, s), n) in &partition_counts {
        println!("    charge={} seg={}: {} partitions", c, s, n);
    }
    if std::env::var_os("MSGF_TRACE_DUMP_PARTITIONS").is_some() {
        println!("  ALL partitions (idx, c, pm, seg):");
        for (i, part) in param.partitions.iter().enumerate() {
            println!("    [{}] c={} pm={} seg={}", i, part.charge, part.parent_mass, part.seg_num);
        }
    }
    // Show distinct ion-type-list sizes across all partitions in (charge=2, seg=0).
    use std::collections::HashSet;
    for (c, s) in [(2_i32, 0_i32), (2, 1)] {
        let mut sizes: Vec<usize> = Vec::new();
        let mut union: HashSet<scoring_crate::param_model::IonType> = HashSet::new();
        for p in &param.partitions {
            if p.charge != c || p.seg_num != s { continue; }
            if let Some(frag_list) = param.frag_off_table.get(p) {
                let n = frag_list.iter()
                    .filter(|f| !matches!(f.ion_type, scoring_crate::param_model::IonType::Noise))
                    .count();
                sizes.push(n);
                for f in frag_list {
                    if !matches!(f.ion_type, scoring_crate::param_model::IonType::Noise) {
                        union.insert(f.ion_type);
                    }
                }
            }
        }
        sizes.sort();
        let len = sizes.len();
        let min_n = sizes.first().copied().unwrap_or(0);
        let max_n = sizes.last().copied().unwrap_or(0);
        let median = if len > 0 { sizes[len / 2] } else { 0 };
        println!(
            "    charge={} seg={}: per-partition ion-list sizes min={} median={} max={}, union={}",
            c, s, min_n, median, max_n, union.len()
        );
    }

    // Load just the requested scan. Auto-detect format by file extension:
    // `.mzML`/`.mzml` → MzMLReader; anything else (e.g. `.mgf`) → MgfReader.
    // For MGF specifically, fall back to extracting `scan=N` from the TITLE
    // line when the reader did not populate `Spectrum::scan` (the BSA parity
    // fixture `test.mgf` has no `SCANS=` field — scan is only encoded in
    // TITLE, matching what `gf_java_parity.rs` does).
    let ext = cli
        .spectrum
        .extension()
        .and_then(|e| e.to_str())
        .map(|s| s.to_lowercase());
    let mut spectra = Vec::new();
    match ext.as_deref() {
        Some("mzml") => {
            let reader = MzMLReader::new(BufReader::new(File::open(&cli.spectrum)?));
            for r in reader {
                let s = r?;
                if s.scan == Some(cli.scan) {
                    spectra.push(s);
                    break;
                }
            }
        }
        _ => {
            // MGF (default / backwards-compatible)
            let reader = MgfReader::new(BufReader::new(File::open(&cli.spectrum)?));
            for r in reader {
                let s = r?;
                let resolved_scan = s
                    .scan
                    .or_else(|| extract_scan_from_title(&s.title));
                if resolved_scan == Some(cli.scan) {
                    spectra.push(s);
                    break;
                }
            }
        }
    }
    if spectra.is_empty() {
        return Err(format!("scan {} not found in {}", cli.scan, cli.spectrum.display()).into());
    }
    let spec = &spectra[0];
    println!(
        "\n=== Spectrum: scan={} precursor_mz={} charge={:?} peaks={} ===",
        cli.scan,
        spec.precursor_mz,
        spec.precursor_charge,
        spec.peaks.len()
    );
    // Per-spectrum partition diagnostic: which partition (and ion list)
    // does THIS spectrum hit for each segment?
    if let Some(z_raw) = spec.precursor_charge {
        let z = z_raw.max(1) as u8;
        let pm = (spec.precursor_mz - PROTON) * z as f64;
        for s in 0..param.num_segments as usize {
            let ion_list = param.ion_types_for_partition(z, pm, s);
            let selected = param.partition_for(z, pm, s);
            println!(
                "  spectrum partition target=(c={} pm={:.2} seg={}) selected=(c={} pm={:.2} seg={}): {} ion types — {:?}",
                z, pm, s,
                selected.charge, selected.parent_mass, selected.seg_num,
                ion_list.len(),
                ion_list.iter().map(|i| match i {
                    scoring_crate::param_model::IonType::Prefix { charge, offset_bits } => format!("P(c={},off={:.3})", charge, f32::from_bits(*offset_bits)),
                    scoring_crate::param_model::IonType::Suffix { charge, offset_bits } => format!("S(c={},off={:.3})", charge, f32::from_bits(*offset_bits)),
                    scoring_crate::param_model::IonType::Noise => "Noise".to_string(),
                }).collect::<Vec<_>>()
            );
        }

        // Hypothesis #1 diagnostic: how many peaks does Rust filter for this
        // spectrum, and what filter m/z values does it use? Java filters by
        // SETTING INTENSITY=0 (peak survives ranking but ranks last), Rust
        // EXCLUDES filtered peaks from ranking entirely. If Rust filters more
        // peaks, ranks shift downward more for the survivors, lowering the
        // log_score lookups for matched ions on long peptides.
        let filter_entries = param.precursor_off_map.get(&(z as i32))
            .map(Vec::as_slice).unwrap_or(&[]);
        let neutral_mass = (spec.precursor_mz - PROTON) * z as f64;
        let mut filter_mzs: Vec<(f64, f64)> = Vec::new();
        for pof in filter_entries {
            let c = (z as i32 - pof.reduced_charge) as f64;
            if c <= 0.0 { continue; }
            let filter_mz = (neutral_mass + c * PROTON) / c + (pof.offset as f64);
            let tol_da = pof.tolerance.as_da(filter_mz);
            filter_mzs.push((filter_mz, tol_da));
        }
        // Determine which peaks would be filtered by Rust's logic.
        let mut n_filtered = 0;
        let mut max_filtered_intensity: f32 = 0.0;
        let mut filtered_examples: Vec<(f64, f32)> = Vec::new();
        for &(mz, intensity) in &spec.peaks {
            let filtered = filter_mzs.iter().any(|&(fmz, tol)| (mz - fmz).abs() <= tol);
            if filtered {
                n_filtered += 1;
                if intensity > max_filtered_intensity {
                    max_filtered_intensity = intensity;
                }
                if filtered_examples.len() < 5 {
                    filtered_examples.push((mz, intensity));
                }
            }
        }
        println!(
            "  Rust filtering: {} of {} peaks filtered ({:.1}%); max filtered intensity={:.1}",
            n_filtered, spec.peaks.len(),
            100.0 * n_filtered as f64 / spec.peaks.len() as f64,
            max_filtered_intensity
        );
        println!("  Filter m/z values (count={}):", filter_mzs.len());
        for (fmz, tol) in &filter_mzs {
            println!("    {:.4} ± {:.4}", fmz, tol);
        }
        if !filtered_examples.is_empty() {
            println!("  First 5 filtered peaks:");
            for (mz, intensity) in &filtered_examples {
                println!("    mz={:.4} intensity={:.1}", mz, intensity);
            }
        }
    }

    // Build search params (same as production harness).
    let mut params = SearchParams::default_tryptic(aa);
    params.precursor_tolerance = PrecursorTolerance::symmetric(Tolerance::Ppm(cli.precursor_tol_ppm));
    params.charge_range = cli.charge_min..=cli.charge_max;
    params.isotope_error_range = cli.isotope_error_min..=cli.isotope_error_max;
    params.top_n_psms_per_spectrum = cli.top_n;
    params.num_tolerable_termini = cli.ntt;
    params.max_missed_cleavages = cli.max_missed_cleavages;
    params.min_peaks = cli.min_peaks;
    params.min_length = cli.min_length;
    params.max_length = cli.max_length;

    // Charges to try.
    let charges_to_try: Vec<u8> = match spec.precursor_charge {
        Some(z) if z > 0 => vec![z as u8],
        _ => params.charge_range.clone().collect(),
    };

    // Print candidate-window bounds per charge, mirroring match_engine.rs.
    println!("\n--- Candidate windows ---");
    for &z in &charges_to_try {
        let charge_f = z as f64;
        let neutral_mass = (spec.precursor_mz - PROTON) * charge_f - H2O;
        let nominal_center = nominal_from(neutral_mass);
        let iso_min = *params.isotope_error_range.start() as i32;
        let iso_max = *params.isotope_error_range.end() as i32;
        let tol_da_left = params.precursor_tolerance.left.as_da(neutral_mass);
        let tol_da_right = params.precursor_tolerance.right.as_da(neutral_mass);
        let widen_left = (tol_da_left - 0.4999_f64).round() as i32;
        let widen_right = (tol_da_right - 0.4999_f64).round() as i32;
        let min_nominal = nominal_center - iso_max - widen_right;
        let max_nominal = nominal_center - iso_min + widen_left;
        println!(
            "  charge={}: neutral_mass={:.4} nominal_center={} window=[{}..={}] (iso_range=[{}..={}], tol_da_left={:.4}, tol_da_right={:.4})",
            z, neutral_mass, nominal_center, min_nominal, max_nominal,
            iso_min, iso_max, tol_da_left, tol_da_right
        );
    }

    // Run the full search on this single spectrum.
    let (queues, run_candidates) = match_spectra(&spectra, &idx, &params, &scorer, 0.5, &cli.decoy_prefix);
    let queue = &queues[0];
    let psms: Vec<_> = queue.iter_psms().collect();

    // Print top-K Rust PSMs.
    println!("\n--- Rust top-{} PSMs ---", psms.len());
    let mut sorted: Vec<&_> = psms.iter().collect();
    sorted.sort_by(|a, b| b.score.partial_cmp(&a.score).unwrap_or(std::cmp::Ordering::Equal));
    for (i, psm) in sorted.iter().enumerate() {
        let cand = &run_candidates[psm.primary_candidate_idx() as usize];
        let prot = idx.protein_at(cand.protein_index);
        let prot_acc = prot.map(|p| p.accession.as_str()).unwrap_or("?");
        let is_decoy = cand.is_decoy;
        let pep_str: String = cand.peptide.residues.iter()
            .map(|aa| aa.residue as char)
            .collect();
        println!(
            "  #{}: peptide={} charge={} score={:.2} rank_score={:.2} iso_off={} prot_idx={} prot={} is_decoy={}",
            i + 1, pep_str, psm.charge_used, psm.score, psm.rank_score,
            psm.isotope_offset, cand.protein_index, prot_acc, is_decoy
        );
    }

    // Set up optional structured JSON trace output.
    let mut trace_json: Option<TraceJson<std::io::BufWriter<File>>> = match cli.trace_json {
        Some(ref path) => {
            let file = File::create(path).map_err(|e| {
                eprintln!("Failed to create --trace-json output {}: {}", path.display(), e);
                e
            })?;
            Some(TraceJson::new(std::io::BufWriter::new(file))?)
        }
        None => None,
    };

    // If user supplied Java top-1, search for it in Rust's enumerated set.
    if let Some(java_str) = &cli.java_top1 {
        let java_pep = parse_flanking(java_str)?;
        println!("\n--- Java top-1 trace: {} ---", java_str);

        // Enumerate all candidates (Rust's view) and search for an exact-residue match.
        let java_residues: Vec<u8> = java_pep.residues.iter().map(|aa| aa.residue).collect();
        let mut found_indices: Vec<usize> = Vec::new();
        let cands: Vec<_> = enumerate_candidates(&idx, &params, &cli.decoy_prefix).collect();
        for (i, c) in cands.iter().enumerate() {
            let cand_residues: Vec<u8> = c.peptide.residues.iter().map(|aa| aa.residue).collect();
            if cand_residues == java_residues {
                found_indices.push(i);
            }
        }
        println!("  Enumerator: {} matches for residue sequence", found_indices.len());
        for &i in found_indices.iter().take(5) {
            let c = &cands[i];
            let prot = idx.protein_at(c.protein_index);
            let prot_acc = prot.map(|p| p.accession.as_str()).unwrap_or("?");
            println!(
                "    cand_idx={} prot_idx={} prot={} is_decoy={} pep_mass={:.4} nominal={}",
                i, c.protein_index, prot_acc, c.is_decoy, c.peptide.mass(),
                c.peptide.nominal_residue_mass()
            );
        }
        if found_indices.is_empty() {
            println!("  WARNING: Java top-1 NOT in Rust's enumerated candidate set (window or enumeration gap)");
        }

        // Check if any of these enumerated candidates are in Rust's top-N queue.
        let in_queue: usize = psms.iter().filter(|psm| {
            let cand = &run_candidates[psm.primary_candidate_idx() as usize];
            let pep_residues: Vec<u8> = cand.peptide.residues.iter()
                .map(|aa| aa.residue).collect();
            pep_residues == java_residues
        }).count();
        println!("  In Rust's top-{} queue: {}", psms.len(), in_queue);

        // Per-split node_score breakdown for Java's peptide.
        // Use the first found candidate to get correct flanking.
        if !found_indices.is_empty() {
            let java_cand_pep = &cands[found_indices[0]].peptide;
            for &z in &charges_to_try {
                println!("\n  Per-split node_score breakdown — Java pep ({}+{}) ---", java_str, z);
                let scored = ScoredSpectrum::new(spec, &scorer, z);
                let total = score_psm(&scored, java_cand_pep, &scorer, z, 0.5);
                print_split_breakdown(
                    &scored,
                    java_cand_pep,
                    &scorer,
                    z,
                    trace_json.as_mut(),
                    cli.scan,
                    java_str,
                    total.round() as i32,
                )?;
                println!("    score_psm total = {}", total);
            }
        }
    }

    // Per-split node_score breakdown for Rust's top-1.
    if let Some(top1) = sorted.first() {
        let rust_top1_pep = &run_candidates[top1.primary_candidate_idx() as usize].peptide;
        let pep_str: String = rust_top1_pep.residues.iter().map(|aa| aa.residue as char).collect();
        println!("\n  Per-split node_score breakdown — Rust top-1 ({} +{}) ---", pep_str, top1.charge_used);
        let scored = ScoredSpectrum::new(spec, &scorer, top1.charge_used);
        let rust_rank_score = top1.score.round() as i32;
        print_split_breakdown(
            &scored,
            rust_top1_pep,
            &scorer,
            top1.charge_used,
            trace_json.as_mut(),
            cli.scan,
            &pep_str,
            rust_rank_score,
        )?;
        println!("    PSM.score (from queue) = {}", top1.score);
    }

    // ---------------------------------------------------------------------
    // Diagnostic: dump the active (post-filter, post-deconvolution) peak list
    // sorted by rank ascending. Read-only; uses the SAME peak/rank set the
    // scorer consumes. Lets us compare Rust's kept-peak ranks against Java's.
    // ---------------------------------------------------------------------
    if cli.dump_peaks {
        let dump_charge = charges_to_try.first().copied().unwrap_or(cli.charge_min);
        let scored = ScoredSpectrum::new(spec, &scorer, dump_charge);
        let active = scored.dump_active_peaks();
        println!(
            "DUMP_PEAKS scan={} charge={} precursor_mz={} active_peaks={}",
            cli.scan,
            dump_charge,
            spec.precursor_mz,
            active.len()
        );
        for (rank, mz, intensity) in &active {
            println!("{}\t{:.5}\t{:.2}", rank, mz, intensity);
        }
    }

    // Quick view of the spectrum's top-10 peaks by intensity.
    println!("\n--- Spectrum top-10 peaks by intensity ---");
    let mut peaks_by_int: Vec<_> = spec.peaks.iter().enumerate().collect();
    peaks_by_int.sort_by(|a, b| b.1.1.partial_cmp(&a.1.1).unwrap_or(std::cmp::Ordering::Equal));
    for (rank, (_idx, &(mz, intensity))) in peaks_by_int.iter().take(10).enumerate() {
        println!("  rank={} mz={:.4} intensity={}", rank + 1, mz, intensity);
    }

    if let Some(tj) = trace_json {
        tj.finish().map_err(|e| {
            eprintln!("Failed to finalize --trace-json output: {}", e);
            e
        })?;
    }

    Ok(())
}

/// Extract `scan=N` from an MGF TITLE string (e.g. mzML
/// `controllerType=0 controllerNumber=1 scan=3416`). Mirrors the helper in
/// `crates/search/tests/gf_java_parity.rs` — required because the BSA parity
/// fixture `test.mgf` has no `SCANS=` line, so `Spectrum::scan` is `None`.
fn extract_scan_from_title(title: &str) -> Option<i32> {
    title
        .split_ascii_whitespace()
        .find_map(|tok| tok.strip_prefix("scan=")?.parse::<i32>().ok())
}

/// Parse a peptide string in `K.PEPTIDE.D` form.
fn parse_flanking(s: &str) -> Result<Peptide, Box<dyn std::error::Error>> {
    let parts: Vec<&str> = s.split('.').collect();
    if parts.len() != 3 {
        return Err(format!("expected K.PEPTIDE.D form, got: {s}").into());
    }
    // Empty flanking residues (e.g. ".PEPTIDE.") fall back to the "no flank"
    // marker `b'-'` instead of indexing an empty slice (which would panic).
    let pre = parts[0].bytes().next().unwrap_or(b'-');
    let post = parts[2].bytes().next().unwrap_or(b'-');
    let body = parts[1];
    // Strip mod annotations like "C+57.021" → "C". Simple heuristic: keep only A-Z.
    let residues: Vec<AminoAcid> = body
        .bytes()
        .filter(|&b| b.is_ascii_uppercase())
        .map(|b| {
            AminoAcid::standard(b)
                .ok_or_else(|| format!("unknown residue: {}", b as char))
        })
        .collect::<Result<_, _>>()?;
    Ok(Peptide::new(residues, pre, post))
}

/// Print per-split node_score: prefix nominal, suffix nominal, score per split,
/// and which ions matched peaks.
///
/// When `trace_json` is `Some`, emits a structured JSON record for this PSM
/// alongside the existing human-readable output.
#[allow(clippy::too_many_arguments)]
fn print_split_breakdown(
    scored: &ScoredSpectrum<'_>,
    peptide: &Peptide,
    scorer: &RankScorer,
    charge: u8,
    mut trace_json: Option<&mut TraceJson<std::io::BufWriter<File>>>,
    scan: i32,
    peptide_label: &str,
    rank_score: i32,
) -> Result<(), Box<dyn std::error::Error>> {
    let n = peptide.length();
    if n < 2 { return Ok(()); }
    // Use SPECTRUM's parent mass for partition lookup (matching score_psm fix).
    let spectrum_parent_mass = scored.parent_mass();
    let peptide_mass = peptide.mass();
    let peptide_nominal = peptide.nominal_residue_mass();
    let mut prefix_acc = 0.0_f64;
    let mut total: i32 = 0;
    let mme = &scorer.param().mme;
    let max_rank = scorer.max_rank();

    // Begin JSON PSM record if a writer is present.
    if let Some(ref mut tj) = trace_json {
        tj.begin_psm(scan, peptide_label, charge, rank_score)?;
    }
    let mut first_json_ion = true;

    println!("    spectrum_parent_mass={:.4}, peptide_mass={:.4}, peptide_nominal={}",
        spectrum_parent_mass, peptide_mass, peptide_nominal);
    for s in 1..n {
        let aa = &peptide.residues[s - 1];
        let residue_mass = aa.mass + aa.mod_.as_ref().map_or(0.0, |m| m.mass_delta);
        prefix_acc += residue_mass;
        let prefix_nominal = nominal_from(prefix_acc);
        let suffix_nominal = peptide_nominal - prefix_nominal;

        // Collect detailed per-ion contributions to compare against Java.
        let mut ion_details: Vec<String> = Vec::new();
        let mut matched_sum: f32 = 0.0;
        let mut missing_sum: f32 = 0.0;
        let mut n_matched = 0;
        let mut n_missing = 0;
        for is_prefix in [true, false] {
            let nom = if is_prefix { prefix_nominal as f64 } else { suffix_nominal as f64 };
            for (ion, theo_mz) in ions_for_node(nom, is_prefix, scorer.param(), spectrum_parent_mass, charge) {
                let seg = scorer.param().segment_num(theo_mz, spectrum_parent_mass);
                let part = scorer.param().partition_for(charge, spectrum_parent_mass, seg);
                let tol_da = mme.as_da(theo_mz);
                let peak_rank = scored.nearest_peak_rank(theo_mz, tol_da);
                let (score_str, contribution, log_prob) = match peak_rank {
                    Some(rank) => {
                        let s = scorer.node_score(part, ion, rank);
                        n_matched += 1;
                        matched_sum += s;
                        (format!("rk{}={:.2}", rank, s), s, s)
                    }
                    None => {
                        let s = scorer.missing_ion_score(part, ion);
                        n_missing += 1;
                        missing_sum += s;
                        (format!("MISS={:.2}", s), s, s)
                    }
                };
                // Emit JSON ion record if writer is present.
                if let Some(ref mut tj) = trace_json {
                    tj.ion(
                        first_json_ion,
                        &format!("{:?}", ion),
                        theo_mz,
                        peak_rank,
                        max_rank,
                        log_prob,
                        contribution,
                    )?;
                    first_json_ion = false;
                }
                let kind = if is_prefix { "P" } else { "S" };
                let off = match ion {
                    scoring_crate::param_model::IonType::Prefix { offset_bits, .. } |
                    scoring_crate::param_model::IonType::Suffix { offset_bits, .. } => f32::from_bits(offset_bits),
                    _ => 0.0,
                };
                ion_details.push(format!("{}{:.1}@{:.1}={}", kind, off, theo_mz, score_str));
            }
        }
        let split_score = (matched_sum + missing_sum).round() as i32;
        total += split_score;

        let resi_char = aa.residue as char;
        println!(
            "    split={} aa[{}]={} pref_nom={} suf_nom={} score={} (matched={} sum={:.2}, missing={} sum={:.2})",
            s, s - 1, resi_char, prefix_nominal, suffix_nominal, split_score,
            n_matched, matched_sum, n_missing, missing_sum
        );
        if s == 4 || s == 1 {
            // Show full per-ion breakdown for selected splits.
            println!("      ions: {}", ion_details.join(" | "));
        }
    }
    println!("    breakdown_total = {}", total);

    // Close JSON PSM record if a writer is present.
    if let Some(ref mut tj) = trace_json {
        tj.end_psm()?;
    }

    Ok(())
}
