//! Value types of the command line: the enum flags, comma lists, modification lists, crosslinker
//! sets, and the value parsers that read them.

use clap::builder::{PossibleValue, TypedValueParser};
use clap::ValueEnum;
use search::decoy::DecoyStrategy;
use search::PrecursorCalMode;

/// Fragmentation method of MGF spectra.
#[derive(Clone, Copy, Debug, PartialEq, Eq, ValueEnum)]
pub(crate) enum Fragmentation {
    #[value(name = "auto")]
    Auto,
    #[value(name = "CID")]
    Cid,
    #[value(name = "ETD")]
    Etd,
    #[value(name = "HCD")]
    Hcd,
    #[value(name = "UVPD")]
    Uvpd,
}

/// Search protocol: the labelling or enrichment of the sample, which selects
/// the protocol model.
#[derive(Clone, Copy, Debug, PartialEq, Eq, ValueEnum)]
pub(crate) enum Protocol {
    #[value(name = "auto")]
    Auto,
    #[value(name = "phospho")]
    Phospho,
    #[value(name = "iTRAQ")]
    Itraq,
    #[value(name = "iTRAQ-phospho")]
    ItraqPhospho,
    #[value(name = "TMT")]
    Tmt,
    #[value(name = "standard")]
    Standard,
}

/// `--tmt-unmix`: whether and how strongly the co-isolation interference of
/// high-resolution MS2 reporter ions is removed.
#[derive(Clone, Copy, Debug, Default, PartialEq, Eq, ValueEnum)]
pub(crate) enum TmtUnmix {
    /// Remove it at a strength learned per run (a no-op without --tmt, and
    /// for low-resolution MS2 or MS3 reporters)
    #[default]
    #[value(name = "auto")]
    Auto,
    /// Write no interference estimate
    #[value(name = "off")]
    Off,
    /// Remove all of the estimated interference
    #[value(name = "on")]
    On,
}

/// `--enzyme-specificity`: the peptide ends that must be cleavage sites.
#[derive(Clone, Copy, Debug, PartialEq, Eq, ValueEnum)]
pub(crate) enum EnzymeSpecificity {
    /// Neither end
    #[value(name = "non-specific")]
    NonSpecific,
    /// At least one end
    #[value(name = "semi")]
    Semi,
    /// Both ends
    #[value(name = "fully")]
    Fully,
}

/// `--score`: what ranks the candidates and fills the PIN `RawScore`.
#[derive(Clone, Copy, Debug, Default, PartialEq, Eq, ValueEnum)]
pub(crate) enum ScoreFlag {
    /// strong for a high-resolution model, rank for a low-resolution one
    #[default]
    Auto,
    /// The rank score of the model
    Rank,
    /// The fused intensity and competition score
    Strong,
}

/// `--fragment-index`: candidate retrieval.
#[derive(Clone, Copy, Debug, Default, PartialEq, Eq, ValueEnum)]
pub(crate) enum FragmentIndexFlag {
    /// The fragment-ion index (--glyco looks candidates up by backbone mass)
    #[default]
    Auto,
    /// As auto; with --glyco, a warning and per-spectrum enumeration
    On,
    /// Score every candidate in each precursor window
    Off,
}

/// One `--predict` entry: a predicted-feature block added to the PIN.
#[derive(Clone, Copy, Debug, PartialEq, Eq, Hash, PartialOrd, Ord, ValueEnum)]
pub(crate) enum PredictFeature {
    /// Retention-time deviation from the bundled andes-rt model
    #[value(name = "rt")]
    Rt,
    /// Spectral similarity to the fragment intensities of the bundled andes-ms2 model
    #[value(name = "ms2")]
    Ms2,
    /// Spectral similarity to the fragment intensities of AlphaPeptDeep
    #[value(name = "peptdeep-ms2")]
    PeptDeepMs2,
    /// Retention-time deviation from AlphaPeptDeep
    #[value(name = "peptdeep-rt")]
    PeptDeepRt,
    /// Retention-time deviation from DeepLC
    #[value(name = "deeplc-rt")]
    DeepLcRt,
    /// No predicted features
    #[value(name = "none")]
    None,
}

/// `--predictor-backend`: ONNX engine of the third-party predictors.
#[derive(Clone, Copy, Debug, PartialEq, Eq, ValueEnum)]
pub(crate) enum PredictorBackend {
    /// tract (pure Rust, built in)
    Tract,
    /// ONNX Runtime (a build with --features ort, plus the runtime library)
    Ort,
}

#[derive(Clone)]
pub(crate) struct Listed<T> {
    pub(super) values: &'static [(&'static str, &'static str)],
    pub(super) parse: fn(&str) -> Result<T, String>,
}

impl<T> Listed<T> {
    pub(crate) const fn new(
        values: &'static [(&'static str, &'static str)],
        parse: fn(&str) -> Result<T, String>,
    ) -> Self {
        Self { values, parse }
    }
}

impl<T: Clone + Send + Sync + 'static> TypedValueParser for Listed<T> {
    type Value = T;

    fn parse_ref(
        &self,
        cmd: &clap::Command,
        arg: Option<&clap::Arg>,
        value: &std::ffi::OsStr,
    ) -> Result<T, clap::Error> {
        self.parse.parse_ref(cmd, arg, value)
    }

    fn possible_values(&self) -> Option<Box<dyn Iterator<Item = PossibleValue> + '_>> {
        Some(Box::new(self.values.iter().map(|&(name, help)| {
            let value = PossibleValue::new(name);
            if help.is_empty() {
                value
            } else {
                value.help(help)
            }
        })))
    }
}

#[derive(Debug, Clone, Default)]
pub(crate) struct ModList(pub(crate) Vec<model::Modification>);

pub(crate) fn parse_fixed_mods(s: &str) -> Result<ModList, String> {
    model::modification::spec::parse_mod_list(s, true)
        .map(ModList)
        .map_err(|e| e.to_string())
}

pub(crate) fn parse_variable_mods(s: &str) -> Result<ModList, String> {
    model::modification::spec::parse_mod_list(s, false)
        .map(ModList)
        .map_err(|e| e.to_string())
}

#[derive(Debug, Clone)]
pub(crate) struct XlSet(pub(crate) Vec<xl::crosslinker::Crosslinker>);

impl XlSet {
    /// One or more crosslinkers with distinct names; `key` names the option in an error.
    pub(crate) fn new(
        key: &str,
        linkers: Vec<xl::crosslinker::Crosslinker>,
    ) -> Result<Self, String> {
        if linkers.is_empty() {
            return Err(format!("{key}: no crosslinker"));
        }
        for (i, a) in linkers.iter().enumerate() {
            if linkers[..i].iter().any(|b| b.name == a.name) {
                return Err(format!("{key}: {} is given twice", a.name));
            }
        }
        Ok(Self(linkers))
    }
}

#[derive(Copy, Clone, Debug, PartialEq, Eq, clap::ValueEnum)]
pub(crate) enum EthcdActivationFlag {
    /// Label EThcD/ETciD as EThcD (default; routes to the EThcD model).
    Ethcd,
    /// Label them ETD so the c/z scoring path engages.
    Etd,
}

/// Bundled species-specific N-glycan databases for `--glyco-species`.
#[derive(Copy, Clone, Debug, PartialEq, Eq, clap::ValueEnum)]
pub(crate) enum GlycoSpeciesFlag {
    /// Human N-glycans.
    Human,
    /// Mouse N-glycans.
    Mouse,
    /// Mouse, extended.
    MouseLarge,
    /// High-mannose N-glycans (species-independent).
    HighMannose,
}

impl GlycoSpeciesFlag {
    pub(crate) fn key(self) -> &'static str {
        match self {
            GlycoSpeciesFlag::Human => "human",
            GlycoSpeciesFlag::Mouse => "mouse",
            GlycoSpeciesFlag::MouseLarge => "mouse-large",
            GlycoSpeciesFlag::HighMannose => "high-mannose",
        }
    }
}

/// Glycan core type for `--glyco-core`: N- vs O-glycan search mode.
#[derive(Copy, Clone, Debug, PartialEq, Eq, clap::ValueEnum)]
pub(crate) enum GlycoCoreFlag {
    /// N-glycan: N-X-S/T sequon gate, trimannosyl core.
    N,
    /// O-glycan: Ser/Thr gate, single reducing-end GalNAc core.
    O,
}

impl From<GlycoCoreFlag> for glyco::glycan::first::GlycanCore {
    fn from(value: GlycoCoreFlag) -> Self {
        match value {
            GlycoCoreFlag::N => glyco::glycan::first::GlycanCore::NGlycan,
            GlycoCoreFlag::O => glyco::glycan::first::GlycanCore::OGlycan,
        }
    }
}

/// `--precursor-mono`: MS1 isotope-envelope precursor correction.
#[derive(Clone, Copy, Debug, PartialEq, Eq, ValueEnum)]
pub(crate) enum PrecursorMonoFlag {
    /// Leave every precursor m/z as recorded.
    Off,
    /// Move a precursor down to the monoisotope its MS1 envelope supports (default).
    Auto,
}

pub(crate) fn parse_enum<E: ValueEnum>(
    what: &str,
    s: &str,
    ignore_case: bool,
) -> Result<E, String> {
    <E as ValueEnum>::from_str(s, ignore_case).map_err(|_| {
        let names: Vec<String> = E::value_variants()
            .iter()
            .filter_map(|v| v.to_possible_value())
            .map(|p| p.get_name().to_string())
            .collect();
        format!("invalid {what} `{s}`: expected {}", names.join("|"))
    })
}

pub(crate) fn parse_decoy_strategy(s: &str) -> Result<DecoyStrategy, String> {
    DecoyStrategy::from_name(s).ok_or_else(|| {
        format!("unknown --decoy-strategy '{s}' (expected shuffle/reverse/sequon/mutate/none)")
    })
}

pub(crate) fn decoy_strategy_name(strategy: DecoyStrategy) -> &'static str {
    match strategy {
        DecoyStrategy::Shuffle => "shuffle",
        DecoyStrategy::Reverse => "reverse",
        DecoyStrategy::Sequon => "sequon",
        DecoyStrategy::Mutate => "mutate",
        DecoyStrategy::None => "none",
    }
}

pub(crate) fn parse_enzyme_list(s: &str) -> Result<String, String> {
    crate::models::parse_enzymes(s)
        .map(|_| s.to_string())
        .map_err(|e| e.to_string())
}

pub(crate) fn parse_plex(s: &str) -> Result<quant::Plex, String> {
    quant::Plex::parse(s)
}

pub(crate) fn parse_precursor_cal(s: &str) -> Result<PrecursorCalMode, String> {
    match s.to_ascii_lowercase().as_str() {
        "auto" => Ok(PrecursorCalMode::Auto),
        "on" => Ok(PrecursorCalMode::On),
        "off" => Ok(PrecursorCalMode::Off),
        _ => Err(format!("invalid precursor-cal `{s}`: expected off|auto|on")),
    }
}

pub(crate) fn parse_uniprot_evidence(s: &str) -> Result<input::EvidenceLevel, String> {
    input::EvidenceLevel::from_name(s).ok_or_else(|| {
        format!("invalid UniProt evidence level `{s}`: expected experimental|curated|any")
    })
}

pub(super) fn parse_transfer_models(s: &str) -> Result<usize, String> {
    match s.parse::<usize>() {
        Ok(n) if (1..=64).contains(&n) => Ok(n),
        _ => Err(format!("{s}: expected 1 to 64 transfer classifier members")),
    }
}
