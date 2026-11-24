//! PSI-MS and UO controlled-vocabulary accessions the mzML reader recognises, and the
//! encodings it refuses rather than decoding into wrong peaks.

pub(super) const CV_ANALYZER_ION_TRAP: &str = "MS:1000264";
pub(super) const CV_ANALYZER_QUAD_ION_TRAP: &str = "MS:1000082";
pub(super) const CV_ANALYZER_RADIAL_LIT: &str = "MS:1000083";
pub(super) const CV_ANALYZER_LINEAR_ION_TRAP: &str = "MS:1000291";
pub(super) const CV_ANALYZER_ORBITRAP: &str = "MS:1000484";
pub(super) const CV_ANALYZER_FTICR: &str = "MS:1000079";
pub(super) const CV_ANALYZER_TOF: &str = "MS:1000084";

pub(super) const CV_MODEL_Q_EXACTIVE: &str = "MS:1001911";
pub(super) const CV_MODEL_Q_EXACTIVE_HF: &str = "MS:1002523";
pub(super) const CV_MODEL_Q_EXACTIVE_HF_X: &str = "MS:1002634";
pub(super) const CV_MODEL_Q_EXACTIVE_PLUS: &str = "MS:1002877";
pub(super) const CV_MODEL_ORBITRAP_FUSION: &str = "MS:1002416";
pub(super) const CV_MODEL_ORBITRAP_ASTRAL: &str = "MS:1003378";

pub(super) const CV_MS_LEVEL: &str = "MS:1000511";
pub(super) const CV_SCAN_TIME: &str = "MS:1000016";
pub(super) const CV_PROFILE: &str = "MS:1000128";
pub(super) const CV_SELECTED_ION_MZ: &str = "MS:1000744";
pub(super) const CV_ISOLATION_LOWER_OFFSET: &str = "MS:1000828";
pub(super) const CV_ISOLATION_TARGET_MZ: &str = "MS:1000827";
pub(super) const CV_ISOLATION_UPPER_OFFSET: &str = "MS:1000829";
pub(super) const CV_SCAN_WINDOW_LOWER: &str = "MS:1000501";
pub(super) const CV_SCAN_WINDOW_UPPER: &str = "MS:1000500";
pub(super) const CV_FAIMS_CV: &str = "MS:1001581";
pub(super) const CV_FILTER_STRING: &str = "MS:1000512";
pub(super) const CV_MZ_PLAIN: &str = "MS:1000040";
pub(super) const CV_CHARGE_STATE: &str = "MS:1000041";
pub(super) const CV_PEAK_INTENSITY: &str = "MS:1000042";
pub(super) const CV_MZ_ARRAY: &str = "MS:1000514";
pub(super) const CV_INTENSITY_ARRAY: &str = "MS:1000515";
pub(super) const CV_NOISE_MZ_ARRAY: &str = "MS:1002743";
pub(super) const CV_NOISE_ARRAY: &str = "MS:1002744";
pub(super) const CV_64BIT: &str = "MS:1000523";
pub(super) const CV_32BIT: &str = "MS:1000521";
pub(super) const CV_ZLIB: &str = "MS:1000574";
pub(super) const CV_TRUNCATION_ZLIB: &str = "MS:1003088";

pub(super) const UNSUPPORTED_COMPRESSION: &[(&str, &str)] = &[
    ("MS:1002312", "MS-Numpress linear prediction"),
    ("MS:1002313", "MS-Numpress positive integer"),
    ("MS:1002314", "MS-Numpress short logged float"),
    ("MS:1002746", "MS-Numpress linear prediction + zlib"),
    ("MS:1002747", "MS-Numpress positive integer + zlib"),
    ("MS:1002748", "MS-Numpress short logged float + zlib"),
    ("MS:1003089", "truncation, delta prediction + zlib"),
    ("MS:1003090", "truncation, linear prediction + zlib"),
    ("MS:1003780", "zstd"),
    ("MS:1003781", "byte-shuffled zstd"),
    ("MS:1003782", "dictionary-encoded zstd"),
    ("MS:1003783", "MS-Numpress linear prediction + zstd"),
    ("MS:1003784", "MS-Numpress positive integer + zstd"),
    ("MS:1003785", "MS-Numpress short logged float + zstd"),
    ("MS:1003826", "coordinate grid encoding"),
];

pub(super) const CV_CID: &str = "MS:1000133";
pub(super) const CV_HCD: &str = "MS:1000422";
pub(super) const CV_ETD: &str = "MS:1000598";
pub(super) const CV_PQD: &str = "MS:1000599";
pub(super) const CV_UVPD: &str = "MS:1000435";
pub(super) const CV_ECD: &str = "MS:1000250";

pub(super) const CV_TRAP_CID: &str = "MS:1002472";
pub(super) const CV_HIGH_ENERGY_BEAM_CID: &str = "MS:1002481";
pub(super) const CV_UVPD_CHILD: &str = "MS:1003246";

pub(super) const CV_SUPPLEMENTAL_HCD: &str = "MS:1002678";
pub(super) const CV_SUPPLEMENTAL_CID: &str = "MS:1002679";

pub(super) const CV_ETHCD: &str = "MS:1002631";
pub(super) const CV_ETCID: &str = "MS:1003182";

const CV_UNIT_MINUTE: &str = "UO:0000031";
const CV_UNIT_MINUTE_LEGACY: &str = "MS:1000038";

pub(super) fn is_minute_unit(unit_accession: &str, unit_name: &str) -> bool {
    match unit_accession {
        CV_UNIT_MINUTE | CV_UNIT_MINUTE_LEGACY => true,
        "" => matches!(
            unit_name.trim().to_ascii_lowercase().as_str(),
            "minute" | "min" | "minutes"
        ),
        _ => false,
    }
}
