//! Arrow schema definitions for the Parquet model store.
//!
//! The store uses a single Parquet file with a top-level `record_kind` discriminator
//! column (values `"manifest"` | `"table"`) so manifest rows and bulk-table rows share
//! one file/schema.  Adding nullable columns for future `sources`/`stats` kinds is safe
//! because Arrow schemas treat new nullable fields as backward-compatible extensions.

use std::sync::Arc;

use arrow::datatypes::{DataType, Field, Fields, Schema, SchemaRef};

// ── helpers ─────────────────────────────────────────────────────────────────

/// Shorthand: nullable field.
#[inline]
fn nf(name: &str, dt: DataType) -> Field {
    Field::new(name, dt, true)
}

/// Shorthand: non-null field.
#[inline]
fn rf(name: &str, dt: DataType) -> Field {
    Field::new(name, dt, false)
}

/// Build `List<item: <inner_dt> (nullable)>` (nullable list).
#[inline]
fn list_of(inner_dt: DataType) -> DataType {
    DataType::List(Arc::new(Field::new("item", inner_dt, true)))
}

/// Build `List<item: Struct<fields…> (non-null)>` (nullable list of structs).
#[inline]
fn list_of_struct(struct_fields: impl IntoIterator<Item = Field>) -> DataType {
    let fields: Fields = struct_fields.into_iter().collect();
    DataType::List(Arc::new(Field::new(
        "item",
        DataType::Struct(fields),
        false, // struct entries themselves are non-null
    )))
}

// ── public schemas ───────────────────────────────────────────────────────────

/// Manifest schema — one row per trained scoring model.
///
/// `model_id` is the only non-null column; all others are nullable so that
/// future schema extensions only need to add new nullable fields.
pub fn models_schema() -> SchemaRef {
    // `charge_hist`: List<Struct<charge: Int32, count: Int32>>
    let charge_hist_dt = list_of_struct([
        Field::new("charge", DataType::Int32, false),
        Field::new("count", DataType::Int32, false),
    ]);

    Arc::new(Schema::new(vec![
        // --- discriminator (shared with tables schema) ---
        nf("record_kind", DataType::Utf8),
        // --- identity ---
        rf("model_id", DataType::Utf8),
        // --- model selector keys ---
        nf("activation", DataType::Utf8),
        nf("instrument", DataType::Utf8),
        nf("enzyme", DataType::Utf8),
        nf("experiment_class", DataType::Utf8), // sorted slug set, e.g. "phospho+tmt"
        // --- scalar model parameters ---
        nf("version", DataType::Int32),
        nf("mme_val", DataType::Float32),
        nf("mme_is_ppm", DataType::Boolean),
        nf("apply_deconv", DataType::Boolean),
        nf("deconv_tol", DataType::Float32),
        nf("num_segments", DataType::Int32),
        nf("max_rank", DataType::Int32),
        nf("error_scaling_factor", DataType::Int32),
        nf("min_charge", DataType::Int32),
        nf("max_charge", DataType::Int32),
        // --- charge histogram ---
        nf("charge_hist", charge_hist_dt),
        // --- provenance ---
        nf("dataset", DataType::Utf8),
        nf("n_psms", DataType::Int64),
        nf("date", DataType::Utf8),
        nf("seed_model", DataType::Utf8),
    ]))
}

/// Bulk-tables schema — one row per (model_id, partition, ion, table_kind).
///
/// `values` is populated for distribution kinds; `offsets` for offset kinds.
/// Both are nullable so that each row carries only the relevant payload.
pub fn tables_schema() -> SchemaRef {
    // `offsets`: List<Struct<offset: Float32, freq: Float32>>
    let offsets_dt = list_of_struct([
        Field::new("offset", DataType::Float32, false),
        Field::new("freq", DataType::Float32, false),
    ]);

    Arc::new(Schema::new(vec![
        // --- discriminator ---
        nf("record_kind", DataType::Utf8),
        // --- identity ---
        rf("model_id", DataType::Utf8),
        // --- partition axes ---
        nf("part_charge", DataType::Int32),
        nf("part_mass_bits", DataType::Int32), // f32::to_bits() as i32 (matches combined_schema)
        nf("part_seg", DataType::Int32),
        // --- ion descriptor ---
        nf("ion_kind", DataType::Utf8),   // "prefix" | "suffix" | "noise" | "-"
        nf("ion_charge", DataType::Int32),
        // --- table kind ---
        nf("table_kind", DataType::Utf8), // "rank_dist"|"ion_err"|"noise_err"|
                                          // "ion_existence"|"frag_off"|"precursor_off"
        // --- payload (one of the two is non-null per row) ---
        nf("values", list_of(DataType::Float32)),
        nf("offsets", offsets_dt),
    ]))
}

/// Combined schema for a single-file store: manifest rows and table rows
/// share one Parquet file.
///
/// Layout:
/// - `record_kind` (non-null) — `"manifest"` | `"table"` | `"source"` | `"stat"`
/// - `model_id` (non-null)
/// - every manifest-only column (nullable) — table rows leave these null
/// - every table-only column (nullable) — manifest rows leave these null
/// - every source-only column (nullable) — other rows leave these null
/// - `source_id` (nullable) — shared by `"source"` and `"stat"` rows
///
/// # Source ledger (`"source"` rows)
///
/// One row per `(model_id, source_id)` that was used during training.
/// Columns: `source_id`, `dataset`, `n_psms`, `date`, `weight`,
/// `train_fdr`, `src_instrument`, `src_experiment_class`.
///
/// # Per-source sufficient statistics (`"stat"` rows)
///
/// One row per `(model_id, source_id, partition, ion_kind, ion_charge, table_kind)`.
/// Reuses `part_charge`, `part_mass_bits`, `part_seg`, `ion_kind`,
/// `ion_charge`, `table_kind` from the table column set.
///
/// `table_kind` values for `"stat"` rows:
/// - `"rank"` — keyed by `(Partition, IonType)`; `ion_kind`/`ion_charge`/
///   `ion_offset_bits` carry the IonType; `counts` carries `Vec<u64>` as
///   `List<Int64>`.
/// - `"error"` / `"noise_error"` — keyed by `Partition`; `ion_kind = "-"`,
///   `ion_charge = 0`, `ion_offset_bits = 0`; `counts = Vec<u64>`.
/// - `"existence"` — keyed by `Partition`; one row per partition; the
///   existence map `{(Partition, u32) → u64}` is serialised as a flat
///   `counts` list where `counts[idx]` = the count for existence index
///   `idx`.  The maximum `u32` idx seen determines the list length.
/// - `"charge"` — one row with `part_charge = 0`, `part_mass_bits = 0`,
///   `part_seg = 0`, `ion_kind = "-"`, `ion_charge = 0`,
///   `ion_offset_bits = 0`.  The charge map `{i32 → u64}` is serialised
///   as two parallel `counts` + `charge_keys: List<Int32>` (same row).
///   `counts[j]` corresponds to `charge_keys[j]`.
pub fn combined_schema() -> SchemaRef {
    let charge_hist_dt = list_of_struct([
        Field::new("charge", DataType::Int32, false),
        Field::new("count", DataType::Int32, false),
    ]);

    // Precursor offsets carry extra fields beyond plain (offset, freq).
    let precursor_off_dt = list_of_struct([
        Field::new("reduced_charge", DataType::Int32, false),
        Field::new("offset", DataType::Float32, false),
        Field::new("tol_is_ppm", DataType::Boolean, false),
        Field::new("tol_val", DataType::Float32, false),
        Field::new("frequency", DataType::Float32, false),
    ]);

    Arc::new(Schema::new(vec![
        // ── shared ──────────────────────────────────────────────────────────
        rf("record_kind", DataType::Utf8),
        rf("model_id", DataType::Utf8),
        // ── manifest-only ───────────────────────────────────────────────────
        nf("activation", DataType::Utf8),
        nf("instrument", DataType::Utf8),
        nf("enzyme", DataType::Utf8),
        nf("protocol", DataType::Utf8),
        nf("version", DataType::Int32),
        nf("mme_val", DataType::Float32),
        nf("mme_is_ppm", DataType::Boolean),
        nf("apply_deconv", DataType::Boolean),
        nf("deconv_tol", DataType::Float32),
        nf("num_segments", DataType::Int32),
        nf("max_rank", DataType::Int32),
        nf("error_scaling_factor", DataType::Int32),
        nf("min_charge", DataType::Int32),
        nf("max_charge", DataType::Int32),
        nf("num_precursor_off", DataType::Int32),
        nf("charge_hist", charge_hist_dt),
        // ── table-only ──────────────────────────────────────────────────────
        nf("part_charge", DataType::Int32),
        nf("part_mass_bits", DataType::Int32), // f32::to_bits() as i32 for bit-exact round-trip
        nf("part_seg", DataType::Int32),
        nf("ion_kind", DataType::Utf8),        // "prefix"|"suffix"|"noise"|"-"
        nf("ion_charge", DataType::Int32),
        nf("ion_offset_bits", DataType::Int32), // f32::to_bits() as i32; 0 for noise/dist rows
        nf("table_kind", DataType::Utf8),
        // "rank_dist", "ion_err", "noise_err", "ion_existence", "frag_off" → values
        nf("values", list_of(DataType::Float32)),
        // "precursor_off" → precursor_offsets (full struct)
        nf("precursor_offsets", precursor_off_dt),
        // ── source-only ─────────────────────────────────────────────────────
        nf("source_id", DataType::Utf8),
        nf("dataset", DataType::Utf8),
        nf("n_psms", DataType::Int64),
        nf("date", DataType::Utf8),
        nf("weight", DataType::Float32),
        nf("train_fdr", DataType::Float32),
        nf("src_instrument", DataType::Utf8),
        nf("src_experiment_class", DataType::Utf8),
        // ── stat-only ────────────────────────────────────────────────────────
        // "stat" rows reuse part_charge/part_mass_bits/part_seg/ion_kind/
        // ion_charge/ion_offset_bits/table_kind from the table section.
        // They carry counts as List<Int64> (u64 histogram bins).
        // For "charge" table_kind, charge_keys carries the i32 charge values
        // parallel to counts.
        nf("counts", list_of(DataType::Int64)),
        nf("charge_keys", list_of(DataType::Int32)),
    ]))
}

// ── tests ────────────────────────────────────────────────────────────────────

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn schemas_have_expected_columns() {
        let m = models_schema();
        assert!(m.field_with_name("model_id").is_ok());
        assert!(m.field_with_name("experiment_class").is_ok());
        assert!(m.field_with_name("charge_hist").is_ok());
        let t = tables_schema();
        assert!(t.field_with_name("table_kind").is_ok());
        assert!(t.field_with_name("values").is_ok());
        assert!(t.field_with_name("offsets").is_ok());
    }

    #[test]
    fn model_id_is_non_null() {
        for schema in [models_schema(), tables_schema()] {
            let f = schema.field_with_name("model_id").unwrap();
            assert!(!f.is_nullable(), "model_id must be non-null in {:?}", schema);
        }
    }

    #[test]
    fn charge_hist_is_list_of_struct() {
        let m = models_schema();
        let f = m.field_with_name("charge_hist").unwrap();
        match f.data_type() {
            DataType::List(inner) => {
                assert!(
                    matches!(inner.data_type(), DataType::Struct(_)),
                    "charge_hist inner must be Struct"
                );
                if let DataType::Struct(fields) = inner.data_type() {
                    assert!(fields.find("charge").is_some());
                    assert!(fields.find("count").is_some());
                }
            }
            other => panic!("charge_hist must be List, got {:?}", other),
        }
    }

    #[test]
    fn offsets_is_list_of_struct() {
        let t = tables_schema();
        let f = t.field_with_name("offsets").unwrap();
        match f.data_type() {
            DataType::List(inner) => {
                assert!(
                    matches!(inner.data_type(), DataType::Struct(_)),
                    "offsets inner must be Struct"
                );
                if let DataType::Struct(fields) = inner.data_type() {
                    assert!(fields.find("offset").is_some());
                    assert!(fields.find("freq").is_some());
                }
            }
            other => panic!("offsets must be List, got {:?}", other),
        }
    }

    #[test]
    fn values_is_list_of_float32() {
        let t = tables_schema();
        let f = t.field_with_name("values").unwrap();
        match f.data_type() {
            DataType::List(inner) => {
                assert_eq!(inner.data_type(), &DataType::Float32);
            }
            other => panic!("values must be List<Float32>, got {:?}", other),
        }
    }
}
