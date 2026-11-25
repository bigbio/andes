use std::sync::Arc;

use arrow::datatypes::{DataType, Field, Fields, Schema, SchemaRef};

#[inline]
fn nf(name: &str, dt: DataType) -> Field {
    Field::new(name, dt, true)
}

#[inline]
fn rf(name: &str, dt: DataType) -> Field {
    Field::new(name, dt, false)
}

#[inline]
fn list_of(inner_dt: DataType) -> DataType {
    DataType::List(Arc::new(Field::new("item", inner_dt, true)))
}

#[inline]
fn list_of_struct(struct_fields: impl IntoIterator<Item = Field>) -> DataType {
    let fields: Fields = struct_fields.into_iter().collect();
    DataType::List(Arc::new(Field::new(
        "item",
        DataType::Struct(fields),
        false,
    )))
}

pub fn models_schema() -> SchemaRef {
    let charge_hist_dt = list_of_struct([
        Field::new("charge", DataType::Int32, false),
        Field::new("count", DataType::Int32, false),
    ]);

    Arc::new(Schema::new(vec![
        nf("record_kind", DataType::Utf8),
        rf("model_id", DataType::Utf8),
        nf("activation", DataType::Utf8),
        nf("instrument", DataType::Utf8),
        nf("enzyme", DataType::Utf8),
        nf("experiment_class", DataType::Utf8),
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
        nf("charge_hist", charge_hist_dt),
        nf("dataset", DataType::Utf8),
        nf("n_psms", DataType::Int64),
        nf("date", DataType::Utf8),
        nf("seed_model", DataType::Utf8),
    ]))
}

pub fn tables_schema() -> SchemaRef {
    let offsets_dt = list_of_struct([
        Field::new("offset", DataType::Float32, false),
        Field::new("freq", DataType::Float32, false),
    ]);

    Arc::new(Schema::new(vec![
        nf("record_kind", DataType::Utf8),
        rf("model_id", DataType::Utf8),
        nf("part_charge", DataType::Int32),
        nf("part_mass_bits", DataType::Int32),
        nf("part_seg", DataType::Int32),
        nf("ion_kind", DataType::Utf8),
        nf("ion_charge", DataType::Int32),
        nf("table_kind", DataType::Utf8),
        nf("values", list_of(DataType::Float32)),
        nf("offsets", offsets_dt),
    ]))
}

pub fn combined_schema() -> SchemaRef {
    let charge_hist_dt = list_of_struct([
        Field::new("charge", DataType::Int32, false),
        Field::new("count", DataType::Int32, false),
    ]);

    let precursor_off_dt = list_of_struct([
        Field::new("reduced_charge", DataType::Int32, false),
        Field::new("offset", DataType::Float32, false),
        Field::new("tol_is_ppm", DataType::Boolean, false),
        Field::new("tol_val", DataType::Float32, false),
        Field::new("frequency", DataType::Float32, false),
    ]);

    Arc::new(Schema::new(vec![
        rf("record_kind", DataType::Utf8),
        rf("model_id", DataType::Utf8),
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
        nf("gbdt_model_bytes", DataType::Binary),
        nf("frag_intensity_model_bytes", DataType::Binary),
        nf("rich_ion_model_bytes", DataType::Binary),
        nf("part_charge", DataType::Int32),
        nf("part_mass_bits", DataType::Int32),
        nf("part_seg", DataType::Int32),
        nf("ion_kind", DataType::Utf8),
        nf("ion_charge", DataType::Int32),
        nf("ion_offset_bits", DataType::Int32),
        nf("ion_loss_class", DataType::Int32),
        nf("table_kind", DataType::Utf8),
        nf("values", list_of(DataType::Float32)),
        nf("precursor_offsets", precursor_off_dt),
        nf("frag_off_loss_classes", list_of(DataType::Int32)),
        nf("source_id", DataType::Utf8),
        nf("dataset", DataType::Utf8),
        nf("n_psms", DataType::Int64),
        nf("date", DataType::Utf8),
        nf("weight", DataType::Float32),
        nf("train_fdr", DataType::Float32),
        nf("decay_half_life_days", DataType::Float32),
        nf("src_instrument", DataType::Utf8),
        nf("src_experiment_class", DataType::Utf8),
        nf("counts", list_of(DataType::Int64)),
        nf("charge_keys", list_of(DataType::Int32)),
    ]))
}

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
            assert!(
                !f.is_nullable(),
                "model_id must be non-null in {schema:?}",
            );
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
            other => panic!("charge_hist must be List, got {other:?}"),
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
            other => panic!("offsets must be List, got {other:?}"),
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
            other => panic!("values must be List<Float32>, got {other:?}"),
        }
    }
}
