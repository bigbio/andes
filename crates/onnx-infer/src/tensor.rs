//! Named dense `f32` tensors exchanged with ONNX models.

use crate::error::{Error, Result};

/// A named, row-major `f32` tensor.
///
/// Shape uses ONNX convention: leftmost axis is typically the batch dimension.
#[derive(Debug, Clone, PartialEq)]
pub struct NamedTensor {
    /// ONNX graph input/output name.
    pub name: String,
    /// Dimension sizes (row-major).
    pub shape: Vec<usize>,
    /// Contiguous `f32` payload; `len == shape.product()`.
    pub data: Vec<f32>,
}

impl NamedTensor {
    /// Build a tensor; `data.len()` must equal the product of `shape`.
    pub fn new(name: impl Into<String>, shape: Vec<usize>, data: Vec<f32>) -> Result<Self> {
        let name = name.into();
        let expected = shape.iter().copied().product::<usize>();
        if data.len() != expected {
            return Err(Error::ShapeMismatch {
                name,
                got: shape,
                got_nelem: data.len(),
                expected_nelem: expected,
            });
        }
        Ok(Self { name, shape, data })
    }

    /// Number of elements.
    pub fn len(&self) -> usize {
        self.data.len()
    }

    /// True when there are no elements.
    pub fn is_empty(&self) -> bool {
        self.data.is_empty()
    }

    /// Batch size when the leading axis is the batch (empty shape → 1).
    pub fn batch_size(&self) -> usize {
        self.shape.first().copied().unwrap_or(1)
    }

    /// Split a batched tensor into per-row views (`shape[0]` rows).
    ///
    /// Useful when a search engine scores top-N PSMs as one batch then wants
    /// per-PSM feature columns.
    pub fn rows(&self) -> Result<Vec<NamedTensor>> {
        if self.shape.is_empty() {
            return Ok(vec![self.clone()]);
        }
        let batch = self.shape[0];
        if batch == 0 {
            return Ok(Vec::new());
        }
        let row_nelem = self.data.len() / batch;
        let row_shape: Vec<usize> = if self.shape.len() == 1 {
            vec![1]
        } else {
            self.shape[1..].to_vec()
        };
        let mut out = Vec::with_capacity(batch);
        for i in 0..batch {
            let start = i * row_nelem;
            let end = start + row_nelem;
            out.push(NamedTensor::new(
                format!("{}[{i}]", self.name),
                row_shape.clone(),
                self.data[start..end].to_vec(),
            )?);
        }
        Ok(out)
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn rejects_mismatched_nelem() {
        let err = NamedTensor::new("x", vec![2, 3], vec![1.0; 5]).unwrap_err();
        assert!(matches!(err, Error::ShapeMismatch { .. }));
    }

    #[test]
    fn splits_batch_rows() {
        let t = NamedTensor::new("x", vec![2, 3], vec![1., 2., 3., 4., 5., 6.]).unwrap();
        let rows = t.rows().unwrap();
        assert_eq!(rows.len(), 2);
        assert_eq!(rows[0].data, vec![1., 2., 3.]);
        assert_eq!(rows[1].data, vec![4., 5., 6.]);
    }
}
