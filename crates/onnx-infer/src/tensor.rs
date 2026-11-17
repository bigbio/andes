use crate::error::{Error, Result};

#[derive(Debug, Clone, PartialEq)]
pub struct NamedTensor {
    pub name: String,
    pub shape: Vec<usize>,
    pub data: Vec<f32>,
}

impl NamedTensor {
    pub fn new(name: impl Into<String>, shape: Vec<usize>, data: Vec<f32>) -> Result<Self> {
        let name = name.into();
        let expected = element_count(&shape).ok_or_else(|| {
            Error::InvalidArg(format!("tensor `{name}`: shape {shape:?} overflows usize"))
        })?;
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

    pub fn len(&self) -> usize {
        self.data.len()
    }

    pub fn is_empty(&self) -> bool {
        self.data.is_empty()
    }

    pub fn batch_size(&self) -> usize {
        self.shape.first().copied().unwrap_or(1)
    }

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

pub fn element_count(dims: &[usize]) -> Option<usize> {
    dims.iter().try_fold(1usize, |acc, &d| acc.checked_mul(d))
}

pub(crate) fn checked_element_count(dims: &[usize], what: &str) -> Result<usize> {
    element_count(dims)
        .ok_or_else(|| Error::InvalidArg(format!("{what}: shape {dims:?} overflows usize")))
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn rejects_overflowing_shape() {
        let err = NamedTensor::new("x", vec![usize::MAX, 2], Vec::new()).unwrap_err();
        assert!(matches!(err, Error::InvalidArg(_)), "{err:?}");
        assert_eq!(element_count(&[]), Some(1));
        assert_eq!(element_count(&[3, 0, usize::MAX]), Some(0));
        assert_eq!(element_count(&[usize::MAX, 2]), None);
    }

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
