//! SplitMix64 and FNV-1a, the two small non-cryptographic hashes used to derive seeds, folds and
//! stable identifiers across the workspace.

use std::borrow::Borrow;

/// Golden-ratio increment of SplitMix64.
pub const SPLITMIX_GOLDEN: u64 = 0x9E37_79B9_7F4A_7C15;

/// FNV-1a 64-bit offset basis.
pub const FNV1A64_OFFSET: u64 = 0xCBF2_9CE4_8422_2325;

const FNV1A64_PRIME: u64 = 0x0000_0100_0000_01B3;

/// Stateless SplitMix64: the output for the generator state `x` before it advances.
#[inline]
pub fn splitmix64(x: u64) -> u64 {
    let mut z = x.wrapping_add(SPLITMIX_GOLDEN);
    z = (z ^ (z >> 30)).wrapping_mul(0xBF58_476D_1CE4_E5B9);
    z = (z ^ (z >> 27)).wrapping_mul(0x94D0_49BB_1331_11EB);
    z ^ (z >> 31)
}

/// Next output of the SplitMix64 stream held in `state`.
#[inline]
pub fn splitmix64_next(state: &mut u64) -> u64 {
    let z = splitmix64(*state);
    *state = state.wrapping_add(SPLITMIX_GOLDEN);
    z
}

/// Continues an FNV-1a 64-bit hash from `h` over `bytes`.
pub fn fnv1a64_extend(h: u64, bytes: impl IntoIterator<Item = impl Borrow<u8>>) -> u64 {
    bytes.into_iter().fold(h, |h, b| {
        (h ^ u64::from(*b.borrow())).wrapping_mul(FNV1A64_PRIME)
    })
}

/// FNV-1a 64-bit hash of `bytes`.
pub fn fnv1a64(bytes: impl IntoIterator<Item = impl Borrow<u8>>) -> u64 {
    fnv1a64_extend(FNV1A64_OFFSET, bytes)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn splitmix64_matches_the_reference_stream() {
        let mut state = 0u64;
        assert_eq!(splitmix64_next(&mut state), 0xE220_A839_7B1D_CDAF);
        assert_eq!(splitmix64(0), 0xE220_A839_7B1D_CDAF);
        assert_eq!(splitmix64_next(&mut state), splitmix64(SPLITMIX_GOLDEN));
    }

    #[test]
    fn fnv1a64_matches_the_reference_vectors() {
        assert_eq!(fnv1a64(&[] as &[u8]), FNV1A64_OFFSET);
        assert_eq!(fnv1a64(*b"a"), 0xAF63_DC4C_8601_EC8C);
        assert_eq!(fnv1a64(*b"foobar"), 0x8594_4171_F739_67E8);
        assert_eq!(
            fnv1a64_extend(fnv1a64(*b"foo"), *b"bar"),
            fnv1a64(*b"foobar")
        );
    }
}
