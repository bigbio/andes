//! Centroiding of profile-mode spectra on the fly: local maxima with an
//! intensity-weighted m/z.

pub(super) fn log_profile_centroid_once() {
    use std::sync::atomic::{AtomicBool, Ordering};
    static WARNED: AtomicBool = AtomicBool::new(false);
    if WARNED
        .compare_exchange(false, true, Ordering::Relaxed, Ordering::Relaxed)
        .is_ok()
    {
        log::info!(
            "profile-mode spectra detected (MS:1000128) — centroiding peaks on the fly \
             (local-max + intensity-weighted m/z). For best results pre-centroid during \
             conversion (e.g. the converter's peak-picking option)."
        );
    }
}

pub(super) fn centroid_profile(profile: &[(f64, f32)]) -> Vec<(f64, f32)> {
    let n = profile.len();
    if n < 3 {
        return profile.to_vec();
    }
    let mut out: Vec<(f64, f32)> = Vec::new();
    let mut i = 1usize;
    while i < n - 1 {
        let (prev, cur, next) = (profile[i - 1].1, profile[i].1, profile[i + 1].1);
        if cur > 0.0 && cur >= prev && cur > next {
            let mut l = i;
            while l > 0 && profile[l - 1].1 <= profile[l].1 && profile[l - 1].1 > 0.0 {
                l -= 1;
            }
            let mut r = i;
            while r + 1 < n && profile[r + 1].1 <= profile[r].1 && profile[r + 1].1 > 0.0 {
                r += 1;
            }
            let mut wsum = 0f64;
            let mut msum = 0f64;
            let mut apex = 0f32;
            for &(m, it) in &profile[l..=r] {
                wsum += it as f64;
                msum += m * it as f64;
                if it > apex {
                    apex = it;
                }
            }
            if wsum > 0.0 {
                out.push((msum / wsum, apex));
            }
            i = r + 1;
        } else {
            i += 1;
        }
    }
    out
}
