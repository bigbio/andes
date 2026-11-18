//! Which MS1 scan an MS2 scan takes its precursor from.
//!
//! An MS2 scan links to the latest MS1 scan acquired at its own FAIMS compensation voltage,
//! else to the latest MS1 scan. With several compensation voltages the MS1 scan just before
//! an MS2 scan is often at another voltage and does not hold its precursor.

/// The latest MS1 scan overall and per compensation voltage, as the caller's MS1 indices.
#[derive(Debug, Default, Clone)]
pub struct Ms1Tracker {
    latest: Option<usize>,
    by_cv: Vec<(u32, usize)>,
}

impl Ms1Tracker {
    /// MS1 scan `idx` was acquired, at compensation voltage `cv` if any.
    pub fn push(&mut self, idx: usize, cv: Option<f32>) {
        self.latest = Some(idx);
        if let Some(cv) = cv {
            let key = cv.to_bits();
            match self.by_cv.iter_mut().find(|(k, _)| *k == key) {
                Some(e) => e.1 = idx,
                None => self.by_cv.push((key, idx)),
            }
        }
    }

    /// The MS1 scan an MS2 scan at compensation voltage `cv` takes its precursor from.
    pub fn link(&self, cv: Option<f32>) -> Option<usize> {
        cv.and_then(|cv| {
            let key = cv.to_bits();
            self.by_cv.iter().find(|(k, _)| *k == key).map(|e| e.1)
        })
        .or(self.latest)
    }

    /// The MS1 indices later MS2 scans can still link to, ascending.
    pub fn live(&self) -> Vec<usize> {
        let mut v: Vec<usize> = self.latest.into_iter().collect();
        v.extend(self.by_cv.iter().map(|e| e.1));
        v.sort_unstable();
        v.dedup();
        v
    }

    /// Renumbers the tracked MS1 scans after the caller kept only `live()`, in order.
    pub fn renumber_to_live(&mut self) {
        let live = self.live();
        let at = |old: usize| live.binary_search(&old).expect("tracked index is live");
        self.latest = self.latest.map(at);
        for e in self.by_cv.iter_mut() {
            e.1 = at(e.1);
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn an_ms2_scan_links_to_the_latest_ms1_at_its_compensation_voltage() {
        let mut t = Ms1Tracker::default();
        assert_eq!(t.link(Some(-50.0)), None);
        t.push(0, Some(-50.0));
        t.push(1, Some(-60.0));
        t.push(2, Some(-75.0));
        assert_eq!(t.link(Some(-60.0)), Some(1));
        assert_eq!(t.link(Some(-50.0)), Some(0));
        // No MS1 at that voltage yet, or no voltage: the latest MS1.
        assert_eq!(t.link(Some(-45.0)), Some(2));
        assert_eq!(t.link(None), Some(2));
        t.push(3, Some(-60.0));
        assert_eq!(t.link(Some(-60.0)), Some(3));
        assert_eq!(t.live(), vec![0, 2, 3]);
        t.renumber_to_live();
        assert_eq!(
            (
                t.link(Some(-50.0)),
                t.link(Some(-75.0)),
                t.link(Some(-60.0))
            ),
            (Some(0), Some(1), Some(2))
        );
    }

    #[test]
    fn without_compensation_voltages_the_latest_ms1_is_linked() {
        let mut t = Ms1Tracker::default();
        t.push(0, None);
        t.push(1, None);
        assert_eq!(t.link(None), Some(1));
        assert_eq!(t.live(), vec![1]);
        t.renumber_to_live();
        assert_eq!(t.link(None), Some(0));
    }
}
