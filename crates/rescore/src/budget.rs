use std::path::{Path, PathBuf};
use std::sync::atomic::{AtomicU64, Ordering};
use std::sync::Arc;

pub struct Budget {
    total: AtomicU64,
    remaining: AtomicU64,
    peak: AtomicU64,
    spill_dir: PathBuf,
}

impl Budget {
    pub fn new(bytes: u64, spill_dir: &Path) -> Arc<Budget> {
        Arc::new(Budget {
            total: AtomicU64::new(bytes),
            remaining: AtomicU64::new(bytes),
            peak: AtomicU64::new(0),
            spill_dir: spill_dir.to_path_buf(),
        })
    }

    pub fn total(&self) -> u64 {
        self.total.load(Ordering::Relaxed)
    }

    /// Lowers the budget to `bytes`, or to what is reserved when that is more. Call it while
    /// no other thread reserves from the budget.
    pub fn shrink_to(&self, bytes: u64) {
        let total = self.total();
        let remaining = self.available();
        let bytes = bytes.max(total - remaining);
        if bytes >= total {
            return;
        }
        self.total.store(bytes, Ordering::Relaxed);
        self.remaining
            .store(remaining - (total - bytes), Ordering::Release);
    }

    pub fn available(&self) -> u64 {
        self.remaining.load(Ordering::Relaxed)
    }

    /// The most bytes reserved at once so far.
    pub fn peak(&self) -> u64 {
        self.peak.load(Ordering::Relaxed)
    }

    fn note_reserved(&self, remaining: u64) {
        self.peak
            .fetch_max(self.total().saturating_sub(remaining), Ordering::Relaxed);
    }

    pub fn spill_dir(&self) -> &Path {
        &self.spill_dir
    }

    pub fn reserve_keeping(&self, bytes: u64, keep: u64) -> bool {
        match self
            .remaining
            .fetch_update(Ordering::AcqRel, Ordering::Acquire, |r| {
                r.checked_sub(bytes).filter(|&left| left >= keep)
            }) {
            Ok(before) => {
                self.note_reserved(before - bytes);
                true
            }
            Err(_) => false,
        }
    }

    pub fn reserve(&self, bytes: u64) -> bool {
        self.reserve_keeping(bytes, 0)
    }

    pub fn take_up_to(&self, bytes: u64) -> u64 {
        let mut got = 0;
        if let Ok(before) = self
            .remaining
            .fetch_update(Ordering::AcqRel, Ordering::Acquire, |r| {
                got = r.min(bytes);
                Some(r - got)
            })
        {
            self.note_reserved(before - got);
        }
        got
    }

    pub fn release(&self, bytes: u64) {
        self.remaining.fetch_add(bytes, Ordering::AcqRel);
    }
}

pub struct Reservation {
    budget: Option<Arc<Budget>>,
    bytes: u64,
}

impl Reservation {
    pub fn none() -> Reservation {
        Reservation {
            budget: None,
            bytes: 0,
        }
    }

    pub fn take(budget: &Arc<Budget>, bytes: u64) -> Option<Reservation> {
        budget.reserve(bytes).then(|| Reservation {
            budget: Some(budget.clone()),
            bytes,
        })
    }

    pub fn grow_up_to(&mut self, budget: &Arc<Budget>, bytes: u64) {
        let got = budget.take_up_to(bytes);
        self.budget = Some(budget.clone());
        self.bytes += got;
    }

    pub fn bytes(&self) -> u64 {
        self.bytes
    }

    pub fn grow(&mut self, budget: &Arc<Budget>, bytes: u64, keep: u64) -> bool {
        if !budget.reserve_keeping(bytes, keep) {
            return false;
        }
        self.budget = Some(budget.clone());
        self.bytes += bytes;
        true
    }
}

impl Drop for Reservation {
    fn drop(&mut self) {
        if let Some(b) = &self.budget {
            b.release(self.bytes);
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn shrinking_keeps_what_is_reserved_and_never_goes_below_it() {
        let b = Budget::new(100, Path::new("."));
        let held = Reservation::take(&b, 30).unwrap();
        b.shrink_to(50);
        assert_eq!((b.total(), b.available()), (50, 20));
        drop(held);
        assert_eq!(b.available(), 50);
        let held = Reservation::take(&b, 40).unwrap();
        b.shrink_to(10);
        assert_eq!((b.total(), b.available()), (40, 0));
        drop(held);
        assert_eq!(b.available(), 40);
        b.shrink_to(1000);
        assert_eq!(b.total(), 40);
        assert_eq!(b.peak(), 40);
    }
}
