use std::path::{Path, PathBuf};
use std::sync::atomic::{AtomicU64, Ordering};
use std::sync::Arc;

pub struct Budget {
    total: u64,
    remaining: AtomicU64,
    spill_dir: PathBuf,
}

impl Budget {
    pub fn new(bytes: u64, spill_dir: &Path) -> Arc<Budget> {
        Arc::new(Budget {
            total: bytes,
            remaining: AtomicU64::new(bytes),
            spill_dir: spill_dir.to_path_buf(),
        })
    }

    pub fn total(&self) -> u64 {
        self.total
    }

    pub fn available(&self) -> u64 {
        self.remaining.load(Ordering::Relaxed)
    }

    pub fn spill_dir(&self) -> &Path {
        &self.spill_dir
    }

    pub fn reserve_keeping(&self, bytes: u64, keep: u64) -> bool {
        self.remaining
            .fetch_update(Ordering::AcqRel, Ordering::Acquire, |r| {
                r.checked_sub(bytes).filter(|&left| left >= keep)
            })
            .is_ok()
    }

    pub fn reserve(&self, bytes: u64) -> bool {
        self.reserve_keeping(bytes, 0)
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
