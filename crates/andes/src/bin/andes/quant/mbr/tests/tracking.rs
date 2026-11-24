use std::cell::Cell;

thread_local! {
    static LIVE: Cell<usize> = const { Cell::new(0) };
    static PEAK: Cell<usize> = const { Cell::new(0) };
    static LOADS: Cell<usize> = const { Cell::new(0) };
}

pub(super) fn acquire() {
    let live = LIVE.get() + 1;
    LIVE.set(live);
    PEAK.set(PEAK.get().max(live));
}

pub(super) fn release() {
    LIVE.set(LIVE.get() - 1);
}

pub(super) fn loaded() {
    LOADS.set(LOADS.get() + 1);
}

pub(super) fn reset() {
    PEAK.set(LIVE.get());
    LOADS.set(0);
}

pub(super) fn stats() -> (usize, usize) {
    (PEAK.get(), LOADS.get())
}
