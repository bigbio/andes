use crate::pin::RowMeta;

const INSERTION_THRESHOLD: usize = 16;

struct ParkMiller(u64);

impl ParkMiller {
    fn next(&mut self) -> u64 {
        self.0 = (self.0 * 279_470_273) % 4_294_967_291;
        self.0
    }
}

fn fast_uint_hash(x: u32) -> u32 {
    let x = ((x >> 16) ^ x).wrapping_mul(0x045d_9f3b);
    let x = ((x >> 16) ^ x).wrapping_mul(0x045d_9f3b);
    (x >> 16) ^ x
}

fn scan_hash(m: &RowMeta) -> u32 {
    fast_uint_hash(m.file) ^ fast_uint_hash(m.scan as u32).wrapping_shl(1)
}

pub fn libstdcxx_sort_by<T: Copy, F: Fn(&T, &T) -> bool>(v: &mut [T], less: F) {
    let n = v.len();
    if n < 2 {
        return;
    }
    let depth = 2 * (usize::BITS - 1 - n.leading_zeros()) as usize;
    introsort_loop(v, 0, n, depth, &less);
    final_insertion_sort(v, 0, n, &less);
}

fn introsort_loop<T: Copy, F: Fn(&T, &T) -> bool>(
    v: &mut [T],
    first: usize,
    mut last: usize,
    mut depth: usize,
    less: &F,
) {
    while last - first > INSERTION_THRESHOLD {
        if depth == 0 {
            make_heap(v, first, last, less);
            sort_heap(v, first, last, less);
            return;
        }
        depth -= 1;
        let cut = partition_pivot(v, first, last, less);
        introsort_loop(v, cut, last, depth, less);
        last = cut;
    }
}

fn partition_pivot<T: Copy, F: Fn(&T, &T) -> bool>(
    v: &mut [T],
    first: usize,
    last: usize,
    less: &F,
) -> usize {
    let mid = first + (last - first) / 2;
    move_median_to_first(v, first, first + 1, mid, last - 1, less);
    unguarded_partition(v, first + 1, last, first, less)
}

fn move_median_to_first<T: Copy, F: Fn(&T, &T) -> bool>(
    v: &mut [T],
    result: usize,
    a: usize,
    b: usize,
    c: usize,
    less: &F,
) {
    let pick = if less(&v[a], &v[b]) {
        if less(&v[b], &v[c]) {
            b
        } else if less(&v[a], &v[c]) {
            c
        } else {
            a
        }
    } else if less(&v[a], &v[c]) {
        a
    } else if less(&v[b], &v[c]) {
        c
    } else {
        b
    };
    v.swap(result, pick);
}

fn unguarded_partition<T: Copy, F: Fn(&T, &T) -> bool>(
    v: &mut [T],
    mut first: usize,
    mut last: usize,
    pivot: usize,
    less: &F,
) -> usize {
    loop {
        while less(&v[first], &v[pivot]) {
            first += 1;
        }
        last -= 1;
        while less(&v[pivot], &v[last]) {
            last -= 1;
        }
        if first >= last {
            return first;
        }
        v.swap(first, last);
        first += 1;
    }
}

fn final_insertion_sort<T: Copy, F: Fn(&T, &T) -> bool>(
    v: &mut [T],
    first: usize,
    last: usize,
    less: &F,
) {
    if last - first > INSERTION_THRESHOLD {
        insertion_sort(v, first, first + INSERTION_THRESHOLD, less);
        for i in first + INSERTION_THRESHOLD..last {
            unguarded_linear_insert(v, i, less);
        }
    } else {
        insertion_sort(v, first, last, less);
    }
}

fn insertion_sort<T: Copy, F: Fn(&T, &T) -> bool>(
    v: &mut [T],
    first: usize,
    last: usize,
    less: &F,
) {
    if first == last {
        return;
    }
    for i in first + 1..last {
        if less(&v[i], &v[first]) {
            let val = v[i];
            v.copy_within(first..i, first + 1);
            v[first] = val;
        } else {
            unguarded_linear_insert(v, i, less);
        }
    }
}

fn unguarded_linear_insert<T: Copy, F: Fn(&T, &T) -> bool>(v: &mut [T], mut last: usize, less: &F) {
    let val = v[last];
    while less(&val, &v[last - 1]) {
        v[last] = v[last - 1];
        last -= 1;
    }
    v[last] = val;
}

fn make_heap<T: Copy, F: Fn(&T, &T) -> bool>(v: &mut [T], first: usize, last: usize, less: &F) {
    let len = last - first;
    if len < 2 {
        return;
    }
    let mut parent = (len - 2) / 2;
    loop {
        let value = v[first + parent];
        adjust_heap(v, first, parent, len, value, less);
        if parent == 0 {
            return;
        }
        parent -= 1;
    }
}

fn adjust_heap<T: Copy, F: Fn(&T, &T) -> bool>(
    v: &mut [T],
    first: usize,
    mut hole: usize,
    len: usize,
    value: T,
    less: &F,
) {
    let top = hole;
    let mut second = hole;
    while second < (len - 1) / 2 {
        second = 2 * (second + 1);
        if less(&v[first + second], &v[first + second - 1]) {
            second -= 1;
        }
        v[first + hole] = v[first + second];
        hole = second;
    }
    if len & 1 == 0 && second == (len - 2) / 2 {
        second = 2 * (second + 1);
        v[first + hole] = v[first + second - 1];
        hole = second - 1;
    }
    while hole > top {
        let parent = (hole - 1) / 2;
        if !less(&v[first + parent], &value) {
            break;
        }
        v[first + hole] = v[first + parent];
        hole = parent;
    }
    v[first + hole] = value;
}

fn sort_heap<T: Copy, F: Fn(&T, &T) -> bool>(v: &mut [T], first: usize, mut last: usize, less: &F) {
    while last - first > 1 {
        last -= 1;
        let value = v[last];
        v[last] = v[first];
        adjust_heap(v, first, 0, last - first, value, less);
    }
}

pub fn folds_by_spectrum(meta: &[RowMeta], k: usize, seed: u64) -> Vec<u8> {
    let k = k.max(1);
    let n = meta.len();
    let mut fold = vec![0u8; n];
    if n == 0 {
        return fold;
    }
    let mut order: Vec<(u32, u32)> = (0..n)
        .filter(|&i| meta[i].label > 0)
        .chain((0..n).filter(|&i| meta[i].label <= 0))
        .map(|i| (scan_hash(&meta[i]), i as u32))
        .collect();
    libstdcxx_sort_by(&mut order, |a, b| a.0 < b.0);
    let mut remain = vec![0i64; k];
    let mut left = n as i64;
    for f in (0..k).rev() {
        remain[f] = left / (f as i64 + 1);
        left -= remain[f];
    }
    let mut rng = ParkMiller(seed);
    let scan = |i: u32| meta[i as usize].scan as u32;
    let mut previous = scan(order[0].1);
    let mut chosen = (rng.next() % k as u64) as usize;
    for &(_, i) in &order {
        let current = scan(i);
        if previous != current {
            chosen = (rng.next() % k as u64) as usize;
            while remain[chosen] <= 0 {
                chosen = (rng.next() % k as u64) as usize;
            }
        }
        fold[i as usize] = chosen as u8;
        remain[chosen] -= 1;
        previous = current;
    }
    fold
}
