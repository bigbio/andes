fn main() {
    let a: Vec<String> = std::env::args().collect();
    let n: usize = a.get(2).and_then(|v| v.parse().ok()).unwrap_or(1);
    let limit: usize = a.get(3).and_then(|v| v.parse().ok()).unwrap_or(3000);
    let reader = input::TimsTofReader::open(&a[1])
        .expect("open")
        .with_dia_mobility_slices(n);
    let total = reader.len();
    let t0 = std::time::Instant::now();
    let mut counts: Vec<usize> = Vec::new();
    for s in reader.take(limit) {
        counts.push(s.expect("spectrum").peaks.len());
    }
    counts.sort_unstable();
    let m = counts.len();
    let mean = counts.iter().sum::<usize>() as f64 / m.max(1) as f64;
    println!(
        "slices={n} spectra={total} read={m} in {:.1}s peaks/spectrum mean={mean:.0} median={} p90={} max={}",
        t0.elapsed().as_secs_f64(),
        counts[m / 2],
        counts[m * 9 / 10],
        counts[m - 1]
    );
}
