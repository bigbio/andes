#[test]
fn probe_sizes() {
    use std::mem::size_of;
    println!(
        "FullGlycoPsm      = {}",
        size_of::<search::glycopeptides::FullGlycoPsm>()
    );
    println!(
        "GlycoPsmKey       = {}",
        size_of::<glyco::glyco_psm::GlycoPsmKey>()
    );
    println!("PsmMatch          = {}", size_of::<search::psm::PsmMatch>());
    println!(
        "GlycoSpectrumResult = {}",
        size_of::<search::glycopeptides::GlycoSpectrumResult>()
    );
}
