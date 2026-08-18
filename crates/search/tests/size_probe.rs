#[test]
fn probe_sizes() {
    use std::mem::size_of;
    println!("FullGlycoPsm      = {}", size_of::<search::glyco_search::FullGlycoPsm>());
    println!("GlycoPsmKey       = {}", size_of::<andes_glyco::glyco_psm::GlycoPsmKey>());
    println!("PsmMatch          = {}", size_of::<search::psm::PsmMatch>());
    println!("GlycoSpectrumResult = {}", size_of::<search::glyco_search::GlycoSpectrumResult>());
}
