//! Explicit runtime initialisation with a file that is not a loadable
//! library. A separate test binary: a failed load leaves `ort`'s runtime
//! unusable for the rest of the process.

#![cfg(feature = "ort")]

use onnx_infer::backend::ort_backend::init_runtime;

#[test]
fn init_runtime_with_a_non_library_file_is_an_error() {
    let dir = tempfile::tempdir().unwrap();
    let path = dir.path().join("libonnxruntime.so");
    std::fs::write(&path, "not a shared library\n").unwrap();
    let r = std::panic::catch_unwind(|| init_runtime(&path));
    let e = r.expect("init_runtime panicked").unwrap_err();
    assert!(e.to_string().contains("libonnxruntime.so"), "{e}");
}
