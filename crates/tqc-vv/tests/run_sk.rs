//! The Solovay–Kitaev density decision at the Atlas instance, asserted in both directions:
//! the probe (which errors on the refuted side by contract) must return the refutation,
//! and the decision witness must pin every exact value. This test can fail — a probe that
//! unexpectedly certifies density, or any drift in the pinned decision, is an error.
#![allow(clippy::unwrap_used)]

use tqc_core::UseCaseParams;
use tqc_vv::witness;

#[test]
fn run_sk_atlas() {
    let p = UseCaseParams::new(4, 3, 8); // atlas

    // The probe reports the refuted side as an error naming the refutation.
    match witness::solovay_kitaev_probe(&p) {
        Ok(m) => panic!(
            "single-qubit density unexpectedly certified (is_dense={}): {}",
            m.is_dense, m.description
        ),
        Err(e) => assert!(
            e.contains("refutes"),
            "probe error must state the exact refutation, got: {e}"
        ),
    }

    // The decision witness pins the exact decided values in both directions.
    witness::solovay_kitaev_decision_witness(&p).unwrap();
}
