//! The open probes (universality, advantage), run over `features/targets/s4_open/`.
//!
//! These are MEASUREMENTS, not assertions: each scenario runs a probe, prints the measured
//! number, and asserts only that a valid measurement was produced — never that the open claim
//! holds. They are non-gating in that sense (they can never go green by asserting universality
//! or advantage).

#![allow(
    clippy::unwrap_used,
    clippy::expect_used,
    clippy::panic,
    clippy::missing_panics_doc
)]

use cucumber::{then, World};
use tqc_atlas::canonical;
use tqc_model::Model;
use tqc_vv::witness;

#[derive(Debug, Default, cucumber::World)]
struct ProbeWorld;

fn atlas() -> tqc_core::UseCaseParams {
    canonical(&Model::load().unwrap()).unwrap()
}

#[then("the generated-subgroup density is measured and universality remains open and unasserted")]
async fn universality(_w: &mut ProbeWorld) {
    let result = witness::universality_probe(&atlas()).unwrap();
    println!("[probe] universality: {result}");
    assert!(
        !result.is_empty(),
        "a measurement or obstruction must be produced"
    );
}

#[then("the topological degeneracy is measured and advantage remains open and unasserted")]
async fn advantage(_w: &mut ProbeWorld) {
    let degeneracy = witness::advantage_probe(&atlas()).unwrap();
    println!(
        "[probe] advantage: topological degeneracy = {degeneracy:.3} braid paths per distinct \
         result κ (isotopy-class collapse via content-addressing; measurement only; advantage \
         OPEN and unasserted)"
    );
    assert!(
        degeneracy >= 1.0,
        "the degeneracy must be a valid measurement (>= 1)"
    );
}

#[then("the Atlas-native MTC construction returns an obstruction")]
async fn t_atlas_native_mtc_obstruction(_w: &mut ProbeWorld) {
    let res = tqc_mtc::native::construct_atlas_native();
    assert!(res.is_err());
}

#[tokio::main]
async fn main() {
    let targets = std::path::Path::new(env!("CARGO_MANIFEST_DIR")).join("../../features/targets");
    ProbeWorld::cucumber()
        .fail_on_skipped()
        .run_and_exit(targets)
        .await;
}
