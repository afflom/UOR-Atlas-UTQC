//! The BDD suite. Each Gherkin step binds to a `tqc-vv` witness (DRY); the runner fails on any
//! failed, skipped, or undefined step, so "fully implemented — no narrowing" is mechanical.
//!
//! Run with `just bdd` (or `cargo test -p tqc-conformance --test bdd`).

// A custom-harness cucumber binary signals failure by panicking in steps and exiting non-zero.
#![allow(
    clippy::unwrap_used,
    clippy::expect_used,
    clippy::panic,
    clippy::missing_panics_doc
)]

use cucumber::{given, then, World};
use tqc_atlas::canonical;
use tqc_core::generators::Generators;
use tqc_core::{labels, UseCaseParams};
use tqc_model::Model;
use tqc_vv::{witness, F1Constants};

#[derive(Debug, Default, cucumber::World)]
struct TqcWorld {
    model: Option<Model>,
    f1: Option<F1Constants>,
    params: Option<UseCaseParams>,
}

impl TqcWorld {
    fn f1(&self) -> &F1Constants {
        self.f1
            .as_ref()
            .expect("the F1 oracle constants step must run first")
    }
    fn params(&self) -> UseCaseParams {
        self.params.expect("a use-case step must run first")
    }
}

#[given("the F1 oracle constants")]
async fn f1_oracle(w: &mut TqcWorld) {
    let model = Model::load().unwrap();
    let f1 = F1Constants::load().unwrap();
    // Loading the oracle also asserts its provenance (sha256 == manifest pin).
    witness::oracle_provenance(&model, &f1).unwrap();
    w.model = Some(model);
    w.f1 = Some(f1);
}

#[given("the UOR Atlas use-case")]
async fn atlas_use_case(w: &mut TqcWorld) {
    let model = w.model.take().unwrap_or_else(|| Model::load().unwrap());
    w.params = Some(canonical(&model).unwrap());
    w.model = Some(model);
}

#[given(expr = "an arbitrary use-case with scope {int} modality {int} context {int}")]
async fn arbitrary_use_case(w: &mut TqcWorld, q: u32, t: u32, o: u32) {
    w.params = Some(UseCaseParams::checked(q, t, o).unwrap());
}

#[then("the objects-labels witness reproduces the F1 Atlas")]
async fn t_objects(w: &mut TqcWorld) {
    witness::objects_labels(&w.params(), w.f1()).unwrap();
}

#[then("classIndex is a bijection over the whole class space")]
async fn t_bijection(w: &mut TqcWorld) {
    assert!(labels::class_index_is_bijection(&w.params()));
}

#[then("the label-space belt witness reproduces the F1 Atlas")]
async fn t_belt(w: &mut TqcWorld) {
    witness::label_space_belt(&w.params(), w.f1()).unwrap();
}

#[then("the inner product is the definite Euclidean companion")]
async fn t_inner(w: &mut TqcWorld) {
    witness::inner_product(&w.params()).unwrap();
}

#[then("the generators have the F1 orders and preserve the inner product")]
async fn t_generators(w: &mut TqcWorld) {
    witness::reflection_generators(&w.params(), w.f1()).unwrap();
}

#[then("the generators have orders scope, context and two")]
async fn t_generator_orders(w: &mut TqcWorld) {
    let p = w.params();
    let g = Generators::new(&p);
    assert_eq!(g.sigma.order(), u64::from(p.sigma_order()));
    assert_eq!(g.tau.order(), u64::from(p.tau_order()));
    assert_eq!(g.mu.order(), u64::from(p.mu_order()));
}

#[then(
    "the absolute structural quotient of the composition algebra forms an associative fusion ring"
)]
async fn structural_quotient_associative(w: &mut TqcWorld) {
    let p = w.params();
    assert!(tqc_core::octonion::absolute_quotient_is_associative(
        p.context as usize
    ));
}

#[then("the spectrum reconciles with the F1 multiplicities and signature")]
async fn t_spectrum(w: &mut TqcWorld) {
    witness::spectrum(&w.params(), w.f1()).unwrap();
}

#[then("the Coxeter rank equals phi of the Coxeter number and the context")]
async fn t_coxeter(w: &mut TqcWorld) {
    witness::coxeter_weyl(&w.params(), w.f1()).unwrap();
}

#[then("the modular identity holds on the F1 coefficients")]
async fn t_modular(w: &mut TqcWorld) {
    witness::modular_identities(&w.params(), w.f1()).unwrap();
}

#[then("the E8 definite anchor reproduces the F1 Atlas")]
async fn t_e8_anchor(w: &mut TqcWorld) {
    witness::definite_anchor_e8(w.f1()).unwrap();
}

#[then("the definite anchor is positive-definite")]
async fn t_anchor_pd(w: &mut TqcWorld) {
    witness::definite_anchor(&w.params()).unwrap();
}

#[then("the composition reduces to the realized g2 product on every sigma-axis")]
async fn t_fusion(w: &mut TqcWorld) {
    witness::fusion_g2(&w.params()).unwrap();
}

#[then("the dual reduces to the realized f4 mirror involution")]
async fn t_dual(w: &mut TqcWorld) {
    witness::dual_f4(&w.params()).unwrap();
}

#[then("the e6/e7/e8 operations reduce to the realized operations")]
async fn t_categorical(w: &mut TqcWorld) {
    witness::categorical_structure(&w.params()).unwrap();
}

#[then("the ground space round-trips with no loss")]
async fn t_ground_space(w: &mut TqcWorld) {
    witness::ground_space_protection(&w.params()).unwrap();
}

#[then("the complex amplitude encoding round-trips and preserves the norm")]
async fn t_amplitude(w: &mut TqcWorld) {
    witness::complex_amplitude_encoding(&w.params()).unwrap();
}

#[then("the modular S and T satisfy the SL(2,Z) relations")]
async fn t_modular_st(w: &mut TqcWorld) {
    witness::modular_s_t(&w.params()).unwrap();
}

#[then("the braiding R satisfies the hexagon and Yang-Baxter")]
async fn t_braiding(w: &mut TqcWorld) {
    witness::braiding_r_matrix(&w.params()).unwrap();
}

#[then("the braid-fuse-read holospace cycle runs and round-trips")]
async fn t_holospace_cycle(w: &mut TqcWorld) {
    witness::holospace_cycle(&w.params()).unwrap();
}

#[tokio::main]
async fn main() {
    let features = std::path::Path::new(env!("CARGO_MANIFEST_DIR")).join("../../features/suites");
    TqcWorld::cucumber()
        .fail_on_skipped()
        .run_and_exit(features)
        .await;
}
