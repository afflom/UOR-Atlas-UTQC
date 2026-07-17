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

use cucumber::{given, then, when, World};
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
    whitepaper_source: Option<String>,
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
    let f1 = w.f1.clone().unwrap_or_else(|| F1Constants::load().unwrap());
    witness::complex_amplitude_encoding(&w.params(), &f1).unwrap();
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

#[then("the quantum realization is unitary and exhibits destructive interference")]
async fn t_quantum_realization(w: &mut TqcWorld) {
    witness::quantum_realization(&w.params()).unwrap();
}

#[then("the generator closure partitions the class space into the derived mirror orbits covering every label")]
async fn t_generative_closure(w: &mut TqcWorld) {
    witness::generative_closure_probe(&w.params()).unwrap();
}

#[then("every UTQC pillar witness passes in the roll-up conjunction")]
async fn t_utqc_proven(w: &mut TqcWorld) {
    let model = w.model.take().unwrap_or_else(|| Model::load().unwrap());
    let f1 = w.f1.clone().unwrap_or_else(|| F1Constants::load().unwrap());
    witness::utqc_proven_probe(&model, &f1, &w.params()).unwrap();
    w.model = Some(model);
}

#[then("the discrete execution replays identical braid words to identical states and kappa")]
async fn t_fault_tolerance(w: &mut TqcWorld) {
    witness::deterministic_replay_witness(&w.params()).unwrap();
}

#[then("the execution cost is exactly the operation count linear in braid depth with no exponential state")]
async fn t_complexity_bound(w: &mut TqcWorld) {
    witness::complexity_bound_witness(&w.params()).unwrap();
}

#[then("any validator can perfectly mathematically reconstruct the final state and identical kappa from the genesis configuration and braid word")]
async fn t_reconstructability(w: &mut TqcWorld) {
    witness::reconstruction_witness(&w.params()).unwrap();
}

#[then("the Solovay-Kitaev density question is exactly decided")]
async fn t_solovay_kitaev_density(w: &mut TqcWorld) {
    // The decision, asserted as a theorem in both directions: unique 2-dim block, confined
    // to the (-1) eigenspace, tr(P1 G_S) = 0 identically, finite projective Clifford image
    // of exact order 24, density refuted. Decided over Q(zeta_24); no f64 in the verdict.
    witness::solovay_kitaev_decision_witness(&w.params()).unwrap();
}

#[then("the archimedean continuity is exactly located on the 22-dim block and saturates PU(22)")]
async fn t_archimedean_continuity(w: &mut TqcWorld) {
    // The positive half of the density decision, saturated: the projective closure on the
    // 22-dim irreducible block is dense in PU(22) (spectral-flow seed, division-free Lie
    // closure, sound mod-p lower bound 483 forcing su(22)). Universal quantum computation
    // on a 22-dim qudit carrier. Decided over Q(zeta_24); no f64 in the verdict.
    witness::archimedean_continuity_witness(&w.params()).unwrap();
}

#[then("the pair carrier is irreducible and its closure is dense in PU(576)")]
async fn t_pair_carrier(w: &mut TqcWorld) {
    // Pinned theorems: pair irreducibility (commutant 1), the diagonal-sector separation
    // (no monodromy power preserves the block tensor code), the native continuous
    // entangling flow (Lie lower bound exceeds the local subalgebra), and PU(576) density
    // (T1 adj-tensor-adj certificate + T2 reachability rank 92 + classical closure T3),
    // with the n-handle corollary by the two-local composition lemma.
    witness::pair_carrier_witness(&w.params()).unwrap();
}

#[then("the S4 modal logic frame satisfies reflexivity and transitivity")]
async fn t_s4_modal_logic(w: &mut TqcWorld) {
    witness::s4_frame_witness(&w.params()).unwrap();
}

#[then("the Mac Lane pentagon and hexagon identities are verified phase-exactly")]
async fn t_mac_lane_pentagon(w: &mut TqcWorld) {
    witness::mac_lane_coherence(&w.params()).unwrap();
}

#[then(
    "the same topological operator resolves to identical κ across the two independent realizations"
)]
async fn t_universality(w: &mut TqcWorld) {
    witness::equivalency_universality_probe(&w.params()).unwrap();
}

#[then("the Atlas-native MTC construction successfully resolves topological obstructions")]
async fn t_atlas_native_mtc_obstruction(w: &mut TqcWorld) {
    let res = tqc_mtc::native::construct_atlas_native(&w.params());
    assert!(res.is_ok(), "The native MTC construction should now mathematically resolve all prior topological obstructions");
}

#[then("the Grover amplitude recurrence is evaluated exactly at the fixed reference instance")]
async fn t_grover_search(_w: &mut TqcWorld) {
    tqc_algorithms::checks::grover_check().unwrap();
}

#[then("the QFT circuit schedules onto a bounded braid word and executes on the class space")]
async fn t_qft_algorithm(w: &mut TqcWorld) {
    let p = w.params();
    let state = tqc_algorithms::checks::qft_word_check(&p).unwrap();
    // Braid words act by permutation, so the evaluated state must be a permutation of the
    // genesis vector 0..n — a real conservation property of the execution.
    let mut sorted = state;
    sorted.sort_unstable();
    let expected: Vec<i64> = (0..p.class_count() as i64).collect();
    assert_eq!(
        sorted, expected,
        "the QFT braid execution must permute the genesis state, not alter it"
    );
}

#[then(
    "the QPE readout peak is computed by exact integer minimization and meets the exact guarantee"
)]
async fn t_qpe_algorithm(_w: &mut TqcWorld) {
    tqc_algorithms::checks::qpe_check().unwrap();
}

#[then(
    "the Shor period is recovered by exact cyclotomic evaluation and verified against the orbit"
)]
async fn t_shor_algorithm(_w: &mut TqcWorld) {
    tqc_algorithms::checks::shor_check().unwrap();
}

#[tokio::main]
async fn main() {
    let features = std::path::Path::new(env!("CARGO_MANIFEST_DIR")).join("../../features/suites");
    TqcWorld::cucumber()
        .fail_on_skipped()
        .run_and_exit(features)
        .await;
}

#[then("the topological execution manifold bounds non-local entanglement entropy")]
async fn t_topological_entanglement(w: &mut TqcWorld) {
    let result = witness::topological_entanglement_probe(&w.params()).unwrap();
    assert!(
        !result.depth_profile.is_empty() && result.max_schmidt_rank >= 1,
        "the measured Schmidt-rank profile must be non-trivial"
    );
    assert!(
        result.entropy_bound > 0.0,
        "measured entanglement entropy must be positive for the braided states"
    );
    assert!(
        result.is_logarithmic_scaling,
        "the measured profile must saturate within the log2 bound (derived verdict): {:?}",
        result.depth_profile
    );
}

#[given(expr = "the whitepaper source in {string}")]
async fn given_whitepaper_source(w: &mut TqcWorld, path: String) {
    let root = std::path::Path::new(env!("CARGO_MANIFEST_DIR")).join("../../");
    let full_path = root.join(&path);
    let content = std::fs::read_to_string(&full_path).expect("failed to read whitepaper source");
    w.whitepaper_source = Some(content);
}

#[then(expr = "it must use the {string} class")]
async fn then_it_must_use_class(w: &mut TqcWorld, class: String) {
    let content = w.whitepaper_source.as_ref().unwrap();
    assert!(
        content.contains(&class),
        "Whitepaper does not use the required document class: {}",
        class
    );
}

#[then(expr = "it must include tikz diagrams for mathematical visual aids")]
async fn then_it_must_include_tikz(w: &mut TqcWorld) {
    let content = w.whitepaper_source.as_ref().unwrap();
    assert!(
        content.contains("\\usepackage{tikz}"),
        "Whitepaper does not include tikz"
    );
}

#[then("isomorphic topological braid operations naturally collide on identical kappa forms via classical equivalence")]
async fn t_tensor_contraction_bypass(w: &mut TqcWorld) {
    witness::isotopy_collision_witness(&w.params()).unwrap();
}

#[then("the generated braiding subgroup is measured as mathematically finite")]
async fn t_finite_closure(w: &mut TqcWorld) {
    let result = witness::finite_closure_probe(&w.params()).unwrap();
    assert!(
        !result.is_dense,
        "The combinatorial execution manifold requires a finite subgroup to map exponential states to cacheable polynomials"
    );
}

#[when("a two-qubit entangling gate is natively constructed from the abelian category")]
async fn w_two_qubit_entangling_gate_natively_constructed(_w: &mut TqcWorld) {
    // The gate is constructed and validated within the subsequent then-clauses.
}

#[then("the native entangling gate is established without a gate-set density claim")]
async fn t_multi_qubit_universality(w: &mut TqcWorld) {
    let result = witness::two_qubit_universality_probe(&w.params()).unwrap();
    assert!(
        result.is_entangling,
        "A native entangling phase gate must be established (no gate-set density is claimed)"
    );
}

#[then("it does not induce a theory collision with the non-abelian construction")]
async fn t_theory_collision_avoided(w: &mut TqcWorld) {
    let result = witness::two_qubit_universality_probe(&w.params()).unwrap();
    assert!(
        result.is_coherent,
        "The entangling gate must reside strictly in the coherent abelian substrate, avoiding collision"
    );
}
