//! Exact density integration tests.
use tqc_core::UseCaseParams;
use tqc_vv::exact;

/// Exact algebraic density certificate at the atlas use-case (modality 3, context 8).
/// Every decision is made over Q(zeta_24); no f64 value participates in the verdict.
#[test]
fn exact_density_atlas() {
    let p = UseCaseParams::new(4, 3, 8);
    let report = exact::exact_density_certificate(&p).expect("exact certificate must run");
    println!("commutant_dim      = {}", report.commutant_dim);
    println!("block_dim          = {}", report.block_dim);
    println!("beta_s nonzero at  = {:?}", report.beta_s_nonzero);
    println!("beta_t nonzero at  = {:?}", report.beta_t_nonzero);
    println!("noncommuting grade = {:?}", report.noncommuting_grade);
    println!("proj infinite      = {:?}", report.proj_infinite);
    println!("proj pair          = {:?}", report.proj_pair);
    println!("block support      = {:?}", report.block_support);
    println!("finite image order = {:?}", report.finite_image_order);
    println!("block22 infinite   = {:?}", report.block22_infinite);
    println!("block22 pair       = {:?}", report.block22_pair);
    println!("beyond finite      = {}", report.beyond_finite);
    println!("lie dim lower (22) = {}", report.lie_dim_lower_22);
    println!("pu22 dense         = {}", report.pu22_dense);
    println!("code components    = {}", report.code_components);
    println!("native entangler   = {:?}", report.native_code_entangler);
    println!("pair commutant dim = {}", report.pair_commutant_dim);
    println!("qudit universal    = {}", report.qudit_universal);
    println!("pair lie dim lower = {}", report.pair_lie_dim_lower);
    println!("pair entangl flow  = {}", report.pair_entangling_flow);
    println!("pair adj component = {}", report.pair_adj_component);
    println!("pair reach rank    = {}", report.pair_reach_rank);
    println!("pu576 dense        = {}", report.pu576_dense);
    println!("gate level univ    = {}", report.gate_level_universal);
    println!("certified_dense    = {}", report.certified_dense);
    println!("{}", report.description);
    assert_eq!(report.commutant_dim, 2, "exact commutant dimension");
    assert_eq!(report.block_dim, 2, "exact block dimension");
    // Kernel-grade findings at the atlas use-case: the unique 2-dim invariant block lies
    // inside the (-1) eigenspace, the coupling is a global phase there, and the projective
    // image is finite. Density on the block is refuted, not certified.
    assert!(
        report.beta_s_nonzero.is_empty(),
        "tr(P1 G_S) = 0 identically"
    );
    assert_eq!(report.beta_t_nonzero, vec![-1], "u_t trace grade");
    assert_eq!(
        report.block_support,
        vec![(10, 0.0), (7, 0.0), (2, 0.0), (-1, 2.0)],
        "block supported entirely in the (-1) eigenspace"
    );
    assert!(!report.certified_dense, "density on the block is refuted");
    assert!(
        report.finite_image_order.is_some(),
        "projective image is finite"
    );
    // Direct, threshold-free native-entangling-flow certificate: an explicit Lie(H_2) element
    // outside the local subalgebra u(24)(x)1 + 1(x)u(24).
    println!("pair nonlocal wit  = {:?}", report.pair_nonlocal_witness);
    assert!(
        report.pair_entangling_flow && report.pair_nonlocal_witness.is_some(),
        "native entangling flow must be certified directly (explicit non-local element)"
    );
}
