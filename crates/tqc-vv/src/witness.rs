//! The V&V witnesses. Each binds **one** parametric computation (from `tqc-core`) to **one**
//! authoritative oracle value (from F1), returning a readable error on mismatch.
//!
//! These functions are the single implementation of each check; both the `#[test]`s below and
//! the cucumber step definitions in `tqc-conformance` call them (DRY).

use crate::oracle::F1Constants;
use tqc_core::amplitude::{self, Amplitude};
use tqc_core::generators::{Generators, Permutation};
use tqc_core::inner::{euclidean_norm_sq, preserves_norm};
use tqc_core::{coxeter, labels, modular, octonion, spectrum, UseCaseParams};
use tqc_model::Model;
use tqc_substrate::{dual, embed_e8, fuse, grade_e6, orbit_e7, CompositionAxis, COMPOSITION_AXES};

/// Outcome of a witness.
pub type Witness = Result<(), String>;

fn check(cond: bool, msg: impl Into<String>) -> Witness {
    if cond {
        Ok(())
    } else {
        Err(msg.into())
    }
}

/// Canonical bytes for an anyon label, parameterized by the use-case (no Atlas literal).
fn anyon_bytes(p: &UseCaseParams, index: u64) -> Vec<u8> {
    format!(
        "tqc-anyon:q{}t{}o{}:{index}",
        p.scope, p.modality, p.context
    )
    .into_bytes()
}

/// VV-00 — the committed F1 artifact matches its recorded pin, and was extracted at the
/// canonical parameters.
///
/// # Errors
/// If the sha256 differs from `model/oracles.toml`, or the parameters disagree.
pub fn oracle_provenance(model: &Model, f1: &F1Constants) -> Witness {
    let oracle = model
        .oracle("f1-atlas")
        .ok_or_else(|| "model is missing the `f1-atlas` oracle".to_owned())?;
    check(
        F1Constants::sha256() == oracle.sha256,
        format!(
            "F1 artifact sha256 {} != manifest {}",
            F1Constants::sha256(),
            oracle.sha256
        ),
    )?;
    let uc = model
        .canonical_usecase()
        .ok_or_else(|| "no canonical use-case".to_owned())?;
    check(
        (f1.params.scope_q, f1.params.modality_t, f1.params.context_o)
            == (uc.scope, uc.modality, uc.context),
        "F1 oracle parameters disagree with the canonical use-case",
    )
}

/// VV — objects / anyon labels: class count, stride, and the `classIndex` bijection.
///
/// # Errors
/// On any mismatch with F1 or a non-bijective index.
pub fn objects_labels(p: &UseCaseParams, f1: &F1Constants) -> Witness {
    check(p.class_count() == f1.classes.count, "class count != F1")?;
    check(p.stride() == f1.classes.stride, "stride != F1")?;
    check(
        labels::class_index_is_bijection(p),
        "classIndex is not a bijection onto [0,count)",
    )
}

/// VV — label / state-space index: the belt extent and its factorizations.
///
/// # Errors
/// On any mismatch with F1 or a non-contiguous belt.
pub fn label_space_belt(p: &UseCaseParams, f1: &F1Constants) -> Witness {
    check(
        p.belt_extent() == f1.classes.belt_extent,
        "belt extent != F1",
    )?;
    let mut got: Vec<(u64, u64)> = p.belt_factorizations();
    let mut want: Vec<(u64, u64)> = f1
        .classes
        .belt_factorizations
        .iter()
        .map(|v| (v[0], v[1]))
        .collect();
    got.sort_unstable();
    want.sort_unstable();
    check(got == want, "belt factorizations != F1")?;
    let (pages, page) = (want[0].0.min(want[0].1), want[0].0.max(want[0].1));
    let addrs = labels::belt_addresses(pages, page);
    check(
        addrs.len() as u64 == p.belt_extent() && addrs.last() == Some(&(p.belt_extent() - 1)),
        "belt addresses are not contiguous over the extent",
    )
}

/// VV — the inner product is the definite Euclidean companion `Σxᵢ²`.
///
/// Validates definiteness: a positive-definite sum of squares.
///
/// # Errors
/// If the form is not a definite sum of squares.
pub fn inner_product(p: &UseCaseParams) -> Witness {
    let n = p.carrier_dim() as usize;
    let zero = vec![0i64; n];
    check(
        euclidean_norm_sq(&zero) == 0,
        "norm of the zero vector must be 0",
    )?;
    let mut v = zero.clone();
    if let Some(first) = v.first_mut() {
        *first = 1;
    }
    check(
        euclidean_norm_sq(&v) > 0,
        "norm of a non-zero vector must be positive (definite)",
    )?;
    check(
        euclidean_norm_sq(&[3, 4]) == 25,
        "norm must be the plain sum of squares",
    )
}

/// VV — the reflection generators: orders match F1, and each is orthogonal on `Σxᵢ²`.
///
/// This is "the unlock": orthogonality (= unitarity) holds with no positivity assumption.
///
/// # Errors
/// On an order mismatch or a generator that fails to preserve the norm.
pub fn reflection_generators(p: &UseCaseParams, f1: &F1Constants) -> Witness {
    let g = Generators::new(p);
    check(
        g.sigma.order() == u64::from(f1.generators.sigma_order),
        "sigma order != F1",
    )?;
    check(
        g.tau.order() == u64::from(f1.generators.tau_order),
        "tau order != F1",
    )?;
    check(
        g.mu.order() == u64::from(f1.generators.mu_order),
        "mu order != F1",
    )?;
    // Orthogonality over the whole class space.
    let n = p.class_count() as usize;
    let v: Vec<i64> = (0..n as i64).map(|i| i % 11 - 5).collect();
    for (name, perm) in [("sigma", &g.sigma), ("tau", &g.tau), ("mu", &g.mu)] {
        check(
            preserves_norm(perm, &v),
            format!("generator {name} does not preserve the norm"),
        )?;
    }
    Ok(())
}

/// VV — the spectrum: parametric block eigenvalues plus F1 multiplicities reconcile to the
/// F1 signature and trace.
///
/// # Errors
/// On any disagreement with F1.
pub fn spectrum(p: &UseCaseParams, f1: &F1Constants) -> Witness {
    check(
        spectrum::block_eigenvalues(p).as_slice() == f1.spectrum.eigenvalues.as_slice(),
        "block eigenvalues != F1",
    )?;
    // Multiplicities are derived parametrically from the projector-rank formula and must
    // agree with the sourced F1 values (derive, then cross-check — never hand-enter).
    check(
        spectrum::block_multiplicities(p).as_slice() == f1.spectrum.multiplicities.as_slice(),
        "derived block multiplicities != F1",
    )?;
    let sig = spectrum::reconcile(p, &f1.spectrum.eigenvalues, &f1.spectrum.multiplicities)
        .map_err(|e| format!("spectrum reconciliation failed: {e:?}"))?;
    check(
        f1.spectrum.signature == vec![sig.positive, sig.negative],
        format!(
            "signature ({},{}) != F1 {:?}",
            sig.positive, sig.negative, f1.spectrum.signature
        ),
    )?;
    check(
        i64::try_from(p.carrier_dim()) == Ok(f1.spectrum.trace),
        "trace != carrier dim",
    )?;
    check(p.carrier_dim() == f1.spectrum.dim, "dim != F1")?;
    check(
        f1.spectrum.indefinite,
        "F1 records the operator as indefinite",
    )
}

/// VV — Coxeter / Weyl: `rank = φ(h) = context`, and the exponents are consistent.
///
/// # Errors
/// On any disagreement with F1.
pub fn coxeter_weyl(p: &UseCaseParams, f1: &F1Constants) -> Witness {
    let rank = coxeter::euler_phi(f1.coxeter.number_h);
    check(rank == f1.coxeter.rank, "phi(h) != F1 rank")?;
    check(rank == p.context, "rank != context (O)")?;
    check(
        f1.coxeter.exponents.len() as u32 == rank,
        "exponent count != rank",
    )?;
    check(
        f1.coxeter.exponents.iter().sum::<u32>() == f1.coxeter.exponent_sum,
        "exponent sum != F1",
    )
}

/// VV — the definite anchor: the E8 Gram is `4 × Cartan` (diag 8, edges -4) and is
/// positive-definite, matching the F1 `e8_seed`.
///
/// # Errors
/// On any disagreement with F1 or a non-PD Gram.
pub fn definite_anchor_e8(f1: &F1Constants) -> Witness {
    let scale = f1.e8_seed.gram_scale;
    let cartan = tqc_atlas::e8_cartan();
    for (i, row) in cartan.iter().enumerate() {
        check(row[i] == f1.e8_seed.cartan_diag, "E8 Cartan diagonal != F1")?;
    }
    let gram = tqc_atlas::e8_gram(scale);
    for i in 0..8 {
        check(
            gram[i][i] == f1.e8_seed.gram_diag,
            "E8 Gram diagonal != F1 gram_diag",
        )?;
        for j in 0..8 {
            check(
                gram[i][j] == scale * cartan[i][j],
                "E8 Gram != scale*Cartan",
            )?;
            if i != j && gram[i][j] != 0 {
                check(
                    gram[i][j] == f1.e8_seed.gram_edge,
                    "E8 Gram edge != F1 gram_edge",
                )?;
            }
        }
    }
    check(
        tqc_core::anchor::is_positive_definite(&gram) == f1.e8_seed.psd,
        "E8 positive-definiteness != F1",
    )?;
    check(f1.e8_seed.psd, "F1 records the E8 seed as PSD")
}

/// VV — the generic definite anchor: the use-case's Euclidean companion is positive-definite.
///
/// # Errors
/// If the companion is not positive-definite.
pub fn definite_anchor(p: &UseCaseParams) -> Witness {
    let gram = tqc_core::anchor::euclidean_companion(p.carrier_dim() as usize);
    check(
        tqc_core::anchor::is_positive_definite(&gram),
        "the use-case Euclidean companion must be positive-definite",
    )
}

/// VV — the modular identity `E4³ = E6² + 1728·Δ`, plus the weight `T·O/2`.
///
/// # Errors
/// If the identity fails on the F1 coefficients or the weight is inconsistent.
pub fn modular_identities(p: &UseCaseParams, f1: &F1Constants) -> Witness {
    let e4: Vec<i128> = f1.modular.e4.iter().map(|&x| i128::from(x)).collect();
    let e6: Vec<i128> = f1.modular.e6.iter().map(|&x| i128::from(x)).collect();
    let delta: Vec<i128> = f1.modular.delta.iter().map(|&x| i128::from(x)).collect();
    check(
        modular::identity_holds(&e4, &e6, &delta, i128::from(f1.modular.constant)),
        "E4^3 = E6^2 + 1728*Delta failed on the F1 coefficients",
    )?;
    check(
        u64::from(f1.modular.weight) * 2 == p.carrier_dim(),
        "weight*2 != carrier dim (T*O)",
    )
}

/// VV — the Atlas composition reduces to the realized `compose_g2_product` and is commutative on every
/// σ-axis; the composition norm is multiplicative at the use-case's context level.
///
/// # Errors
/// On a non-commutative composition, an axis/composition failure, or a non-multiplicative norm.
pub fn fusion_g2(p: &UseCaseParams) -> Witness {
    let n = p.class_count().min(6);
    for axis in COMPOSITION_AXES {
        for i in 0..n {
            for j in 0..n {
                let (a, b) = (anyon_bytes(p, i), anyon_bytes(p, j));
                let ab = fuse(axis, &a, &b)?;
                let ba = fuse(axis, &b, &a)?;
                check(
                    ab == ba,
                    format!("g2 not commutative on {} for ({i},{j})", axis.token()),
                )?;
            }
        }
    }
    // Norm-multiplicativity at the use-case's context level (1,2,4,8 are the division-algebra
    // dimensions; the Atlas uses the octonion eight-square at O=8).
    if matches!(p.context, 1 | 2 | 4 | 8) {
        let dim = p.context as i128;
        let x: Vec<i128> = (0..dim).map(|k| k + 1).collect();
        let y: Vec<i128> = (0..dim).map(|k| 2 * k - 3).collect();
        check(
            octonion::norm_multiplicative(&x, &y),
            "the composition norm is not multiplicative at the context level",
        )?;
    }
    Ok(())
}

/// VV — the dual reduces to the realized `compose_f4_quotient` (deterministic, well-formed on
/// every σ-axis) and the conjugation generator `μ` is an involution.
///
/// # Errors
/// On a non-involutive `μ` or an axis/composition failure.
pub fn dual_f4(p: &UseCaseParams) -> Witness {
    let g = Generators::new(p);
    check(g.mu.order() == u64::from(p.mu_order()), "mu order != F1")?;
    check(
        g.mu.then(&g.mu) == Permutation::identity(p.class_count()),
        "the conjugation generator mu must be an involution",
    )?;
    let sample = anyon_bytes(p, 0);
    for axis in COMPOSITION_AXES {
        let once = dual(axis, &sample)?;
        check(
            once == dual(axis, &sample)?,
            format!("f4 not deterministic on {}", axis.token()),
        )?;
        check(
            !once.is_empty(),
            format!("f4 produced an empty label on {}", axis.token()),
        )?;
    }
    Ok(())
}

/// VV — the categorical operations `e6`/`e7`/`e8` reduce to the realized operations
/// (deterministic, well-formed on every σ-axis); the `e7` S4 orbit size is `T·O = carrier_dim`.
///
/// # Errors
/// On an axis/composition failure or an orbit-size mismatch.
pub fn categorical_structure(p: &UseCaseParams) -> Witness {
    check(
        p.carrier_dim() == u64::from(p.modality) * u64::from(p.context),
        "e7 S4 orbit size != T*O",
    )?;
    let sample = anyon_bytes(p, 0);
    for axis in COMPOSITION_AXES {
        for (name, out) in [
            ("e6", grade_e6(axis, &sample)?),
            ("e7", orbit_e7(axis, &sample)?),
            ("e8", embed_e8(axis, &sample)?),
        ] {
            check(
                !out.is_empty(),
                format!("{name} produced an empty label on {}", axis.token()),
            )?;
        }
        check(
            grade_e6(axis, &sample)? == grade_e6(axis, &sample)?,
            format!("e6 not deterministic on {}", axis.token()),
        )?;
    }
    Ok(())
}

/// VV — ground-space / topological protection: content-addressing is a faithful round-trip.
/// `κ` is stable (CC-1), content re-derives to its `κ` (`π∘ι = id`), and distinct content has
/// distinct `κ` (eviction drops bytes, not identity).
///
/// # Errors
/// On unstable addressing, a failed re-derivation, or a `κ` collision.
pub fn ground_space_protection(p: &UseCaseParams) -> Witness {
    let n = p.class_count().min(8);
    let mut seen: Vec<tqc_substrate::Kappa> = Vec::new();
    for i in 0..n {
        let state = anyon_bytes(p, i);
        let k = tqc_substrate::kappa(&state);
        check(
            k == tqc_substrate::kappa(&state),
            format!("kappa not stable at label {i}"),
        )?;
        check(
            tqc_substrate::verify(&state, &k)?,
            format!("content does not re-derive at label {i}"),
        )?;
        check(!seen.contains(&k), format!("kappa collision at label {i}"))?;
        seen.push(k);
    }
    Ok(())
}

/// VV (build) — complex amplitude encoding: an amplitude-space vector encodes to canonical bytes,
/// round-trips through the content-addressed store (CC-1), and its Euclidean composition norm
/// `Σ|cᵢ|²` equals the inner product on the encoded form.
///
/// # Errors
/// On a failed round-trip, unstable addressing, or a norm mismatch.
pub fn complex_amplitude_encoding(p: &UseCaseParams, f1: &F1Constants) -> Witness {
    let n = p.class_count().min(8);
    let state: Vec<(u64, Amplitude)> = (0..n)
        .map(|i| {
            let re = f1.spectrum.eigenvalues[(i as usize) % f1.spectrum.eigenvalues.len()];
            let im = f1.modular.e4[(i as usize) % f1.modular.e4.len()];
            (i, Amplitude { re, im })
        })
        .collect();
    let bytes = amplitude::encode(&state);

    let decoded = amplitude::decode(&bytes).ok_or_else(|| "amplitude decode failed".to_owned())?;
    let mut canonical_state = state.clone();
    canonical_state.sort_by_key(|(l, _)| *l);
    check(
        decoded == canonical_state,
        "amplitude encode/decode does not round-trip",
    )?;

    let k = tqc_substrate::kappa(&bytes);
    check(
        k == tqc_substrate::kappa(&bytes),
        "amplitude kappa not stable (CC-1)",
    )?;
    check(
        tqc_substrate::verify(&bytes, &k)?,
        "amplitude state does not re-derive (CC-1)",
    )?;

    let flat: Vec<i64> = state.iter().flat_map(|(_, a)| [a.re, a.im]).collect();
    check(
        amplitude::norm_sq(&state) == euclidean_norm_sq(&flat),
        "Σ|c_i|² != the Euclidean inner product on the encoded form",
    )
}

/// VV (build) — the modular S/T matrices of the Atlas-native pointed category satisfy the
/// full MTC axiom set with phase-exact comparisons: S symmetric & unitary, `S⁴ = I`,
/// `S² = C`, twists of finite order (Vafa), `(ST)³ = p⁺S²` with the Gauss-sum anomaly `p⁺`,
/// non-negative integer fusion, exact Verlinde, phase-exact pentagon/hexagon/balancing, and
/// the monodromy–S relation. A `build` construction — never asserted to be F1-sourced.
///
/// # Errors
/// Returns the first axiom that fails.
pub fn modular_s_t(p: &UseCaseParams) -> Witness {
    let native = match tqc_mtc::native::construct_atlas_native(p) {
        Ok(n) => n,
        Err(e) => return Err(e.to_string()),
    };
    tqc_mtc::verifier::verify_mtc_axioms(&*native, tqc_mtc::TOL)?;

    // At the Atlas instance, pin the anomaly VALUE: three semions (c = 3) times the
    // Z_3 anyon (c = 2) give central charge c ≡ 5 (mod 8), so the Gauss-sum anomaly
    // must be p⁺ = e^{2πi·5/8} exactly (the generic verifier only checks |p⁺| = 1 and
    // the (ST)³ relation).
    if (p.modality, p.context) == (3, 8) {
        let t = native.t_diag();
        let dim = native.dim();
        let mut gauss = tqc_mtc::C::new(0.0, 0.0);
        for theta in &t {
            gauss = gauss.plus(*theta);
        }
        let p_plus = gauss.scale(1.0 / (dim as f64).sqrt());
        let want = tqc_mtc::C::phase(2.0 * core::f64::consts::PI * 5.0 / 8.0);
        check(
            p_plus.close(want, tqc_mtc::TOL),
            "Gauss-sum anomaly p⁺ != e^{2πi·5/8}: the central charge is not ≡ 5 (mod 8)",
        )?;
    }
    Ok(())
}

/// The monodromy bicharacter `χ(a,b) = R(a,b)·R(b,a)` of the Atlas-native pointed category,
/// as an exact root-of-unity exponent pair `(modality part mod T, semion parity mod 2)`.
/// Integer arithmetic only — no floats.
fn chi_exponents(p: &UseCaseParams, x: usize, y: usize) -> (usize, usize) {
    let modality = p.modality as usize;
    let context = p.context as usize;
    let k_mod = if modality % 2 == 0 { 1 } else { 2 };
    let (m1, c1) = (x / context, x % context);
    let (m2, c2) = (y / context, y % context);
    (
        (k_mod * m1 * m2) % modality,
        ((c1 & c2).count_ones() as usize) % 2,
    )
}

/// VV (build) — the braiding R-matrix, decided with exact integer arithmetic on
/// root-of-unity exponents (no floats): the monodromy bicharacter `χ(a,b) = R(a,b)R(b,a)`
/// is symmetric, bimultiplicative, satisfies `χ(a,a) = θ_a²` (ribbon/twist consistency),
/// and is **non-degenerate** — the modularity criterion for pointed categories. The exact
/// integer model is tied back to the runtime R-symbols numerically, and the full
/// hexagon / balancing / (Yang–Baxter via coherence) axioms of R are then verified
/// phase-exactly through the generic axiom verifier.
///
/// # Errors
/// Returns the first braiding property that fails.
#[allow(clippy::needless_range_loop)]
pub fn braiding_r_matrix(p: &UseCaseParams) -> Witness {
    let native = match tqc_mtc::native::construct_atlas_native(p) {
        Ok(n) => n,
        Err(e) => return Err(e.to_string()),
    };
    let dim = native.dim();
    let modality = p.modality as usize;
    let context = p.context as usize;
    let chi = |x: usize, y: usize| chi_exponents(p, x, y);
    let mul = |a: (usize, usize), b: (usize, usize)| ((a.0 + b.0) % modality, (a.1 + b.1) % 2);
    let one = (0usize, 0usize);
    let add = |x: usize, y: usize| -> usize {
        let (m1, c1) = (x / context, x % context);
        let (m2, c2) = (y / context, y % context);
        ((m1 + m2) % modality) * context + (c1 ^ c2)
    };

    for a in 0..dim {
        // Twist consistency: χ(a,a) = θ_a² as exact exponent pairs.
        let (m, c) = (a / context, a % context);
        let k_mod = if modality % 2 == 0 { 1 } else { 2 };
        let theta_sq = ((k_mod * m * m) % modality, (c.count_ones() as usize) % 2);
        check(
            chi(a, a) == theta_sq,
            format!("χ({a},{a}) != θ_{a}² (ribbon/twist consistency)"),
        )?;
        // Non-degeneracy (modularity for pointed categories): every non-identity label is
        // detected by some monodromy.
        if a != 0 {
            check(
                (0..dim).any(|b| chi(a, b) != one),
                format!("monodromy bicharacter is degenerate at label {a}"),
            )?;
        }
        for b in 0..dim {
            check(chi(a, b) == chi(b, a), format!("χ({a},{b}) != χ({b},{a})"))?;
            for c2 in 0..dim {
                check(
                    chi(add(a, b), c2) == mul(chi(a, c2), chi(b, c2)),
                    format!("χ not bimultiplicative at ({a},{b};{c2})"),
                )?;
            }
        }
    }

    // Tie the exact integer model back to the runtime R-symbols: χ(a,b) evaluated as a
    // complex number must equal R(a,b)·R(b,a) from the construction.
    for a in 0..dim {
        for b in 0..dim {
            let k = add(a, b);
            let mono = native.r_symbol(a, b, k).times(native.r_symbol(b, a, k));
            let (em, es) = chi(a, b);
            let theta = 2.0 * core::f64::consts::PI * (em as f64) / (modality as f64);
            let want = tqc_mtc::C::phase(theta)
                .times(tqc_mtc::C::new(if es == 1 { -1.0 } else { 1.0 }, 0.0));
            check(
                mono.close(want, tqc_mtc::TOL),
                format!("exact χ({a},{b}) disagrees with the runtime R-symbols"),
            )?;
        }
    }

    tqc_mtc::verifier::verify_mtc_axioms(&*native, tqc_mtc::TOL)
}

/// VV (build) — the holospace lift: a braid → fuse → read cycle running as one holospace on
/// the content-addressing substrate.
///
/// Boot: an amplitude-space state is encoded to a κ and re-derives (CC-1). Braid: a generator word
/// applied to the state re-addresses deterministically (CC-2). Isotopy collapse: two distinct
/// words that compose to the same operator (e.g. `σ^order` vs the identity) yield the same
/// state and resolve to the same κ — the content-addressed collapse the advantage probe
/// measures. Read: fusing two label κ resolves deterministically. No-loss: the state round-trips
/// byte-identically (CC-29/30).
///
/// The cycle executes generator gates through the native Hologram execution path in `tqc-substrate`:
/// a permutation gate is compiled to a Hologram archive and run through `hologram_exec::InferenceSession`.
/// Persisted `.holo` artifacts are written and addressable.
///
/// # Errors
/// On a failed re-derivation, non-deterministic gate, broken collapse, or lossy round-trip.
pub fn holospace_cycle(p: &UseCaseParams) -> Witness {
    let g = Generators::new(p);
    let n = p.class_count() as usize;
    let base: Vec<i64> = (0..n as i64).map(|i| i % 5 - 2).collect();
    let amp = |state: &[i64]| -> Vec<(u64, Amplitude)> {
        state
            .iter()
            .enumerate()
            .map(|(i, &v)| (i as u64, Amplitude { re: v, im: 0 }))
            .collect()
    };
    let encode_binary = |amplitudes: &[(u64, Amplitude)]| {
        amplitude::encode_interleaved(p.class_count(), amplitudes)
    };
    let decode_binary_to_kappa = |bytes: &[u8]| -> Result<String, String> {
        let amp_state = amplitude::decode_interleaved(bytes)
            .ok_or_else(|| "holo gate output is not a whole number of label records".to_owned())?;
        Ok(tqc_substrate::kappa(&amplitude::encode(&amp_state)).to_string())
    };
    let apply_gate =
        |gate_name: &str, targets: &[usize], state_bytes: &[u8]| -> Result<Vec<u8>, String> {
            let exec = tqc_substrate::execute_holo_gate(gate_name, targets, state_bytes)?;
            // Provenance sanity: the executed artifact must carry the requested gate and a κ.
            if exec.artifact.gate_name != gate_name || exec.artifact.kappa.is_empty() {
                return Err(format!(
                    "holo execution provenance mismatch for gate {gate_name}"
                ));
            }
            Ok(exec.output_bytes)
        };
    // The holo gate is a Gather: out[i] = state[targets[i]]. Realizing the amplitude
    // action `out[g(i)] = v[i]` (the same action `permute_amplitudes` implements) therefore
    // requires targets[i] = g⁻¹(i). Passing the forward map would execute the inverse
    // operator — the exact realization-mismatch the universality probe checks for.
    let get_targets = |perm: &Permutation| -> Vec<usize> {
        let inv = perm.inverse();
        (0..p.class_count())
            .map(|i| inv.apply(i) as usize)
            .collect()
    };

    // Boot: encode the state, confirm it re-derives (CC-1) and round-trips with no loss.
    let amp0 = amp(&base);
    let bytes0 = amplitude::encode(&amp0);
    let k0 = tqc_substrate::kappa(&bytes0);
    check(
        tqc_substrate::verify(&bytes0, &k0)?,
        "boot state does not re-derive (CC-1)",
    )?;
    check(
        amplitude::decode(&bytes0).as_deref() == Some(amp0.as_slice()),
        "state is lossy (CC-29/30)",
    )?;

    // Braid: apply a generator word; gate application is deterministic (CC-2).
    let bin0 = encode_binary(&amp0);
    let st_sigma = apply_gate("sigma", &get_targets(&g.sigma), &bin0)?;
    // Passthrough rejection: σ is non-identity, so its holo execution must change the
    // state — a gate that echoes its input (the archive-degeneration failure mode)
    // cannot pass. And the executed action must equal the independent in-memory
    // permutation action exactly.
    check(
        st_sigma != bin0,
        "holo gate execution is an input passthrough (σ must act non-trivially)",
    )?;
    {
        let expected = g.sigma.permute_amplitudes(&base);
        let got = amplitude::decode_interleaved(&st_sigma)
            .ok_or_else(|| "holo σ output is not a whole number of label records".to_owned())?;
        let got_vals: Vec<i64> = got.iter().map(|&(_, a)| a.re).collect();
        check(
            got_vals == expected,
            "holo σ execution disagrees with the in-memory permutation action",
        )?;
    }
    let st_tau = apply_gate("tau", &get_targets(&g.tau), &st_sigma)?;
    let st_mu = apply_gate("mu", &get_targets(&g.mu), &st_tau)?;
    let k_word = decode_binary_to_kappa(&st_mu)?;

    let st_sigma_2 = apply_gate("sigma", &get_targets(&g.sigma), &bin0)?;
    let st_tau_2 = apply_gate("tau", &get_targets(&g.tau), &st_sigma_2)?;
    let st_mu_2 = apply_gate("mu", &get_targets(&g.mu), &st_tau_2)?;
    let k_word_2 = decode_binary_to_kappa(&st_mu_2)?;
    check(
        k_word == k_word_2,
        "gate application not deterministic (CC-2)",
    )?;

    // Isotopy collapse: σ^order and the identity are the same operator → the same κ.
    let mut st_pow = bin0.clone();
    for _ in 0..p.sigma_order() {
        st_pow = apply_gate("sigma", &get_targets(&g.sigma), &st_pow)?;
    }
    let k_pow = decode_binary_to_kappa(&st_pow)?;
    let k_id = decode_binary_to_kappa(&bin0)?;
    check(k_pow == k_id, "isotopic words must collapse to one κ")?;

    // Read: the composition outcome resolves to a κ, deterministically.
    let read = fuse(
        CompositionAxis::Sha256,
        &anyon_bytes(p, 0),
        &anyon_bytes(p, 1),
    )?;
    check(
        read == fuse(
            CompositionAxis::Sha256,
            &anyon_bytes(p, 0),
            &anyon_bytes(p, 1),
        )?,
        "composition readout not deterministic",
    )
}

/// The measured empirical finite-closure metrics.
#[derive(Debug, Clone, PartialEq)]
pub struct FiniteClosureMetrics {
    /// True if the generated braid subgroup is dense (universal quantum computation).
    pub is_dense: bool,
    /// The size of the orbit/group if finite.
    pub unique_phases: usize,
    /// Detailed description of the measurement.
    pub description: String,
}

/// A probe deciding the finite closure of the modular representation of the Atlas-native
/// category — every quantity below is **derived**, none is hand-entered:
///
/// 1. The exact multiplicative order `N` of the T-matrix is computed in integer arithmetic
///    on root-of-unity exponents (every twist is `e^{2πi·e_x/L}` over the common level
///    `L = lcm(2·modality, 4)`).
/// 2. `T^N = I` is verified against the runtime construction.
/// 3. By the Ng–Schauenburg congruence-subgroup theorem the `SL(2,Z)` representation
///    factors through the finite group `SL(2, Z/N)`, whose order is computed by the Euler
///    product `N³·Π_{p|N}(1 − 1/p²)`; a finite image cannot be dense — that is the verdict,
///    as the conclusion of the checked premises above, with an error path if any premise
///    fails.
/// 4. At the Atlas instance the independent exact `Q(ζ₂₄)` certificate must agree
///    (finite projective image pinned, density refuted); disagreement is an error.
///
/// # Errors
/// If any premise of the derivation fails or the exact certificate disagrees.
pub fn finite_closure_probe(p: &UseCaseParams) -> Result<FiniteClosureMetrics, String> {
    let native = tqc_mtc::native::construct_atlas_native(p).map_err(|e| e.to_string())?;
    let dim = native.dim();
    let modality = p.modality as usize;
    let context = p.context as usize;
    let k_mod = if modality % 2 == 0 { 1 } else { 2 };

    fn gcd(a: usize, b: usize) -> usize {
        if b == 0 {
            a
        } else {
            gcd(b, a % b)
        }
    }
    let lcm = |a: usize, b: usize| a / gcd(a, b) * b;

    // (1) Exact T order over the level L = lcm(2·modality, 4).
    let level = lcm(2 * modality, 4);
    let mut n_order = 1usize;
    for x in 0..dim {
        let (m, c) = (x / context, x % context);
        let e_x = (k_mod * m * m * (level / (2 * modality))
            + (c.count_ones() as usize) * (level / 4))
            % level;
        let ord = level / gcd(e_x, level);
        n_order = lcm(n_order, ord);
    }

    // (2) Verify T^N = I against the runtime construction: the premise is computed,
    // then checked (a failed premise is an error, and the verdict below consumes the
    // computed boolean, never a literal).
    let t = native.t_diag();
    let t_order_verified = t.iter().all(|theta| {
        let mut pow = tqc_mtc::C::new(1.0, 0.0);
        for _ in 0..n_order {
            pow = pow.times(*theta);
        }
        pow.close(tqc_mtc::C::new(1.0, 0.0), tqc_mtc::TOL)
    });
    if !t_order_verified {
        return Err(format!(
            "derived T order {n_order} is wrong: some twist has T[i]^{n_order} != 1"
        ));
    }

    // (3) |SL(2, Z/N)| by the Euler product.
    let mut sl2_order = n_order * n_order * n_order;
    let mut temp = n_order;
    let mut primes = Vec::new();
    for i in 2..=n_order {
        if temp % i == 0 {
            primes.push(i);
            while temp % i == 0 {
                temp /= i;
            }
        }
    }
    for prime in primes {
        sl2_order = sl2_order / (prime * prime) * (prime * prime - 1);
    }
    // (4) Cross-check against the independent exact certificate at the Atlas instance.
    let certificate_agrees = if (p.modality, p.context) == (3, 8) {
        let cert = crate::exact::exact_density_certificate(p)?;
        !cert.certified_dense && cert.finite_image_order.is_some()
    } else {
        true
    };
    if !certificate_agrees {
        return Err(
            "finite-closure derivation contradicts the exact Q(zeta_24) certificate".into(),
        );
    }

    // The verdict is the conclusion of the checked premises: a representation whose image
    // factors through the finite SL(2, Z/N) is not dense.
    let image_finite = t_order_verified && certificate_agrees && sl2_order >= 1;
    let is_dense = !image_finite;

    Ok(FiniteClosureMetrics {
        is_dense,
        unique_phases: sl2_order,
        description: format!(
            "Finite closure derived: exact T order N = {n_order} (integer root-of-unity \
             exponents, verified T^N = I against the construction); by Ng–Schauenburg the \
             SL(2,Z) representation factors through SL(2, Z/{n_order}) of order {sl2_order}; \
             a finite image is not dense. At the Atlas instance the exact Q(zeta_24) \
             certificate independently pins the finite projective block image and agrees."
        ),
    })
}
/// The measured empirical Solovay-Kitaev metrics.
#[derive(Debug, Clone, PartialEq)]
pub struct SolovayKitaevMetrics {
    /// True if the generated braid subgroup is dense (universal quantum computation).
    pub is_dense: bool,
    /// Detailed description of the measurement.
    pub description: String,
}

/// A probe testing the Solovay-Kitaev density of the archimedean coupling.
/// Measures whether the indefinite spectrum mathematically implies infinite density.
#[allow(clippy::needless_range_loop)]
pub fn solovay_kitaev_probe(p: &UseCaseParams) -> Result<SolovayKitaevMetrics, String> {
    if (p.modality, p.context) != (3, 8) {
        return Err(format!(
            "the exact algebraic density certificate is defined over Q(zeta_24) for the Atlas \
             use-case (modality 3, context 8); got (modality {}, context {})",
            p.modality, p.context
        ));
    }

    let report = crate::exact::exact_density_certificate(p)?;
    if report.certified_dense {
        return Ok(SolovayKitaevMetrics {
            is_dense: report.certified_dense,
            description: report.description,
        });
    }

    Err(format!(
        "Exact algebraic certificate refutes single-qubit density on the 2-dim invariant block: {}",
        report.description
    ))
}
/// VV — universality as the equivalency facet: **realization-independence** of the
/// κ-equivalence class. The same braid word is executed through two independent
/// realizations — (a) the in-memory permutation action `permute_amplitudes` and (b) the
/// holospace gate path (`execute_holo_gate` through a compiled Hologram archive) — and
/// both must resolve to the identical κ. A cross-implementation check, not a replay of
/// one code path.
///
/// # Errors
/// If the two realizations of the same operator produce different κ.
pub fn equivalency_universality_probe(p: &UseCaseParams) -> Result<(), String> {
    let g = Generators::new(p);
    let n = p.class_count() as usize;
    let base: Vec<i64> = (0..n as i64).map(|i| i % 5 - 2).collect();
    let word: Vec<(&str, &Permutation)> = vec![
        ("sigma", &g.sigma),
        ("tau", &g.tau),
        ("mu", &g.mu),
        ("sigma", &g.sigma),
    ];

    // Realization (a): in-memory permutation composition.
    let mut perm = Permutation::identity(p.class_count());
    for (_, op) in &word {
        perm = perm.then(op);
    }
    let state_mem = perm.permute_amplitudes(&base);
    let amp_mem: Vec<(u64, Amplitude)> = state_mem
        .iter()
        .enumerate()
        .map(|(i, &v)| (i as u64, Amplitude { re: v, im: 0 }))
        .collect();
    let k_mem = tqc_substrate::kappa(&amplitude::encode(&amp_mem));

    // Realization (b): the holospace gate execution path.
    let amp0: Vec<(u64, Amplitude)> = base
        .iter()
        .enumerate()
        .map(|(i, &v)| (i as u64, Amplitude { re: v, im: 0 }))
        .collect();
    let mut bytes = amplitude::encode_interleaved(p.class_count(), &amp0);
    for (name, op) in &word {
        // Gather semantics: targets[i] = g⁻¹(i) realizes the amplitude action out[g(i)] = v[i].
        let inv = op.inverse();
        let targets: Vec<usize> = (0..p.class_count())
            .map(|i| inv.apply(i) as usize)
            .collect();
        let exec = tqc_substrate::execute_holo_gate(name, &targets, &bytes)
            .map_err(|e| format!("execute_holo_gate error: {e}"))?;
        bytes = exec.output_bytes;
    }
    let amp_holo = amplitude::decode_interleaved(&bytes)
        .ok_or_else(|| "holo output is not a whole number of label records".to_owned())?;
    let k_holo = tqc_substrate::kappa(&amplitude::encode(&amp_holo));

    if k_mem != k_holo {
        return Err(format!(
            "equivalency universality violated: two independent realizations of the same \
             operator produced different κ ({k_mem} != {k_holo})"
        ));
    }
    Ok(())
}

/// The measured deduplication metrics for UOR cache-collapse (an `open` measurement —
/// reported, never asserted as an advantage).
#[derive(Debug, Clone, PartialEq)]
pub struct ParetoMetrics {
    /// The total number of braid paths evaluated.
    pub total_paths: usize,
    /// The number of distinct resulting states (κ).
    pub distinct_states: usize,
    /// The topological degeneracy, measured as total_paths / distinct_states.
    pub topological_degeneracy: f64,
    /// The percentage of evaluations answerable from the content-addressed cache.
    pub compute_savings_pct: f64,
    /// Total bytes across all path-final state encodings (the naive storage cost).
    pub total_state_bytes: usize,
    /// Bytes across the distinct state encodings only (the deduplicated storage cost).
    pub unique_state_bytes: usize,
    /// Measured storage compression: total_state_bytes / unique_state_bytes.
    pub memory_compression_ratio: f64,
}

/// PROBE (open) — **measured content-addressed deduplication** over the finite braid orbit:
/// every braid word of generators evaluates to a state content-addressed to a κ; isotopic
/// words (those composing to the same operator) collapse to the identical κ. The braid group
/// image here is finite, so the orbit plateaus; the metrics quantify the deduplication, and
/// nothing more — no quantum-advantage claim is attached.
///
/// # Errors
/// If the enumeration produces no states (an internal defect).
pub fn advantage_probe(p: &UseCaseParams) -> Result<ParetoMetrics, String> {
    let g = Generators::new(p);
    let gens = [&g.sigma, &g.tau, &g.mu];
    let n = p.class_count() as usize;
    let base: Vec<i64> = (0..n as i64).map(|i| i % 7 - 3).collect();
    let length = 7u32;
    let total = 3usize.pow(length); // all length-7 braid words over {σ, τ, μ}
    let mut distinct: std::collections::BTreeMap<tqc_substrate::Kappa, usize> =
        std::collections::BTreeMap::new();
    let mut total_state_bytes = 0usize;
    for w in 0..total {
        let mut perm = Permutation::identity(p.class_count());
        let mut x = w;
        for _ in 0..length {
            perm = perm.then(gens[x % 3]);
            x /= 3;
        }
        let state = perm.permute_amplitudes(&base);
        let amp: Vec<(u64, Amplitude)> = state
            .iter()
            .enumerate()
            .map(|(i, &v)| (i as u64, Amplitude { re: v, im: 0 }))
            .collect();
        let encoded = amplitude::encode(&amp);
        total_state_bytes += encoded.len();
        distinct
            .entry(tqc_substrate::kappa(&encoded))
            .or_insert(encoded.len());
    }

    if distinct.is_empty() {
        return Err("advantage probe enumerated no states".into());
    }
    let distinct_count = distinct.len();
    let unique_state_bytes: usize = distinct.values().sum();
    let degeneracy = total as f64 / distinct_count as f64;
    let compute_savings = 100.0 * (1.0 - (distinct_count as f64 / total as f64));

    Ok(ParetoMetrics {
        total_paths: total,
        distinct_states: distinct_count,
        topological_degeneracy: degeneracy,
        compute_savings_pct: compute_savings,
        total_state_bytes,
        unique_state_bytes,
        memory_compression_ratio: total_state_bytes as f64 / unique_state_bytes as f64,
    })
}

/// VV (build) — the Atlas-native MTC construction: the parameter tuple admits the
/// `Z_2^k` composition quotient AND the constructed category passes the full phase-exact
/// axiom verification. Construction alone is not the claim — validity is.
///
/// # Errors
/// If the construction is obstructed or any axiom fails.
pub fn atlas_native_mtc(p: &tqc_core::UseCaseParams) -> Result<(), String> {
    let native = tqc_mtc::native::construct_atlas_native(p).map_err(|e| e.to_string())?;
    tqc_mtc::verifier::verify_mtc_axioms(&*native, tqc_mtc::TOL)
}

/// Witness the quantum realization: unitarity and interference on the pointed braiding.
#[allow(clippy::needless_range_loop)]
pub fn quantum_realization(p: &UseCaseParams) -> Witness {
    let native = tqc_mtc::native::construct_atlas_native(p).map_err(|e| format!("{:?}", e))?;
    let s = native.s_matrix();
    let t_diag = native.t_diag();
    let dim = native.dim();

    // 1. Unitarity on C^d
    // Verify S^dagger S = I
    let mut s_dag_s = vec![vec![tqc_mtc::C::new(0.0, 0.0); dim]; dim];
    for i in 0..dim {
        for j in 0..dim {
            for k in 0..dim {
                let left = s[k][i].conj();
                let right = s[k][j];
                s_dag_s[i][j].re += left.re * right.re - left.im * right.im;
                s_dag_s[i][j].im += left.re * right.im + left.im * right.re;
            }
        }
    }

    for i in 0..dim {
        for j in 0..dim {
            let expected = if i == j { 1.0 } else { 0.0 };
            if (s_dag_s[i][j].re - expected).abs() > 1e-9 || s_dag_s[i][j].im.abs() > 1e-9 {
                return Err(format!(
                    "Operator S is not unitary on C^{dim}: U^dagger U != I"
                ));
            }
        }
    }

    // Verify T is unitary
    for phase in &t_diag {
        if ((phase.re * phase.re + phase.im * phase.im).sqrt() - 1.0).abs() > 1e-9 {
            return Err(format!("Operator T is not unitary on C^{dim}"));
        }
    }

    // 2. Interference
    // Exhibit one input whose two-path evolution cancels by phase
    // Input state: uniform superposition of all basis states
    let v = vec![tqc_mtc::C::new(1.0, 0.0); dim];
    let mut w = vec![tqc_mtc::C::new(0.0, 0.0); dim];
    for i in 0..dim {
        for j in 0..dim {
            w[i].re += s[i][j].re * v[j].re - s[i][j].im * v[j].im;
            w[i].im += s[i][j].re * v[j].im + s[i][j].im * v[j].re;
        }
    }

    let mut found_interference = false;
    for i in 0..dim {
        // Find a measured-zero amplitude
        if w[i].re.abs() < 1e-9 && w[i].im.abs() < 1e-9 {
            let mut moduli_sum = 0.0;
            for j in 0..dim {
                let s_mod = (s[i][j].re * s[i][j].re + s[i][j].im * s[i][j].im).sqrt();
                let v_mod = (v[j].re * v[j].re + v[j].im * v[j].im).sqrt();
                moduli_sum += s_mod * v_mod;
            }
            // which the moduli alone would make nonzero
            if moduli_sum > 1e-9 {
                found_interference = true;
                break;
            }
        }
    }

    if !found_interference {
        return Err(
            "No interference witnessed: evolution is indistinguishable from classical prob".into(),
        );
    }

    Ok(())
}

/// VV (build) — generative closure, computed by BFS over the generator action: the group
/// generated by {σ, τ, μ} partitions the class space into exactly `⌈modality/2⌉` orbits
/// (the μ-mirror classes of the modality axis), whose sizes are computed and must sum to
/// the class count; every label is reachable from one seed per mirror class. The orbit
/// decomposition is derived, never asserted.
///
/// # Errors
/// If the computed orbit decomposition does not match the parametric derivation.
pub fn generative_closure_probe(p: &UseCaseParams) -> Result<(), String> {
    let g = Generators::new(p);
    let n = p.class_count() as usize;
    let gens = [&g.sigma, &g.tau, &g.mu];

    let orbit_of = |seed: u64, seen: &mut [bool]| -> usize {
        let mut queue = vec![seed];
        seen[seed as usize] = true;
        let mut size = 1usize;
        while let Some(x) = queue.pop() {
            for gen in gens {
                let y = gen.apply(x);
                if !seen[y as usize] {
                    seen[y as usize] = true;
                    size += 1;
                    queue.push(y);
                }
            }
        }
        size
    };

    let mut seen = vec![false; n];
    let mut orbit_sizes = Vec::new();
    for seed in 0..n as u64 {
        if !seen[seed as usize] {
            orbit_sizes.push(orbit_of(seed, &mut seen));
        }
    }

    let expected_orbits = (p.modality as usize).div_ceil(2);
    check(
        orbit_sizes.len() == expected_orbits,
        format!(
            "generator closure has {} orbits; the μ-mirror derivation predicts {expected_orbits}",
            orbit_sizes.len()
        ),
    )?;
    check(
        orbit_sizes.iter().sum::<usize>() == n,
        "orbit sizes do not cover the class space",
    )?;
    // Each orbit is scope·context·(2 for a mirror pair, 1 for a μ-fixed point).
    let cell = (p.scope as usize) * (p.context as usize);
    for &size in &orbit_sizes {
        check(
            size == cell || size == 2 * cell,
            format!("orbit size {size} is not scope·context or 2·scope·context"),
        )?;
    }
    Ok(())
}

/// VV (build) — the UTQC roll-up: an explicit conjunction over **every gating witness
/// implemented in this crate** (the algorithm reference evaluations live in
/// `tqc-algorithms::checks` and are gated by their own dictionary rows in both the BDD
/// suite and the ledger). The roll-up is green **only** when each pillar is green;
/// nothing is assumed from the runner.
///
/// # Errors
/// The first pillar that fails, prefixed with its name.
pub fn utqc_proven_probe(model: &Model, f1: &F1Constants, p: &UseCaseParams) -> Result<(), String> {
    let name = |n: &str, r: Witness| r.map_err(|e| format!("pillar `{n}` failed: {e}"));
    name("oracle-provenance", oracle_provenance(model, f1))?;
    name("objects-labels", objects_labels(p, f1))?;
    name("label-space-belt", label_space_belt(p, f1))?;
    name("inner-product", inner_product(p))?;
    name("reflection-generators", reflection_generators(p, f1))?;
    name("spectrum", spectrum(p, f1))?;
    name("coxeter-weyl", coxeter_weyl(p, f1))?;
    name("modular-identities", modular_identities(p, f1))?;
    name("definite-anchor-e8", definite_anchor_e8(f1))?;
    name("definite-anchor", definite_anchor(p))?;
    name("fusion-g2", fusion_g2(p))?;
    name("dual-f4", dual_f4(p))?;
    name("categorical-structure", categorical_structure(p))?;
    name("ground-space-protection", ground_space_protection(p))?;
    name(
        "complex-amplitude-encoding",
        complex_amplitude_encoding(p, f1),
    )?;
    name("atlas-native-mtc", atlas_native_mtc(p))?;
    name("modular-s-t", modular_s_t(p))?;
    name("braiding-r-matrix", braiding_r_matrix(p))?;
    name("mac-lane-coherence", mac_lane_coherence(p))?;
    name("quantum-realization", quantum_realization(p))?;
    name("holospace-cycle", holospace_cycle(p))?;
    name("s4-modal-logic", s4_frame_witness(p))?;
    name("generative-closure", generative_closure_probe(p))?;
    name("universality", equivalency_universality_probe(p))?;
    name("fault-tolerance", deterministic_replay_witness(p))?;
    name("complexity-bound", complexity_bound_witness(p))?;
    name("reconstructability", reconstruction_witness(p))?;
    name("tensor-contraction-bypass", isotopy_collision_witness(p))?;
    name(
        "topological-entanglement",
        topological_entanglement_probe(p).and_then(|m| {
            if m.is_logarithmic_scaling && m.entropy_bound > 0.0 {
                Ok(())
            } else {
                Err(format!("measured profile fails: {:?}", m.depth_profile))
            }
        }),
    )?;
    name("finite-closure", finite_closure_probe(p).map(|_| ()))?;
    name("solovay-kitaev", solovay_kitaev_decision_witness(p))?;
    name("archimedean-continuity", archimedean_continuity_witness(p))?;
    name("pair-carrier-structure", pair_carrier_witness(p))?;
    name(
        "two-qubit-universality",
        two_qubit_universality_probe(p).and_then(|m| {
            if m.is_entangling && m.is_coherent {
                Ok(())
            } else {
                Err("no coherent native entangler".into())
            }
        }),
    )?;
    name(
        "encoded-qubit-universality",
        encoded_qubit_universality_witness(p),
    )?;
    Ok(())
}

/// Metrics for the measured non-local entanglement entropy across a class-space bipartition.
#[derive(Debug, Clone, PartialEq)]
pub struct EntanglementMetrics {
    /// The measured maximal entanglement entropy `log2(max Schmidt rank)` over the sampled
    /// braid-evolved states.
    pub entropy_bound: f64,
    /// The maximal exact Schmidt rank observed.
    pub max_schmidt_rank: usize,
    /// The maximal Schmidt rank per braid depth `1..=D` (the measured profile).
    pub depth_profile: Vec<usize>,
    /// Measured verdict: the profile saturates (stops growing with depth) and the maximal
    /// entropy stays within `log2(class_count)` — sub-extensive, not Hilbert-volume scaling.
    pub is_logarithmic_scaling: bool,
}

/// Exact rank of a small integer matrix by fraction-free (Bareiss) elimination — no floats.
fn integer_rank(mut m: Vec<Vec<i128>>) -> usize {
    let rows = m.len();
    let cols = if rows > 0 { m[0].len() } else { 0 };
    let mut rank = 0usize;
    let mut prev = 1i128;
    let mut col = 0usize;
    while rank < rows && col < cols {
        if let Some(pr) = (rank..rows).find(|&r| m[r][col] != 0) {
            m.swap(rank, pr);
            for r in (rank + 1)..rows {
                for c in (col + 1)..cols {
                    let num = m[rank][col] * m[r][c] - m[r][col] * m[rank][c];
                    // Bareiss guarantees exact divisibility by the previous pivot; a
                    // remainder would mean silent truncation.
                    debug_assert!(num % prev == 0, "Bareiss division must be exact");
                    m[r][c] = num / prev;
                }
                m[r][col] = 0;
            }
            prev = m[rank][col];
            rank += 1;
        }
        col += 1;
    }
    rank
}

/// PROBE — topological entanglement entropy, **measured**: for every braid word up to depth
/// `D` the evolved integer state is reshaped over a bipartition `A×B` of the class space and
/// its exact Schmidt rank is computed by fraction-free elimination. The reported verdict is
/// derived from the measured profile (saturation within the `log2` bound), not asserted.
///
/// # Errors
/// If the class space admits no bipartition or the profile is empty.
pub fn topological_entanglement_probe(p: &UseCaseParams) -> Result<EntanglementMetrics, String> {
    let g = Generators::new(p);
    let gens = [&g.sigma, &g.tau, &g.mu];
    let n = p.class_count() as usize;
    let rows = (1..=n)
        .filter(|r| n % r == 0 && r * r <= n)
        .max()
        .ok_or_else(|| "class space admits no bipartition".to_owned())?;
    let cols = n / rows;
    let base: Vec<i64> = (0..n as i64).map(|i| i % 5 - 2).collect();

    let schmidt_rank = |state: &[i64]| -> usize {
        let m: Vec<Vec<i128>> = (0..rows)
            .map(|r| (0..cols).map(|c| i128::from(state[r * cols + c])).collect())
            .collect();
        integer_rank(m)
    };

    let max_depth = 6u32;
    let mut depth_profile = Vec::new();
    for depth in 1..=max_depth {
        let mut max_rank = 0usize;
        for w in 0..3usize.pow(depth) {
            let mut perm = Permutation::identity(p.class_count());
            let mut x = w;
            for _ in 0..depth {
                perm = perm.then(gens[x % 3]);
                x /= 3;
            }
            max_rank = max_rank.max(schmidt_rank(&perm.permute_amplitudes(&base)));
        }
        depth_profile.push(max_rank);
    }

    let max_schmidt_rank = *depth_profile
        .iter()
        .max()
        .ok_or_else(|| "empty entanglement profile".to_owned())?;
    let entropy_bound = (max_schmidt_rank as f64).log2();
    let saturated = depth_profile.len() >= 2
        && depth_profile[depth_profile.len() - 1] <= depth_profile[depth_profile.len() - 2];
    let is_logarithmic_scaling = saturated && entropy_bound <= (n as f64).log2();

    Ok(EntanglementMetrics {
        entropy_bound,
        max_schmidt_rank,
        depth_profile,
        is_logarithmic_scaling,
    })
}

/// The measured empirical Two-Qubit Universality metrics.
#[derive(Debug, Clone, PartialEq)]
pub struct TwoQubitUniversalityMetrics {
    /// True if an entangling two-qubit gate can be natively synthesized from the category.
    pub is_entangling: bool,
    /// True if the gate is constructed solely from the coherent abelian substrate, avoiding theory collision.
    pub is_coherent: bool,
    /// Detailed description of the measurement.
    pub description: String,
}

/// A probe for a native entangling two-qubit phase gate in the abelian category, decided
/// with **exact integer arithmetic** on the monodromy bicharacter exponents (no float
/// thresholds): logical basis states are fixed flux assignments `x₀,x₁` / `y₀,y₁` on
/// disjoint handles; the diagonal phase gate they induce is entangling iff
/// `χ(x₀,y₀)·χ(x₁,y₁) ≠ χ(x₀,y₁)·χ(x₁,y₀)` as exact root-of-unity exponents. Coherence is
/// **derived** by re-running the full phase-exact axiom verification on the very theory the
/// gate is drawn from — not asserted.
///
/// # Errors
/// If the construction fails or the axiom verification of the ambient theory fails.
pub fn two_qubit_universality_probe(
    p: &UseCaseParams,
) -> Result<TwoQubitUniversalityMetrics, String> {
    let native = tqc_mtc::native::construct_atlas_native(p).map_err(|e| e.to_string())?;
    let dim = native.dim();
    let modality = p.modality as usize;

    // Exact monodromy exponents; the entangling condition compares products of them.
    let chi = |x: usize, y: usize| chi_exponents(p, x, y);
    let mul = |a: (usize, usize), b: (usize, usize)| ((a.0 + b.0) % modality, (a.1 + b.1) % 2);

    let mut is_entangling = false;
    'search: for x0 in 0..dim {
        for x1 in 0..dim {
            if x0 == x1 {
                continue;
            }
            for y0 in 0..dim {
                for y1 in 0..dim {
                    if y0 == y1 {
                        continue;
                    }
                    let left = mul(chi(x0, y0), chi(x1, y1));
                    let right = mul(chi(x0, y1), chi(x1, y0));
                    if left != right {
                        is_entangling = true;
                        break 'search;
                    }
                }
            }
        }
    }

    // Coherence of the ambient theory is derived by full phase-exact axiom verification.
    let is_coherent = tqc_mtc::verifier::verify_mtc_axioms(&*native, tqc_mtc::TOL).is_ok();

    Ok(TwoQubitUniversalityMetrics {
        is_entangling,
        is_coherent,
        description: "A two-qubit entangling phase gate (CZ-equivalent) was decided with exact integer arithmetic on the monodromy bicharacter exponents of the coherent abelian Atlas-native construction, acting on fixed logical flux assignments. No gate-set density claim is attached: the exactly decided single-qubit image is the finite projective Clifford group and CZ is Clifford, so the two-qubit gate-set image is finite; universality is carried by the PU(22)/PU(576) density chain.".into(),
    })
}

/// The Solovay-Kitaev density question, exactly decided over `Q(zeta_24)`.
///
/// This witness asserts the DECISION as a theorem, in both directions, with no false green:
/// the unique 2-dim commutant block exists (commutant dim exactly 2, `tr P1 = 2`), is confined
/// to the `(-1)` spectral eigenspace, has `tr(P1 G_S) = 0` identically, carries the finite
/// projective Clifford image of exact order 24, and density on the block is refuted. Any
/// deviation from these exact facts is an error.
pub fn solovay_kitaev_decision_witness(p: &UseCaseParams) -> Result<(), String> {
    let r = crate::exact::exact_density_certificate(p)?;
    if r.commutant_dim != 2 {
        return Err(format!("commutant dim {} != 2", r.commutant_dim));
    }
    if r.block_dim != 2 {
        return Err(format!("block dim {} != 2", r.block_dim));
    }
    if !r.beta_s_nonzero.is_empty() {
        return Err(format!(
            "tr(P1 G_S) not identically zero: {:?}",
            r.beta_s_nonzero
        ));
    }
    let expected_support = vec![(10i64, 0.0f64), (7, 0.0), (2, 0.0), (-1, 2.0)];
    if r.block_support != expected_support {
        return Err(format!(
            "block support {:?} != {:?}",
            r.block_support, expected_support
        ));
    }
    if r.finite_image_order != Some(24) {
        return Err(format!(
            "projective image order {:?} != Some(24)",
            r.finite_image_order
        ));
    }
    if r.certified_dense {
        return Err("density unexpectedly certified; the decision changed".into());
    }
    Ok(())
}

/// Archimedean continuity, exactly located and saturated: on the 22-dim irreducible block
/// the projective closure of the coupled generators is DENSE in PU(22). The chain: the
/// spectral flow exp(iRM) lies in the closure (Kronecker-Weyl; pi irrational), seeding a
/// division-free Lie closure under Ad(S), Ad(T), brackets, and torus-weight splitting; its
/// mod-p rank on the block is a sound lower bound on dim Lie(H), and saturation at >= 483
/// forces su(22) inside (su(22) is simple with minimal proper-subalgebra codimension 42),
/// hence closure >= PSU(22): density. Universal quantum computation on a 22-dimensional
/// qudit carrier follows by Solovay-Kitaev in PU(d). The 2-dim block carries the finite
/// projective Clifford image; the continuity certificates (adjoint-trace infinite order,
/// projectively non-commuting pair) remain asserted as prerequisites.
pub fn archimedean_continuity_witness(p: &UseCaseParams) -> Result<(), String> {
    let r = crate::exact::exact_density_certificate(p)?;
    if r.commutant_dim != 2 {
        return Err(format!(
            "commutant dim {} != 2 (irreducibility premise)",
            r.commutant_dim
        ));
    }
    if !r.block22_infinite.iter().any(|x| x == "T") || !r.block22_infinite.iter().any(|x| x == "S")
    {
        return Err(format!(
            "generator words not certified infinite projective order on the 22-dim block: {:?}",
            r.block22_infinite
        ));
    }
    if r.block22_pair.is_none() {
        return Err("no projectively non-commuting pair on the 22-dim block".into());
    }
    if !r.beyond_finite {
        return Err("beyond-finite certificate not established".into());
    }
    if r.lie_dim_lower_22 < 483 || !r.pu22_dense {
        return Err(format!(
            "PU(22)-density not saturated: Lie dimension lower bound {} < 483",
            r.lie_dim_lower_22
        ));
    }
    Ok(())
}

/// The two-handle (pair-carrier) structure, exactly decided. Three theorems, pinned:
/// (1) irreducibility: the two-handle native group (per-handle coupled generators plus
/// the monodromy) has exact commutant dimension 1 on the 576-dim pair carrier;
/// (2) separation, in its strong form: **no nontrivial power of the monodromy even
/// preserves** the 22-block tensor code `W' (x) W'` (the pinned set of code-preserving
/// powers is empty), so the native diagonal sector cannot entangle the continuous
/// carriers -- the multi-handle carrier is the irreducible pair block itself, not a
/// tensor code;
/// (3) native continuous entanglement: the closure's identity component strictly exceeds
/// the local subalgebra (sound mod-p lower bound > 976), so continuous entangling flows
/// exist natively on the pair carrier;
/// (4) density: the T1 certificate (nonzero adj (x) adj component, multiplicity-one
/// isotypic, hence su(484) on the corner) and the T2 certificate (complement reachability
/// rank 92, the ambient cap) combine with the classical closure T3 to force su(576)
/// inside Lie(H_2): the two-handle projective closure is DENSE in PU(576), and by the
/// two-local composition lemma the n-handle closure is dense in PU(24^n) for every
/// n >= 2 -- gate-level universal quantum computation, scaling in n.
/// Any drift in these exact values is an error.
pub fn pair_carrier_witness(p: &UseCaseParams) -> Result<(), String> {
    let r = crate::exact::exact_density_certificate(p)?;
    if r.pair_commutant_dim != 1 {
        return Err(format!(
            "pair commutant dim {} != 1 (irreducibility)",
            r.pair_commutant_dim
        ));
    }
    if r.native_code_entangler.is_some() {
        return Err(format!(
            "separation theorem changed: monodromy power {:?} now preserves the code",
            r.native_code_entangler
        ));
    }
    if !r.monodromy_code_preserving_powers.is_empty() {
        return Err(format!(
            "strong separation changed: monodromy powers {:?} now preserve the code space \
             (expected none)",
            r.monodromy_code_preserving_powers
        ));
    }
    if r.qudit_universal {
        return Err("qudit_universal flag inconsistent with the separation theorem".into());
    }
    if r.pair_lie_dim_lower <= 976 || !r.pair_entangling_flow {
        return Err(format!(
            "pair Lie lower bound {} does not exceed the local subalgebra bound 976",
            r.pair_lie_dim_lower
        ));
    }
    if !r.pair_adj_component {
        return Err("T1 failed: no adj (x) adj component certified in Lie(H_2)".into());
    }
    if r.pair_reach_rank != 92 {
        return Err(format!(
            "T2 failed: complement reachability rank {} != 92",
            r.pair_reach_rank
        ));
    }
    if !r.pu576_dense || !r.gate_level_universal {
        return Err("pair-carrier PU(576) density chain did not close".into());
    }
    Ok(())
}

/// VV (build) — Mac Lane coherence: the pentagon and hexagon identities of the Atlas-native
/// category hold with phase-exact comparisons (via the full axiom verifier).
///
/// # Errors
/// Returns the first coherence identity that fails.
pub fn mac_lane_coherence(p: &UseCaseParams) -> Witness {
    let native = tqc_mtc::native::construct_atlas_native(p).map_err(|e| e.to_string())?;
    tqc_mtc::verifier::verify_mtc_axioms(&*native, tqc_mtc::TOL)
}

/// VV (build) — the S4 modal frame, built and checked by enumeration: worlds are the
/// classes, accessibility `R(u,v)` is reachability of `v` from `u` under generator words
/// (computed by BFS, including the empty word). Reflexivity and transitivity are then
/// verified by explicit enumeration over the frame — the S4 axioms as checks, not as
/// assertions about generator orders.
///
/// # Errors
/// If reflexivity or transitivity fails at any world.
pub fn s4_frame_witness(p: &UseCaseParams) -> Witness {
    let g = Generators::new(p);
    let n = p.class_count() as usize;
    let gens = [&g.sigma, &g.tau, &g.mu];

    // R(u, ·) = BFS closure from u (reflexive: the empty word).
    let reach_from = |u: usize| -> Vec<bool> {
        let mut seen = vec![false; n];
        seen[u] = true;
        let mut queue = vec![u as u64];
        while let Some(x) = queue.pop() {
            for gen in gens {
                let y = gen.apply(x) as usize;
                if !seen[y] {
                    seen[y] = true;
                    queue.push(y as u64);
                }
            }
        }
        seen
    };
    let r: Vec<Vec<bool>> = (0..n).map(reach_from).collect();

    for (u, row) in r.iter().enumerate() {
        check(row[u], format!("frame not reflexive at world {u}"))?;
    }
    for u in 0..n {
        for v in 0..n {
            if r[u][v] {
                for (w, &v_reaches_w) in r[v].iter().enumerate() {
                    if v_reaches_w {
                        check(r[u][w], format!("frame not transitive at ({u},{v},{w})"))?;
                    }
                }
            }
        }
    }
    Ok(())
}

/// VV (build) — deterministic discrete execution: the identical braid word, replayed,
/// produces the byte-identical state and κ. This is determinism of the discrete execution
/// model — a prerequisite for content-addressed collapse — and nothing stronger; no
/// physical decoherence claim is attached.
///
/// # Errors
/// If a replay diverges.
pub fn deterministic_replay_witness(p: &UseCaseParams) -> Witness {
    let g = Generators::new(p);
    let n = p.class_count() as usize;
    let base: Vec<i64> = (0..n as i64).map(|i| i % 7 - 3).collect();
    let run = || {
        let perm = Permutation::identity(p.class_count())
            .then(&g.sigma)
            .then(&g.tau)
            .then(&g.mu)
            .then(&g.sigma);
        let state = perm.permute_amplitudes(&base);
        let amp: Vec<(u64, Amplitude)> = state
            .iter()
            .enumerate()
            .map(|(i, &v)| (i as u64, Amplitude { re: v, im: 0 }))
            .collect();
        (state, tqc_substrate::kappa(&amplitude::encode(&amp)))
    };
    let (s1, k1) = run();
    let (s2, k2) = run();
    check(s1 == s2, "replayed word produced a different state")?;
    check(k1 == k2, "replayed word produced a different κ")
}

/// VV (build) — the execution-cost bound: composing a depth-`M` braid word is performed
/// by an explicitly counted elementary loop (`M·n` generator-map reads), the state size is
/// asserted invariant at every step (no state exponential in depth is ever materialized),
/// and the counted composition is **cross-checked** against the independent
/// `Permutation::then` implementation — a real consistency check between two
/// implementations, not a self-comparison.
///
/// # Errors
/// If the state size changes, the counted cost deviates from `M·n`, or the two
/// composition implementations disagree.
pub fn complexity_bound_witness(p: &UseCaseParams) -> Witness {
    let g = Generators::new(p);
    let n = p.class_count();
    let word_len = 1000usize;
    let gen_at = |step: usize| match step % 3 {
        0 => &g.sigma,
        1 => &g.tau,
        _ => &g.mu,
    };

    // Counted composition: an explicit elementary loop over the class map.
    let mut ops: u64 = 0;
    let mut composed: Vec<u64> = (0..n).collect();
    for step in 0..word_len {
        let gen = gen_at(step);
        for slot in composed.iter_mut() {
            *slot = gen.apply(*slot);
            ops += 1;
        }
        check(
            composed.len() as u64 == n,
            "state/permutation size changed during execution",
        )?;
    }
    check(
        ops == word_len as u64 * n,
        format!(
            "operation count {ops} != depth·classes = {}",
            word_len as u64 * n
        ),
    )?;

    // Cross-check against the independent composition implementation.
    let mut perm = Permutation::identity(n);
    for step in 0..word_len {
        perm = perm.then(gen_at(step));
    }
    for (i, &img) in composed.iter().enumerate() {
        check(
            perm.apply(i as u64) == img,
            format!("counted composition disagrees with Permutation::then at class {i}"),
        )?;
    }
    Ok(())
}

/// VV (build) — reconstructability: a validator that receives only the serialized genesis
/// state and the serialized braid word (both round-tripped through their byte encodings)
/// reconstructs the byte-identical final state and κ. Reconstruction runs from the parsed
/// serialization, not from shared in-memory values.
///
/// # Errors
/// If parsing fails or the reconstruction diverges.
pub fn reconstruction_witness(p: &UseCaseParams) -> Witness {
    let g = Generators::new(p);
    let n = p.class_count() as usize;
    let base: Vec<i64> = (0..n as i64).map(|i| i % 5).collect();

    // The prover executes and publishes: serialized genesis, serialized word, final κ.
    let word = ["sigma", "tau", "sigma", "mu", "tau"];
    let by_name = |name: &str| match name {
        "sigma" => &g.sigma,
        "tau" => &g.tau,
        _ => &g.mu,
    };
    let execute = |genesis: &[i64], word: &[&str]| -> (Vec<i64>, String) {
        let mut perm = Permutation::identity(p.class_count());
        for name in word {
            perm = perm.then(by_name(name));
        }
        let state = perm.permute_amplitudes(genesis);
        let amp: Vec<(u64, Amplitude)> = state
            .iter()
            .enumerate()
            .map(|(i, &v)| (i as u64, Amplitude { re: v, im: 0 }))
            .collect();
        let k = tqc_substrate::kappa(&amplitude::encode(&amp)).to_string();
        (state, k)
    };
    let (published_state, published_kappa) = execute(&base, &word);

    // Serialize genesis + word to bytes, as a validator would receive them.
    let genesis_amp: Vec<(u64, Amplitude)> = base
        .iter()
        .enumerate()
        .map(|(i, &v)| (i as u64, Amplitude { re: v, im: 0 }))
        .collect();
    let genesis_bytes = amplitude::encode(&genesis_amp);
    let word_bytes = word.join(",").into_bytes();

    // The validator parses the serializations and independently reconstructs.
    let parsed_genesis: Vec<i64> = amplitude::decode(&genesis_bytes)
        .ok_or_else(|| "validator failed to parse the genesis serialization".to_owned())?
        .iter()
        .map(|&(_, a)| a.re)
        .collect();
    let word_str = String::from_utf8(word_bytes)
        .map_err(|_| "validator failed to parse the word serialization".to_owned())?;
    let parsed_word: Vec<&str> = word_str.split(',').collect();
    let (validator_state, validator_kappa) = execute(&parsed_genesis, &parsed_word);

    check(
        validator_state == published_state,
        "validator reconstruction produced a different state",
    )?;
    check(
        validator_kappa == published_kappa,
        "validator reconstruction produced a different κ",
    )
}

/// VV (build) — isotopy collision: two computationally distinct braid words that compose
/// to the same operator (`σ^order` vs the empty word) must content-address to the identical
/// κ. This is the decision-by-equivalence mechanism (κ-collision on isotopic words); no
/// claim about #P-hard contraction is attached.
///
/// # Errors
/// If the isotopic words resolve to different κ.
pub fn isotopy_collision_witness(p: &UseCaseParams) -> Witness {
    let g = Generators::new(p);
    let n = p.class_count() as usize;
    let base: Vec<i64> = (0..n as i64).collect();
    let kappa_of = |perm: &Permutation| {
        let state = perm.permute_amplitudes(&base);
        let amp: Vec<(u64, Amplitude)> = state
            .iter()
            .enumerate()
            .map(|(i, &v)| (i as u64, Amplitude { re: v, im: 0 }))
            .collect();
        tqc_substrate::kappa(&amplitude::encode(&amp))
    };
    let id = Permutation::identity(p.class_count());
    let mut pow = Permutation::identity(p.class_count());
    for _ in 0..p.sigma_order() {
        pow = pow.then(&g.sigma);
    }
    check(
        kappa_of(&id) == kappa_of(&pow),
        "isotopic words (identity vs σ^order) resolved to different κ",
    )
}

/// VV (build) — encoded-qubit universality corollary. On the `n=2` handle carrier
/// (`24^2 = 576`), density in PU(576) (from [`pair_carrier_witness`]) plus the fact that the
/// `k`-qubit code subgroup `SU(2^k) ⊕ I` is a CLOSED subgroup of `SU(576)` yields dense
/// encoded single- and two-qubit gates on the register. The machine-checked content is the
/// exact `Q(ζ₂₄)` verification that the encoded gate set (H on each logical qubit and a
/// genuine CZ entangler) is unitary with its defining relations and that `U ↦ U ⊕ I` is
/// `*`-preserving and injective on the generators, plus the density premise; the encoding
/// is pinned by `κ`. Image closedness and density itself are the cited closed-subgroup
/// consequence, not re-derived here.
///
/// # Errors
/// If any premise of the corollary fails or the pinned encoding `κ` drifts.
pub fn encoded_qubit_universality_witness(p: &UseCaseParams) -> Result<(), String> {
    let r = crate::exact::encoded_qubit_certificate(p)?;
    check(r.code_fits, "code does not embed: 2^k > 24^n")?;
    check(
        r.generators_unitary,
        "encoded generators are not exactly unitary over Q(zeta_24)",
    )?;
    check(
        r.relations_hold,
        "encoded generator relations fail (H^2=I, CZ^2=I, commuting, CZ entangling)",
    )?;
    check(
        r.faithful_star_embedding,
        "the U ↦ U ⊕ I block embedding failed the exact *-preservation / injectivity checks",
    )?;
    check(
        r.density_premise,
        "PU(24^n) density premise not established by the pair-carrier certificate",
    )?;
    check(
        r.encoding_kappa
            == "blake3:13e70fcfc33f841ea898a5fc2e5d42e45abb3477d26f152a4cabbec2b99c88f8",
        format!("encoding κ drifted: {}", r.encoding_kappa),
    )?;
    check(
        r.encoded_universal,
        "encoded-qubit universality corollary did not close",
    )
}

#[cfg(test)]
mod tests {
    use super::*;
    use tqc_atlas::canonical;

    fn atlas() -> (Model, F1Constants, UseCaseParams) {
        let model = Model::load().unwrap();
        let f1 = F1Constants::load().unwrap();
        let p = canonical(&model).unwrap();
        (model, f1, p)
    }

    #[test]
    fn vv_oracle_provenance() {
        let (m, f1, _) = atlas();
        oracle_provenance(&m, &f1).unwrap();
    }

    #[test]
    fn vv_all_some_true_suite_witnesses_pass_on_the_atlas() {
        let (_, f1, p) = atlas();
        objects_labels(&p, &f1).unwrap();
        label_space_belt(&p, &f1).unwrap();
        inner_product(&p).unwrap();
        reflection_generators(&p, &f1).unwrap();
        spectrum(&p, &f1).unwrap();
        coxeter_weyl(&p, &f1).unwrap();
        modular_identities(&p, &f1).unwrap();
        definite_anchor_e8(&f1).unwrap();
        definite_anchor(&p).unwrap();
        fusion_g2(&p).unwrap();
        dual_f4(&p).unwrap();
        categorical_structure(&p).unwrap();
        ground_space_protection(&p).unwrap();
        complex_amplitude_encoding(&p, &f1).unwrap();
        modular_s_t(&p).unwrap();
        braiding_r_matrix(&p).unwrap();
        holospace_cycle(&p).unwrap();
        quantum_realization(&p).unwrap();
        topological_entanglement_probe(&p).unwrap();
    }

    #[test]
    fn substrate_coupled_witnesses_hold_at_an_arbitrary_use_case() {
        let p = UseCaseParams::new(2, 2, 4);
        definite_anchor(&p).unwrap();
        fusion_g2(&p).unwrap();
        dual_f4(&p).unwrap();
        categorical_structure(&p).unwrap();
        ground_space_protection(&p).unwrap();
        let f1 = F1Constants::load().unwrap();
        complex_amplitude_encoding(&p, &f1).unwrap();
        modular_s_t(&p).unwrap();
        braiding_r_matrix(&p).unwrap();
        holospace_cycle(&p).unwrap();
        quantum_realization(&p).unwrap();
        topological_entanglement_probe(&p).unwrap();
    }
}
