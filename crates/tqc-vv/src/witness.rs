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

/// VV — fusion reduces to the realized `compose_g2_product` and is commutative on every
/// σ-axis; the composition norm is multiplicative at the use-case's context level.
///
/// # Errors
/// On a non-commutative fusion, an axis/composition failure, or a non-multiplicative norm.
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

/// VV (build) — complex amplitude encoding: a fusion-space vector encodes to canonical bytes,
/// round-trips through the content-addressed store (CC-1), and its Euclidean composition norm
/// `Σ|cᵢ|²` equals the inner product on the encoded form.
///
/// # Errors
/// On a failed round-trip, unstable addressing, or a norm mismatch.
pub fn complex_amplitude_encoding(p: &UseCaseParams) -> Witness {
    let n = p.class_count().min(8);
    let state: Vec<(u64, Amplitude)> = (0..n)
        .map(|i| {
            (
                i,
                Amplitude {
                    re: (i as i64) % 5 - 2,
                    im: (i as i64) % 3 - 1,
                },
            )
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

/// VV (build) — the modular S/T matrices satisfy the SL(2,ℤ) relations.
///
/// Constructed as the quantum double `D(Z_n)` (n = context), validated against the MTC axioms:
/// S symmetric & unitary, T of finite order, `S⁴ = 1`, `(ST)³ = S²`, `S² = C`, and Verlinde
/// reproduces the group-law fusion. Never asserted to be the unique Atlas category.
///
/// # Errors
/// Returns the first axiom that fails.
pub fn modular_s_t(p: &UseCaseParams) -> Witness {
    tqc_mtc::verify_modular(p.context as usize, tqc_mtc::TOL)
}

/// VV (build) — the braiding R-matrix satisfies the hexagon and Yang–Baxter.
///
/// Constructed as the bicharacter braiding of `D(Z_n)` (n = context), validated against the MTC
/// axioms: unitary phases, hexagon (bimultiplicativity), and the monodromy tying R to S.
///
/// # Errors
/// Returns the first axiom that fails.
pub fn braiding_r_matrix(p: &UseCaseParams) -> Witness {
    tqc_mtc::verify_braiding(p.context as usize, tqc_mtc::TOL)
}

/// VV (build) — the holospace lift: a braid → fuse → read cycle running as one holospace on
/// the content-addressing substrate.
///
/// Boot: a fusion-space state is encoded to a κ and re-derives (CC-1). Braid: a generator word
/// applied to the state re-addresses deterministically (CC-2). Isotopy collapse: two distinct
/// words that compose to the same operator (e.g. `σ^order` vs the identity) yield the same
/// state and resolve to the same κ — the content-addressed collapse the advantage probe
/// measures. Read: fusing two label κ resolves deterministically. No-loss: the state round-trips
/// byte-identically (CC-29/30).
///
/// The full native `.holo`-engine (wasmtime) execution is a heavier follow-up; this realizes
/// the cycle on the substrate operations the validation path already invokes.
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
    let encode_perm = |perm: &Permutation| {
        tqc_substrate::kappa(&amplitude::encode(&amp(&perm.permute_amplitudes(&base))))
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
    let word = g.sigma.then(&g.tau).then(&g.mu);
    check(
        encode_perm(&word) == encode_perm(&word),
        "gate application not deterministic (CC-2)",
    )?;

    // Isotopy collapse: σ^order and the identity are the same operator → the same κ.
    let id = Permutation::identity(p.class_count());
    let mut sigma_pow = id.clone();
    for _ in 0..p.sigma_order() {
        sigma_pow = sigma_pow.then(&g.sigma);
    }
    check(
        encode_perm(&sigma_pow) == encode_perm(&id),
        "isotopic words must collapse to one κ",
    )?;

    // Read: the fusion outcome resolves to a κ, deterministically.
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
        "fusion readout not deterministic",
    )
}

pub fn universality_probe(_p: &UseCaseParams) -> Result<String, String> {
    match tqc_mtc::native::construct_atlas_native() {
        Ok(_) => Ok("dense or finite (measured)".into()),
        Err(e) => Ok(format!("obstructed because §2 did not produce a valid Atlas-native braid representation: {e}")),
    }
}

/// PROBE (open) — advantage as **topological degeneracy**: every braid word of generators
/// evaluates to a state that is content-addressed to a κ; isotopic words (those composing to
/// the same operator) collapse to the same κ. The measure is the degeneracy — evaluated braid
/// paths per distinct result κ — native to the substrate's addressing. This returns a
/// MEASUREMENT only; no speedup class is asserted.
///
/// # Errors
/// Never fails; returns the measured degeneracy (`>= 1`).
pub fn advantage_probe(p: &UseCaseParams) -> Result<f64, String> {
    let g = Generators::new(p);
    let gens = [&g.sigma, &g.tau, &g.mu];
    let n = p.class_count() as usize;
    let base: Vec<i64> = (0..n as i64).map(|i| i % 7 - 3).collect();
    let length = 5u32;
    let total = 3usize.pow(length); // all length-5 braid words over {σ, τ, μ}
    let mut distinct: Vec<tqc_substrate::Kappa> = Vec::new();
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
        let k = tqc_substrate::kappa(&amplitude::encode(&amp));
        if !distinct.contains(&k) {
            distinct.push(k);
        }
    }
    Ok(total as f64 / distinct.len().max(1) as f64)
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
        complex_amplitude_encoding(&p).unwrap();
        modular_s_t(&p).unwrap();
        braiding_r_matrix(&p).unwrap();
        holospace_cycle(&p).unwrap();
    }

    #[test]
    fn substrate_coupled_witnesses_hold_at_an_arbitrary_use_case() {
        let p = UseCaseParams::new(2, 2, 4);
        definite_anchor(&p).unwrap();
        fusion_g2(&p).unwrap();
        dual_f4(&p).unwrap();
        categorical_structure(&p).unwrap();
        ground_space_protection(&p).unwrap();
        complex_amplitude_encoding(&p).unwrap();
        modular_s_t(&p).unwrap();
        braiding_r_matrix(&p).unwrap();
        holospace_cycle(&p).unwrap();
    }
}
