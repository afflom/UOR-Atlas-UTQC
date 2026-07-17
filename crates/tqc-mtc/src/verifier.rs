use crate::{close_mat, identity, is_symmetric, is_unitary, mat_pow, zeros, Matrix, C};

/// Data for a generic Modular Tensor Category (MTC).
pub trait ModularData {
    /// The number of simple objects.
    fn dim(&self) -> usize;
    /// The modular S matrix.
    fn s_matrix(&self) -> Matrix;
    /// The diagonal topological spins T.
    fn t_diag(&self) -> Vec<C>;
    /// The charge-conjugation permutation matrix.
    fn charge_conjugation(&self) -> Matrix;
    /// Fusion coefficient N_{ij}^k.
    fn n_ijk(&self, i: usize, j: usize, k: usize) -> f64;
    /// The F-symbol (associator).
    fn f_symbol(&self, i: usize, j: usize, k: usize, l: usize, m: usize, n: usize) -> C;
    /// The R-symbol (braiding).
    fn r_symbol(&self, i: usize, j: usize, k: usize) -> C;
}

/// The search cap for the multiplicative order of a twist. Vafa's theorem guarantees the
/// twists of a genuine MTC are roots of unity; a data set whose twists exceed this order
/// is rejected rather than silently accepted.
const TWIST_ORDER_CAP: usize = 1 << 16;

/// Verify the universal MTC axioms for the supplied modular data.
///
/// **Implemented checks** (each compares full complex values within `tol` — no
/// magnitude-only comparisons):
/// - `S` is symmetric and unitary; `S⁴ = I`; `S² = C` (charge conjugation).
/// - Every twist `θᵢ` is a phase of finite multiplicative order (Vafa's theorem).
/// - `(ST)³ = p⁺ S²` where the anomaly `p⁺` is a phase, and `p⁺` equals the Gauss sum
///   `(Σᵢ dᵢ² θᵢ)/D` — it is cross-derived, not merely read off `(ST)³`.
/// - Fusion coefficients `N_{ij}^k` are non-negative integers with a unique identity object.
/// - The Verlinde formula reproduces `N_{ij}^k` exactly (real part = N, imaginary part = 0).
/// - `R`-symbols are unit phases on every admissible fusion channel.
/// - The balancing (ribbon) equation `R^k_{ij} R^k_{ji} = θ_k θᵢ⁻¹ θⱼ⁻¹` holds exactly.
/// - The pentagon equation holds exactly for the `F`-symbols.
/// - Both hexagon equations hold exactly for `F` and `R` (by Joyal–Street coherence,
///   pentagon + hexagons imply the Yang–Baxter relation for the induced braiding).
/// - The monodromy–S relation `D·S_{ab} = Σ_c N_{āb}^c d_c θ_c θ_ā⁻¹ θ_b⁻¹` ties the
///   braiding data back to the modular data.
///
/// # Errors
/// Returns a description of the first axiom that fails within `tol`.
#[allow(clippy::needless_range_loop, clippy::too_many_lines)]
pub fn verify_mtc_axioms(m: &dyn ModularData, tol: f64) -> Result<(), String> {
    let dim = m.dim();
    let s = m.s_matrix();
    let t = m.t_diag();

    if !is_symmetric(&s, tol) {
        return Err("S is not symmetric".into());
    }
    if !is_unitary(&s, tol) {
        return Err("S is not unitary".into());
    }

    // Twists are phases of finite multiplicative order (Vafa).
    for (i, theta) in t.iter().enumerate() {
        if (theta.abs2() - 1.0).abs() > tol {
            return Err(format!("T[{i}] is not a phase"));
        }
        let mut pow = *theta;
        let mut order = 1usize;
        while !pow.close(C::new(1.0, 0.0), tol) {
            pow = pow.times(*theta);
            order += 1;
            if order > TWIST_ORDER_CAP {
                return Err(format!(
                    "T[{i}] has no finite order <= {TWIST_ORDER_CAP} (Vafa violated)"
                ));
            }
        }
    }

    // Fusion coefficients are non-negative integers.
    for i in 0..dim {
        for j in 0..dim {
            for k in 0..dim {
                let n_val = m.n_ijk(i, j, k);
                if n_val < -tol || (n_val - n_val.round()).abs() > tol {
                    return Err(format!(
                        "N_{{{i},{j}}}^{k} = {n_val} is not a non-negative integer"
                    ));
                }
            }
        }
    }

    // Unique identity object: N_{e,j}^k = δ_{jk}.
    let mut identity_obj = None;
    for i in 0..dim {
        let is_id = (0..dim).all(|j| {
            (0..dim).all(|k| {
                let expected = if j == k { 1.0 } else { 0.0 };
                (m.n_ijk(i, j, k) - expected).abs() <= tol
            })
        });
        if is_id {
            if let Some(prev) = identity_obj {
                return Err(format!("identity object is not unique: {prev} and {i}"));
            }
            identity_obj = Some(i);
        }
    }
    let e =
        identity_obj.ok_or_else(|| "No identity object found (N_{e,i}^j != δ_{ij})".to_string())?;

    // Quantum dimensions d_i = S_{ei}/S_{ee} and the total dimension D = 1/S_{ee}.
    let s_ee = s[e][e];
    if s_ee.abs2() < tol * tol {
        return Err("S_{ee} is zero".into());
    }
    if s_ee.im.abs() > tol {
        return Err("S_{ee} is not real".into());
    }
    let d_total = 1.0 / s_ee.re;
    let dims: Vec<f64> = (0..dim)
        .map(|i| {
            let d = s[e][i];
            d.re / s_ee.re
        })
        .collect();
    for i in 0..dim {
        if s[e][i].im.abs() > tol {
            return Err(format!("quantum dimension d_{i} is not real"));
        }
    }

    let s2 = mat_pow(&s, 2);
    if !close_mat(&mat_pow(&s, 4), &identity(dim), tol) {
        return Err("S^4 != I".into());
    }

    // (ST)^3 = p^+ S^2 where p^+ is the anomaly phase e^{2πi c/8}. The anomaly is
    // cross-derived from the Gauss sum (Σ d_a² θ_a)/D and must agree with (ST)³ S⁻².
    let mut st = zeros(dim);
    for (i, row) in st.iter_mut().enumerate() {
        for (j, cell) in row.iter_mut().enumerate() {
            *cell = s[i][j].times(t[j]);
        }
    }
    let st3 = mat_pow(&st, 3);
    let mut gauss = C::new(0.0, 0.0);
    for i in 0..dim {
        gauss = gauss.plus(t[i].scale(dims[i] * dims[i]));
    }
    let p_plus = gauss.scale(1.0 / d_total);
    if (p_plus.abs2() - 1.0).abs() > tol {
        return Err(format!(
            "Anomaly phase p^+ (Gauss sum / D) has invalid magnitude: {}",
            p_plus.abs2().sqrt()
        ));
    }
    let mut scaled_s2 = zeros(dim);
    for (i, row) in scaled_s2.iter_mut().enumerate() {
        for (j, cell) in row.iter_mut().enumerate() {
            *cell = s2[i][j].times(p_plus);
        }
    }
    if !close_mat(&st3, &scaled_s2, tol) {
        return Err("(ST)^3 != p^+ S^2 for the Gauss-sum anomaly p^+".into());
    }

    if !close_mat(&s2, &m.charge_conjugation(), tol) {
        return Err("S^2 != charge conjugation".into());
    }

    // Verlinde formula, exact: the sum must equal N_{ij}^k (already known to be a
    // non-negative integer) with vanishing imaginary part.
    for i in 0..dim {
        for j in 0..dim {
            for k in 0..dim {
                let n_val = m.n_ijk(i, j, k);
                let mut sum = C::new(0.0, 0.0);
                for l in 0..dim {
                    let s_0l = s[e][l];
                    if s_0l.abs2() < 1e-14 {
                        return Err(format!("S_{{0,{l}}} is zero, obstructing Verlinde"));
                    }
                    let inv = C::new(s_0l.re / s_0l.abs2(), -s_0l.im / s_0l.abs2());
                    let term = s[i][l].times(s[j][l]).times(s[k][l].conj());
                    sum = sum.plus(term.times(inv));
                }
                if (sum.re - n_val).abs() > tol || sum.im.abs() > tol {
                    return Err(format!("Verlinde formula fails at N_{{{i},{j}}}^{k}"));
                }
            }
        }
    }

    // R-symbols are unit phases on admissible channels, and the balancing (ribbon)
    // equation holds exactly: R^k_{ij} R^k_{ji} = θ_k θ_i⁻¹ θ_j⁻¹.
    for i in 0..dim {
        for j in 0..dim {
            for k in 0..dim {
                if m.n_ijk(i, j, k).abs() > tol {
                    let r1 = m.r_symbol(i, j, k);
                    if (r1.abs2() - 1.0).abs() > tol {
                        return Err(format!("R_{{{i},{j}}}^{k} is not a unit phase"));
                    }
                    if m.n_ijk(j, i, k).abs() > tol {
                        let r2 = m.r_symbol(j, i, k);
                        let lhs = r1.times(r2);
                        let rhs = t[k].times(t[i].conj()).times(t[j].conj());
                        if !lhs.close(rhs, tol) {
                            return Err(format!("Balancing equation fails at N_{{{i},{j}}}^{k}"));
                        }
                    }
                }
            }
        }
    }

    // Monodromy–S relation: D·S_{ab} = Σ_c N_{āb}^c d_c θ_c θ_ā⁻¹ θ_b⁻¹, where ā is the
    // charge conjugate of a. This ties the braiding/twist data back to the S matrix.
    let cc = m.charge_conjugation();
    for a in 0..dim {
        let a_bar = (0..dim)
            .find(|&x| (cc[a][x].re - 1.0).abs() <= tol && cc[a][x].im.abs() <= tol)
            .ok_or_else(|| format!("charge conjugation row {a} is not a permutation row"))?;
        for b in 0..dim {
            let mut sum = C::new(0.0, 0.0);
            for c in 0..dim {
                let n_val = m.n_ijk(a_bar, b, c);
                if n_val.abs() > tol {
                    let term = t[c].times(t[a_bar].conj()).times(t[b].conj());
                    sum = sum.plus(term.scale(n_val * dims[c]));
                }
            }
            let want = s[a][b].scale(d_total);
            if !sum.close(want, tol) {
                return Err(format!("monodromy–S relation fails at ({a},{b})"));
            }
        }
    }

    // Pentagon equation, exact.
    for i1 in 0..dim {
        for i2 in 0..dim {
            for a in 0..dim {
                if m.n_ijk(i1, i2, a).abs() < tol {
                    continue;
                }
                for i3 in 0..dim {
                    for b in 0..dim {
                        if m.n_ijk(a, i3, b).abs() < tol {
                            continue;
                        }
                        for i4 in 0..dim {
                            for c in 0..dim {
                                if m.n_ijk(i3, i4, c).abs() < tol {
                                    continue;
                                }
                                for i5 in 0..dim {
                                    if m.n_ijk(b, i4, i5).abs() < tol {
                                        continue;
                                    }
                                    if m.n_ijk(a, c, i5).abs() < tol {
                                        continue;
                                    }
                                    for d in 0..dim {
                                        if m.n_ijk(i2, c, d).abs() < tol {
                                            continue;
                                        }
                                        if m.n_ijk(i1, d, i5).abs() < tol {
                                            continue;
                                        }

                                        let lhs = m
                                            .f_symbol(i1, i2, c, i5, a, d)
                                            .times(m.f_symbol(a, i3, i4, i5, b, c));

                                        let mut rhs = C::new(0.0, 0.0);
                                        for e_idx in 0..dim {
                                            if m.n_ijk(i2, i3, e_idx).abs() < tol {
                                                continue;
                                            }
                                            if m.n_ijk(i1, e_idx, b).abs() < tol {
                                                continue;
                                            }
                                            if m.n_ijk(e_idx, i4, d).abs() < tol {
                                                continue;
                                            }

                                            let term = m
                                                .f_symbol(i2, i3, i4, d, e_idx, c)
                                                .times(m.f_symbol(i1, e_idx, i4, i5, b, d))
                                                .times(m.f_symbol(i1, i2, i3, b, a, e_idx));
                                            rhs = rhs.plus(term);
                                        }

                                        if !lhs.close(rhs, tol) {
                                            return Err(format!("Pentagon equation fails at 1={i1} 2={i2} 3={i3} 4={i4} 5={i5} a={a} b={b} c={c} d={d}"));
                                        }
                                    }
                                }
                            }
                        }
                    }
                }
            }
        }
    }

    // Hexagon equations, exact.
    for i1 in 0..dim {
        for i2 in 0..dim {
            for a in 0..dim {
                if m.n_ijk(i1, i2, a).abs() < tol {
                    continue;
                }
                if m.n_ijk(i2, i1, a).abs() < tol {
                    continue;
                }
                for i3 in 0..dim {
                    for d in 0..dim {
                        if m.n_ijk(a, i3, d).abs() < tol {
                            continue;
                        }

                        for c in 0..dim {
                            if m.n_ijk(i1, i3, c).abs() < tol {
                                continue;
                            }
                            if m.n_ijk(i2, c, d).abs() < tol {
                                continue;
                            }

                            let lhs1 = m
                                .r_symbol(i1, i3, c)
                                .times(m.f_symbol(i2, i1, i3, d, a, c))
                                .times(m.r_symbol(i1, i2, a));

                            let lhs2 = m
                                .r_symbol(i3, i1, c)
                                .conj()
                                .times(m.f_symbol(i2, i1, i3, d, a, c))
                                .times(m.r_symbol(i2, i1, a).conj());

                            let mut rhs1 = C::new(0.0, 0.0);
                            let mut rhs2 = C::new(0.0, 0.0);

                            for b in 0..dim {
                                if m.n_ijk(i2, i3, b).abs() < tol {
                                    continue;
                                }
                                if m.n_ijk(i1, b, d).abs() < tol {
                                    continue;
                                }
                                if m.n_ijk(b, i1, d).abs() < tol {
                                    continue;
                                }

                                let f1 = m.f_symbol(i1, i2, i3, d, a, b);
                                let f2 = m.f_symbol(i2, i3, i1, d, b, c);

                                rhs1 = rhs1.plus(f2.times(m.r_symbol(i1, b, d)).times(f1));
                                rhs2 = rhs2.plus(f2.times(m.r_symbol(b, i1, d).conj()).times(f1));
                            }

                            if !lhs1.close(rhs1, tol) {
                                return Err(format!(
                                    "Hexagon 1 fails at 1={i1} 2={i2} 3={i3} d={d} a={a} c={c}"
                                ));
                            }
                            if !lhs2.close(rhs2, tol) {
                                return Err(format!(
                                    "Hexagon 2 fails at 1={i1} 2={i2} 3={i3} d={d} a={a} c={c}"
                                ));
                            }
                        }
                    }
                }
            }
        }
    }

    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::native::construct_atlas_native;
    use crate::DoubleZn;
    use tqc_core::params::UseCaseParams;

    #[test]
    fn verify_generalized_mtc_axioms_for_double_zn() {
        let tol = 1e-9;
        // Verify for a small double category
        let d2 = DoubleZn { n: 2 };
        verify_mtc_axioms(&d2, tol)
            .unwrap_or_else(|e| panic!("verify_mtc_axioms failed for D(Z_2): {e}"));

        let d3 = DoubleZn { n: 3 };
        verify_mtc_axioms(&d3, tol)
            .unwrap_or_else(|e| panic!("verify_mtc_axioms failed for D(Z_3): {e}"));
    }

    #[test]
    fn verify_full_phase_exact_axioms_for_atlas_native() {
        // The Atlas-native pointed category C(Z_3 × Z_2^3, q) must pass every axiom
        // with full complex-valued (phase-sensitive) comparisons — and so must the
        // second (even-modality) use-case instance, which exercises the level-1
        // Eilenberg–MacLane cocycle branch.
        // (2,4,4) exercises an even modality >= 4, where the S-matrix bicharacter level
        // differs from the naive +2·m1m2/n exponent.
        for (scope, modality, context) in [(4u32, 3u32, 8u32), (2, 2, 4), (2, 4, 4), (3, 5, 2)] {
            let p = UseCaseParams::new(scope, modality, context);
            let native = construct_atlas_native(&p).expect("construction");
            verify_mtc_axioms(&*native, 1e-9).unwrap_or_else(|e| {
                panic!("phase-exact axioms failed at ({scope},{modality},{context}): {e}")
            });
        }
    }

    #[test]
    fn phase_blind_data_is_rejected() {
        // A mutant of D(Z_2) whose R-symbol carries a wrong sign must be caught. The
        // magnitude of every symbol is unchanged, so a magnitude-only checker would
        // accept it; the phase-exact checker must reject it.
        struct WrongSign(DoubleZn);
        impl ModularData for WrongSign {
            fn dim(&self) -> usize {
                ModularData::dim(&self.0)
            }
            fn s_matrix(&self) -> Matrix {
                ModularData::s_matrix(&self.0)
            }
            fn t_diag(&self) -> Vec<C> {
                ModularData::t_diag(&self.0)
            }
            fn charge_conjugation(&self) -> Matrix {
                ModularData::charge_conjugation(&self.0)
            }
            fn n_ijk(&self, i: usize, j: usize, k: usize) -> f64 {
                self.0.n_ijk(i, j, k)
            }
            fn f_symbol(&self, i: usize, j: usize, k: usize, l: usize, m: usize, n: usize) -> C {
                self.0.f_symbol(i, j, k, l, m, n)
            }
            fn r_symbol(&self, i: usize, j: usize, k: usize) -> C {
                // Flip the sign on one non-identity channel.
                let r = self.0.r_symbol(i, j, k);
                if i == 1 && j == 1 {
                    r.scale(-1.0)
                } else {
                    r
                }
            }
        }
        let mutant = WrongSign(DoubleZn { n: 2 });
        assert!(
            verify_mtc_axioms(&mutant, 1e-9).is_err(),
            "a phase error must be detected"
        );
    }
}
