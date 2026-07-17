//! Pins the holospace gate execution semantics: `execute_holo_gate` is a **Gather**
//! (`out[i] = state[targets[i]]`), so realizing the amplitude action `out[g(i)] = v[i]`
//! requires `targets[i] = g⁻¹(i)`. This test also rejects the passthrough failure mode
//! (an archive that echoes its input), which the determinism-only checks cannot see.
#![allow(clippy::unwrap_used)]

use tqc_core::amplitude;
use tqc_core::generators::Generators;
use tqc_core::UseCaseParams;

#[test]
fn holo_gate_is_a_real_gather_not_a_passthrough() {
    let p = UseCaseParams::new(4, 3, 8);
    let g = Generators::new(&p);
    let n = p.class_count() as usize;
    let base: Vec<i64> = (0..n as i64).map(|i| i % 5 - 2).collect();
    let amp0: Vec<(u64, amplitude::Amplitude)> = base
        .iter()
        .enumerate()
        .map(|(i, &v)| (i as u64, amplitude::Amplitude { re: v, im: 0 }))
        .collect();
    let bytes = amplitude::encode_interleaved(p.class_count(), &amp0);

    // Forward targets: the gather convention applies g to positions, i.e. out[i] = v[g(i)].
    let fwd: Vec<usize> = (0..p.class_count())
        .map(|i| g.sigma.apply(i) as usize)
        .collect();
    let exec = tqc_substrate::execute_holo_gate("sigma", &fwd, &bytes).unwrap();
    let out: Vec<i64> = amplitude::decode_interleaved(&exec.output_bytes)
        .unwrap()
        .iter()
        .map(|&(_, a)| a.re)
        .collect();
    let gather_expected: Vec<i64> = (0..n)
        .map(|i| base[g.sigma.apply(i as u64) as usize])
        .collect();
    assert_eq!(
        out, gather_expected,
        "holo gate must implement the Gather semantics"
    );
    assert_ne!(out, base, "holo gate must not be an input passthrough");

    // Inverse targets realize the amplitude action out[g(i)] = v[i] — the same action
    // `permute_amplitudes` implements. This is the binding the witnesses rely on.
    let inv = g.sigma.inverse();
    let inv_targets: Vec<usize> = (0..p.class_count())
        .map(|i| inv.apply(i) as usize)
        .collect();
    let exec = tqc_substrate::execute_holo_gate("sigma", &inv_targets, &bytes).unwrap();
    let out: Vec<i64> = amplitude::decode_interleaved(&exec.output_bytes)
        .unwrap()
        .iter()
        .map(|&(_, a)| a.re)
        .collect();
    assert_eq!(
        out,
        g.sigma.permute_amplitudes(&base),
        "inverse targets must realize the permutation amplitude action"
    );
}
