//! The complex amplitude encoding: a fusion-space vector as a content-addressable map.
//!
//! Realizes the pure side of the `complex-amplitude-encoding` dictionary row. A state is a map
//! `{ label index → complex amplitude }`; this module gives its canonical byte encoding (which
//! the substrate then content-addresses) and the Euclidean composition norm `Σ|cᵢ|²`. Exact
//! Gaussian-integer amplitudes keep it no_std-clean and reproducible.

use alloc::format;
use alloc::string::String;
use alloc::vec::Vec;

/// A complex amplitude as a Gaussian integer.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub struct Amplitude {
    /// Real part.
    pub re: i64,
    /// Imaginary part.
    pub im: i64,
}

/// Canonical byte encoding of a state (sorted by label; `label,re,im` joined by `;`).
#[must_use]
pub fn encode(state: &[(u64, Amplitude)]) -> Vec<u8> {
    let mut v: Vec<(u64, Amplitude)> = state.to_vec();
    v.sort_by_key(|(l, _)| *l);
    let parts: Vec<String> = v
        .iter()
        .map(|(l, a)| format!("{l},{},{}", a.re, a.im))
        .collect();
    parts.join(";").into_bytes()
}

/// Decode a canonical encoding back into a state (in canonical, label-sorted order).
#[must_use]
pub fn decode(bytes: &[u8]) -> Option<Vec<(u64, Amplitude)>> {
    let s = core::str::from_utf8(bytes).ok()?;
    if s.is_empty() {
        return Some(Vec::new());
    }
    let mut out = Vec::new();
    for part in s.split(';') {
        let mut it = part.split(',');
        let l = it.next()?.parse().ok()?;
        let re = it.next()?.parse().ok()?;
        let im = it.next()?.parse().ok()?;
        if it.next().is_some() {
            return None;
        }
        out.push((l, Amplitude { re, im }));
    }
    Some(out)
}

/// Fixed-width interleaved binary encoding of a dense state over `n` labels:
/// `[re₀, im₀, re₁, im₁, …]` as little-endian `i64`s (16 bytes per label). This is the wire
/// format the holospace gate path executes on; [`decode_interleaved`] is its exact inverse.
#[must_use]
pub fn encode_interleaved(n: u64, state: &[(u64, Amplitude)]) -> Vec<u8> {
    let mut v = alloc::vec![0i64; (n * 2) as usize];
    for &(l, a) in state {
        let l = l as usize;
        v[l * 2] = a.re;
        v[l * 2 + 1] = a.im;
    }
    v.iter().flat_map(|x| x.to_le_bytes()).collect()
}

/// Decode the fixed-width interleaved binary encoding. Trailing bytes that do not fill a
/// full 16-byte label record are rejected (never silently zeroed).
#[must_use]
pub fn decode_interleaved(bytes: &[u8]) -> Option<Vec<(u64, Amplitude)>> {
    if bytes.len() % 16 != 0 {
        return None;
    }
    let mut out = Vec::with_capacity(bytes.len() / 16);
    for (i, chunk) in bytes.chunks_exact(16).enumerate() {
        let mut re_b = [0u8; 8];
        let mut im_b = [0u8; 8];
        re_b.copy_from_slice(&chunk[0..8]);
        im_b.copy_from_slice(&chunk[8..16]);
        out.push((
            i as u64,
            Amplitude {
                re: i64::from_le_bytes(re_b),
                im: i64::from_le_bytes(im_b),
            },
        ));
    }
    Some(out)
}

/// The Euclidean composition norm `Σ|cᵢ|² = Σ(reᵢ² + imᵢ²)`.
#[must_use]
pub fn norm_sq(state: &[(u64, Amplitude)]) -> u128 {
    state.iter().fold(0u128, |acc, (_, a)| {
        let re_sq = (a.re.unsigned_abs() as u128).saturating_pow(2);
        let im_sq = (a.im.unsigned_abs() as u128).saturating_pow(2);
        acc.saturating_add(re_sq.saturating_add(im_sq))
    })
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn encode_decode_round_trips() {
        let state = [
            (2, Amplitude { re: 1, im: -1 }),
            (0, Amplitude { re: 3, im: 0 }),
        ];
        let decoded = decode(&encode(&state)).unwrap();
        // Canonical order is label-sorted.
        assert_eq!(
            decoded,
            [
                (0, Amplitude { re: 3, im: 0 }),
                (2, Amplitude { re: 1, im: -1 })
            ]
        );
    }

    #[test]
    fn norm_is_sum_of_squared_moduli() {
        let state = [
            (0, Amplitude { re: 3, im: 4 }),
            (1, Amplitude { re: 0, im: 2 }),
        ];
        assert_eq!(norm_sq(&state), 25 + 4);
    }
}
