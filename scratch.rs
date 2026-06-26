use std::collections::HashSet;

fn cd_mul(a: &[i128], b: &[i128]) -> Vec<i128> {
    let n = a.len();
    if n == 1 { return vec![a[0] * b[0]]; }
    let h = n / 2;
    let (a1, a2) = a.split_at(h);
    let (b1, b2) = b.split_at(h);
    
    let conj = |x: &[i128]| -> Vec<i128> {
        let mut out = x.to_vec();
        for i in 1..x.len() { out[i] = -out[i]; }
        out
    };
    
    let first = cd_mul(a1, b1).iter().zip(cd_mul(&conj(b2), a2)).map(|(x,y)| x - y).collect::<Vec<_>>();
    let second = cd_mul(b2, a1).iter().zip(cd_mul(a2, &conj(b1))).map(|(x,y)| x + y).collect::<Vec<_>>();
    
    let mut out = first;
    out.extend(second);
    out
}

fn main() {
    let dim = 8;
    for i in 0..dim {
        for j in 0..dim {
            for k in 0..dim {
                // Check (i * j) * k == i * (j * k) with absolute values
                let mut ei = vec![0; dim]; ei[i] = 1;
                let mut ej = vec![0; dim]; ej[j] = 1;
                let mut ek = vec![0; dim]; ek[k] = 1;
                
                let ij = cd_mul(&ei, &ej);
                let mut ij_abs = vec![0; dim];
                for x in 0..dim { ij_abs[x] = ij[x].abs(); }
                
                let ij_k = cd_mul(&ij_abs, &ek);
                let mut ij_k_abs = vec![0; dim];
                for x in 0..dim { ij_k_abs[x] = ij_k[x].abs(); }
                
                let jk = cd_mul(&ej, &ek);
                let mut jk_abs = vec![0; dim];
                for x in 0..dim { jk_abs[x] = jk[x].abs(); }
                
                let i_jk = cd_mul(&ei, &jk_abs);
                let mut i_jk_abs = vec![0; dim];
                for x in 0..dim { i_jk_abs[x] = i_jk[x].abs(); }
                
                if ij_k_abs != i_jk_abs {
                    println!("Not associative! {} {} {}", i, j, k);
                    return;
                }
            }
        }
    }
    println!("Associative!");
}
