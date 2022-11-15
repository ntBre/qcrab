#![allow(unused)]

use super::{contraction::Contraction, shell::Shell};
use crate::{
    hf_libint::dfact,
    molecule::{Atom, Molecule},
    Dmat,
};
use itertools::Itertools;
use std::{
    cmp::{max, min},
    f64::consts::PI,
    fs::read_to_string,
    ops::Index,
    path::Path,
};

#[derive(PartialEq, Debug)]
pub(crate) struct Basis(pub(crate) Vec<Shell>);

mod json {
    use std::collections::HashMap;

    use serde::Deserialize;

    #[derive(Debug, Deserialize)]
    pub(crate) struct JsonBasis {
        pub(crate) elements: HashMap<usize, Shells>,
    }

    #[derive(Debug, Deserialize)]
    pub(crate) struct Shells {
        pub(crate) electron_shells: Vec<RawShell>,
    }

    #[derive(Debug, Deserialize)]
    pub(crate) struct RawShell {
        pub(crate) angular_momentum: Vec<usize>,
        pub(crate) exponents: Vec<String>,
        pub(crate) coefficients: Vec<Vec<String>>,
        pub(crate) function_type: String,
    }
}

#[test]
fn load() {
    use nalgebra::vector;

    let mol = Molecule::load("testfiles/h2o/STO-3G/geom.dat");
    let got = Basis::load("basis_sets/sto-3g.json", &mol);
    let want = Basis(vec![
        Shell {
            alpha: vec![130.7093214, 23.80886605, 6.443608313],
            contr: vec![Contraction {
                l: 0,
                pure: false,
                coeff: vec![
                    4.251943277787213,
                    4.11229442420408,
                    1.2816225514343584,
                ],
            }],
            origin: vector![0.0, -0.143225816552, 0.0],
        },
        Shell {
            alpha: vec![5.033151319, 1.169596125, 0.38038896],
            contr: vec![Contraction {
                l: 0,
                pure: false,
                coeff: vec![
                    -0.2394130049456894,
                    0.32023423543952656,
                    0.24168555455632082,
                ],
            }],
            origin: vector![0.0, -0.143225816552, 0.0],
        },
        Shell {
            alpha: vec![5.033151319, 1.169596125, 0.38038896],
            contr: vec![Contraction {
                l: 1,
                pure: false,
                coeff: vec![
                    1.6754501961195145,
                    1.0535680440115096,
                    0.1669028790880833,
                ],
            }],
            origin: vector![0.0, -0.143225816552, 0.0],
        },
        Shell {
            alpha: vec![3.425250914, 0.6239137298, 0.168855404],
            contr: vec![Contraction {
                l: 0,
                pure: false,
                coeff: vec![
                    0.2769343550790519,
                    0.26783885160947885,
                    0.08347367112984118,
                ],
            }],
            origin: vector![1.638036840407, 1.136548822547, -0.0],
        },
        Shell {
            alpha: vec![3.425250914, 0.6239137298, 0.168855404],
            contr: vec![Contraction {
                l: 0,
                pure: false,
                coeff: vec![
                    0.2769343550790519,
                    0.26783885160947885,
                    0.08347367112984118,
                ],
            }],
            origin: vector![-1.638036840407, 1.136548822547, -0.0],
        },
    ]);
    assert!(got.len() == want.len());
    assert_eq!(got, want, "got\n{:#?}, want\n{:#?}", got, want);
}

impl Basis {
    /// load the basis set from `path` and combine it with the origins in `mol`.
    /// panics if `path` cannot be read, if its data cannot be deserialized, and
    /// if it does not contain one of the elements in `mol`.
    #[allow(unused)]
    pub(crate) fn load(path: impl AsRef<Path>, mol: &Molecule) -> Self {
        let data = read_to_string(&path).unwrap();
        let raw: json::JsonBasis = serde_json::from_str(&data).unwrap();
        let mut shells = Vec::new();
        for Atom {
            atomic_number,
            coord,
        } in &mol.atoms
        {
            for shell in &raw.elements[atomic_number].electron_shells {
                let alpha: Vec<f64> = shell
                    .exponents
                    .iter()
                    .map(|s| s.parse().unwrap())
                    .collect();
                for (i, l) in shell.angular_momentum.iter().enumerate() {
                    let coeff: Vec<f64> = shell.coefficients[i]
                        .iter()
                        .map(|s| s.parse().unwrap())
                        .collect();
                    shells.push(Shell::new(
                        alpha.clone(),
                        vec![Contraction::new(
                            *l,
                            shell.function_type == "gto_spherical",
                            coeff,
                        )],
                        *coord,
                    ));
                }
            }
        }
        Self(shells)
    }

    pub(crate) fn nbasis(&self) -> usize {
        self.0.iter().map(|s| s.size()).sum()
    }

    /// the name pretty much says it don't you think?
    pub(crate) fn map_shell_to_basis_function(&self) -> Vec<usize> {
        let mut ret = Vec::with_capacity(self.len());
        let mut n = 0;
        for shell in &self.0 {
            ret.push(n);
            n += shell.size();
        }
        ret
    }

    /// STO-3G basis set
    ///
    /// cite: W. J. Hehre, R. F. Stewart, and J. A. Pople, The Journal of Chemical
    /// Physics 51, 2657 (1969) doi: 10.1063/1.1672392
    ///
    /// obtained from https://bse.pnl.gov/bse/portal via libint
    pub(crate) fn sto3g(mol: &Molecule) -> Self {
        let mut shells = Vec::new();
        for Atom {
            atomic_number,
            coord,
        } in &mol.atoms
        {
            match atomic_number {
                // hydrogen
                1 => shells.push(Shell::new(
                    vec![3.425250910, 0.623913730, 0.168855400],
                    vec![Contraction::new(
                        0,
                        false,
                        vec![0.15432897, 0.53532814, 0.44463454],
                    )],
                    *coord,
                )),
                // carbon
                6 => shells.extend([
                    Shell::new(
                        vec![71.616837000, 13.045096000, 3.530512200],
                        vec![Contraction::new(
                            0,
                            false,
                            vec![0.15432897, 0.53532814, 0.44463454],
                        )],
                        *coord,
                    ),
                    Shell::new(
                        vec![2.941249400, 0.683483100, 0.222289900],
                        vec![Contraction::new(
                            0,
                            false,
                            vec![-0.09996723, 0.39951283, 0.70011547],
                        )],
                        *coord,
                    ),
                    Shell::new(
                        vec![2.941249400, 0.683483100, 0.222289900],
                        vec![Contraction::new(
                            1,
                            false,
                            vec![0.15591627, 0.60768372, 0.39195739],
                        )],
                        *coord,
                    ),
                ]),
                // oxygen
                8 => shells.extend([
                    Shell::new(
                        vec![130.709320000, 23.808861000, 6.443608300],
                        vec![Contraction::new(
                            0,
                            false,
                            vec![0.15432897, 0.53532814, 0.44463454],
                        )],
                        *coord,
                    ),
                    Shell::new(
                        vec![5.033151300, 1.169596100, 0.380389000],
                        vec![Contraction::new(
                            0,
                            false,
                            vec![-0.09996723, 0.39951283, 0.70011547],
                        )],
                        *coord,
                    ),
                    Shell::new(
                        vec![5.033151300, 1.169596100, 0.380389000],
                        vec![Contraction::new(
                            1,
                            false,
                            vec![0.15591627, 0.60768372, 0.39195739],
                        )],
                        *coord,
                    ),
                ]),
                _ => panic!("unsupported element {atomic_number}"),
            }
            // shells.push(
        }
        Self(shells)
    }

    pub(crate) fn len(&self) -> usize {
        self.0.len()
    }

    pub(crate) fn overlap_ints(&self) -> Dmat {
        let n = self.nbasis();
        let mut result = Dmat::zeros(n, n);

        let shell2bf = self.map_shell_to_basis_function();

        for s1 in 0..self.len() {
            let bf1 = shell2bf[s1];
            let n1 = self[s1].cartesian_size();
            for s2 in 0..=s1 {
                let bf2 = shell2bf[s2];
                let n2 = self[s2].cartesian_size();
                let s1: &Shell = &self[s1];
                let s2: &Shell = &self[s2];
                assert!(
                    s1.ncontr() == 1 && s2.ncontr() == 1,
                    "generally-contracted shells not yet supported"
                );

                let l1 = s1.contr[0].l;
                let l2 = s2.contr[0].l;
                let mut buf = vec![0.0; n1 * n2];
                // somewhere in here I have to loop over the actual
                // orbitals, not just the shells. that's where the missing
                // entries are coming from - 3*1 for pxs overlap and 3*3 for
                // the pxp overlap
                for (p1, alpha) in s1.alpha.iter().enumerate() {
                    for (p2, beta) in s2.alpha.iter().enumerate() {
                        let a = s1.origin;
                        let b = s2.origin;
                        let c1 = s1.contr[0].coeff[p1];
                        let c2 = s2.contr[0].coeff[p2];
                        let gammap = alpha + beta;
                        let oogammap = 1.0 / gammap;
                        let rhop_over_alpha1 = beta * oogammap;
                        let rhop = alpha * rhop_over_alpha1;
                        let ab_x = a[0] - b[0];
                        let ab_y = a[1] - b[1];
                        let ab_z = a[2] - b[2];
                        let ab2_x = ab_x * ab_x;
                        let ab2_y = ab_y * ab_y;
                        let ab2_z = ab_z * ab_z;

                        let ovlp_ss_x = f64::exp(-rhop * ab2_x) * c1 * c2;
                        let ovlp_ss_y = f64::exp(-rhop * ab2_y);
                        let ovlp_ss_z = f64::exp(-rhop * ab2_z);

                        let p = (*alpha * a + *beta * b) / gammap;
                        let pa = p - a;
                        let pb = p - b;

                        let ls1 = l_perm(l1);
                        let ls2 = l_perm(l2);

                        let pfac = (PI / (alpha + beta)).sqrt();
                        let combos = ls1.iter().cartesian_product(ls2);
                        for (i, ((l1x, l1y, l1z), (l2x, l2y, l2z))) in
                            combos.enumerate()
                        {
                            let ix = pfac
                                * s_xyz(*l1x, l2x, pa.x, pb.x, alpha, beta);
                            let iy = pfac
                                * s_xyz(*l1y, l2y, pa.y, pb.y, alpha, beta);
                            let iz = pfac
                                * s_xyz(*l1z, l2z, pa.z, pb.z, alpha, beta);
                            buf[i] += ovlp_ss_x
                                * ovlp_ss_y
                                * ovlp_ss_z
                                * ix
                                * iy
                                * iz;
                        }
                    }
                }
                let mut ij = 0;
                for i in bf1..bf1 + n1 {
                    for j in bf2..bf2 + n2 {
                        // println!("{} {} {:.8}", i + 1, j + 1, buf[ij]);
                        result[(i, j)] = buf[ij];
                        result[(j, i)] = buf[ij];
                        ij += 1;
                    }
                }
            }
        }
        result
    }
}

fn tform_cols(n1: usize, l: usize, src: Vec<f64>) -> Vec<f64> {
    let coefs = SolidHarmonicCoeffs::new(l);

    let nc = (l + 1) * (l + 2) / 2;
    let n = 2 * l + 1;
    let mut tgt = vec![0.0; n1 * n];

    // loop over shg
    for s in 0..n {
        let nc_s = coefs.nnz(s); // # of cartesians contributing to shg s
        let c_idxs = coefs.row_idx(s); // indices of cartesians contributing to shg s
        let c_vals = coefs.row_values(s); // coefficients of cartesians contributing to shg s

        let tgt_blk_s_offset = s;

        for ic in 0..nc_s {
            let c = c_idxs[ic];
            let s_c_coeff = c_vals[ic];

            let src_blk_s = &src[c..];
            let tgt_blk_s = &mut tgt[tgt_blk_s_offset..];

            // loop over other dims
            let mut si = 0;
            let mut ti = 0;
            (0..n1).for_each(|_| {
                tgt_blk_s[ti] += s_c_coeff * src_blk_s[si];
                si += nc;
                ti += n;
            });
        }
    }

    tgt
}

fn tform_rows(l: usize, n2: usize, src: Vec<f64>) -> Vec<f64> {
    let coefs = SolidHarmonicCoeffs::new(l);
    let n = 2 * l + 1;
    let mut tgt = vec![0.0; n * n2];
    for s in 0..n {
        let nc_s = coefs.nnz(s);
        let c_idxs = coefs.row_idx(s);
        let c_vals = coefs.row_values(s);
        let tgt_blk_s_offset = s * n2;
        for ic in 0..nc_s {
            let c = c_idxs[ic];
            let s_c_coeff = c_vals[ic];
            let src_blk_s = &src[c * n2..];
            let tgt_blk_s = &mut tgt[tgt_blk_s_offset..];

            // loop over other dims
            (0..n2).for_each(|i2| {
                tgt_blk_s[i2] += s_c_coeff * src_blk_s[ic];
            });
        }
    }
    tgt
}

struct SolidHarmonicCoeffs {
    /// elements
    values: Vec<f64>,

    /// column indices
    colidx: Vec<usize>,

    /// "pointer" to the beginning of each row
    row_offset: Vec<usize>,
}

impl SolidHarmonicCoeffs {
    fn new(l: usize) -> Self {
        let npure = 2 * l + 1;
        let ncart = (l + 1) * (l + 2) / 2;
        let mut full_coeff = vec![0.0; npure * ncart];
        let mut m = -(l as isize);
        // assuming LIBINT_SHGSHELL_ORDERING_STANDARD, not GAUSSIAN ordering
        for pure_idx in 0..npure {
            let mut cart_idx = 0;
            // this is the FOR_CART macro
            for lx in (0..=l).rev() {
                for ly in (0..l - lx).rev() {
                    let lz = l - lx - ly;
                    full_coeff[pure_idx * ncart + cart_idx] =
                        coeff(l, m, lx, ly, lz);
                    cart_idx += 1;
                }
            }
            m += 1;
        }

        let nnz = full_coeff.iter().map(|c| usize::from(*c != 0.0)).sum();
        let mut values = vec![0.0; nnz];
        let mut colidx = vec![0; nnz];
        let mut row_offset = vec![0; npure + 1];

        let mut pc = 0;
        let mut cnt = 0;
        (0..npure).for_each(|p| {
            row_offset[p] = cnt;
            for c in 0..ncart {
                if full_coeff[pc] != 0.0 {
                    values[cnt] = full_coeff[pc];
                    colidx[cnt] = c;
                    cnt += 1;
                }
                pc += 1;
            }
        });
        row_offset[npure] = cnt;

        Self {
            values,
            colidx,
            row_offset,
        }
    }

    fn nnz(&self, r: usize) -> usize {
        self.row_offset[r + 1] - self.row_offset[r]
    }

    fn row_idx(&self, r: usize) -> &[usize] {
        &self.colidx[self.row_offset[r]..]
    }

    fn row_values(&self, r: usize) -> &[f64] {
        &self.values[self.row_offset[r]..]
    }
}

fn parity(i: isize) -> isize {
    if i % 2 != 0 {
        -1
    } else {
        1
    }
}

fn coeff(l: usize, m: isize, lx: usize, ly: usize, lz: usize) -> f64 {
    let l = l as isize;
    let lx = lx as isize;
    let ly = ly as isize;
    let lz = lz as isize;
    let abs_m = m.abs();
    if (lx + ly - abs_m) % 2 != 0 {
        return 0.0;
    }

    let j = (lx + ly - abs_m) / 2;
    if j < 0 {
        return 0.0;
    }

    /*----------------------------------------------------------------------------------------
     Checking whether the cartesian polynomial contributes to the requested component of Ylm
    ----------------------------------------------------------------------------------------*/
    let comp = if m >= 0 { 1 } else { -1 };
    /*  if (comp != ((abs_m-lx)%2 ? -1 : 1))*/
    let i = abs_m - lx;
    if comp != parity(i.abs()) {
        return 0.0;
    }

    let mut pfac = f64::sqrt(
        (((fac(2 * lx)) * (fac(2 * ly)) * (fac(2 * lz))) / fac(2 * l))
            * ((fac(l - abs_m)) / (fac(l)))
            * ((1.0) / fac(l + abs_m))
            * ((1.0) / (fac(lx) * fac(ly) * fac(lz))),
    );
    /*  pfac = sqrt(fac[l-abs_m]/(fac[l]*fac[l]*fac[l+abs_m]));*/
    pfac /= (1 << l) as f64;
    if m < 0 {
        pfac *= parity((i - 1) / 2) as f64;
    } else {
        pfac *= parity(i / 2) as f64;
    }

    let i_min = j;
    let i_max = (l - abs_m) / 2;
    let mut sum = 0.0;
    for i in i_min..=i_max {
        let mut pfac1 = binom(l, i) * binom(i, j);
        pfac1 *= (parity(i) * factorial(2 * (l - i))) as f64
            / fac(l - abs_m - 2 * i);
        let mut sum1 = 0.0;
        let k_min = max((lx - abs_m) / 2, 0);
        let k_max = min(j, lx / 2);
        for k in k_min..=k_max {
            if lx - 2 * k <= abs_m {
                sum1 +=
                    binom(j, k) * binom(abs_m, lx - 2 * k) * parity(k) as f64;
            }
        }
        sum += pfac1 * sum1;
    }
    sum *= f64::sqrt(
        (dfact(2 * l)) / (dfact(2 * lx) * dfact(2 * ly) * dfact(2 * lz)),
    );

    if m == 0 {
        pfac * sum
    } else {
        std::f64::consts::SQRT_2 * pfac * sum
    }
}

fn tform(l_row: usize, l_col: usize, source_blk: Vec<f64>) -> Vec<f64> {
    // in libint these are built by calling `instance` but as far as I can tell
    // `instance` just returns the lth element of an iterator of
    // SolidHarmoniccoeffs, which should just be the same as calling my `new`
    // directly with `l`
    let coefs_row = SolidHarmonicCoeffs::new(l_row);
    let coefs_col = SolidHarmonicCoeffs::new(l_col);

    let ncart_col = (l_col + 1) * (l_col + 2) / 2;
    let npure_row = 2 * l_row + 1;
    let npure_col = 2 * l_col + 1;
    let mut target_blk = vec![0.0; npure_row * npure_col];

    // loop over row shg
    for s1 in 0..npure_row {
        let nc1 = coefs_row.nnz(s1); // # of cartesians contributing to shg s1
        let c1_idxs = coefs_row.row_idx(s1); // indices of cartesians contributing to shg s1
        let c1_vals = coefs_row.row_values(s1); // coefficients of cartesians contributing to shg s1

        let target_blk_s1 = &mut target_blk[s1 * npure_col..];

        // loop over col shg
        (0..npure_col).for_each(|s2| {
            let nc2 = coefs_col.nnz(s2); // # of cartesians contributing to shg s2
            let c2_idxs = coefs_col.row_idx(s2); // indices of cartesians contributing to shg s2
            let c2_vals = coefs_col.row_values(s2); // coefficients of cartesians contributing to shg s2

            for ic1 in 0..nc1 {
                let c1 = c1_idxs[ic1];
                let s1_c1_coeff = c1_vals[ic1];

                let source_blk_c1 = &source_blk[c1 * ncart_col..];

                for ic2 in 0..nc2 {
                    let c2 = c2_idxs[ic2];
                    let s2_c2_coeff = c2_vals[ic2];

                    target_blk_s1[s2] +=
                        source_blk_c1[c2] * s1_c1_coeff * s2_c2_coeff;
                }
            }
        });
    }

    target_blk
}

/// return the permutations of angular momentum values associated with `l`.
/// panics if `l` greater than 2
fn l_perm(l: usize) -> Vec<(isize, isize, isize)> {
    match l {
        0 => vec![(0, 0, 0)],
        1 => vec![(1, 0, 0), (0, 1, 0), (0, 0, 1)],
        2 => vec![
            (2, 0, 0),
            (0, 1, 1),
            (1, 1, 0),
            (1, 0, 1),
            (0, 0, 2),
            (0, 2, 0),
        ],
        _ => panic!("unmatched l value {l}"),
    }
}

/// compute individual S_[xyz] values
fn s_xyz(
    ax: isize,
    bx: isize,
    pa: f64,
    pb: f64,
    alpha: &f64,
    beta: &f64,
) -> f64 {
    let mut acc = 0.0;
    for ix in 0..=ax {
        for jx in 0..=bx {
            if (ix + jx) % 2 != 1 {
                acc += binom(ax, ix)
                    * binom(bx, jx)
                    * dfact(ix + jx - 1)
                    * pa.powi((ax - ix) as i32)
                    * pb.powi((bx - jx) as i32)
                    / ((2.0 * (alpha + beta)).powi(((ix + jx) / 2) as i32));
            }
        }
    }
    acc
}

/// return the factorial of i
fn factorial(mut i: isize) -> isize {
    let mut prod = 1;
    while i > 1 {
        prod *= i;
        i -= 1;
    }
    prod
}

fn fac(i: isize) -> f64 {
    factorial(i) as f64
}

/// return the binomial coefficient n choose k
fn binom(n: isize, k: isize) -> f64 {
    (factorial(n) / (factorial(k) * factorial(n - k))) as f64
}

impl Index<usize> for Basis {
    type Output = Shell;

    fn index(&self, index: usize) -> &Self::Output {
        self.0.index(index)
    }
}
