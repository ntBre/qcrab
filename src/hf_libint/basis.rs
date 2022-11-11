use super::{contraction::Contraction, shell::Shell};
use crate::{
    hf_libint::dfact,
    molecule::{Atom, Molecule},
    Dmat,
};
use itertools::Itertools;
use std::{f64::consts::PI, ops::Index};

pub(crate) struct Basis(pub(crate) Vec<Shell>);

impl Basis {
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
            let n1 = self[s1].size();
            for s2 in 0..=s1 {
                let bf2 = shell2bf[s2];
                let n2 = self[s2].size();
                let buf = {
                    let s1: &Shell = &self[s1];
                    let s2: &Shell = &self[s2];
                    assert!(
                        s1.ncontr() == 1 && s2.ncontr() == 1,
                        "generally-contracted shells not yet supported"
                    );

                    let l1 = s1.contr[0].l;
                    let l2 = s2.contr[0].l;
                    let mut results = vec![0.0; n1 * n2];
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

                            let ls1 = match l1 {
                                0 => vec![(0, 0, 0)],
                                1 => vec![(1, 0, 0), (0, 1, 0), (0, 0, 1)],
                                _ => panic!("unmatched l value {l1}"),
                            };
                            let ls2 = match l2 {
                                0 => vec![(0, 0, 0)],
                                1 => vec![(1, 0, 0), (0, 1, 0), (0, 0, 1)],
                                _ => panic!("unmatched l value {l1}"),
                            };

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
                                results[i] += ovlp_ss_x
                                    * ovlp_ss_y
                                    * ovlp_ss_z
                                    * ix
                                    * iy
                                    * iz;
                            }
                        }
                    }
                    results
                };
                let mut ij = 0;
                for i in bf1..bf1 + n1 {
                    for j in bf2..bf2 + n2 {
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
