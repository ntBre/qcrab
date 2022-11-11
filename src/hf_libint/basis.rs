use super::{contraction::Contraction, shell::Shell};
use crate::{
    molecule::{Atom, Molecule},
    Dmat,
};
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
                    let mut result = 0.0;
                    for (p1, alpha1) in s1.alpha.iter().enumerate() {
                        for (p2, alpha2) in s2.alpha.iter().enumerate() {
                            let a = s1.origin;
                            let b = s2.origin;
                            let c1 = s1.contr[0].coeff[p1];
                            let c2 = s2.contr[0].coeff[p2];
                            let gammap = alpha1 + alpha2;
                            let oogammap = 1.0 / gammap;
                            let rhop_over_alpha1 = alpha2 * oogammap;
                            let rhop = alpha1 * rhop_over_alpha1;
                            let ab_x = a[0] - b[0];
                            let ab_y = a[1] - b[1];
                            let ab_z = a[2] - b[2];
                            let ab2_x = ab_x * ab_x;
                            let ab2_y = ab_y * ab_y;
                            let ab2_z = ab_z * ab_z;

                            let xyz_pfac: f64 = 1.0;
                            let ovlp_ss_x =
                                f64::exp(-rhop * ab2_x) * xyz_pfac * c1 * c2;
                            let ovlp_ss_y = f64::exp(-rhop * ab2_y) * xyz_pfac;
                            let ovlp_ss_z = f64::exp(-rhop * ab2_z) * xyz_pfac;

                            let p = (*alpha1 * a + *alpha2 * b) / gammap;
                            let pa = p - a;
                            let pb = p - b;

                            let ix = i_helper(l1, l2, pa.x, pb.x, gammap);
                            let iy = i_helper(l1, l2, pa.y, pb.y, gammap);
                            let iz = i_helper(l1, l2, pa.z, pb.z, gammap);

                            if l1 != 0 || l2 != 0 {
                                // println!("result_{p1}{p2}={:.8}", result);
                            }
                            result += ovlp_ss_x
                                * ovlp_ss_y
                                * ovlp_ss_z
                                * ix
                                * iy
                                * iz;
                        }
                    }
                    if l1 != 0 || l2 != 0 {
                        // println!("result={:.8}", result);
                        vec![f64::NAN; n1 * n2]
                    } else {
                        vec![result]
                    }
                };
                // println!(
                //     "setting {},{} to {},{}",
                //     bf1,
                //     bf2,
                //     bf1 + n1 - 1,
                //     bf2 + n2 - 1
                // );
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

fn i_helper(l1: usize, l2: usize, pax: f64, pbx: f64, gammap: f64) -> f64 {
    fn f(k: usize, l1: usize, l2: usize, pax: f64, pbx: f64) -> f64 {
        fn binom(n: isize, k: isize) -> f64 {
            fn fact(mut i: isize) -> isize {
                let mut prod = 1;
                while i > 1 {
                    prod *= i;
                    i -= 1;
                }
                prod
            }
            (fact(n) / (fact(k) * fact(n - k))) as f64
        }
        let k = k as isize;
        let l1 = l1 as isize;
        let l2 = l2 as isize;
        let mut res = 0.0;
        for i in 0..=l1 {
            for j in 0..=l2 {
                if i + j <= k {
                    res += pax.powi((l1 - i) as i32)
                        * binom(l1, i)
                        * pbx.powi((l2 - j) as i32)
                        * binom(l2, j);
                }
            }
        }
        res
    }
    let mut ret = 0.0;
    for i in 0..=l1 + l2 {
        let v = if i < 1 {
            1.0
        } else {
            super::DF_KMINUS1[2 * i - 1] as f64
        };
        ret += f(2 * i, l1, l2, pax, pbx) * v / (2.0 * gammap).powi(i as i32)
            * (PI / gammap).sqrt();
    }
    ret
}

impl Index<usize> for Basis {
    type Output = Shell;

    fn index(&self, index: usize) -> &Self::Output {
        self.0.index(index)
    }
}
