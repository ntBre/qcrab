use super::{contraction::Contraction, shell::Shell};
use crate::{
    hf_libint::dfact,
    molecule::{Atom, Molecule},
    Dmat,
};
use std::{f64::consts::PI, fs::read_to_string, ops::Index, path::Path};

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
            for s2 in 0..=s1 {
                let bf2 = shell2bf[s2];
                let s1 = &self[s1];
                let s2 = &self[s2];
                assert!(
                    s1.ncontr() == 1 && s2.ncontr() == 1,
                    "generally-contracted shells not yet supported"
                );

                let l1 = s1.contr[0].l;
                let l2 = s2.contr[0].l;
                let ls1 = l_perm(l1);
                let ls2 = l_perm(l2);
                for (p1, alpha) in s1.alpha.iter().enumerate() {
                    for (p2, beta) in s2.alpha.iter().enumerate() {
                        let a = s1.origin;
                        let b = s2.origin;
                        let c1 = s1.contr[0].coeff[p1];
                        let c2 = s2.contr[0].coeff[p2];
                        let gamma = alpha + beta;
                        // in-lining this fails the test...
                        let oogammap = 1.0 / gamma;
                        let rhop = alpha * beta * oogammap;
                        let ab_x = a[0] - b[0];
                        let ab_y = a[1] - b[1];
                        let ab_z = a[2] - b[2];
                        let ab2_x = ab_x * ab_x;
                        let ab2_y = ab_y * ab_y;
                        let ab2_z = ab_z * ab_z;

                        let eabx = f64::exp(-rhop * ab2_x) * c1 * c2;
                        let eaby = f64::exp(-rhop * ab2_y);
                        let eabz = f64::exp(-rhop * ab2_z);

                        let p = (*alpha * a + *beta * b) / gamma;
                        let pa = p - a;
                        let pb = p - b;

                        if s1.contr[0].pure || s2.contr[0].pure {
                            panic!("can't handle spherical harmonics yet");
                        }

                        let pfac = (PI / (alpha + beta)).sqrt();
                        // actually best to keep these here, not outer loops
                        for (i, (l1x, l1y, l1z)) in ls1.iter().enumerate() {
                            for (j, (l2x, l2y, l2z)) in ls2.iter().enumerate() {
                                let ix = pfac
                                    * s_xyz(
                                        *l1x, *l2x, pa.x, pb.x, alpha, beta,
                                    );
                                let iy = pfac
                                    * s_xyz(
                                        *l1y, *l2y, pa.y, pb.y, alpha, beta,
                                    );
                                let iz = pfac
                                    * s_xyz(
                                        *l1z, *l2z, pa.z, pb.z, alpha, beta,
                                    );
                                result[(bf1 + i, bf2 + j)] +=
                                    eabx * eaby * eabz * ix * iy * iz;
                            }
                        }
                    }
                }
            }
        }
        result.fill_upper_triangle_with_lower_triangle();
        result
    }
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
