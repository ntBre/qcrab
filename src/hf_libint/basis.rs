use super::{contraction::Contraction, shell::Shell};
use crate::{
    hf_libint::dfact,
    molecule::{Atom, Molecule},
    Dmat,
};
use std::{f64::consts::PI, fs::read_to_string, ops::Index, path::Path};

#[derive(PartialEq, Debug)]
pub(crate) struct Basis {
    pub(crate) shells: Vec<Shell>,
}

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
    let want = Basis {
        shells: vec![
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
        ],
    };
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
        Self { shells }
    }

    pub(crate) fn nbasis(&self) -> usize {
        self.shells.iter().map(|s| s.size()).sum()
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
        Self { shells }
    }

    #[allow(unused)]
    pub(crate) fn len(&self) -> usize {
        self.shells.len()
    }

    pub(crate) fn overlap_ints(&self) -> Dmat {
        let n = self.nbasis();
        let mut result = Dmat::zeros(n, n);

        let mut ls = Vec::new();
        let mut ss = Vec::new();
        for (i, shell) in self.shells.iter().enumerate() {
            let l = match shell.contr[0].l {
                0 => vec![(0, 0, 0)],
                1 => vec![(1, 0, 0), (0, 1, 0), (0, 0, 1)],
                2 => vec![
                    (2, 0, 0),
                    (0, 2, 0),
                    (0, 0, 2),
                    (0, 1, 1),
                    (1, 1, 0),
                    (1, 0, 1),
                ],
                _ => panic!("unmatched l value {}", shell.contr[0].l),
            };
            let s = vec![i; l.len()];
            ls.extend(l);
            ss.extend(s);
        }

        // loop over orbitals
        for (i, (l1x, l1y, l1z)) in ls.iter().enumerate() {
            for (j, (l2x, l2y, l2z)) in ls[..=i].iter().enumerate() {
                let s1 = &self.shells[ss[i]];
                let s2 = &self.shells[ss[j]];
                let a = s1.origin;
                let b = s2.origin;
                let ab = a - b;
                let ab2 = ab.dot(&ab);
                let mut res = 0.0;
                // loop over primitives within shells
                for (p1, alpha) in s1.alpha.iter().enumerate() {
                    for (p2, beta) in s2.alpha.iter().enumerate() {
                        let c1 = s1.contr[0].coeff[p1];
                        let c2 = s2.contr[0].coeff[p2];
                        let gamma = alpha + beta;
                        let rhop = alpha * beta / gamma;
                        let eabx = f64::exp(-rhop * ab2) * c1 * c2;
                        let p = (*alpha * a + *beta * b) / gamma;
                        let pa = p - a;
                        let pb = p - b;

                        if s1.contr[0].pure || s2.contr[0].pure {
                            panic!("can't handle spherical harmonics yet");
                        }

                        let pfac = (PI / (alpha + beta)).powf(1.5);
                        let sx = s_xyz(*l1x, *l2x, pa.x, pb.x, alpha, beta);
                        let sy = s_xyz(*l1y, *l2y, pa.y, pb.y, alpha, beta);
                        let sz = s_xyz(*l1z, *l2z, pa.z, pb.z, alpha, beta);
                        res += pfac * eabx * sx * sy * sz;
                    }
                }
                result[(i, j)] = res;
                result[(j, i)] = res;
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
        self.shells.index(index)
    }
}
