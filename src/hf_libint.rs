use crate::{
    hf::nuclear_repulsion,
    molecule::{Atom, Molecule},
    Dmat, Vec3,
};

struct Contraction {
    /// angular momentum
    l: usize,

    /// spherical(?)
    pure: bool,

    /// contraction coefficients
    coeff: Vec<f64>,
}

impl Contraction {
    fn new(l: usize, pure: bool, coeff: Vec<f64>) -> Self {
        Self { l, pure, coeff }
    }

    const fn cartesian_size(&self) -> usize {
        (self.l + 1) * (self.l + 2) / 2
    }

    const fn size(&self) -> usize {
        if self.pure {
            2 * self.l + 1
        } else {
            (self.l + 1) * (self.l + 2) / 2
        }
    }
}

/// double factorial of k-1 = (k-1)!!
const DF_KMINUS1: [i64; 31] = [
    1,
    1,
    1,
    2,
    3,
    8,
    15,
    48,
    105,
    384,
    945,
    3840,
    10395,
    46080,
    135135,
    645120,
    2027025,
    10321920,
    34459425,
    185794560,
    654729075,
    3715891200,
    13749310575,
    81749606400,
    316234143225,
    1961990553600,
    7905853580625,
    51011754393600,
    213458046676875,
    1428329123020800,
    6190283353629375,
];

/// A shell in a basis set
struct Shell {
    /// exponents of primitive Gaussians
    alpha: Vec<f64>,

    /// see [Contraction]
    contr: Vec<Contraction>,

    /// origin of the shell
    origin: Vec3,

    /// maximum ln of (absolute) contraction coefficient for each primitive
    max_ln_coeff: Vec<f64>,
}

impl Shell {
    fn new(
        alpha: Vec<f64>,
        mut contr: Vec<Contraction>,
        origin: Vec3,
        embed_norm: bool,
    ) -> Self {
        Self {
            max_ln_coeff: if embed_norm {
                Self::renorm(&alpha, &mut contr)
            } else {
                Self::update_max_ln_coeff(&alpha, &contr)
            },
            alpha,
            contr,
            origin,
        }
    }

    /// embeds normalization constants into contraction coefficients.
    fn renorm(alpha: &[f64], contr: &mut [Contraction]) -> Vec<f64> {
        const SQRT_PI_CUBED: f64 = 5.568_327_996_831_708;
        let np = alpha.len();
        for c in contr.iter_mut() {
            // size of df_kminus1 array
            assert!(c.l <= 15);
            for (i, p) in alpha.iter().enumerate() {
                assert!(*p >= 0.0);
                if *p != 0.0 {
                    let two_alpha = 2.0 * p;
                    let two_alpha_to_am32 =
                        two_alpha.powi(c.l as i32 + 1) * two_alpha.sqrt();
                    let normalization_factor = f64::sqrt(
                        f64::powi(2.0, c.l as i32) * two_alpha_to_am32
                            / (SQRT_PI_CUBED * DF_KMINUS1[2 * c.l] as f64),
                    );
                    c.coeff[i] *= normalization_factor;
                }
            }

            // this should be do_enforce_unit_normalization(), but it's set by
            // some kind of flags, and the default is to do it
            if true {
                let mut norm = 0.0;
                for p in 0..np {
                    for q in 0..=p {
                        let gamma = alpha[p] + alpha[q];
                        let a = if p == q { 1.0 } else { 2.0 };
                        norm += a
                            * DF_KMINUS1[2 * c.l] as f64
                            * SQRT_PI_CUBED
                            * c.coeff[p]
                            * c.coeff[q]
                            / (f64::powi(2.0, c.l as i32)
                                * f64::powi(gamma, c.l as i32 + 1)
                                * gamma.sqrt());
                    }
                }
                let fac = 1.0 / norm.sqrt();
                for p in 0..np {
                    c.coeff[p] *= fac;
                }
            }
        }

        Self::update_max_ln_coeff(alpha, contr)
    }

    // NOTE: if there can be zero contractions, we might want to initialize
    // max_ln_c to -std::numeric_limits<real_t>::max() in C++ terms
    fn update_max_ln_coeff(alpha: &[f64], contr: &[Contraction]) -> Vec<f64> {
        let nprim = alpha.len();
        let mut ret = vec![0.0; nprim];
        for (p, v) in ret.iter_mut().enumerate() {
            let mut max_ln_c = contr[0].coeff[p].abs().ln();
            for c in contr.iter().skip(1) {
                max_ln_c = max_ln_c.max(c.coeff[p].abs().ln());
            }
            *v = max_ln_c;
        }
        ret
    }
}

/// STO-3G basis set
///
/// cite: W. J. Hehre, R. F. Stewart, and J. A. Pople, The Journal of Chemical
/// Physics 51, 2657 (1969) doi: 10.1063/1.1672392
///
/// obtained from https://bse.pnl.gov/bse/portal via libint
fn make_sto3g_basis(mol: &Molecule) -> Vec<Shell> {
    let mut shells = Vec::new();
    for Atom {
        atomic_number,
        coord,
    } in &mol.atoms
    {
        match atomic_number {
            1 => shells.push(Shell::new(
                vec![3.425250910, 0.623913730, 0.168855400],
                vec![Contraction::new(
                    0,
                    false,
                    vec![0.15432897, 0.53532814, 0.44463454],
                )],
                *coord,
                true,
            )),
            // oxygen
            8 => shells.extend([
                //
                Shell::new(
                    vec![130.709320000, 23.808861000, 6.443608300],
                    vec![Contraction::new(
                        0,
                        false,
                        vec![0.15432897, 0.53532814, 0.44463454],
                    )],
                    *coord,
                    true,
                ),
                Shell::new(
                    vec![5.033151300, 1.169596100, 0.380389000],
                    vec![Contraction::new(
                        0,
                        false,
                        vec![-0.09996723, 0.39951283, 0.70011547],
                    )],
                    *coord,
                    true,
                ),
                Shell::new(
                    vec![5.033151300, 1.169596100, 0.380389000],
                    vec![Contraction::new(
                        1,
                        false,
                        vec![0.15591627, 0.60768372, 0.39195739],
                    )],
                    *coord,
                    true,
                ),
            ]),
            _ => panic!("unsupported element {atomic_number}"),
        }
        // shells.push(
    }
    shells
}

/// types of Operators supported by Engine
enum Operator {
    /// overlap
    Overlap,

    /// electronic kinetic energy, -1/2∇²
    Kinetic,

    /// Coulomb potential due to point charges
    Nuclear,

    /// 2-body Coulomb operator, 1/r₁₂
    Coulomb,
}

fn compute_1body_ints(shells: &[Shell], overlap: Operator) -> Dmat {
    todo!()
}

#[test]
fn hf_libint() {
    let mol = Molecule::load("testfiles/h2o/STO-3G/geom.dat");
    let enuc = nuclear_repulsion(&mol);
    let shells = make_sto3g_basis(&mol);

    let s = compute_1body_ints(&shells, Operator::Overlap);
    let t = compute_1body_ints(&shells, Operator::Kinetic);
    let v = compute_1body_ints(&shells, Operator::Kinetic);
}
