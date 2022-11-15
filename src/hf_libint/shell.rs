use super::contraction::Contraction;
use super::DF_KMINUS1;
use crate::Vec3;

/// A shell in a basis set
#[derive(PartialEq, Debug)]
pub(crate) struct Shell {
    /// exponents of primitive Gaussians
    pub(crate) alpha: Vec<f64>,

    /// see [Contraction]
    pub(crate) contr: Vec<Contraction>,

    /// origin of the shell
    pub(crate) origin: Vec3,
}

impl Shell {
    pub(crate) fn new(
        alpha: Vec<f64>,
        mut contr: Vec<Contraction>,
        origin: Vec3,
    ) -> Self {
        Self::renorm(&alpha, &mut contr);
        Self {
            alpha,
            contr,
            origin,
        }
    }

    pub(crate) fn ncontr(&self) -> usize {
        self.contr.len()
    }

    #[allow(unused)]
    pub(crate) fn nprim(&self) -> usize {
        self.alpha.len()
    }

    /// embeds normalization constants into contraction coefficients.
    pub(crate) fn renorm(alpha: &[f64], contr: &mut [Contraction]) -> Vec<f64> {
        pub(crate) const SQRT_PI_CUBED: f64 = 5.568_327_996_831_708;
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
    pub(crate) fn update_max_ln_coeff(
        alpha: &[f64],
        contr: &[Contraction],
    ) -> Vec<f64> {
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

    pub(crate) fn cartesian_size(&self) -> usize {
        self.contr.iter().map(|c| c.cartesian_size()).sum()
    }

    pub(crate) fn size(&self) -> usize {
        self.contr.iter().map(|c| c.size()).sum()
    }
}
