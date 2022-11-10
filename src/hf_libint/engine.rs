use super::{shell::Shell, Operator};

struct Params;

/// an enum eventually I think
struct BraKet;

struct ScreeningMethod;

/// `Engine` computes integrals of operators (or operator sets) specified by
/// combination of `Operator` and `BraKet`. NOTE the C++ docs say that only
/// one-contraction Shells are supported, but I don't think this is a
/// problem for normal basis sets
pub(crate) struct Engine {
    /// the operator
    oper: Operator,

    /// the maximum number of primitives per contracted Gaussian (must be
    /// greater than 0)
    max_nprim: usize,

    /// maximum angular momentum of a Gaussian shell
    max_l: usize,

    /// if not 0, compute geometric derivatives of Gaussian integrals of
    /// order `deriv_order`
    deriv_order: usize,

    /// numerical precision, defaults to f64::EPSILON. only used for
    /// two-body integrals
    precision: f64,

    /// a value of type Engine::operator_traits<oper>::oper_params_type (!)
    /// specifying the parameters of the operator set, e.g. position and
    /// magnitude of charges creating the Coulomb potential for
    /// [Operator::Nuclear]. For most values of `oper`, this is not needed.
    /// NOTE: in Rust these will probably fields on the Operator enum
    params: Params,

    /// a value of the BraKet type
    braket: BraKet,

    /// not documented in libint
    screening_method: ScreeningMethod,

    /// scale the target integrals by this factor. defaults to 1.0
    scale: f64,
}

impl Engine {
    pub(crate) fn new(
        oper: Operator,
        max_nprim: usize,
        max_l: usize,
        deriv_order: usize,
    ) -> Self {
        Self {
            oper,
            max_nprim,
            max_l,
            deriv_order,
            precision: f64::EPSILON,
            params: Params,
            braket: BraKet,
            screening_method: ScreeningMethod,
            scale: 1.0,
        }
    }

    fn nparams(&self) -> usize {
        match &self.oper {
            Operator::Overlap => 1,
            Operator::Kinetic => 1,
            Operator::Nuclear { q } => q.len(),
            Operator::Coulomb => 1,
        }
    }

    /// libint lumps this and [compute2] together, but they have different
    /// arities, so I'll just keep them separate
    pub(crate) fn compute1(&self, s1: &Shell, s2: &Shell) -> Vec<f64> {
        let mut ret = Vec::new();
        assert!(
            s1.ncontr() == 1 && s2.ncontr() == 1,
            "generally-contracted shells not yet supported"
        );

        let is_nuc = matches!(self.oper, Operator::Nuclear { q: _ });
        let l1 = s1.contr[0].l;
        let l2 = s2.contr[0].l;

        // we can skip the check on params because the nuclear charges are
        // required in the Nuclear operator

        let n1 = s1.size();
        let n2 = s2.size();
        let ncart1 = s1.cartesian_size();
        let ncart2 = s2.cartesian_size();
        let ncart12 = ncart1 * ncart2;

        // there are a bunch of checks here I think we can skip
        // (engine.impl.h:201) also skipping the scratch stuff

        let lmax = l1.max(l2);

        // N.B. for l=0 no need to transform to solid harmonics this is a
        // workaround for the corner case of oper_ == Operator::*nuclear, and
        // solid harmonics (s|s) integral ... beware the integral storage state
        // machine
        let tform_to_solids =
            (s1.contr[0].pure || s2.contr[0].pure) && lmax != 0;

        // simple (s|s) ints will be computed directly and accumulated in the first
        // element of stack
        let compute_directly = lmax == 0
            && self.deriv_order == 0
            && (matches!(self.oper, Operator::Overlap) || is_nuc);

        // TODO implement primdata, absolutely insane definition in libint
        let mut primdata = Vec::new();

        let nprim1 = s1.nprim();
        let nprim2 = s2.nprim();

        // loop over accumulation batches. might not need this since I don't
        // really want to batch anything for now
        for pset in 0..self.nparams() {
            let mut p12 = 0;
            for p1 in 0..nprim1 {
                for p2 in 0..nprim2 {
                    let p = self.compute_primdata(s1, s2, p1, p2, pset);
                    primdata.push(p);
                    p12 += 1;
                }
            }
            let contrdepth = p12;

            if compute_directly {
                let mut result = 0.0;
                match &self.oper {
                    Operator::Overlap => {
                        for p in &primdata {
                            result += p.ovlp_ss_x * p.ovlp_ss_y * p.ovlp_ss_z;
                        }
                        ret.push(result);
                    }
                    Operator::Nuclear { q } => todo!(),
                    _ => {}
                }
            } else {
                // default braket for oper of rank 1 is x_x, for 2 is xx_xx
                // let buildfnidx = s1.contr[0].l * hard_lmax + s2.contr[0].l;
                // todo!()
                ret.push(f64::NAN);
            }
        }

        if tform_to_solids {
            // todo!();
        }
        // todo!();
        ret
    }

    pub(crate) fn compute_primdata(
        &self,
        s1: &Shell,
        s2: &Shell,
        p1: usize,
        p2: usize,
        pset: usize,
    ) -> Primdata {
        // return a libint_t which depends on something to be defined
        let a = s1.origin;
        let b = s2.origin;
        let alpha1 = s1.alpha[p1];
        let alpha2 = s2.alpha[p2];
        let c1 = s1.contr[0].coeff[p1];
        let c2 = s2.contr[0].coeff[p2];
        let gammap = alpha1 + alpha2;
        let oogammap = 1.0 / gammap;
        let rhop_over_alpha1 = alpha2 * oogammap;
        let rhop = alpha1 * rhop_over_alpha1;
        let px = (alpha1 * a[0] + alpha2 * b[0]) * oogammap;
        let py = (alpha1 * a[1] + alpha2 * b[1]) * oogammap;
        let pz = (alpha1 * a[2] + alpha2 * b[2]) * oogammap;
        let ab_x = a[0] - b[0];
        let ab_y = a[1] - b[1];
        let ab_z = a[2] - b[2];
        let ab2_x = ab_x * ab_x;
        let ab2_y = ab_y * ab_y;
        let ab2_z = ab_z * ab_z;

        let is_nuc = self.oper.is_nuclear();

        let l1 = s1.contr[0].l;
        let l2 = s2.contr[0].l;

        // NOTE this is a nother operator variant that we don't have yet
        let is_sphemultipole = false;
        let use_hrr = (is_nuc || is_sphemultipole) && l1 > 0 && l2 > 0;
        let hrr_ket_to_bra = l1 >= l2;

        let mut ret = Primdata::default();

        if use_hrr {
            if hrr_ket_to_bra {
                ret.ab_x = ab_x;
                ret.ab_y = ab_y;
                ret.ab_z = ab_z;
            } else {
                ret.ba_x = -ab_x;
                ret.ba_y = -ab_y;
                ret.ba_z = -ab_z;
            }
        }

        ret.pa_x = px - a[0];
        ret.pa_y = py - a[1];
        ret.pa_z = pz - a[2];

        ret.pb_x = px - b[0];
        ret.pb_y = py - b[1];
        ret.pb_z = pz - b[2];

        // if oper == emultipole[123] ...

        // if oper == sphemultipole ...

        ret.oo2z = oogammap;
        const SQRT_PI: f64 = 1.772_453_850_905_516;
        let xyz_pfac: f64 = SQRT_PI * f64::sqrt(oogammap);
        ret.ovlp_ss_x =
            f64::exp(-rhop * ab2_x) * xyz_pfac * c1 * c2 * self.scale;
        ret.ovlp_ss_y = f64::exp(-rhop * ab2_y) * xyz_pfac;
        ret.ovlp_ss_z = f64::exp(-rhop * ab2_z) * xyz_pfac;

        if self.oper.is_kinetic() || self.deriv_order > 0 {
            ret.two_alpha0_bra = 2.0 * alpha1;
            ret.two_alpha0_ket = 2.0 * alpha2;
        }

        if self.oper.is_nuclear() {
            todo!("engine.impl.h:1025")
        }
        ret
    }
}

// I think some of these are enums/mutually exclusive, but I'm not sure yet
#[derive(Default)]
pub(crate) struct Primdata {
    ab_x: f64,
    ab_y: f64,
    ab_z: f64,
    //
    ba_x: f64,
    ba_y: f64,
    ba_z: f64,
    //
    pa_x: f64,
    pa_y: f64,
    pa_z: f64,
    //
    pb_x: f64,
    pb_y: f64,
    pb_z: f64,
    //
    oo2z: f64,
    //
    ovlp_ss_x: f64,
    ovlp_ss_y: f64,
    ovlp_ss_z: f64,
    //
    two_alpha0_bra: f64,
    two_alpha0_ket: f64,
}
