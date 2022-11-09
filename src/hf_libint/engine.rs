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
        }
    }

    fn nparams(&self) -> usize {
        match self.oper {
            Operator::Overlap => 1,
            Operator::Kinetic => 1,
            Operator::Nuclear { q } => q.len(),
            Operator::Coulomb => 1,
        }
    }

    /// libint lumps this and [compute2] together, but they have different
    /// arities, so I'll just keep them separate
    pub(crate) fn compute1(&self, s1: &Shell, s2: &Shell) -> Vec<f64> {
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
            let mut p12 = 0.0;
            for p1 in 0..nprim1 {
                for p2 in 0..nprim2 {
                    primdata.push(self.compute_primdata(s1, s2, p1, p2, pset));
                }
            }
        }

        // TODO fix this after implementing primdata

        // if compute_directly {
        //     let mut result = primdata[0];
        //     match self.oper {
        //         Operator::Overlap => for p12 in 0..contrdepth {},
        //         Operator::Kinetic => todo!(),
        //         Operator::Nuclear { q } => todo!(),
        //         Operator::Coulomb => todo!(),
        //     }
        // }
        todo!()
    }

    pub(crate) fn compute_primdata(
        &self,
        s1: &Shell,
        s2: &Shell,
        p1: usize,
        p2: usize,
        pset: usize,
    ) -> _ {
        todo!()
    }
}
