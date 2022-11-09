use super::{
    shell::{self, Shell},
    Operator,
};

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

    pub(crate) fn compute(&self, s1: &Shell, s2: &Shell) -> Vec<f64> {
        todo!()
    }
}
