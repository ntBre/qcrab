use crate::{
    hf::nuclear_repulsion,
    molecule::{Atom, Molecule},
    Vec3,
};

pub(crate) mod basis;
mod contraction;
mod engine;
mod shell;

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

/// types of Operators supported by Engine. TODO might not be copy if these get
/// extensive fields
#[derive(Clone)]
pub(crate) enum Operator {
    /// overlap
    Overlap,

    /// electronic kinetic energy, -1/2∇²
    Kinetic,

    /// Coulomb potential due to point charges
    Nuclear { q: Vec<(usize, Vec3)> },

    /// 2-body Coulomb operator, 1/r₁₂
    Coulomb,
}

impl Operator {
    pub(crate) fn rank(&self) -> usize {
        match self {
            Operator::Overlap | Operator::Kinetic => 1,
            Operator::Nuclear { q: _ } => 1,
            Operator::Coulomb => 2,
        }
    }

    /// return an Operator::Nuclear with the charges generated from `mol`
    pub(crate) fn nuclear(mol: &Molecule) -> Self {
        Self::Nuclear {
            q: mol
                .atoms
                .iter()
                .map(
                    |&Atom {
                         atomic_number,
                         coord,
                     }| (atomic_number, coord),
                )
                .collect(),
        }
    }

    /// Returns `true` if the operator is [`Nuclear`].
    ///
    /// [`Nuclear`]: Operator::Nuclear
    #[must_use]
    pub(crate) fn is_nuclear(&self) -> bool {
        matches!(self, Self::Nuclear { .. })
    }

    /// Returns `true` if the operator is [`Kinetic`].
    ///
    /// [`Kinetic`]: Operator::Kinetic
    #[must_use]
    pub(crate) fn is_kinetic(&self) -> bool {
        matches!(self, Self::Kinetic)
    }
}

fn libint() {
    let mol = Molecule::load("testfiles/h2o/STO-3G/geom.dat");
    let enuc = nuclear_repulsion(&mol);
    let shells = basis::Basis::sto3g(&mol);
    let nao: usize = shells.0.iter().map(|s| s.size()).sum();

    let s = shells.compute_1body_ints(Operator::Overlap);
    let t = shells.compute_1body_ints(Operator::Kinetic);
    let v = shells.compute_1body_ints(Operator::nuclear(&mol));
}

#[test]
fn hf_libint() {
    libint()
}
