use std::ops::Index;

use crate::{
    hf::nuclear_repulsion,
    molecule::{Atom, Molecule},
    Vec3,
};

use self::{contraction::Contraction, engine::Engine, shell::Shell};

mod basis;
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

#[test]
fn hf_libint() {
    let mol = Molecule::load("testfiles/h2o/STO-3G/geom.dat");
    let enuc = nuclear_repulsion(&mol);
    let shells = basis::Basis::sto3g(&mol);
    let nao: usize = shells.0.iter().map(|s| s.size()).sum();

    let q = mol
        .atoms
        .iter()
        .map(
            |&Atom {
                 atomic_number,
                 coord,
             }| (atomic_number, coord),
        )
        .collect();

    let s = shells.compute_1body_ints(Operator::Overlap);
    let t = shells.compute_1body_ints(Operator::Kinetic);
    let v = shells.compute_1body_ints(Operator::Nuclear { q });
}
