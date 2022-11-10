use crate::{
    hf::{nuclear_repulsion, overlap_integrals},
    molecule::Molecule,
};

pub(crate) mod basis;
mod contraction;
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

#[allow(unused)]
fn libint() {
    let mol = Molecule::load("testfiles/h2o/STO-3G/geom.dat");
    let _enuc = nuclear_repulsion(&mol);
    let shells = basis::Basis::sto3g(&mol);
    let _nao: usize = shells.0.iter().map(|s| s.size()).sum();

    let s = overlap_integrals("testfiles/h2o/STO-3G/s.dat");
    println!("loaded s={:.8}", s);
    let s = shells.overlap_ints();
    println!("computed s={:.8}", s);
    // let t = shells.compute_1body_ints(Operator::Kinetic);
    // let v = shells.compute_1body_ints(Operator::nuclear(&mol));
}

#[test]
fn hf_libint() {
    libint()
}
