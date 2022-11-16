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

fn dfact(n: isize) -> f64 {
    if n == -1 {
        1.0
    } else {
        DF_KMINUS1[n as usize] as f64
    }
}

#[allow(unused)]
fn libint() {
    let mol = Molecule::load("testfiles/h2o/STO-3G/geom.dat");
    let _enuc = nuclear_repulsion(&mol);
    let shells = basis::Basis::sto3g(&mol);
    let _nao: usize = shells.shells.iter().map(|s| s.size()).sum();

    let want = overlap_integrals("testfiles/h2o/STO-3G/s.dat");
    let got = shells.overlap_ints();
    approx::assert_abs_diff_eq!(got, want, epsilon = 1e-13);
    // let t = shells.compute_1body_ints(Operator::Kinetic);
    // let v = shells.compute_1body_ints(Operator::nuclear(&mol));
}

#[test]
fn hf_libint() {
    libint()
}

#[test]
fn overlap() {
    let mol = Molecule::load("testfiles/h2o/STO-3G/geom.dat");
    let shells = basis::Basis::sto3g(&mol);
    let want = overlap_integrals("testfiles/h2o/STO-3G/s.dat");
    let got = shells.overlap_ints();
    if approx::abs_diff_ne!(got, want, epsilon = 1e-13) {
        println!("got={:.8}", got);
        println!("want={:.8}", want);
        panic!("differ by {:.2e}", (got - want).abs().max());
    }
}

#[test]
fn kinetic() {
    let mol = Molecule::load("testfiles/h2o/STO-3G/geom.dat");
    let shells = basis::Basis::sto3g(&mol);
    let want = crate::hf::kinetic_integrals("testfiles/h2o/STO-3G/t.dat");
    let got = shells.kinetic_ints();
    if approx::abs_diff_ne!(got, want, epsilon = 1e-13) {
        println!("got={:.8}", got);
        println!("want={:.8}", want);
        println!("diff={:.8}", &got - &want);
        panic!("differ by {:.2e}", (got - want).abs().max());
    }
}

#[test]
fn kinetic_dz() {
    let mol = Molecule::load("testfiles/h2o/DZ/geom.dat");
    let shells = basis::Basis::load("basis_sets/dz.json", &mol);
    let want = crate::hf::kinetic_integrals("testfiles/h2o/DZ/t.dat");
    let got = shells.kinetic_ints();
    if approx::abs_diff_ne!(got, want, epsilon = 1e-13) {
        println!("got={:.8}", got);
        println!("want={:.8}", want);
        println!("diff={:.8}", &got - &want);
        panic!("differ by {:.2e}", (got - want).abs().max());
    }
}

#[cfg(test)]
mod benches {
    extern crate test;

    use test::Bencher;

    use crate::molecule::Molecule;

    use super::basis;

    #[bench]
    fn overlap(b: &mut Bencher) {
        let mol = Molecule::load("testfiles/h2o/STO-3G/geom.dat");
        let shells = basis::Basis::sto3g(&mol);
        b.iter(|| shells.overlap_ints());
    }
}
