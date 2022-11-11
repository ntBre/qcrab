use nalgebra::dmatrix;

use crate::{
    hf::{nuclear_repulsion, overlap_integrals},
    molecule::Molecule,
    Dmat,
};

use std::f64::NAN;

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

fn nan_cmp(a: &Dmat, b: &Dmat, eps: f64) -> bool {
    let (ar, ac) = a.shape();
    let (br, bc) = b.shape();
    if ar != br || ac != bc {
        return false;
    }
    for i in 0..ar {
        for j in 0..ac {
            // if they're not equal and not both NaN
            if (a[(i, j)] - b[(i, j)]).abs() > eps {
                if !(a[(i, j)].is_nan() && b[(i, j)].is_nan()) {
                    return false;
                }
            }
        }
    }
    true
}

#[allow(unused)]
fn libint() {
    let mol = Molecule::load("testfiles/h2o/STO-3G/geom.dat");
    let _enuc = nuclear_repulsion(&mol);
    let shells = basis::Basis::sto3g(&mol);
    let _nao: usize = shells.0.iter().map(|s| s.size()).sum();

    let s = overlap_integrals("testfiles/h2o/STO-3G/s.dat");
    println!("loaded s={:.8}", s);
    let s = shells.overlap_ints();
    let want_s = dmatrix![
    1.00000000, 0.23670394,        NAN,        NAN,        NAN, 0.03840560, 0.03840560;
    0.23670394, 1.00000000,        NAN,        NAN,        NAN, 0.38613884, 0.38613884;
           NAN,        NAN,        NAN,        NAN,        NAN,        NAN,        NAN;
           NAN,        NAN,        NAN,        NAN,        NAN,        NAN,        NAN;
           NAN,        NAN,        NAN,        NAN,        NAN,        NAN,        NAN;
    0.03840560, 0.38613884,        NAN,        NAN,        NAN, 1.00000000, 0.18175989;
    0.03840560, 0.38613884,        NAN,        NAN,        NAN, 0.18175989, 1.00000000;
     ];
    assert!(nan_cmp(&s, &want_s, 1e-8));
    println!("computed s={:.8}", s);
    // let t = shells.compute_1body_ints(Operator::Kinetic);
    // let v = shells.compute_1body_ints(Operator::nuclear(&mol));
}

#[test]
fn hf_libint() {
    libint()
}
