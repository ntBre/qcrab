use std::path::Path;

use approx::assert_abs_diff_eq;
use nalgebra::dmatrix;

use crate::eri::Eri;
use crate::hf::*;
use crate::hf_libint::basis::Basis;
use crate::molecule::*;

#[test]
fn test_enuc() {
    let mol = Molecule::load("testfiles/h2o/STO-3G/geom.dat");
    let got = nuclear_repulsion(&mol);
    let want = 8.002_367_061_810_45;
    assert_abs_diff_eq!(got, want, epsilon = 1e-12);
}

#[test]
fn test_core_hamiltonian() {
    // this tests T and V loading
    let t = kinetic_integrals("testfiles/h2o/STO-3G/t.dat");
    let v = attraction_integrals("testfiles/h2o/STO-3G/v.dat");
    let got = t + v;
    let want = dmatrix![
        -32.5773954,  -7.5788328,   0.0000000,  -0.0144738,   0.0000000,
    -1.2401023,  -1.2401023;
        -7.5788328,  -9.2009433,   0.0000000,  -0.1768902,   0.0000000,
    -2.9067098,  -2.9067098;
        0.0000000,   0.0000000,  -7.4588193,   0.0000000,   0.0000000,
    -1.6751501,   1.6751501;
        -0.0144738,  -0.1768902,   0.0000000,  -7.4153118,   0.0000000,
    -1.3568683,  -1.3568683;
        0.0000000,   0.0000000,   0.0000000,   0.0000000,  -7.3471449,
    0.0000000,   0.0000000;
        -1.2401023,  -2.9067098,  -1.6751501,  -1.3568683,   0.0000000,
    -4.5401711,  -1.0711459;
        -1.2401023,  -2.9067098,   1.6751501,  -1.3568683,   0.0000000,
    -1.0711459,  -4.5401711;
       ];
    assert_abs_diff_eq!(got, want, epsilon = 1e-7);
}

#[test]
fn test_orthog() {
    let got =
        build_orthog_matrix(overlap_integrals("testfiles/h2o/STO-3G/s.dat"));
    let want = dmatrix![
        1.0236346,  -0.1368547,  -0.0000000,  -0.0074873,  -0.0000000,
    0.0190279,   0.0190279;
        -0.1368547,   1.1578632,   0.0000000,   0.0721601,   0.0000000,
    -0.2223326,  -0.2223326;
        -0.0000000,   0.0000000,   1.0733148,   0.0000000,  -0.0000000,
    -0.1757583,   0.1757583;
        -0.0074873,   0.0721601,   0.0000000,   1.0383050,   0.0000000,
    -0.1184626,  -0.1184626;
        -0.0000000,   0.0000000,  -0.0000000,   0.0000000,   1.0000000,
    -0.0000000,  -0.0000000;
        0.0190279,  -0.2223326,  -0.1757583,  -0.1184626,  -0.0000000,
    1.1297234,  -0.0625975;
        0.0190279,  -0.2223326,   0.1757583,  -0.1184626,  -0.0000000,
    -0.0625976,   1.1297234;
      ];
    assert_abs_diff_eq!(got, want, epsilon = 1e-7);
}

#[test]
fn test_density() {
    let fock = kinetic_integrals("testfiles/h2o/STO-3G/t.dat")
        + attraction_integrals("testfiles/h2o/STO-3G/v.dat");
    let s12 =
        build_orthog_matrix(overlap_integrals("testfiles/h2o/STO-3G/s.dat"));
    let mol = Molecule::load("testfiles/h2o/STO-3G/geom.dat");
    let got = density(&fock, &s12, mol.nelec());
    let want = dmatrix![
        1.0650117,  -0.2852166,  -0.0000000,  -0.0195534,  -0.0000000,
    0.0334496,   0.0334496;
        -0.2852166,   1.2489657,   0.0000000,   0.1135594,   0.0000000,
    -0.1442809,  -0.1442809;
        -0.0000000,   0.0000000,   1.1258701,  -0.0000000,  -0.0000000,
    -0.1461317,   0.1461317;
        -0.0195534,   0.1135594,  -0.0000000,   1.0660638,   0.0000000,
    -0.0993583,  -0.0993583;
        -0.0000000,   0.0000000,  -0.0000000,   0.0000000,   1.0000000,
    -0.0000000,  -0.0000000;
        0.0334496,  -0.1442809,  -0.1461317,  -0.0993583,  -0.0000000,
    0.0426802,   0.0047460;
        0.0334496,  -0.1442809,   0.1461317,  -0.0993583,  -0.0000000,
    0.0047460,   0.0426802;
      ];
    assert_abs_diff_eq!(got, want, epsilon = 1e-7);
}

#[test]
fn test_scf_init() {
    let hcore = kinetic_integrals("testfiles/h2o/STO-3G/t.dat")
        + attraction_integrals("testfiles/h2o/STO-3G/v.dat");
    let s12 =
        build_orthog_matrix(overlap_integrals("testfiles/h2o/STO-3G/s.dat"));
    let mol = Molecule::load("testfiles/h2o/STO-3G/geom.dat");
    let d = density(&hcore, &s12, mol.nelec());
    let got = energy(&d, &hcore, &hcore);
    let want = -125.842077437699;
    assert_abs_diff_eq!(got, want, epsilon = 1e-12);
}

#[test]
fn test_fock() {
    let hcore = kinetic_integrals("testfiles/h2o/STO-3G/t.dat")
        + attraction_integrals("testfiles/h2o/STO-3G/v.dat");
    let s12 =
        build_orthog_matrix(overlap_integrals("testfiles/h2o/STO-3G/s.dat"));
    let mol = Molecule::load("testfiles/h2o/STO-3G/geom.dat");
    let d = density(&hcore, &s12, mol.nelec());
    let eri = Eri::new("testfiles/h2o/STO-3G/eri.dat");
    let got = fock(&d, &hcore, &eri);
    let want = dmatrix![
        -18.8132695,  -4.8726875,  -0.0000000,  -0.0115290,   0.0000000,
    -0.8067323,  -0.8067323;
        -4.8726875,  -1.7909029,  -0.0000000,  -0.1808692,   0.0000000,
    -0.5790557,  -0.5790557;
        -0.0000000,  -0.0000000,   0.1939644,   0.0000000,   0.0000000,
    -0.1708886,   0.1708886;
        -0.0115290,  -0.1808692,   0.0000000,   0.2391247,   0.0000000,
    -0.1828683,  -0.1828683;
        0.0000000,   0.0000000,   0.0000000,   0.0000000,   0.3091071,
    0.0000000,   0.0000000;
        -0.8067323,  -0.5790557,  -0.1708886,  -0.1828683,   0.0000000,
    -0.1450338,  -0.1846675;
        -0.8067323,  -0.5790557,   0.1708886,  -0.1828683,   0.0000000,
    -0.1846675,  -0.1450338;
       ];
    assert_abs_diff_eq!(got, want, epsilon = 1e-7);
}

#[test]
fn test_do_scf() {
    struct Test {
        dir: &'static str,
        want: f64,
        eps: f64,
    }
    let tests = [
        Test {
            dir: "testfiles/h2o/STO-3G",
            want: -82.944446990003,
            eps: 1e-7,
        },
        Test {
            dir: "testfiles/h2o/DZ",
            want: -83.980246037187,
            eps: 1e-12,
        },
        Test {
            dir: "testfiles/h2o/DZP",
            want: -84.011188854711,
            eps: 1e-12,
        },
        Test {
            dir: "testfiles/ch4/STO-3G",
            want: -53.224154786383,
            eps: 1e-7,
        },
    ];
    for test in &tests[..] {
        let dir = Path::new(test.dir);
        let s = if test.dir.contains("STO-3G") {
            let mol = Molecule::load(dir.join("geom.dat"));
            let shells = Basis::load("basis_sets/sto-3g.json", &mol);
            shells.overlap_ints()
        } else if test.dir.ends_with("DZ") {
            let mol = Molecule::load(dir.join("geom.dat"));
            let shells = Basis::load("basis_sets/dz.json", &mol);
            shells.overlap_ints()
        // } else if test.dir.ends_with("DZP") {
        //     let mol = Molecule::load(dir.join("geom.dat"));
        //     let shells = Basis::load("basis_sets/dzp.json", &mol);
        //     // let want = overlap_integrals(dir.join("s.dat"));
        //     // println!("want={:.4}", want);
        //     let got = shells.overlap_ints();
        //     // println!("got={:.4}", got);
        //     got
        } else {
            overlap_integrals(dir.join("s.dat"))
        };
        // let t = shells.compute_1body_ints(Operator::Kinetic);
        // let v = shells.compute_1body_ints(Operator::nuclear(&mol));
        let t = kinetic_integrals(dir.join("t.dat"));
        let v = attraction_integrals(dir.join("v.dat"));
        let hcore = t + v;
        let s12 = build_orthog_matrix(s);
        let mol = Molecule::load(dir.join("geom.dat"));
        let eri = Eri::new(dir.join("eri.dat"));
        let got = do_scf(&hcore, &s12, &eri, mol.nelec(), 1e-12, 7e-12);
        assert_abs_diff_eq!(got, test.want, epsilon = test.eps);
    }
}
