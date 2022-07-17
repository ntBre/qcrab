use std::{
    fs::read_to_string,
    io::{BufRead, BufReader},
};

use nalgebra::SymmetricEigen;

use crate::{eri::Eri, Dmat, Mat3};

// TODO eventually take a Molecule, not a filename
pub fn nuclear_repulsion(filename: &str) -> f64 {
    let data = read_to_string(filename).unwrap();
    data.trim().parse().unwrap()
}

pub fn load_sym_matrix(filename: &str) -> Dmat {
    let f = std::fs::File::open(filename).unwrap();
    let f = BufReader::new(f).lines();
    let mut rows = 1;
    let mut cols = 1;
    let mut ret = Dmat::zeros(rows, cols);
    for line in f.flatten() {
        let split: Vec<_> = line.split_whitespace().collect();
        if split.len() != 3 {
            continue;
        }
        let i: usize = split[0].parse().unwrap();
        let j: usize = split[1].parse().unwrap();
        let v: f64 = split[2].parse().unwrap();
        while i > rows {
            ret = ret.insert_row(rows, 0.0);
            rows += 1;
        }
        while j > cols {
            ret = ret.insert_column(cols, 0.0);
            cols += 1;
        }
        ret[(i - 1, j - 1)] = v;
    }
    ret.fill_upper_triangle_with_lower_triangle();
    ret
}

pub fn overlap_integrals(filename: &str) -> Dmat {
    load_sym_matrix(filename)
}

pub fn kinetic_integrals(filename: &str) -> Dmat {
    load_sym_matrix(filename)
}

pub fn attraction_integrals(filename: &str) -> Dmat {
    load_sym_matrix(filename)
}

pub fn build_orthog_matrix(s: Dmat) -> Dmat {
    let (r, c) = s.shape();
    let sym = SymmetricEigen::new(s);
    let ls = sym.eigenvectors;
    let mut lambda = Dmat::zeros(r, c);
    assert_eq!(r, c);
    for i in 0..r {
        lambda[(i, i)] = 1.0 / sym.eigenvalues[i].sqrt();
    }
    ls.clone() * lambda * ls.transpose()
}

pub fn density(fock: &Dmat, s12: &Dmat, nelec: usize) -> Dmat {
    let f0 = s12.transpose() * fock * s12;
    let sym = SymmetricEigen::new(f0.clone());
    let vecs = sym.eigenvectors;
    let vals = sym.eigenvalues;
    // this is wrong for open shells
    let nocc = nelec / 2;
    // sort the eigenvalues and then the eigenvectors correspondingly
    let mut pairs: Vec<_> = vals.iter().enumerate().collect();
    pairs.sort_by(|(_, a), (_, b)| a.partial_cmp(&b).unwrap());
    let (rows, cols) = f0.shape();
    let mut c0p = Dmat::zeros(rows, cols);
    for i in 0..cols {
        c0p.set_column(i, &vecs.column(pairs[i].0));
    }
    let c0 = s12 * c0p.clone();
    let mut ret = Dmat::zeros(rows, cols);
    for mu in 0..rows {
        for nu in 0..cols {
            for m in 0..nocc {
                ret[(mu, nu)] += c0[(mu, m)] * c0[(nu, m)];
            }
        }
    }
    ret
}

pub fn energy(dens: &Dmat, hcore: &Dmat, fock: &Dmat) -> f64 {
    let mut energy = 0.0;
    let (r, c) = dens.shape();
    for i in 0..r {
        for j in 0..c {
            energy += dens[(i, j)] * (hcore[(i, j)] + fock[(i, j)]);
        }
    }
    energy
}

fn fock(dens: &Dmat, hcore: &Dmat, eri: &Eri) -> Dmat {
    let (r, c) = hcore.shape();
    assert_eq!(r, c);
    let mut fock = hcore.clone();
    for m in 0..r {
        for n in 0..r {
            for l in 0..r {
                for s in 0..r {
                    fock[(m, n)] += dens[(l, s)]
                        * (2.0 * eri[(m, n, l, s)] - eri[(m, l, n, s)]);
                }
            }
        }
    }
    fock
}

#[cfg(test)]
mod tests {
    use approx::{abs_diff_eq, assert_abs_diff_eq};
    use nalgebra::dmatrix;

    use crate::{eri::Eri, molecule::Molecule};

    use super::*;

    #[test]
    fn test_enuc() {
        let got = nuclear_repulsion("testfiles/h2o/STO-3G/enuc.dat");
        let want = 8.002367061810450;
        assert_eq!(got, want);
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
    fn test_new_fock() {
        let _eri = Eri::new("testfiles/h2o/STO-3G/eri.dat");
    }

    #[test]
    fn test_orthog() {
        let got = build_orthog_matrix(overlap_integrals(
            "testfiles/h2o/STO-3G/s.dat",
        ));
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
        let s12 = build_orthog_matrix(overlap_integrals(
            "testfiles/h2o/STO-3G/s.dat",
        ));
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
        let s12 = build_orthog_matrix(overlap_integrals(
            "testfiles/h2o/STO-3G/s.dat",
        ));
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
        let s12 = build_orthog_matrix(overlap_integrals(
            "testfiles/h2o/STO-3G/s.dat",
        ));
        let mol = Molecule::load("testfiles/h2o/STO-3G/geom.dat");
        let d = density(&hcore, &s12, mol.nelec());
        let eri = Eri::new("testfiles/h2o/STO-3G/eri.dat");
        let got = fock(&d, &hcore, &eri);
        let want = dmatrix![
        -18.8132695,  -4.8726875,  -0.0000000,  -0.0115290,   0.0000000,  -0.8067323,  -0.8067323;
         -4.8726875,  -1.7909029,  -0.0000000,  -0.1808692,   0.0000000,  -0.5790557,  -0.5790557;
         -0.0000000,  -0.0000000,   0.1939644,   0.0000000,   0.0000000,  -0.1708886,   0.1708886;
         -0.0115290,  -0.1808692,   0.0000000,   0.2391247,   0.0000000,  -0.1828683,  -0.1828683;
          0.0000000,   0.0000000,   0.0000000,   0.0000000,   0.3091071,   0.0000000,   0.0000000;
         -0.8067323,  -0.5790557,  -0.1708886,  -0.1828683,   0.0000000,  -0.1450338,  -0.1846675;
         -0.8067323,  -0.5790557,   0.1708886,  -0.1828683,   0.0000000,  -0.1846675,  -0.1450338;
           ];
        assert_abs_diff_eq!(got, want, epsilon = 1e-7);
    }
}
