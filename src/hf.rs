use std::{
    fs::read_to_string,
    io::{BufRead, BufReader},
};

use crate::Dmat;

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

#[cfg(test)]
mod tests {
    use approx::assert_abs_diff_eq;
    use nalgebra::dmatrix;

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
}
