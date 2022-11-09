use std::{
    io::{BufRead, BufReader},
    path::Path,
};

use nalgebra::SymmetricEigen;

use crate::{eri::Eri, molecule::Molecule, Dmat};

/// compute the nuclear repulsion energy for `mol`
pub fn nuclear_repulsion(mol: &Molecule) -> f64 {
    let mut ret = 0.0;
    for (i, ai) in mol.atoms.iter().enumerate() {
        for aj in mol.atoms.iter().skip(i + 1) {
            let r = (ai.coord - aj.coord).norm();
            ret += ai.atomic_number as f64 * aj.atomic_number as f64 / r;
        }
    }
    ret
}

pub fn load_sym_matrix<P: AsRef<Path>>(filename: P) -> Dmat {
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

pub fn overlap_integrals<P: AsRef<Path>>(filename: P) -> Dmat {
    load_sym_matrix(filename)
}

pub fn kinetic_integrals<P: AsRef<Path>>(filename: P) -> Dmat {
    load_sym_matrix(filename)
}

pub fn attraction_integrals<P: AsRef<Path>>(filename: P) -> Dmat {
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
    pairs.sort_by(|(_, a), (_, b)| a.partial_cmp(b).unwrap());
    let (rows, cols) = f0.shape();
    let mut c0p = Dmat::zeros(rows, cols);
    (0..cols).for_each(|i| {
        c0p.set_column(i, &vecs.column(pairs[i].0));
    });
    let c0 = s12 * c0p;
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

pub fn fock(dens: &Dmat, hcore: &Dmat, eri: &Eri) -> Dmat {
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

fn rmsd(dnew: &Dmat, dold: &Dmat) -> f64 {
    (dnew - dold).norm()
}

/// run the scf procedure and return the final energy
pub fn do_scf(
    hcore: &Dmat,
    s12: &Dmat,
    eri: &Eri,
    nelec: usize,
    d1: f64,
    d2: f64,
) -> f64 {
    // compute the initial guess density
    println!(
        "\n{:<5}{:^21}{:^21}{:^21}",
        "Iter", "E(elec)", "Delta(E)", "RMSD"
    );
    let mut e_old = 0.0;
    let mut r_old = 0.0;
    let mut f = hcore.clone();
    let mut d = density(&f, s12, nelec);
    let mut d_old = d.clone();
    let mut de = 2.0 * d1;
    let mut drmsd = 2.0 * d2;
    let mut iter = 0;
    while de.abs() > d1 || drmsd.abs() > d2 {
        let e = energy(&d, hcore, &f);
        f = fock(&d, hcore, eri);
        let r = rmsd(&d, &d_old);
        de = e - e_old;
        drmsd = r - r_old;
        println!("{:5}{:21.12}{:21.12}{:21.12}", iter, e, de, r);
        d_old = d;
        d = density(&f, s12, nelec);
        e_old = e;
        r_old = r;
        iter += 1;
    }
    e_old
}
