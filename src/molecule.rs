use na::Vector3;
use nalgebra as na;
use std::{
    f64::consts::PI, fs::File, io::BufRead, io::BufReader, iter::zip,
    path::Path,
};

use crate::{Angle, Bond, Mat3, Tors, Vec3};

const H: f64 = 6.626_070_15e-34; // J * s
const C: f64 = 299_792_458.; // m / s
const AMU: f64 = 1.660_539_066_60e-27; // kg
const BOHR: f64 = 5.291_772_109_03e-11; // m
const MOI: f64 = H / (8.0 * PI * PI * C);

const PTABLE: [f64; 9] = [
    0.0,
    1.007_825_032,
    4.002_603_254,
    7.016_003_436,
    9.012_183_0,
    11.009_305,
    12.00,
    14.003_074_004,
    15.994_914_619,
];

#[derive(Debug, Default, PartialEq)]
pub struct Atom {
    pub(crate) atomic_number: usize,
    pub(crate) coord: Vec3,
}

impl Atom {
    pub fn mass(&self) -> f64 {
        PTABLE[self.atomic_number]
    }
}

#[derive(Debug, Default, PartialEq)]
pub struct Molecule {
    pub(crate) atoms: Vec<Atom>,
}

/// return the unit vector in the direction of v
fn unit(v: Vector3<f64>) -> Vector3<f64> {
    v / v.magnitude()
}

impl Molecule {
    pub fn new(atoms: Vec<Atom>) -> Self {
        Self { atoms }
    }

    pub fn nelec(&self) -> usize {
        self.atoms.iter().map(|a| a.atomic_number).sum()
    }

    pub fn from_vecs(zs: Vec<usize>, coords: Vec<Vec3>) -> Self {
        let mut ret = Self::default();
        assert_eq!(zs.len(), coords.len());
        for (atomic_number, coord) in zip(zs, coords) {
            ret.atoms.push(Atom {
                atomic_number,
                coord,
            })
        }
        ret
    }

    pub fn load<P: AsRef<Path>>(geomfile: P) -> Self {
        let f = File::open(geomfile).expect("failed to open geomfile");
        let lines = BufReader::new(f).lines().flatten().skip(1);
        let mut mol = Molecule::default();
        for line in lines {
            let vec: Vec<&str> = line.split_whitespace().collect();
            if vec.len() != 4 {
                continue;
            }
            let coords: Vec<f64> = vec[1..4]
                .iter()
                .map(|c| c.parse().expect("failed to parse coord"))
                .collect();
            mol.atoms.push(Atom {
                atomic_number: vec[0]
                    .parse()
                    .expect("failed to parse atomic number"),
                coord: na::vector![coords[0], coords[1], coords[2]],
            });
        }
        mol
    }

    pub fn dist(&self, i: usize, j: usize) -> f64 {
        (self.atoms[i].coord - self.atoms[j].coord).magnitude()
    }

    /// compute the angle between atoms `i`, `j`, and `k`, assuming that `j` is
    /// the central atom. panics if any of these is out of range
    pub fn angle(&self, i: usize, j: usize, k: usize) -> f64 {
        let atom = self.atoms[i].coord;
        let btom = self.atoms[j].coord;
        let ctom = self.atoms[k].coord;
        let eji = unit(atom - btom);
        let ejk = unit(ctom - btom);
        eji.dot(&ejk).acos()
    }

    pub fn center_of_mass(&self) -> Vector3<f64> {
        let mut ret = na::vector![0.0, 0.0, 0.0];
        let mut m: f64 = 0.0;
        for atom in &self.atoms {
            let mi = atom.mass();
            m += mi;
            ret += mi * atom.coord;
        }
        ret / m
    }

    pub fn translate(&mut self, vec: Vector3<f64>) {
        for atom in &mut self.atoms {
            atom.coord -= vec;
        }
    }

    pub fn bond_lengths(&self) -> Vec<Bond> {
        let mol = self;
        let mut ret = Vec::new();
        let len = mol.atoms.len();
        for i in 0..len {
            for j in 0..i {
                ret.push(Bond::new(i, j, mol.dist(i, j)));
            }
        }
        ret
    }

    /// return all of the bond angles in mol for which the legs of the angle are
    /// less than 4.0 in degrees
    pub fn bond_angles(&self) -> Vec<Angle> {
        let mut ret = Vec::new();
        let len = self.atoms.len();
        for i in 0..len {
            for j in 0..i {
                for k in 0..j {
                    if self.dist(i, j) < 4.0 && self.dist(j, k) < 4.0 {
                        ret.push(Angle::new(
                            k,
                            j,
                            i,
                            self.angle(i, j, k).to_degrees(),
                        ));
                    }
                }
            }
        }
        ret
    }

    pub fn oop(&self, i: usize, j: usize, k: usize, l: usize) -> Tors {
        let atom = self.atoms[i].coord;
        let btom = self.atoms[j].coord;
        let ctom = self.atoms[k].coord;
        let dtom = self.atoms[l].coord;
        let ekj = unit(btom - ctom);
        let ekl = unit(dtom - ctom);
        let eki = unit(atom - ctom);
        let pjkl = self.angle(j, k, l).sin();
        let mut tmp = (ekj.cross(&ekl) / pjkl).dot(&eki);
        tmp = tmp.clamp(-1.0, 1.0);
        Tors::new(i, j, k, l, tmp.asin().to_degrees())
    }

    /// Compute the out-of-plane angles for `mol`
    pub fn oop_angles(&self) -> Vec<Tors> {
        let mol = self;
        let mut ret = Vec::new();
        let len = mol.atoms.len();
        for i in 0..len {
            for k in 0..len {
                for j in 0..len {
                    for l in 0..j {
                        if i == j || i == k || i == l || j == k || k == l {
                            continue;
                        }
                        if mol.dist(i, k) >= 4.0
                            || mol.dist(j, k) >= 4.0
                            || mol.dist(k, l) >= 4.0
                        {
                            continue;
                        }
                        ret.push(self.oop(i, j, k, l));
                    }
                }
            }
        }
        ret
    }

    pub fn torsional_angles(&self) -> Vec<Tors> {
        let mol = self;
        let mut ret = Vec::new();
        let len = mol.atoms.len();
        for i in 0..len {
            for j in 0..i {
                for k in 0..j {
                    for l in 0..k {
                        if mol.dist(i, j) >= 4.0
                            || mol.dist(j, k) >= 4.0
                            || mol.dist(k, l) >= 4.0
                        {
                            continue;
                        }
                        let atom = mol.atoms[i].coord;
                        let btom = mol.atoms[j].coord;
                        let ctom = mol.atoms[k].coord;
                        let dtom = mol.atoms[l].coord;
                        let eij = unit(btom - atom);
                        let ejk = unit(ctom - btom);
                        let ekl = unit(dtom - ctom);
                        let p_ijk = mol.angle(i, j, k).sin();
                        let p_jkl = mol.angle(j, k, l).sin();
                        let cross1 = eij.cross(&ejk);
                        let cross2 = ejk.cross(&ekl);
                        let dot = cross1.dot(&cross2);
                        let mut div = dot / (p_ijk * p_jkl);
                        div = div.clamp(-1.0, 1.0);
                        ret.push(Tors::new(
                            i,
                            j,
                            k,
                            l,
                            (div.acos()).to_degrees(),
                        ));
                    }
                }
            }
        }
        ret
    }

    pub fn moi(&self) -> Mat3 {
        let mut i = Mat3::zeros();
        const X: usize = 0;
        const Y: usize = 1;
        const Z: usize = 2;
        for atom in &self.atoms {
            let mi = atom.mass();
            // diagonal
            i[(X, X)] += mi * (atom.coord.y.powi(2) + atom.coord.z.powi(2));
            i[(Y, Y)] += mi * (atom.coord.x.powi(2) + atom.coord.z.powi(2));
            i[(Z, Z)] += mi * (atom.coord.x.powi(2) + atom.coord.y.powi(2));
            // off-diagonal
            i[(X, Y)] -= mi * atom.coord.x * atom.coord.y;
            i[(X, Z)] -= mi * atom.coord.x * atom.coord.z;
            i[(Y, Z)] -= mi * atom.coord.y * atom.coord.z;
        }
        i[(Y, X)] = i[(X, Y)];
        i[(Z, X)] = i[(X, Z)];
        i[(Z, Y)] = i[(Y, Z)];
        i
    }

    /// compute the rotational constants (in MHz) for `self` by diagonlizing the
    /// moment of inertia matrix in `moi`
    pub fn rots(&self, moi: &Mat3) -> Vec<f64> {
        let mut evals = moi.eigenvalues().expect("failed to diagonalize moi");
        let evals = evals.as_mut_slice();
        evals.sort_by(|a, b| a.partial_cmp(b).unwrap());
        evals
            .iter()
            .map(|e| 1.0e-6 * C * MOI / (e * AMU * BOHR * BOHR))
            .collect()
    }

    pub(crate) fn len(&self) -> usize {
        self.atoms.len()
    }
}

impl approx::AbsDiffEq for Atom {
    type Epsilon = f64;

    fn default_epsilon() -> Self::Epsilon {
        1e-7
    }

    fn abs_diff_eq(&self, other: &Self, epsilon: Self::Epsilon) -> bool {
        if self.atomic_number != other.atomic_number {
            return false;
        }
        if (self.coord - other.coord).norm() > epsilon {
            return false;
        }
        true
    }
}

impl approx::AbsDiffEq for Molecule {
    type Epsilon = f64;

    fn default_epsilon() -> Self::Epsilon {
        1e-7
    }

    fn abs_diff_eq(&self, other: &Self, epsilon: Self::Epsilon) -> bool {
        self.atoms.len() == other.atoms.len()
            && zip(&self.atoms, &other.atoms)
                .all(|(a, b)| a.abs_diff_eq(b, epsilon))
    }
}
