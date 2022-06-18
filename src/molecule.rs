use na::Vector3;
use nalgebra as na;
use std::{fs::File, io::BufRead, io::BufReader, iter::zip};

use crate::{Angle, Bond, Tors};

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
pub struct Molecule {
    /// atomic charges
    zs: Vec<i32>,

    /// coordinates
    coords: Vec<Vector3<f64>>,
}

/// return the unit vector in the direction of v
fn unit(v: Vector3<f64>) -> Vector3<f64> {
    v / v.magnitude()
}

impl Molecule {
    pub fn new(zs: Vec<i32>, coords: Vec<Vector3<f64>>) -> Self {
        Self { zs, coords }
    }

    pub fn load(geomfile: &str) -> Self {
        let f = File::open(geomfile).expect("failed to open geomfile");
        let lines = BufReader::new(f).lines().flatten().skip(1);
        let mut mol = Molecule::default();
        for line in lines {
            let vec: Vec<&str> = line.split_whitespace().collect();
            if vec.len() != 4 {
                continue;
            }
            mol.zs
                .push(vec[0].parse().expect("failed to parse atomic number"));
            let coords: Vec<f64> = vec[1..4]
                .iter()
                .map(|c| c.parse().expect("failed to parse coord"))
                .collect();
            mol.coords
                .push(na::vector![coords[0], coords[1], coords[2]]);
        }
        mol
    }

    pub fn dist(&self, i: usize, j: usize) -> f64 {
        (self.coords[i] - self.coords[j]).magnitude()
    }

    pub fn angle(&self, i: usize, j: usize, k: usize) -> f64 {
        let atom = self.coords[i];
        let btom = self.coords[j];
        let ctom = self.coords[k];
        let eji = unit(atom - btom);
        let ejk = unit(ctom - btom);
        return eji.dot(&ejk).acos();
    }

    pub fn center_of_mass(&self) -> Vector3<f64> {
        let mut ret = na::vector![0.0, 0.0, 0.0];
        let mut m: f64 = 0.0;
        for (z, c) in zip(&self.zs, &self.coords) {
            let mi = PTABLE[*z as usize];
            m += mi;
            ret += mi * c;
        }
        ret / m
    }

    pub fn translate(&mut self, vec: Vector3<f64>) {
        for c in &mut self.coords {
            *c -= vec;
        }
    }

    pub fn bond_lengths(&self) -> Vec<Bond> {
        let mol = self;
        let mut ret = Vec::new();
        let len = mol.coords.len();
        for i in 0..len {
            for j in i + 1..len {
                ret.push(Bond::new(j, i, mol.dist(i, j)));
            }
        }
        ret
    }

    /// return all of the bond angles in mol for which the legs of the angle are
    /// less than 4.0 in degrees
    pub fn bond_angles(&self) -> Vec<Angle> {
        let mol = self;
        let mut ret = Vec::new();
        let len = mol.coords.len();
        for i in 0..len {
            for j in i + 1..len {
                for k in j + 1..len {
                    if mol.dist(i, j) > 4.0 || mol.dist(j, k) > 4.0 {
                        continue;
                    }
                    ret.push(Angle::new(
                        k,
                        j,
                        i,
                        mol.angle(i, j, k).to_degrees(),
                    ));
                }
            }
        }
        ret
    }

    /// Compute the out-of-plane angles for `mol`
    pub fn oop_angles(&self) -> Vec<Tors> {
        let mol = self;
        let mut ret = Vec::new();
        let len = mol.coords.len();
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
                        let atom = mol.coords[i];
                        let btom = mol.coords[j];
                        let ctom = mol.coords[k];
                        let dtom = mol.coords[l];
                        let ekl = unit(dtom - ctom);
                        let eki = unit(atom - ctom);
                        let ekj = unit(btom - ctom);
                        let p_jkl = mol.angle(j, k, l).sin();
                        let mut tmp = (ekj.cross(&ekl) / p_jkl).dot(&eki);
                        if tmp < -1.0 {
                            tmp = -1.0;
                        } else if tmp > 1.0 {
                            tmp = 1.0;
                        }
                        ret.push(Tors::new(
                            i,
                            j,
                            k,
                            l,
                            tmp.asin().to_degrees(),
                        ));
                    }
                }
            }
        }
        ret
    }

    pub fn torsional_angles(&self) -> Vec<Tors> {
        let mol = self;
        let mut ret = Vec::new();
        let len = mol.coords.len();
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
                        let atom = mol.coords[i];
                        let btom = mol.coords[j];
                        let ctom = mol.coords[k];
                        let dtom = mol.coords[l];
                        let eij = unit(btom - atom);
                        let ejk = unit(ctom - btom);
                        let ekl = unit(dtom - ctom);
                        let p_ijk = mol.angle(i, j, k).sin();
                        let p_jkl = mol.angle(j, k, l).sin();
                        let cross1 = eij.cross(&ejk);
                        let cross2 = ejk.cross(&ekl);
                        let dot = cross1.dot(&cross2);
                        let mut div = dot / (p_ijk * p_jkl);
                        if div < -1.0 {
                            div = -1.0;
                        } else if div > 1.0 {
                            div = 1.0;
                        }
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
}

impl approx::AbsDiffEq for Molecule {
    type Epsilon = f64;

    fn default_epsilon() -> Self::Epsilon {
        1e-7
    }

    fn abs_diff_eq(&self, other: &Self, epsilon: Self::Epsilon) -> bool {
        if self.zs != other.zs {
            return false;
        }
        if self.coords.len() != other.coords.len() {
            return false;
        }
        for (i, c) in self.coords.iter().enumerate() {
            if (c - other.coords[i]).norm() > epsilon {
                return false;
            }
        }
        return true;
    }
}