use std::fmt::Display;

use nalgebra as na;

#[cfg(test)]
mod mol_tests;

pub mod hf;
pub mod molecule;

type Vec3 = na::Vector3<f64>;
type Mat3 = na::Matrix3<f64>;
type Dvec = na::DVector<f64>;
type Dmat = na::DMatrix<f64>;

#[derive(Clone, Debug, PartialEq)]
pub struct Bond {
    i: usize,
    j: usize,
    val: f64,
}

impl Display for Bond {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        write!(f, "{:5}{:5}{:8.5}", self.i, self.j, self.val)
    }
}

impl Bond {
    pub fn new(i: usize, j: usize, val: f64) -> Self {
        Self { i, j, val }
    }
}

#[derive(Clone, Debug, PartialEq)]
pub struct Angle {
    i: usize,
    j: usize,
    k: usize,
    val: f64,
}

impl Display for Angle {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        write!(f, "{:5}{:5}{:5}{:11.6}", self.i, self.j, self.k, self.val)
    }
}

impl Angle {
    pub fn new(i: usize, j: usize, k: usize, val: f64) -> Self {
        Self { i, j, k, val }
    }
}

#[derive(Clone, Debug, PartialEq)]
pub struct Tors {
    i: usize,
    j: usize,
    k: usize,
    l: usize,
    val: f64,
}

impl Display for Tors {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        write!(
            f,
            "{:5}{:5}{:5}{:5}{:11.6}",
            self.i, self.j, self.k, self.l, self.val
        )
    }
}

impl Tors {
    pub fn new(i: usize, j: usize, k: usize, l: usize, val: f64) -> Tors {
        Tors { i, j, k, l, val }
    }
}
