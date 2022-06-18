use nalgebra as na;

#[cfg(test)]
mod mol_tests;

pub mod molecule;

type Vec3 = na::Vector3<f64>;
type Mat3 = na::Matrix3<f64>;

#[derive(Clone, Debug, PartialEq)]
pub struct Bond {
    i: usize,
    j: usize,
    val: f64,
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

impl Tors {
    pub fn new(i: usize, j: usize, k: usize, l: usize, val: f64) -> Tors {
        Tors { i, j, k, l, val }
    }
}
