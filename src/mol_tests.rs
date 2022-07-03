use std::fs::read_to_string;

use approx::{assert_abs_diff_eq, AbsDiffEq};
use na::{dvector, DVector};
use nalgebra as na;

use crate::molecule::Molecule;
use crate::{Angle, Bond, Dvec, Tors};

impl AbsDiffEq for Bond {
    type Epsilon = f64;

    fn default_epsilon() -> Self::Epsilon {
        1e-5
    }

    fn abs_diff_eq(&self, other: &Self, epsilon: Self::Epsilon) -> bool {
        self.i == other.i
            && self.j == other.j
            && (self.val - other.val).abs() < epsilon
    }
}

impl AbsDiffEq for Angle {
    type Epsilon = f64;

    fn default_epsilon() -> Self::Epsilon {
        1e-6
    }

    fn abs_diff_eq(&self, other: &Self, epsilon: Self::Epsilon) -> bool {
        self.i == other.i
            && self.j == other.j
            && self.k == other.k
            && (self.val - other.val).abs() < epsilon
    }
}

impl AbsDiffEq for Tors {
    type Epsilon = f64;

    fn default_epsilon() -> Self::Epsilon {
        1e-6
    }

    fn abs_diff_eq(&self, other: &Self, epsilon: Self::Epsilon) -> bool {
        self.i == other.i
            && self.j == other.j
            && self.k == other.k
            && self.l == other.l
            && (self.val - other.val).abs() < epsilon
    }
}

#[test]
fn test_read_geom() {
    let got = Molecule::load("testfiles/acetaldehyde.dat");
    let want = Molecule::from_vecs(
        vec![6, 6, 8, 1, 1, 1, 1],
        vec![
            na::vector![0.000000000000, 0.000000000000, 0.000000000000],
            na::vector![0.000000000000, 0.000000000000, 2.845112131228],
            na::vector![1.899115961744, 0.000000000000, 4.139062527233],
            na::vector![-1.894048308506, 0.000000000000, 3.747688672216],
            na::vector![1.942500819960, 0.000000000000, -0.701145981971],
            na::vector![-1.007295466862, -1.669971842687, -0.705916966833],
            na::vector![-1.007295466862, 1.669971842687, -0.705916966833],
        ],
    );
    assert_eq!(got, want);
}

fn load_bonds(filename: &str) -> na::DVector<Bond> {
    let mut ret = Vec::new();
    let data = read_to_string(filename).unwrap();
    for line in data.lines() {
        let split: Vec<_> = line.trim().split_whitespace().collect();
        ret.push(Bond::new(
            split[0].parse().unwrap(),
            split[1].parse().unwrap(),
            split[2].parse().unwrap(),
        ))
    }
    DVector::from(ret)
}

#[test]
fn test_bond_lengths() {
    struct Test {
        infile: &'static str,
        want: na::DVector<Bond>,
    }
    let tests = vec![
        Test {
            infile: "testfiles/acetaldehyde.dat",
            want: load_bonds("testfiles/acetaldehyde_bonds.dat"),
        },
        Test {
            infile: "testfiles/benzene.dat",
            want: load_bonds("testfiles/benzene_bonds.dat"),
        },
        Test {
            infile: "testfiles/allene.dat",
            want: load_bonds("testfiles/allene_bonds.dat"),
        },
    ];
    for test in tests {
        let got = Molecule::load(test.infile).bond_lengths();
        assert_abs_diff_eq!(na::DVector::from_row_slice(&got), test.want);
    }
}

#[test]
fn test_bond_angles() {
    let got = Molecule::load("testfiles/acetaldehyde.dat").bond_angles();
    let want = vec![
        Angle::new(2, 1, 0, 124.268308),
        Angle::new(3, 1, 0, 115.479341),
        Angle::new(5, 4, 0, 35.109529),
        Angle::new(6, 4, 0, 35.109529),
        Angle::new(6, 5, 0, 36.373677),
        Angle::new(3, 2, 1, 28.377448),
        Angle::new(6, 5, 4, 60.484476),
    ];
    assert_abs_diff_eq!(
        na::DVector::from_row_slice(&got),
        na::DVector::from_row_slice(&want)
    );
}

#[test]
fn test_op_angles() {
    let got = Molecule::load("testfiles/acetaldehyde.dat").oop_angles();
    let want = vec![
        Tors::new(0, 3, 1, 2, -0.000000),
        Tors::new(0, 6, 4, 5, 19.939726),
        Tors::new(0, 6, 5, 4, -19.850523),
        Tors::new(0, 5, 6, 4, 19.850523),
        Tors::new(1, 5, 0, 4, 53.678778),
        Tors::new(1, 6, 0, 4, -53.678778),
        Tors::new(1, 6, 0, 5, 54.977064),
        Tors::new(2, 3, 1, 0, 0.000000),
        Tors::new(3, 2, 1, 0, -0.000000),
        Tors::new(4, 5, 0, 1, -53.651534),
        Tors::new(4, 6, 0, 1, 53.651534),
        Tors::new(4, 6, 0, 5, -54.869992),
        Tors::new(4, 6, 5, 0, 29.885677),
        Tors::new(4, 5, 6, 0, -29.885677),
        Tors::new(5, 4, 0, 1, 53.626323),
        Tors::new(5, 6, 0, 1, -56.277112),
        Tors::new(5, 6, 0, 4, 56.194621),
        Tors::new(5, 6, 4, 0, -30.558964),
        Tors::new(5, 4, 6, 0, 31.064344),
        Tors::new(6, 4, 0, 1, -53.626323),
        Tors::new(6, 5, 0, 1, 56.277112),
        Tors::new(6, 5, 0, 4, -56.194621),
        Tors::new(6, 5, 4, 0, 30.558964),
        Tors::new(6, 4, 5, 0, -31.064344),
    ];
    assert_abs_diff_eq!(
        na::DVector::from_row_slice(&got),
        na::DVector::from_row_slice(&want)
    );
}

#[test]
fn test_torsional_angles() {
    let got = Molecule::load("testfiles/acetaldehyde.dat").torsional_angles();
    let want = vec![
        Tors::new(3, 2, 1, 0, 180.000000),
        Tors::new(6, 5, 4, 0, 36.366799),
    ];
    assert_abs_diff_eq!(
        na::DVector::from_row_slice(&got),
        na::DVector::from_row_slice(&want)
    );
}

#[test]
fn test_com() {
    let mol = Molecule::load("testfiles/acetaldehyde.dat");
    let got = mol.center_of_mass();
    let want = na::vector![0.64494926, 0.00000000, 2.31663792];
    assert_abs_diff_eq!(got, want, epsilon = 2e-8);
}

#[test]
fn test_translate() {
    let mut mol = Molecule::load("testfiles/acetaldehyde.dat");
    let com = mol.center_of_mass();
    mol.translate(com);
    let want = Molecule::from_vecs(
        vec![6, 6, 8, 1, 1, 1, 1],
        vec![
            na::vector![-0.644949260000, 0.000000000000, -2.316637920000],
            na::vector![-0.644949260000, 0.000000000000, 0.528474211228],
            na::vector![1.254166701744, 0.000000000000, 1.822424607233],
            na::vector![-2.538997568506, 0.000000000000, 1.431050752216],
            na::vector![1.297551559960, 0.000000000000, -3.017783901971],
            na::vector![-1.652244726862, -1.669971842687, -3.022554886833],
            na::vector![-1.652244726862, 1.669971842687, -3.022554886833],
        ],
    );
    assert_abs_diff_eq!(mol, want);
}

#[test]
fn test_moi() {
    let mut mol = Molecule::load("testfiles/acetaldehyde.dat");
    let com = mol.center_of_mass();
    mol.translate(com);
    let got = mol.moi();
    #[rustfmt::skip]
    let want = na::Matrix3::new(
        156.154091561645, 0.000000000000, -52.855584120568,
        0.000000000000, 199.371126996236, 0.000000000000,
        -52.855584120568, 0.000000000000, 54.459548882464,
    );
    assert_abs_diff_eq!(got, want, epsilon = 8e-7);
}

#[test]
fn test_rots() {
    struct Test {
        infile: &'static str,
        want: Dvec,
        eps: f64,
    }
    let tests = vec![
        Test {
            infile: "testfiles/acetaldehyde.dat",
            want: dvector![56461.542, 10102.130, 9052.169],
            eps: 1e-3,
        },
        Test {
            infile: "testfiles/benzene.dat",
            want: dvector![5791.637, 5791.636, 2895.818],
            eps: 1e-3,
        },
        Test {
            infile: "testfiles/allene.dat",
            want: dvector![167151.728, 8558.659, 8558.659],
            eps: 2e-3,
        },
    ];
    for test in tests {
        let mut mol = Molecule::load(test.infile);
        let com = mol.center_of_mass();
        mol.translate(com);
        let moi = mol.moi();
        let got = mol.rots(&moi);
        assert_abs_diff_eq!(Dvec::from(got), test.want, epsilon = test.eps);
    }
}
