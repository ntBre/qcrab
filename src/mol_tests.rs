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

fn load_angles(filename: &str) -> na::DVector<Angle> {
    let mut ret = Vec::new();
    let data = read_to_string(filename).unwrap();
    for line in data.lines() {
        let split: Vec<_> = line.trim().split_whitespace().collect();
        ret.push(Angle::new(
            split[0].parse().unwrap(),
            split[1].parse().unwrap(),
            split[2].parse().unwrap(),
            split[3].parse().unwrap(),
        ))
    }
    DVector::from(ret)
}

#[test]
fn test_bond_angles() {
    struct Test {
        infile: &'static str,
        want: na::DVector<Angle>,
    }

    impl Test {
        fn new(infile: &'static str, want: na::DVector<Angle>) -> Self {
            Self { infile, want }
        }
    }
    let tests = [
        Test::new(
            "testfiles/acetaldehyde.dat",
            load_angles("testfiles/acetaldehyde_angles.dat"),
        ),
        Test::new(
            "testfiles/benzene.dat",
            load_angles("testfiles/benzene_angles.dat"),
        ),
        Test::new(
            "testfiles/allene.dat",
            load_angles("testfiles/allene_angles.dat"),
        ),
    ];
    for test in tests {
        let got = Molecule::load(test.infile).bond_angles();
        assert_eq!(got.len(), test.want.len());
        assert_abs_diff_eq!(
            na::DVector::from_row_slice(&got),
            test.want,
            epsilon = 1e-6
        );
    }
}

fn load_oop(filename: &str) -> na::DVector<Tors> {
    let mut ret = Vec::new();
    let data = read_to_string(filename).unwrap();
    for line in data.lines() {
        let split: Vec<_> = line.trim().split_whitespace().collect();
        ret.push(Tors::new(
            split[0].parse().unwrap(),
            split[1].parse().unwrap(),
            split[2].parse().unwrap(),
            split[3].parse().unwrap(),
            split[4].parse().unwrap(),
        ))
    }
    DVector::from(ret)
}

struct OopTest {
    infile: &'static str,
    want: DVector<Tors>,
}

impl OopTest {
    fn new(infile: &'static str, want: DVector<Tors>) -> Self {
        Self { infile, want }
    }
}

#[test]
fn test_op_angles() {
    let tests = [
        OopTest::new(
            "testfiles/acetaldehyde.dat",
            load_oop("testfiles/acetaldehyde_oop.dat"),
        ),
        OopTest::new(
            "testfiles/benzene.dat",
            load_oop("testfiles/benzene_oop.dat"),
        ),
        OopTest::new(
            "testfiles/allene.dat",
            load_oop("testfiles/allene_oop.dat"),
        ),
    ];
    for test in tests {
        let got = Molecule::load(test.infile).oop_angles();
        assert_eq!(got.len(), test.want.len());
        assert_abs_diff_eq!(na::DVector::from_row_slice(&got), test.want);
    }
}

#[test]
fn test_torsional_angles() {
    let tests = [
        OopTest::new(
            "testfiles/acetaldehyde.dat",
            load_oop("testfiles/acetaldehyde_tors.dat"),
        ),
        OopTest::new(
            "testfiles/benzene.dat",
            load_oop("testfiles/benzene_tors.dat"),
        ),
        OopTest::new(
            "testfiles/allene.dat",
            load_oop("testfiles/allene_tors.dat"),
        ),
    ];
    for test in tests {
        let got = na::DVector::from_row_slice(
            &Molecule::load(test.infile).torsional_angles(),
        );
        assert_abs_diff_eq!(got, test.want, epsilon = 3e-6);
    }
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
