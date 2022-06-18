use approx::{assert_abs_diff_eq, AbsDiffEq};
use nalgebra as na;

use crate::molecule::Molecule;
use crate::{Angle, Bond, Tors};

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
    let got = Molecule::load("inp/geom.xyz");
    let want = Molecule::new(
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

#[test]
fn test_bond_lengths() {
    let got = Molecule::load("inp/geom.xyz").bond_lengths();
    let want = na::dvector![
        Bond::new(1, 0, 2.84511),
        Bond::new(2, 0, 4.55395),
        Bond::new(3, 0, 4.19912),
        Bond::new(4, 0, 2.06517),
        Bond::new(5, 0, 2.07407),
        Bond::new(6, 0, 2.07407),
        Bond::new(2, 1, 2.29803),
        Bond::new(3, 1, 2.09811),
        Bond::new(4, 1, 4.04342),
        Bond::new(5, 1, 4.05133),
        Bond::new(6, 1, 4.05133),
        Bond::new(3, 2, 3.81330),
        Bond::new(4, 2, 4.84040),
        Bond::new(5, 2, 5.89151),
        Bond::new(6, 2, 5.89151),
        Bond::new(4, 3, 5.87463),
        Bond::new(5, 3, 4.83836),
        Bond::new(6, 3, 4.83836),
        Bond::new(5, 4, 3.38971),
        Bond::new(6, 4, 3.38971),
        Bond::new(6, 5, 3.33994)
    ];
    assert_abs_diff_eq!(na::DVector::from_row_slice(&got), want);
}

#[test]
fn test_bond_angles() {
    let got = Molecule::load("inp/geom.xyz").bond_angles();
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
    let got = Molecule::load("inp/geom.xyz").oop_angles();
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
    let got = Molecule::load("inp/geom.xyz").torsional_angles();
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
    let mol = Molecule::load("inp/geom.xyz");
    let got = mol.center_of_mass();
    let want = na::vector![0.64494926, 0.00000000, 2.31663792];
    assert_abs_diff_eq!(got, want, epsilon = 2e-8);
}

#[test]
fn test_translate() {
    let mut mol = Molecule::load("inp/geom.xyz");
    let com = mol.center_of_mass();
    mol.translate(com);
    let want = Molecule::new(
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
