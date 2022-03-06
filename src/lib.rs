use na::Vector3;
use nalgebra as na;
use std::{fs::File, io::BufRead, io::BufReader};

#[derive(Debug)]
pub struct Molecule {
    // atomic charges
    zs: Vec<i32>,
    // coordinates
    coords: Vec<Vector3<f64>>,
}

impl Molecule {
    pub fn new() -> Molecule {
        Molecule {
            zs: Vec::new(),
            coords: Vec::new(),
        }
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
}

pub fn read_geom(geomfile: &str) -> Molecule {
    let f = File::open(geomfile).expect("failed to open geomfile");
    let mut lines = BufReader::new(f).lines();
    lines.next(); // skip the atom count
    let mut mol = Molecule::new();
    for line in lines.map(|x| x.unwrap()) {
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

#[derive(Debug)]
pub struct Bond {
    i: usize,
    j: usize,
    val: f64,
}

impl Bond {
    fn new(i: usize, j: usize, val: f64) -> Self {
        Self { i, j, val }
    }
}

impl PartialEq for Bond {
    fn eq(&self, other: &Self) -> bool {
        self.i == other.i && self.j == other.j && (self.val - other.val).abs() < 1e-5
    }
}

pub fn bond_lengths(mol: &Molecule) -> Vec<Bond> {
    let mut ret = Vec::new();
    let len = mol.coords.len();
    for i in 0..len {
        for j in i + 1..len {
            ret.push(Bond::new(j, i, mol.dist(i, j)));
        }
    }
    ret
}

/// return the unit vector in the direction of v
fn unit(v: Vector3<f64>) -> Vector3<f64> {
    v / v.magnitude()
}

#[derive(Debug)]
pub struct Angle {
    i: usize,
    j: usize,
    k: usize,
    val: f64,
}

impl Angle {
    fn new(i: usize, j: usize, k: usize, val: f64) -> Self {
        Self { i, j, k, val }
    }
}

impl PartialEq for Angle {
    fn eq(&self, other: &Self) -> bool {
        self.i == other.i
            && self.j == other.j
            && self.k == other.k
            && (self.val - other.val).abs() < 1e-5
    }
}

/// return all of the bond angles in mol for which the legs of the angle are
/// less than 4.0 in degrees
pub fn bond_angles(mol: &Molecule) -> Vec<Angle> {
    let mut ret = Vec::new();
    let len = mol.coords.len();
    for i in 0..len {
        for j in i + 1..len {
            for k in j + 1..len {
                if mol.dist(i, j) > 4.0 || mol.dist(j, k) > 4.0 {
                    continue;
                }
                ret.push(Angle::new(k, j, i, mol.angle(i, j, k).to_degrees()));
            }
        }
    }
    ret
}

/// Compute the out-of-plane angles for `mol`
pub fn oop_angles(mol: &Molecule) -> Vec<Tors> {
    let mut ret = Vec::new();
    let len = mol.coords.len();
    for i in 0..len {
        for k in 0..len {
            for j in 0..len {
                for l in 0..j {
                    if i == j || i == k || i == l || j == k || k == l {
                        continue;
                    }
                    if mol.dist(i, k) >= 4.0 || mol.dist(j, k) >= 4.0 || mol.dist(k, l) >= 4.0 {
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
                    ret.push(Tors::new(i, j, k, l, tmp.asin().to_degrees()));
                }
            }
        }
    }
    ret
}

#[derive(Debug)]
pub struct Tors {
    i: usize,
    j: usize,
    k: usize,
    l: usize,
    val: f64,
}

impl Tors {
    fn new(i: usize, j: usize, k: usize, l: usize, val: f64) -> Tors {
        Tors { i, j, k, l, val }
    }
}

impl PartialEq for Tors {
    fn eq(&self, other: &Self) -> bool {
        self.i == other.i
            && self.j == other.j
            && self.k == other.k
            && self.l == other.l
            && (self.val - other.val).abs() < 1e-6
    }
}

pub fn torsional_angles(mol: &Molecule) -> Vec<Tors> {
    let mut ret = Vec::new();
    let len = mol.coords.len();
    for i in 0..len {
        for j in 0..i {
            for k in 0..j {
                for l in 0..k {
                    if mol.dist(i, j) >= 4.0 || mol.dist(j, k) >= 4.0 || mol.dist(k, l) >= 4.0 {
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
                    ret.push(Tors::new(i, j, k, l, (div.acos()).to_degrees()));
                }
            }
        }
    }
    ret
}

#[cfg(test)]
mod tests {
    use super::*;
    #[test]
    fn test_read_geom() {
        let got = read_geom("inp/geom.xyz");
        let want = Molecule {
            zs: vec![6, 6, 8, 1, 1, 1, 1],
            coords: vec![
                na::vector![0.000000000000, 0.000000000000, 0.000000000000],
                na::vector![0.000000000000, 0.000000000000, 2.845112131228],
                na::vector![1.899115961744, 0.000000000000, 4.139062527233],
                na::vector![-1.894048308506, 0.000000000000, 3.747688672216],
                na::vector![1.942500819960, 0.000000000000, -0.701145981971],
                na::vector![-1.007295466862, -1.669971842687, -0.705916966833],
                na::vector![-1.007295466862, 1.669971842687, -0.705916966833],
            ],
        };
        assert_eq!(got.zs, want.zs, "got {:?}, wanted {:?}", got.zs, want.zs);
        assert_eq!(got.coords, want.coords);
    }

    #[test]
    fn test_bond_lengths() {
        let got = bond_lengths(&read_geom("inp/geom.xyz"));
        let want = vec![
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
            Bond::new(6, 5, 3.33994),
        ];
        assert_eq!(got, want);
    }

    #[test]
    fn test_bond_angles() {
        let got = bond_angles(&read_geom("inp/geom.xyz"));
        let want = vec![
            Angle::new(2, 1, 0, 124.268308),
            Angle::new(3, 1, 0, 115.479341),
            Angle::new(5, 4, 0, 35.109529),
            Angle::new(6, 4, 0, 35.109529),
            Angle::new(6, 5, 0, 36.373677),
            Angle::new(3, 2, 1, 28.377448),
            Angle::new(6, 5, 4, 60.484476),
        ];
        assert_eq!(got, want);
    }

    #[test]
    fn test_op_angles() {
        let got = oop_angles(&read_geom("inp/geom.xyz"));
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
        assert_eq!(got, want);
    }

    #[test]
    fn test_torsional_angles() {
        let got = torsional_angles(&read_geom("inp/geom.xyz"));
        let want = vec![
            Tors::new(3, 2, 1, 0, 180.000000),
            Tors::new(6, 5, 4, 0, 36.366799),
        ];
        assert_eq!(got, want);
    }
}
