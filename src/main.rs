use na::Vector3;
use nalgebra as na;
use std::{fs::File, io::BufRead, io::BufReader};

#[derive(Debug)]
struct Molecule {
    // atomic charges
    zs: Vec<i32>,
    // coordinates
    coords: Vec<Vector3<f64>>,
}

impl Molecule {
    fn new() -> Molecule {
        Molecule {
            zs: Vec::new(),
            coords: Vec::new(),
        }
    }
}

fn read_geom(geomfile: &str) -> Molecule {
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

fn bond_lengths(mol: Molecule) -> Vec<(usize, usize, f64)> {
    let mut ret = Vec::new();
    for (i, atom) in mol.coords.iter().enumerate() {
        for (j, btom) in mol.coords[i + 1..].iter().enumerate() {
            let diff = (atom - btom).magnitude();
            ret.push((j + i + 1, i, diff));
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
        let got = bond_lengths(read_geom("inp/geom.xyz"));
        let want = vec![
            (1, 0, 2.84511),
            (2, 0, 4.55395),
            (3, 0, 4.19912),
            (4, 0, 2.06517),
            (5, 0, 2.07407),
            (6, 0, 2.07407),
            (2, 1, 2.29803),
            (3, 1, 2.09811),
            (4, 1, 4.04342),
            (5, 1, 4.05133),
            (6, 1, 4.05133),
            (3, 2, 3.81330),
            (4, 2, 4.84040),
            (5, 2, 5.89151),
            (6, 2, 5.89151),
            (4, 3, 5.87463),
            (5, 3, 4.83836),
            (6, 3, 4.83836),
            (5, 4, 3.38971),
            (6, 4, 3.38971),
            (6, 5, 3.33994),
        ];
        assert_eq!(got.len(), want.len());
        assert_eq!(
            got.iter().map(|x| x.0).collect::<Vec<usize>>(),
            want.iter().map(|x| x.0).collect::<Vec<usize>>(),
        );
        assert_eq!(
            got.iter().map(|x| x.1).collect::<Vec<usize>>(),
            want.iter().map(|x| x.1).collect::<Vec<usize>>(),
        );
        let eps = 1e-5;
        for i in 0..got.len() {
            assert!(
                got[i].2 - want[i].2 < eps,
                "got {}, wanted {}",
                got[i].2,
                want[i].2
            );
        }
    }
}

fn main() {
    let geomfile = "../inp/geom.xyz";
    let mol = read_geom(geomfile);
    let lens = bond_lengths(mol);
    println!("{:?}", lens);
}
