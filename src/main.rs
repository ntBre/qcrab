use std::{fs::File, io::BufRead, io::BufReader};

#[derive(Debug)]
struct Molecule {
    // atomic charges
    zs: Vec<i32>,
    // coordinates
    coords: Vec<f64>,
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
        mol.zs
            .push(vec[0].parse().expect("failed to parse atomic number"));
        for coord in &vec[1..4] {
            let coord: f64 = coord.parse().expect("failed to parse coord");
            mol.coords.push(coord);
        }
    }
    mol
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
                0.000000000000,
                0.000000000000,
                0.000000000000,
                0.000000000000,
                0.000000000000,
                2.845112131228,
                1.899115961744,
                0.000000000000,
                4.139062527233,
                -1.894048308506,
                0.000000000000,
                3.747688672216,
                1.942500819960,
                0.000000000000,
                -0.701145981971,
                -1.007295466862,
                -1.669971842687,
                -0.705916966833,
                -1.007295466862,
                1.669971842687,
                -0.705916966833,
            ],
        };
        assert_eq!(got.zs, want.zs, "got {:?}, wanted {:?}", got.zs, want.zs);
        assert_eq!(got.coords, want.coords);
    }
}

fn main() {
    let geomfile = "../inp/geom.xyz";
    let mol = read_geom(geomfile);
    println!("{:?}", mol);
}
