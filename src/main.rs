use qcrab::molecule::Molecule;
use qcrab::{self};

fn main() {
    let geomfile = "inp/geom.xyz";
    let mol = Molecule::load(geomfile);
    let lens = mol.bond_lengths();
    let angs = mol.bond_angles();
    let opbs = mol.oop_angles();
    let tors = mol.torsional_angles();
    println!("{:?}, {:?}, {:?}, {:?}", lens, angs, opbs, tors);
}
