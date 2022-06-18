use qcrab::{
    self, bond_angles, bond_lengths, oop_angles, torsional_angles, Molecule,
};

fn main() {
    let geomfile = "inp/geom.xyz";
    let mol = Molecule::load(geomfile);
    let lens = bond_lengths(&mol);
    let angs = bond_angles(&mol);
    let opbs = oop_angles(&mol);
    let tors = torsional_angles(&mol);
    println!("{:?}, {:?}, {:?}, {:?}", lens, angs, opbs, tors);
}
