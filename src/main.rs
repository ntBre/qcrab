use qcrab;

fn main() {
    let geomfile = "../inp/geom.xyz";
    let mol = qcrab::read_geom(geomfile);
    let lens = qcrab::bond_lengths(&mol);
    let angs = qcrab::bond_angles(&mol);
    let opbs = qcrab::op_angles(&mol);
    let tors = qcrab::torsional_angles(&mol);
    println!("{:?}, {:?}, {:?}, {:?}", lens, angs, opbs, tors);
}
