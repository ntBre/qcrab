use super::{contraction::Contraction, engine::Engine, shell::Shell, Operator};
use crate::{
    molecule::{Atom, Molecule},
    Dmat,
};
use std::ops::Index;

pub(crate) struct Basis(pub(crate) Vec<Shell>);

impl Basis {
    /// return the maximum number of primitives (alpha length) of `shells`. panics
    /// if called with empty `shells`
    pub(crate) fn max_nprim(&self) -> usize {
        self.0.iter().map(|s| s.alpha.len()).max().unwrap()
    }

    /// return the maximum angular momentum across the contractions of `shells`
    pub(crate) fn max_l(&self) -> usize {
        self.0
            .iter()
            .flat_map(|s| s.contr.iter().map(|c| c.l))
            .max()
            .unwrap()
    }

    pub(crate) fn nbasis(&self) -> usize {
        self.0.iter().map(|s| s.size()).sum()
    }

    /// the name pretty much says it don't you think?
    pub(crate) fn map_shell_to_basis_function(&self) -> Vec<usize> {
        let mut ret = Vec::with_capacity(self.len());
        let mut n = 0;
        for shell in &self.0 {
            ret.push(n);
            n += shell.size();
        }
        ret
    }

    /// STO-3G basis set
    ///
    /// cite: W. J. Hehre, R. F. Stewart, and J. A. Pople, The Journal of Chemical
    /// Physics 51, 2657 (1969) doi: 10.1063/1.1672392
    ///
    /// obtained from https://bse.pnl.gov/bse/portal via libint
    pub(crate) fn sto3g(mol: &Molecule) -> Self {
        let mut shells = Vec::new();
        for Atom {
            atomic_number,
            coord,
        } in &mol.atoms
        {
            match atomic_number {
                // hydrogen
                1 => shells.push(Shell::new(
                    vec![3.425250910, 0.623913730, 0.168855400],
                    vec![Contraction::new(
                        0,
                        false,
                        vec![0.15432897, 0.53532814, 0.44463454],
                    )],
                    *coord,
                    true,
                )),
                // carbon
                6 => shells.extend([
                    Shell::new(
                        vec![71.616837000, 13.045096000, 3.530512200],
                        vec![Contraction::new(
                            0,
                            false,
                            vec![0.15432897, 0.53532814, 0.44463454],
                        )],
                        *coord,
                        true,
                    ),
                    Shell::new(
                        vec![2.941249400, 0.683483100, 0.222289900],
                        vec![Contraction::new(
                            0,
                            false,
                            vec![-0.09996723, 0.39951283, 0.70011547],
                        )],
                        *coord,
                        true,
                    ),
                    Shell::new(
                        vec![2.941249400, 0.683483100, 0.222289900],
                        vec![Contraction::new(
                            1,
                            false,
                            vec![0.15591627, 0.60768372, 0.39195739],
                        )],
                        *coord,
                        true,
                    ),
                ]),
                // oxygen
                8 => shells.extend([
                    Shell::new(
                        vec![130.709320000, 23.808861000, 6.443608300],
                        vec![Contraction::new(
                            0,
                            false,
                            vec![0.15432897, 0.53532814, 0.44463454],
                        )],
                        *coord,
                        true,
                    ),
                    Shell::new(
                        vec![5.033151300, 1.169596100, 0.380389000],
                        vec![Contraction::new(
                            0,
                            false,
                            vec![-0.09996723, 0.39951283, 0.70011547],
                        )],
                        *coord,
                        true,
                    ),
                    Shell::new(
                        vec![5.033151300, 1.169596100, 0.380389000],
                        vec![Contraction::new(
                            1,
                            false,
                            vec![0.15591627, 0.60768372, 0.39195739],
                        )],
                        *coord,
                        true,
                    ),
                ]),
                _ => panic!("unsupported element {atomic_number}"),
            }
            // shells.push(
        }
        Self(shells)
    }

    pub(crate) fn len(&self) -> usize {
        self.0.len()
    }

    pub(crate) fn compute_1body_ints(&self, obtype: Operator) -> Dmat {
        let n = self.nbasis();
        let mut result = Dmat::zeros(n, n);
        let engine = Engine::new(obtype, self.max_nprim(), self.max_l(), 0);

        let shell2bf = self.map_shell_to_basis_function();

        for s1 in 0..self.len() {
            let bf1 = shell2bf[s1];
            let n1 = self[s1].size();

            for s2 in 0..=s1 {
                let bf2 = shell2bf[s2];
                let n2 = self[s2].size();

                println!("{} {}", s1, s2);
                // they do some weird buffer stuff where we look inside the buffer
                // previously returned by engine.result(), but I'll just return the
                // result each time, at least for now. I can see why you wouldn't
                // want to allocate a new vector on every call to compute, though
                let buf = engine.compute1(&self[s1], &self[s2]);
                if !buf.is_empty() {
                    result[(bf1, bf2)] = buf[0];
                    result[(bf2, bf1)] = buf[0];
                }
            }
        }
        result
    }

    pub(crate) fn overlap_ints(&self) -> Dmat {
        let n = self.nbasis();
        let mut result = Dmat::zeros(n, n);

        let shell2bf = self.map_shell_to_basis_function();

        for s1 in 0..self.len() {
            let bf1 = shell2bf[s1];
            let n1 = self[s1].size();

            for s2 in 0..=s1 {
                let bf2 = shell2bf[s2];
                let n2 = self[s2].size();

                println!("{} {}", s1, s2);
                // they do some weird buffer stuff where we look inside the buffer
                // previously returned by engine.result(), but I'll just return the
                // result each time, at least for now. I can see why you wouldn't
                // want to allocate a new vector on every call to compute, though
                let buf = compute1(&self[s1], &self[s2]);
                result[(bf1, bf2)] = buf;
                result[(bf2, bf1)] = buf;
            }
        }
        result
    }
}

pub(crate) fn compute1(s1: &Shell, s2: &Shell) -> f64 {
    assert!(
        s1.ncontr() == 1 && s2.ncontr() == 1,
        "generally-contracted shells not yet supported"
    );

    let l1 = s1.contr[0].l;
    let l2 = s2.contr[0].l;

    // we can skip the check on params because the nuclear charges are
    // required in the Nuclear operator

    let n1 = s1.size();
    let n2 = s2.size();
    let ncart1 = s1.cartesian_size();
    let ncart2 = s2.cartesian_size();
    let ncart12 = ncart1 * ncart2;

    // there are a bunch of checks here I think we can skip
    // (engine.impl.h:201) also skipping the scratch stuff

    let lmax = l1.max(l2);

    // N.B. for l=0 no need to transform to solid harmonics this is a
    // workaround for the corner case of oper_ == Operator::*nuclear, and
    // solid harmonics (s|s) integral ... beware the integral storage state
    // machine
    let tform_to_solids = (s1.contr[0].pure || s2.contr[0].pure) && lmax != 0;

    // TODO implement primdata, absolutely insane definition in libint
    let nprim1 = s1.nprim();
    let nprim2 = s2.nprim();

    let mut result = 0.0;
    for p1 in 0..nprim1 {
        for p2 in 0..nprim2 {
            let p = compute_primdata(s1, s2, p1, p2);
            result += p;
        }
    }

    result
}

/// this is only for overlap right now, need to inline it
pub(crate) fn compute_primdata(
    s1: &Shell,
    s2: &Shell,
    p1: usize,
    p2: usize,
) -> f64 {
    // return a libint_t which depends on something to be defined
    let a = s1.origin;
    let b = s2.origin;
    let alpha1 = s1.alpha[p1];
    let alpha2 = s2.alpha[p2];
    let c1 = s1.contr[0].coeff[p1];
    let c2 = s2.contr[0].coeff[p2];
    let gammap = alpha1 + alpha2;
    let oogammap = 1.0 / gammap;
    let rhop_over_alpha1 = alpha2 * oogammap;
    let rhop = alpha1 * rhop_over_alpha1;
    let px = (alpha1 * a[0] + alpha2 * b[0]) * oogammap;
    let py = (alpha1 * a[1] + alpha2 * b[1]) * oogammap;
    let pz = (alpha1 * a[2] + alpha2 * b[2]) * oogammap;
    let ab_x = a[0] - b[0];
    let ab_y = a[1] - b[1];
    let ab_z = a[2] - b[2];
    let ab2_x = ab_x * ab_x;
    let ab2_y = ab_y * ab_y;
    let ab2_z = ab_z * ab_z;

    let l1 = s1.contr[0].l;
    let l2 = s2.contr[0].l;

    // NOTE this is a nother operator variant that we don't have yet
    let is_sphemultipole = false;
    let hrr_ket_to_bra = l1 >= l2;

    if l1 != 0 || l2 != 0 {
        return f64::NAN;
    }

    const SQRT_PI: f64 = 1.772_453_850_905_516;
    let xyz_pfac: f64 = SQRT_PI * f64::sqrt(oogammap);
    let ovlp_ss_x = f64::exp(-rhop * ab2_x) * xyz_pfac * c1 * c2;
    let ovlp_ss_y = f64::exp(-rhop * ab2_y) * xyz_pfac;
    let ovlp_ss_z = f64::exp(-rhop * ab2_z) * xyz_pfac;

    ovlp_ss_x * ovlp_ss_y * ovlp_ss_z
}

impl Index<usize> for Basis {
    type Output = Shell;

    fn index(&self, index: usize) -> &Self::Output {
        self.0.index(index)
    }
}
