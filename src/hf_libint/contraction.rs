#[derive(PartialEq, Debug)]
pub(crate) struct Contraction {
    /// angular momentum
    pub(crate) l: usize,

    /// spherical(?)
    pub(crate) pure: bool,

    /// contraction coefficients
    pub(crate) coeff: Vec<f64>,
}

impl Contraction {
    pub(crate) fn new(l: usize, pure: bool, coeff: Vec<f64>) -> Self {
        Self { l, pure, coeff }
    }

    pub(crate) const fn cartesian_size(&self) -> usize {
        (self.l + 1) * (self.l + 2) / 2
    }

    pub(crate) const fn size(&self) -> usize {
        if self.pure {
            2 * self.l + 1
        } else {
            self.cartesian_size()
        }
    }
}
