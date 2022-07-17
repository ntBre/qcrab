use std::{
    io::{BufRead, BufReader},
    ops::Index,
};

pub(crate) struct Eri {
    inner: Vec<f64>,
}

impl Index<(usize, usize, usize, usize)> for Eri {
    type Output = f64;

    fn index(&self, index: (usize, usize, usize, usize)) -> &Self::Output {
        let (m, n, l, s) = index;
        &self.inner[Self::index_inner(m, n, l, s)]
    }
}

impl Eri {
    /// returns the compound index associated with `m`, `n`, `l`, and `s`
    fn index_inner(m: usize, n: usize, l: usize, s: usize) -> usize {
        let compound = |a: usize, b: usize| a * (a + 1) / 2 + b;
        let compare = |a: usize, b: usize| if a < b { (b, a) } else { (a, b) };
        let (m, n) = compare(m, n);
        let (l, s) = compare(l, s);
        let mn = compound(m, n);
        let ls = compound(l, s);
        let (mn, ls) = compare(mn, ls);
        compound(mn, ls)
    }

    pub fn new(filename: &str) -> Self {
        let mut inner = Vec::new();
        let f = std::fs::File::open(filename).unwrap();
        let lines = BufReader::new(f).lines();
        for line in lines.flatten() {
            let split: Vec<_> = line.split_whitespace().collect();
            if split.len() != 5 {
                continue;
            }
            let m: usize = split[0].parse().unwrap();
            let n: usize = split[1].parse().unwrap();
            let l: usize = split[2].parse().unwrap();
            let s: usize = split[3].parse().unwrap();
            let v: f64 = split[4].parse().unwrap();
            let index = Self::index_inner(m - 1, n - 1, l - 1, s - 1);
            if index >= inner.len() {
                inner.resize(index + 1, 0.0);
            }
            inner[index] = v;
        }
        Self { inner }
    }
}
