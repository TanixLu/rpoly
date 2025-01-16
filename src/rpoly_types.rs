#[derive(Debug, Clone, Copy, Default)]
pub struct RpolyComplex {
    pub re: f64,
    pub im: f64,
}

impl RpolyComplex {
    pub fn is_real(&self) -> bool {
        self.im == 0.0
    }

    pub fn is_imaginary(&self) -> bool {
        !self.is_real()
    }

    pub fn is_positive(&self) -> bool {
        self.is_real() && self.re > 0.0
    }

    pub fn is_negtive(&self) -> bool {
        self.is_real() && self.re < 0.0
    }

    pub fn is_non_positive(&self) -> bool {
        self.is_real() && self.re <= 0.0
    }

    pub fn is_non_negtive(&self) -> bool {
        self.is_real() && self.re >= 0.0
    }

    pub fn is_zero(&self) -> bool {
        self.re == 0.0 && self.im == 0.0
    }

    pub fn abs(&self) -> f64 {
        if self.is_real() {
            self.re
        } else if self.re == 0.0 {
            self.im
        } else {
            (self.re * self.re + self.im * self.im).sqrt()
        }
    }
}

#[derive(Debug)]
pub enum RpolyError {
    RpolyLeadingCoefficientZero,
    RpolyNotConvergent,
}

pub struct RpolyRoots<const MDP1: usize>(pub [RpolyComplex; MDP1]);

impl<const MDP1: usize> RpolyRoots<MDP1> {
    pub fn root_count(&self) -> usize {
        MDP1 - 1
    }

    pub fn to_vec(&self) -> Vec<RpolyComplex> {
        self.0[..MDP1 - 1].to_vec()
    }
}

pub struct RpolyRootsIntoIterator<const MDP1: usize> {
    data: [RpolyComplex; MDP1],
    index: usize,
}

impl<const MDP1: usize> Iterator for RpolyRootsIntoIterator<MDP1> {
    type Item = RpolyComplex;

    fn next(&mut self) -> Option<Self::Item> {
        if self.index < MDP1 {
            let result = self.data[self.index];
            self.index += 1;
            Some(result)
        } else {
            None
        }
    }
}

impl<const MDP1: usize> IntoIterator for RpolyRoots<MDP1> {
    type Item = RpolyComplex;
    type IntoIter = RpolyRootsIntoIterator<MDP1>;

    fn into_iter(self) -> Self::IntoIter {
        Self::IntoIter {
            data: self.0,
            index: 0,
        }
    }
}
