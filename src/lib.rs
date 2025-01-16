//! # Rpoly - Polynomial Root Finding in Rust
//!
//! This crate provides a Rust interface for finding the roots of real-coefficient polynomials using
//! an adapted version of the RPOLY algorithm. The core logic in [`rpoly_ak1.rs`] is manually converted
//! from the C++ code found at <https://www.akiti.ca/rpoly_ak1_Intro.html>.
//!
//! The primary entry point is the [`rpoly`] function, which returns the roots of a polynomial of
//! degree \(n - 1\). The polynomial is expressed in the form
//!
//! a\[0\]*x^(n-1) + a\[1\]*x^(n-2) + ... + a\[n-2\]*x + a\[n-1\] = 0
//!
//! where `n == MDP1`. Thus, `a.len()` must be `MDP1`, and the degree is `MDP1 - 1`. It is essential
//! that `a[0] != 0.0` (the leading coefficient must not be zero).
//!
//! ## Usage
//!
//! 1. Include this crate in your `Cargo.toml`.
//! 2. Call [`rpoly`] with an array of type `[f64; MDP1]` representing the polynomial coefficients.
//! 3. If successful, [`rpoly`] returns an [`RpolyRoots`] struct, which you can iterate to get
//!    [`RpolyComplex`] roots (both real and imaginary parts).
//!
//! ## Examples
//!
//! ```rust
//! use rpoly::{rpoly, RpolyComplex};
//!
//! // Solve x^2 - 5x + 6 = 0, i.e., roots should be x = 2.0 and x = 3.0
//! let coefficients = [1.0, -5.0, 6.0];
//! let roots = rpoly(&coefficients).unwrap();
//! assert_eq!(roots.root_count(), 2);
//!
//! for root in roots {
//!     println!("Root = {} + {}i", root.re, root.im);
//! }
//! ```
//!
//! ## Error Handling
//!
//! If the leading coefficient `a[0]` is zero, you will get an error of type [`RpolyError::RpolyLeadingCoefficientZero`].
//! If the method fails to converge, it returns [`RpolyError::RpolyNotConvergent`]. Handle these cases
//! by matching against the `Result` type returned by [`rpoly`].

mod rpoly_ak1;

use rpoly_ak1::{rpoly_ak1, LEADING_COEFFICIENT_ZERO_NUM, NOT_CONVERGENT_NUM};

#[derive(Debug, Clone, Copy, Default)]
pub struct RpolyComplex {
    pub re: f64,
    pub im: f64,
}

pub struct RpolyRoots<const MDP1: usize>(pub [RpolyComplex; MDP1]);

pub struct RpolyRootsIntoIterator<const MDP1: usize> {
    data: [RpolyComplex; MDP1],
    index: usize,
}

#[derive(Debug)]
pub enum RpolyError {
    RpolyLeadingCoefficientZero,
    RpolyNotConvergent,
}

pub use RpolyError::RpolyLeadingCoefficientZero;
pub use RpolyError::RpolyNotConvergent;

impl<const MDP1: usize> Iterator for RpolyRootsIntoIterator<MDP1> {
    type Item = RpolyComplex;

    fn next(&mut self) -> Option<Self::Item> {
        let result = if self.index + 1 < MDP1 {
            self.data[self.index]
        } else {
            return None;
        };
        self.index += 1;
        Some(result)
    }
}

impl<const MDP1: usize> RpolyRoots<MDP1> {
    pub fn root_count(&self) -> usize {
        MDP1 - 1
    }

    pub fn to_vec(&self) -> Vec<RpolyComplex> {
        self.0[..MDP1 - 1].to_vec()
    }
}

impl<const MDP1: usize> IntoIterator for RpolyRoots<MDP1> {
    type Item = RpolyComplex;

    type IntoIter = RpolyRootsIntoIterator<MDP1>;

    fn into_iter(self) -> Self::IntoIter {
        RpolyRootsIntoIterator {
            data: self.0,
            index: 0,
        }
    }
}

/// Solve the polynomial given by:
///
/// a\[0\]*x^(n-1) + a\[1\]*x^(n-2) + ... + a\[n-2\]*x + a\[n-1\] = 0
///
/// where `n == a.len() == MDP1`.
///
/// # Errors
///
/// - [`RpolyLeadingCoefficientZero`]: If `a[0] == 0.0`.
/// - [`RpolyNotConvergent`]: If the method fails to converge.
///
/// # Examples
///
/// ```
/// # use rpoly::rpoly;
/// #
/// let roots = rpoly(&[1.0, -3.0, 2.0]).unwrap();  // x^2 - 3x + 2 = 0
/// assert_eq!(roots.root_count(), 2);
/// ```
///
pub fn rpoly<const MDP1: usize>(a: &[f64; MDP1]) -> Result<RpolyRoots<MDP1>, RpolyError> {
    let mut degree = MDP1 - 1;
    let mut zeror = [0.0; MDP1];
    let mut zeroi = [0.0; MDP1];
    rpoly_ak1(a, &mut degree, &mut zeror, &mut zeroi);
    match degree {
        LEADING_COEFFICIENT_ZERO_NUM => Err(RpolyLeadingCoefficientZero),
        NOT_CONVERGENT_NUM => Err(RpolyNotConvergent),
        _ => {
            let mut roots = [RpolyComplex::default(); MDP1];
            for i in 0..degree {
                roots[i].re = zeror[i];
                roots[i].im = zeroi[i];
            }
            Ok(RpolyRoots(roots))
        }
    }
}

#[test]
fn test_rpoly() {
    let roots = rpoly(&[
        1.0,
        -687712.0,
        -3.12295507665e11,
        7.503861976776e12,
        3.803023224e15,
        0.0,
    ])
    .unwrap();

    assert!(roots.root_count() == 5);
    for root in roots {
        assert!(root.im == 0.0);
        // dbg!(root);
    }
}
