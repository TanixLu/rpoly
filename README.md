# rpoly

[![crates.io](https://img.shields.io/crates/v/rpoly.svg)](https://crates.io/crates/rpoly)
[![docs.rs](https://img.shields.io/docsrs/rpoly)](https://docs.rs/rpoly)
[![GitHub](https://img.shields.io/badge/GitHub-TanixLu/rpoly-blue)](https://github.com/TanixLu/rpoly)

`rpoly` is a Rust implementation of the [RPOLY algorithm](https://en.wikipedia.org/wiki/Jenkins%E2%80%93Traub_algorithm#:~:text=known%20as%20the%20%22-,RPOLY,-%22%20algorithm.%20The%20latter), also known as the Jenkins-Traub method. This algorithm is widely regarded for its robust and efficient computation of all roots (both real and complex) of a real-coefficient univariate polynomial.

---

## Features

- **Polynomial Solver**: Efficiently finds all roots of a univariate polynomial with real coefficients.
- **Iterative Interface**: Retrieve computed roots as a simple iterable collection.
- **Error Handling**: Provides clear error messages for cases like leading coefficient zero or non-convergence.

---

## Getting Started

### Installation

Add `rpoly` to your `Cargo.toml`:

```toml
[dependencies]
rpoly = "*"
```

### Example Usage

Solve the polynomial \(x^2 - 5x + 6 = 0\):

```rust
use rpoly::rpoly;

fn main() {
    let coefficients = [1.0, -5.0, 6.0]; // Represents x^2 - 5x + 6
    let roots = rpoly(&coefficients).unwrap();

    println!("Number of roots: {}", roots.root_count());
    for root in roots {
        println!("Root = {} + {}i", root.re, root.im);
    }
}
```

Output:
```
Number of roots: 2
Root = 1.9999999999999998 + 0i
Root = 3.0000000000000004 + 0i
```

---

## API Overview

### Main Function

#### `rpoly`

```rust
pub fn rpoly<const MDP1: usize>(
    a: &[f64; MDP1],
) -> Result<RpolyRoots<MDP1>, RpolyError>
```

Solves the polynomial:

a\[0\]*x^(n-1) + a\[1\]*x^(n-2) + ... + a\[n-2\]*x + a\[n-1\] = 0

**Parameters**:
- `a`: Array of coefficients, where `a.len() == MDP1` and the degree of the polynomial is \(MDP1 - 1\).

**Returns**:
- `Ok(RpolyRoots)`: Contains all roots of the polynomial as `RpolyComplex` values.
- `Err(RpolyError)`: Possible errors:
  - `RpolyLeadingCoefficientZero`: Leading coefficient `a[0]` is zero.
  - `RpolyNotConvergent`: The method failed to converge.

### Supporting Structures

- `RpolyRoots`: A collection of roots returned by the algorithm. It implements `IntoIterator` for easy iteration.
- `RpolyComplex`: Represents a complex number with fields `re` (real part) and `im` (imaginary part).

---

## Error Handling

The crate uses the `Result` type for error handling, returning `RpolyError` variants in cases where the computation cannot proceed:

- **`RpolyLeadingCoefficientZero`**: Raised if the leading coefficient (`a[0]`) is zero.
- **`RpolyNotConvergent`**: Raised if the algorithm fails to converge on a solution.

Example:

```rust
use rpoly::{rpoly, RpolyError};

let result = rpoly(&[0.0, 1.0, 2.0]); // Invalid: leading coefficient is zero
match result {
    Err(RpolyError::RpolyLeadingCoefficientZero) => {
        println!("Error: Leading coefficient cannot be zero.");
    }
    _ => {}
}
```

---

## Contributing

Contributions are welcome! Feel free to submit pull requests or issues on [GitHub](https://github.com/TanixLu/rpoly).

---

## License

The program was initially published as a FORTRAN program at [Netlib TOMS](https://www.netlib.org/toms/index.html#:~:text=ACM%2018%20200-,file%3A%20493.gz,-keywords%3A%20polynomial) and adhered to the [ACM software copyright notice](https://www.acm.org/publications/policies/software-copyright-notice). Later, it was translated into a C++ program available at [Akiti](https://www.akiti.ca/rpoly_ak1_Intro.html). I have translated this C++ program into Rust, and therefore, this program continues to comply with the ACM copyright policy. See LICENSE-ACM for more details.
