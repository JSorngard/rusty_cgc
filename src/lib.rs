#[cfg(test)]
mod truths;

use num_complex::Complex;

use std::cmp::Ordering;
use std::f64::consts::PI;

#[derive(Debug, Copy, Clone, PartialEq, Eq)]
enum Sign {
    Plus,
    Minus,
}

impl Sign {
    fn into_flipped(self) -> Self {
        match self {
            Self::Plus => Self::Minus,
            Self::Minus => Self::Plus,
        }
    }

    fn flip(&mut self) {
        *self = match self {
            Self::Plus => Self::Minus,
            Self::Minus => Self::Plus,
        };
    }
}

impl From<Sign> for f64 {
    fn from(s: Sign) -> Self {
        match s {
            Sign::Plus => 1.0,
            Sign::Minus => -1.0,
        }
    }
}

macro_rules! impl_from_sign_for_int {
    ($($t:ty),+) => {
        $(
            impl From<Sign> for $t {
                fn from(s: Sign) -> Self {
                    match s {
                        Sign::Plus => 1,
                        Sign::Minus => -1,
                    }
                }
            }
        )+
    };
}
impl_from_sign_for_int! {i8, i16, i32, i64, i128}

/// Returns the value of the Wigner 3j symbol for the given integer inputs.
/// The first three inputs are the angular momentum quantum
/// numbers, while the last three are the magnetic quantum numbers.
/// # Example
/// Basic usage
/// ```
/// # use rusty_cgc::{Error, wigner_3j};
/// assert_eq!(wigner_3j(1, 1, 1, 1, -1, 0)?, 1.0/f64::sqrt(6.0));
/// # Ok::<(), Error>(())
/// ```
/// There can be floating point errors for larger inputs, so the following assertion will panic:
/// ```should_panic
/// # use rusty_cgc::{Error, wigner_3j};
/// assert_eq!(wigner_3j(10, 10, 10, 8, 2, -10)?, 2.0 * f64::sqrt(561.0 / 723695.0));
/// # Ok::<(), Error>(())
/// ```
/// If we use float apropriate comparison instead the assert succeeds:
/// ```
/// use approx::assert_relative_eq;
/// # use rusty_cgc::{Error, wigner_3j};
/// assert_relative_eq!(wigner_3j(10, 10, 10, 8, 2, -10)?, 2.0 * f64::sqrt(561.0 / 723695.0));
/// # Ok::<(), Error>(())
/// ```
/// For larger inputs this error grows,
/// but for all valid 3j symbols with j1 and j2 at most equal to 10, this error is less than 100 * f64::EPSILON.
pub fn wigner_3j(j1: u32, j2: u32, j3: u32, m1: i32, m2: i32, m3: i32) -> Result<f64, Error> {
    is_unphysical(j1, j2, j3, m1, m2, -m3)?;

    let (j1, j2, j3, m1, m2, m3, mut sign) =
        reorder_3j_arguments(j1, j2, j3, m1, m2, m3, Sign::Plus);

    if (i64::from(j1) - i64::from(j2) - i64::from(m3)) % 2 != 0 {
        sign.flip();
    }

    let cg = clebsch_gordan(j1, j2, j3, m1, m2, -m3)?;

    Ok(f64::from(sign) * cg / (2.0 * f64::from(j3) + 1.0).sqrt())
}

/// Reorder j1/m1, j2/m2, j3/m3 such that j1 >= j2 >= j3 and m1 >= 0 or m1 == 0 && m2 >= 0
fn reorder_3j_arguments(
    j1: u32,
    j2: u32,
    j3: u32,
    m1: i32,
    m2: i32,
    m3: i32,
    sign: Sign,
) -> (u32, u32, u32, i32, i32, i32, Sign) {
    // An odd permutation of the columns or a
    // sign change of the m-quantum values (time reversal)
    // give a phase factor of (-1)^(j1+j2+j3).
    // If we assume that this phase factor is -1, we are only wrong if j1+j2+j3 is even,
    // which we correct for at the end.
    if j1 < j2 {
        reorder_3j_arguments(j2, j1, j3, m2, m1, m3, sign.into_flipped())
    } else if j2 < j3 {
        reorder_3j_arguments(j1, j3, j2, m1, m3, m2, sign.into_flipped())
    } else if m1 < 0 || (m1 == 0 && m2 < 0) {
        reorder_3j_arguments(j1, j2, j3, -m1, -m2, -m3, sign.into_flipped())
    } else {
        (
            j1,
            j2,
            j3,
            m1,
            m2,
            m3,
            if (j1 + j2 + j3) % 2 == 0 {
                // Sign doesn't matter if total J = j1 + j2 + j3 is even
                Sign::Plus
            } else {
                sign
            },
        )
    }
}

/// Returns the value of the Wigner 6j-symbol.
/// # Example
/// ```
/// # use rusty_cgc::{wigner_6j, Error};
/// use approx::assert_relative_eq;
/// assert_relative_eq!(wigner_6j(1, 1, 2, 1, 1, 0)?, 1.0 / 3.0);
///
/// // As the arguments increase the function becomes less accurate
/// assert_relative_eq!(
///     wigner_6j(5, 6, 7, 6, 7, 8)?,
///     1327.0 / (92378.0 * f64::sqrt(10.0)),
///     // Not accurate to within 1e-13
///     epsilon=1e-12
/// );
/// # Ok::<(), Error>(())
/// ```
pub fn wigner_6j(j1: u32, j2: u32, j3: u32, j4: u32, j5: u32, j6: u32) -> Result<f64, Error> {
    if !is_triad(j1, j2, j3)
        || !is_triad(j1, j5, j6)
        || !is_triad(j4, j2, j6)
        || !is_triad(j4, j5, j3)
    {
        return Err(Error::NotTriangular);
    }

    //The arguments to the factorials should always be positive
    let fac = delta(j2, j4, j6)?
        * delta(j2, j1, j3)?
        * delta(j6, j5, j1)?
        * delta(j4, j5, j3)?
        * ratio_of_factorials(
            &mut [j2 + j4 + j6 + 1, j4 + j5 + j3 + 1],
            &mut [
                j2 + j4 - j6,
                j6 + j1 - j5,
                j6 + j5 - j1,
                j2 + j3 - j1,
                j1 + j3 - j2,
                j4 + j5 - j3,
            ],
        )
        * phase(j4 + j6 + j1 + j3);
    let sum: f64 = (0..=(2 * j4).min(j4 + j6 - j2).min(j4 + j3 - j5))
        .map(|z| {
            ratio_of_factorials(
                &mut [2 * j4 - z, j4 + j6 + j3 - z - j1, j4 + j6 + j1 + j3 + 1 - z],
                &mut [
                    z,
                    j4 + j6 - z - j2,
                    j4 + j3 - z - j5,
                    j2 + j4 + j6 + 1 - z,
                    j4 + j5 + j3 + 1 - z,
                ],
            ) * phase(z)
        })
        .sum();
    Ok(sum * fac)
}

/// Returns the value of the Wigner 9j symbol.
/// # Note
/// Quickly becomes overwhelmed by floating point errors for inputs around 10.
pub fn wigner_9j(
    j1: u32,
    j2: u32,
    j3: u32,
    j4: u32,
    j5: u32,
    j6: u32,
    j7: u32,
    j8: u32,
    j9: u32,
) -> Result<f64, Error> {
    //Check that all rows are triads
    if !is_triad(j1, j2, j3) || !is_triad(j4, j5, j6) || !is_triad(j7, j8, j9) {
        return Err(Error::NotTriangular);
    }

    //Check that all columns are triads
    if !is_triad(j1, j4, j7) || !is_triad(j2, j5, j8) || !is_triad(j3, j6, j9) {
        return Err(Error::NotTriangular);
    }

    let prefactor = phase(j7 + j8 - j9) * nabla(j2, j1, j3)? / nabla(j2, j5, j8)?
        * nabla(j4, j5, j6)?
        / nabla(j4, j1, j7)?
        * nabla(j9, j3, j6)?
        / nabla(j9, j7, j8)?;

    let mut sum: f64 = 0.0;
    for x in 0..=(j5 + j8 - j2).min(j7 + j8 - j9) {
        for y in 0..=(j3 + j9 - j6).min(j4 + j5 - j6) {
            for z in 0..=(j1 + j7 - j4).min(j1 + j2 - j3) {
                if j4 + j9 + x + z < j1 + j8
                    || j2 + j6 + x + y < j4 + j8
                    || j1 + j2 + j9 < j6 + y + z
                {
                    continue;
                }
                sum += phase(x + y + z)
                    * ratio_of_factorials(
                        &mut [
                            2 * j8 - x,
                            j2 + j5 + x - j8,
                            j7 + j9 + x - j8,
                            j5 + j6 + y - j4,
                            j3 + j6 + y - j9,
                            j1 + j2 + j9 - y - z - j6,
                            2 * j1 - z,
                            j4 + j7 + z - j1,
                        ],
                        &mut [
                            x,
                            j5 + j8 - x - j2,
                            j7 + j8 - j9 - x,
                            j2 + j6 + x + y - j8 - j4,
                            j4 + j9 + x + z - j1 - j8,
                            y,
                            j4 + j5 - j6 - y,
                            j3 + j9 - y - j6,
                            2 * j6 + 1 + y,
                            z,
                            j1 + j2 - j3 - z,
                            j1 + j7 - z - j4,
                            j1 + j2 + j3 + 1 - z,
                        ],
                    );
            }
        }
    }
    Ok(prefactor * sum)
}

/// Returns the value of the Racah W coefficient.
pub fn racah_w(j1: u32, j2: u32, j: u32, j3: u32, j12: u32, j23: u32) -> Result<f64, Error> {
    Ok(phase(j1 + j2 + j3 + j) * wigner_6j(j1, j2, j12, j3, j, j23)?)
}

/// Returns the Gaunt coefficient for the input angular momenta.
/// The Gaunt coefficient is defined as the integral over three spherical harmonics:  
/// Y(l1, m1, θ, φ) * Y(l2, m2, θ, φ) * Y(l3, m3, θ, φ)
pub fn gaunt(l1: u32, l2: u32, l3: u32, m1: i32, m2: i32, m3: i32) -> Result<f64, Error> {
    is_unphysical(l1, l2, l3, m1, m2, -m3)?;

    Ok(
        ((2.0 * f64::from(l1) + 1.0) * (2.0 * f64::from(l2) + 1.0) * (2.0 * f64::from(l3) + 1.0)
            / (4.0 * PI))
            .sqrt()
            * wigner_3j(l1, l2, l3, 0, 0, 0)?
            * wigner_3j(l1, l2, l3, m1, m2, m3)?,
    )
}

/// Returns the value of the triangle coefficient \Delta(abc) used in the computation of 6j and 9j symbols.
fn delta(a: u32, b: u32, c: u32) -> Result<f64, Error> {
    if a + c >= b && a + b >= c && b + c >= a {
        Ok(
            ratio_of_factorials(&mut [a + c - b, a + b - c, c + b - a], &mut [a + c + b + 1])
                .sqrt(),
        )
    } else {
        Err(Error::NotTriangular)
    }
}

/// Returns the value of the \nabla(abc) triangle coefficient used in the computation of 9j symbols.
fn nabla(a: u32, b: u32, c: u32) -> Result<f64, Error> {
    if a + c >= b && a + b >= c && b + c >= a {
        Ok(
            ratio_of_factorials(&mut [a + c - b, a + b - c, a + b + c + 1], &mut [b + c - a])
                .sqrt(),
        )
    } else {
        Err(Error::NotTriangular)
    }
}

#[derive(Debug)]
pub struct WignerError;

impl std::fmt::Display for WignerError {
    fn fmt(&self, f: &mut std::fmt::Formatter) -> std::fmt::Result {
        write!(f, "m-values can not be larger than j")
    }
}

impl std::error::Error for WignerError {}

/// Returns the value of the small Wigner d-matrix in the z-y-z convention.
pub fn wigner_small_d(j: u32, mp: i32, m: i32, beta: f64) -> Result<f64, WignerError> {
    if mp.unsigned_abs() > j || m.unsigned_abs() > j {
        return Err(WignerError {});
    }

    let prefactor = f64::sqrt(
        factorial(j_plus_m(j, mp))
            * factorial(j_plus_m(j, -mp))
            * factorial(j_plus_m(j, m))
            * factorial(j_plus_m(j, -m)),
    );
    let mut sum: f64 = 0.0;
    for s in (m - mp).max(0).unsigned_abs()..=j_plus_m(j, m).min(j_plus_m(j, -mp)) {
        sum += f64::powi((beta / 2.0).cos(), (j_plus_m(j, m) + j_plus_m(j, -mp) - 2 * s) as i32)//(2 * j + m - mp)
            * f64::powi((beta / 2.0).sin(), 2 * s as i32 + mp - m)
            / (factorial(j_plus_m(j, m) - s)
                * factorial(s)
                * factorial((s as i32 + mp - m) as u32)
                * factorial(j_plus_m(j, -mp) - s))
            * phase((s as i32 + mp - m) as u32); //(-1)^x = (-1)^(-x) if x is real, and a positive i32 always fits in a u32.
    }
    Ok(sum * prefactor)
}

/// Returns the value of the Wigner D-matrix in the z-y-z convention.
pub fn wigner_d(
    j: u32,
    mp: i32,
    m: i32,
    alpha: f64,
    beta: f64,
    gamma: f64,
) -> Result<Complex<f64>, WignerError> {
    Ok(
        (wigner_small_d(j, mp, m, beta)? * -1.0 * Complex::<f64>::i() * f64::from(mp) * alpha)
            .exp()
            * (-1.0 * Complex::<f64>::i() * f64::from(m) * gamma).exp(),
    )
}

/// Returns the value of the Clebsch-Gordan coefficient for
/// the given integer inputs.
/// The first three inputs are the angular momentum quantum numbers,
/// while the last three are their projections.
pub fn clebsch_gordan(j1: u32, j2: u32, j3: u32, m1: i32, m2: i32, m3: i32) -> Result<f64, Error> {
    //This code is simply ported Fortran code,
    //as such it is not completely idiomatic rust.

    is_unphysical(j1, j2, j3, m1, m2, m3)?;

    if m1.unsigned_abs() + m2.unsigned_abs() == 0 && (j1 + j2 + j3) % 2 == 1 {
        return Ok(0.0);
    }

    let ia1 = j3 + j2 - j1;
    let ia2 = j_plus_m(j3, m3);
    let ia3: i64 = i64::from(j2) + i64::from(m3) - i64::from(j1);

    let ni = if ia3 > 0 {
        u32::try_from(ia3).unwrap()
    } else {
        0
    };
    let nm = if ia2 <= ia1 { ia2 } else { ia1 };

    let mut ip1 = if m1 >= 0 {
        j2 + j3 + (m1 as u32)
    } else {
        j2 + j3 - (-m1 as u32)
    } - ni; //j2 + j3 + m1 - ni
    let mut ip2 = j_plus_m(j1, -m1) + ni + 1;
    let mut ir1 = ni + 1;
    let mut ir2 = j3 + j2 - ni - j1;
    let mut ir3 = j_plus_m(j3, m3) - ni;
    let mut ir4 = if m3 < 0 {
        j1 + ni + 1 + (-m3 as u32)
    } else {
        j1 + ni + 1 - (m3 as u32)
    } - j2; //j1 + ni + 1 - j2 - m3

    let mut s1 = phase(ni + j_plus_m(j2, m2))
        * ratio_of_factorials(&mut [ip1, ip2 - 1], &mut [ir2, ni, ir3, ir4 - 1]);

    let n = nm - ni;
    let mut fa;
    if n != 0 {
        fa = s1;
        for _ in 1..=n {
            fa = -fa * f64::from(ip2) * f64::from(ir2) / f64::from(ip1) * f64::from(ir3)
                / (f64::from(ir1) * f64::from(ir4));
            s1 += fa;
            ip1 -= 1;
            ip2 += 1;
            ir1 += 1;
            ir2 -= 1;
            ir3 -= 1;
            ir4 += 1;
        }
    }

    let res = (f64::from(2 * j3 + 1)
        * ratio_of_factorials(
            &mut [j3 + j1 - j2, ia1, j1 + j2 - j3, ia2, j_plus_m(j3, -m3)],
            &mut [
                j1 + j2 + j3 + 1,
                j_plus_m(j1, -m1),
                j_plus_m(j2, -m2),
                j_plus_m(j2, m2),
                j_plus_m(j1, m1),
            ],
        ))
    .sqrt()
        * s1;

    if res.abs() < 1e-14 {
        //guard against floating point errors making a zero result non-zero
        Ok(0.0)
    } else {
        Ok(res)
    }
}

/// Returns whether the given triplet of angular momenta form a triad.
fn is_triad(j1: u32, j2: u32, j3: u32) -> bool {
    j3 >= j1.abs_diff(j2) && j3 <= j1 + j2
}

/// Returns the result of adding an angular momentum to its projection.
/// Debug_asserts that |m| <= j.
fn j_plus_m(j: u32, m: i32) -> u32 {
    debug_assert!(m.unsigned_abs() <= j);
    if m >= 0 {
        j + m.unsigned_abs()
    } else {
        j - m.unsigned_abs()
    }
}

///Returns the factorial of the input integer.
fn factorial(n: u32) -> f64 {
    (2..=n).map(f64::from).product()
}

/// Returns an f64 with a value of 1.0 if the input is even, and -1.0 if it is odd
fn phase(x: u32) -> f64 {
    if x % 2 == 0 {
        1.0
    } else {
        -1.0
    }
}

/// Returns whether the given quantum numbers represent something unphysical in a CG-coeff or 3j symbol.
fn is_unphysical(j1: u32, j2: u32, j3: u32, m1: i32, m2: i32, m3: i32) -> Result<(), Error> {
    if m1.unsigned_abs() > j1 && m2.unsigned_abs() > j2 && m3.unsigned_abs() > j3 {
        Err(Error::TooLargeM)
    } else if !is_triad(j1, j2, j3) {
        Err(Error::NotTriangular)
    } else if m1 + m2 != m3 {
        Err(Error::IncorrectMSum)
    } else {
        Ok(())
    }
}

#[derive(Debug)]
pub enum Error {
    TooLargeM,
    NotTriangular,
    IncorrectMSum,
}

impl std::fmt::Display for Error {
    fn fmt(&self, f: &mut std::fmt::Formatter) -> std::fmt::Result {
        match self {
            Self::TooLargeM => write!(f, "|m| is larger than its corresponding j"),
            Self::NotTriangular => {
                write!(f, "arguments do not fulfill the needed triangle conditions")
            }
            Self::IncorrectMSum => write!(f, "m1 + m2 do not equal m3"),
        }
    }
}

impl std::error::Error for Error {}

/// Takes in two lists of integers representing a ratio of two sets of factorials
/// and returns the value of that ratio as an `f64`.
/// E.g. an input of [5, 2, 7] and [3, 9] represents the equation (5!*2!*7!)/(3!*9!)
/// and would give a value of 5/9 = 0.5555...
///
/// Can handle large factorials accurately as long as both the numerator and denominator
/// have factorials of similar size.
///
/// Note that an empty vector is treated as if it contains a single 1.
fn ratio_of_factorials(numerators: &mut [u32], denominators: &mut [u32]) -> f64 {
    // In this function we pair up the arguments in the numerator and denominator
    // in order to find pairs of similar values.

    let number_of_numerators = numerators.len();
    let number_of_denominators = denominators.len();

    // We begin by ordering the terms in descending order
    // so the input [5, 2, 7], [3, 9] becomes [7, 5, 2], [9, 3]

    if number_of_numerators > 1 {
        numerators.sort_unstable();
        numerators.reverse();
    }

    if number_of_denominators > 1 {
        denominators.sort_unstable();
        denominators.reverse();
    }

    if numerators == denominators {
        return 1.0;
    }

    // Split into (numerator, denominator) pairs.
    // So [7, 5, 2], [9, 3] becomes the pairs (7, 9) and (5, 3).
    let res = numerators
        .iter()
        .zip(denominators.iter())
        .fold(1.0, |res, (n, d)| {
            match n.cmp(d) {
                // n > d => multiply by all the terms in n! that are not cancelled by dividing by d!.
                // The second pair in the example goes to this branch and we multiply the result by 5 * 4.
                Ordering::Greater => res * (d + 1..n + 1).map(f64::from).product::<f64>(),
                // n < d => divide by all the terms in d! that are not cancelled by multiplying by n!.
                // The first pair in the example goes to this branch and we divide the result by 9 * 8.
                Ordering::Less => res / (n + 1..d + 1).map(f64::from).product::<f64>(),
                // if n == d the terms cancel out completely, so we just leave the result as it is.
                Ordering::Equal => res,
            }
        });
    // After this the example result is 5/18 = 0.27777...

    // Deal with the unpaired numbers.
    // The example input has an unpaired 2 in the numerator.
    match number_of_numerators.cmp(&number_of_denominators) {
        Ordering::Greater => {
            // Multiply the result by the factorials of the remaining numeratrs.
            // For the example input this results in multiplying the result by 2.
            res * numerators
                .iter()
                .skip(number_of_denominators)
                .map(|n| factorial(*n))
                .product::<f64>()
        }
        Ordering::Less => {
            // Divide the result by the factorials of the remaining denominators.
            res / denominators
                .iter()
                .skip(number_of_numerators)
                .map(|d| factorial(*d))
                .product::<f64>()
        }
        // if the lengths are equal we have dealt with all terms in the loop, and can just return the result.
        Ordering::Equal => res,
    }
    // After this the result is 5/18 * 2 = 5/9 = 0.555555...
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::truths::return_3j_truths;
    use approx::assert_relative_eq;

    use std::f64::consts::PI;

    //We allow floating point errors on the scale of TOL.
    //The answers we compare against are exact,
    //but we use floats to compute the answers, so we will have
    //errors due to the formulas for 3j symbols containing alternating
    //products and divisions by large factorials.
    const TOL: f64 = 100.0 * f64::EPSILON;

    #[test]
    fn test_ratio_of_factorials() {
        assert_relative_eq!(ratio_of_factorials(&mut [5, 2, 7], &mut [3, 6]), 280.0);
        assert_relative_eq!(
            ratio_of_factorials(&mut [1, 2, 3], &mut [0, 5, 8]),
            1.0 / 403200.0
        );
        assert_relative_eq!(ratio_of_factorials(&mut [200], &mut [197]), 7880400.0);
        assert_relative_eq!(
            ratio_of_factorials(&mut [1000000], &mut [999999, 8]),
            3125.0 / 126.0
        );
        assert_eq!(ratio_of_factorials(&mut [], &mut []), 1.0);
    }

    #[test]
    fn test_bad_cgc_input() {
        assert!(clebsch_gordan(1, 1, 0, 1, 1, 0).is_err()); //m3 != m1 + m2
        assert!(clebsch_gordan(1, 2, 0, 0, 0, 0).is_err()); //|m1| > j1
        assert!(clebsch_gordan(2, 2, 5, 0, -1, -1).is_err()); //non-triangular
        assert_eq!(clebsch_gordan(1, 1, 1, 0, 0, 0).unwrap(), 0.0); //all ms 0 => j3 must be even
    }

    #[test]
    fn test_good_cgc_input() {
        assert_relative_eq!(
            clebsch_gordan(5, 5, 0, 1, -1, 0).unwrap(),
            f64::sqrt(1.0 / 11.0)
        );
        assert_relative_eq!(
            clebsch_gordan(1, 2, 1, 1, -1, 0).unwrap(),
            f64::sqrt(3.0 / 10.0)
        );
        assert_relative_eq!(
            clebsch_gordan(1, 1, 0, 1, -1, 0).unwrap(),
            f64::sqrt(1.0 / 3.0)
        );
        assert_relative_eq!(clebsch_gordan(1, 2, 3, 1, 2, 3).unwrap(), 1.0);
    }

    #[test]
    fn test_bad_3j_input() {
        assert!(wigner_3j(1, 1, 0, 1, 1, 0).is_err()); //m3 != m1 + m2
        assert!(wigner_3j(1, 2, 0, 0, 0, 0).is_err()); //|m| > j
        assert!(wigner_3j(2, 2, 5, 0, 0, 0).is_err()); //non-triangular
        assert_eq!(wigner_3j(1, 1, 1, 0, 0, 0).unwrap(), 0.0); //all ms = 0 => j3 must be even
    }

    #[test]
    fn test_good_3j_input() {
        println!("1");
        assert_relative_eq!(wigner_3j(1, 1, 0, 0, 0, 0).unwrap(), -f64::sqrt(1.0 / 3.0));
        println!("2");
        assert_relative_eq!(wigner_3j(1, 1, 2, -1, 1, 0).unwrap(), f64::sqrt(1.0 / 30.0));
        println!("3");
        assert_relative_eq!(wigner_3j(1, 2, 3, 1, 2, -3).unwrap(), f64::sqrt(1.0 / 7.0));
        println!("4");
        assert_relative_eq!(
            wigner_3j(10, 10, 9, 5, 4, -9).unwrap(),
            11.0 / 3.0 * f64::sqrt(91.0 / 126730.0)
        );
    }

    #[test]
    fn test_many_3js_to_within_relaxed_tolerance() {
        //Tests every valid combination of inputs with j1 and j2 <= 10

        //Load the right answers, computed with Mathematica.
        let threej_truths = return_3j_truths();

        //The maximum values of the first two angular momenta.
        const MAXJ: u32 = 10;

        for j1 in 0..=MAXJ {
            for j2 in 0..=MAXJ {
                for j3 in ((j1 as i32) - (j2 as i32)).unsigned_abs()..=j1 + j2 {
                    for m1 in -(j1 as i32)..=j1 as i32 {
                        for m2 in -(j2 as i32)..=j2 as i32 {
                            if (m1 + m2).unsigned_abs() > j3 {
                                continue;
                            }

                            let s = format!(
                                "(({}, {}), ({}, {}), ({}, {}))",
                                j1,
                                m1,
                                j2,
                                m2,
                                j3,
                                -m1 - m2
                            );

                            println!("{}", s);
                            assert_relative_eq!(
                                wigner_3j(j1, j2, j3, m1, m2, -m1 - m2).unwrap(),
                                threej_truths[s.as_str()],
                                epsilon = TOL,
                            );
                        }
                    }
                }
            }
        }
    }

    #[test]
    fn test_good_wigner_small_d_inputs() {
        let points = 12;
        for i in 0..=points {
            let arg = i as f64 * 2.0 * PI / points as f64;
            assert_relative_eq!(wigner_small_d(1, 0, 0, arg).unwrap(), arg.cos());
            assert_relative_eq!(
                wigner_small_d(1, 1, 1, arg).unwrap(),
                ((1.0 + arg.cos()) / 2.0)
            );
            assert_relative_eq!(
                wigner_small_d(2, 2, 0, arg).unwrap(),
                f64::sqrt(3.0 / 8.0) * arg.sin() * arg.sin()
            );
            assert_relative_eq!(
                wigner_small_d(2, 0, 2, arg).unwrap(),
                f64::sqrt(3.0 / 8.0) * arg.sin() * arg.sin()
            );
            assert_relative_eq!(
                wigner_small_d(2, 1, -1, arg).unwrap(),
                (-2.0 * arg.cos() * arg.cos() + arg.cos() + 1.0) / 2.0
            );
        }
    }

    #[test]
    fn test_good_wigner_d_inputs() {
        assert_relative_eq!(wigner_d(2, 1, -1, PI / 2.0, PI, -PI / 2.0).unwrap().re, 1.0);
        assert_relative_eq!(wigner_d(2, 1, -1, PI / 2.0, PI, -PI / 2.0).unwrap().im, 0.0);
    }

    #[test]
    fn test_good_6j_inputs() {
        assert_relative_eq!(
            wigner_6j(1, 2, 3, 4, 5, 6).unwrap(),
            f64::sqrt(2.0 / 715.0) / 3.0,
            epsilon = TOL
        );
        assert_relative_eq!(wigner_6j(1, 1, 2, 1, 1, 0).unwrap(), 1.0 / 3.0);
    }

    #[test]
    fn test_bad_9j_inputs() {
        assert!(wigner_9j(0, 0, 0, 0, 9, 0, 0, 0, 0).is_err());
        assert!(wigner_9j(100, 4, 6, 4, 6, 8, 6, 8, 10).is_err());
        assert!(wigner_9j(2, 4, 6, 4, 6, 8, 6, 8, 100).is_err());
    }

    #[test]
    fn test_good_9j_inputs() {
        assert_relative_eq!(
            wigner_9j(2, 4, 6, 4, 6, 8, 6, 8, 10).unwrap(),
            -37903.0 / 97274034.0
        );
        assert_relative_eq!(
            wigner_9j(2, 4, 6, 4, 6, 8, 6, 8, 9).unwrap(),
            199.0 / 1907334.0
        );
        assert_relative_eq!(
            wigner_9j(1, 2, 3, 2, 3, 4, 3, 4, 5).unwrap(),
            -61.0 / 66150.0
        );
        assert_relative_eq!(
            wigner_9j(1, 3, 2, 2, 3, 1, 2, 1, 3).unwrap(),
            -4.0 / 315.0 * f64::sqrt(2.0 / 35.0)
        );
        assert_relative_eq!(wigner_9j(1, 1, 2, 1, 1, 1, 2, 1, 1).unwrap(), -1.0 / 45.0);
        assert_relative_eq!(
            wigner_9j(1, 2, 3, 1, 2, 3, 2, 4, 6).unwrap(),
            (1.0 / 21.0) * f64::sqrt(0.2)
        );
        assert_relative_eq!(
            wigner_9j(5, 7, 6, 6, 8, 7, 7, 9, 8).unwrap(),
            8129.0 / 3136350672.0,
            // This value is harder to compute
            epsilon = 1e-11,
        );
        // These two are wildly inaccurate
        // assert_relative_eq!(
        //     wigner_9j(10, 11, 12, 11, 12, 13, 12, 13, 14).unwrap(),
        //     1.0 / 31648909588721100.0,
        // );
        // assert_relative_eq!(
        //     wigner_9j(11, 12, 13, 12, 13, 14, 13, 14, 15).unwrap(),
        //     -26976629723.0 / 599200846745523750.0
        // );
    }

    #[test]
    fn test_good_gaunt_inputs() {
        assert_relative_eq!(gaunt(1, 0, 1, 1, 0, -1).unwrap(), -1.0 / (2.0 * PI.sqrt()));
        assert_relative_eq!(
            gaunt(2, 1, 1, 1, 0, -1).unwrap(),
            -0.5 * (3.0 / (5.0 * PI)).sqrt()
        );
        assert_relative_eq!(
            gaunt(10, 4, 10, 2, 3, -5).unwrap(),
            1323.0 * (91.0 / (2.0 * PI)).sqrt() / 37145.0
        );
    }
}
