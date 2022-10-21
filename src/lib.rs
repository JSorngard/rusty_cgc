#[cfg(test)]
#[macro_use]
extern crate approx;

#[cfg(test)]
mod truths;

use itertools::Itertools;
use num::complex::Complex;
use std::f64::consts::PI;

/// Returns the value of the Wigner 3j symbol for the given integer inputs.
/// The first three inputs are the angular momentum quantum
/// numbers, while the last three are the magnetic quantum numbers.
/// # Example
/// Basic usage
/// ```
/// # use rusty_cgc::wigner_3j;
/// assert_eq!(wigner_3j(1, 1, 1, 1, -1, 0).unwrap(), 1.0/f64::sqrt(6.0));
/// ```
/// There can be floating point errors for larger inputs, so the following assertion will panic:
/// ```should_panic
/// # use rusty_cgc::wigner_3j;
/// assert_eq!(wigner_3j(10, 10, 10, 8, 2, -10).unwrap(), 2.0*f64::sqrt(561.0/723695.0));
/// ```
/// If we use float apropriate comparison instead the assert succeeds:
/// ```
/// #[macro_use]
/// extern crate approx;
/// # use rusty_cgc::wigner_3j;
/// # fn main() {
/// assert_relative_eq!(wigner_3j(10, 10, 10, 8, 2, -10).unwrap(), 2.0*f64::sqrt(561.0/723695.0));
/// # }
/// ```
/// For certain inputs this error can compound during computation.
/// For all valid 3j symbols with j1 and j2 at most equal to 10, this error is less than 100 * f64::EPSILON.
pub fn wigner_3j(j1: u32, j2: u32, j3: u32, m1: i32, m2: i32, m3: i32) -> Result<f64, String> {
    if let Some(e) = is_unphysical(j1, j2, j3, m1, m2, -m3) {
        return Err(e);
    };

    let (j1, j2, j3, m1, m2, m3, mut sign) = reorder3j(j1, j2, j3, m1, m2, m3, 1);

    if (i64::from(j1) - i64::from(j2) - i64::from(m3)) % 2 != 0 {
        sign *= -1
    }

    let cg = clebsch_gordan(j1, j2, j3, m1, m2, -m3)?;

    Ok(f64::from(sign) * cg / (2.0 * f64::from(j3) + 1.0).sqrt())
}

/// Reorder j1/m1, j2/m2, j3/m3 such that j1 >= j2 >= j3 and m1 >= 0 or m1 == 0 && m2 >= 0
fn reorder3j(
    j1: u32,
    j2: u32,
    j3: u32,
    m1: i32,
    m2: i32,
    m3: i32,
    sign: i8,
) -> (u32, u32, u32, i32, i32, i32, i8) {
    //An odd permutation of the columns or a
    //sign change of the m-quantum values (time reversal)
    //give a phase factor of (-1)^(j1+j2+j3).
    //If we assume that this phase factor is -1, we are only wrong if j1+j2+j3 is even,
    //which we correct for at the end.
    if j1 < j2 {
        reorder3j(j2, j1, j3, m2, m1, m3, -sign)
    } else if j2 < j3 {
        reorder3j(j1, j3, j2, m1, m3, m2, -sign)
    } else if m1 < 0 || (m1 == 0 && m2 < 0) {
        reorder3j(j1, j2, j3, -m1, -m2, -m3, -sign)
    } else {
        (
            j1,
            j2,
            j3,
            m1,
            m2,
            m3, //Sign doesn't matter if total J = j1 + j2 + j3 is even
            if (j1 + j2 + j3) % 2 == 0 { 1 } else { sign },
        )
    }
}

/// Returns the value of the Wigner 6j-symbol.
pub fn wigner_6j(j1: u32, j2: u32, j3: u32, j4: u32, j5: u32, j6: u32) -> Result<f64, String> {
    if !is_triad(j1, j2, j3)
        || !is_triad(j1, j5, j6)
        || !is_triad(j4, j2, j6)
        || !is_triad(j4, j5, j3)
    {
        return Err("the inputs do not fulfill the required triangle conditions".to_owned());
    }

    //The arguments to the factorials should always be positive
    let fac = delta(j2, j4, j6)
        * delta(j2, j1, j3)
        * delta(j6, j5, j1)
        * delta(j4, j5, j3)
        * ratio_of_factorials(
            vec![j2 + j4 + j6 + 1, j4 + j5 + j3 + 1],
            vec![
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
                vec![2 * j4 - z, j4 + j6 + j3 - z - j1, j4 + j6 + j1 + j3 + 1 - z],
                vec![
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

/// This function fails for some inputs, and I have not figured out why yet
fn wigner_9j(
    j1: u32,
    j2: u32,
    j3: u32,
    j4: u32,
    j5: u32,
    j6: u32,
    j7: u32,
    j8: u32,
    j9: u32,
) -> Result<f64, String> {
    //Check that all rows are triads
    if !is_triad(j1, j2, j3) || !is_triad(j4, j5, j6) || !is_triad(j7, j8, j9) {
        return Err("A row does not fulfill the triangle conditions".to_owned());
    }

    //Check that all columns are triads
    if !is_triad(j1, j4, j7) || !is_triad(j2, j5, j8) || !is_triad(j3, j6, j9) {
        return Err("A column does not fulfill the triangle conditions".to_owned());
    }

    let prefactor = phase(j7 + j8 - j9) * nabla(j2, j1, j3) / nabla(j2, j5, j8) * nabla(j4, j5, j6)
        / nabla(j4, j1, j7)
        * nabla(j9, j3, j6)
        / nabla(j9, j7, j8);

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
                        vec![
                            2 * j8 - x,
                            j2 + j5 + x - j8,
                            j7 + j9 + x - j8,
                            j5 + j6 + y - j4,
                            j3 + j6 + y - j9,
                            j1 + j2 + j9 - y - z - j6,
                            2 * j1 - z,
                            j4 + j7 + z - j1,
                        ],
                        vec![
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
pub fn racah_w(j1: u32, j2: u32, j: u32, j3: u32, j12: u32, j23: u32) -> Result<f64, String> {
    Ok(phase(j1 + j2 + j3 + j) * wigner_6j(j1, j2, j12, j3, j, j23)?)
}

/// Returns the Gaunt coefficient for the input angular momenta.
/// The Gaunt coefficient is defined as the integral over three spherical harmonics.
pub fn gaunt(l1: u32, l2: u32, l3: u32, m1: i32, m2: i32, m3: i32) -> Result<f64, String> {
    if let Some(e) = is_unphysical(l1, l2, l3, m1, m2, -m3) {
        return Err(e);
    };

    let w1 = wigner_3j(l1, l2, l3, 0, 0, 0)?;
    let w2 = wigner_3j(l1, l2, l3, m1, m2, m3)?;

    Ok(
        ((2.0 * f64::from(l1) + 1.0) * (2.0 * f64::from(l2) + 1.0) * (2.0 * f64::from(l3) + 1.0)
            / (4.0 * PI))
            .sqrt()
            * w1
            * w2,
    )
}

/// Returns the value of the triangle coefficient \Delta(abc) used in the computation of 6j and 9j symbols.
/// The inputs are debug_asserted to be triangular.
fn delta(a: u32, b: u32, c: u32) -> f64 {
    debug_assert!(a + c >= b && a + b >= c && b + c >= a);
    ratio_of_factorials(vec![a + c - b, a + b - c, c + b - a], vec![a + c + b + 1]).sqrt()
}

/// Returns the value of the \nabla(abc) triangle coefficient used in the computation of 9j symbols.
/// The inputs are debug_asserted to be triangular.
fn nabla(a: u32, b: u32, c: u32) -> f64 {
    debug_assert!(a + c >= b && a + b >= c && b + c >= a);
    ratio_of_factorials(vec![a + c - b, a + b - c, a + b + c + 1], vec![b + c - a]).sqrt()
    // (factorial((a + c - b).into()) * factorial((a + b - c).into()) / factorial((b + c - a).into())
    //     * factorial((a + b + c + 1).into()))
    // .sqrt()
}

/// Returns the value of the small Wigner d-matrix in the z-y-z convention.
pub fn wigner_small_d(j: u32, mp: i32, m: i32, beta: f64) -> Result<f64, String> {
    //abs(i32) always fits in a u32.
    if u32::try_from(mp.abs()).unwrap() > j || u32::try_from(m.abs()).unwrap() > j {
        return Err("m-values can not be larger than j".to_owned());
    }

    let prefactor = f64::sqrt(
        factorial(j_plus_m(j, mp).into())
            * factorial(j_plus_m(j, -mp).into())
            * factorial(j_plus_m(j, m).into())
            * factorial(j_plus_m(j, -m).into()),
    );
    let mut sum: f64 = 0.0;
    for s in (i32::max(0, m - mp) as u32)..=u32::min(j_plus_m(j, m), j_plus_m(j, -mp)) {
        sum += f64::powf((beta / 2.0).cos(), (j_plus_m(j, m) + j_plus_m(j, -mp) - 2 * s).into())//(2 * j + m - mp)
            * f64::powf((beta / 2.0).sin(), (2 * s + u32::try_from(mp - m).unwrap()).into())
            / (factorial((j_plus_m(j, m) - s).try_into().unwrap())
                * factorial(s.try_into().unwrap())
                * factorial((s + u32::try_from(mp - m).unwrap()).try_into().unwrap())
                * factorial((j_plus_m(j, -mp) - s).try_into().unwrap()))
            * phase(s + u32::try_from(mp - m).unwrap()); //(-1)^x = (-1)^(-x) if x is real, and a positive i32 always fits in a u32.
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
) -> Result<Complex<f64>, String> {
    Ok(
        (wigner_small_d(j, mp, m, beta)? * -1.0 * Complex::<f64>::i() * (mp as f64) * alpha).exp()
            * (-1.0 * Complex::<f64>::i() * (m as f64) * gamma).exp(),
    )
}

/// Returns the value of the Clebsch-Gordan coefficient for
/// the given integer inputs.
/// The first three inputs are the angular momentum quantum numbers,
/// while the last three are the magnetic quantum numbers.
pub fn clebsch_gordan(j1: u32, j2: u32, j3: u32, m1: i32, m2: i32, m3: i32) -> Result<f64, String> {
    //This code is simply ported Fortran code,
    //as such it is not completely idiomatic rust.

    if let Some(e) = is_unphysical(j1, j2, j3, m1, m2, m3) {
        return Err(e);
    };

    if m1.abs() + m2.abs() == 0 && (j1 + j2 + j3) % 2 == 1 {
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

    let cc = (f64::from(2 * j3 + 1)
        * ratio_of_factorials(
            vec![j3 + j1 - j2, ia1, j1 + j2 - j3, ia2, j_plus_m(j3, -m3)],
            vec![
                j1 + j2 + j3 + 1,
                j_plus_m(j1, -m1),
                j_plus_m(j2, -m2),
                j_plus_m(j2, m2),
                j_plus_m(j1, m1),
            ],
        ))
    .sqrt();

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
        * ratio_of_factorials(vec![ip1, ip2 - 1], vec![ir2, ni, ir3, ir4 - 1]);

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

    let res = cc * s1;
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
    debug_assert!(u32::try_from(m.abs()).unwrap() <= j);
    if m >= 0 {
        j + u32::try_from(m).unwrap()
    } else {
        j - u32::try_from(-m).unwrap()
    }
}

///Returns whether the given triplet of angular momenta form a triad.
fn is_float_triad(j1: f32, j2: f32, j3: f32) -> bool {
    //normal triangle condition                j1 + j2 + j3 must be an integer
    j3 >= (j1 - j2).abs() && j3 <= j1 + j2 && (j1 + j1 + j3).fract() == 0.0
}

///Returns the factorial of the input integer as a float
fn factorial(n: u64) -> f64 {
    if n == 0 {
        1.0
    } else {
        let mut res = 1.0;
        for i in 2..=n {
            res *= i as f64;
        }
        res
    }
}

///Returns an f64 with a value of 1.0 if the input is even, and -1.0 if it is odd
fn phase(x: u32) -> f64 {
    if x % 2 == 0 {
        1.0
    } else {
        -1.0
    }
}

///Returns whether the given quantum numbers represent something unphysical in a CG-coeff or 3j symbol.
fn is_unphysical(j1: u32, j2: u32, j3: u32, m1: i32, m2: i32, m3: i32) -> Option<String> {
    if m1.unsigned_abs() > j1 && m2.unsigned_abs() > j2 && m3.unsigned_abs() > j3 {
        Some("|m| is larger than its corresponding j".to_owned())
    } else if !is_triad(j1, j2, j3) {
        Some("j1, j2, and j3 do not fulfill the triangle condition".to_owned())
    } else if m1 + m2 != m3 {
        Some("m1 + m2 do not equal m3".to_owned())
    } else {
        None
    }
}

/// Takes in two lists of integers representing a ratio of two sets of factorials
/// and returns the value of that ratio as an `f64`.
/// E.g. an input of [5, 2, 7] and [3, 6] represents the equation (5!*2!*7!)/(3!*6!)
/// and would give a value of 280.0
///
/// Can handle large factorials, as long as both the numerator and denominator
/// have factorials of similar size.
///
/// # Example
/// ```
/// # use rusty_cgc::ratio_of_factorials;
/// assert_eq!(ratio_of_factorials(vec![1000000], vec![999999, 8]), 3125.0 / 126.0);
/// ```
pub fn ratio_of_factorials(mut numerators: Vec<u32>, mut denominators: Vec<u32>) -> f64 {
    // In this function we pair up the arguments in the numerator and denominator
    // in order to find pairs of similar values.

    // We begin by ordering the terms in descending order
    numerators.sort_unstable();
    numerators = numerators.into_iter().rev().collect();

    denominators.sort_unstable();
    denominators = denominators.into_iter().rev().collect();

    // These will keep track of which terms remain unpaired
    let mut unpaired_numerators = vec![true; numerators.len()];
    let mut unpaired_denominators = vec![true; denominators.len()];

    // Then we loop through the numerators (largest first)
    let mut candidate_pairs: Vec<(usize, usize, i64)> = Vec::new();
    for (i, n) in numerators.iter().enumerate() {
        match denominators
            .iter()
            .enumerate()
            .filter(|(j, d)| unpaired_denominators[*j])
            .next()
        {
            // and if we can find a denominator (largest first)
            Some((j, d)) => {
                unpaired_numerators[i] = false;
                unpaired_denominators[j] = false;
                // we pair them up
                candidate_pairs.push((i, j, i64::from(*n) - i64::from(*d)));
            }
            None => break,
        }
    }

    let mut result = 1.0;
    for pair in candidate_pairs {
        let num = numerators[pair.0];
        let den = denominators[pair.1];

        // Cancel each pair member against the other,
        // so 7!/5! becomes just 6*7.
        // This shrinks the largest numbers
        // that show up during the computation.
        result *= if pair.2 >= 0 {
            ((den + 1)..=num).map(f64::from).product()
        } else {
            1.0 / ((num + 1)..=den).map(f64::from).product::<f64>()
        };
    }

    // Then we multiply in the factorials of the remaining numerators
    result *= numerators
        .iter()
        .enumerate()
        .filter_map(|(i, e)| {
            if unpaired_numerators[i] {
                Some(factorial((*e).into()))
            } else {
                None
            }
        })
        .product::<f64>();

    // and divide away the factorials of the remaining denominators
    result /= denominators
        .iter()
        .enumerate()
        .filter_map(|(i, e)| {
            if unpaired_denominators[i] {
                Some(factorial((*e).into()))
            } else {
                None
            }
        })
        .product::<f64>();

    result
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::truths::return_3j_truths;
    use std::f64::consts::PI;

    //We allow floating point errors on the scale of TOL.
    //The answers we compare against are exact,
    //but we use floats to compute the answers, so we will have
    //errors due to the formulas for 3j symbols containing alternating
    //products and divisions by large factorials.
    const TOL: f64 = 100.0 * f64::EPSILON;

    #[test]
    fn test_ratio_of_factorials() {
        assert_relative_eq!(ratio_of_factorials(vec![5, 2, 7], vec![3, 6]), 280.0);
        assert_relative_eq!(
            ratio_of_factorials(vec![1, 2, 3], vec![0, 5, 8]),
            1.0 / 403200.0
        );
        assert_relative_eq!(ratio_of_factorials(vec![200], vec![197]), 7880400.0);
        assert_relative_eq!(
            ratio_of_factorials(vec![1000000], vec![999999, 8]),
            3125.0 / 126.0
        );
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
                for j3 in ((j1 as i32) - (j2 as i32)).abs() as u32..=j1 + j2 {
                    for m1 in -(j1 as i32)..=j1 as i32 {
                        for m2 in -(j2 as i32)..=j2 as i32 {
                            if (m1 + m2).abs() as u32 > j3 {
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
