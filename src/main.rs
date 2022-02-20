#![deny(clippy::all)]

fn main() {
    const MAXJ: i32 = 5;
    let start_time = std::time::Instant::now();
    let mut acc: f64 = 0.0;
    for j1 in 0..=MAXJ {
        for m1 in -j1..=j1 {
            for j2 in 0..=MAXJ {
                for m2 in -j2..=j2 {
                    for j3 in 0..=MAXJ {
                        for m3 in -j3..=j3 {
                            acc += wigner_3j(j1, j2, j3, m1, m2, m3);
                        }
                    }
                }
            }
        }
    }
    let elapsed = start_time.elapsed();
    let num_symbols = i32::pow(MAXJ * (2 * MAXJ + 1), 3) as u32;
    println!(
        "Took {:.2?} to compute all {} 3j-symbols with js of at most {}. Their sum is {}",
        elapsed, num_symbols, MAXJ, acc,
    );
    println!(
        "This gives an average speed of {:.2?} per function call",
        elapsed / num_symbols
    );
}

///Returns the value of the Wigner 3j symbol for the given integer inputs.
///Returns 0.0 if the arguments are invalid. The first three inputs are the
///angular momentum quantum numbers, while the last three are the magnetic
///quantum numbers belonging to the first three angular momentum quantum numbers.
fn wigner_3j(j1: i32, j2: i32, j3: i32, m1: i32, m2: i32, m3: i32) -> f64 {
    let sign: f64 = if (j1 - j2 - m3) % 2 == 0 { 1.0 } else { -1.0 };
    sign * clebsch_gordan_coefficient(j1, j2, j3, m1, m2, -m3)
        / f64::sqrt(2.0 * f64::from(j3) + 1.0)
}

///Returns the value of the Clebsch-Gordan coefficient for
///the given integer inputs. Returns 0.0 if the arguments are invalid.
///The first three inputs are the angular momentum quantum numbers,
///while the last three are the magnetic quantum numbers belonging to the first
///three angular momentum quantum numbers.
fn clebsch_gordan_coefficient(j1: i32, j2: i32, j3: i32, m1: i32, m2: i32, m3: i32) -> f64 {
    //Normal Fortran rules: variables beginning with
    //i,j,...,n are i32 and everything else is f32
    //Original Fortran code says: IMPLICIT REAL*8(A-H,O-Z)
    //=> all variables that begin with A-H or O-Z must be f64
    //This code is simply ported Fortran code,
    //as such it is not idiomatic rust.

    if m3 != m1 + m2 {
        return 0.0;
    }

    if j1 < i32::abs(m1) || j2 < i32::abs(m2) || j3 < i32::abs(m3) {
        return 0.0;
    }

    if j3 > (j1 + j2) || j3 < i32::abs(j1 - j2) {
        return 0.0;
    }

    if i32::abs(m1) + i32::abs(m2) == 0 && (j1 + j2 + j3) % 2 == 1 {
        return 0.0;
    }

    let ia1 = j3 + j2 - j1;
    let ia2 = j3 + m3;
    let ia3 = j2 + m3 - j1;

    let ni = if ia3 > 0 { ia3 } else { 0 };
    let nm = if ia2 <= ia1 { ia2 } else { ia1 };

    let cc = f64::sqrt(
        //All inputs to the factorials will be >= 0, so casting to u64 loses no sign information
        f64::from(2 * j3 + 1) * (factorial((j3 + j1 - j2) as u64)
            / factorial((j1 + j2 + j3 + 1) as u64))
            * factorial(ia1 as u64)
            * factorial((j1 + j2 - j3) as u64)
            / factorial((j1 - m1) as u64)
            / factorial((j2 - m2) as u64)
            * factorial(ia2 as u64)
            / factorial((j2 + m2) as u64)
            * factorial((j3 - m3) as u64)
            / factorial((j1 + m1) as u64),
    );

    let mut ip1 = j2 + j3 + m1 - ni;
    let mut ip2 = j1 - m1 + ni + 1;
    let mut ir1 = ni + 1;
    let mut ir2 = j3 - j1 + j2 - ni;
    let mut ir3 = j3 + m3 - ni;
    let mut ir4 = j1 - j2 - m3 + ni + 1;
    //Same here: all inputs to the factorials will be >= 0, so casting to u64 loses no sign information
    let mut s1 = factorial(ip1 as u64) / factorial(ir2 as u64) * factorial((ip2 - 1) as u64)
        / (factorial(ni as u64) * factorial(ir3 as u64) * factorial((ir4 - 1) as u64))
        * if (ni + j2 + m2) % 2 == 1 { -1.0 } else { 1.0 };
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
    if f64::abs(res) < 1e-14 {
        //guard against floating point errors making a zero result non-zero
        0.0
    } else {
        res
    }
}

///Returns the factorial of the input integer as a float
fn factorial(n: u64) -> f64 {
    if n == 0 {
        1.0
    } else {
        let mut res = 1.0;
        for i in 1..=n {
            res *= i as f64;
        }
        res
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_bad_cgc_input() {
        assert_eq!(clebsch_gordan_coefficient(1, 1, 0, 1, 1, 0), 0.0); //m3 != m1 + m2
        assert_eq!(clebsch_gordan_coefficient(1, 2, 0, 0, 0, 0), 0.0); //|m1| > j1
        assert_eq!(clebsch_gordan_coefficient(2, 2, 5, 0, -1, -1), 0.0); //non-triangular
        assert_eq!(clebsch_gordan_coefficient(1, 1, 1, 0, 0, 0), 0.0); //all ms 0 => j3 must be even
    }

    #[test]
    fn test_good_cgc_input() {
        assert_eq!(
            clebsch_gordan_coefficient(5, 5, 0, 1, -1, 0),
            0.30151134457776363
        );
        assert_eq!(
            clebsch_gordan_coefficient(1, 2, 1, 1, -1, 0),
            0.5477225575051661
        );
        /*//Fails on the last decimal due to floating point errors
        assert_eq!(
            clebsch_gordan_coefficient(1, 1, 0, 1, -1),
            0.5773502691896258//last digit of _8... has been rounded up from _76..
        );*/
        assert_eq!(clebsch_gordan_coefficient(1, 2, 3, 1, 2, 3), 1.0);
    }

    #[test]
    fn test_bad_3j_input() {
        assert_eq!(wigner_3j(1, 1, 0, 1, 1, 0), 0.0); //m3 != m1 + m2
        assert_eq!(wigner_3j(1, 2, 0, 0, 0, 0), 0.0); //|m| > j
        assert_eq!(wigner_3j(2, 2, 5, 0, 0, 0), 0.0); //non-triangular
        assert_eq!(wigner_3j(1, 1, 1, 0, 0, 0), 0.0); //all ms = 0 => j3 must be even
    }

    #[test]
    fn test_good_3j_input() {
        /*//Fails due to floating point errors on the last decimal
        assert_eq!(wigner_3j(1, 1, 0, 0, 0, 0), -0.5773502691896258);*/
        assert_eq!(wigner_3j(1, 1, 2, -1, 1, 0), 0.18257418583505536);
        assert_eq!(wigner_3j(1, 2, 3, 1, 2, -3), 0.3779644730092272);
    }

    #[test]
    fn test_sum_of_many_3j() {
        const MAXJ: i32 = 10;
        let mut acc: f64 = 0.0;
        for j1 in 0..=MAXJ {
            for m1 in -j1..=j1 {
                for j2 in 0..=MAXJ {
                    for m2 in -j2..=j2 {
                        for j3 in 0..=MAXJ {
                            for m3 in -j3..=j3 {
                                acc += wigner_3j(j1, j2, j3, m1, m2, m3);
                            }
                        }
                    }
                }
            }
        }
        assert_eq!(acc, 48.03155288822479);
    }
}
