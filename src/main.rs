#![deny(clippy::all)]

fn main() {
    const TRIALS: u32 = 100;
    let now = std::time::Instant::now();
    let mut acc: f64 = 0.0;
    for _ in 0..TRIALS {
        acc += wigner_3j(5, 5, 0, 1, -1, 0);
    }
    let elapsed = now.elapsed();
    println!(
        "Took {:.2?} to get wigner_3j(5,5,0,1,-1,0) = {} {} times.",
        elapsed,
        acc / f64::from(TRIALS),
        TRIALS,
    );
    println!(
        "This given an average speed of {:.2?} per function call",
        elapsed / TRIALS
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
    //The type of all the l<1-8> variables have been changed
    //from i32 to u64 and l9 changed to f64 to reduce type casting

    if m3 != m1 + m2 {
        return 0.0;
    }

    if j1 < i32::abs(m1) || j2 < i32::abs(m2) || j3 < i32::abs(m3) {
        return 0.0;
    }

    if j3 > (j1 + j2) || j3 < i32::abs(j1 - j2) {
        return 0.0;
    }

    if i32::abs(m1) + i32::abs(m2) == 0 && u32_bit_at((j1 + j2 + j3) as u32, 0) {
        return 0.0;
    }

    let ia1 = j3 + j2 - j1;
    let ia2 = j3 + m3;
    let ia3 = j2 + m3 - j1;

    let ni = if ia3 > 0 { ia3 } else { 0 };
    let nm = if ia2 <= ia1 { ia2 } else { ia1 };

    let l1 = (j3 + j1 - j2) as u64;
    let l2 = (j1 + j2 + j3 + 1) as u64;
    let l3 = (j1 + j2 - j3) as u64;
    let l4 = (j1 - m1) as u64;
    let l5 = (j2 - m2) as u64;
    let l6 = (j2 + m2) as u64;
    let l7 = (j3 - m3) as u64;
    let l8 = (j1 + m1) as u64;
    let l9 = f64::from(2 * j3 + 1);

    let cc: f64 = f64::sqrt(
        l9 * factorial(l1) / factorial(l2) * factorial(ia1.try_into().unwrap()) * factorial(l3)
            / factorial(l4)
            / factorial(l5)
            * factorial(ia2.try_into().unwrap())
            / factorial(l6)
            * factorial(l7)
            / factorial(l8),
    );

    let mut ip1 = j2 + j3 + m1 - ni;
    let b1 = factorial(ip1.try_into().unwrap());
    let mut ip2 = j1 - m1 + ni;
    let b2 = factorial(ip2.try_into().unwrap());
    ip2 += 1;
    let d1 = factorial(ni.try_into().unwrap());
    let mut ir1 = ni + 1;
    let mut ir2 = j3 - j1 + j2 - ni;
    let d2 = factorial(ir2.try_into().unwrap());
    let mut ir3 = j3 + m3 - ni;
    let d3 = factorial(ir3.try_into().unwrap());
    let mut ir4 = j1 - j2 - m3 + ni;
    let d4 = factorial(ir4.try_into().unwrap());
    ir4 += 1;
    let mut fac: f64 = 1.0;
    if u32_bit_at((ni + j2 + m2) as u32, 0) {
        fac *= -1.0;
    }
    let mut s1 = b1 / d2 * b2 / (d1 * d3 * d4) * fac;
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

///Gets the bit at position n. Bits are numbered from 0 (least significant) to 31 (most significant).
fn u32_bit_at(input: u32, n: u8) -> bool {
    if n < 32 {
        input & (1 << n) != 0
    } else {
        false
    }
}

///Returns the factorial of the input integer as a float
fn factorial(n: u64) -> f64 {
    (match n {
        0 => 1,
        1.. => (1..=n).product(),
    }) as f64
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
}
