use num::complex::Complex;
use std::f64::consts::PI;

///Returns the value of the Wigner 3j symbol for the given integer inputs. Returns 0.0
///if the arguments are invalid. The first three inputs are the angular momentum quantum
///numbers, while the last three are the magnetic quantum numbers.
pub fn wigner_3j(j1: u32, j2: u32, j3: u32, m1: i32, m2: i32, m3: i32) -> Result<f64, String> {
    match is_unphysical(j1, j2, j3, m1, m2, -m3) {
        Some(e) => return Err(e),
        None => (),
    }

    let (j1, j2, j3, m1, m2, m3, mut sign) = reorder3j(j1, j2, j3, m1, m2, m3, 1.0);

    if (i64::from(j1) - i64::from(j2) - i64::from(m3)) % 2 != 0 {
        sign *= -1.0
    }

    Ok(sign
        * match clebsch_gordan(j1, j2, j3, m1, m2, -m3) {
            Ok(x) => x,
            Err(e) => return Err(e),
        }
        / (2.0 * f64::from(j3) + 1.0).sqrt())
}

///Reorder j1/m1, j2/m2, j3/m3 such that j1 >= j2 >= j3 and m1 >= 0 or m1 == 0 && m2 >= 0
fn reorder3j(
    j1: u32,
    j2: u32,
    j3: u32,
    m1: i32,
    m2: i32,
    m3: i32,
    mut sign: f64,
) -> (u32, u32, u32, i32, i32, i32, f64) {
    //An odd permutation of the columns or a
    //sign change of the m-quantum values (time reversal)
    //give a phase factor of (-1)^(j1+j2+j3).
    //If we assume that this phase factor is -1, we are only wrong if j1+j2+j3 is even,
    //which we correct for at the end.
    if j1 < j2 {
        return reorder3j(j2, j1, j3, m2, m1, m3, -sign);
    } else if j2 < j3 {
        return reorder3j(j1, j3, j2, m1, m3, m2, -sign);
    } else if m1 < 0 || (m1 == 0 && m2 < 0) {
        return reorder3j(j1, j2, j3, -m1, -m2, -m3, -sign);
    } else {
        //Sign doesn't matter if total J = j1 + j2 + j3 is even
        if (j1 + j2 + j3) % 2 == 0 {
            sign = 1.0;
        }
        return (j1, j2, j3, m1, m2, m3, sign);
    }
}

pub fn wigner_6j(j1: u32, j2: u32, j3: u32, j4: u32, j5: u32, j6: u32) -> f64 {
    if !is_triad(j1, j2, j3)
        || !is_triad(j1, j5, j6)
        || !is_triad(j4, j2, j6)
        || !is_triad(j4, j5, j3)
    {
        return 0.0;
    }

    //The arguments to the factorials should always be positive
    let fac = delta(j2, j4, j6) * delta(j2, j1, j3) / factorial((j2 + j4 - j6).into())
        * delta(j6, j5, j1)
        / factorial((j6 + j1 - j5).into())
        * factorial((j2 + j4 + j6 + 1).into())
        / factorial((j6 + j5 - j1).into())
        * delta(j4, j5, j3)
        / factorial((j2 + j3 - j1).into())
        * factorial((j4 + j5 + j3 + 1).into())
        / factorial((j1 + j3 - j2).into())
        / factorial((j4 + j5 - j3).into())
        * phase(j4 + j6 + j1 + j3);
    let mut sum: f64 = 0.0;
    for z in 0..=u32::min(u32::min(2 * j4, j4 + j6 - j2), j4 + j3 - j5) {
        sum += factorial((2 * j4 - z).into())
            * factorial((j4 + j6 + j3 - z - j1).into())
            * factorial((j4 + j6 + j1 + j3 + 1 - z).into())
            / factorial(z.into())
            / factorial((j4 + j6 - z - j2).into())
            / factorial((j4 + j3 - z - j5).into())
            / factorial((j2 + j4 + j6 + 1 - z).into())
            / factorial((j4 + j5 + j3 + 1 - z).into())
            * phase(z);
    }
    sum * fac
}

pub fn wigner_9j(
    j11: u32,
    j21: u32,
    j31: u32,
    j12: u32,
    j22: u32,
    j32: u32,
    j13: u32,
    j23: u32,
    j33: u32,
) -> f64 {
    println!("In function");
    //Check that all rows are triads
    if !is_triad(j11, j21, j31) || !is_triad(j12, j22, j32) || !is_triad(j13, j23, j33) {
        return 0.0;
    }

    //Check that all columns are triads
    if !is_triad(j11, j12, j13) || !is_triad(j21, j22, j23) || !is_triad(j31, j32, j33) {
        return 0.0;
    }

    println!("Passed triad checks");

    let prefactor = phase(j13 + j23 - j33) * nabla(j21, j11, j31) / nabla(j21, j22, j23)
        * nabla(j12, j22, j32)
        / nabla(j12, j11, j13)
        * nabla(j33, j31, j32)
        / nabla(j33, j13, j23);

    println!("prefactor is {}", prefactor);
    let mut sum: f64 = 0.0;
    for x in 0..=*[2 * j33, j22 + j23 - j21, j13 + j23 - j33]
        .iter()
        .min()
        .unwrap()
    //array is never empty
    {
        println!("x = {}", x);
        for y in 0..=j31 + j33 - j32 {
            println!(" y = {}", y);
            for z in u32::try_from((i64::from(j11) - i64::from(j12) - i64::from(j13)).max(0))
                .unwrap()..=(j11 + j13 - j12)
            {
                println!("  z = {}", z);
                if j12 + j33 + x + z < j11 + j23 {
                    continue;
                }
                println!("   Computing for that");
                println!("    fact({} - {})", j11 + j21, j31 + z);
                sum += phase(x + y + z) * factorial((2 * j23 - x).into())
                    / factorial(x.into())
                    / factorial((j22 + j23 - x - j21).into())
                    * factorial((j21 + j22 + x - j23).into())
                    / factorial((j13 + j23 - j33 - x).into())
                    / factorial((j21 + j32 + x + y - j23 - j12).into())
                    * factorial((j13 + j33 + x - j23).into())
                    / factorial((j12 + j33 + x + z - j11 - j23).into())
                    / factorial(y.into())
                    * factorial((j22 + j32 + y - j12).into())
                    / factorial((j12 + j22 - j32 - y).into())
                    * factorial((j31 + j32 + y - j33).into())
                    / factorial((j31 + j33 - y - j32).into())
                    / factorial((2 * j32 + 1 + y).into())
                    * factorial((j11 + j21 + j33 - y - z - j32).into())
                    / factorial(z.into())
                    / factorial((j11 + j21 - j31 - z).into())
                    * factorial((2 * j11 - z).into())
                    / factorial((j11 + j13 - z - j12).into())
                    * factorial((j12 + j13 + z - j11).into())
                    / factorial((j11 + j21 + j31 + 1 - z).into());
            }
        }
    }
    println!("sum is {}", sum);
    prefactor * sum
}

///Returns the value of the Racah W coefficient.
pub fn racah_w(j1: u32, j2: u32, j: u32, j3: u32, j12: u32, j23: u32) -> f64 {
    phase(j1 + j2 + j3 + j) * wigner_6j(j1, j2, j12, j3, j, j23)
}

///Returns the Gaunt coefficient for the input angular momenta.
///The Gaunt coefficient is defined as the integral over three spherical harmonics.
pub fn gaunt(l1: u32, l2: u32, l3: u32, m1: i32, m2: i32, m3: i32) -> Result<f64, String> {
    match is_unphysical(l1, l2, l3, m1, m2, -m3) {
        Some(e) => return Err(e),
        None => (),
    }

    Ok(f64::sqrt(
        (2.0 * (l1 as f64) + 1.0) * (2.0 * (l2 as f64) + 1.0) * (2.0 * (l3 as f64) + 1.0)
            / (4.0 * PI),
    ) * wigner_3j(l1, l2, l3, 0, 0, 0).unwrap()//Have already checked for unphysicality.
        * wigner_3j(l1, l2, l3, m1, m2, m3).unwrap())
}

///Returns the value of the triangle coefficient \Delta(abc) used in the computation of 6j and 9j symbols.
///The inputs are debug_asserted to be triangular.
fn delta(a: u32, b: u32, c: u32) -> f64 {
    debug_assert!(a + c >= b && a + b >= c && b + c >= a);
    (factorial((a + c - b).into()) * factorial((a + b - c).into())
        / factorial((a + c + b + 1).into())
        * factorial((c + b - a).into()))
    .sqrt()
}

///Returns the value of the \nabla(abc) triangle coefficient used in the computation of 9j symbols.
///The inputs are debug_asserted to be triangular.
fn nabla(a: u32, b: u32, c: u32) -> f64 {
    debug_assert!(a + c >= b && a + b >= c && b + c >= a);
    (factorial((a + c - b).into()) * factorial((a + b - c).into()) / factorial((b + c - a).into())
        * factorial((a + b + c + 1).into()))
    .sqrt()
}

///Returns the value of the small Wigner d-matrix in the z-y-z convention.
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

///Returns the value of the Wigner D-matrix in the z-y-z convention.
pub fn wigner_d(
    j: u32,
    mp: i32,
    m: i32,
    alpha: f64,
    beta: f64,
    gamma: f64,
) -> Result<Complex<f64>, String> {
    let d = match wigner_small_d(j, mp, m, beta) {
        Ok(x) => x,
        Err(x) => return Err(x),
    };
    Ok((-1.0 * Complex::<f64>::i() * (mp as f64) * alpha).exp()
        * d
        * (-1.0 * Complex::<f64>::i() * (m as f64) * gamma).exp())
}

///Returns the value of the Clebsch-Gordan coefficient for
///the given integer inputs. Returns 0.0 if the arguments are invalid.
///The first three inputs are the angular momentum quantum numbers,
///while the last three are the magnetic quantum numbers.
pub fn clebsch_gordan(j1: u32, j2: u32, j3: u32, m1: i32, m2: i32, m3: i32) -> Result<f64, String> {
    //This code is simply ported Fortran code,
    //as such it is not completely idiomatic rust.

    match is_unphysical(j1, j2, j3, m1, m2, m3) {
        Some(e) => return Err(e),
        None => (),
    }

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

    let cc = f64::sqrt(
        //All inputs to the factorials will be >= 0, so casting to u64 loses no sign information
        f64::from(2 * j3 + 1) * factorial((j3 + j1 - j2).into())
            / factorial((j1 + j2 + j3 + 1).into())
            * factorial(ia1.into())
            * factorial((j1 + j2 - j3).into())
            / factorial(j_plus_m(j1, -m1).into())
            / factorial(j_plus_m(j2, -m2).into())
            * factorial(ia2.into())
            / factorial(j_plus_m(j2, m2).into())
            * factorial(j_plus_m(j3, -m3).into())
            / factorial(j_plus_m(j1, m1).into()),
    );

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
            //Same here: all inputs to the factorials will be >= 0, so casting to u64 loses no sign information
    let mut s1 = phase(ni + j_plus_m(j2, m2)) * factorial(ip1.into()) / factorial(ir2.into())
        * factorial((ip2 - 1).into())
        / (factorial(ni.into()) * factorial(ir3.into()) * factorial((ir4 - 1).into()));

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

///Returns whether the given triplet of angular momenta form a triad.
fn is_triad(j1: u32, j2: u32, j3: u32) -> bool {
    j3 >= abs_diff(j1, j2) && j3 <= j1 + j2
}

///Returns the result of adding an angular momentum to its projection.
///Debug_asserts that |m| <= j.
fn j_plus_m(j: u32, m: i32) -> u32 {
    debug_assert!(u32::try_from(m.abs()).unwrap() <= j);
    if m >= 0 {
        j + u32::try_from(m).unwrap()
    } else {
        j - u32::try_from(-m).unwrap()
    }
}

///Returns the absolute value of the difference of the unsigned integers x and y.
fn abs_diff(x: u32, y: u32) -> u32 {
    if x > y {
        x - y
    } else {
        y - x
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
    if (m1.abs() as u32) > j1 && (m2.abs() as u32) > j2 && (m3.abs() as u32) > j3 {
        Some("|m| is larger than its corresponding j".to_owned())
    } else if !is_triad(j1, j2, j3) {
        Some("j1, j2, and j3 do not fulfill the triangle condition".to_owned())
    } else if m1 + m2 != m3 {
        Some("m1 + m2 do not equal m3".to_owned())
    } else {
        None
    }
}
