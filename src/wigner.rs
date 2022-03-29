use num::complex::Complex;

///Returns the value of the Wigner 3j symbol for the given integer inputs.
///Returns 0.0 if the arguments are invalid. The first three inputs are the
///angular momentum quantum numbers, while the last three are the magnetic
///quantum numbers belonging to the first three angular momentum quantum numbers.
pub fn wigner_3j(j1: u32, j2: u32, j3: u32, m1: i32, m2: i32, m3: i32) -> f64 {
    let (j1, j2, j3, m1, m2, m3, mut sign) = reorder3j(j1, j2, j3, m1, m2, m3, 1.0);

    if (i64::from(j1) - i64::from(j2) - i64::from(m3)) % 2 != 0 {
        sign *= -1.0
    }

    sign * clebsch_gordan(j1, j2, j3, m1, m2, -m3) / (2.0 * f64::from(j3) + 1.0).sqrt()
}

///Reorder j1/m1, j2/m2, j3/m3 such that j1 >= j2 >= j3 and m1 >= 0 or m1 == 0 && m2 >= 0
#[allow(dead_code)]
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

pub fn wigner_6j(j1: i32, j2: i32, j3: i32, j4: i32, j5: i32, j6: i32) -> f64 {
    if j4 < 0 {
        return 0.0;
    }

    if !is_triad(j1, j2, j3)
        || !is_triad(j1, j5, j6)
        || !is_triad(j4, j2, j6)
        || !is_triad(j4, j5, j3)
    {
        return 0.0;
    }

    let fac = delta(j2, j4, j6) * delta(j2, j1, j3) / factorial((j2 + j4 - j6).try_into().unwrap())
        * delta(j6, j5, j1)
        / factorial((j6 - j5 + j1).try_into().unwrap())
        * factorial((j2 + j4 + j6 + 1).try_into().unwrap())
        / factorial((j6 + j5 - j1).try_into().unwrap())
        * delta(j4, j5, j3)
        / factorial((j2 - j1 + j3).try_into().unwrap())
        * factorial((j4 + j5 + j3 + 1).try_into().unwrap())
        / factorial((-j2 + j1 + j3).try_into().unwrap())
        / factorial((j4 + j5 - j3).try_into().unwrap())
        * phase(j4 + j6 + j1 + j3);
    let mut sum: f64 = 0.0;
    for z in 0..=i32::min(i32::min(2 * j4, -j2 + j4 + j6), j4 - j5 + j3) {
        sum += factorial((2 * j4 - z).try_into().unwrap())
            * factorial((j4 + j6 - j1 + j3 - z).try_into().unwrap())
            * factorial((j4 + j6 + j1 + j3 + 1 - z).try_into().unwrap())
            / factorial(z.try_into().unwrap())
            / factorial((-j2 + j4 + j6 - z).try_into().unwrap())
            / factorial((j4 - j5 + j3 - z).try_into().unwrap())
            / factorial((j2 + j4 + j6 + 1 - z).try_into().unwrap())
            / factorial((j4 + j5 + j3 + 1 - z).try_into().unwrap())
            * phase(z);
    }
    sum * fac
}

pub fn wigner_9j(
    j11: i32,
    j21: i32,
    j31: i32,
    j12: i32,
    j22: i32,
    j32: i32,
    j13: i32,
    j23: i32,
    j33: i32,
) -> f64 {
    //Check that all rows are triads
    if !is_triad(j11, j21, j31) || !is_triad(j12, j22, j32) || !is_triad(j13, j23, j33) {
        return 0.0;
    }

    //Check that all columns are triads
    if !is_triad(j11, j12, j13) || !is_triad(j21, j22, j23) || !is_triad(j31, j32, j33) {
        return 0.0;
    }

    let prefactor = phase(j13 + j23 - j33) * nabla(j21, j11, j31)
        / nabla(j21, j22, j23)
        * nabla(j12, j22, j32)
        / nabla(j12, j11, j13)
        * nabla(j33, j31, j32)
        / nabla(j33, j13, j23);
    let mut sum: f64 = 0.0;
    for x in 0..=i32::min(i32::min(2 * j33, j22 - j21 + j23), j13 + j23 - j33) {
        for y in 0..=j31 - j32 + j33 {
            for z in i32::max(j11 - j12 - j13, 0)..=(j11 - j12 + j13) {
                if -j11 + j12 + j33 - j23 + x + z < 0 {
                    continue;
                }
                // let numerator: f64 = vec![
                //     2 * j23 - x,
                //     j21 + j22 - j23 + x,
                //     j13 - j23 + j33 + x,
                //     j22 - j12 + j32 + y,
                //     j31 + j32 - j33 + y,
                //     j11 + j21 - j32 + j33 - y - z,
                //     2 * j11 - z,
                //     j12 - j11 + j13 + z,
                // ]
                // .into_iter()
                // .map(|x| factorial(x.try_into().unwrap()))
                // .product();
                // let denominator: f64 = vec![
                //     x,
                //     j22 - j21 + j23 - x,
                //     j13 + j23 - j33 - x,
                //     j21 - j12 + j32 - j23 + x + y,
                //     j12 - j11 - j23 + j33 + x + z,
                //     y,
                //     j12 + j22 - j32 - y,
                //     j31 - j32 + j33 - y,
                //     2 * j32 + 1 + y,
                //     z,
                //     j11 + j21 - j31 - z,
                //     j11 - j12 + j13 - z,
                //     j11 + j21 + j31 + 1 - z,
                // ]
                // .into_iter()
                // .map(|x| factorial(x.try_into().unwrap()))
                // .product();
                // sum += if (x + y + z) % 2 == 1 { -1.0 } else { 1.0 } * numerator / denominator
                /*println!("x={}, y={}, z={}", x, y, z);
                println!("gives j11 + j21 - j31 - z = {}", j11 + j21 - j31 - z);*/

                sum += phase(x + y + z)
                    * factorial((2 * j23 - x).try_into().unwrap())
                    / factorial(x.try_into().unwrap())
                    / factorial((j22 - j21 + j23 - x).try_into().unwrap())
                    * factorial((j21 + j22 - j23 + x).try_into().unwrap())
                    / factorial((j13 + j23 - j33 - x).try_into().unwrap())
                    / factorial((j21 - j12 + j32 - j23 + x + y).try_into().unwrap())
                    * factorial((j13 - j23 + j33 + x).try_into().unwrap())
                    / factorial((j12 - j11 - j23 + j33 + x + z).try_into().unwrap())
                    / factorial(y.try_into().unwrap())
                    * factorial((j22 - j12 + j32 + y).try_into().unwrap())
                    / factorial((j12 + j22 - j32 - y).try_into().unwrap())
                    * factorial((j31 + j32 - j33 + y).try_into().unwrap())
                    / factorial((j31 - j32 + j33 - y).try_into().unwrap())
                    / factorial((2 * j32 + 1 + y).try_into().unwrap())
                    * factorial((j11 + j21 - j32 + j33 - y - z).try_into().unwrap())
                    / factorial(z.try_into().unwrap())
                    / factorial((j11 + j21 - j31 - z).try_into().unwrap())
                    * factorial((2 * j11 - z).try_into().unwrap())
                    / factorial((j11 - j12 + j13 - z).try_into().unwrap())
                    * factorial((j12 - j11 + j13 + z).try_into().unwrap())
                    / factorial((j11 + j21 + j31 + 1 - z).try_into().unwrap());
            }
        }
    }
    prefactor * sum
}

fn delta(a: i32, b: i32, c: i32) -> f64 {
    (factorial((a + c - b).try_into().unwrap()) * factorial((a - c + b).try_into().unwrap())
        / factorial((a + c + b + 1).try_into().unwrap())
        * factorial((-a + c + b).try_into().unwrap()))
    .sqrt()
}

fn nabla(a: i32, b: i32, c: i32) -> f64 {
    (factorial((a - b + c).try_into().unwrap()) * factorial((a + b - c).try_into().unwrap())
        / factorial((b + c - a).try_into().unwrap())
        * factorial((a + b + c + 1).try_into().unwrap()))
    .sqrt()
}

///Returns the value of the small Wigner d-matrix in the z-y-z convention.
///If the input is not valid it returns an error instead.
pub fn wigner_small_d(j: i32, mp: i32, m: i32, beta: f64) -> Result<f64, String> {
    if mp.abs() > j || m.abs() > j {
        return Err("m-values can not be larger than j".to_owned());
    }

    let prefactor = f64::sqrt(
        factorial((j + mp).try_into().unwrap())
            * factorial((j - mp).try_into().unwrap())
            * factorial((j + m).try_into().unwrap())
            * factorial((j - m).try_into().unwrap()),
    );
    let mut sum: f64 = 0.0;
    for s in i32::max(0, m - mp)..=i32::min(j + m, j - mp) {
        sum += f64::powf((beta / 2.0).cos(), (2 * j + m - mp - 2 * s).into())
            * f64::powf((beta / 2.0).sin(), (mp - m + 2 * s).into())
            / (factorial((j + m - s).try_into().unwrap())
                * factorial(s.try_into().unwrap())
                * factorial((mp - m + s).try_into().unwrap())
                * factorial((j - mp - s).try_into().unwrap()))
            * phase(mp - m + s);
    }
    Ok(sum * prefactor)
}

///Returns the value of the Wigner D-matrix in the z-y-z convention.
///If the input is not valid it returns an error instead.
pub fn wigner_d(
    j: i32,
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
///while the last three are the magnetic quantum numbers belonging to the first
///three angular momentum quantum numbers.
pub fn clebsch_gordan(uj1: u32, uj2: u32, uj3: u32, m1: i32, m2: i32, m3: i32) -> f64 {
    //Normal Fortran rules: variables beginning with
    //i,j,...,n are i32 and everything else is f32
    //Original Fortran code says: IMPLICIT REAL*8(A-H,O-Z)
    //=> all variables that begin with A-H or O-Z must be f64
    //This code is simply ported Fortran code,
    //as such it is not idiomatic rust.

    let j1: i32 = uj1.try_into().unwrap();
    let j2: i32 = uj2.try_into().unwrap();
    let j3: i32 = uj3.try_into().unwrap();

    if m3 != m1 + m2 {
        return 0.0;
    }

    if j1 < m1.abs() || j2 < m2.abs() || j3 < m3.abs() {
        return 0.0;
    }

    if !is_triad(j1, j2, j3) {
        return 0.0;
    }

    if m1.abs() + m2.abs() == 0 && (j1 + j2 + j3) % 2 == 1 {
        return 0.0;
    }

    let ia1 = j3 + j2 - j1;
    let ia2 = j3 + m3;
    let ia3 = j2 + m3 - j1;

    let ni = if ia3 > 0 { ia3 } else { 0 };
    let nm = if ia2 <= ia1 { ia2 } else { ia1 };

    let cc = f64::sqrt(
        //All inputs to the factorials will be >= 0, so casting to u64 loses no sign information
        f64::from(2 * j3 + 1) * factorial((j3 + j1 - j2) as u64)
            / factorial((j1 + j2 + j3 + 1) as u64)
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
        / (factorial(ni as u64) * factorial(ir3 as u64) * factorial((ir4 - 1) as u64));
    if (ni + j2 + m2) % 2 == 1 {
        s1 = -s1;
    }
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
        0.0
    } else {
        res
    }
}

///Returns whether the given triplet of angular momenta form a triad
fn is_triad(j1: i32, j2: i32, j3: i32) -> bool {
    j3 >= (j1 - j2).abs() && j3 <= j1 + j2
}

#[allow(dead_code)]
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

///Returns a f64 with a value of 1.0 if the input is even, and -1.0 if it is odd
fn phase(x: i32) -> f64 {
    if x % 2 == 0 {
        1.0
    } else {
        -1.0
    }
}