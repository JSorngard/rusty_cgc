use std::time::Instant;
fn main() {
    let now = Instant::now();
    let mut acc: f64 = 0.0;
    for _ in 0..1000 {
        acc += clebsch_gordan_coefficient(5, 5, 0, 1, -1);
    }
    let elapsed = now.elapsed();
    println!(
        "Took {:.2?} to get CGC(5,5,0,1,-1) = {}",
        elapsed,
        acc / 1000.0
    );
}

///Returns the value of the Clebsch-Gordan coefficient for
///the given inputs. Returns 0 if the arguments are invalid.
fn clebsch_gordan_coefficient(j1: i32, j2: i32, j3: i32, m1: i32, m2: i32) -> f64 {
    //Normal Fortran rules: variables beginning with
    //i,j,...,n are i32 and everything else is f32
    //Original Fortran code says: IMPLICIT REAL*8(A-H,O-Z)
    //=> all variables that begin with A-H or O-Z must be f64
    //This code is simply ported Fortran code,
    //as such it is not ideomatic rust.

    let m3 = m1 + m2;
    if j1 < i32::abs(m1) || j2 < i32::abs(m2) || j3 < i32::abs(m3) {
        return 0.0;
    }

    if j3 > (j1 + j2) || j3 < i32::abs(j1 - j2) {
        return 0.0;
    }

    if i32::abs(m1) + i32::abs(m2) == 0 && u32_bit_at((j1 + j2 + j3) as u32, 0) {
        return 0.0;
    }

    //Precompute table of factorials
    /*const FACTS: usize = 99;
    let mut fact: Vec<f64> = vec![1.0; FACTS];
    for i in 1..FACTS {
        fact[i] = i as f64 * fact[i - 1];
    }*/

    let ia1: i32 = j3 + j2 - j1;
    let ia2: i32 = j3 + m3;
    let ia3: i32 = j2 + m3 - j1;

    let ni: i32 = if ia3 > 0 { ia3 } else { 0 };
    let nm: i32 = if ia2 <= ia1 { ia2 } else { ia1 };

    let l1: i32 = j3 + j1 - j2;
    let l2: i32 = j1 + j2 + j3 + 1;
    let l3: i32 = j1 + j2 - j3;
    let l4: i32 = j1 - m1;
    let l5: i32 = j2 - m2;
    let l6: i32 = j2 + m2;
    let l7: i32 = j3 - m3;
    let l8: i32 = j1 + m1;
    let l9: i32 = 2 * j3 + 1;

    let cc: f64 = f64::sqrt(
        l9 as f64 * factorial(l1.try_into().unwrap()) as f64
            / factorial(l2.try_into().unwrap()) as f64
            * factorial(ia1.try_into().unwrap()) as f64
            * factorial(l3.try_into().unwrap()) as f64
            / factorial(l4.try_into().unwrap()) as f64
            / factorial(l5.try_into().unwrap()) as f64
            * factorial(ia2.try_into().unwrap()) as f64
            / factorial(l6.try_into().unwrap()) as f64
            * factorial(l7.try_into().unwrap()) as f64
            / factorial(l8.try_into().unwrap()) as f64,
    );

    let mut ip1: i32 = j2 + j3 + m1 - ni;
    let b1: f64 = factorial(ip1.try_into().unwrap()) as f64;
    let mut ip2: i32 = j1 - m1 + ni;
    let b2: f64 = factorial(ip2.try_into().unwrap()) as f64;
    ip2 += 1;
    let d1: f64 = factorial(ni.try_into().unwrap()) as f64;
    let mut ir1 = ni + 1;
    let mut ir2 = j3 - j1 + j2 - ni;
    let d2 = factorial(ir2.try_into().unwrap()) as f64;
    let mut ir3 = j3 + m3 - ni;
    let d3 = factorial(ir3.try_into().unwrap()) as f64;
    let mut ir4 = j1 - j2 - m3 + ni;
    let d4 = factorial(ir4.try_into().unwrap()) as f64;
    ir4 += 1;
    let mut fac: f64 = 1.0;
    if u32_bit_at((ni + j2 + m2) as u32, 0) {
        fac *= -1.0;
    }
    let mut s1: f64 = b1 / d2 * b2 / (d1 * d3 * d4) * fac;
    let n = nm - ni;

    if n != 0 {
        let mut fa = s1;
        for _ in 1..=n {
            fa =
                -fa * ip2 as f64 * ir2 as f64 / ip1 as f64 * ir3 as f64 / (ir1 as f64 * ir4 as f64);
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
    //guard against floating point errors making a zero result non-zero
    if f64::abs(res) < 1e-14 {
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

///Returns the factorial of the input number
fn factorial(n: u64) -> u64 {
    match n {
        0 => 1,
        1.. => (1..=n).product(),
    }
}
