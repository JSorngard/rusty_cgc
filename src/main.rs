use crate::wigner::wigner_3j;
use crate::wigner::wigner_6j;
mod wigner;

fn main() {
    const MAXJ: i32 = 10;
    let start_time = std::time::Instant::now();
    let mut acc: f64 = 0.0;
    for j1 in 0..=MAXJ {
        for j2 in 0..=MAXJ {
            for j3 in i32::abs(j1 - j2)..=j1 + j2 {
                for m1 in -j1..=j1 {
                    for m2 in -j2..=j2 {
                        acc += wigner_3j(j1, j2, j3, m1, m2, -m1 - m2);
                    }
                }
            }
        }
    }
    println!("{}", wigner_6j(1, 2, 3, 4, 5, 6));
    let elapsed = start_time.elapsed();
    let num_symbols = i32::pow(MAXJ * (2 * MAXJ + 1), 3) as u32;
    println!(
        "Took {:.2?} to compute all {} 3j-symbols with j1 and j2 of at most {}. Their sum is {}",
        elapsed, num_symbols, MAXJ, acc,
    );
    println!(
        "This gives an average speed of {:.2?} per function call",
        elapsed / num_symbols
    );
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::wigner::clebsch_gordan;

    #[test]
    fn test_bad_cgc_input() {
        assert_eq!(clebsch_gordan(1, 1, 0, 1, 1, 0), 0.0); //m3 != m1 + m2
        assert_eq!(clebsch_gordan(1, 2, 0, 0, 0, 0), 0.0); //|m1| > j1
        assert_eq!(clebsch_gordan(2, 2, 5, 0, -1, -1), 0.0); //non-triangular
        assert_eq!(clebsch_gordan(1, 1, 1, 0, 0, 0), 0.0); //all ms 0 => j3 must be even
    }

    #[test]
    fn test_good_cgc_input() {
        assert_eq!(clebsch_gordan(5, 5, 0, 1, -1, 0), 0.30151134457776363);
        assert_eq!(clebsch_gordan(1, 2, 1, 1, -1, 0), 0.5477225575051661);
        assert_eq!(clebsch_gordan(1, 1, 0, 1, -1, 0), 0.5773502691896257);
        assert_eq!(clebsch_gordan(1, 2, 3, 1, 2, 3), 1.0);
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
        assert_eq!(wigner_3j(1, 1, 0, 0, 0, 0), -0.5773502691896257);
        assert_eq!(wigner_3j(1, 1, 2, -1, 1, 0), 0.18257418583505536);
        assert_eq!(wigner_3j(1, 2, 3, 1, 2, -3), 0.3779644730092272);
        assert_eq!(wigner_3j(10, 10, 9, 5, 4, -9), 0.09825449077444699);
    }

    #[test]
    fn test_sum_of_many_3j() {
        const MAXJ: i32 = 10;
        let mut acc: f64 = 0.0;
        for j1 in 0..=MAXJ {
            for j2 in 0..=MAXJ {
                for j3 in i32::abs(j1 - j2)..=j1 + j2 {
                    for m1 in -j1..=j1 {
                        for m2 in -j2..=j2 {
                            acc += wigner_3j(j1, j2, j3, m1, m2, -m1 - m2);
                        }
                    }
                }
            }
        }
        assert_eq!(acc, 57.84244636528764);
    }
}
