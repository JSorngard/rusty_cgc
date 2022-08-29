pub mod wigner;

#[cfg(test)]
#[macro_use]
extern crate approx;

#[cfg(test)]
mod truths;

#[cfg(test)]
mod tests {
    use crate::truths::return_3j_truths;
    use crate::wigner::*;
    use std::f64::consts::PI;

    //We allow floating point errors on the scale of TOL.
    //The answers we compare against are exact,
    //but we use floats to compute the answers, so we will have
    //errors due to the formulas for 3j symbols containing alternating
    //products and divisions by large factorials.
    const TOL: f64 = 100.0 * f64::EPSILON;

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
            8129.0 / 3136350672.0
        );
        assert_relative_eq!(
            wigner_9j(10, 11, 12, 11, 12, 13, 12, 13, 14).unwrap(),
            1.0 / 31648909588721100.0
        );
        assert_relative_eq!(
            wigner_9j(11, 12, 13, 12, 13, 14, 13, 14, 15).unwrap(),
            -26976629723.0 / 599200846745523750.0
        );
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
