use std::collections::LinkedList;

use landau::dyadic::Interval;
use landau::dyadic::{ComplexDyadic, Dyadic};
use landau::dyadic::{add_complex_vec, add_vec, psi};

#[test]
fn addition_dyadic() {
    let result = Dyadic::new(1, 1) + Dyadic::new(1, 1);
    assert_eq!(result, Dyadic::new(2, 1));
}

#[test]
fn multiplication_dyadic() {
    let result = Dyadic::new(5, 2) * Dyadic::new(3, 1);
    assert_eq!(result, Dyadic::new(15, 3));
}

#[test]
fn division_with_precision_real() {
    // 3/32 ÷ 5/4 = 3/40
    let a = Dyadic {
        numerator: 3,
        exponent: -5,
    };
    let b = Dyadic {
        numerator: 5,
        exponent: -2,
    };
    let q = a.div_with_precision(b, 30);
    let (qa, qb) = (q.to_f64(), 3.0 / 40.0);
    assert!((qa - qb).abs() < 1e-9, "got {qa}, expected {qb}");
}

#[test]
fn test_div_with_precision_complex() {
    let a = ComplexDyadic {
        re: Dyadic {
            numerator: 3,
            exponent: -5,
        }, // 3/32
        im: Dyadic {
            numerator: 1,
            exponent: -4,
        }, // 1/16
    };
    let b = ComplexDyadic {
        re: Dyadic {
            numerator: 5,
            exponent: -2,
        }, // 5/4
        im: Dyadic {
            numerator: 1,
            exponent: -3,
        }, // 1/8
    };
    let q = a.div_with_precision(b, 40);
    let (approx_re, approx_im) = q.to_f64();

    let a_re = 3.0 / 32.0;
    let a_im = 1.0 / 16.0;
    let b_re = 5.0 / 4.0;
    let b_im = 1.0 / 8.0;

    let denom = b_re * b_re + b_im * b_im;
    let ref_re = (a_re * b_re + a_im * b_im) / denom;
    let ref_im = (a_im * b_re - a_re * b_im) / denom;
    // ----------------------------------------------------------

    assert!((approx_re - ref_re).abs() < 1e-9);
    assert!((approx_im - ref_im).abs() < 1e-9);
}

#[test]
fn dyadic_reduction_and_equality() {
    // 4 * 2^(-2) == 1
    assert_eq!(Dyadic::new(4, -2).reduce(), Dyadic::new(1, 0));
    // Unreduced equal values compare equal
    assert_eq!(Dyadic::new(8, -3), Dyadic::new(1, 0));
    // Negative values, different (n,e) same product
    assert_eq!(Dyadic::new(-6, 1), Dyadic::new(-3, 2));
}

#[test]
fn dyadic_approximate_bound_and_shape() {
    // Should choose an exponent within ±max and be finite
    let d = Dyadic::approximate(0.3, 6);
    assert!(d.exponent.unsigned_abs() <= 6);
    let v = d.to_f64();
    assert!(v.is_finite());
}

#[test]
fn complex_division_multiplication_inverse() {
    // (z1 / z2) * z2 ≈ z1 for z2 ≠ 0
    let z1 = ComplexDyadic::new(Dyadic::new(3, -2), Dyadic::new(5, -3));
    let z2 = ComplexDyadic::new(Dyadic::new(7, -3), Dyadic::new(-1, -2));
    let back = (z1 / z2) * z2;
    let (br, bi) = back.to_f64();
    let (r, i) = z1.to_f64();
    assert!((br - r).abs() < 1e-8 && (bi - i).abs() < 1e-8);
}

#[test]
fn complex_multiplication_matches_f64() {
    // (1/2 + i/4) * (3/4 - i/2)
    let z1 = ComplexDyadic::new(Dyadic::new(1, -1), Dyadic::new(1, -2));
    let z2 = ComplexDyadic::new(Dyadic::new(3, -2), Dyadic::new(-1, -1));
    let z = z1 * z2;

    let (zr, zi) = z.to_f64();
    let (z1r, z1i) = z1.to_f64();
    let (z2r, z2i) = z2.to_f64();
    let ref_r = z1r * z2r - z1i * z2i;
    let ref_i = z1r * z2i + z1i * z2r;

    assert!((zr - ref_r).abs() < 1e-9);
    assert!((zi - ref_i).abs() < 1e-9);
}

#[test]
fn ring_identities_small_grid() {
    let nums: [i128; 5] = [-2, -1, 0, 1, 2];
    let exps: [i32; 3] = [-2, 0, 2];

    for &an in &nums {
        for &ae in &exps {
            let a = Dyadic::new(an, ae);
            // identity + annihilator
            assert_eq!(a + Dyadic::zero(), a);
            assert_eq!(a * Dyadic::new(1, 0), a);
            assert_eq!(a * Dyadic::zero(), Dyadic::zero());

            // inverse-ish sanity: (a + a) - a = a
            assert_eq!((a + a) - a, a);

            // distributivity on a tiny grid
            for &bn in &nums {
                for &be in &exps {
                    let b = Dyadic::new(bn, be);
                    for &cn in &nums {
                        for &ce in &exps {
                            let c = Dyadic::new(cn, ce);
                            assert_eq!(a * (b + c), a * b + a * c);
                        }
                    }
                }
            }
        }
    }
}

#[test]
fn complex_abs_matches_pythagoras() {
    let zs = [
        ComplexDyadic::new(Dyadic::new(1, -1), Dyadic::new(1, -2)),
        ComplexDyadic::new(Dyadic::new(-3, 1), Dyadic::new(5, -3)),
        ComplexDyadic::new(Dyadic::new(0, 0), Dyadic::new(0, 0)),
        ComplexDyadic::new(Dyadic::new(7, -4), Dyadic::new(-2, -4)),
    ];
    for z in zs {
        let (re, im) = z.to_f64();
        let abs = z.abs();
        assert!((abs * abs - (re * re + im * im)).abs() < 1e-12);
    }
}

#[test]
fn vector_ops_handle_mismatched_lengths() {
    let a = vec![Dyadic::new(1, 0), Dyadic::new(1, -1)]; // [1, 0.5]
    let b = vec![Dyadic::new(3, 0)]; // [3]
    let sum = add_vec(&a, &b); // [4, 0.5]
    assert_eq!(sum.len(), 2);
    assert_eq!(sum[0], Dyadic::new(4, 0));
    assert_eq!(sum[1], Dyadic::new(1, -1));

    let ac = vec![
        ComplexDyadic::new(Dyadic::new(1, 0), Dyadic::new(1, -1)),
        ComplexDyadic::new(Dyadic::new(0, 0), Dyadic::new(1, -2)),
    ];
    let bc = vec![ComplexDyadic::new(Dyadic::new(2, 0), Dyadic::new(0, 0))];
    let sumc = add_complex_vec(&ac, &bc);
    assert_eq!(sumc.len(), 2);
    let (r0, i0) = sumc[0].to_f64();
    assert!((r0 - 3.0).abs() < 1e-12 && (i0 - 0.5).abs() < 1e-12);
    let (r1, i1) = sumc[1].to_f64();
    assert!((r1 - 0.0).abs() < 1e-12 && (i1 - 0.25).abs() < 1e-12);
}

#[test]
fn interval_split_and_psi_midpoint() {
    // Start with unit square [0,1] x [0,1]
    let unit = Interval::new(
        Dyadic::new(0, 0),
        Dyadic::new(1, 0),
        Dyadic::new(0, 0),
        Dyadic::new(1, 0),
    );
    // Choose quadrant 4 (upper-right): x∈[0.5,1], y∈[0.5,1]
    let mut word = LinkedList::new();
    word.push_back(4);
    let sub = psi(unit, &word);
    let mid = sub.midpoint().expect("not empty");
    let (x, y) = mid.to_f64();
    assert!(x >= 0.5 && x <= 1.0 && y >= 0.5 && y <= 1.0);
}
