// use std::intrinsics::sqrtf64;
use std::cmp::Ordering;
use std::collections::LinkedList;
//use num_complex::Complex;
use std::fmt::{self, UpperExp, write};
use std::hash::Hash;
use std::io::Empty;
use std::ops::{Add, Div, Mul, Sub};

use plotters::element::ComposedElement; //division not yet implemented

#[derive(Debug, Clone, Copy)]

// Define dyadic type
pub struct Dyadic {
    pub numerator: i128,
    pub exponent: i32,
}

impl Dyadic {
    pub fn new(num: i128, exp: i32) -> Self {
        Dyadic {
            numerator: num,
            exponent: exp,
        }
    }
    pub fn to_f64(&self) -> f64 {
        (self.numerator as f64) * (2.0f64).powi(self.exponent)
    }
    // This is used to determine if a number is zero, which need to be checked before division
    pub fn zero() -> Dyadic {
        Dyadic {
            numerator: 0,
            exponent: 0,
        }
    }

    pub fn from_i64(n: i64) -> Dyadic {
        Dyadic {
            numerator: n as i128,
            exponent: 0,
        }
    }

    pub fn reduce(self) -> Dyadic {
        let mut numerator = self.numerator;
        let mut exponent = self.exponent;
        while numerator % 2 == 0 && numerator != 0 {
            numerator /= 2;
            exponent += 1;
        }
        Dyadic::new(numerator, exponent)
    }

    pub fn powi(self, power: i32) -> Self {
        if power == 0 {
            return Dyadic::new(1, 0);
        }

        let new_numerator = if power > 0 {
            self.numerator.pow(power as u32)
        } else {
            // Negative power: numerator will be in denominator
            1 // simplified: we will adjust exponent instead
        };

        let new_exponent = self.exponent * power
            - if power < 0 {
                (self.numerator.abs().ilog2() as i32) * power.abs()
            } else {
                0
            };

        Dyadic {
            numerator: new_numerator,
            exponent: new_exponent,
        }
    }

    pub fn approximate(x: f64, max_exponent: u32) -> Dyadic {
        let mut best = Dyadic {
            numerator: 0,
            exponent: 0,
        };
        let mut min_error = f64::MAX;

        for e in -(max_exponent as i32)..=(max_exponent as i32) {
            let factor = 2f64.powi(e);
            let approx_coeff = (x / factor).round() as i128;
            let approx_val = approx_coeff as f64 * factor;
            let error = (x - approx_val).abs();

            if error < min_error {
                min_error = error;
                best = Dyadic {
                    numerator: approx_coeff,
                    exponent: e,
                };
            }
        }

        best
    }

    pub fn div_with_precision(&self, denom: Dyadic, bits: u32) -> Dyadic {
        assert!(denom.numerator != 0, "division by zero dyadic");

        // integer scaling: n_a * 2^bits
        let scaled = self.numerator << bits;

        // rounded integer division (nearest, ties to even)
        let div = {
            let num = scaled;
            let den = denom.numerator;
            let q = num / den;
            let r = num % den;
            // rounding decision
            if 2 * r.abs() > den.abs() {
                q + q.signum()
            } else if 2 * r.abs() == den.abs() && (q & 1) != 0 {
                q + q.signum()
            } else {
                q
            }
        };

        // exponent: (k_a - k_b - bits)
        Dyadic {
            numerator: div,
            exponent: self.exponent - denom.exponent - (bits as i32),
        }
        .reduce()
    }
}

impl Add for Dyadic {
    type Output = Dyadic;
    fn add(self, other: Dyadic) -> Dyadic {
        let exp_diff = self.exponent - other.exponent;

        let result_numerator;
        let result_exponent;

        if exp_diff > 0 {
            result_numerator = other.numerator + (self.numerator << exp_diff);
            result_exponent = other.exponent;
        } else {
            result_numerator = (other.numerator << -exp_diff) + self.numerator;
            result_exponent = self.exponent;
        }
        Dyadic::new(result_numerator, result_exponent).reduce()
    }
}

impl Sub for Dyadic {
    type Output = Dyadic;

    fn sub(self, other: Dyadic) -> Dyadic {
        let exp_diff = self.exponent - other.exponent;

        if exp_diff > 0 {
            Dyadic::new(
                -other.numerator + (self.numerator << exp_diff),
                other.exponent,
            )
            .reduce()
        } else {
            Dyadic::new(
                -(other.numerator << -exp_diff) + self.numerator,
                self.exponent,
            )
            .reduce()
        }
    }
}

impl Mul for Dyadic {
    type Output = Dyadic;
    fn mul(self, other: Dyadic) -> Dyadic {
        // println!("{} {}", self.numerator, other.numerator);
        Dyadic::new(
            self.numerator * other.numerator,
            self.exponent + other.exponent,
        )
        .reduce()
    }
}

impl Div for Dyadic {
    type Output = Dyadic;

    fn div(self, other: Dyadic) -> Dyadic {
        if other.numerator == 0 {
            panic!("Division with zero Dyadic!");
        }

        let result_f64 = self.to_f64() / other.to_f64();
        Dyadic::approximate(result_f64, 30).reduce() // adjust precision as needed
    }
}

impl PartialOrd for Dyadic {
    fn partial_cmp(&self, other: &Dyadic) -> Option<Ordering> {
        Some(self.to_f64().partial_cmp(&other.to_f64()).unwrap())
    }
}

impl PartialEq for Dyadic {
    fn eq(&self, other: &Self) -> bool {
        // Reduce both dyadics and compare
        let self_reduced = self.reduce();
        let other_reduced = other.reduce();

        self_reduced.numerator == other_reduced.numerator
            && self_reduced.exponent == other_reduced.exponent
    }
}

impl Eq for Dyadic {}

impl Hash for Dyadic {
    fn hash<H: std::hash::Hasher>(&self, state: &mut H) {
        self.numerator.hash(state);
        self.exponent.hash(state)
    }
}

pub fn add_vec(a: &Vec<Dyadic>, b: &Vec<Dyadic>) -> Vec<Dyadic> {
    let len = a.len().max(b.len()); // Use the length of the longer vector
    let mut result = Vec::with_capacity(len);

    for i in 0..len {
        // If index i is within bounds of self, use the value; otherwise, use Dyadic(0)
        let val_self = if i < a.len() { a[i] } else { Dyadic::zero() };
        let val_other = if i < b.len() { b[i] } else { Dyadic::zero() };
        result.push(val_self + val_other);
    }

    result
}

// Define type ComplexDyadic
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
pub struct ComplexDyadic {
    pub re: Dyadic,
    pub im: Dyadic,
}

impl ComplexDyadic {
    pub fn new(re: Dyadic, im: Dyadic) -> ComplexDyadic {
        ComplexDyadic { re, im }
    }
    pub fn abs(&self) -> f64 {
        let real = self.re.to_f64();
        let im = self.im.to_f64();
        return (real * real + im * im).sqrt();
    }

    // Add fast powering
    pub fn powi(self, power: i32) -> Self {
        if power == 0 {
            return ComplexDyadic::one();
        } else if power < 0 {
            unimplemented!("Negative powers not implemented");
        }

        fn fast_pow(base: ComplexDyadic, exp: i32) -> ComplexDyadic {
            if exp == 0 {
                ComplexDyadic::one()
            } else if exp % 2 == 0 {
                let half = fast_pow(base * base, exp / 2);
                half
            } else {
                base * fast_pow(base, exp - 1)
            }
        }

        fast_pow(self, power)
    }

    pub fn one() -> ComplexDyadic {
        ComplexDyadic::new(Dyadic::new(1, 0), Dyadic::new(0, 0))
    }

    pub fn zero() -> ComplexDyadic {
        ComplexDyadic::new(Dyadic::zero(), Dyadic::zero())
    }

    pub fn absolute_value(self) -> f64 {
        let re = self.re.powi(2);
        let im = self.im.powi(2);
        (re.to_f64() + im.to_f64()).powf(0.5)
    }

    pub fn from_i64(n: i64) -> Self {
        let d = Dyadic::from_i64(n);
        Self {
            re: d.clone(),
            im: Dyadic::zero(),
        }
    }
    /// Divide by a small integer with dyadic precision control (applies to both parts).
    pub fn div_i64(self, n: i64, bits: u32) -> Self {
        let denom = Dyadic::from_i64(n);
        Self {
            re: self.re.div_with_precision(denom.clone(), bits),
            im: self.im.div_with_precision(denom, bits),
        }
    }

    pub fn div_with_precision(&self, other: ComplexDyadic, bits: u32) -> ComplexDyadic {
        // denominator = c^2 + d^2
        let denom = other.re.clone() * other.re.clone() + other.im.clone() * other.im.clone();

        assert!(
            denom.numerator != 0,
            "division by zero complex dyadic"
        );

        // real = (ac + bd) / (c^2 + d^2)
        let num_re = self.re.clone() * other.re.clone() + self.im.clone() * other.im.clone();
        let re = num_re.div_with_precision(denom.clone(), bits);

        // imag = (bc - ad) / (c^2 + d^2)
        let num_im = self.im.clone() * other.re.clone() - self.re.clone() * other.im.clone();
        let im = num_im.div_with_precision(denom, bits);

        ComplexDyadic { re, im }
    }

    pub fn to_f64(&self) -> (f64, f64) {
        (self.re.to_f64(), self.im.to_f64())
    }
}

impl Add for ComplexDyadic {
    type Output = ComplexDyadic;
    fn add(self, other: ComplexDyadic) -> ComplexDyadic {
        println!("Adding Dyadic: {:?} + {:?}", self, other);
        ComplexDyadic::new(self.re + other.re, self.im + other.im) // Does this compile?
    }
}

impl Sub for ComplexDyadic {
    type Output = ComplexDyadic;
    fn sub(self, other: ComplexDyadic) -> ComplexDyadic {
        ComplexDyadic::new(self.re - other.re, self.im - other.im) // Is this in line with definitions of add and sub for Dyadic numbers?
    }
}

impl Mul for ComplexDyadic {
    type Output = ComplexDyadic;
    fn mul(self, other: ComplexDyadic) -> ComplexDyadic {
        return ComplexDyadic::new(
            self.re * other.re - self.im * other.im,
            self.re * other.im + self.im * other.re,
        ); // Again, we need to make sure this runs in terms of operations
    }
}

impl Div for ComplexDyadic {
    type Output = ComplexDyadic;
    fn div(self, other: ComplexDyadic) -> ComplexDyadic {
        if other != ComplexDyadic::zero() {
            ComplexDyadic::new(
                (self.re * other.re + self.im * other.im)
                    / (other.re * other.re + other.im * other.im),
                (self.im * other.re - self.re * other.im)
                    / (other.re * other.re + other.im * other.im),
            )
        } else {
            panic!("Division with ComplexDyadic zero!")
        }
    }
}

pub fn add_complex_vec(a: &Vec<ComplexDyadic>, b: &Vec<ComplexDyadic>) -> Vec<ComplexDyadic> {
    let len = a.len().max(b.len()); // Use the length of the longer vector
    let mut result = Vec::with_capacity(len);

    for i in 0..len {
        let val_self = if i < a.len() {
            a[i]
        } else {
            ComplexDyadic::zero()
        };
        let val_other = if i < b.len() {
            b[i]
        } else {
            ComplexDyadic::zero()
        };
        result.push(val_self + val_other);
    }

    result
}

pub fn sub_complex_vec(a: &Vec<ComplexDyadic>, b: &Vec<ComplexDyadic>) -> Vec<ComplexDyadic> {
    let len = a.len().max(b.len()); // Use the length of the longer vector
    let mut result = Vec::with_capacity(len);

    for i in 0..len {
        let val_self = if i < a.len() {
            a[i]
        } else {
            ComplexDyadic::zero()
        };
        let val_other = if i < b.len() {
            b[i]
        } else {
            ComplexDyadic::zero()
        };
        result.push(val_self - val_other);
    }

    result
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]

// Define interval type
pub enum Interval {
    Empty,
    Bounded(Dyadic, Dyadic, Dyadic, Dyadic),
}

// Implement constructor for type interval
impl Interval {
    pub fn new(x_lower: Dyadic, x_upper: Dyadic, y_lower: Dyadic, y_upper: Dyadic) -> Interval {
        if x_lower > x_upper {
            Interval::Empty
        } else if y_lower > y_upper {
            Interval::Empty
        } else {
            Interval::Bounded(x_lower, x_upper, y_lower, y_upper)
        }
    }
    pub fn midpoint(&self) -> Option<ComplexDyadic> {
        match self {
            Interval::Empty => None,
            Interval::Bounded(x_lower, x_upper, y_lower, y_upper) => Some(ComplexDyadic::new(
                (*x_lower + *x_upper) * Dyadic::new(1, -1),
                (*y_lower + *y_upper) * Dyadic::new(1, -1),
            )),
        }
    }
}
#[derive(Debug, Clone, Copy, PartialEq, Eq)]

// Define alphabet of words
enum Letter {
    On,
    Two,
    Three,
    Four,
}

struct Word {
    length: i32,
    // Luka : Words could be lists of vectors (Luka probably meant or)
}

// Implement display methods for defined types
impl fmt::Display for Dyadic {
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        write!(f, "{:.15}", self.to_f64())
        // Luka: This could be problematic with stricter accuracy requirements
    }
}

impl fmt::Display for Interval {
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        match self {
            Interval::Empty => write!(f, "Empty"),
            Interval::Bounded(x_min, x_max, y_min, y_max) => {
                write!(
                    f,
                    "[{}, {}] x [{}, {}]",
                    x_min, x_max, y_min, y_max
                )
            }
        }
    }
}

impl fmt::Display for ComplexDyadic {
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        if self.im.numerator >= 0 {
            write!(f, "{} + {}i", self.re, self.im)
        } else {
            write!(
                f,
                "{} - {}i",
                self.re,
                Dyadic::new(-self.im.numerator, self.im.exponent)
            )
        }
    }
}

fn split(rect: Interval, n: u8) -> Interval {
    match rect {
        Interval::Empty => Interval::Empty,
        Interval::Bounded(x_lower, x_upper, y_lower, y_upper) => {
            let x_mid = (x_lower + x_upper) * (Dyadic::new(1, -1));
            let y_mid = (y_lower + y_upper) * (Dyadic::new(1, -1));

            match n {
                1 => Interval::new(x_lower, x_mid, y_lower, y_mid),
                2 => Interval::new(x_mid, x_upper, y_lower, y_mid),
                3 => Interval::new(x_lower, x_mid, y_mid, y_upper),
                4 => Interval::new(x_mid, x_upper, y_mid, y_upper),
                _ => Interval::Empty,
            }
        }
    }
}

pub fn psi(int1: Interval, lst: &LinkedList<u8>) -> Interval {
    match lst.front() {
        None => int1,
        Some(&fst) => psi(
            split(int1, fst),
            &lst.iter().skip(1).cloned().collect(),
        ),
    }
}


#[cfg(test)]
mod tests {
    use super::*;
    use crate::holomorphic::ExpansionCoefficients;

    #[test]
    fn test_div_with_precision_real() {
        let a = Dyadic {
            numerator: 3,
            exponent: -5,
        }; // 3 * 2^-5 = 3/32
        let b = Dyadic {
            numerator: 5,
            exponent: -2,
        }; // 5 * 2^-2 = 5/4
        let q = a.div_with_precision(b, 30);

        // Expected: (3/32) / (5/4) = 3/40 = 0.075
        let approx = q.to_f64();
        assert!((approx - 0.075).abs() < 1e-9);
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
    fn derivative_basic() {
        // f(z) = z + 2z^2 + 3z^3  =>  f'(z) = 1 + 4z + 9z^2
        let f = ExpansionCoefficients {
            vector: vec![
                ComplexDyadic::zero(),      // c0
                ComplexDyadic::from_i64(1), // c1
                ComplexDyadic::from_i64(2), // c2
                ComplexDyadic::from_i64(3), // c3
            ],
        };
        let fp = f.derivative();
        assert_eq!(
            fp.vector,
            vec![
                ComplexDyadic::from_i64(1), // 1
                ComplexDyadic::from_i64(4), // 4
                ComplexDyadic::from_i64(9), // 9
            ]
        );
    }

    #[test]
    fn antiderivative_basic() {
        // f'(z) = 1 + 2z + 3z^2  =>  f(z) = z + (2/2)z^2 + (3/3)z^3 = z + z^2 + z^3
        let fprime = ExpansionCoefficients {
            vector: vec![
                ComplexDyadic::from_i64(1),
                ComplexDyadic::from_i64(2),
                ComplexDyadic::from_i64(3),
            ],
        };
        let f = fprime.antiderivative(60); // 60 dyadic bits for safety
        assert_eq!(
            f.vector,
            vec![
                ComplexDyadic::zero(),
                ComplexDyadic::from_i64(1),
                ComplexDyadic::from_i64(1),
                ComplexDyadic::from_i64(1),
            ]
        );
        // And derivative gets us back:
        assert_eq!(f.derivative().vector, fprime.vector);
    }
}
