// use std::intrinsics::sqrtf64;
use std::cmp::Ordering;
use std::collections::LinkedList;
use std::fmt::{self};
use std::hash::Hash;
use std::ops::{Add, Div, Mul, Sub};

#[derive(Debug, Clone, Copy)]

// Define dyadic type as numerator * 2 ^ exponent
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

    // used to get a float from a dyadic number
    pub fn to_f64(&self) -> f64 {
        (self.numerator as f64) * (2.0f64).powi(self.exponent)
    }

    // This is used to determine if a number is zero, which need to be checked before division. Also helpful when testing
    pub fn zero() -> Dyadic {
        Dyadic {
            numerator: 0,
            exponent: 0,
        }
    }

    // Create a Dyadic with integer value
    pub fn from_i64(n: i64) -> Dyadic {
        Dyadic {
            numerator: n as i128,
            exponent: 0,
        }
    }

    // Reduce a dyadic number by dividing numerator and increasing exponent. Two dyadics which are equal when reduced, are viewed as the same dyadics, forming an equivalence class. This should also be called after every operation to ensure we keep the numerator as small as possible.
    pub fn reduce(self) -> Dyadic {
        if self.numerator == 0 {
            return Dyadic::zero();
        }
        let mut n = self.numerator;
        let mut e = self.exponent;
        while n % 2 == 0 {
            n /= 2;
            e += 1;
        }
        Dyadic::new(n, e)
    }

    // Fast power implementation to reduce computational time
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

    // Returns a dyadic number which approximates x while setting the maximum exponent
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

    // Division where we consider our precision
    pub fn div_with_precision(&self, denom: Dyadic, bits: u32) -> Dyadic {
        assert!(denom.numerator != 0, "division by zero dyadic"); // Error

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
        Dyadic::approximate(result_f64, 30).reduce() // adjust precision as needed. This may also be a bad way of doing this to be honest
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

    // Absolute value of complex dyadic number can be calculated by converting to f64
    pub fn abs(&self) -> f64 {
        let real = self.re.to_f64();
        let im = self.im.to_f64();
        (real * real + im * im).sqrt()
    }

    // Add fast powering for better performance
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
                fast_pow(base * base, exp / 2) //return half
                
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

    // Create new ComplexDyadic from integer, very useful when operating with functions
    pub fn from_i64(n: i64) -> Self {
        let d = Dyadic::from_i64(n);
        Self {
            re: d,
            im: Dyadic::zero(),
        }
    }
    // Divide by a small integer with dyadic precision control
    pub fn div_i64(self, n: i64, bits: u32) -> Self {
        let denom = Dyadic::from_i64(n);
        Self {
            re: self.re.div_with_precision(denom, bits),
            im: self.im.div_with_precision(denom, bits),
        }
    }

    pub fn div_with_precision(&self, other: ComplexDyadic, bits: u32) -> ComplexDyadic {
        // denominator = c^2 + d^2
        let denom = other.re * other.re + other.im * other.im;

        assert!(
            denom.numerator != 0,
            "division by zero complex dyadic"
        );

        // real = (ac + bd) / (c^2 + d^2)
        let num_re = self.re * other.re + self.im * other.im;
        let re = num_re.div_with_precision(denom, bits);

        // im = (bc - ad) / (c^2 + d^2)
        let num_im = self.im * other.re - self.re * other.im;
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
        ComplexDyadic::new(
            self.re * other.re - self.im * other.im,
            self.re * other.im + self.im * other.re,
        ) // Again, we need to make sure this runs in terms of operations
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

// Adding to vectors of dyadic numbers which is then used when operating with BoundingSequences and ExpansionCoefficients
pub fn add_vec(a: &[Dyadic], b: &[Dyadic]) -> Vec<Dyadic> {
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

// Adding and subtracting complex vectors
pub fn add_complex_vec(a: &[ComplexDyadic], b: &[ComplexDyadic]) -> Vec<ComplexDyadic> {
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

pub fn sub_complex_vec(a: &[ComplexDyadic], b: &[ComplexDyadic]) -> Vec<ComplexDyadic> {
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
        if x_lower > x_upper || y_lower > y_upper {
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
// Implement display methods for defined types
impl fmt::Display for Dyadic {
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        write!(f, "{:.15}", self.to_f64())
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





