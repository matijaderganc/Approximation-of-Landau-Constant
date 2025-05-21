use std::collections::LinkedList;
//use num_complex::Complex;
use std::fmt::{self, write, UpperExp};

// use std::intrinsics::sqrtf64;
use std::cmp::Ordering;
use std::io::Empty;
use std::ops::{Add, Div, Mul, Sub};

use plotters::element::ComposedElement; //division not yet implemented

#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]

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

    pub fn reduce(self) -> Dyadic {
        if self.numerator % 2 == 0 {
            return Dyadic::new(self.numerator / 2, self.exponent + 1).reduce();
        }
        self
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
}

impl Add for Dyadic {
    type Output = Dyadic;
    fn add(self, other: Dyadic) -> Dyadic {
        let exp_diff = self.exponent - other.exponent;
        if exp_diff > 0 {
            Dyadic::new(
                other.numerator + (self.numerator << exp_diff),
                other.exponent,
            )
        } else {
            Dyadic::new(
                (other.numerator << -exp_diff) + self.numerator,
                self.exponent,
            )
        }
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
        } else {
            Dyadic::new(
                -(other.numerator << -exp_diff) + self.numerator,
                self.exponent,
            )
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
    }
}

impl Div for Dyadic {
    type Output = Dyadic;
    fn div(self, other: Dyadic) -> Dyadic {
        if other.numerator != 0 {
            Dyadic::new(
                self.numerator / other.numerator,
                self.exponent - other.exponent,
            )
        } else {
            panic!("Division with zero Dyadic!")
        }
    }
}

impl PartialOrd for Dyadic {
    fn partial_cmp(&self, other: &Dyadic) -> Option<Ordering> {
        Some(self.to_f64().partial_cmp(&other.to_f64()).unwrap())
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
        ComplexDyadic { re: re, im: im }
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
        }

        let mut result = ComplexDyadic::one();
        for _ in 0..power {
            // println!("{result}");
            result = result * self;
        }
        result
        } else if power < 0 {
            // If you want to support negative powers, handle here.
            // For example, use self.inverse().powi(-power)
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
        return (re.to_f64() + im.to_f64()).powf(0.5);
    }
}

impl Add for ComplexDyadic {
    type Output = ComplexDyadic;
    fn add(self, other: ComplexDyadic) -> ComplexDyadic {
        return ComplexDyadic::new(self.re + other.re, self.im + other.im); // Does this compile?
    }
}

impl Sub for ComplexDyadic {
    type Output = ComplexDyadic;
    fn sub(self, other: ComplexDyadic) -> ComplexDyadic {
        return ComplexDyadic::new(self.re - other.re, self.im - other.im); // Is this in line with definitions of add and sub for Dyadic numbers?
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
                write!(f, "[{}, {}] x [{}, {}]", x_min, x_max, y_min, y_max)
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
        Some(&fst) => psi(split(int1, fst), &lst.iter().skip(1).cloned().collect()),
    }
}
