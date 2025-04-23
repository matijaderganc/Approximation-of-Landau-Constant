use std::collections::LinkedList;
//use num_complex::Complex;
use std::fmt::{self, UpperExp};

use std::io::Empty;
use std::ops::{Add, Mul, Sub, Div}; //division not yet implemented
use std::cmp::Ordering;


#[derive(Debug, Clone, Copy, PartialEq, Eq)]

// Define dyadic type
struct Dyadic {
    numerator : i128 ,
    exponent : i32 ,
}

impl Dyadic {
    fn new(num : i128, exp : i32) -> Self {
        Dyadic {numerator : num, exponent : exp}
    }
    fn to_f64(&self) -> f64 {
        (self.numerator as f64) / (2.0f64).powi(self.exponent)
    }
    // This is used to determine if a number is zero, which need to be checked before division 
    fn zero() -> Dyadic {
        Dyadic {numerator : 0, exponent : 0} 
    }
}

impl Add for Dyadic {  
    type Output = Dyadic; 
    fn add(self, other : Dyadic) -> Dyadic {
        let exp_diff = self.exponent - other.exponent;
        if exp_diff > 0 {
            Dyadic::new(self.numerator + (other.numerator << exp_diff), self.exponent)
        }
        else {
            Dyadic::new((self.numerator << -exp_diff) + other.numerator, other.exponent)
        }
    }
}
impl Sub for Dyadic {
    type Output = Dyadic;

    fn sub(self, other: Dyadic) -> Dyadic {
        let exp_diff = self.exponent - other.exponent;

        if exp_diff > 0 {
            Dyadic::new(self.numerator - (other.numerator << exp_diff), self.exponent)
        } else {
            Dyadic::new((self.numerator << -exp_diff) - other.numerator, other.exponent)
        }
    }
}

impl Mul for Dyadic {
    type Output = Dyadic;
    fn mul(self, other: Dyadic) -> Dyadic {
        Dyadic::new(self.numerator * other.numerator, self.exponent + other.exponent)
    }
}

impl Div for Dyadic {
    type Output = Dyadic
    fn div(self, other : Dyadic) -> Dyadic {
        if other != Dyadic::zero() {

        }
    }
}

impl PartialOrd for Dyadic {
    fn partial_cmp(&self, other: &Dyadic) -> Option<Ordering> {
        Some(self.to_f64().partial_cmp(&other.to_f64()).unwrap())
    }
}

// Define type ComplexDyadic
struct ComplexDyadic {
    re : Dyadic 
    im : Dyadic 
}

impl ComplexDyadic {
    fn new(re : Dyadic, im : Dyadic) -> ComplexDyadic {
        ComplexDyadic{re : re, im : im}
    }
}

impl Add for ComplexDyadic {
    type Output = ComplexDyadic
    fn add(self, other : ComplexDyadic) -> ComplexDyadic {
        return ComplexDyadic::new(self.re + other.re, self.im + other.im) // Does this compile?
    }
}

impl Sub for ComplexDyadic {
    type Output = ComplexDyadic
    fn sub(self, other : ComplexDyadic) -> ComplexDyadic {
        return ComplexDyadic::new(self.re - other.re, self.im - other.im) // Is this in line with definitions of add and sub for Dyadic numbers? 
    }
}

impl Mul for ComplexDyadic {
    type Output = ComplexDyadic
    fn mul(self, other : ComplexDyadic) -> ComplexDyadic {
        return ComplexDyadic::new(self.re * other.re - self.im * other.im, self.re * other.im + self.im * other.re) // Again, we need to make sure this runs in terms of operations
    }
}

impl Div for ComplexDyadic {
    type Output = ComplexDyadic
    fn div(self, other : ComplexDyadic) -> ComplexDyadic {
        if other.re != Dyadic::zero() && other.im != Dyadic::zero() {
            ComplexDyadic::new((self.re * other.re + self.im * other.im) / (other.re * other.re + other.im * other.im), (self.im * other.re - self.re * other.im) / (other.re * other.re + other.im * other.im)) // Check if operations work as needed 
        }
    }
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]

// Define interval type 
enum Interval {
    Empty , 
    Bounded(Dyadic, Dyadic, Dyadic, Dyadic),
}

// Implement constructor for type interval 
impl Interval {
    fn new(x_lower : Dyadic, x_upper : Dyadic, y_lower : Dyadic, y_upper : Dyadic) -> Interval {
        if x_lower > x_upper {
            Interval::Empty
        }
        else if y_lower > y_upper  {
            Interval::Empty
        } 
        else {
            Interval::Bounded(x_lower, x_upper, y_lower, y_upper)
        }

    }

  
}
#[derive(Debug, Clone, Copy, PartialEq, Eq)]

// Define alphabet of words 
enum Letter {
    On,
    Two,
    Three, 
    Four
}

struct Word {
    length : i32, 
    // Luka : Words could be lists of vectors
}

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

fn split(rect : Interval, n : u8) -> Interval {
    match rect {
        Interval::Empty => Interval::Empty,
        Interval::Bounded(x_lower,x_upper, y_lower, y_upper ) => {
            let x_mid = (x_lower + x_upper) * (Dyadic::new(1, -1));
            let y_mid = (y_lower + y_upper) * (Dyadic::new(1, -1));

            match n {
                1 => Interval::new(x_lower, x_mid, y_lower, y_mid),
                2 => Interval::new(x_mid, x_upper, y_lower, y_mid),
                3 => Interval::new(x_lower, x_mid, y_mid, y_upper),
                4 => Interval::new(x_mid, x_upper, y_mid, y_upper),
                _ => Interval::Empty
            }
        }
    }
} 

fn psi(int1 : Interval, lst : &LinkedList<u8>) -> Interval {
    // Kaj je tukaj misljeno z lst.front()?
    match lst.front() {
        None => int1 ,
        Some(&fst) => psi(split(int1, fst), &lst.iter().skip(1).cloned().collect())
    }
}

fn main() {
    let y = Dyadic::new(10, 5) ;
    let x = Dyadic::new(3, 5);
    let z = Dyadic::new(4, 6);
    let w = Dyadic::new(7, 6);
    println!("{}", x.to_f64()) ;
    let i1 = Interval::new(x, y, z, w);
    let mut word = LinkedList::new();
    for _ in 0..5 {
        word.push_back(1);
    }
    let i3 = psi(i1, &word) ;
    // print!("{}", i3)
   
}

