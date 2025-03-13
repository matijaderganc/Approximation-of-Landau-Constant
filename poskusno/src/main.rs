//use num_complex::Complex;
use std::fmt::{self, UpperExp};

use std::io::Empty;
use std::ops::{Add, Mul, Sub, Div}; //division not yet implemented
use std::cmp::Ordering;


#[derive(Debug, Clone, Copy, PartialEq, Eq)]

struct Dyadic {
    numerator : i64 ,
    exponent : i32 ,
}

impl Dyadic {
    fn new(num : i64, exp : i32) -> Self {
        Dyadic {numerator : num, exponent : exp}
    }
    fn to_f64(&self) -> f64 {
        (self.numerator as f64) / (2.0f64).powi(self.exponent)
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


impl PartialOrd for Dyadic {
    fn partial_cmp(&self, other: &Dyadic) -> Option<Ordering> {
        Some(self.to_f64().partial_cmp(&other.to_f64()).unwrap())
    }
}
#[derive(Debug, Clone, Copy, PartialEq, Eq)]

enum Interval {
    Empty , 
    Bounded(Dyadic, Dyadic),
}

impl Interval {
    fn new(lower : Dyadic, upper : Dyadic) -> Interval {
        if lower > upper {
            Interval::Empty
        }
        else {
            Interval::Bounded(lower, upper)
        }

    }

    fn midpoint(&self) -> Option<Dyadic> {
        match self {
            Interval::Empty => None, 
            Interval::Bounded(a,b ) => Some(Dyadic::new(
                (a.numerator + b.numerator) / 2,
                a.exponent.min(b.exponent))), 
        }
        
    }
}
enum Letter {
    One = 1 ,
    Two = 2 ,
    Three = 3, 
    Four = 4
}
impl fmt::Display for Dyadic {
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        write!(f, "{:.6}", self.to_f64()) // Prints with 6 decimal places
    }
}
impl fmt::Display for Interval {
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        match self {
            Interval::Empty => write!(f, "Empty"),
            Interval::Bounded(x_min, x_max) => {
                write!(f, "({}, {})", x_min, x_max)
            }
        }
    }
}



fn split(first: Interval, second : Interval, n : Letter) -> (Interval, Interval) {
    match (first, second) {
        (Interval::Empty, _) => (Interval::Empty, Interval::Empty) ,
        (_, Interval::Empty) => (Interval::Empty, Interval::Empty) ,
        (Interval::Bounded(x,y ), Interval::Bounded(z,w )) => {
        let x_mid = Dyadic::new((x.numerator + y.numerator) / 2, y.exponent.min(x.exponent));
        let y_mid = Dyadic::new((z.numerator + w.numerator) / 2, w.exponent.min(z.exponent));


         match n {
             Letter::One => (Interval::new(x, x_mid), Interval::new(z, y_mid)),
             Letter::Two => (Interval::new(x_mid, y), Interval::new(z, y_mid)),
             Letter::Three => (Interval::new(x, x_mid), Interval::new(y_mid, w)),
             Letter::Four => (Interval::new(x_mid, y), Interval::new(y_mid, w)),
         }
        }
    }
} 

fn main() {
    let x = Dyadic::new(10, 5) ;
    let y = Dyadic::new(3, 5);
    let z = x + y;
    let i1 = Interval::new(y, x);
    let i2 = Interval::new(Dyadic::new(3,5), Dyadic::new(4, 5));
    let i3 = split(i1, i2, Letter::One);
    println!("{:?}", i3);
}

