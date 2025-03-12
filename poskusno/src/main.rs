//use num_complex::Complex;

use std::ops::Add;

#[derive(Debug, Clone, Copy, PartialEq)]

struct Dyadic {
    numerator : i64 ,
    exponent : i32 ,
}

impl Dyadic {
    fn new(num : i64, exp : i32) -> Self {
        Dyadic {numerator : num, exponent : exp}
    }
    fn to_float(&self) -> f64 {
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

fn main() {
    let mut x = Dyadic::new(10, 5) ;
    let y = Dyadic::new(3, 5);
    let z = x + y;
    println!("{}", x.to_float());
    println!("{}", y.to_float());
    println!("{}", z.to_float());


}

