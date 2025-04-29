use crate::dyadic::{Dyadic, ComplexDyadic, Interval, psi} ;

pub struct Holomorphic {
    pub coef : Vec<ComplexDyadic> 
}
impl Holomorphic {
    pub fn new(v : Vec<ComplexDyadic>) -> Self {
        return Holomorphic{coef : v}
    }

    pub fn eval(&self, z : ComplexDyadic) -> ComplexDyadic {
        let mut result = ComplexDyadic::zero();

        for (k, coeff) in self.coef.iter().enumerate() {
            let zk = z.powi(k as i32);
            result = result + (*coeff * zk);
        }
        result
    }
}
