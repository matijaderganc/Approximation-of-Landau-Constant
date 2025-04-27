// TODO: 
// implement holomorphic funtions as vectors, have addition, subtraction, mult (not sure if needed)
// not yet sure if we should implement them just for dyadic complex numbers
// implement derivate, integral (choose +C so that it maps 0 to 0)
// implement methods eval, absolute value

// Call the Dyadic type from dyadic.rs: 
use landau::dyadic::{ComplexDyadic, Dyadic} ;
// Call operations 
use std::ops::{Add, Mul, Sub, Div}; //division not yet implemented
use std::sync::Arc;

// Define the bounding sequences, used to limit the set of all sequences. By definition, the series (m_i)p^i converges for all p from [0, 1).
// It may be necessary to expand this implementation, depending on needs. 
#[derive(Clone)]
struct BoundingSequence {
    n_th: Arc<dyn Fn(u32) -> Dyadic>,
}

impl Add for BoundingSequence {
    type Output = Self;

    fn add(self, other: Self) -> Self {
        let left  = self.n_th;
        let right = other.n_th;

        Self {
            n_th: Arc::new(move |n| {
                let a = left(n);
                let b = right(n);
                if a > b { a } else { b }
            }),
        }
    }
}


impl BoundingSequence {
    pub fn new<F>(f: F) -> Self
    where
        F: Fn(u32) -> Dyadic + Send + Sync + 'static,
    {
        Self { n_th: Arc::new(f) }
    }
}

// Implement the sequences, used to identify complex holomorphic functions. 
// By definition of any given sequence, the real an imaginary part of a_i are both bounded by m_i for a sequence m_n of the type BoundingSequence. 
#[derive(Clone)]
struct ExpansionCoefficients {
    n_th : Arc<dyn Fn(u32) -> ComplexDyadic>,
}

impl Add for ExpansionCoefficients {
    type Output = Self;

    fn add(self, other: Self) -> Self {
        let left  = self.n_th;
        let right = other.n_th;

        Self {
            n_th: Arc::new(move |n| left(n) + right(n)),
        }
    }
}

impl ExpansionCoefficients {
    pub fn new<F>(f: F) -> Self
    where
        F: Fn(u32) -> ComplexDyadic + Send + Sync + 'static,
    {
        Self { n_th: Arc::new(f) }
    }
}

// A given holomorphic function consists of a bounding sequence (type BoundingSequence) and a bounded sequence, which represents the coefficients in its expansion (type ExpansionCoefficients)
// Functions are (for now) defined only for ComplexDyadic numbers, not arbitrary complex numbers. This could need to be changed later on. 
struct ComplexFunction {
    bounding_sequence: BoundingSequence,
    expansion_coefficients: ExpansionCoefficients,
    function: Box<dyn Fn(ComplexDyadic) -> ComplexDyadic>,
    // Upper limit defaults to 100, but a separate constructor if defined 
    upper_limit_of_summation: u32,
}

impl ComplexFunction {
    /// default = 100 terms
    pub fn new(
        bounding_sequence: BoundingSequence,
        expansion_coefficients: ExpansionCoefficients,
    ) -> Self {
        Self::new_with_upper_limit(bounding_sequence, expansion_coefficients, 100)
    }

    pub fn new_with_upper_limit(
        bounding_sequence: BoundingSequence,
        expansion_coefficients: ExpansionCoefficients,
        upper_limit_of_summation: u32,
    ) -> Self {
        let nth_fn = expansion_coefficients.n_th.clone();

        let func = move |z: ComplexDyadic| {
            let mut sum = ComplexDyadic::zero();
            for i in 1..=upper_limit_of_summation {
                sum = sum + nth_fn(i) * z.powi(i as i32);
            }
            ComplexDyadic::one() + sum
        };

        Self {
            bounding_sequence,
            expansion_coefficients,          // still intact – not moved
            function: Box::new(func),
            upper_limit_of_summation,
        }
    }
}

impl Add for ComplexFunction { 
    type Output = ComplexFunction;
    fn add(self, other : ComplexFunction) -> ComplexFunction {
        ComplexFunction::new_with_upper_limit(
            self.bounding_sequence + other.bounding_sequence,
            self.expansion_coefficients + other.expansion_coefficients,
            std::cmp::max(self.upper_limit_of_summation, other.upper_limit_of_summation)
            )    
    }
}
