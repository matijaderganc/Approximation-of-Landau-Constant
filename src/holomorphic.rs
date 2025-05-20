// TODO:
// - implement holomorphic mult (not sure if needed)
// - not yet sure if we should implement them just for dyadic complex numbers. This may need to be revised later on. Functions now work only on ComplexDyadic
// - implement derivate, integral (choose +C so that it maps 0 to 0) -> This is the most important aspect of this struct, as it will be used later on
// - do some basic testing to see function return expected values

// Call the Dyadic type from dyadic.rs:
use crate::dyadic::{add_complex_vec, add_vec, sub_complex_vec, ComplexDyadic, Dyadic};
use std::ops::{Add, Bound, Div, Mul, Sub}; //division not yet implemented
use std::sync::Arc;
use std::vec;

// Function to map vector to function which returns n-th element.
pub fn vec_to_sequence(vec: &Vec<Dyadic>) -> impl Fn(u32) -> Dyadic + '_ {
    move |n: u32| {
        vec.get(n as usize) // safe, no panic
            .cloned() // use `copied()` if Dyadic implements Copy
            .unwrap_or_else(Dyadic::zero)
    }
}

pub fn comp_vec_to_sequence(vec: &Vec<ComplexDyadic>) -> impl Fn(u32) -> ComplexDyadic + '_ {
    move |n: u32| {
        vec.get(n as usize)
            .cloned()
            .unwrap_or_else(ComplexDyadic::zero)
    }
}
// Define the bounding sequences, used to limit the set of all sequences. By definition, the series (m_i)p^i converges for all p from [0, 1).
// It may be necessary to expand this implementation, depending on needs.
#[derive(Clone)]
pub struct BoundingSequence {
    vector: Vec<Dyadic>,
}

impl Add for BoundingSequence {
    type Output = Self;

    fn add(self, other: Self) -> Self {
        Self {
            vector: add_vec(&self.vector, &other.vector),
        }
    }
}

impl BoundingSequence {
    pub fn new(vector: Vec<Dyadic>) -> BoundingSequence {
        BoundingSequence { vector: vector }
    }
}
// Implement the sequences, used to identify complex holomorphic functions.
// By definition of any given sequence, the real an imaginary part of a_i are both bounded by m_i for a sequence m_n of the type BoundingSequence.
#[derive(Clone, Debug)]
// a_1, the first element should always be zero here.
pub struct ExpansionCoefficients {
    vector: Vec<ComplexDyadic>,
}

impl Add for ExpansionCoefficients {
    type Output = Self;

    fn add(self, other: Self) -> Self {
        ExpansionCoefficients {
            vector: add_complex_vec(&self.vector, &other.vector),
        }
    }
}

impl Sub for ExpansionCoefficients {
    type Output = Self;

    fn sub(self, other: Self) -> ExpansionCoefficients {
        ExpansionCoefficients {
            vector: sub_complex_vec(&self.vector, &other.vector),
        }
    }
}

impl ExpansionCoefficients {
    pub fn new(vector: Vec<ComplexDyadic>) -> ExpansionCoefficients {
        ExpansionCoefficients { vector: vector }
    }
}

// A given holomorphic function consists of a bounding sequence (type BoundingSequence) and a bounded sequence, which represents the coefficients in its expansion (type ExpansionCoefficients)
// Functions are (for now) defined only for ComplexDyadic numbers, not arbitrary complex numbers. This could need to be changed later on.
pub struct ComplexFunction {
    bounding_sequence: BoundingSequence,
    pub expansion_coefficients: ExpansionCoefficients,
    // Upper limit defaults to 100, but a separate constructor if defined
    upper_limit_of_summation: u32,
}

impl ComplexFunction {
    /// default = 100 terms
    pub fn new(
        bounding_sequence: BoundingSequence,
        expansion_coefficients: ExpansionCoefficients,
    ) -> ComplexFunction {
        let len = expansion_coefficients.vector.len();
        let len = len as u32;

        ComplexFunction {
            bounding_sequence: bounding_sequence,
            expansion_coefficients: expansion_coefficients,
            upper_limit_of_summation: len,
        }
    }

    pub fn eval(&self, z: ComplexDyadic) -> ComplexDyadic {
        let mut sum = ComplexDyadic::zero();
        for i in 0..=self.expansion_coefficients.vector.len() - 1 {
            println!("{}, {}", sum, i);
            sum = sum + self.expansion_coefficients.vector[i as usize] * (z.powi((i + 1) as i32));
            println!("{}", sum)
        }
        return sum + ComplexDyadic::one();
    }

    // This could need to be checked if it works properly.
    pub fn derivative(&self) -> ComplexFunction {
        let coefficients = &self.expansion_coefficients.vector;
        let length = coefficients.len();
        let mut derivative_coefficients: Vec<ComplexDyadic> = Vec::with_capacity(length - 1);

        // Calculate the derivative coefficients

        for i in 1..length {
            derivative_coefficients.push(
                coefficients[i] * ComplexDyadic::new(Dyadic::new(i as i128, 0), Dyadic::zero()),
            );
        }
        ComplexFunction::new(
            self.bounding_sequence.clone(),
            ExpansionCoefficients::new(derivative_coefficients),
        )
    }

    // This still needs to be fixed so it maps zero correctly
    pub fn antiderivative(&self) -> ComplexFunction {
        let sequence = &self.expansion_coefficients.vector;
        let mut antiderivative_sequence = Vec::with_capacity(sequence.len() + 1);
        antiderivative_sequence.push(ComplexDyadic::one());
        for i in 0..sequence.len() {
            antiderivative_sequence.push(
                sequence[i] / ComplexDyadic::new(Dyadic::new(i as i128 + 2, 0), Dyadic::zero()),
            )
        }
        ComplexFunction {
            bounding_sequence: self.bounding_sequence.clone(),
            expansion_coefficients: ExpansionCoefficients::new(antiderivative_sequence),
            upper_limit_of_summation: self.upper_limit_of_summation + 1,
        }
    }
}

impl Add for ComplexFunction {
    type Output = ComplexFunction;
    fn add(self, other: ComplexFunction) -> ComplexFunction {
        ComplexFunction::new(
            self.bounding_sequence + other.bounding_sequence,
            self.expansion_coefficients + other.expansion_coefficients,
        )
    }
}

impl Sub for ComplexFunction {
    type Output = Self;

    fn sub(self, other: Self) -> Self {
        Self::new(
            self.bounding_sequence + other.bounding_sequence,
            self.expansion_coefficients - other.expansion_coefficients,
        )
    }
}
