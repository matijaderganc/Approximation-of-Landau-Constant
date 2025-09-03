use std::ops::{Add, Sub};

use crate::dyadic::{ComplexDyadic, Dyadic, add_complex_vec, add_vec, sub_complex_vec};

// Function to map vector to function which returns n-th element.
pub fn vec_to_sequence(vec: &Vec<Dyadic>) -> impl Fn(u32) -> Dyadic + '_ {
    move |n: u32| {
        vec.get(n as usize) // safe, no panic
            .cloned() // use `copied()` if Dyadic implements Copy
            .unwrap_or_else(Dyadic::zero)
    }
}

pub fn comp_vec_to_sequence(vec: Vec<ComplexDyadic>) -> impl Fn(u32) -> ComplexDyadic {
    move |n: u32| {
        vec.get(n as usize) // safe, no panic
            .cloned() // use `copied()` if Dyadic: Copy
            .unwrap_or_else(ComplexDyadic::zero)
    }
}
// Define the bounding sequences, used to limit the set of all sequences. By definition, the series (m_i)p^i converges for all p from [0, 1).
// It may be necessary to expand this implementation, depending on needs.
#[derive(Clone)]
pub struct BoundingSequence {
    vector: Vec<Dyadic>,
}

// Used to add two BoundingSequences together
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
        BoundingSequence { vector }
    }
}

// Implement the sequences, used to identify complex holomorphic functions, as are listed in the paper
// By definition for any given sequence, the real an imaginary part of a_i are both bounded by m_i for a sequence m_n of the type BoundingSequence.
// Different names are kept to keep things clear, despite them being the same struct under the hood.
#[derive(Clone)]
pub struct ExpansionCoefficients {
    pub vector: Vec<ComplexDyadic>,
}

// Similar to the type BoundingSequence, we define the addition of them
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
        ExpansionCoefficients { vector }
    }

    // f'(z) where f(z) = sum c_k z^k
    pub fn derivative(&self) -> ExpansionCoefficients {
        let n = self.vector.len();
        if n <= 1 {
            return ExpansionCoefficients {
                vector: vec![ComplexDyadic::zero()],
            };
        }
        let mut out = Vec::with_capacity(n - 1);
        for k in 1..n {
            let factor = ComplexDyadic::from_i64(k as i64);
            out.push(self.vector[k].clone() * factor);
        }

        ExpansionCoefficients { vector: out }
    }

    //  integrate f(z)dz with constant of integration = 0
    pub fn antiderivative(&self, bits: u32) -> ExpansionCoefficients {
        let n = self.vector.len();
        let mut out = Vec::with_capacity(n + 1);

        // constant term = 0
        out.push(ComplexDyadic::zero());

        for k in 0..n {
            out.push(self.vector[k].clone().div_i64((k as i64) + 1, bits));
        }

        ExpansionCoefficients { vector: out }
    }

    // Evaluate the function, represented with its Taylor expansion in a given point using Horners algorithm
    pub fn eval(&self, z: &ComplexDyadic) -> ComplexDyadic {
        let mut acc = ComplexDyadic::zero();
        for coeff in self.vector.iter().rev() {
            acc = acc * z.clone() + coeff.clone();
        }
        acc
    }

    // Define mu_prime as in source paper
    pub fn mu_prime(&self, r: f64) -> f64 {
        self.vector.iter().enumerate().skip(1) // skip constant term
            .map(|(n, c)| (n as f64) * c.abs() * r.powi((n-1) as i32))
            .sum()
    }

    // Define mu_double_prime as in the source paper
    pub fn mu_double_prime(&self, r: f64) -> f64 {
        self.vector
            .iter()
            .enumerate()
            .skip(2)
            .map(|(n, c)| (n as f64) * (n - 1) as f64 * c.abs() * r.powi((n - 2) as i32))
            .sum()
    }
}

// A given holomorphic function consists of a bounding sequence (type BoundingSequence) and a bounded sequence, which represents the coefficients in its expansion (type ExpansionCoefficients)
// Functions are (for now) defined only for ComplexDyadic numbers, not arbitrary complex numbers. It may be useful to change the definition of a ComplexFunction to take only its expansion coefficients
// Upper limit of summation is also not useful probably and was initially used to define the number of terms to consider when evaluating
pub struct ComplexFunction {
    bounding_sequence: BoundingSequence,
    pub expansion_coefficients: ExpansionCoefficients,
    // Upper limit defaults to 1000, but a separate constructor if defined
    upper_limit_of_summation: u32,
}

impl ComplexFunction {
    /// default = 1000 terms
    pub fn new(
        bounding_sequence: BoundingSequence,
        expansion_coefficients: ExpansionCoefficients,
    ) -> ComplexFunction {
        let len = expansion_coefficients.vector.len();
        let len = len as u32;

        ComplexFunction {
            bounding_sequence,
            expansion_coefficients,
            upper_limit_of_summation: len,
        }
    }

    pub fn eval(&self, z: &ComplexDyadic) -> ComplexDyadic {
        let mut acc = ComplexDyadic::zero();
        for coeff in self.expansion_coefficients.vector.iter().rev() {
            acc = acc * z.clone() + coeff.clone();
        }
        acc
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
        antiderivative_sequence.push(ComplexDyadic::zero());
        for i in 0..sequence.len() {
            antiderivative_sequence.push(
                sequence[i] / ComplexDyadic::new(Dyadic::new(i as i128 + 1, 0), Dyadic::zero()),
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
