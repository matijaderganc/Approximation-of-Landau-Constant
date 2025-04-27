// TODO: 
// implement holomorphic funtions as vectors, have addition, subtraction, mult (not sure if needed)
// not yet sure if we should implement them just for dyadic complex numbers
// implement derivate, integral (choose +C so that it maps 0 to 0)
// implement methods eval, absolute value

// Define the bounding sequences, used to limit the set of all sequences. By definition, the series (m_i)p^i converges for all p from [0, 1).
// It may be necessary to expand this implementation, depending on needs. 

// Call the Dyadic type from dyadic.rs: 
use landau::dyadic::{ComplexDyadic, Dyadic} ;

struct BoundingSequence {
    n_th: fn(u32) -> Dyadic,
}

impl BoundingSequence {
    fn new(function: fn(u32) -> Dyadic) -> Self {
        BoundingSequence {
            n_th: function,
        }
    }
}

// Implement the sequences, used to identify complex holomorphic functions. 
// By definition of any given sequence, the real an imaginary part of a_i are both bounded by m_i for a sequence m_n of the type BoundingSequence. 
struct ExpansionCoefficients {
    n_th : fn(u32) -> ComplexDyadic
}

impl ExpansionCoefficients {
    fn new(function : fn(u32) -> ComplexDyadic) -> Self {
        ExpansionCoefficients {
            n_th : function, 
        }
    }
}

// A given holomorphic function consists of a bounding sequence (type BoundingSequence) and a bounded sequence, which represents the coefficients in its expansion (type ExpansionCoefficients)
// Functions are (for now) defined only for ComplexDyadic numbers, not arbitrary complex numbers. This could need to be changed later on. 
struct ComplexFunction {
    bounding_sequence: BoundingSequence,
    expansion_coefficients: ExpansionCoefficients,
    function: Box<dyn Fn(ComplexDyadic) -> ComplexDyadic + Send + Sync>,
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

    /// caller-supplied term limit
    pub fn new_with_upper_limit(
        bounding_sequence: BoundingSequence,
        expansion_coefficients: ExpansionCoefficients,
        upper_limit_of_summation: u32,
    ) -> Self {
        // grab only the raw function pointer (Copy, no Clone needed)
        let nth_fn = expansion_coefficients.n_th;

        let func = move |z: ComplexDyadic| {
            let mut sum = ComplexDyadic::zero();
            for i in 1..=upper_limit_of_summation {
                let coeff = (nth_fn)(i);
                sum = sum + coeff * z.powi(i as i32);
            }
            ComplexDyadic::one() + sum
        };

        Self {
            bounding_sequence,
            expansion_coefficients, // stored intact
            function: Box::new(func),
            upper_limit_of_summation,
        }
    }
}
