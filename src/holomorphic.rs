// TODO: 
// implement holomorphic funtions as vectors, have addition, subtraction, mult (not sure if needed)
// not yet sure if we should implement them just for dyadic complex numbers
// implement derivate, integral (choose +C so that it maps 0 to 0)
// implement methods eval, absolute value

// Define the bounding sequences, used to limit the set of all sequences. By definition, the series (m_i)p^i converges for all p from [0, 1).
// It may be necessary to expand this implementation, depending on needs. 

// Call the Dyadic type from dyadic.rs: 
use landau::dyadic::{Dyadic} ;

struct BoundingSequences {
    n_th: fn(u32) -> Dyadic,
}

impl BoundingSequences {
    fn new(function: fn(u32) -> Dyadic) -> Self {
        BoundingSequences {
            n_th: function,
        }
    }
}

// Implement the sequences, used to identify complex holomorphic functions. 
// By definition of any given sequence, the real an imaginary part of a_i are both bounded by m_i for from sequence m_n of the type BoundingSequence. 
struct HolomorphicSequence {
    n_th : fn(u32) -> Dyadic
}

struct ComplexFunction {

}