use landau::covering_grids::{
    self, create_covering_grid, extreme_points, grid_approx, grid_complement, unit_disk_n,
};
use landau::dyadic::{psi, ComplexDyadic, Dyadic, Interval};
use landau::plot::{self, plot_covering_grid, plot_set};

use landau::holomorphic::{
    comp_vec_to_sequence, vec_to_sequence, BoundingSequence, ComplexFunction, ExpansionCoefficients,
};
use landau::psi::{generate_all_words, generate_word, mu_first, mu_second, psi_infinity};
use plotters::prelude::*;

use std::collections::LinkedList;

fn main() -> Result<(), Box<dyn std::error::Error>> {
    let x = Dyadic::new(12, -2);

    println!("{}", Dyadic::approximate(0.3, 20));
    println!("{}", x.reduce());
    let m_seq = vec![x.clone(), x.clone(), x.clone(), x.clone()];
    let t_seq = vec![0, 1, 0, 1, 2, 0, 1, 2, 2, 0, 1, 2, 3, 0, 1, 2, 3, 0];
    let word = vec![1, 1, 1, 3, 2, 1, 3, 1, 2, 3, 1, 2, 3, 1, 3, 2, 1, 2, 3, 2];
    let holo = psi_infinity(&m_seq, &t_seq, &word);
    for t in &holo {
        println!("{}", t)
    }
    println!("done");
    let points1 = unit_disk_n(-6);

    let f = ComplexFunction::new(
        BoundingSequence::new(m_seq),
        ExpansionCoefficients::new(holo),
    );
    println!("{}", f.eval(ComplexDyadic::new(x.clone(), x.clone())));
    let mut points2 = Vec::new();
    for x in &points1 {
        points2.push(f.eval(*x))
    }
    let grid = create_covering_grid(&points2, Dyadic::new(1, -2), Dyadic::new(1, -4));
    plot_set(&points2, "image2.png");
    plot_covering_grid(&grid, "covering_grid.png")?;
    plot_covering_grid(
        &grid_complement(&grid, Dyadic::new(1, -4)),
        "complement.png",
    )?;
    println!("{:?}", generate_all_words(10).len());
    println!(
        "approximation is {}",
        grid_approx(&grid, Dyadic::new(1, -4)).unwrap()
    );
    Ok(())
}
// add these to tests!!!

// let y = Dyadic::new(5, -2);
// let z = Dyadic::new(4, -1);
// let w = Dyadic::new(7, -1);
// let i1 = Interval::new(x, y, z, w);

// let mut word = LinkedList::new();
// for _ in 0..5 {
//     word.push_back(1);
// }
// let i3 = psi(i1, &word) ;
// print!("{}", i3) ;

// let alpha = ComplexDyadic::new(x, y);
// let beta = ComplexDyadic::new(z, w) ;
// let gamma = i1.midpoint().unwrap() ;
// println!("{}", alpha * beta) ;
// println!("{}, {}", i1, gamma)

// for a in holo{
//     println!("{}, {}", a, a.abs())
// }
// println!("{}", mu_first(&2.0, &x)) ;
// println!("{}", mu_second(&2.0, &x))

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_addition_dyadic() {
        let result = Dyadic::new(1, 1) + Dyadic::new(1, 1);
        assert_eq!(result, Dyadic::new(2, 1));
    }

    #[test]
    fn test_multiplication_dyadic() {
        let result = Dyadic::new(5, 2) * Dyadic::new(3, 1);
        assert_eq!(result, Dyadic::new(15, 3));
    }

    #[test]
    fn test_addition_complex_dyadic() {
        let result = ComplexDyadic::new(Dyadic::new(1, 0), Dyadic::new(2, 0))
            + ComplexDyadic::new(Dyadic::new(2, 0), Dyadic::new(3, 0));
        assert_eq!(
            result,
            ComplexDyadic::new(Dyadic::new(3, 0), Dyadic::new(5, 0))
        );
    }

    #[test]
    fn test_multiplication_complex_dyadic() {
        let result = ComplexDyadic::new(Dyadic::new(1, 2), Dyadic::new(2, 0))
            * ComplexDyadic::new(Dyadic::new(2, -1), Dyadic::new(3, 0));
        assert_eq!(
            result,
            ComplexDyadic::new(Dyadic::new(-2, 0), Dyadic::new(14, 0))
        );
    }

    #[test]
    fn test_eval() {
        let x = ComplexDyadic::new(Dyadic::new(1, 0), Dyadic::zero());
        let m_seq = vec![x.clone(), x.clone(), x.clone(), x.clone()];
        let f = ComplexFunction::new(
            BoundingSequence::new(vec![Dyadic::new(2, 0), Dyadic::new(2, 0)]),
            ExpansionCoefficients::new(vec![x.clone(), x.clone(), x.clone()]),
        );
        assert_eq!(
            f.eval(ComplexDyadic::new(Dyadic::new(1, 0), Dyadic::zero())),
            ComplexDyadic::new(Dyadic::new(3, 0), Dyadic::zero())
        );
    }

    // Implement test for derrivative of complex functions
}
