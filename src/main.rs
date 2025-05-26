use landau::covering_grids::{
    self, create_covering_grid, extreme_points, unit_disk_n, grid_complement, grid_approx
};
use landau::dyadic::{psi, ComplexDyadic, Dyadic, Interval};
use landau::plot::{self, plot_covering_grid, plot_set};

use landau::holomorphic::{
    comp_vec_to_sequence, vec_to_sequence, BoundingSequence, ComplexFunction, ExpansionCoefficients,
};
use landau::psi::{generate_all_words, generate_word, mu_first, mu_second, psi_infinity};
use plotters::prelude::*;

use std::collections::LinkedList;
use std::f64::consts::E;
use std::vec;

fn evaluate_function(f : &Vec<ComplexDyadic>, delta : Dyadic, disk : &Vec<ComplexDyadic>) -> f64 {
    let m_seq = vec![Dyadic::new(1, -1), Dyadic::new(1, -1), Dyadic::new(1, -1)] ; // no role for now
    let complex = ComplexFunction::new(
        BoundingSequence::new(m_seq),
        ExpansionCoefficients::new(f.clone()),
    );
    let integral: ComplexFunction = complex.antiderivative(); 
    let mut image = Vec::new() ;
    for x in disk {
        image.push(integral.eval(*x))
    }
    let grid = create_covering_grid(&image, delta.clone()* Dyadic::new(4, 0), delta.clone()) ;
    grid_approx(&grid, delta)
} 
fn main()-> Result<(), Box<dyn std::error::Error>>
{   
    let mut min_word = vec![] ;
    let mut min = 5.0 ;
    let disk1 = unit_disk_n(-5) ;

    for i in 0..150 {
    let m_seq = vec![Dyadic::new(1, 0), Dyadic::new(1, -1), Dyadic::new(1, -1), Dyadic::new(1, -2)] ;
    let t_seq = vec![0, 1, 0, 1, 2, 0, 1, 2, 2, 0, 1, 2, 3, 0, 1, 2, 3, 0, 1, 2, 3, 0, 1, 2, 3];
    let word = generate_word(30);
    let holo1 = psi_infinity(&m_seq, &t_seq, &word) ;
    let val = evaluate_function(&holo1, Dyadic::new(1, -3), &disk1);
    
    if  val < min {
        min = val ;
        min_word = word
    }
}
    print!(" this is min word{:?}", min_word) ;
    println!("\n this is min approx {}", min) ;
    let m_seq1 = vec![Dyadic::new(1, 0), Dyadic::new(1, -1), Dyadic::new(1, -1), Dyadic::new(1, -2)] ;
    let t_seq1 = vec![0, 1, 0, 1, 2, 0, 1, 2, 2, 0, 1, 2, 3, 0, 1, 2, 3, 0];
    let min_func = psi_infinity(&m_seq1, &t_seq1, &min_word) ;
    for z in min_func.clone() {
        println!("{z}")
    }
    println!("done with vector") ;
    let mut min_im: Vec<_> = Vec::new() ;


    let disk = unit_disk_n(-5) ;
    let min_complex = ComplexFunction::new(
        BoundingSequence::new(m_seq1),
        ExpansionCoefficients::new(min_func.clone()),
    );
    let min_integral = min_complex.antiderivative();
    for x in &disk {
        min_im.push(min_integral.eval(*x))
    }    ;
    plot_set(&min_im, "min_image.png") ;
    let min_grid = create_covering_grid(&min_im, Dyadic::new(1, -3), Dyadic::new(1, -5)) ;
    plot_covering_grid(&min_grid, "min_grid.png") ;
    
    Ok(())
}
// let min_func: Vec<ComplexDyadic> = psi_infinity(&m_seq1, &t_seq1, &word_min) ;
    // let min_func: Vec<ComplexDyadic> = vec![ComplexDyadic::zero(), ComplexDyadic::one(), ComplexDyadic::zero(), ComplexDyadic::new(Dyadic::new(-1, -2), Dyadic::zero()), ComplexDyadic::zero(), ComplexDyadic::new(Dyadic::new(-1, -3), Dyadic::zero()), ComplexDyadic::zero(), ComplexDyadic::new(Dyadic::new(-1, -3), Dyadic::zero())] ;
    
// println!("{}", evaluate_function(&min_func, Dyadic::new(1, -3), &disk).unwrap()) ;
    // let x: ComplexDyadic = ComplexDyadic::new(Dyadic::new(1, -2), Dyadic::new(1, -2)) ;

// println!("{}, {}", x, min_complex.eval(x)) ;
    // println!("done evaluating") ;
    // for c in min_complex.expansion_coefficients.vector {
    //         println!("{}", c)
    //     }
        
    // println!("done with first complex function") ;
    // for x in min_integral.expansion_coefficients.vector {
    //     println!("{}", x)
    // }
    // let value = evaluate_function(&min_func, Dyadic::new(1, -5), &disk) ;
    // println!("{}", value) ;
// // add these to tests!!!

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
