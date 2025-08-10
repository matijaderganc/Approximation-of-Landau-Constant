use landau::covering_grids::{
    create_covering_grid, grid_approx_with_edge, unit_disk_n, unit_disk_n_boundary
};
use landau::dyadic::{ComplexDyadic, Dyadic};
use landau::plot::{plot_covering_grid, plot_set};

use landau::holomorphic::{
    BoundingSequence, ComplexFunction, ExpansionCoefficients,
};
use landau::psi::{generate_all_words, psi_infinity, t_vector};
use std::vec;

fn evaluate_function(f : &Vec<ComplexDyadic>, delta : Dyadic, image_acc : i32) -> f64 {
    let m_seq = vec![Dyadic::new(1, -1), Dyadic::new(1, -1), Dyadic::new(1, -1)] ; // no role for now
    let complex = ComplexFunction::new(
        BoundingSequence::new(m_seq),
        ExpansionCoefficients::new(f.clone()),
    );
    let integral: ComplexFunction = complex.antiderivative(); 

    let inside = unit_disk_n(image_acc) ;
    
    let edge = unit_disk_n_boundary(image_acc) ; 
    let mut image_inside: Vec<ComplexDyadic> = Vec::new() ;
    let mut image_outside: Vec<ComplexDyadic> = Vec::new() ;
    

    for x in &inside {
        image_inside.push(integral.eval(*x))
    }
    for y in &edge {
        image_outside.push(integral.eval(*y))
    }
    let grid_inside = create_covering_grid(&image_inside, delta.clone()* Dyadic::new(4, 0), delta.clone()) ;
    let grid_outside = create_covering_grid(&image_outside, delta.clone()* Dyadic::new(4, 0), delta.clone()) ;
    grid_approx_with_edge(&grid_inside, &grid_outside,  delta)
} 
fn main() -> Result<(), Box<dyn std::error::Error>>
{   
    

    let m_seq = vec![Dyadic::new(1, 0), Dyadic::new(1, -1), Dyadic::new(1, -1), Dyadic::new(1, -1)] ;
    let t_seq = t_vector(20) ;
    let all_6 = generate_all_words(6) ;

    let mut min_word = vec![1, 1, 3, 1, 1, 1] ;
    let mut min = 5.0 ;
    let holo1 = psi_infinity(&m_seq, &t_seq, &min_word) ;
    let val = evaluate_function(&holo1, Dyadic::new(1, -5), -5);
    println!("\n this is min approx {}", val) ;
    
    // for word in all_6 { 
    //     let m_seq = vec![Dyadic::new(1, 0), Dyadic::new(1, -1), Dyadic::new(1, -1), Dyadic::new(1, -2)] ;
    //     let t_seq = vec![0, 1, 0, 1, 2, 0, 1, 2, 2, 0, 1, 2, 3, 0, 1, 2, 3, 0, 1, 2, 3, 0, 1, 2, 3];
    //     let holo1 = psi_infinity(&m_seq, &t_seq, &word) ;
    //     let val = evaluate_function(&holo1, Dyadic::new(1, -5), -5);
        
    //     if  val < min {
    //         min = val ;
    //         min_word = word
    //     }
    // } 
    // print!(" this is min word{:?}", min_word) ;
    // println!("\n this is min approx {}", min) ;
   




    let m_seq1 = vec![Dyadic::new(1, 0), Dyadic::new(1, -1), Dyadic::new(1, -1), Dyadic::new(1, -1)] ;
    let t_seq1 = vec![0, 1, 0, 1, 2, 0, 1, 2, 2, 0, 1, 2, 3, 0, 1, 2, 3, 0];
    let min_func = psi_infinity(&m_seq1, &t_seq1, &min_word) ;

    let test_func = vec![ComplexDyadic::one(), ComplexDyadic::zero(), ComplexDyadic::zero(),ComplexDyadic::new(Dyadic::new(81, -5), Dyadic::zero()), ComplexDyadic::new(Dyadic::new(-81 * 5, -7), Dyadic::zero())];
    let disk1 = unit_disk_n(-5) ;
    let disk_boundary1 = unit_disk_n_boundary(-5) ;
    let min_complex = ComplexFunction::new(
        BoundingSequence::new(m_seq1),
        ExpansionCoefficients::new(test_func.clone()),
    );
    let mut min_im: Vec<ComplexDyadic> = Vec::new() ;

    let min_integral = min_complex.antiderivative();
    for x in &disk1 {
        min_im.push(min_integral.eval(*x))
    }    ;
    plot_set(&min_im, "min_image.png") ;
    let min_grid = create_covering_grid(&min_im, Dyadic::new(1, -4), Dyadic::new(1, -5)) ;
    plot_covering_grid(&min_grid, "min_grid.png") ;
    Ok(())
}


// let mut min_word: Vec<_> = vec![] ;
    // let mut min = 5.0 ;
    // let disk1 = unit_disk_n(-3) ;
    // let disk_boundary = unit_disk_n_boundary(-5) ;
    // plot_set(&disk_boundary, "boundary.png") ;

    // for i in 0..100 {
    //     let m_seq = vec![Dyadic::new(1, 0), Dyadic::new(1, -1), Dyadic::new(1, -1), Dyadic::new(1, -2)] ;
    //     let t_seq = vec![0, 1, 0, 1, 2, 0, 1, 2, 2, 0, 1, 2, 3, 0, 1, 2, 3, 0, 1, 2, 3, 0, 1, 2, 3];
    //     let word = generate_word(30);
    //     let holo1 = psi_infinity(&m_seq, &t_seq, &word) ;
    //     let val = evaluate_function(&holo1, Dyadic::new(1, -3), &disk1);
        
    //     if  val < min {
    //         min = val ;
    //         min_word = word
    //     }
    // }
    // print!(" this is min word{:?}", min_word) ;
    // println!("\n this is min approx {}", min) ;
    // min = 5.0 ;

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
