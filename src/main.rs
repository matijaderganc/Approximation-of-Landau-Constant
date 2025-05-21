use landau::covering_grids::{
    self, create_covering_grid, extreme_points, unit_disk_n, grid_complement, grid_approx
};
use landau::dyadic::{psi, ComplexDyadic, Dyadic, Interval};
use landau::holomorphic::{
    comp_vec_to_sequence, vec_to_sequence, BoundingSequence, ComplexFunction, ExpansionCoefficients,
};
use landau::psi::{mu_first, mu_second, psi_infinity};
use landau::plot::{self, plot_set, plot_covering_grid} ;

use plotters::prelude::*;
use std::collections::LinkedList;

fn evaluate_function(f : &Vec<ComplexDyadic>, delta : Dyadic, disk : &Vec<ComplexDyadic>) -> Option<f64> {
    let m_seq = vec![Dyadic::new(1, -1), Dyadic::new(1, -1)] ; // no role for now
    let complex = ComplexFunction::new(
        BoundingSequence::new(m_seq),
        ExpansionCoefficients::new(f.clone()),
    );
    let mut image = Vec::new() ;
    for x in disk {
        image.push(complex.eval(*x))
    }
    let grid = create_covering_grid(&image, delta.clone()* Dyadic::new(4, 0), delta.clone()) ;
    grid_approx(&grid, delta)
} 
fn main()-> Result<(), Box<dyn std::error::Error>>
{   
    let x = Dyadic::new(1, -1) ;
    let m_seq = vec![x.clone(), x.clone()] ;
    let t_seq = vec![0, 1, 0, 1, 2, 0, 1, 2, 2, 0, 1, 2, 3, 0, 1, 2, 3, 0];
    let word = vec![1, 1, 1, 3, 2, 1, 3, 1, 2, 3, 1, 2, 3, 1, 3, 2, 1, 2, 3, 2];
    let holo = psi_infinity(&m_seq, &t_seq, &word);
    for t in &holo{
        println!("{}", t)
    }
    println!("done") ;
    let points1 = unit_disk_n(-6);

    let f = ComplexFunction::new(
        BoundingSequence::new(m_seq),
        ExpansionCoefficients::new(holo),
    );
    let mut points2 = Vec::new() ;
    for x in &points1 {
        points2.push(f.eval(*x))
    }    ;
    let grid = create_covering_grid(&points2, Dyadic::new(1, -2), Dyadic::new(1, -4)) ;
    plot_set(&points2, "image2.png") ;
    plot_covering_grid(&grid, "covering_grid.png")?;
    plot_covering_grid(&grid_complement(&grid, Dyadic::new(1, -4)), "complement.png")?;
    println!("approximation is {}", grid_approx(&grid, Dyadic::new(1, -4)).unwrap()) ;
    let m_seq1 = vec![x.clone(), x.clone()] ;

    let holo1 = psi_infinity(&m_seq1, &t_seq, &word);

    println!("new approx {}", evaluate_function(&holo1, Dyadic::new(1, -4), &points1).unwrap()) ;
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
