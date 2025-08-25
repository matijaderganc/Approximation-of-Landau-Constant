use std::collections::LinkedList;

use landau::covering_grids::{
    covering_grid_bitmap, create_covering_grid, extreme_points, grid_approx_with_edge, unit_disk_n, unit_disk_n_boundary
};
use landau::dyadic::{ComplexDyadic, Dyadic};
use landau::plot::{plot_covering_grid, plot_set};

use landau::holomorphic::{
    BoundingSequence, ComplexFunction, ExpansionCoefficients,
};
use landau::psi::{generate_all_words, psi_infinity, t_vector, m_vec, generate_word};
use landau::edt::{edt_1d_squared, edt_2d_squared, landau_l_via_edt_from_bitmap, landau_l_with_edt_from_complex_set, print_grid} ;

use std::vec;
use std::time::{Duration, Instant};
fn ms(d: Duration) -> f64 { d.as_secs_f64() * 1000.0 }


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
fn main() -> Result<(), Box<dyn std::error::Error>> {
    // --- parameters (tweak if you like) ---
    let eps_over_4 = Dyadic::new(1, -5); // ε/4 = 2^-3
    let delta      = Dyadic::new(1, -7); // δ   = 2^-5  (δ ≤ ε/4)
    let radius_n   = -6;                 // sample D_{1-2^{-6}}
    // --- sequences used consistently everywhere ---
    // Use the same m_seq & t_seq for the search and the final plot.
    // (Here I keep your small hand-picked m_seq; swap for m_vec(N) if you want.)
    // let m_seq = vec![Dyadic::new(1, -1), Dyadic::new(1, -1), Dyadic::new(1, -2), Dyadic::new(1, -3), Dyadic::new(1, -3)];
    let m_seq = m_vec(15) ;
    let t_seq: Vec<usize> = t_vector(64); // dispatcher; make it long enough

    // --- precompute the domain sample once ---
    let domain = unit_disk_n(radius_n);
    

    // --- search over words of length 6 ---
    let mut best_word: Vec<u8> = vec![];
    let mut best_l = 100.0;

    // Reuse buffers across iterations
    let mut img: Vec<ComplexDyadic> = Vec::with_capacity(domain.len());
    //let all_8 = generate_all_words(1) ;
    for i  in 0..10 {
        // let t1 = Instant::now();

        let word = generate_word(30);
        // coefficients from ψ∞
        let coeffs = psi_infinity(&m_seq, &t_seq, &word);

        // build f from f' (coeffs are for f'), then integrate
        let fprime = ComplexFunction::new(
            BoundingSequence::new(m_seq.clone()),
            ExpansionCoefficients::new(coeffs),
        );
        let f = fprime.antiderivative();
        // let dt_1 = t1.elapsed() ;
        // println!("time spent creating func: {}", ms(dt_1)) ;
        // evaluate on the domain
        // let t2: Instant = Instant::now();

        img.clear();
        for &z in &domain {
            img.push(f.eval(z));
        }
        // let dt_2 = t2.elapsed() ;
        // println!("time spent evaluating: {}", ms(dt_2)) ;
        // covering grid of the image (ε/4, δ)  — ensure this matches your function’s arg order
        // let t3: Instant = Instant::now();

        // let grid = create_covering_grid(&img, eps_over_4, delta);
        // let dt_3 = t3.elapsed() ;
        // println!("time spent creating grid: {}", ms(dt_3)) ;

        // let t5: Instant = Instant::now();
        let grid_bitmap = covering_grid_bitmap(&img, eps_over_4, delta);
        // let dt_5 = t5.elapsed() ;
        // println!("time spent creating grid: {}", ms(dt_5)) ;
        // linear-time Landau l via EDT
        // let t4: Instant = Instant::now();
        let l_val = landau_l_via_edt_from_bitmap(&grid_bitmap, delta);
        // let l_val = landau_l_with_edt_from_complex_set(&grid, delta);
        // let dt_4 = t4.elapsed() ;
        // println!("time spent edt: {}", ms(dt_4)) ;
        if l_val < best_l {
            best_l = l_val;
            best_word = word; // `word` is moved from the iterator; fine
        }
    }

    println!("min word = {:?}\nmin approx = {}", best_word, best_l);
    // best_word = vec![1, 1, 1, 1, 1, 1] ;

    // --- render/plot the best one ---
    let coeffs = psi_infinity(&m_seq, &t_seq, &best_word);
    let fprime = ComplexFunction::new(
        BoundingSequence::new(m_seq.clone()),
        ExpansionCoefficients::new(coeffs.clone()),
    );
    let f = fprime.antiderivative();
    for z in &f.expansion_coefficients.vector {
        println!("{}+{}i", z.re.to_f64(), z.im.to_f64())
    }
    img.clear();
    for &z in &domain {
        img.push(f.eval(z));
    }
    // for &z in &img {
    //     println!("{:?}", z) 
    // }
    let grid = create_covering_grid(&img, eps_over_4, delta);
    plot_covering_grid(&grid, "min_grid.png");
    let ext = extreme_points(&img) ;
    println!("{:?}", ext) ;
    plot_set(&img, "min_image.png") ;
    let l_check = landau_l_with_edt_from_complex_set(&grid, delta);
    println!("l(best) = {}", l_check);

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


// Test to ensure basic operation and functions perform properly. More test can be added as needed.
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
            f.eval(ComplexDyadic::new(
                Dyadic::new(1, 0),
                Dyadic::zero()
            )),
            ComplexDyadic::new(Dyadic::new(3, 0), Dyadic::zero())
        );
    }
}
