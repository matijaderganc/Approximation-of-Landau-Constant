use landau::covering_grids::{
    create_covering_grid, grid_approx_with_edge, unit_disk_n, unit_disk_n_boundary
};
use landau::dyadic::{ComplexDyadic, Dyadic};
use landau::plot::{plot_covering_grid, plot_set};

use landau::holomorphic::{
    BoundingSequence, ComplexFunction, ExpansionCoefficients,
};
use landau::psi::{generate_all_words, psi_infinity, t_vector, m_n, m_vec};
use landau::edt::{edt_1d_squared, edt_2d_squared, print_grid, landau_l_with_edt_from_complex_set} ;

use std::vec;

const INF: f64 = 1e20;


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
    let eps_over_4 = Dyadic::new(1, -3); // ε/4 = 2^-3
    let delta      = Dyadic::new(1, -5); // δ   = 2^-5  (δ ≤ ε/4)
    let radius_n   = -4;                 // sample D_{1-2^{-6}}

    // --- sequences used consistently everywhere ---
    // Use the same m_seq & t_seq for the search and the final plot.
    // (Here I keep your small hand-picked m_seq; swap for m_vec(N) if you want.)
    let m_seq = m_vec(4) ;
    let t_seq: Vec<usize> = t_vector(64); // dispatcher; make it long enough

    // --- precompute the domain sample once ---
    let domain = unit_disk_n(radius_n);

    // --- search over words of length 6 ---
    let mut best_word: Vec<u8> = vec![1, 1, 3, 1];
    let mut best_l = f64::INFINITY;

    // Reuse buffers across iterations
    let mut img: Vec<ComplexDyadic> = Vec::with_capacity(domain.len());

    for word in generate_all_words(8) {
        // coefficients from ψ∞
        let coeffs = psi_infinity(&m_seq, &t_seq, &word);

        // build f from f' (coeffs are for f'), then integrate
        let fprime = ComplexFunction::new(
            BoundingSequence::new(m_seq.clone()),
            ExpansionCoefficients::new(coeffs),
        );
        let f = fprime.antiderivative();

        // evaluate on the domain
        img.clear();
        for &z in &domain {
            img.push(f.eval(z));
        }

        // covering grid of the image (ε/4, δ)  — ensure this matches your function’s arg order
        let grid = create_covering_grid(&img, eps_over_4, delta);

        // linear-time Landau l via EDT
        let l_val = landau_l_with_edt_from_complex_set(&grid, delta);
        if l_val < best_l {
            best_l = l_val;
            best_word = word; // `word` is moved from the iterator; fine
        }
    }

    println!("min word = {:?}\nmin approx = {}", best_word, best_l);

    // --- render/plot the best one ---
    let coeffs = psi_infinity(&m_seq, &t_seq, &best_word);
    let fprime = ComplexFunction::new(
        BoundingSequence::new(m_seq.clone()),
        ExpansionCoefficients::new(coeffs.clone()),
    );
    let f = fprime.antiderivative();

    img.clear();
    for &z in &domain {
        img.push(f.eval(z));
    }
    let grid = create_covering_grid(&img, eps_over_4, delta);
    plot_covering_grid(&grid, "min_grid.png");

    let l_check = landau_l_with_edt_from_complex_set(&grid, delta);
    println!("l(best) = {}", l_check);

    Ok(())
}

    // for z in &holo1 {
    //     println!("{:?}", z)
    // }
    // let mut min = 5.0 ;
   

    // let min_func = psi_infinity(&m_seq, &t_seq, &min_word) ;

    // // let val = evaluate_function(&holo1, Dyadic::new(1, -4), -5);
    // // println!("\n this is min approx {}", val) ;

    

    // let disk1 = unit_disk_n(-6) ;
    // for x in &disk1 {
    //     min_im.push(min_integral.eval(*x))
    // }    ;
    // plot_set(&min_im, "min_image.png") ;
    // let min_grid = create_covering_grid(&min_im, Dyadic::new(1, -2), Dyadic::new(1, -4)) ;
    // println!("{}", min_grid.len()) ;

    // plot_covering_grid(&min_grid, "min_grid.png") ;


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



