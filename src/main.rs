use std::sync::Arc;
use std::vec;

use futures::stream::{self, StreamExt};
use landau::corollary_2::{calculate_epsilon, certify_r_hat_rho};
use landau::covering_grids::{covering_grid_bitmap, create_covering_grid, unit_disk_n};
use landau::dyadic::Dyadic;
use landau::edt::{landau_l_via_edt_from_bitmap, landau_l_with_edt_from_complex_set, print_grid};
use landau::evaluation::{
    calculate_for_all_words_updated,
    calculate_for_length,
    calculate_for_word_updated,
    sweep_mseq_len3,
};
use landau::holomorphic::{BoundingSequence, ComplexFunction, ExpansionCoefficients};
use landau::plot::{plot_covering_grid, plot_set};
use landau::psi::{generate_all_words, generate_word, m_vec, mu_second, psi_infinity, t_vector};


pub async fn calculate_random(
    number: usize,
    word_length: i32,
    delta: Dyadic,
    disk_radius: i32,
) -> (Vec<u8>, f64) {
    let epsilon = delta * Dyadic::new(1, 2);
    let m_seq = vec![
        Dyadic::new(1, -1),
        Dyadic::new(1, -1),
        Dyadic::new(1, -2),
        Dyadic::new(1, -3),
        Dyadic::new(1, -3),
    ];
    let m_seq1 = Arc::new(m_seq);
    let t_seq = Arc::new(t_vector(1000));
    let domain = Arc::new(unit_disk_n(disk_radius));
    let concurrency = num_cpus::get().max(1);
    println!("Working on {} CPU cores", concurrency);

    let mut results = stream::iter(0..number)
        .map(|_| {
            let domain = Arc::clone(&domain);
            let m_seq = Arc::clone(&m_seq1);
            let t_seq = Arc::clone(&t_seq);

            // Heavy CPU work => spawn_blocking keeps the async scheduler happy
            tokio::task::spawn_blocking(move || {
                // --- build a random function ---
                let word = generate_word(word_length);
                let coeffs = psi_infinity(&m_seq, &t_seq, &word);

                let fprime = ComplexFunction::new(
                    BoundingSequence::new((*m_seq).clone()),
                    ExpansionCoefficients::new(coeffs),
                );
                let f = fprime.antiderivative();

                // --- eval on domain ---
                let mut img = Vec::with_capacity(domain.len());
                for &z in domain.iter() {
                    img.push(f.eval(&z));
                }

                // --- EDT path you asked for ---
                let grid_bitmap = covering_grid_bitmap(&img, epsilon, delta);
                let l_val = landau_l_via_edt_from_bitmap(&grid_bitmap, delta);

                (l_val, word)
            })
        })
        .buffer_unordered(concurrency);
    let mut best_l = f64::INFINITY;
    let mut best_word: Vec<u8> = vec![];

    while let Some(res) = results.next().await {
        let (l_val, word) = res.expect("worker panicked");
        if l_val < best_l {
            best_l = l_val;
            best_word = word;
        }
    }
    (best_word, best_l)
}



#[tokio::main(flavor = "multi_thread")]
async fn main() -> anyhow::Result<()> {
    landau::ui::run_server().await
}
// async fn main() -> Result<(), Box<dyn std::error::Error>> {

//     let delta      = Dyadic::new(1, -7);
//     let m_seq = vec![
//         Dyadic::new(1, -2),
//         Dyadic::new(1, -2),
//         Dyadic::new(1, -2)
//     ];
//     // let (best_word_random, approx_random) = calculate_random(1000, 30, delta, -8).await ;
//     // println!("best approx for random words is: {}, best word is : {:?}", approx_random, best_word_random) ;
//     let (best_word_length, approx_length) = calculate_for_length(3, delta, -5, m_seq).await ;
//     println!("best approx for all words with this length is: {}, best word is : {:?}", approx_length, best_word_length) ;
//     let t_seq1 = t_vector(100);
//     let m_seq1 = vec![
//         Dyadic::new(1, -2),
//         Dyadic::new(1, -2),
//         Dyadic::new(1, -2),
//         Dyadic::new(1, -3),
//         Dyadic::new(1, -3),
//     ];
//     let m_seq2 = m_vec(5) ;
//     let coeffs = psi_infinity(&m_seq2, &t_seq1, &best_word_length);
//     let fprime = ComplexFunction::new(
//         BoundingSequence::new(m_seq1.clone()),
//         ExpansionCoefficients::new(coeffs.clone()),
//     );
//     let f = fprime.antiderivative();
//     for z in &f.expansion_coefficients.vector {
//         println!("{}+{}i", z.re.to_f64(), z.im.to_f64())
//     }
//     let mut img = vec![] ;
//     for &z in &unit_disk_n(-7) {
//         img.push(f.eval(&z));
//     }

//     plot_set(&img, "min_image.png") ;
//     let grid = create_covering_grid(&img, delta * Dyadic::new(1, 2), delta);
//     plot_covering_grid(&grid, "min_grid.png");

//     let possible_m = vec![
//         Dyadic::new(3, -1),
//         Dyadic::new(1, 0),
//     ];
//     // let _top3 = sweep_mseq_len3(&possible_m, 6, delta, -7).await;
//     let val = calculate_for_all_words_updated(3, -1, m_seq2).await;
//     println!("{}", val) ;

//     println!("{}", mu_second(&1.0, &Dyadic::new(7, -3)));
//     Ok(())
// }
// Test to ensure basic operation and functions perform properly. More test can be added as needed.
#[cfg(test)]
mod tests {
    use landau::dyadic::ComplexDyadic;

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
            f.eval(&ComplexDyadic::new(
                Dyadic::new(1, 0),
                Dyadic::zero()
            )),
            ComplexDyadic::new(Dyadic::new(3, 0), Dyadic::zero())
        );
    }
    #[test]
    fn test_epsilon() {
        let r = 0.95;
        let r_hat = Dyadic::new(127, -7);
        println!("{}", r_hat.to_f64());
        let rho = Dyadic::new(1, -4);
        let eps = calculate_epsilon(r, r_hat, rho);
        assert_eq!(eps, Dyadic::new(1, -37));
    }
}
