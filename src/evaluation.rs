use std::cmp::Ordering;
use std::sync::Arc;

use futures::stream::{self, StreamExt};
use num_cpus;

use crate::corollary_2::{certify_r_hat_rho, calculate_epsilon};
use crate::covering_grids::{covering_grid_bitmap, create_covering_grid, unit_disk_n, unit_disk_radius};
use crate::dyadic::{Dyadic, ComplexDyadic};
use crate::edt::{landau_l_via_edt_from_bitmap};
use crate::holomorphic::{BoundingSequence, ComplexFunction, ExpansionCoefficients};
use crate::psi::{psi_infinity, t_vector, generate_all_words, generate_word};
use crate::plot::{plot_covering_grid, plot_set};

/// calculate_random approximates value of lambda for a "number" number of words on length word_length, on
/// given delta for grid generation and disk_accuracy for a unit disk domain
pub async fn calculate_random(
    number: usize,
    word_length: i32,
    delta: Dyadic,
    disk_accuracy : i32,
) -> (Vec<u8>, f64) {
    let epsilon = delta * Dyadic::new(1, 2) ;
    let m_seq = vec![
        Dyadic::new(1, -1),
        Dyadic::new(1, -1),
        Dyadic::new(1, -2),
        Dyadic::new(1, -3),
        Dyadic::new(1, -3),
    ];
    let m_seq1  = Arc::new(m_seq);
    let t_seq = Arc::new(t_vector(1000));
    let domain = Arc::new(unit_disk_n(disk_accuracy));
    let concurrency = num_cpus::get().max(1);

    let mut results = stream::iter(0..number)
        .map(|_| {
            let domain = Arc::clone(&domain);
            let m_seq  = Arc::clone(&m_seq1);
            let t_seq  = Arc::clone(&t_seq);

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



/// calculate_for_length evaluates lambda_f values for grids size delta using EDT algorithm. For the domain
/// it takes the standard disk with size 1. This is not exactly as the article suggests, but it is a good approximation none the less.
pub async fn calculate_for_length(
    length: usize,
    delta: Dyadic,
    disk_accuracy: i32,
    m_seq : Vec<Dyadic>,
) -> (Vec<u8>, f64) {
    let epsilon = delta * Dyadic::new(1, 2) ;
    
    let m_seq1  = Arc::new(m_seq);
    let t_seq = Arc::new(t_vector(100));
    let domain = Arc::new(unit_disk_n(disk_accuracy));
    let concurrency = num_cpus::get().max(1);
    let words: Vec<Vec<u8>> = generate_all_words(length);

    let mut results = stream::iter(words.into_iter())
        .map(|word| {
            let domain = Arc::clone(&domain);
            let m_seq2  = Arc::clone(&m_seq1);
            let t_seq  = Arc::clone(&t_seq);

            tokio::task::spawn_blocking(move || {
                let coeffs = psi_infinity(&m_seq2, &t_seq, &word);
                let fprime = ComplexFunction::new(
                    BoundingSequence::new((*m_seq2).clone()),
                    ExpansionCoefficients::new(coeffs),
                );
                let f = fprime.antiderivative();

                // Evaluate on domain
                let mut img: Vec<ComplexDyadic> = Vec::with_capacity(domain.len());
                for &z in domain.iter() {
                    img.push(f.eval(&z));
                }

                // standard evaluation using GridBitmap and EDT
                let grid_bitmap = covering_grid_bitmap(&img, epsilon, delta);
                let l_val = landau_l_via_edt_from_bitmap(&grid_bitmap, delta);

                (l_val, word)
            })
        })
        .buffer_unordered(concurrency);

    let mut best_l = f64::INFINITY;
    let mut best_word: Vec<u8> = vec![];

    // we have results for all words, we then find the one with the best approximation
    while let Some(res) = results.next().await {
        let (l_val, word) = res.expect("worker panicked");
        if l_val < best_l {
            best_l = l_val;
            best_word = word;
        }
    }

    (best_word, best_l)

 }

 /**  
This part was only used partially in order to evaluate performance of different m_sequences. It evaluates the calculate_for_length
on all possible m_sequences formed by given dyadics (seq of finite length) and returns ones with best approximations. Currently not in use.
*/
struct MSeqResult {
    m_seq: Vec<Dyadic>,
    word: Vec<u8>,
    approx: f64,
}
/// Saves top_k smallest values in a vector of MSeq, helper function
fn keep_top_k(top: &mut Vec<MSeqResult>, entry: MSeqResult, k: usize) {
    top.push(entry);
    top.sort_by(|a, b| a.approx.partial_cmp(&b.approx).unwrap_or(Ordering::Equal));
    if top.len() > k {
        top.truncate(k);
    }
}

/// Used to check different m sequences, sweeps through all combinations of length 3 vectors, where each value is one of possible_m elements (Dyadics).
pub async fn sweep_mseq_len3(
    possible_m: &[Dyadic],
    length: usize,
    delta: Dyadic,
    disk_radius: i32,
) -> Vec<(Vec<Dyadic>, Vec<u8>, f64)> {
    let n = possible_m.len();
    let mut top3: Vec<MSeqResult> = Vec::new();

    for i in 0..n {
        for j in 0..n {
            for k in 0..n {
                let m_seq = vec![
                    possible_m[i],
                    possible_m[j],
                    possible_m[k],
                ];

                let (best_word, best_l) =
                    calculate_for_length(length, delta, disk_radius, m_seq.clone()).await;

                keep_top_k(
                    &mut top3,
                    MSeqResult { m_seq, word: best_word, approx: best_l },
                    3,
                );
            }
        }
    }

    println!("=== Top 3 m_seq by approx (ascending) ===");
    for (rank, e) in top3.iter().enumerate() {
        let ms = e
            .m_seq
            .iter()
            .map(|d| format!("{:.6}", d.to_f64()))
            .collect::<Vec<_>>()
            .join(", ");
        println!(
            "#{:>2}: m_seq = [{}],  approx = {:.8},  word = {:?}",
            rank + 1,
            ms,
            e.approx,
            e.word
        );
    }

    top3
        .into_iter()
        .map(|e| (e.m_seq, e.word, e.approx))
        .collect()
}
// we just use this to calculate disk_accuracy needed, epsilon/2 is fine as points arent spread out too far apart
fn minimal_pow2_step_below(eps: Dyadic) -> (i32, Dyadic) {
    assert!(eps > Dyadic::zero(), "eps must be > 0");
    
    let mut n = (-eps.to_f64().log2()).ceil() as i32;
    // assure strict inequality
    let mut step = Dyadic::new(1, -n).reduce(); // 2^{-n}
    if step >= eps {
        n += 1;
        step = Dyadic::new(1, -n).reduce();
    }
    (-n-1, step)
}


/// calculates the approximation for a certain word on a disk with radius 1-2^(-n), as described in corollary 2. We currently use this
/// this one only evaluates one word, disk_decrease tells us what radius we use for domain. We can also choose to plot results. We do that
/// inside of the function because visualisation helps us understand if the sets are being evaluated properly.
pub async fn calculate_for_word_updated(word : Vec<u8>, disk_decrease : i32, m_seq : Vec<Dyadic>, plot : bool) -> f64 {
    let r_dy = (Dyadic::new(1, 0) - Dyadic::new(1, disk_decrease)).reduce();
    let r_f = r_dy.to_f64() ; 
    let m_seq1  = Arc::new(m_seq);
    let t_seq = Arc::new(t_vector(100)); // Could be bigger, but for now 100 is MORE than enough
    let n_samples = 500 ; //subject to change, used to evaluate rho as described in Lemma 6
    tokio::task::spawn_blocking(move || {
        // Build f' from the word, then integrate to f
        let coeffs = psi_infinity(&m_seq1, &t_seq, &word);
        let fprime = ComplexFunction::new(
            BoundingSequence::new((*m_seq1).clone()),
            ExpansionCoefficients::new(coeffs),
        );
        let f = fprime.antiderivative();
        
        // we calculate r_hat, rho as described in article. This later gives us epsilon, so we can complute sufficient epsilon covering grid
        let (r_hat, rho) = certify_r_hat_rho(r_dy, &fprime, n_samples, 30, 30) ;
        let mut eps = calculate_epsilon(r_f, r_hat, rho) * Dyadic::new(1, 2);
        let min_eps = Dyadic::new(1, -8); // a custom fix for now
        if eps < min_eps {
            eps = min_eps;
        }
        assert!(eps > Dyadic::zero(), "ε came out zero; cannot proceed");
        let delta = (eps * Dyadic::new(1, -2)).reduce();
        let (acc_n, _step) = minimal_pow2_step_below(eps);
        let domain = unit_disk_radius(acc_n, r_f) ;
        let mut img: Vec<ComplexDyadic> = Vec::with_capacity(domain.len());
        for &z in domain.iter() {
            img.push(f.eval(&z));
        } 
        let grid_bitmap = covering_grid_bitmap(&img, eps, delta);
        // println!("{}", grid_bitmap.height * grid_bitmap.width);
        let l_val = landau_l_via_edt_from_bitmap(&grid_bitmap, delta);
        if plot {
            println!("{}, {} is eps, delta", eps.to_f64(), delta.to_f64());
            plot_set(&img, "test_image.png");
            let grid_set = create_covering_grid(&img, eps, delta);
            plot_covering_grid(&grid_set, "test_grid.png");
        }
        l_val
    }) 
    .await
    .expect("worker panicked")
}


/// approximates Landau's constant as minimum of all values on words specific lenght, evaluated on domain of disk with radius 1 - disk_decrease
pub async fn approximate_all_words_length(
    length: usize,
    disk_decrease: i32,
    m_seq: Vec<Dyadic>
) -> (f64, Vec<u8>) {
    let words: Vec<Vec<u8>> = generate_all_words(length);

    let mut best_l = f64::INFINITY;
    let mut best_word: Vec<u8> = Vec::new();

    for w in &words {
        let l = calculate_for_word_updated(w.clone(), disk_decrease, m_seq.clone(), false).await;
        if l < best_l {
            best_l = l;
            best_word = w.clone();
        }
    }
    (best_l, best_word)
}

/// approximates Landau's constant as minimum of all values on words with length up or equal to max_length, evaluated on domain of disk with radius 1 - disk_decrease (step)
pub async fn approximate_all_words(
    max_length: usize,
    disk_decrease: i32,
    m_seq: Vec<Dyadic>,
) -> f64 {
    let parallelism = std::thread::available_parallelism()
        .map(|n| n.get())
        .unwrap_or(4);

    let m_seq = Arc::new(m_seq);

    let mut best_l = f64::INFINITY;
    let mut best_word: Vec<u8> = Vec::new();

    for n in 1..=max_length {
        let words: Vec<Vec<u8>> = generate_all_words(n);

        let mut len_best = f64::INFINITY;
        let mut len_best_word: Vec<u8> = Vec::new();

        // Kick off all evaluations for this length, but bound how many run at once.
        let mut futs = stream::iter(words.into_iter().map(|w| {
            let m_seq = m_seq.clone();
            async move {
                let l = calculate_for_word_updated(
                    w.clone(),
                    disk_decrease,
                    (*m_seq).clone(),
                    /*plot:*/ false,
                ).await;
                (w, l)
            }
        }))
        .buffer_unordered(parallelism);

        // Reduce to the best word for this length
        while let Some((w, l)) = futs.next().await {
            if l < len_best {
                len_best = l;
                len_best_word = w;
            }
        }

        // Update global best
        if len_best < best_l {
            best_l = len_best;
            best_word = len_best_word;
        }
    }

    // Plot the global best (same as your original)
    let _ = calculate_for_word_updated(
        best_word.clone(),
        disk_decrease,
        (*m_seq).clone(),
        /*plot:*/ true,
    ).await;

    println!("Best word: {:?}\nMin approx = {}", best_word, best_l);
    best_l
}






