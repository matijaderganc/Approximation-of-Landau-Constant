use std::cmp::Ordering;
use std::sync::Arc;

use futures::stream::{self, StreamExt};
use num_cpus;

use crate::corollary_2::{certify_r_hat_rho, calculate_epsilon};
use crate::covering_grids::{covering_grid_bitmap, unit_disk_n, unit_disk_radius, create_covering_grid};
use crate::dyadic::{Dyadic, ComplexDyadic};
use crate::edt::landau_l_via_edt_from_bitmap;
use crate::holomorphic::{BoundingSequence, ComplexFunction, ExpansionCoefficients};
use crate::psi::{psi_infinity, t_vector, generate_all_words};
use crate::plot::{plot_covering_grid, plot_set};

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
    // println!("Working on {} CPU cores", concurrency) ;
    let words: Vec<Vec<u8>> = generate_all_words(length);

    let mut results = stream::iter(words.into_iter())
        .map(|word| {
            let domain = Arc::clone(&domain);
            let m_seq2  = Arc::clone(&m_seq1);
            let t_seq  = Arc::clone(&t_seq);

            tokio::task::spawn_blocking(move || {
                // Build f' from the word, then integrate to f
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

                // EDT path (bitmap → Landau l)
                let grid_bitmap = covering_grid_bitmap(&img, epsilon, delta);
                let l_val = landau_l_via_edt_from_bitmap(&grid_bitmap, delta);

                (l_val, word)
            })
        })
        .buffer_unordered(concurrency);

    // Reduce to the best (min ℓ)
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

struct MSeqResult {
    m_seq: Vec<Dyadic>,
    word: Vec<u8>,
    approx: f64,
}

fn keep_top_k(top: &mut Vec<MSeqResult>, entry: MSeqResult, k: usize) {
    top.push(entry);
    top.sort_by(|a, b| a.approx.partial_cmp(&b.approx).unwrap_or(Ordering::Equal));
    if top.len() > k {
        top.truncate(k);
    }
}

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
                    possible_m[i].clone(),
                    possible_m[j].clone(),
                    possible_m[k].clone(),
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
// we just use this to calculate disk_accuracy
fn minimal_pow2_step_below(eps: Dyadic) -> (i32, Dyadic) {
    assert!(eps > Dyadic::zero(), "eps must be > 0");
    // Start with n = ceil(-log2(eps)) using f64 to get a candidate…
    let mut n = (-eps.to_f64().log2()).ceil() as i32;
    // …and then check/adjust in exact dyadics to ensure strict inequality.
    let mut step = Dyadic::new(1, -n).reduce(); // 2^{-n}
    if step >= eps {
        n += 1;
        step = Dyadic::new(1, -n).reduce();
    }
    (-n-1, step)
}
// calculates the approximation for a certain word on a disk with radius 1-2^(-n), as described in corollary 2.
pub async fn calculate_for_word_updated(word : Vec<u8>, disk_decrease : i32, m_seq : Vec<Dyadic>, plot : bool) -> f64 {
    let r_dy = (Dyadic::new(1, 0) - Dyadic::new(1, disk_decrease)).reduce();
    let r_f = r_dy.to_f64() ; 
    let m_seq1  = Arc::new(m_seq);
    let t_seq = Arc::new(t_vector(100));
    let n_samples = 500 ; //subject to change!
    tokio::task::spawn_blocking(move || {
        // Build f' from the word, then integrate to f
        let coeffs = psi_infinity(&m_seq1, &t_seq, &word);
        let fprime = ComplexFunction::new(
            BoundingSequence::new((*m_seq1).clone()),
            ExpansionCoefficients::new(coeffs),
        );
        let f = fprime.antiderivative();;
        
        let (r_hat, rho) = certify_r_hat_rho(r_dy, &fprime, n_samples, 30, 30) ;
        let mut eps = calculate_epsilon(r_f, r_hat, rho);
        let min_eps = Dyadic::new(1, -6); // THIS IS CUSTOM FIX!!!
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
        let l_val = landau_l_via_edt_from_bitmap(&grid_bitmap, delta);
        if plot {
            plot_set(&img, "test_image.png");
            let grid_set = create_covering_grid(&img, eps, delta);
            plot_covering_grid(&grid_set, "test_grid.png");
        }
        l_val
    }) 
    .await
    .expect("worker panicked")
}

pub async fn calculate_for_all_words_updated(
    length: usize,
    disk_decrease: i32,
    m_seq: Vec<Dyadic>,
) -> f64 {
    let words: Vec<Vec<u8>> = generate_all_words(length);

    let mut best_l = f64::INFINITY;
    let mut best_word: Vec<u8> = Vec::new();

    // Evaluate each word (no plots during the sweep)
    for w in &words {
        let l = calculate_for_word_updated(w.clone(), disk_decrease, m_seq.clone(), false).await;
        if l < best_l {
            best_l = l;
            best_word = w.clone();
        }
    }

    // Re-run the best one with plotting enabled
    let _ = calculate_for_word_updated(best_word.clone(), disk_decrease, m_seq.clone(), true).await;

    println!("Best word (len {}): {:?}\nMin approx = {}", length, best_word, best_l);
    best_l
}