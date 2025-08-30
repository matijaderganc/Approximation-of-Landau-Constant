use std::cmp::Ordering;
use std::sync::Arc;

use futures::stream::{self, StreamExt};
use num_cpus;

use crate::covering_grids::{unit_disk_n, covering_grid_bitmap};
use crate::dyadic::{Dyadic, ComplexDyadic};
use crate::edt::landau_l_via_edt_from_bitmap;
use crate::holomorphic::{BoundingSequence, ComplexFunction, ExpansionCoefficients};
use crate::psi::{psi_infinity, t_vector, generate_all_words};

pub async fn calculate_for_length(
    length: usize,
    delta: Dyadic,
    disk_radius: i32,
    m_seq : Vec<Dyadic>,
) -> (Vec<u8>, f64) {
    let epsilon = delta * Dyadic::new(1, 2) ;
    
    let m_seq1  = Arc::new(m_seq);
    let t_seq = Arc::new(t_vector(100));
    let domain = Arc::new(unit_disk_n(disk_radius));
    let concurrency = num_cpus::get().max(1);
    // println!("Working on {} CPU cores", concurrency) ;
    let words: Vec<Vec<u8>> = generate_all_words(length);

    let mut results = stream::iter(words.into_iter())
        .map(|word| {
            let domain = Arc::clone(&domain);
            let m_seq  = Arc::clone(&m_seq1);
            let t_seq  = Arc::clone(&t_seq);

            tokio::task::spawn_blocking(move || {
                // Build f' from the word, then integrate to f
                let coeffs = psi_infinity(&m_seq, &t_seq, &word);
                let fprime = ComplexFunction::new(
                    BoundingSequence::new((*m_seq).clone()),
                    ExpansionCoefficients::new(coeffs),
                );
                let f = fprime.antiderivative();

                // Evaluate on domain
                let mut img: Vec<ComplexDyadic> = Vec::with_capacity(domain.len());
                for &z in domain.iter() {
                    img.push(f.eval(z));
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

    // Print nicely
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