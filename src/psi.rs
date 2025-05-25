use crate::dyadic::{psi, ComplexDyadic, Dyadic, Interval};
use rand::Rng;
use std::collections::LinkedList;

fn vec_to_linked_list(vec: Vec<u8>) -> LinkedList<u8> {
    let mut list = LinkedList::new();

    for elem in vec {
        list.push_back(elem);
    }

    list
}

// Generates a random word of a given length
pub fn generate_word(length: i32) -> Vec<u8> {
    let mut rng = rand::thread_rng();
    let mut word = Vec::with_capacity(length as usize);
    for _ in 0..length {
        word.push(rng.gen_range(1..5))
    }
    return word;
}

// Generates all possible vectors of given length with values ranging from 1 to 4 inclusive.
// Used to get all possible words
pub fn generate_all_words(length: usize) -> Vec<Vec<u8>> {
    if length == 0 {
        return vec![vec![]]; // base case: one empty vector
    }

    let smaller_vectors = generate_all_words(length - 1);
    let mut result = Vec::with_capacity(smaller_vectors.len() * 4);

    for vec in smaller_vectors {
        for val in 1..=4 {
            let mut new_vec = vec.clone();
            new_vec.push(val);
            result.push(new_vec);
        }
    }

    result
}

fn select_matching(n: u32, vec1: &Vec<u32>, vec2: &Vec<u8>) -> Vec<u8> {
    //words will be u8
    let mut output = Vec::new();

    for (v1, v2) in vec1.iter().zip(vec2.iter()) {
        if *v1 == n {
            output.push(*v2);
        }
    }
    output
}

pub fn psi_infinity(m_seq: &Vec<Dyadic>, t_seq: &Vec<u32>, word: &Vec<u8>) -> Vec<ComplexDyadic> {
    let mut out = Vec::new();
    out.push(ComplexDyadic::new(Dyadic::new(1, 0), Dyadic::new(0, 0)));
    for n in 0..m_seq.len() {
        let word_on_step = select_matching(n as u32, t_seq, word);
        let interval1 = Interval::new(
            m_seq[n] * Dyadic::new(-1, 0),
            m_seq[n],
            m_seq[n] * Dyadic::new(-1, 0),
            m_seq[n],
        );
        let word_as_list = vec_to_linked_list(word_on_step);

        let interval2 = psi(interval1, &word_as_list);
        let mid = interval2.midpoint();
        match mid {
            Some(x) => out.push(x),
            None => panic!("Empty interval"),
        }
    }

    return out;
}

pub fn mu_first(c: &f64, r: &Dyadic) -> f64 {
    let e = std::f64::consts::E;
    let a = (2.0_f64).sqrt()
        * (c * e
            * (1.0 / 2.0)
            * (1.0 / (1.0 - r.to_f64()) + (1.0 / ((1.0 - r.to_f64()).powf(2.0))))
            + (2.0 - c * e));
    return a;
}

pub fn mu_second(c: &f64, r: &Dyadic) -> f64 {
    let e = std::f64::consts::E;
    let a = (2.0_f64).sqrt()
        * (((c * e) / ((1.0 - r.to_f64()).powf(3.0)))
            + (c * e) / (2.0 * (1.0 - r.to_f64()).powf(2.0))
            + 2.0);
    return a;
}
