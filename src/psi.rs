use crate::dyadic::{Dyadic, ComplexDyadic, Interval, psi} ;
use std::collections::LinkedList;

fn vec_to_linked_list(vec: Vec<u8>) -> LinkedList<u8> {
    let mut list = LinkedList::new();

    for elem in vec {
        list.push_back(elem);
    }

    list
}

fn select_matching(n: u32, vec1: &Vec<u32>, vec2: &Vec<u8>) -> Vec<u8> { //words will be u8
    let mut output = Vec::new();

    for (v1, v2) in vec1.iter().zip(vec2.iter()) {
        if *v1 == n {
            output.push(*v2);
        }
    }
    output
}

pub fn psi_infinity(m_seq: &Vec<Dyadic>, t_seq : &Vec<u32>, word : &Vec<u8>) -> Vec<ComplexDyadic> {
    let mut out = Vec::new() ;
    out.push(ComplexDyadic::new(Dyadic::new(1,0), Dyadic::new(0, 0)));
    for n in 0..m_seq.len() {
        let word_on_step = select_matching(n as u32, t_seq, word) ;
        let interval1 = Interval::new(m_seq[n]*Dyadic::new(-1,0), m_seq[n], m_seq[n]*Dyadic::new(-1,0), m_seq[n]);
        let word_as_list = vec_to_linked_list(word_on_step) ;

        let interval2 = psi(interval1, &word_as_list) ;

        let mid = interval2.midpoint();
        match mid {
            Some(x) => out.push(x), 
            None => panic!("Empty interval")
        }    
    }
    return out
}


