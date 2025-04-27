use landau::dyadic::{Dyadic, ComplexDyadic, Interval, psi} ;
use landau::psi::{psi_infinity} ;
use std::collections::LinkedList;
fn main() {
    let x = Dyadic::new(3, -1) ;
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
    let m_seq = vec![x.clone(), x.clone(), x.clone()] ;
    let t_seq = vec![0, 1, 0, 1, 2] ;
    let word = vec![2, 3, 1, 4, 1] ;
    let holo = psi_infinity(&m_seq, &t_seq, &word) ;
    for a in holo{
        println!("{}, {}", a, a.abs())
    }
}


