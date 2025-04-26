use landau::dyadic::{Dyadic, ComplexDyadic, Interval, psi} ;
use std::collections::LinkedList;

fn main() {
    let x = Dyadic::new(3, -2) ;
    let y = Dyadic::new(5, -2);
    let z = Dyadic::new(4, -1);
    let w = Dyadic::new(7, -1);
    // let i1 = Interval::new(x, y, z, w);
    // let mut word = LinkedList::new();
    // for _ in 0..5 {
    //     word.push_back(1);
    // }
    // let i3 = psi(i1, &word) ;
    // print!("{}", i3) ;
    let alpha = ComplexDyadic::new(x, y);
    let beta = ComplexDyadic::new(z, w) ;
    println!("{}", alpha * beta)
}


