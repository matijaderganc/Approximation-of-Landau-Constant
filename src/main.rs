use landau::covering_grids::{
    self, create_covering_grid, extreme_points, plot_covering_grid, unit_disk_n,
};
use landau::dyadic::{psi, ComplexDyadic, Dyadic, Interval};
use landau::holomorphic::{
    comp_vec_to_sequence, vec_to_sequence, BoundingSequence, ComplexFunction, ExpansionCoefficients,
};
use landau::psi::{mu_first, mu_second, psi_infinity};

use plotters::prelude::*;
use std::collections::LinkedList;

fn main() // Result<(), Box<dyn std::error::Error>>
{
    // let x = Dyadic::new(3, -1);
    // let m_seq = vec![x.clone(), x.clone(), x.clone()];
    // let x = Dyadic::new(3, -1) ;
    // let m_seq = vec![x.clone(), x.clone(), x.clone()] ;
    // let t_seq = vec![0, 1, 0, 1, 2, 1, 4, 1, 2, 4, 4, 3, 1, 4, 5];
    // let word = vec![1, 3, 2, 1, 2, 4, 3, 1, 2, 3, 1, 2, 3, 1, 3, 2];
    // let holo = psi_infinity(&m_seq, &t_seq, &word);
    // let points1 = unit_disk_n(-1);

    // let f = ComplexFunction::new(
    //     BoundingSequence::new(m_seq),
    //     ExpansionCoefficients::new(holo),
    // );
    // println!("Complex function compiled");
    // TODO use holomprphic.rs to represent the f above, then run the following code
    // let mut points2 = Vec::new() ;
    // println!("{}", f.eval(ComplexDyadic::new(x.clone(), x.clone())));
    let x = Dyadic::new(1, 0);
    let m_seq = vec![x.clone(), x.clone(), x.clone()];
    // This is a polynomial with factors 1, 1, 1
    let holo = vec![
        ComplexDyadic::new(Dyadic::new(1, 0), Dyadic::new(0, 0)),
        ComplexDyadic::new(Dyadic::new(1, 0), Dyadic::new(0, 0)),
        ComplexDyadic::new(Dyadic::new(1, 0), Dyadic::new(0, 0)),
    ];

    let f = ComplexFunction::new(
        BoundingSequence::new(m_seq),
        ExpansionCoefficients::new(holo),
    );
    println!("Complex function compiled");

    println!(
        "Function result: {}, {:?}",
        f.eval(ComplexDyadic::new(Dyadic::new(1, 0), Dyadic::new(1, 0))),
        f.expansion_coefficients
    );
    println!(
        "Derivative result: {}, {:?}",
        f.derivative()
            .eval(ComplexDyadic::new(Dyadic::new(1, 0), Dyadic::new(1, 0))),
        f.derivative().expansion_coefficients
    )
    // let root = BitMapBackend::new("dyadic_disk.png", (600, 600)).into_drawing_area();
    // root.fill(&WHITE)?;

    let mut chart = ChartBuilder::on(&root)
        .caption("Dyadic Points in Unit Disk", ("sans-serif", 30))
        .margin(20)
        .x_label_area_size(30)
        .y_label_area_size(30)
        .build_cartesian_2d(-1.1..1.1, -1.1..1.1)?;

    chart.configure_mesh().draw()?;

    // chart.draw_series(
    //     points1.iter().map(|p| {
    //         let (x, y) = (p.re.to_f64(), p.im.to_f64());
    //         Circle::new((x, y), 3, RED.filled())
    //     })
    // )?;
    // println!("Plot saved as 'dyadic_disk.png'");

    let root2 = BitMapBackend::new("dyadic_image.png", (600, 600)).into_drawing_area();
    root2.fill(&WHITE)?;

    let mut chart2 = ChartBuilder::on(&root2)
        .caption("Dyadic Image of Unit Disk", ("sans-serif", 30))
        .margin(20)
        .x_label_area_size(30)
        .y_label_area_size(30)
        .build_cartesian_2d(-5.1..5.1, -5.1..5.1)?;

    chart2.configure_mesh().draw()?;

    chart2.draw_series(
        points2.iter().map(|p| {
            let (x, y) = (p.re.to_f64(), p.im.to_f64());
            Circle::new((x, y), 3, RED.filled())
        })
    )?;
    println!("Plot saved as 'dyadic_image.png'");
    let grid = create_covering_grid(&points2, Dyadic::new(1, -2), Dyadic::new(1, -5)) ;
    plot_covering_grid(&grid, "covering_grid.png")?;
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
