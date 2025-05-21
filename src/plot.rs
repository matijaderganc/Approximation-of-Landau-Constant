use crate::dyadic::{Dyadic, ComplexDyadic} ;
use crate::covering_grids::{extreme_points};
use plotters::prelude::*;
use std::collections::HashSet;



pub fn plot_set(points: &Vec<ComplexDyadic>, filename : &str) -> Result<(), Box<dyn std::error::Error>>{
    let root = BitMapBackend::new(filename, (600, 600)).into_drawing_area();
    root.fill(&WHITE)?;
    let extremes = extreme_points(points).unwrap() ;
    let x_min = &extremes[0].re.to_f64() - 1.0 ; 
    let x_max = &extremes[1].re.to_f64() + 1.0 ; 
    let y_min = &extremes[2].im.to_f64() - 1.0 ; 
    let y_max = &extremes[3].im.to_f64() + 1.0 ; 
    let mut chart = ChartBuilder::on(&root)
        .caption("Complex Dyadic Points", ("sans-serif", 25))
        .margin(20)
        .x_label_area_size(30)
        .y_label_area_size(30)
        .build_cartesian_2d(x_min..x_max, y_min..y_max)?;

    chart.configure_mesh().draw()?;

    chart.draw_series(
        points.iter().map(|z| {
            let (x, y) = (z.re.to_f64(), z.im.to_f64());
            Circle::new((x, y), 2, RED.filled())
        })
    )?;

    println!("Plot saved to '{}'", filename);
    Ok(())
}

pub fn plot_covering_grid(grid: &HashSet<ComplexDyadic>, filename: &str) -> Result<(), Box<dyn std::error::Error>> {
    let extremes = extreme_points(&grid.iter().copied().collect()).unwrap();
    let min_real = extremes[0].re.to_f64() - 1.0 ;
    let max_real = extremes[1].re.to_f64() + 1.0 ;
    let min_imag = extremes[2].im.to_f64() - 1.0 ;
    let max_imag = extremes[3].im.to_f64() + 1.0;
    let root = BitMapBackend::new(filename, (600, 600)).into_drawing_area();
    root.fill(&WHITE)?;

    let range = -5.1f64..5.1f64;

    let mut chart = ChartBuilder::on(&root)
        .caption("ε-Covering Grid", ("sans-serif", 30))
        .margin(20)
        .x_label_area_size(30)
        .y_label_area_size(30)
        .build_cartesian_2d(min_real..max_real, min_imag..max_imag)?;

    chart.configure_mesh().draw()?;

    chart.draw_series(
        grid.iter().map(|z| {
            let (x, y) = (z.re.to_f64(), z.im.to_f64());
            Circle::new((x, y), 1, RED.filled())
        })
    )?;
    println!("{}", grid.len()) ;
    println!("Grid plotted to {}", filename);
    Ok(())
}