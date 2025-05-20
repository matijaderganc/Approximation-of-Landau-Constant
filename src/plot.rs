use crate::dyadic::{Dyadic, ComplexDyadic} ;
use crate::covering_grids::{extreme_points};
use plotters::prelude::*;


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
