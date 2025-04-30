use crate::dyadic::{Dyadic, ComplexDyadic, Interval, psi} ;
use std::collections::HashSet;
use plotters::prelude::*;

pub struct CoveringGrid {
    pub grid : Vec<(i32, i32)> ,
    pub epsilon : f64 ,
    pub delta : Dyadic , 
}

pub fn unit_disk_n( n : i32) -> Vec<ComplexDyadic> { //could be improved, we check the whole square
    let mut points = Vec::new()  ;
    let max = (2i128).pow(-n as u32) ;
    for real in -max..=max {
        for imaginary in -max..=max {
            let re = Dyadic::new(real, n) ;
            let im = Dyadic::new(imaginary, n) ;
            let z = ComplexDyadic::new(re, im) ;
            if z.abs() <= 1.0 {
                points.push(z)
            }
        }
    }
    return points
}

pub fn extreme_points(vec : &Vec<ComplexDyadic>) -> Option<Vec<ComplexDyadic>> {
    if vec.is_empty() {
        return None; // No points, return None
    }
    let mut lowest_real = vec[0];
    let mut highest_real = vec[0];
    let mut lowest_im = vec[0];
    let mut highest_im = vec[0];

    for &point in vec {
        if point.re.to_f64() < lowest_real.re.to_f64() {
            lowest_real = point
        }
        if point.re.to_f64() > highest_real.re.to_f64() {
            highest_real = point
        }
        if point.im.to_f64() < lowest_im.im.to_f64() {
            lowest_im = point
        }
        if point.im.to_f64() > highest_im.im.to_f64() {
            highest_im = point
        }
    }
    return Some(vec![lowest_real, highest_real, lowest_im, highest_im])
}

pub fn create_covering_grid(set : &Vec<ComplexDyadic>, epsilon : Dyadic, delta : Dyadic) -> HashSet<ComplexDyadic> {
    let mut grid = HashSet::new();
    let delta_f64 = delta.to_f64();
    let epsilon_f64 = epsilon.to_f64();
    if set.is_empty() {
        return grid;
    }

    let grid_radius = (epsilon.to_f64() / delta.to_f64()).ceil() as i128;
    for a in set {
        let real_a = a.re.to_f64();
        let imag_a = a.im.to_f64();
        let center_real = (real_a / delta.to_f64()).round() as i128;
        let center_imag = (imag_a / delta.to_f64()).round() as i128;
        let center = ComplexDyadic::new(Dyadic::new(center_real, 0) * delta, Dyadic::new(center_imag, 0) * delta) ;
        for dr in -grid_radius..=grid_radius {
            for di in -grid_radius..=grid_radius {
                let diff = ComplexDyadic::new(Dyadic::new(dr, 0) * delta, Dyadic::new(di, 0) * delta) ;

                let candidate = center + diff ;

                let sub = *a - candidate.clone() ;
                let dist2 = sub.abs() ;
                if dist2 <= epsilon_f64 {
                    grid.insert(candidate);
                }
            }
        }
    }
    grid
}

pub fn plot_covering_grid(grid: &HashSet<ComplexDyadic>, filename: &str) -> Result<(), Box<dyn std::error::Error>> {
    let root = BitMapBackend::new(filename, (600, 600)).into_drawing_area();
    root.fill(&WHITE)?;

    let range = -5.1f64..5.1f64;

    let mut chart = ChartBuilder::on(&root)
        .caption("ε-Covering Grid", ("sans-serif", 30))
        .margin(20)
        .x_label_area_size(30)
        .y_label_area_size(30)
        .build_cartesian_2d(range.clone(), range.clone())?;

    chart.configure_mesh().draw()?;

    chart.draw_series(
        grid.iter().map(|z| {
            let (x, y) = (z.re.to_f64(), z.im.to_f64());
            Circle::new((x, y), 1, RED.filled())
        })
    )?;

    println!("Grid plotted to {}", filename);
    Ok(())
}