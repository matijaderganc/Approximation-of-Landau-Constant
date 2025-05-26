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

pub fn grid_complement(grid : &HashSet<ComplexDyadic>, delta : Dyadic) -> HashSet<ComplexDyadic> {
    let extremes = extreme_points(&grid.iter().copied().collect()).unwrap();
    let min_real = extremes[0].re.to_f64();
    let max_real = extremes[1].re.to_f64();
    let min_imag = extremes[2].im.to_f64();
    let max_imag = extremes[3].im.to_f64();

    let delta_f64 = delta.to_f64();

    let r_start = (min_real / delta_f64).floor() as i64 - 1;
    let r_end   = (max_real / delta_f64).ceil() as i64 + 1;
    let i_start = (min_imag / delta_f64).floor() as i64 - 1;
    let i_end   = (max_imag / delta_f64).ceil() as i64 + 1;

    let mut full_grid = HashSet::new();
    for dr in r_start..=r_end {
        for di in i_start..=i_end {
            let point: ComplexDyadic = ComplexDyadic::new(Dyadic::new(dr as i128, 0) * delta, Dyadic::new(di as i128, 0) * delta) ;
            full_grid.insert(point);
        }
    } ;
    for g in grid {
        full_grid.remove(g);
    }
    full_grid
}

pub fn grid_approx(grid : &HashSet<ComplexDyadic>, delta : Dyadic) -> f64 {
    let mut max : Option<f64> = Some(delta.to_f64()) ;
    let comp = grid_complement(grid, delta) ;
    for point1 in grid{
        let mut min : Option<f64> = None ;
        for point2 in &comp {
            let dist = (*point1 - *point2).abs() ;
            if let Some(current_min) = min {
                if dist < current_min {
                    min = Some(dist)
                }
            }
            else {
                min = Some(dist)
            }
        }
        if let Some(current_max) = max {
            if min.unwrap() > current_max {
                max = Some(min.unwrap())
            }
        }
        else {
            max = Some(min.unwrap())
        }
    }
    max.unwrap() - (delta.to_f64()*4.0)
}


