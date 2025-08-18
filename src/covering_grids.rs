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
pub fn unit_disk_n_boundary(n : i32) -> Vec<ComplexDyadic> {
    let mut points = Vec::new()  ;
    let max = (2i128).pow(-n as u32) ;
    for real in -max..=max {
        for imaginary in -max..=max {
            let re = Dyadic::new(real, n) ;
            let im = Dyadic::new(imaginary, n) ;
            let z = ComplexDyadic::new(re, im) ;
            if z.abs() <= 1.0 + Dyadic::new(1,n).to_f64() && z.abs() >= 1.0 {
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
    let mut out: HashSet<ComplexDyadic> = HashSet::new();
    let delta_f64 = delta.to_f64();
    if set.is_empty() { return out; }

    let delta_f = delta.to_f64();
    let eps_f   = epsilon.to_f64();
    let inv_d   = 1.0 / delta_f;
    
    let r_float = eps_f * inv_d;           // ε / δ
    let r       = r_float.ceil() as i32;   // search square radius
    let r2      = r_float * r_float; 
    let mut lattice: HashSet<(i32, i32)> = HashSet::with_capacity(set.len() * ((2*r + 1) as usize));

    for a in set {
        let ax = a.re.to_f64();
        let ay = a.im.to_f64();
        
        let cx_f = ax * inv_d;
        let cy_f = ay * inv_d;
        let cx   = cx_f.round() as i32;
        let cy   = cy_f.round() as i32;

        for di in -r..=r {
            let dx = (cx + di) as f64 - cx_f;    // lattice delta in x (units of lattice steps)
            let dx2 = dx * dx;
            // quick reject if already beyond radius horizontally
            if dx2 > r2 { continue; }

            for dj in -r..=r {
                let dy = (cy + dj) as f64 - cy_f;
                if dx2 + dy*dy <= r2 + 1e-12 {   // small epsilon to be robust against FP ties
                    lattice.insert((cx + di, cy + dj));
                }
            }}
    } 
    out.reserve(lattice.len());
    for (i, j) in lattice {
        // i*δ, j*δ  (use your dyadic ops; this happens once per kept cell)
        let re = Dyadic::new(i as i128, 0) * delta;
        let im = Dyadic::new(j as i128, 0) * delta;
        out.insert(ComplexDyadic::new(re, im));
    }

    out
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

pub fn grid_approx_with_edge(inside : &HashSet<ComplexDyadic>, edge: &HashSet<ComplexDyadic>, delta : Dyadic) -> f64 {
    let mut max : Option<f64> = Some(delta.to_f64()) ;
    for point1 in inside{
        let mut min : Option<f64> = None ;
        for point2 in edge {
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
    max.unwrap() 
}

