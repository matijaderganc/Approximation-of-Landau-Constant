use std::collections::HashSet;

use crate::dyadic::{ComplexDyadic, Dyadic};


//create a unit disk domain as a grid, disk_accuracy means points will be 2^{-n} apart
pub fn unit_disk_n(disk_accuracy: i32) -> Vec<ComplexDyadic> {
    let mut points = Vec::new();
    let max = (2i128).pow(-disk_accuracy as u32);
    for real in -max..=max {
        for imaginary in -max..=max {
            let re = Dyadic::new(real, disk_accuracy);
            let im = Dyadic::new(imaginary, disk_accuracy);
            let z = ComplexDyadic::new(re, im);
            if z.abs() <= 1.0 {
                points.push(z)
            }
        }
    }
    points
}
// unit_disk_radius return a unit disk wiht radius r, and disk_accuracy means points will be 2^{-n} apart
pub fn unit_disk_radius(disk_accuracy: i32, r: f64) -> Vec<ComplexDyadic> {
    let mut points = Vec::new();
    let max = (2i128).pow(-disk_accuracy as u32);
    for real in -max..=max {
        for imaginary in -max..=max {
            let re = Dyadic::new(real, disk_accuracy);
            let im = Dyadic::new(imaginary, disk_accuracy);
            let z = ComplexDyadic::new(re, im);
            if z.abs() <= r {
                points.push(z)
            }
        }
    }
    points
}

//returns just the boundary of unit disk
pub fn unit_disk_n_boundary(disk_accuracy: i32) -> Vec<ComplexDyadic> {
    let mut points = Vec::new();
    let max = (2i128).pow(-disk_accuracy as u32);
    for real in -max..=max {
        for imaginary in -max..=max {
            let re = Dyadic::new(real, disk_accuracy);
            let im = Dyadic::new(imaginary, disk_accuracy);
            let z = ComplexDyadic::new(re, im);
            if z.abs() <= 1.0 + Dyadic::new(1, disk_accuracy).to_f64() && z.abs() >= 1.0 {
                points.push(z)
            }
        }
    }
    points
}

//determines extreme points (lowest/highest real/complex part)
pub fn extreme_points(vec: &Vec<ComplexDyadic>) -> Option<Vec<ComplexDyadic>> {
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
    Some(vec![
        lowest_real,
        highest_real,
        lowest_im,
        highest_im,
    ])
}

//create an epsilon covering grid from an image: a grid point is inside if it is epsilon or closer away from one point in the image
pub fn create_covering_grid(
    set: &Vec<ComplexDyadic>,
    epsilon: Dyadic,
    delta: Dyadic,
) -> HashSet<ComplexDyadic> {
    let mut out: HashSet<ComplexDyadic> = HashSet::new();
    if set.is_empty() {
        return out;
    }

    let delta_f = delta.to_f64();
    let eps_f = epsilon.to_f64();
    let inv_d = 1.0 / delta_f; //precompute

    let r_float = eps_f * inv_d;
    let r = r_float.ceil() as i32; // search square radius
    let r2 = r_float * r_float;
    let mut lattice: HashSet<(i32, i32)> =
        HashSet::with_capacity(set.len() * ((2 * r + 1) as usize));

    for a in set {
        let ax = a.re.to_f64();
        let ay = a.im.to_f64();

        let cx_f = ax * inv_d;
        let cy_f = ay * inv_d;
        let cx = cx_f.round() as i32;
        let cy = cy_f.round() as i32;

        for di in -r..=r {
            let dx = (cx + di) as f64 - cx_f; // lattice delta in x (units of lattice steps)
            let dx2 = dx * dx;
            // quick reject if already beyond radius horizontally
            if dx2 > r2 {
                continue;
            }

            for dj in -r..=r {
                let dy = (cy + dj) as f64 - cy_f;
                if dx2 + dy * dy <= r2 + 1e-14 {
                    // small epsilon to be robust against floating point ties
                    lattice.insert((cx + di, cy + dj));
                }
            }
        }
    }
    out.reserve(lattice.len());
    for (i, j) in lattice {
        let re = Dyadic::new(i as i128, 0) * delta;
        let im = Dyadic::new(j as i128, 0) * delta;
        out.insert(ComplexDyadic::new(re, im));
    }

    out
}

//creates complement of a grid on a scale of its extremes
pub fn grid_complement(grid: &HashSet<ComplexDyadic>, delta: Dyadic) -> HashSet<ComplexDyadic> {
    let extremes = extreme_points(&grid.iter().copied().collect()).unwrap();
    let min_real = extremes[0].re.to_f64();
    let max_real = extremes[1].re.to_f64();
    let min_imag = extremes[2].im.to_f64();
    let max_imag = extremes[3].im.to_f64();

    let delta_f64 = delta.to_f64();

    let r_start = (min_real / delta_f64).floor() as i64 - 1;
    let r_end = (max_real / delta_f64).ceil() as i64 + 1;
    let i_start = (min_imag / delta_f64).floor() as i64 - 1;
    let i_end = (max_imag / delta_f64).ceil() as i64 + 1;

    let mut full_grid = HashSet::new();
    for dr in r_start..=r_end {
        for di in i_start..=i_end {
            let point: ComplexDyadic = ComplexDyadic::new(
                Dyadic::new(dr as i128, 0) * delta,
                Dyadic::new(di as i128, 0) * delta,
            );
            full_grid.insert(point);
        }
    }
    for g in grid {
        full_grid.remove(g);
    }
    full_grid
}

/// approximates the value of a largest disk inside the grid, but this is a slow, quadratic algorithm
/// Not used probably as we have faster algorithms now
pub fn grid_approx(grid: &HashSet<ComplexDyadic>, delta: Dyadic) -> f64 {
    let mut max: Option<f64> = Some(delta.to_f64());
    let comp = grid_complement(grid, delta);
    for point1 in grid {
        let mut min: Option<f64> = None;
        for point2 in &comp {
            let dist = (*point1 - *point2).abs();
            if let Some(current_min) = min {
                if dist < current_min {
                    min = Some(dist)
                }
            } else {
                min = Some(dist)
            }
        }
        if let Some(current_max) = max {
            if min.unwrap() > current_max {
                max = Some(min.unwrap())
            }
        } else {
            max = Some(min.unwrap())
        }
    }
    max.unwrap() - (delta.to_f64() * 4.0)
}

/// We will now write a faster version of grid creation that creates a hash map directly and doesnt bother with
/// complex dyadics, as the type is not important for creating grids and calculating distances. We can compute everything
/// in u8 grid, and then just multiply the result by delta. This saves us a lot of time compared to previous approach
pub struct GridBitmap {
    pub width: usize,
    pub height: usize,
    pub origin_i: i32,
    pub origin_j: i32,
    pub data: Vec<u8>, // 0 = complement, 1 = inside
}

pub fn covering_grid_bitmap(
    points: &[ComplexDyadic],
    epsilon: Dyadic,
    delta: Dyadic,
) -> GridBitmap {
    assert!(
        !points.is_empty(),
        "covering_grid_bitmap: empty image set"
    );

    let delta_f = delta.to_f64();
    let eps_f = epsilon.to_f64();
    let inv_d = 1.0 / delta_f;

    let r_f = eps_f * inv_d; // ε / δ
    let r = r_f.ceil() as i32; // integer square radius
    let r2 = r_f * r_f;

    #[derive(Copy, Clone)]
    struct Center {
        cx: i32,
        cy: i32,
        fx: f64,
        fy: f64,
    }
    let mut centers: Vec<Center> = Vec::with_capacity(points.len());
    let (mut min_i, mut max_i, mut min_j, mut max_j) = (i32::MAX, i32::MIN, i32::MAX, i32::MIN);

    for z in points {
        let ax = z.re.to_f64();
        let ay = z.im.to_f64();

        let cx_f = ax * inv_d;
        let cy_f = ay * inv_d;

        let cx = cx_f.round() as i32;
        let cy = cy_f.round() as i32;

        // sub-pixel offsets (in lattice units), in [-0.5, 0.5]
        let fx = cx as f64 - cx_f;
        let fy = cy as f64 - cy_f;

        centers.push(Center { cx, cy, fx, fy });

        // bbox grows by ±r around each center
        if cx - r < min_i {
            min_i = cx - r;
        }
        if cx + r > max_i {
            max_i = cx + r;
        }
        if cy - r < min_j {
            min_j = cy - r;
        }
        if cy + r > max_j {
            max_j = cy + r;
        }
    }

    min_i -= 1;
    max_i += 1;
    min_j -= 1;
    max_j += 1;

    let width = (max_i - min_i + 1) as usize;
    let height = (max_j - min_j + 1) as usize;

    // Dense bitmap (0 = complement, 1 = inside)
    let mut data = vec![0u8; width * height];

    for c in &centers {
        // sweep around the center within the square [-r, r]^2
        for di in -r..=r {
            let dx = di as f64 + c.fx;
            let dx2 = dx * dx;
            if dx2 > r2 {
                continue;
            } // early skip if horizontal distance already too large
            let x = (c.cx + di - min_i) as usize;

            for dj in -r..=r {
                let dy = dj as f64 + c.fy;
                if dx2 + dy * dy <= r2 + 1e-12 {
                    let y = (c.cy + dj - min_j) as usize;
                    data[y * width + x] = 1;
                }
            }
        }
    }

    GridBitmap {
        width,
        height,
        origin_i: min_i,
        origin_j: min_j,
        data,
    }
}

/// we use this function just so we can plot our covering grids as HashSets of ComplexDyadics
pub fn bitmap_to_complex_set(g: &GridBitmap, delta: Dyadic) -> HashSet<ComplexDyadic> {
    let mut out = HashSet::new();
    if g.data.is_empty() {
        return out;
    }

    // Precompute dyadic coordinates for each column (x) and row (y).
    // re_vals[x] = (origin_i + x) * delta, im_vals[y] = (origin_j + y) * delta
    let mut re_vals: Vec<Dyadic> = Vec::with_capacity(g.width);
    for x in 0..g.width {
        let i = g.origin_i as i128 + x as i128;
        re_vals.push(Dyadic::new(i, 0) * delta);
    }
    let mut im_vals: Vec<Dyadic> = Vec::with_capacity(g.height);
    for y in 0..g.height {
        let j = g.origin_j as i128 + y as i128;
        im_vals.push(Dyadic::new(j, 0) * delta);
    }

    out.reserve(g.width * g.height / 8 + 1);
    #[allow(clippy::needless_range_loop)]
    for y in 0..g.height {
        for x in 0..g.width {
            if g.data[y * g.width + x] != 0 {
                let re = re_vals[x];
                let im = im_vals[y];
                out.insert(ComplexDyadic::new(re, im));
            }
        }
    }

    out
}
