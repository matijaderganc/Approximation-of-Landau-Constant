use std::collections::HashSet;
use std::f64;

use crate::covering_grids::GridBitmap;
use crate::dyadic::{ComplexDyadic, Dyadic};

// Felzenszwalb & Huttenlocher (2004), "Distance Transforms of Sampled Functions"
// Given f[q] where f[p]=0 for "sources" and f[p]=∞ elsewhere, returns
// g[q] = min_p ( (q-p)^2 + f[p] ).
const INF: f64 = 1e20;

pub fn edt_1d_squared(f: &[f64]) -> Vec<f64> {
    let n = f.len();
    let mut v = vec![0usize; n]; // positions of parabolas in envelope
    let mut z = vec![f64::NEG_INFINITY; n + 1]; // breakpoints between parabolas
    let mut g = vec![0f64; n]; // output

    let mut k: usize = 0; // number of parabolas in envelope - 1
    v[0] = 0;
    z[0] = f64::NEG_INFINITY;
    z[1] = f64::INFINITY;

    // Build lower envelope
    for q in 1..n {
        // intersection with current parabola v[k]
        let mut s = ((f[q] + (q * q) as f64) - (f[v[k]] + (v[k] * v[k]) as f64))
            / (2.0 * (q as f64 - v[k] as f64));
        // pop while new parabola takes over earlier than current breakpoint
        while k > 0 && s <= z[k] {
            k -= 1;
            s = ((f[q] + (q * q) as f64) - (f[v[k]] + (v[k] * v[k]) as f64))
                / (2.0 * (q as f64 - v[k] as f64));
        }
        k += 1;
        v[k] = q;
        z[k] = s;
        z[k + 1] = f64::INFINITY;
    }

    // Evaluate envelope
    k = 0;
    #[allow(clippy::needless_range_loop)]
    //clippy is an annoying bastard and doesnt allow needless range loops :(
    for q in 0..n {
        while z[k + 1] < q as f64 {
            k += 1;
        }
        let p = v[k];
        let dq = q as f64 - p as f64;
        g[q] = dq * dq + f[p];
    }

    g
}


pub fn edt_2d_squared(img: &mut [f64], width: usize, height: usize) {
    // first pass: columns
    {
        let mut col = vec![0.0; height];
        let mut out = vec![0.0; height];
        for x in 0..width {
            for y in 0..height {
                col[y] = img[y * width + x];
            }
            let tmp = edt_1d_squared(&col);
            out.copy_from_slice(&tmp);
            for y in 0..height {
                img[y * width + x] = out[y];
            }
        }
    }
    // second pass: rows
    {
        let mut row = vec![0.0; width];
        let mut out = vec![0.0; width];
        for y in 0..height {
            for x in 0..width {
                row[x] = img[y * width + x];
            }
            let tmp = edt_1d_squared(&row);
            out.copy_from_slice(&tmp);
            for x in 0..width {
                img[y * width + x] = out[x];
            }
        }
    }
}

pub fn print_grid(label: &str, img: &[f64], w: usize, h: usize, sqrt: bool) {
    println!("{label}:");
    for y in 0..h {
        for x in 0..w {
            let v = img[y * w + x];
            let v = if sqrt { v.sqrt() } else { v };
            print!("{:4.0}", v);
        }
        println!();
    }
    println!();
}

/// Compute l(ε,δ,G) in linear time using your EDT, from a grid of ComplexDyadic.
/// - `grid`: your image set G as lattice points (ComplexDyadic snapped by δ when plotting).
/// - `delta`: the lattice step δ (Dyadic).
///
/// Steps:
/// 1) snap ComplexDyadic -> (i,j) integer lattice via δ
/// 2) build a dense field over a 1-cell padded bbox: 0.0 on complement, INF on inside
/// 3) run edt_2d_squared in-place (squared distances to nearest complement)
/// 4) read max distance over INSIDE sites, convert to metric: δ + δ * sqrt(max_sq)
pub fn landau_l_with_edt_from_complex_set(grid: &HashSet<ComplexDyadic>, delta: Dyadic) -> f64 {
    assert!(!grid.is_empty(), "grid is empty");

    // --- 1) snap to δ-lattice (i,j) ---
    let step = delta.to_f64();
    let inv = 1.0 / step;
    let mut lattice: HashSet<(i32, i32)> = HashSet::with_capacity(grid.len());
    for z in grid {
        let i = (z.re.to_f64() * inv).round() as i32;
        let j = (z.im.to_f64() * inv).round() as i32;
        lattice.insert((i, j));
    }

    // --- 2) padded bbox so complement exists on border ---
    let (mut min_i, mut max_i, mut min_j, mut max_j) = (i32::MAX, i32::MIN, i32::MAX, i32::MIN);
    for &(i, j) in &lattice {
        if i < min_i {
            min_i = i;
        }
        if i > max_i {
            max_i = i;
        }
        if j < min_j {
            min_j = j;
        }
        if j > max_j {
            max_j = j;
        }
    }
    min_i -= 1;
    max_i += 1;
    min_j -= 1;
    max_j += 1;

    let width = (max_i - min_i + 1) as usize;
    let height = (max_j - min_j + 1) as usize;

    // --- 3) dense field: 0.0 at complement, INF at inside ---
    let mut field = vec![0.0f64; width * height]; // row-major
    let big = INF; // if not pub: let big = 1e20;
    for &(i, j) in &lattice {
        let x = (i - min_i) as usize;
        let y = (j - min_j) as usize;
        field[y * width + x] = big;
    }

    // --- 4) exact squared Euclidean distances (in-place) ---
    edt_2d_squared(&mut field, width, height);

    // --- 5) max distance among INSIDE sites only ---
    let mut max_sq = 0.0f64;
    for &(i, j) in &lattice {
        let x = (i - min_i) as usize;
        let y = (j - min_j) as usize;
        let d2 = field[y * width + x];
        if d2.is_finite() && d2 > max_sq {
            max_sq = d2;
        }
    }

    // --- 6) convert pixels -> metric and add the leading +δ ---
    step + step * max_sq.sqrt()
}

pub fn landau_l_via_edt_from_bitmap(g: &GridBitmap, delta: Dyadic) -> f64 {
    // Build the EDT field: 0 at complement, INF at inside
    let mut field = vec![0.0f64; g.width * g.height];
    for (idx, &pix) in g.data.iter().enumerate() {
        if pix != 0 {
            field[idx] = INF;
        }
    }

    // Exact squared Euclidean distances, in-place
    edt_2d_squared(&mut field, g.width, g.height);

    // Max over inside cells only
    let mut max_sq = 0.0f64;
    for (idx, &pix) in g.data.iter().enumerate() {
        if pix != 0 {
            let d2 = field[idx];
            if d2.is_finite() && d2 > max_sq {
                max_sq = d2;
            }
        }
    }

    let d = delta.to_f64();
    d + d * max_sq.sqrt()
}
