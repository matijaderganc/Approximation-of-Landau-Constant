use std::collections::HashSet;
use std::f64::consts::PI;

use plotters::prelude::*;

use crate::covering_grids::extreme_points;
use crate::dyadic::ComplexDyadic;
use std::path::{Path, PathBuf};


fn plots_path(filename: &str) -> PathBuf {
    let out_dir = Path::new("plots");
    // Create the folder if missing
    let _ = std::fs::create_dir_all(out_dir);

    let p = Path::new(filename);
    if p.components().count() > 1 {
        // Caller already gave a path; use as-is
        p.to_path_buf()
    } else {
        out_dir.join(p)
    }
}

pub fn plot_set(
    points: &Vec<ComplexDyadic>,
    filename: &str,
) -> Result<(), Box<dyn std::error::Error>> {
    // Render to PNG on disk through Plotters (no image crate needed)
    let out_path = plots_path(filename);
    let root = BitMapBackend::new(out_path.to_str().unwrap(), (600, 600)).into_drawing_area();

    // Dark background to make colors pop
    root.fill(&RGBColor(10, 12, 16))?;

    // Compute padded bounding box from your helper
    let extremes = extreme_points(points).ok_or("extreme_points returned None")?;
    let x_min = extremes[0].re.to_f64() - 1.0;
    let x_max = extremes[1].re.to_f64() + 1.0;
    let y_min = extremes[2].im.to_f64() - 1.0;
    let y_max = extremes[3].im.to_f64() + 1.0;

    // Build chart; you can hide mesh if you want a cleaner image
    let mut chart = ChartBuilder::on(&root)
        .margin(10)
        .x_label_area_size(20)
        .y_label_area_size(20)
        .build_cartesian_2d(x_min..x_max, y_min..y_max)?;

    chart.configure_mesh()
        .light_line_style(RGBAColor(255, 255, 255, 16.0)) // very faint grid
        .label_style(("sans-serif", 12).into_font().color(&RGBColor(170, 176, 190)))
        .draw()?;

    // Draw points with domain coloring
    chart.draw_series(points.iter().map(|z| {
        let x = z.re.to_f64();
        let y = z.im.to_f64();
        let c = color_from_complex_rgb(x, y);
        // A tiny 2-pixel radius helps visibility at dense sampling
        Circle::new((x, y), 2, c.filled())
    }))?;

    root.present()?; // ensure file is written
    println!("Set plotted to {}", out_path.display());
    Ok(())
}

pub fn plot_covering_grid(
    grid: &HashSet<ComplexDyadic>,
    filename: &str,
) -> Result<(), Box<dyn std::error::Error>> {
    let out_path = plots_path(filename);
    let root = BitMapBackend::new(out_path.to_str().unwrap(), (600, 600)).into_drawing_area();

    let extremes = extreme_points(&grid.iter().copied().collect()).unwrap();
    let min_real = extremes[0].re.to_f64() - 1.0;
    let max_real = extremes[1].re.to_f64() + 1.0;
    let min_imag = extremes[2].im.to_f64() - 1.0;
    let max_imag = extremes[3].im.to_f64() + 1.0;
    root.fill(&WHITE)?;
    let mut chart = ChartBuilder::on(&root)
        .caption("ε-Covering Grid", ("sans-serif", 30))
        .margin(20)
        .x_label_area_size(30)
        .y_label_area_size(30)
        .build_cartesian_2d(min_real..max_real, min_imag..max_imag)?;

    chart.configure_mesh().draw()?;

    chart.draw_series(grid.iter().map(|z| {
        let (x, y) = (z.re.to_f64(), z.im.to_f64());
        Circle::new((x, y), 1, RED.filled())
    }))?;
    println!("Grid plotted to {}", out_path.display());
    Ok(())
}

fn hsv_to_rgb_u8(h: f64, s: f64, v: f64) -> RGBColor {
    // h in [0,1], s,v in [0,1]
    let h6 = (h.fract() + 1.0).fract() * 6.0;
    let i = h6.floor() as i32;
    let f = h6 - i as f64;
    let p = v * (1.0 - s);
    let q = v * (1.0 - s * f);
    let t = v * (1.0 - s * (1.0 - f));

    let (r, g, b) = match i {
        0 => (v, t, p),
        1 => (q, v, p),
        2 => (p, v, t),
        3 => (p, q, v),
        4 => (t, p, v),
        _ => (v, p, q),
    };
    RGBColor(
        (r.clamp(0.0, 1.0) * 255.0) as u8,
        (g.clamp(0.0, 1.0) * 255.0) as u8,
        (b.clamp(0.0, 1.0) * 255.0) as u8,
    )
}

fn color_from_complex_rgb(x: f64, y: f64) -> RGBColor {
    // Hue from argument
    let arg = y.atan2(x); // [-π, π]
    let hue = (arg + PI) / (2.0 * PI); // [0,1]
    // Value from log|w| stripes
    let r = (x * x + y * y).sqrt();
    let stripes = if r > 0.0 {
        (r.ln() / (2.0_f64).ln()).rem_euclid(1.0)
    } else {
        0.0
    };
    let v = 0.25 + 0.75 * stripes; // prevent over-dark regions
    hsv_to_rgb_u8(hue, 1.0, v)
}




// pub fn plot_set(
//     points: &Vec<ComplexDyadic>,
//     filename: &str,
// ) -> Result<(), Box<dyn std::error::Error>> {
//     let root = BitMapBackend::new(filename, (600, 600)).into_drawing_area();
//     root.fill(&WHITE)?;
//     let extremes = extreme_points(points).unwrap();
//     let x_min = &extremes[0].re.to_f64() - 1.0;
//     let x_max = &extremes[1].re.to_f64() + 1.0;
//     let y_min = &extremes[2].im.to_f64() - 1.0;
//     let y_max = &extremes[3].im.to_f64() + 1.0;
//     let mut chart = ChartBuilder::on(&root)
//         .caption("Complex Dyadic Points", ("sans-serif", 25))
//         .margin(20)
//         .x_label_area_size(30)
//         .y_label_area_size(30)
//         .build_cartesian_2d(x_min..x_max, y_min..y_max)?;

//     chart.configure_mesh().draw()?;

//     chart.draw_series(points.iter().map(|z| {
//         let (x, y) = (z.re.to_f64(), z.im.to_f64());
//         Circle::new((x, y), 2, RED.filled())
//     }))?;

//     println!("Plot saved to '{}'", filename);
//     Ok(())
// }
