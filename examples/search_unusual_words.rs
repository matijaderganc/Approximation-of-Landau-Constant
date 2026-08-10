use landau::{
    covering_grids::unit_disk_n,
    dyadic::Dyadic,
    holomorphic::{BoundingSequence, ComplexFunction, ExpansionCoefficients},
    psi::{psi_infinity, t_vector},
};
use plotters::prelude::*;

const WORD_LENGTH: usize = 8;
const BOUNDARY_SAMPLES: usize = 192;
const GALLERY_SIZE: usize = 12;

#[derive(Clone, Debug)]
struct Candidate {
    word: Vec<u8>,
    score: f64,
    concavity: f64,
    radial_variation: f64,
    cusp_strength: f64,
}

fn m_sequence() -> Vec<Dyadic> {
    vec![Dyadic::new(1, -1), Dyadic::new(1, -1), Dyadic::new(1, -1)]
}

fn function_for_word(word: &[u8]) -> ComplexFunction {
    let m_seq = m_sequence();
    let t_seq = t_vector(word.len());
    let derivative_coefficients = psi_infinity(&m_seq, &t_seq, word);
    let f_prime = ComplexFunction::new(
        BoundingSequence::new(m_seq),
        ExpansionCoefficients::new(derivative_coefficients),
    );
    f_prime.antiderivative()
}

fn f64_coefficients(word: &[u8]) -> Vec<(f64, f64)> {
    function_for_word(word)
        .expansion_coefficients
        .vector
        .iter()
        .map(|z| (z.re.to_f64(), z.im.to_f64()))
        .collect()
}

fn evaluate(coefficients: &[(f64, f64)], z: (f64, f64)) -> (f64, f64) {
    coefficients.iter().rev().fold((0.0, 0.0), |acc, &c| {
        (
            acc.0 * z.0 - acc.1 * z.1 + c.0,
            acc.0 * z.1 + acc.1 * z.0 + c.1,
        )
    })
}

fn boundary(coefficients: &[(f64, f64)]) -> Vec<(f64, f64)> {
    (0..BOUNDARY_SAMPLES)
        .map(|i| {
            let theta = std::f64::consts::TAU * i as f64 / BOUNDARY_SAMPLES as f64;
            evaluate(coefficients, (theta.cos(), theta.sin()))
        })
        .collect()
}

fn cross(o: (f64, f64), a: (f64, f64), b: (f64, f64)) -> f64 {
    (a.0 - o.0) * (b.1 - o.1) - (a.1 - o.1) * (b.0 - o.0)
}

fn polygon_area(points: &[(f64, f64)]) -> f64 {
    if points.len() < 3 {
        return 0.0;
    }
    points
        .iter()
        .zip(points.iter().cycle().skip(1))
        .take(points.len())
        .map(|(a, b)| a.0 * b.1 - a.1 * b.0)
        .sum::<f64>()
        .abs()
        * 0.5
}

fn convex_hull_area(points: &[(f64, f64)]) -> f64 {
    let mut sorted = points.to_vec();
    sorted.sort_by(|a, b| a.0.total_cmp(&b.0).then(a.1.total_cmp(&b.1)));
    sorted.dedup();
    if sorted.len() < 3 {
        return 0.0;
    }

    let mut lower = Vec::with_capacity(sorted.len());
    for &point in &sorted {
        while lower.len() >= 2
            && cross(
                lower[lower.len() - 2],
                lower[lower.len() - 1],
                point,
            ) <= 0.0
        {
            lower.pop();
        }
        lower.push(point);
    }

    let mut upper = Vec::with_capacity(sorted.len());
    for &point in sorted.iter().rev() {
        while upper.len() >= 2
            && cross(
                upper[upper.len() - 2],
                upper[upper.len() - 1],
                point,
            ) <= 0.0
        {
            upper.pop();
        }
        upper.push(point);
    }

    lower.pop();
    upper.pop();
    lower.extend(upper);
    polygon_area(&lower)
}

fn score_word(word: Vec<u8>) -> Candidate {
    let points = boundary(&f64_coefficients(&word));
    let center = points
        .iter()
        .fold((0.0, 0.0), |sum, p| (sum.0 + p.0, sum.1 + p.1));
    let center = (
        center.0 / points.len() as f64,
        center.1 / points.len() as f64,
    );

    let radii: Vec<f64> = points
        .iter()
        .map(|p| ((p.0 - center.0).powi(2) + (p.1 - center.1).powi(2)).sqrt())
        .collect();
    let mean_radius = radii.iter().sum::<f64>() / radii.len() as f64;
    let radial_variation =
        (radii.iter().map(|r| (r - mean_radius).powi(2)).sum::<f64>() / radii.len() as f64).sqrt()
            / mean_radius;

    let steps: Vec<f64> = points
        .iter()
        .zip(points.iter().cycle().skip(1))
        .take(points.len())
        .map(|(a, b)| ((a.0 - b.0).powi(2) + (a.1 - b.1).powi(2)).sqrt())
        .collect();
    let mean_step = steps.iter().sum::<f64>() / steps.len() as f64;
    let min_step = steps.iter().copied().fold(f64::INFINITY, f64::min);
    let cusp_strength = (1.0 - min_step / mean_step).clamp(0.0, 1.0);

    let area = polygon_area(&points);
    let hull_area = convex_hull_area(&points);
    let concavity = if hull_area > 0.0 {
        (1.0 - area / hull_area).clamp(0.0, 1.0)
    } else {
        0.0
    };

    // Concavity identifies inward dents, radial variation rejects near-circles,
    // and cusp strength favors cardioid/heart-like pointed boundaries.
    let score = 6.0 * concavity + 2.0 * radial_variation + cusp_strength;
    Candidate {
        word,
        score,
        concavity,
        radial_variation,
        cusp_strength,
    }
}

fn word_from_index(mut index: usize) -> Vec<u8> {
    let mut word = vec![1; WORD_LENGTH];
    for digit in &mut word {
        *digit = (index % 4 + 1) as u8;
        index /= 4;
    }
    word
}

fn render_gallery(candidates: &[Candidate]) -> Result<(), Box<dyn std::error::Error>> {
    std::fs::create_dir_all("plots")?;
    let root =
        BitMapBackend::new("plots/unusual_words_gallery.png", (1600, 1200)).into_drawing_area();
    root.fill(&WHITE)?;
    let panels = root.split_evenly((3, 4));
    let domain = unit_disk_n(-6);

    for (panel, candidate) in panels.into_iter().zip(candidates) {
        let f = function_for_word(&candidate.word);
        let image: Vec<_> = domain.iter().map(|z| f.eval(z)).collect();
        let points: Vec<_> = image
            .iter()
            .map(|z| (z.re.to_f64(), z.im.to_f64()))
            .collect();

        let min_x = points.iter().map(|p| p.0).fold(f64::INFINITY, f64::min);
        let max_x = points.iter().map(|p| p.0).fold(f64::NEG_INFINITY, f64::max);
        let min_y = points.iter().map(|p| p.1).fold(f64::INFINITY, f64::min);
        let max_y = points.iter().map(|p| p.1).fold(f64::NEG_INFINITY, f64::max);
        let center = ((min_x + max_x) * 0.5, (min_y + max_y) * 0.5);
        let half_span = ((max_x - min_x).max(max_y - min_y) * 0.58).max(0.1);
        let title = format!(
            "{:?}  score {:.3}",
            candidate.word, candidate.score
        );

        let mut chart = ChartBuilder::on(&panel)
            .caption(title, ("sans-serif", 18))
            .margin(12)
            .build_cartesian_2d(
                center.0 - half_span..center.0 + half_span,
                center.1 - half_span..center.1 + half_span,
            )?;
        chart.configure_mesh().disable_mesh().draw()?;
        chart.draw_series(points.into_iter().map(|p| Circle::new(p, 1, RED.filled())))?;
    }

    root.present()?;
    Ok(())
}

fn main() -> Result<(), Box<dyn std::error::Error>> {
    let word_count = 4usize.pow(WORD_LENGTH as u32);
    let mut candidates: Vec<_> = (0..word_count)
        .map(|index| score_word(word_from_index(index)))
        .collect();
    // Fix the first quadrant choice so the gallery does not fill up with
    // rotations/reflections of essentially the same high-scoring outline.
    candidates.retain(|candidate| candidate.word[0] == 4);
    candidates.sort_by(|a, b| b.score.total_cmp(&a.score));
    candidates.truncate(GALLERY_SIZE);

    for candidate in &candidates {
        println!(
            "word={:?} score={:.4} concavity={:.4} radial={:.4} cusp={:.4}",
            candidate.word,
            candidate.score,
            candidate.concavity,
            candidate.radial_variation,
            candidate.cusp_strength,
        );
    }
    render_gallery(&candidates)?;
    println!("Gallery written to plots/unusual_words_gallery.png");
    Ok(())
}
