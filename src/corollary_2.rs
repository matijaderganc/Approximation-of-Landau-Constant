use crate::dyadic::{ComplexDyadic, Dyadic};
use crate::holomorphic::ComplexFunction;
use crate::psi::mu_second;

/// From given r, r_hat and rho, we calculate epsilon as described in Lemma 6 of the article
pub fn calculate_epsilon(r: f64, r_hat: Dyadic, rho: Dyadic) -> Dyadic {
    let c = 1.0 + (2.0f64).powi(-10);

    let mu = mu_second(&c, &r_hat);
    let r_hat_f = r_hat.to_f64();

    let cap_from_rho_f = rho.to_f64() / (4.0 * mu); // first condition for epsilon

    let cap_from_r_f = (r_hat_f - r) / 2.0; // second condition for epsilon

    let cap_f = cap_from_rho_f.min(cap_from_r_f);
    assert!(
        cap_f.is_finite() && cap_f > 0.0,
        "No positive Δ available; check rho, r_hat, and mu_second."
    );

    let e: i32 = cap_f.log2().floor() as i32;
    let delta = Dyadic::new(1, e).reduce();
    (rho * delta * Dyadic::new(1, -4)).reduce()
    // (rho * delta * Dyadic::new(1, -4)).reduce()

}

///  we floor dyadic numbers
fn dyadic_floor_bits(x: f64, prec_bits: i32) -> Dyadic {
    if !x.is_finite() || x <= 0.0 {
        return Dyadic::zero();
    }
    let scale: i128 = 1i128 << prec_bits;
    let n = (x * (scale as f64)).floor() as i128;
    Dyadic::new(n, -prec_bits).reduce()
}

/// Build a ComplexDyadic point z = r_hat * e^{iθ}, rounding each coordinate
/// down to dyadic with `frac_bits` precision. This is used to calculate rho from r_hat.
fn point_on_circle(r_hat: Dyadic, theta: f64, frac_bits: i32) -> ComplexDyadic {
    let rc = r_hat.to_f64();
    let x = rc * theta.cos();
    let y = rc * theta.sin();
    let sx = if x < 0.0 { -1 } else { 1 };
    let sy = if y < 0.0 { -1 } else { 1 };
    ComplexDyadic::new(
        dyadic_floor_bits(x.abs(), frac_bits) * Dyadic::new(sx, 0),
        dyadic_floor_bits(y.abs(), frac_bits) * Dyadic::new(sy, 0),
    )
}

/// we calculate r_hat and rho from a given input r, we have to justify the conditions described in Lemma 6.
pub fn certify_r_hat_rho(
    r: Dyadic,
    fprime: &ComplexFunction,
    n_samples: usize,     // samples around the circle
    coord_frac_bits: i32, // dyadic rounding for z
    rho_frac_bits: i32,   // dyadic rounding for rho
) -> (Dyadic, Dyadic) {
    assert!(
        r > Dyadic::zero() && r < Dyadic::new(1, 0),
        "r must be in (0,1)"
    );
    let one = Dyadic::new(1, 0);

    // Pick r_hat slightly above r
    let mut step = ((one - r) * Dyadic::new(1, -3)).reduce();
    let fallback = Dyadic::new(1, -6);
    if step <= Dyadic::zero() || step > fallback {
        step = fallback;
    }

    let mut r_hat = (r + step).reduce();
    if r_hat >= one {
        r_hat = (one - Dyadic::new(1, -12)).reduce();
    }

    // sample f' on |z| = r_hat
    let n = n_samples.max(8);
    let dtheta = std::f64::consts::TAU / (n as f64);
    let mut min_mag = f64::INFINITY;

    for j in 0..n {
        let theta = (j as f64) * dtheta;
        let z = point_on_circle(r_hat, theta, coord_frac_bits);
        let val = fprime.eval(&z);
        let mag = val.abs();
        if mag < min_mag {
            min_mag = mag;
        }
    }

    let rho = dyadic_floor_bits(min_mag, rho_frac_bits);
    (r_hat, rho)
}
