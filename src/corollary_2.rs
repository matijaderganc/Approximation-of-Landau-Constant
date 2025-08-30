use crate::covering_grids::covering_grid_bitmap;
use crate::dyadic::{ComplexDyadic, Dyadic};
use crate::edt::landau_l_via_edt_from_bitmap;
use crate::psi::mu_second;
use crate::holomorphic::ComplexFunction; 


pub fn calculate_epsilon(r : f64, r_hat : Dyadic, rho : Dyadic) -> Dyadic {
    let  c = 1.0 + (2.0f64).powi(-10);

    let mu = mu_second(&c, &r_hat) ;
    let r_hat_f = r_hat.to_f64();

    let cap_from_rho_f = rho.to_f64() / (4.0 * mu) ; // Δ ≤ ρ / (4 μ'')
    
    let cap_from_r_f = (r_hat_f - r) / 2.0;                 // > 0
   
    let mut cap_f = cap_from_rho_f.min(cap_from_r_f);
    assert!(cap_f.is_finite() && cap_f > 0.0, "No positive Δ available; check rho, r_hat, and mu_second.");

    let mut e: i32 = cap_f.log2().floor() as i32;
    let mut delta = Dyadic::new(1, e).reduce();


    // let cap_gap_strict_f = cap_from_r_f; // enforce: 2Δ < gap  ⇔  Δ < gap/2
    // if mu > 0.0 {
    //     // ensure both constraints together
    //     while delta.to_f64() > cap_from_rho_f || 2.0 * delta.to_f64() >= cap_gap_strict_f {
    //         e -= 1;
    //         delta = Dyadic::new(1, e).reduce();
    //     }
    // } else {
    //     // only the gap constraint
    //     while 2.0 * delta.to_f64() >= cap_gap_strict_f {
    //         e -= 1;
    //         delta = Dyadic::new(1, e).reduce();
    //     }
    // }
    let eps = (rho * delta * Dyadic::new(1, -4)).reduce();
    eps  
}

fn dyadic_floor_bits(x: f64, prec_bits: i32) -> Dyadic {
    if !x.is_finite() || x <= 0.0 {
        return Dyadic::zero();
    }
    let scale: i128 = 1i128 << prec_bits;
    let n = (x * (scale as f64)).floor() as i128;
    Dyadic::new(n, -prec_bits).reduce()
}

/// Build a ComplexDyadic point z = r_hat * e^{iθ}, rounding each coordinate
/// down to dyadic with `frac_bits` fractional bits.
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

/// NAIVE version:
/// - choose r_hat = r + min((1 - r)/8, 2^{-6})
/// - sample f' at N points on |z| = r_hat
/// - set rho = min_j |f'(z_j)| (rounded down to dyadic)
///
/// Returns (r_hat, rho).
pub fn certify_r_hat_rho(
    r: Dyadic,
    fprime: &ComplexFunction,
    n_samples: usize,   // e.g. 512 or 1024
    coord_frac_bits: i32, // dyadic rounding for z, e.g. 50
    rho_frac_bits: i32,   // dyadic rounding for rho, e.g. 40
) -> (Dyadic, Dyadic) {
    assert!(r > Dyadic::zero() && r < Dyadic::new(1, 0), "r must be in (0,1)");
    let one = Dyadic::new(1, 0);

    // Pick r_hat slightly above r
    let mut step = ((one - r) * Dyadic::new(1, -3)).reduce(); // (1 - r)/8
    let fallback = Dyadic::new(1, -6);                        // 2^-6
    if step <= Dyadic::zero() || step > fallback { step = fallback; }

    let mut r_hat = (r + step).reduce();
    if r_hat >= one { r_hat = (one - Dyadic::new(1, -12)).reduce(); }

    // Sample f' on |z| = r_hat
    let n = n_samples.max(8);
    let dtheta = std::f64::consts::TAU / (n as f64);
    let mut min_mag = f64::INFINITY;

    for j in 0..n {
        let theta = (j as f64) * dtheta;
        let z = point_on_circle(r_hat, theta, coord_frac_bits);
        let val = fprime.eval(&z); // ComplexDyadic
        let mag = val.abs();
        if mag < min_mag { min_mag = mag; }
    }

    // Naive rho = min magnitude (rounded down to dyadic)
    let rho = dyadic_floor_bits(min_mag, rho_frac_bits);
    (r_hat, rho)
}

    