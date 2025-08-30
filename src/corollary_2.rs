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
    if !(x.is_finite()) || x <= 0.0 {
        return Dyadic::zero();
    }
    let scale: i128 = 1i128 << prec_bits;
    let n = (x * (scale as f64)).floor() as i128;
    Dyadic::new(n, -prec_bits).reduce()
}

/// Build a ComplexDyadic from polar (r_hat * e^{iθ}) with dyadic rounding.
fn point_on_circle(r_hat: Dyadic, theta: f64, frac_bits: i32) -> ComplexDyadic {
    let rc = r_hat.to_f64();
    let x = rc * theta.cos();
    let y = rc * theta.sin();
    ComplexDyadic::new(
        dyadic_floor_bits(x.abs(), frac_bits) * Dyadic::new(if x < 0.0 { -1 } else { 1 }, 0),
        dyadic_floor_bits(y.abs(), frac_bits) * Dyadic::new(if y < 0.0 { -1 } else { 1 }, 0),
    )
}

/// Certified lower bound ρ on |f'(z)| over |z|=r_hat, using:
///   ρ  = min_j {|f'(z_j)|} - μ''_{r̂} * r_hat * Δθ
/// where z_j = r_hat * e^{2π i j / N}, and Δθ = 2π/N.
/// Increases N if needed; if still non-positive, moves r_hat closer to r.
pub fn certify_r_hat_rho(
    r: Dyadic,                   // base radius (dyadic in (0,1))
    fprime: &ComplexFunction,    // f' (your existing derivative object)
    target_margin_bits: i32,     // how finely to round ρ down to dyadic (e.g. 40)
) -> (Dyadic, Dyadic) {
    assert!(r > Dyadic::zero() && r < Dyadic::new(1,0), "r must be in (0,1)");
    let one = Dyadic::new(1,0);

    // Start slightly above r, safely below 1.
    let mut step = ((one - r) * Dyadic::new(1, -3)).reduce(); // (1-r)/8
    let mut r_hat = (r + step).reduce();
    if r_hat >= one { r_hat = (one - Dyadic::new(1, -12)).reduce(); }

    // Outer loop: if we fail to get ρ>0, move r̂ closer to r by halving the step.
    let min_step = Dyadic::new(1, -60);

    'outer: loop {
        // μ''_{r̂} from your psi crate (Lemma 6 uses an upper bound on sup|f''| on D_{r̂})
        let c = 1.0 + (2.0f64).powi(-100);
        let mu2 = mu_second(&c, &r_hat); // Dyadic upper bound on sup_{|z|<=r̂} |f''(z)|

        // Try increasing angular resolutions N until Lipschitz-corrected min is > 0
        let mut n: usize = 256;
        while n <= 1 << 16 {
            let delta_theta = std::f64::consts::TAU / (n as f64);
            // Lipschitz radius on the circle between samples: ≤ μ'' * r̂ * Δθ
            let lip = mu2 * r_hat.to_f64() * delta_theta;

            // Sample f' on |z|=r_hat
            let mut min_mag = std::f64::INFINITY;
            for j in 0..n {
                let theta = (j as f64) * delta_theta;
                let z = point_on_circle(r_hat, theta, /*frac_bits*/ 50);
                let val = fprime.eval(&z); // ComplexDyadic
                let mag = (val.re.to_f64().powi(2) + val.im.to_f64().powi(2)).sqrt();
                if mag < min_mag { min_mag = mag; }
            }

            // Certified uniform lower bound over the whole circle
            let rho_f64 = min_mag - lip;
            if rho_f64 > 0.0 {
                // Round down to a dyadic strictly below rho_f64
                let rho = dyadic_floor_bits(rho_f64, target_margin_bits);
                if rho > Dyadic::zero() {
                    return (r_hat, rho);
                }
            }
            n *= 2; // refine sampling
        }

        // If we get here, even with fine sampling and this r̂ we couldn't certify ρ>0.
        // Move r̂ closer to r and try again.
        step = (step * Dyadic::new(1, -1)).reduce();
        assert!(step > min_step, "Could not certify a positive ρ; try smaller r̂ or adjust bounds.");
        r_hat = (r + step).reduce();
        if r_hat >= one { r_hat = (one - Dyadic::new(1, -12)).reduce(); }
    }
}

    