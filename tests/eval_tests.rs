use landau::corollary_2::calculate_epsilon;
use landau::covering_grids::{covering_grid_bitmap, unit_disk_n};
use landau::holomorphic::{BoundingSequence, ComplexFunction, ExpansionCoefficients};
use landau::dyadic::{ComplexDyadic, Dyadic};
use landau::edt::landau_l_via_edt_from_bitmap;
use landau::psi::{m_vec, psi_infinity, t_vector};

#[test]
    fn test_epsilon() {
        let r = 0.95;
        let r_hat = Dyadic::new(127, -7);
        println!("{}", r_hat.to_f64());
        let rho = Dyadic::new(1, -4);
        let eps = calculate_epsilon(r, r_hat, rho);
        assert_eq!(eps, Dyadic::new(1, -37));
    }

#[test]
fn test_approx() {
    let m_seq = m_vec(5);
    let delta = Dyadic::new(1, -7);
    let t_seq = t_vector(100);
    let test_word = [1, 2, 3, 4, 1, 1];
    let coeffs = psi_infinity(&m_seq, &t_seq, &test_word);
    let fprime = ComplexFunction::new(
        BoundingSequence::new(m_seq),
        ExpansionCoefficients::new(coeffs),
    );
    let f = fprime.antiderivative();

    let domain = unit_disk_n(-6);
    let mut img: Vec<ComplexDyadic> = Vec::with_capacity(domain.len());
    for &z in domain.iter() {
        img.push(f.eval(&z));
    }
    let eps = delta * Dyadic::new(1, 2);
    let grid_bitmap = covering_grid_bitmap(&img, eps, delta);
    let l_val = landau_l_via_edt_from_bitmap(&grid_bitmap, delta);
    assert_eq!(l_val, 0.9677649931279152);
}