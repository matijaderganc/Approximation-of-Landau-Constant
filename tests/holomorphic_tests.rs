use landau::dyadic::{ComplexDyadic, Dyadic};
use landau::holomorphic::{BoundingSequence, ComplexFunction, ExpansionCoefficients};

#[test]
fn derivative_basic_polynomial() {
    // f(z) = z + 2z^2 + 3z^3  =>  f'(z) = 1 + 4z + 9z^2
    let f = ExpansionCoefficients {
        vector: vec![
            ComplexDyadic::zero(),      // c0
            ComplexDyadic::from_i64(1), // c1
            ComplexDyadic::from_i64(2), // c2
            ComplexDyadic::from_i64(3), // c3
        ],
    };

    let fprime = f.derivative();

    let expected = [
        ComplexDyadic::from_i64(1), // 1
        ComplexDyadic::from_i64(4), // 4 z
        ComplexDyadic::from_i64(9), // 9 z^2
    ];

    for (got, exp) in fprime.vector.iter().zip(expected.iter()) {
        let (gr, gi) = got.to_f64();
        let (er, ei) = exp.to_f64();
        assert!((gr - er).abs() < 1e-9 && (gi - ei).abs() < 1e-9);
    }
}

#[test]
fn derivative_basic() {
    // f(z) = z + 2z^2 + 3z^3  =>  f'(z) = 1 + 4z + 9z^2
    let f = ExpansionCoefficients {
        vector: vec![
            ComplexDyadic::zero(),      // c0
            ComplexDyadic::from_i64(1), // c1
            ComplexDyadic::from_i64(2), // c2
            ComplexDyadic::from_i64(3), // c3
        ],
    };
    let fp = f.derivative();
    assert_eq!(
        fp.vector,
        vec![
            ComplexDyadic::from_i64(1), // 1
            ComplexDyadic::from_i64(4), // 4
            ComplexDyadic::from_i64(9), // 9
        ]
    );
}

#[test]
fn antiderivative_basic() {
    // f'(z) = 1 + 2z + 3z^2  =>  f(z) = z + z^2 + z^3
    let fprime = ExpansionCoefficients {
        vector: vec![
            ComplexDyadic::from_i64(1),
            ComplexDyadic::from_i64(2),
            ComplexDyadic::from_i64(3),
        ],
    };
    let f = fprime.antiderivative(60); // 60 dyadic bits for safety

    // expected coefficients
    let expected = vec![
        ComplexDyadic::zero(),
        ComplexDyadic::from_i64(1),
        ComplexDyadic::from_i64(1),
        ComplexDyadic::from_i64(1),
    ];

    // compare each coefficient with tolerance
    for (got, exp) in f.vector.iter().zip(expected.iter()) {
        let (gr, gi) = got.to_f64();
        let (er, ei) = exp.to_f64();
        assert!(
            (gr - er).abs() < 1e-9 && (gi - ei).abs() < 1e-9,
            "mismatch: got ({gr}, {gi}), expected ({er}, {ei})"
        );
    }

    // and derivative gets us back (tolerance again)
    let back = f.derivative();
    for (got, exp) in back.vector.iter().zip(fprime.vector.iter()) {
        let (gr, gi) = got.to_f64();
        let (er, ei) = exp.to_f64();
        assert!(
            (gr - er).abs() < 1e-9 && (gi - ei).abs() < 1e-9,
            "mismatch in derivative: got ({gr}, {gi}), expected ({er}, {ei})"
        );
    }
}

#[test]
fn derivative_linearity() {
    // g(z) = (1 + 2z + 3z^2) + (2 + 1z + 0z^2)
    let g1 = ExpansionCoefficients {
        vector: vec![
            ComplexDyadic::from_i64(1),
            ComplexDyadic::from_i64(2),
            ComplexDyadic::from_i64(3),
        ],
    };
    let g2 = ExpansionCoefficients {
        vector: vec![
            ComplexDyadic::from_i64(2),
            ComplexDyadic::from_i64(1),
            ComplexDyadic::from_i64(0),
        ],
    };
    let g = ExpansionCoefficients {
        vector: g1
            .vector
            .iter()
            .zip(g2.vector.iter())
            .map(|(a, b)| *a + *b)
            .collect(),
    };

    let gp = g.derivative();
    let g1p = g1.derivative();
    let g2p = g2.derivative();
    for (got, exp) in gp.vector.iter().zip(
        g1p.vector
            .iter()
            .zip(g2p.vector.iter())
            .map(|(a, b)| *a + *b),
    ) {
        let (gr, gi) = got.to_f64();
        let (er, ei) = exp.to_f64();
        assert!((gr - er).abs() < 1e-9 && (gi - ei).abs() < 1e-9);
    }
}

#[test]
fn eval_constant_polynomial() {
    // f(z) = 3
    let coeffs = ExpansionCoefficients {
        vector: vec![ComplexDyadic::from_i64(3)],
    };
    let one = Dyadic::new(1, 0);
    let bs = BoundingSequence::new(vec![one, one, one, one]);
    let f = ComplexFunction::new(bs, coeffs);

    // f(any z) = 3
    let z = ComplexDyadic::new(Dyadic::new(5, -3), Dyadic::new(-1, -4));
    let (re, im) = f.eval(&z).to_f64();
    assert!((re - 3.0).abs() < 1e-12 && im.abs() < 1e-12);
}

#[test]
fn derivative_and_antiderivative_inverse() {
    // Start with f'(z) = 1 + 2z + 3z^2
    let fprime = ExpansionCoefficients {
        vector: vec![
            ComplexDyadic::from_i64(1),
            ComplexDyadic::from_i64(2),
            ComplexDyadic::from_i64(3),
        ],
    };

    // Integrate then differentiate back. (antiderivative may pick c0=0)
    let one = Dyadic::new(1, 0);
    let bs = BoundingSequence::new(vec![one, one, one, one]);
    let f = ComplexFunction::new(bs, fprime).antiderivative();
    let back = f.derivative();

    for (got, exp) in back.expansion_coefficients.vector.iter().zip(
        [
            ComplexDyadic::from_i64(1),
            ComplexDyadic::from_i64(2),
            ComplexDyadic::from_i64(3),
        ]
        .iter(),
    ) {
        let (gr, gi) = got.to_f64();
        let (er, ei) = exp.to_f64();
        assert!((gr - er).abs() < 1e-9 && (gi - ei).abs() < 1e-9);
    }
}

#[test]
fn complex_function_eval_and_derivative_consistency() {
    // f(z) = 1 + 2z + 3z^2
    let f_coeffs = ExpansionCoefficients {
        vector: vec![
            ComplexDyadic::from_i64(1),
            ComplexDyadic::from_i64(2),
            ComplexDyadic::from_i64(3),
        ],
    };
    // Bounding sequence values are not used by eval/derivative algebra itself here,
    // just provide a small vector of positives.
    let bs = BoundingSequence {
        vector: vec![Dyadic::from_i64(1); 4],
    };

    let f = ComplexFunction::new(bs.clone(), f_coeffs.clone());
    let fp = f.derivative();
    let got: Vec<(f64, f64)> = fp
        .expansion_coefficients
        .vector
        .iter()
        .map(|c| c.to_f64())
        .collect();
    let want = [(2.0, 0.0), (6.0, 0.0)];
    for (g, w) in got.iter().zip(want.iter()) {
        assert!((g.0 - w.0).abs() < 1e-12 && (g.1 - w.1).abs() < 1e-12);
    }

    // Check Horner eval against direct polynomial at z = 1/2 + i/4
    let z = ComplexDyadic::new(Dyadic::new(1, -1), Dyadic::new(1, -2));
    let (zr, zi) = z.to_f64();
    let (fr, fi) = f.eval(&z).to_f64();
    let (er, ei) = {
        // 1 + 2z + 3z^2 in f64
        let z2r = zr * zr - zi * zi;
        let z2i = 2.0 * zr * zi;
        let rr = 1.0 + 2.0 * zr + 3.0 * z2r;
        let ii = 0.0 + 2.0 * zi + 3.0 * z2i;
        (rr, ii)
    };
    assert!((fr - er).abs() < 1e-12 && (fi - ei).abs() < 1e-12);
}
