use landau::edt::edt_1d_squared;

#[test]
fn single_source_parabola() {
    // f has a single zero (source) at index k
    let n = 11usize;
    let k = 4usize;
    let mut f = vec![1.0e20; n];
    f[k] = 0.0;

    let g = edt_1d_squared(&f);
    for i in 0..n {
        let ref_val = (i as f64 - k as f64).powi(2);
        assert!(
            (g[i] - ref_val).abs() < 1e-9,
            "i={i}, got={}, want={ref_val}",
            g[i]
        );
    }
}
