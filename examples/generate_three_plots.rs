use landau::{
    covering_grids::{create_covering_grid, grid_complement, unit_disk_n},
    dyadic::Dyadic,
    holomorphic::{
        BoundingSequence, ComplexFunction, ExpansionCoefficients,
    },
    plot::{plot_covering_grid, plot_grid_complement, plot_set_old},
    psi::{psi_infinity, t_vector},
};

fn main() -> Result<(), Box<dyn std::error::Error>> {
    // 1. Sample D ∩ 2^(-n) Z².
    //
    // Important: unit_disk_n expects the exponent itself, so -6 means
    // that neighboring points are 2^-6 = 1/64 apart.
    let disk_accuracy = -4;
    let domain = unit_disk_n(disk_accuracy);

    plot_set_old(
        &domain,
        "plots/01_unit_disk_grid.png",
    )?;

    // 2. Construct a normalized holomorphic function.
    //
    // psi_infinity constructs the coefficients of f'. Its constant
    // coefficient is 1, and integration therefore gives:
    //     f(0) = 0, f'(0) = 1.
    let m_seq = vec![
        Dyadic::new(1, -1), // 1/2
        Dyadic::new(1, -1), // 1/2
        Dyadic::new(1, -1), // 1/2
    ];

    let word = [1u8, 4, 4, 4, 1, 4, 4, 3];
    let t_seq = t_vector(word.len());

    let derivative_coefficients =
        psi_infinity(&m_seq, &t_seq, &word);

    let f_prime = ComplexFunction::new(
        BoundingSequence::new(m_seq),
        ExpansionCoefficients::new(derivative_coefficients),
    );

    let f = f_prime.antiderivative();

    // Apply f to every point in the original grid.
    let image: Vec<_> = domain
        .iter()
        .map(|z| f.eval(z))
        .collect();

    plot_set_old(
        &image,
        "plots/02_holomorphic_image.png",
    )?;

    // 3. Construct an epsilon-covering grid of f(D).
    //
    // epsilon controls the neighborhood around the sampled image.
    // delta is the spacing of the output lattice. The repository
    // normally uses delta = epsilon / 4.
    let epsilon = Dyadic::new(1, -3); // 1/8
    let delta = epsilon * Dyadic::new(1, -2); // 1/32

    let image_grid =
        create_covering_grid(&image, epsilon, delta);

    plot_covering_grid(
        &image_grid,
        "plots/03_image_epsilon_grid.png",
    )?;

    // 4. Plot the complement of the epsilon-covering grid inside its
    // padded rectangular lattice bounding box.
    let complement = grid_complement(&image_grid, delta);
    plot_grid_complement(
        &complement,
        "plots/04_grid_complement.png",
    )?;

    println!("Domain samples: {}", domain.len());
    println!("Image samples: {}", image.len());
    println!("Covering-grid points: {}", image_grid.len());
    println!("Complement points: {}", complement.len());

    Ok(())
}
