use landau::{covering_grids::{create_covering_grid, unit_disk_n}, dyadic::Dyadic, holomorphic::{BoundingSequence, ComplexFunction, ExpansionCoefficients}, plot::{plot_covering_grid, plot_set, plot_set_old}, psi::{psi_infinity, t_vector}};

#[tokio::main(flavor = "multi_thread", worker_threads = 8)]

async fn main() -> anyhow::Result<()> {
    landau::ui::run_server().await
}

// fn main() -> Result<(), Box<dyn std::error::Error>> {
//     let m_seq = vec![Dyadic::new(1, -1), Dyadic::new(1, -1), Dyadic::new(1, -1)] ;
//     let t_seq = t_vector(10) ;
//     let domain = unit_disk_n(-4) ;
//     // plot_set_old(&domain, "ppt_domain.png") ;
//     let word = [1, 3, 1, 3, 1, 1] ;
//     let holo = psi_infinity(&m_seq, &t_seq, &word);
//     let f_prime = ComplexFunction::new(
//         BoundingSequence::new(m_seq),
//         ExpansionCoefficients::new(holo));
//     let mut im: Vec<_> = Vec::new() ;
//     let f = f_prime.antiderivative();
//     for x in &domain {
//         im.push(f.eval(&x)) }
//     plot_set_old(&im, "ppt_image.png") ;
//     let grid = create_covering_grid(&im, Dyadic::new(1, -3), Dyadic::new(1, -5));
//     plot_covering_grid(&grid, "ppt_grid.png")
// }