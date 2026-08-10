use std::{
    alloc::{GlobalAlloc, Layout, System},
    hint::black_box,
    mem::size_of_val,
    sync::atomic::{AtomicUsize, Ordering},
    time::{Duration, Instant},
};

use landau::{
    covering_grids::{covering_grid_bitmap, create_covering_grid, unit_disk_n},
    dyadic::Dyadic,
    edt::landau_l_via_edt_from_bitmap,
    holomorphic::{BoundingSequence, ComplexFunction, ExpansionCoefficients},
    psi::{psi_infinity, t_vector},
};
use rand::{Rng, SeedableRng, rngs::StdRng};

// These settings mirror the main program's usual delta = epsilon / 4 rule.
// With a unit-disk domain they normally produce dense grids containing roughly
// one to several million lattice cells, depending on the random function.
const SAMPLE_COUNT: usize = 8;
const EDT_REPEATS_PER_SAMPLE: usize = 3;
const WORD_LENGTH: usize = 12;
const DOMAIN_ACCURACY: i32 = -7; // domain spacing = 2^-7
const EPSILON_EXPONENT: i32 = -7; // epsilon = 2^-7
const DELTA_EXPONENT: i32 = -9; // delta = 2^-9 = epsilon / 4
const RANDOM_SEED: u64 = 0x4c41_4e44_4155;

/// Counts requested live heap bytes in this benchmark process. This measures
/// allocations retained by a grid, rather than transient peak allocations or
/// the operating system allocator's page-level overhead.
struct TrackingAllocator;

static LIVE_HEAP_BYTES: AtomicUsize = AtomicUsize::new(0);

#[global_allocator]
static GLOBAL_ALLOCATOR: TrackingAllocator = TrackingAllocator;

unsafe impl GlobalAlloc for TrackingAllocator {
    unsafe fn alloc(&self, layout: Layout) -> *mut u8 {
        let ptr = System.alloc(layout);
        if !ptr.is_null() {
            LIVE_HEAP_BYTES.fetch_add(layout.size(), Ordering::Relaxed);
        }
        ptr
    }

    unsafe fn alloc_zeroed(&self, layout: Layout) -> *mut u8 {
        let ptr = System.alloc_zeroed(layout);
        if !ptr.is_null() {
            LIVE_HEAP_BYTES.fetch_add(layout.size(), Ordering::Relaxed);
        }
        ptr
    }

    unsafe fn dealloc(&self, ptr: *mut u8, layout: Layout) {
        System.dealloc(ptr, layout);
        LIVE_HEAP_BYTES.fetch_sub(layout.size(), Ordering::Relaxed);
    }

    unsafe fn realloc(&self, ptr: *mut u8, old_layout: Layout, new_size: usize) -> *mut u8 {
        let new_ptr = System.realloc(ptr, old_layout, new_size);
        if !new_ptr.is_null() {
            if new_size >= old_layout.size() {
                LIVE_HEAP_BYTES.fetch_add(new_size - old_layout.size(), Ordering::Relaxed);
            } else {
                LIVE_HEAP_BYTES.fetch_sub(old_layout.size() - new_size, Ordering::Relaxed);
            }
        }
        new_ptr
    }
}

#[derive(Debug)]
struct SampleResult {
    complex_points: usize,
    complex_grid_bytes: usize,
    bitmap_width: usize,
    bitmap_height: usize,
    bitmap_cells: usize,
    bitmap_inside_cells: usize,
    bitmap_bytes: usize,
    average_edt_time: Duration,
    edt_cells_per_second: f64,
}

fn random_word(rng: &mut StdRng) -> Vec<u8> {
    (0..WORD_LENGTH).map(|_| rng.gen_range(1..=4)).collect()
}

fn function_for_word(word: &[u8]) -> ComplexFunction {
    // This is the bounded sequence used by the interactive single-word path.
    let m_seq = vec![
        Dyadic::new(1, -1),
        Dyadic::new(1, -1),
        Dyadic::new(1, -2),
        Dyadic::new(1, -2),
    ];
    let coefficients = psi_infinity(&m_seq, &t_vector(word.len()), word);
    let f_prime = ComplexFunction::new(
        BoundingSequence::new(m_seq),
        ExpansionCoefficients::new(coefficients),
    );
    f_prime.antiderivative()
}

fn retained_bytes_since(before: usize, value_stack_bytes: usize) -> usize {
    LIVE_HEAP_BYTES
        .load(Ordering::Relaxed)
        .saturating_sub(before)
        + value_stack_bytes
}

fn benchmark_sample(
    domain: &[landau::dyadic::ComplexDyadic],
    word: &[u8],
    epsilon: Dyadic,
    delta: Dyadic,
) -> SampleResult {
    let f = function_for_word(word);
    let image: Vec<_> = domain.iter().map(|z| f.eval(z)).collect();

    // Measure the retained representation of HashSet<ComplexDyadic>.
    let before_complex_grid = LIVE_HEAP_BYTES.load(Ordering::Relaxed);
    let complex_grid = create_covering_grid(&image, epsilon, delta);
    let complex_grid_bytes = retained_bytes_since(before_complex_grid, size_of_val(&complex_grid));
    let complex_points = complex_grid.len();
    black_box(&complex_grid);
    drop(complex_grid);

    // Also report the retained bitmap size, since this is the representation
    // consumed by the current EDT implementation.
    let before_bitmap = LIVE_HEAP_BYTES.load(Ordering::Relaxed);
    let bitmap = covering_grid_bitmap(&image, epsilon, delta);
    let bitmap_bytes = retained_bytes_since(before_bitmap, size_of_val(&bitmap));
    let bitmap_cells = bitmap.width * bitmap.height;
    let bitmap_inside_cells = bitmap.data.iter().filter(|&&cell| cell != 0).count();

    // Warm allocator and instruction caches before collecting timings.
    black_box(landau_l_via_edt_from_bitmap(&bitmap, delta));

    let mut elapsed = Duration::ZERO;
    for _ in 0..EDT_REPEATS_PER_SAMPLE {
        let started = Instant::now();
        black_box(landau_l_via_edt_from_bitmap(
            black_box(&bitmap),
            delta,
        ));
        elapsed += started.elapsed();
    }
    let average_edt_time = elapsed / EDT_REPEATS_PER_SAMPLE as u32;
    let edt_cells_per_second = bitmap_cells as f64 / average_edt_time.as_secs_f64();

    SampleResult {
        complex_points,
        complex_grid_bytes,
        bitmap_width: bitmap.width,
        bitmap_height: bitmap.height,
        bitmap_cells,
        bitmap_inside_cells,
        bitmap_bytes,
        average_edt_time,
        edt_cells_per_second,
    }
}

fn mean(values: &[f64]) -> f64 {
    values.iter().sum::<f64>() / values.len() as f64
}

fn standard_deviation(values: &[f64]) -> f64 {
    let average = mean(values);
    (values
        .iter()
        .map(|value| (value - average).powi(2))
        .sum::<f64>()
        / values.len() as f64)
        .sqrt()
}

fn mib(bytes: usize) -> f64 {
    bytes as f64 / (1024.0 * 1024.0)
}

fn main() {
    assert!(SAMPLE_COUNT > 0);
    assert!(EDT_REPEATS_PER_SAMPLE > 0);

    let epsilon = Dyadic::new(1, EPSILON_EXPONENT);
    let delta = Dyadic::new(1, DELTA_EXPONENT);
    assert_eq!(delta, epsilon * Dyadic::new(1, -2));

    println!("Landau grid/EDT benchmark");
    println!("Run with --release for meaningful timing results.");
    println!(
        "samples={SAMPLE_COUNT}, EDT repeats/sample={EDT_REPEATS_PER_SAMPLE}, \
         domain step=2^{DOMAIN_ACCURACY}, epsilon=2^{EPSILON_EXPONENT}, \
         delta=2^{DELTA_EXPONENT}"
    );

    let domain = unit_disk_n(DOMAIN_ACCURACY);
    println!("domain points={}\n", domain.len());

    let mut rng = StdRng::seed_from_u64(RANDOM_SEED);
    let mut results = Vec::with_capacity(SAMPLE_COUNT);
    for sample_index in 0..SAMPLE_COUNT {
        let word = random_word(&mut rng);
        let result = benchmark_sample(&domain, &word, epsilon, delta);
        println!(
            "sample {:>2}: word={:?}\n\
             complex grid: {:>9} points, {:>8.2} MiB, {:>7.2} bytes/point\n\
             EDT bitmap:  {:>4} x {:<4} = {:>9} cells ({:>9} inside), \
             {:>6.2} MiB, {:>5.2} bytes/cell\n\
             EDT: {:>8.3} ms/run, {:>10.0} cells/s\n",
            sample_index + 1,
            word,
            result.complex_points,
            mib(result.complex_grid_bytes),
            result.complex_grid_bytes as f64 / result.complex_points as f64,
            result.bitmap_width,
            result.bitmap_height,
            result.bitmap_cells,
            result.bitmap_inside_cells,
            mib(result.bitmap_bytes),
            result.bitmap_bytes as f64 / result.bitmap_cells as f64,
            result.average_edt_time.as_secs_f64() * 1_000.0,
            result.edt_cells_per_second,
        );
        results.push(result);
    }

    let bytes_per_complex_point: Vec<_> = results
        .iter()
        .map(|r| r.complex_grid_bytes as f64 / r.complex_points as f64)
        .collect();
    let edt_times_ms: Vec<_> = results
        .iter()
        .map(|r| r.average_edt_time.as_secs_f64() * 1_000.0)
        .collect();
    let edt_rates: Vec<_> = results.iter().map(|r| r.edt_cells_per_second).collect();
    let bitmap_bytes_per_cell: Vec<_> = results
        .iter()
        .map(|r| r.bitmap_bytes as f64 / r.bitmap_cells as f64)
        .collect();
    let total_complex_bytes: usize = results.iter().map(|r| r.complex_grid_bytes).sum();
    let total_complex_points: usize = results.iter().map(|r| r.complex_points).sum();
    let total_edt_cells: usize = results.iter().map(|r| r.bitmap_cells).sum();
    let total_edt_seconds: f64 = results
        .iter()
        .map(|r| r.average_edt_time.as_secs_f64())
        .sum();

    println!("Summary across {SAMPLE_COUNT} random samples");
    println!(
        "complex-grid RAM: {:.2} bytes/point weighted average \
         ({:.2} mean, {:.2} standard deviation)",
        total_complex_bytes as f64 / total_complex_points as f64,
        mean(&bytes_per_complex_point),
        standard_deviation(&bytes_per_complex_point),
    );
    println!(
        "EDT time:         {:.3} ms/grid average ({:.3} ms standard deviation)",
        mean(&edt_times_ms),
        standard_deviation(&edt_times_ms),
    );
    println!(
        "bitmap RAM:       {:.3} bytes/cell average ({:.3} standard deviation)",
        mean(&bitmap_bytes_per_cell),
        standard_deviation(&bitmap_bytes_per_cell),
    );
    println!(
        "EDT throughput:   {:.0} cells/s aggregate ({:.0} mean, {:.0} standard deviation)",
        total_edt_cells as f64 / total_edt_seconds,
        mean(&edt_rates),
        standard_deviation(&edt_rates),
    );
    println!("Note: one 'cell' means one point in the dense bitmap for one complete 2D EDT.");
}
