use std::sync::{Arc, Mutex};

use axum::{
    Router,
    extract::{Form, State},
    http::{StatusCode, header},
    response::{Html, IntoResponse, Response},
    routing::{get, post},
};
use serde::Deserialize;
use tokio::{net::TcpListener, task};
use std::time::Instant;

use crate::{covering_grids::{covering_grid_bitmap, unit_disk_n}, dyadic::Dyadic};
// NEW: bring the sweep into scope
use crate::evaluation::approximate_all_words;
use crate::holomorphic::{BoundingSequence, ComplexFunction, ExpansionCoefficients};
use crate::plot::plot_set;
use crate::psi::{m_vec, psi_infinity, t_vector};
use crate::edt::{landau_l_via_edt_from_bitmap};

const CONTACT_EMAIL: &str = "landau@constant.com";

// CSS for all pages
fn common_css() -> &'static str {
    r#"
    :root { --bg:#0a0b10; --fg:#e6e9ef; --muted:#a0a7b4; --card:#131522; --brand:#7aa2ff; --accent:#99e2b4 }
    * { box-sizing: border-box }
    html, body { margin:0; padding:0; height:100%; font-family: system-ui, -apple-system, Segoe UI, Roboto, Inter, Arial, sans-serif; color:var(--fg); background:radial-gradient(1200px 800px at 20% -30%, #14172a 0%, #0b0d17 40%, #0a0b10 100%) fixed }
    a { color: var(--brand); text-decoration: none }
    a:hover { text-decoration: underline }
    .container { max-width: 980px; margin: 0 auto; padding: 2rem 1.25rem }
    .card { background: linear-gradient(180deg, rgba(255,255,255,.03), rgba(255,255,255,.01)); border:1px solid rgba(255,255,255,.06); border-radius: 16px; box-shadow: 0 10px 30px rgba(0,0,0,.25); padding: 1.25rem }
    h1, h2 { letter-spacing:.2px; }
    h1 { font-size: 1.6rem; margin: 0 0 1rem }
    h2 { font-size: 1.25rem; margin: 1.25rem 0 .5rem }
    p { color: var(--muted); line-height: 1.55 }
    code { background:#0f1120; padding:.12rem .35rem; border-radius:.35rem }
    label { display:block; margin:.5rem 0 .25rem }
    input[type="text"], input[type="number"] {
      width:100%; padding:.6rem .7rem; background:#0f1120; color:var(--fg);
      border:1px solid rgba(255,255,255,.08); border-radius:.5rem; outline:none;
    }
    input[type="number"] { max-width: 220px }
    button {
      display:inline-flex; align-items:center; gap:.5rem; padding:.6rem .9rem; border-radius:.55rem; border:1px solid rgba(255,255,255,.1);
      color:#0b0d17; background: linear-gradient(180deg, #adcbff, #7aa2ff); font-weight:600; cursor:pointer;
    }
    button:hover { filter: brightness(1.02) }
    hr { margin:1.75rem 0; border:none; border-top:1px solid rgba(255,255,255,.08) }
    .row { display:grid; grid-template-columns: 1fr 220px; gap:1rem; align-items:end }
    .nav { background: rgba(255,255,255,.02); border-bottom: 1px solid rgba(255,255,255,.06) }
    .nav-inner { max-width: 980px; margin:0 auto; padding:.8rem 1.25rem; display:flex; gap:1rem; align-items:center; justify-content:space-between }
    .nav a.brand { font-weight: 700; letter-spacing:.4px; color: var(--fg) }
    .nav .links { display:flex; gap:1rem }
    .footer { margin-top: 2rem; padding: 1rem 1.25rem; border-top: 1px solid rgba(255,255,255,.08); color: var(--muted) }
    .badge { display:inline-block; font-size:.75rem; padding:.2rem .5rem; border-radius:.4rem; background: rgba(153,226,180,.15); color: var(--accent); border:1px solid rgba(153,226,180,.25) }
    ul { margin: .5rem 0 .5rem 1rem }
    li { margin:.25rem 0 }
    figure { margin: 1rem 0; text-align: center }
    figcaption { color: var(--muted); font-size:.9rem; margin-top:.35rem }
    "#
}

// Navigation html code
fn nav_html() -> String {
    format!(
        r#"
    <div class="nav">
      <div class="nav-inner">
        <a class="brand" href="/">Landau Explorer</a>
        <div class="links">
          <a href="/">Home</a>
          <a href="/intro">Introduction</a>
        </div>
      </div>
    </div>
    "#
    )
}

fn escape_html(s: &str) -> String {
  let mut out = String::with_capacity(s.len());
  for ch in s.chars() {
      match ch {
          '&' => out.push_str("&amp;"),
          '<' => out.push_str("&lt;"),
          '>' => out.push_str("&gt;"),
          '"' => out.push_str("&quot;"),
          '\'' => out.push_str("&#39;"),
          _ => out.push(ch),
      }
  }
  out
}

// Footer with contact email
fn footer_html() -> String {
    format!(
        r#"
    <div class="container footer">
      <div>© Landau Explorer · Contact: <a href="mailto:{0}">{0}</a></div>
    </div>
    "#,
        CONTACT_EMAIL
    )
}

// Page layout wrapper
fn page_html(title: &str, inner: &str) -> String {
    format!(
        r#"<!doctype html>
<html>
  <head>
    <meta charset="utf-8" />
    <title>{title}</title>
    <meta name="viewport" content="width=device-width, initial-scale=1" />
    <style>{css}</style>
  </head>
  <body>
    {nav}
    <main class="container">
      {inner}
    </main>
    {footer}
  </body>
</html>"#,
        title = title,
        css = common_css(),
        nav = nav_html(),
        inner = inner,
        footer = footer_html()
    )
}


/// Shared app state: store the most recent PNG bytes
#[derive(Default)]
struct AppState {
    png: Mutex<Vec<u8>>,
}
type SharedState = Arc<AppState>;


#[derive(Deserialize)]
struct PlotForm {
    word: String,
    acc_n: Option<i32>, 
}

const DEFAULT_PLOT_ACC_N: i32 = 7;   // for /plot (step = 2^-n)
const DEFAULT_CALC_ALL_LEN: usize = 3;
const DEFAULT_CALC_ALL_STEP: i32 = 1;

fn default_calc_all_len() -> usize { DEFAULT_CALC_ALL_LEN }
fn default_calc_all_step() -> i32 { DEFAULT_CALC_ALL_STEP }
#[derive(Deserialize)]
struct CalcAllForm {
  #[serde(default = "default_calc_all_len")]
  length: usize,

  #[serde(default = "default_calc_all_step")]
  step: i32,
}

/// Both "1,2,3,4" or "1234" work fine
fn parse_word(input: &str) -> Vec<u8> {
    if input.contains(',') || input.contains(' ') {
        input
            .split(|c: char| c == ',' || c.is_whitespace())
            .filter(|s| !s.is_empty())
            .filter_map(|s| s.parse::<u8>().ok())
            .collect()
    } else {
        input
            .chars()
            .filter_map(|c| c.to_digit(10).map(|d| d as u8))
            .collect()
    }
}

async fn home() -> Html<String> {
    let inner = r#"
    <div class="card">
      <h1>Word → image of D<sub>1</sub></h1>
      <p>Enter a word over <code>{1,2,3,4}</code> to generate the derivative via <code>psi_infinity</code>, integrate to a normalised <code>f</code>, evaluate on a discrete unit disk, create a covering grid of the image and evaluate λ<sub>f</sub> using EDT, then render.</p>

      <form method="post" action="/plot">
        <div class="row">
          <div>
            <label>Word (e.g. <code>1,2,3,4</code> or <code>1234</code>)</label>
            <input type="text" name="word" placeholder="1,2,3,4" required />
          </div>
          <div>
            <label>Accuracy n (step = 2<sup>−n</sup>)</label>
            <input type="number" name="acc_n" min="6" max="18" value="7" />
          </div>
        </div>
        <p><button type="submit">Plot</button></p>
      </form>

      <hr/>

      <h2>Approximate the constant</h2>
      <p>Final algorithm: approximate Landau's constant on all words up to (including) selected length. Recieved approximation is up to λ*2^{-step} away from true value, theoretically proven. :)</p>
      <form method="post" action="/calc-all">
        <div class="row">
          <div>
            <label>Word length</label>
            <input type="number" name="length" min="1" max="10" value="3" />
          </div>
          <div>
            <label>Step (UI) → <code>disk_decrease</code></label>
            <input type="number" name="step" min="1" max="10" value="1" />
          </div>
        </div>
        <p><button type="submit">Run sweep</button></p>
      </form>
    </div>
  "#;

    Html(page_html("Landau Explorer — Home", inner))
}

async fn plot(State(state): State<SharedState>, Form(form): Form<PlotForm>) -> impl IntoResponse {
    // Parse inputs
    let word_vec = parse_word(&form.word);
    let acc_n_ui = form.acc_n.unwrap_or(DEFAULT_PLOT_ACC_N).clamp(3, 18);
    let epsilon = Dyadic::new(1, -acc_n_ui) ;
    let delta = epsilon * Dyadic::new(1, -2);
    let word = form.word.trim().to_string();
    let acc_n_display = match form.acc_n {
        Some(n) => n.to_string(),
        None => "(default)".to_string(),
    };

    // Build inputs for the blocking job
    let m_seq = Arc::new(vec![Dyadic::new(1, -1), Dyadic::new(1, -1), Dyadic::new(1, -2), Dyadic::new(1, -2), ]);
    let t_seq = Arc::new(t_vector(1000));
    let word_job = word_vec.clone();

    let res = task::spawn_blocking(move || -> Result<(Vec<u8>, f64, usize, u128), String> {
        // Build f' and f
        let coeffs = psi_infinity(&m_seq, &t_seq, &word_job);
        let fprime = ComplexFunction::new(
            BoundingSequence::new((*m_seq).clone()),
            ExpansionCoefficients::new(coeffs),
        );
        let f = fprime.antiderivative();

        // Create unit disk of radius 1.0, with step = 2^-acc_n_ui
        let domain = unit_disk_n(-acc_n_ui);

        // Evaluate function
        let mut img = Vec::with_capacity(domain.len());
        for &z in &domain {
            img.push(f.eval(&z));
        }
        let t0 = Instant::now();

        // build covering grid
        let bitmap = covering_grid_bitmap(&img, epsilon, delta);
        // If your function returns (bitmap, meta), split it accordingly.

        let grid_points = bitmap.width * bitmap.height; // or bitmap.width * bitmap.height if you have dims

        // >>> run EDT to estimate lambda
        let lambda_f = landau_l_via_edt_from_bitmap(&bitmap, delta);

        #[allow(clippy::unnecessary_cast)]
        let elapsed_ms = t0.elapsed().as_millis() as u128;

        // Plot to temp path
        let tmpdir = tempfile::tempdir().map_err(|e| e.to_string())?;
        let path = tmpdir.path().join("plot.png");
        plot_set(&img, path.to_str().ok_or("invalid temp path")?)
            .map_err(|e| format!("plot_set failed: {e}"))?;

        // Read file
        let png_bytes = std::fs::read(&path).map_err(|e| e.to_string())?;
        let png_magic: &[u8] = b"\x89PNG\r\n\x1a\n";
        if png_bytes.len() < 8 || &png_bytes[..8] != png_magic {
            return Err(format!("not a PNG ({} bytes)", png_bytes.len()));
        }
        Ok((png_bytes, lambda_f, grid_points, elapsed_ms))
    })
    .await;

    match res {
      Ok(Ok((png_bytes, lambda_hat, grid_points, elapsed_ms))) => {
        {
            let mut guard = state.png.lock().unwrap();
            *guard = png_bytes;
        }

        // Render
        let inner = format!(
                  r#"
        <div class="card">
        <span class="badge">Plot</span>
        <h1>Word <code>{word_display}</code>, step = 2<sup>-{acc_n_display}</sup></h1>

        <p style="margin: .25rem 0 1rem; font-size: .95rem;">
          <strong>λ<sub>f</sub> (EDT):</strong> {lambda:.6}
          <span style="opacity:.7"> · grid points: {grid_pts} · {elapsed_ms} ms</span>
        </p>

        <img src="/img" alt="plot" style="width:100%; border-radius:.5rem; border:1px solid rgba(255,255,255,.08)" />
        </div>
        "#,
        word_display = escape_html(&word),
            acc_n_display = acc_n_display,
            lambda      = lambda_hat,
            grid_pts    = grid_points,
            elapsed_ms  = elapsed_ms,
        );
        Html(page_html("Landau Explorer — Plot", &inner)).into_response()
    }
    Ok(Err(e)) => (StatusCode::INTERNAL_SERVER_ERROR, e).into_response(),
    Err(e) => (StatusCode::INTERNAL_SERVER_ERROR, e.to_string()).into_response(),
    }
}

async fn intro() -> Html<String> {
    let inner = r#"
    <div class="card">
      <span class="badge">Overview</span>
      <h1>Introduction to Landau's constant and this app</h1>
      <p>
        <strong>Landau's constant</strong> λ is the largest universal c &gt; 0 such that for every holomorphic function
        <code>f</code> on a disk <code>D<sub>r</sub>(z₀)</code> with <code>f'(z₀) ≠ 0</code>, the image <code>f(D<sub>r</sub>(z₀))</code> contains
        a Euclidean disk of radius <code>|f'(z₀)| · r · c</code>. The <em>supremum</em> of such c is λ. Known bounds are
        <code>0.5 &lt; λ ≤ 0.54325…</code>. The exact value is unknown, but it is computably approximable.
      </p>

      <h2>Algorithmic idea (Rettinger, 2012)</h2>
      <p>
        The key is to work with <em>normalised</em> functions (f(0)=0, f'(0)=1) represented via a derivative
        power series with bounded coefficients. A word over {1,2,3,4} selects nested intervals to produce
        the coefficients (your <code>psi_infinity</code>), yielding <code>f'</code>. Integrating gives <code>f</code>.
        For a fixed <code>f</code>, we estimate <code>λ<sub>f</sub> = l(f(D))</code>, the radius of the largest disk contained
        in the image <code>f(D)</code>. Taking an infimum over a compact, well-chosen class yields λ.
      </p>

      <h2>What this app does</h2>
      <ul>
        <li><strong>Build</strong> <code>f'</code> from a chosen word via <code>psi_infinity</code>, then compute the antiderivative <code>f</code>.</li>
        <li><strong>Sample</strong> a discrete domain (unit disk at dyadic mesh 2<sup>−n</sup>) and evaluate <code>f</code>.</li>
        <li><strong>Covering grid</strong>: convert the image set into a bitmap / grid for approximation.</li>
        <li><strong>EDT</strong>: run an Euclidean Distance Transform to estimate the maximal inscribed disk radius inside each image</li>
        <li><strong>Approximate</strong>: use many different words to get an accurate approximation of the constant</li>
      </ul>

      <h2>Reading colored plots of holomorphic maps</h2>
      <p>
        A common visualization is <em>domain coloring</em>: hue encodes the argument <code>arg f(z)</code>, while
        brightness/saturation encodes the magnitude (often <code>log |f(z)|</code>). Zeros appear where the hue
        wheel completes a full turn and brightness is near minimum; critical points cause characteristic
        color “pinwheels”. In our current plot we show samples of the image set; domain coloring
        swaps the point plotter for a per-pixel shader over the sampled domain. This provides us with more information that would otherwise be available with a simple plot.
      </p>

      <h2>Current problems</h2>
      <p>
        The algorithm requires very precise covering grids even in the earliest steps (epsilon requirements in Corollary 2). This means
        we very quickly reach grids with billions and billions of points, very difficult to store on computer RAM, while also keeping calculation very expensive. 
        This means algorithm cannot run on steps after 1, at least not on our current laptop computers. 
      </p>  
      
      <h2>Reference</h2>
      <p>
        R. Rettinger (2012), <em>On Computable Approximations of Landau's Constant</em>, Logical Methods
        in Computer Science, 8(4:15), 1–11. <a href="https://lmcs.episciences.org/1189/pdf" target="_blank" rel="noopener">PDF</a>
      </p>  
    </div>
  "#;

    Html(page_html("Landau Explorer — Introduction", inner))
}

async fn calc_all(
    State(state): State<SharedState>,
    Form(form): Form<CalcAllForm>,
) -> impl IntoResponse {
    let length = form.length;        // now defaults to 3 if omitted
    let step   = form.step; // guard
    let disk_decrease = -step;

    // Keep m_seq
    let m_seq = m_vec(5);

    // Run the async sweep (it re-runs the best word with plotting enabled inside)
    let approx = approximate_all_words(length, disk_decrease, m_seq).await;

    // After the sweep re-plots the best word, try to read the PNG it produced
    // If calculate_for_word_updated writes to a known "plot.png" this will pick it up
    let png_bytes = std::fs::read("test_image.png");
    match png_bytes {
        Ok(bytes) => {
            *state.png.lock().unwrap() = bytes;
            let inner = format!(
                r#"
                <div class="card">
                  <p><a href="/">← back</a></p>
                  <h1>Sweep result</h1>
                  <p><b>Length:</b> {length}</p>
                  <p><b>Step:</b> {step} (disk_decrease = {disk_decrease})</p>
                  <p><b>Best approximation:</b> {approx:.12}</p>
                  <figure>
                    <img alt="plot" src="/img" />
                    <figcaption>Best image among sampled words.</figcaption>
                  </figure>
                </div>
              "#
            );
            Html(page_html("Landau Explorer — Sweep", &inner)).into_response()
        }
        Err(e) => (
            StatusCode::INTERNAL_SERVER_ERROR,
            format!("Sweep done (best ≈ {approx:.12}), but couldn't read plot.png: {e}"),
        )
            .into_response(),
    }
}

async fn img(State(state): State<SharedState>) -> Response {
    let bytes = state.png.lock().unwrap().clone();
    if bytes.is_empty() {
        return (
            StatusCode::NOT_FOUND,
            "no image yet — POST /plot or /calc-all first",
        )
            .into_response();
    }
    Response::builder()
        .status(StatusCode::OK)
        .header(header::CONTENT_TYPE, "image/png")
        .body(bytes.into())
        .unwrap()
}

pub async fn run_server() -> anyhow::Result<()> {
    let state: SharedState = Arc::new(AppState::default());

    let app = Router::new()
        .route("/", get(home))
        .route("/intro", get(intro))
        .route("/plot", post(plot))
        .route("/calc-all", post(calc_all))
        .route("/img", get(img))
        .with_state(state);

    let addr = "127.0.0.1:3000";
    let listener = TcpListener::bind(addr).await?;
    println!("Open http://{addr}/");
    axum::serve(listener, app).await?;
    Ok(())
}
