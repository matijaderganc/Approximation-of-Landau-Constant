use axum::{
    extract::{Form, State},
    http::{header, StatusCode},
    response::{Html, IntoResponse, Response},
    routing::{get, post},
    Router,
};
use serde::Deserialize;
use std::sync::{Arc, Mutex};
use tokio::{net::TcpListener, task};

use crate::covering_grids::{unit_disk_n, unit_disk_radius};
use crate::holomorphic::{BoundingSequence, ComplexFunction, ExpansionCoefficients};
use crate::plot::plot_set;
use crate::psi::{m_vec, psi_infinity, t_vector};

// NEW: bring the sweep into scope
use crate::evaluation::calculate_for_all_words_updated;

/// Shared app state: store the most recent PNG bytes (simple demo).
#[derive(Default)]
struct AppState {
    png: Mutex<Vec<u8>>,
}
type SharedState = Arc<AppState>;

#[derive(Deserialize)]
struct PlotForm {
    word: String,
    acc_n: Option<i32>, // step exponent: step = 2^-acc_n
}

// NEW: form for the sweep window
#[derive(Deserialize)]
struct CalcAllForm {
    length: usize,
    step: i32, // UI step 1 -> disk_decrease -1
}

/// Accept either "1,2,3,4" or "1234"
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
    Html(
        r#"
<!doctype html>
<html>
  <head>
    <meta charset="utf-8"/>
    <title>Word → f(D_{1})</title>
    <style>
      body{font:16px/1.4 system-ui, sans-serif; margin:2rem auto; max-width:800px}
      input[type=text], input[type=number]{width:100%; padding:.5rem; font-size:1rem}
      button{padding:.5rem 1rem; font-size:1rem}
      code{background:#f6f8fa; padding:.1rem .3rem; border-radius:.25rem}
      .row{display:grid; grid-template-columns: 1fr 180px; gap:1rem; align-items:end}
      label{display:block; margin:.5rem 0 .25rem}
      hr{margin:1.75rem 0; border:none; border-top:1px solid #eee}
    </style>
  </head>
  <body>
    <h1>Word → image of D<sub>1</sub></h1>

    <!-- Existing window -->
    <form method="post" action="/plot">
      <div class="row">
        <div>
          <label>Word (use only symbols: <code>1, 2, 3, 4</code> or <code>1234</code>):</label>
          <input type="text" name="word" placeholder="1,2,3,4" required />
        </div>
        <div>
          <label>Accuracy n (step = 2<sup>−n</sup>):</label>
          <input type="number" name="acc_n" min="6" max="18" value="10" />
        </div>
      </div>
      <p><button type="submit">Plot</button></p>
    </form>

    <hr/>

    <!-- NEW window: sweep all words of given length -->
    <form method="post" action="/calc-all">
      <div class="row">
        <div>
          <label>Word length:</label>
          <input type="number" name="length" min="1" value="3" required />
        </div>
        <div>
          <label>Step (1 → disk_decrease = −1):</label>
          <input type="number" name="step" min="1" value="1" required />
        </div>
      </div>
      <p><button type="submit">Calculate best & plot</button></p>
    </form>
  </body>
</html>
"#.to_string(),
    )
}

async fn plot(State(state): State<SharedState>, Form(form): Form<PlotForm>) -> impl IntoResponse {
    // Parse inputs
    let word_vec = parse_word(&form.word);
    let word_disp = format!("{:?}", word_vec);
    let acc_n_ui = form.acc_n.unwrap_or(5).clamp(3, 18);

    // Build inputs for the blocking job (as in your code)
    let m_seq = Arc::new(m_vec(5));
    let t_seq = Arc::new(t_vector(1000));
    let word_job = word_vec.clone();

    let res = task::spawn_blocking(move || -> Result<Vec<u8>, String> {
        // 1) Build f' and f
        let coeffs = psi_infinity(&m_seq, &t_seq, &word_job);
        let fprime = ComplexFunction::new(
            BoundingSequence::new((*m_seq).clone()),
            ExpansionCoefficients::new(coeffs),
        );
        let f = fprime.antiderivative();

        // 2) Domain: unit disk of radius 1.0, with step = 2^-acc_n_ui
        let domain = unit_disk_n(-acc_n_ui);

        // 3) Evaluate
        let mut img = Vec::with_capacity(domain.len());
        for &z in &domain {
            img.push(f.eval(&z));
        }

        // 4) Plot to temp path (let plot_set create the file)
        let tmpdir = tempfile::tempdir().map_err(|e| e.to_string())?;
        let path = tmpdir.path().join("plot.png");
        plot_set(&img, path.to_str().ok_or("invalid temp path")?)
            .map_err(|e| format!("plot_set failed: {e}"))?;

        // 5) Read bytes back
        let png_bytes = std::fs::read(&path).map_err(|e| e.to_string())?;
        let png_magic: &[u8] = b"\x89PNG\r\n\x1a\n";
        if png_bytes.len() < 8 || &png_bytes[..8] != png_magic {
            return Err(format!("not a PNG ({} bytes)", png_bytes.len()));
        }
        Ok(png_bytes)
    })
    .await;

    match res {
        Ok(Ok(png_bytes)) => {
            *state.png.lock().unwrap() = png_bytes;
            let html = format!(
                r#"
<!doctype html>
<html>
  <head>
    <meta charset="utf-8"/>
    <title>Result</title>
    <style>body{{font:16px/1.4 system-ui, sans-serif; margin:2rem auto; max-width:900px}}</style>
  </head>
  <body>
    <p><a href="/">← back</a></p>
    <p><b>Word:</b> {word_disp}</p>
    <p><b>Accuracy:</b> n = {acc_n_ui} (step = 2<sup>−{acc_n_ui}</sup>)</p>
    <img alt="plot" src="/img" />
  </body>
</html>
"#
            );
            Html(html).into_response()
        }
        Ok(Err(e)) => (StatusCode::INTERNAL_SERVER_ERROR, e).into_response(),
        Err(e) => (StatusCode::INTERNAL_SERVER_ERROR, e.to_string()).into_response(),
    }
}

// NEW: sweep handler
async fn calc_all(State(state): State<SharedState>, Form(form): Form<CalcAllForm>) -> impl IntoResponse {
    let length = form.length;
    let step = form.step.max(1); // guard
    let disk_decrease = -step;

    // Keep m_seq "as usual"
    let m_seq = m_vec(5);

    // Run the async sweep (it re-runs the best word with plotting enabled inside)
    let approx = calculate_for_all_words_updated(length, disk_decrease, m_seq).await;

    // After the sweep re-plots the best word, try to read the PNG it produced.
    // If your calculate_for_word_updated writes to a known "plot.png", this will pick it up.
    let png_bytes = std::fs::read("test_grid.png");
    match png_bytes {
        Ok(bytes) => {
            *state.png.lock().unwrap() = bytes;
            let html = format!(
                r#"
<!doctype html>
<html>
  <head>
    <meta charset="utf-8"/>
    <title>Sweep result</title>
    <style>body{{font:16px/1.4 system-ui, sans-serif; margin:2rem auto; max-width:900px}}</style>
  </head>
  <body>
    <p><a href="/">← back</a></p>
    <p><b>Length:</b> {length}</p>
    <p><b>Step:</b> {step} (disk_decrease = {disk_decrease})</p>
    <p><b>Best approximation:</b> {approx:.12}</p>
    <img alt="plot" src="/img" />
  </body>
</html>
"#
            );
            Html(html).into_response()
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
        return (StatusCode::NOT_FOUND, "no image yet — POST /plot or /calc-all first").into_response();
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
        .route("/plot", post(plot))
        // NEW route
        .route("/calc-all", post(calc_all))
        .route("/img", get(img))
        .with_state(state);

    let addr = "127.0.0.1:3000";
    let listener = TcpListener::bind(addr).await?;
    println!("Open http://{addr}/");

    axum::serve(listener, app).await?;
    Ok(())
}