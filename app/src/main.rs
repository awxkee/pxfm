use bessel::{bessel_i0, bessel_i1, bessel_k0, bessel_k1};
use num_complex::Complex;
use pxfm::{
    f_cathetusf, f_cosm1, f_cosm1f, f_cotf, f_cotpi, f_cotpif, f_digamma, f_digammaf, f_erfcinv,
    f_erfcinvf, f_erfcx, f_erfcxf, f_erfinv, f_erfinvf, f_expf, f_expm1, f_i0, f_i0f, f_i1, f_i1f,
    f_i2, f_i2f, f_j0, f_j1, f_k0, f_k0f, f_k1, f_k1f, f_k2f, f_lgammaf, f_log1pf, f_log1pmx,
    f_log1pmxf, f_sin, f_sinmx, f_sinmxf, f_tanf, f_tanpi, f_tanpif, f_trigamma, f_trigammaf, f_y0,
    f_y0f, f_y1, f_y1f,
};
use rand::Rng;
use rayon::iter::IntoParallelIterator;
use rayon::iter::ParallelIterator;
use rug::{Assign, Float};
use std::ops::{Add, Mul, Sub};
use std::process::Command;
use std::sync::atomic::{AtomicUsize, Ordering};
use std::sync::{Arc, Mutex};
use std::thread;
use std::time::Duration;
use zbessel_rs::{bessel_i, bessel_k};

fn compute_besselk(x: f64) -> Result<Float, Box<dyn std::error::Error>> {
    let r = x.to_string();
    let output = Command::new("python3")
        .arg("bessel/besselk.py")
        .arg(r)
        .output()?;

    if !output.status.success() {
        let stdout = str::from_utf8(&output.stderr)?;
        return Err(format!("Error: {:?}", stdout).into());
    }

    let stdout = str::from_utf8(&output.stdout)?;
    let result = Float::parse(stdout.trim()).unwrap();
    Ok(Float::with_val(100, result))
}

fn count_ulp(d: f32, c: &Float) -> f32 {
    let c2 = c.to_f32();

    if (c2 == 0. || c2.is_subnormal()) && (d == 0. || d.is_subnormal()) {
        return 0.;
    }

    if (c2 == 0.) && (d != 0.) {
        return 10000.;
    }

    if c2.is_infinite() && d.is_infinite() {
        return 0.;
    }

    if d.is_nan() && c.is_nan() {
        return 0.;
    }

    let prec = c.prec();

    let mut fry = Float::with_val(prec, d);

    let mut frw = Float::new(prec);

    let (_, e) = c.to_f32_exp();

    frw.assign(Float::u_exp(1, e - 24_i32));

    fry -= c;
    fry /= &frw;
    fry.to_f32().abs()
}

fn count_ulp_f64(d: f64, c: &Float) -> f64 {
    let c2 = c.to_f64();

    if (c2 == 0. || c2.is_subnormal()) && (d == 0. || d.is_subnormal()) {
        return 0.;
    }

    if (c2 == 0.) && (d != 0.) {
        return 10000.;
    }

    if c2.is_infinite() && d.is_infinite() {
        return 0.;
    }

    if d.is_nan() && c.is_nan() {
        return 0.;
    }

    let prec = c.prec();

    let mut fry = Float::with_val(prec, d);

    let mut frw = Float::new(prec);

    let (_, e) = c.to_f64_exp();

    frw.assign(Float::u_exp(1, e - 53_i32));

    fry -= c;
    fry /= &frw;
    fry.to_f64().abs()
}

fn mpfr_cosm1f(x: f32) -> Float {
    let r = Float::with_val(100, (x as f64) * 0.5).sin();
    let r = r.clone().mul(&r.clone());
    r.mul(Float::with_val(100, -2))
}

fn log1pmxf(x: f32) -> Float {
    Float::with_val(90, x).ln_1p().sub(&Float::with_val(90, x))
}

fn sinmxf(x: f32) -> Float {
    Float::with_val(150, x).sin().sub(&Float::with_val(150, x))
}

fn erfcx(x: f64) -> Float {
    let dx2 = Float::with_val(150, x).mul(&Float::with_val(150, x));
    Float::with_val(150, x).erfc().mul(&dx2.exp())
}

fn sinmx(x: f64) -> Float {
    Float::with_val(250, x).sin().sub(&Float::with_val(250, x))
}

fn log1pmx(x: f64) -> Float {
    Float::with_val(150, x)
        .ln_1p()
        .sub(&Float::with_val(150, x))
}

fn mpfr_cosm1(x: f64) -> Float {
    let r = Float::with_val(100, (x as f64) * 0.5).sin();
    let r = r.clone().mul(&r.clone());
    r.mul(Float::with_val(100, -2))
}

#[allow(static_mut_refs)]

fn test_f32_against_mpfr_multithreaded() {
    use std::time::Instant;

    let start = Instant::now();
    let executions = Arc::new(AtomicUsize::new(0));
    let failures = Arc::new(AtomicUsize::new(0));

    let exec1 = executions.clone();
    let fail1 = failures.clone();

    thread::spawn(move || {
        loop {
            unsafe {
                thread::sleep(Duration::from_secs(25));
                let elapsed = start.elapsed();
                eprintln!(
                    "[{:?}] Failures so far: {} Executions so far: {}, Percentage {}",
                    elapsed,
                    fail1.load(Ordering::Relaxed),
                    exec1.load(Ordering::Relaxed),
                    fail1.load(Ordering::Relaxed) as f32 / exec1.load(Ordering::Relaxed) as f32,
                );
            }
        }
    });
    let mut exceptions = Arc::new(Mutex::new(Vec::<f32>::new()));

    let start_bits = (3f32).to_bits();
    let end_bits = (150f32).to_bits();
    println!("amount {}", end_bits - start_bits);

    // Exhaustive: 0..=u32::MAX
    (start_bits..end_bits).into_par_iter().for_each(|bits| {
        let x = f32::from_bits(bits);

        if !x.is_finite() {
            return; // skip NaNs and infinities
        }

        // let v = match bessel_k(
        //     Complex {
        //         re: x as f64,
        //         im: 0.,
        //     },
        //     2.,
        //     1,
        //     1,
        // ) {
        //     Ok(v) => v,
        //     Err(_) => return,
        // };

        let y = rand::rng().random_range(0f32..1e10f32);

        let expected_sin_pi = Float::with_val(90, x)
            .mul(Float::with_val(90, x))
            .sub(Float::with_val(90, y).mul(Float::with_val(90, y)))
            .sqrt();
        let actual = f_cathetusf(x, y);
        if actual.is_infinite() {
            return;
        }

        executions.fetch_add(1, Ordering::Relaxed);

        let diff = count_ulp(actual, &Float::with_val(90, expected_sin_pi.clone()));
        // if diff.is_nan() || diff.is_infinite() {
        //     return;
        // }

        if diff > 0.5 {
            failures.fetch_add(1, Ordering::Relaxed);
            exceptions.lock().unwrap().push(x);
            eprintln!(
                "Mismatch: x = {x:?}, expected = {:?}, got = {actual:?}, ULP diff = {diff}",
                expected_sin_pi.to_f32(),
            );
        }
    });

    let start_bits = (-8f64).to_bits();
    let end_bits = (start_bits + 350000);

    // Mismatch: x = 0.9999900000195318, expected = 0.6019174596052772, got = 0.6019174596052773, ULP diff = 0.5242313917684331, correct 10790, wrong 435

    //
    // // // Exhaustive: 0..=u64::MAX
    // (start_bits..=end_bits).into_par_iter().for_each(|bits| {
    //     let x = f64::from_bits(bits);
    //
    //     if !x.is_finite() {
    //         return; // skip NaNs and infinities
    //     }
    //
    //     // let v = match bessel_k(
    //     //     Complex {
    //     //         re: x,
    //     //         im: 0.,
    //     //     },
    //     //     0.,
    //     //     1,
    //     //     1,
    //     // ) {
    //     //     Ok(v) => v,
    //     //     Err(_) => return,
    //     // };
    //
    //     let expected = erfcx(x); //Float::with_val(90, (x as f64).trigamma());
    //     let actual = f_erfcx(x);
    //
    //     let diff = count_ulp_f64(actual, &expected);
    //
    //     let execs = executions.fetch_add(1, Ordering::Relaxed);
    //
    //     if diff > 0.5 {
    //         let f = failures.fetch_add(1, Ordering::Relaxed);
    //         exceptions.lock().unwrap().push(x);
    //         eprintln!(
    //             "Mismatch: x = {x:?}, expected = {:?}, got = {actual:?}, ULP diff = {diff}, correct {}, wrong {}",
    //             expected.to_f64(),
    //             execs - f,
    //             f,
    //         );
    //     }
    // });

    let ex = exceptions.lock().unwrap();
    println!("exceptions count {}: {:?}", ex.len(), ex);

    let total_failures = failures.load(Ordering::Relaxed);
    println!(
        "Done in {:.2?}, total {}, failures: {}",
        start.elapsed(),
        end_bits - start_bits,
        total_failures
    );
    assert_eq!(total_failures, 0, "ULP failures found");
}

fn find_cutoff() {
    let mut scratch = -0.1;
    let mut value = -9f64;
    let mut depth = 0;
    loop {
        let rs = f_erfcx(value);
        if rs.is_infinite() || rs.is_nan() {
            // if rs == 0. {
            println!(
                "found basic cutoff between {}, next {}",
                value - scratch,
                value
            );
            value -= scratch;
            if depth >= 9 {
                // time to refine
                loop {
                    value = f64::from_bits(value.to_bits() + 1);
                    let rs = f_erfcx(value);
                    // if rs == 0. {
                    if rs.is_infinite() || rs.is_nan() {
                        panic!(
                            "found basic cutoff between {}, 0x{:16x}u64, next {}",
                            value - scratch,
                            value.to_bits(),
                            value
                        );
                    }
                }
                break;
            }
            depth += 1;
            scratch /= 10.0;
        }
        value += scratch;
    }
}

fn main() {
    // find_cutoff();
    test_f32_against_mpfr_multithreaded();
}
