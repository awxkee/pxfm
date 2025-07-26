use pxfm::f_asinhf;
use rayon::iter::IntoParallelIterator;
use rayon::iter::ParallelIterator;
use rug::{Assign, Float};
use std::sync::atomic::{AtomicUsize, Ordering};
use std::sync::{Arc, Mutex};
use std::thread;
use std::time::Duration;

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

    // let start_bits = 0f32.to_bits();
    // let end_bits = 1.35f32.to_bits();

    // Exhaustive: 0..=u32::MAX
    (0..=u32::MAX).into_par_iter().for_each(|bits| {
        let x = f32::from_bits(bits);

        if !x.is_finite() {
            return; // skip NaNs and infinities
        }

        let expected = Float::with_val(60, x).asinh();
        let actual = f_asinhf(x);

        executions.fetch_add(1, Ordering::Relaxed);

        let diff = count_ulp(actual, &expected);
        // if diff.is_nan() || diff.is_infinite() {
        //     return;
        // }

        if diff > 0.5 {
            failures.fetch_add(1, Ordering::Relaxed);
            exceptions.lock().unwrap().push(x);
            eprintln!(
                "Mismatch: x = {x:?}, expected = {:?}, got = {actual:?}, ULP diff = {diff}",
                expected.to_f32(),
            );
        }
    });

    // let start_bits = 20.35f64.to_bits();
    // let end_bits = 20.5f64.to_bits();

    // println!("amount {}", (end_bits - start_bits));
    //
    // // Exhaustive: 0..=u64::MAX
    // (start_bits..=end_bits).into_par_iter().for_each(|bits| {
    //     let x = f64::from_bits(bits);
    //
    //     if !x.is_finite() {
    //         return; // skip NaNs and infinities
    //     }
    //
    //     let expected = Float::with_val(70, x).y0();
    //     let actual = f_y0(x);
    //
    //     let diff = count_ulp_f64(actual, &expected);
    //
    //     if diff > 0.5 {
    //         failures.fetch_add(1, Ordering::Relaxed);
    //         exceptions.lock().unwrap().push(x);
    //         eprintln!(
    //             "Mismatch: x = {x:?}, expected = {:?}, got = {actual:?}, ULP diff = {diff}",
    //             expected.to_f64(),
    //         );
    //     }
    // });

    println!("exceptions {:?}", exceptions.lock().unwrap());

    let total_failures = failures.load(Ordering::Relaxed);
    println!(
        "Done in {:.2?}, failures: {}",
        start.elapsed(),
        total_failures
    );
    assert_eq!(total_failures, 0, "ULP failures found");
}

fn main() {
    test_f32_against_mpfr_multithreaded();
}
