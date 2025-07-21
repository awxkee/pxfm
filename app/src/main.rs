use pxfm::{f_j1, f_j1f};
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

#[allow(static_mut_refs)]

fn test_f32_against_mpfr_multithreaded() {
    use std::time::Instant;

    let start = Instant::now();
    let failures = AtomicUsize::new(0);

    // thread::spawn(move || {
    //     loop {
    //         unsafe {
    //             thread::sleep(Duration::from_secs(25));
    //             let elapsed = start.elapsed();
    //             eprintln!(
    //                 "[{:?}] Failures so far: {} Success so far: {}",
    //                 elapsed,
    //                 pxfm::OVERSHOOT.load(Ordering::Relaxed),
    //                 pxfm::SUCCESS.load(Ordering::Relaxed),
    //             );
    //         }
    //     }
    // });

    let mut exceptions = Arc::new(Mutex::new(Vec::<f32>::new()));

    // Exhaustive: 0..=u32::MAX
    (0u32..=u32::MAX).into_par_iter().for_each(|bits| {
        let x = f32::from_bits(bits);

        if !x.is_finite() {
            return; // skip NaNs and infinities
        }

        let expected = Float::with_val(60, x).j1();
        let actual = f_j1f(x);

        let diff = count_ulp(actual, &expected);

        if diff > 0.5 {
            failures.fetch_add(1, Ordering::Relaxed);
            exceptions.lock().unwrap().push(x);
            eprintln!(
                "Mismatch: x = {x:?}, expected = {:?}, got = {actual:?}, ULP diff = {diff}",
                expected.to_f32(),
            );
        }
    });

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
