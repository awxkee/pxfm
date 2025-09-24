use num_complex::Complex;
use pxfm::{
    f_cbrtf, f_cos, f_cospi, f_cospif, f_cotpif, f_erfcx, f_exp2m1f, f_expm1f, f_i0ef, f_i0f,
    f_i1ef, f_i1f, f_j0, f_j0f, f_j1f, f_jincpi, f_jincpif, f_k0ef, f_k0f, f_k1ef, f_k1f,
    f_lgamma_rf, f_lgammaf, f_log10p1f, f_logisticf, f_sin, f_sincpi, f_sincpif, f_sinpif, f_tanf,
    f_tanpif, f_y0f, floorf,
};
use rayon::iter::IntoParallelIterator;
use rayon::iter::ParallelIterator;
use rug::float::Constant;
use rug::ops::Pow;
use rug::{Assign, Float};
use std::ops::{Div, Mul, Sub};
use std::process::Command;
use std::sync::atomic::{AtomicUsize, Ordering};
use std::sync::{Arc, Mutex};
use std::time::Duration;
use std::{cmp, thread};
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

fn cathethusf(x: f32, y: f32) -> Float {
    Float::with_val(90, x)
        .mul(Float::with_val(90, x))
        .sub(Float::with_val(90, y).mul(Float::with_val(90, y)))
        .sqrt()
}

fn sincf(x: f32) -> Float {
    if x == 0. {
        return Float::with_val(100, 1);
    }
    Float::with_val(90, x)
        .sin_pi()
        .div(Float::with_val(90, x).mul(&Float::with_val(90, Constant::Pi)))
}

fn jincf(x: f32) -> Float {
    if x == 0. {
        return Float::with_val(100, 0.5);
    }
    Float::with_val(100, x)
        .mul(Float::with_val(100, Constant::Pi))
        .j1()
        .div(Float::with_val(100, x).mul(&Float::with_val(100, Constant::Pi)))
        .mul(&Float::with_val(100, 2))
}

fn jinc(x: f64) -> Float {
    if x == 0. {
        return Float::with_val(90, 0.5);
    }
    Float::with_val(500, x)
        .mul(Float::with_val(500, Constant::Pi))
        .j1()
        .div(Float::with_val(500, x).mul(&Float::with_val(500, Constant::Pi)))
        .mul(&Float::with_val(500, 2))
}

fn sinc(x: f64) -> Float {
    if x == 0. {
        return Float::with_val(100, 1);
    }
    Float::with_val(110, x)
        .sin_pi()
        .div(Float::with_val(110, x).mul(&Float::with_val(110, Constant::Pi)))
}

fn cathethus(x: f64, y: f64) -> Float {
    Float::with_val(90, x)
        .mul(Float::with_val(90, x))
        .sub(Float::with_val(90, y).mul(Float::with_val(90, y)))
        .sqrt()
}

#[inline(always)]
pub(crate) fn fmlaf(a: f32, b: f32, c: f32) -> f32 {
    #[cfg(any(
        all(
            any(target_arch = "x86", target_arch = "x86_64"),
            target_feature = "fma"
        ),
        target_arch = "aarch64"
    ))]
    {
        f32::mul_add(a, b, c)
    }
    #[cfg(not(any(
        all(
            any(target_arch = "x86", target_arch = "x86_64"),
            target_feature = "fma"
        ),
        target_arch = "aarch64"
    )))]
    {
        a * b + c
    }
}

fn power_1_over_2p4(x: f32) -> f32 {
    // <<FunctionApproximations`
    // ClearAll["Global`*"]
    // f[x_]:=x^(1/2.4)
    // {err,approx}=MiniMaxApproximation[f[x],{x,{0.0031308,1},5,5},WorkingPrecision->120]
    // poly=Numerator[approx][[1]];
    // coeffs=CoefficientList[poly,x];
    // TableForm[Table[Row[{"'",NumberForm[coeffs[[i+1]],{50,50},ExponentFunction->(Null&)],"',"}],{i,0,Length[coeffs]-1}]]
    // poly=Denominator[approx][[1]];
    // coeffs=CoefficientList[poly,x];
    // TableForm[Table[Row[{"'",NumberForm[coeffs[[i+1]],{50,50},ExponentFunction->(Null&)],"',"}],{i,0,Length[coeffs]-1}]]

    // {(0.0289338010459052768667829149000052753343873331610775276582557678250851411234301218465061660352860953295187058758737912162+53.8623763009617827238804917522486374010590522684609844916014721001247959077817057340461311250842991318301379035219977324 x+8461.66353644535730096787657341327383515294183720807166590641979127794213938139114937567679471650463944535822518303170833 x^2+211128.049488272220422036251872130329956832786987819720589556696955655573725068505543545903583186876442851765586421440421 x^3+900313.518840443904876317689759512896411972559937409425974431130627993863616349335036691645384583492243036482496509529713 x^4+492315.307389668369449762240479335126870840793920889759480547361599802163460950922479829770546763836543288961586364795476 x^5)/(1+564.074027620019536981196523902363185281203265251643728891604383764379752662060776624555299302071834364195384627871876281 x+40001.4751754530351018342773752031877898539860304695946652893823712363377779862830093570314799330771927269626977941502951 x^2+481502.165256419970578502904463056896789378927643242864543258820988720081861499820733694542963854146311003501446078552596 x^3+934232.440283741344491840394012151287524709719415465841614934991382211082112980461022109810287512020640854392423290894099 x^4+155984.639070076857760721626042860907787145662041895286887431732348206040358627511992445378871099253349296424713626537253 x^5),-8.40481927933470710194099105449220633672892953521701393192821980736071525155554362532804325719077695781594081505798762639*10^-6}}
    const P: [u32; 6] = [
        0x3ced0694, 0x42577313, 0x460436a7, 0x484e2e03, 0x495bcd98, 0x48f0636a,
    ];
    const Q: [u32; 6] = [
        0x3f800000, 0x440d04bd, 0x471c417a, 0x48eb1bc5, 0x49641587, 0x48185429,
    ];
    #[inline(always)]
    #[allow(clippy::too_many_arguments)]
    pub(crate) fn f_polyeval6(x: f32, a0: f32, a1: f32, a2: f32, a3: f32, a4: f32, a5: f32) -> f32 {
        let x2 = x * x;

        let u0 = fmlaf(x, a5, a4);
        let u1 = fmlaf(x, a3, a2);
        let u2 = fmlaf(x, a1, a0);

        let v0 = fmlaf(x2, u0, u1);

        fmlaf(x2, v0, u2)
    }
    let p = f_polyeval6(
        x,
        f32::from_bits(P[0]),
        f32::from_bits(P[1]),
        f32::from_bits(P[2]),
        f32::from_bits(P[3]),
        f32::from_bits(P[4]),
        f32::from_bits(P[5]),
    );
    let q = f_polyeval6(
        x,
        f32::from_bits(Q[0]),
        f32::from_bits(Q[1]),
        f32::from_bits(Q[2]),
        f32::from_bits(Q[3]),
        f32::from_bits(Q[4]),
        f32::from_bits(Q[5]),
    );
    p / q
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

    let start_bits = 700f32.to_bits();
    let end_bits = 1500f32.to_bits();
    println!("amount {}", end_bits - start_bits);

    // Exhaustive: 0..=u32::MAX
    (0..u32::MAX).into_par_iter().for_each(|bits| {
        let x = f32::from_bits(bits);

        if !x.is_finite() {
            return; // skip NaNs and infinities
        }

        // let v = match bessel_k(
        //     Complex {
        //         re: x as f64,
        //         im: 0.,
        //     },
        //     1.,
        //     2,
        //     1,
        // ) {
        //     Ok(v) => v,
        //     Err(_) => return,
        // };

        let expected_sin_pi = Float::with_val(60, x).log10_1p();
        let actual = f_log10p1f(x);
        // if actual.is_infinite() {
        //     return;
        // }

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

    // let start_bits = (700f64).to_bits();
    // let end_bits = (start_bits + 3500000);
    //
    // // Mismatch: x = 0.9999900000195318, expected = 0.6019174596052772, got = 0.6019174596052773, ULP diff = 0.5242313917684331, correct 10790, wrong 435
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
    //     let expected = jinc(x);
    //     let actual = f_jincpi(x);
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
