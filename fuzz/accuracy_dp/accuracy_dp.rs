#![no_main]

use libfuzzer_sys::fuzz_target;
use pxfm::{
    f_acos, f_asin, f_atan, f_cbrt, f_cos, f_cospi, f_exp, f_exp2, f_exp2m1, f_exp10, f_exp10m1,
    f_expm1, f_log, f_log2, f_log10, f_pow, f_sin, f_sincos, f_sinpi, f_tan, f_tanpi,
};
use rug::ops::Pow;
use rug::{Assign, Float};

pub fn count_ulp_f64(d: f64, c: &Float) -> f64 {
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

fn test_method(
    value: f64,
    method: fn(f64) -> f64,
    mpfr_value: &Float,
    method_name: String,
    max_ulp: f64,
) {
    let xr = method(value);
    let ulp = count_ulp_f64(xr, mpfr_value);
    assert!(
        ulp <= max_ulp,
        "ULP should be less than {max_ulp}, but it was {}, on {} result {}, using {method_name} and MPFR {:.19}",
        ulp,
        value,
        xr,
        mpfr_value.to_f64(),
    );
}

fn test_method_2_outputs(
    value: f64,
    method: fn(f64) -> (f64, f64),
    mpfr_value0: &Float,
    mpfr_value1: &Float,
    method_name: String,
    max_ulp: f64,
) {
    let (xr, yr) = method(value);
    let ulp = count_ulp_f64(xr, mpfr_value0);
    assert!(
        ulp <= max_ulp,
        "SIN ULP should be less than {max_ulp}, but it was {}, on {} using {method_name} and MPFR {}",
        ulp,
        value,
        mpfr_value0.to_f64(),
    );

    let ulp = count_ulp_f64(yr, mpfr_value1);
    assert!(
        ulp <= max_ulp,
        "COS ULP should be less than {max_ulp}, but it was {}, on {} using {method_name} and MPFR {}",
        ulp,
        value,
        mpfr_value1.to_f64(),
    );
}

fn test_method_allow_not_normals(
    value: f64,
    method: fn(f64) -> f64,
    mpfr_value: &Float,
    method_name: String,
    max_ulp: f64,
) {
    let xr = method(value);
    let ulp = count_ulp_f64(xr, mpfr_value);
    if !ulp.is_normal() {
        return;
    }
    assert!(
        ulp <= ulp,
        "ULP should be less than {max_ulp}, but it was {}, on {} using {method_name} on {value} and MPFR {}",
        xr,
        ulp,
        mpfr_value.to_f64(),
    );
}

fn test_method_2vals_ignore_nan(
    value0: f64,
    value1: f64,
    method: fn(f64, f64) -> f64,
    mpfr_value: &Float,
    method_name: String,
    max_ulp: f64,
) {
    let xr = method(value0, value1);
    let ulp = count_ulp_f64(xr, mpfr_value);
    if !ulp.is_normal() {
        return;
    }
    assert!(
        ulp <= max_ulp,
        "ULP should be less than {max_ulp}, but it was {}, using {method_name} on x: {value0}, y: {value1}",
        ulp
    );
}

fuzz_target!(|data: (f64, f64)| {
    let x0 = data.0;
    let x1 = data.0;
    let mpfr_x0 = Float::with_val(100, x0);
    let mpfr_x1 = Float::with_val(100, x1);
    test_method(
        x0,
        f_exp10m1,
        &mpfr_x0.clone().exp10_m1(),
        "f_exp10m1".to_string(),
        0.5001,
    );
    test_method(
        x0,
        f_exp2m1,
        &mpfr_x0.clone().exp2_m1(),
        "f_exp2m1".to_string(),
        0.5000,
    );
    test_method(
        x0,
        f_expm1,
        &mpfr_x0.clone().exp_m1(),
        "f_expm1".to_string(),
        0.5005,
    );
    test_method(
        x0,
        f_tanpi,
        &mpfr_x0.clone().tan_pi(),
        "f_tanpi".to_string(),
        0.5001,
    );
    test_method(
        x0,
        f_cospi,
        &mpfr_x0.clone().cos_pi(),
        "f_cospi".to_string(),
        0.5,
    );
    test_method(
        x0,
        f_sinpi,
        &mpfr_x0.clone().sin_pi(),
        "f_sinpi".to_string(),
        0.5,
    );
    test_method(
        x0,
        f_cbrt,
        &mpfr_x0.clone().cbrt(),
        "f_cbrt".to_string(),
        0.5,
    );
    test_method(x0, f_log, &mpfr_x0.clone().ln(), "f_log".to_string(), 0.5);
    test_method(
        x0,
        f_log2,
        &mpfr_x0.clone().log2(),
        "f_log2".to_string(),
        0.5,
    );
    test_method(
        x0,
        f_log10,
        &mpfr_x0.clone().log10(),
        "f_log10".to_string(),
        0.5,
    );
    test_method(x0, f_exp, &mpfr_x0.clone().exp(), "f_exp".to_string(), 0.5);
    test_method(
        x0,
        f_exp2,
        &mpfr_x0.clone().exp2(),
        "f_exp2".to_string(),
        0.5000,
    );
    test_method(
        x0,
        f_exp10,
        &mpfr_x0.clone().exp10(),
        "f_exp10".to_string(),
        0.5000,
    );
    test_method(
        x0,
        f_sin,
        &mpfr_x0.clone().sin(),
        "f_sin".to_string(),
        0.5005,
    );
    test_method(
        x0,
        f_cos,
        &mpfr_x0.clone().cos(),
        "f_cos".to_string(),
        0.5005,
    );
    test_method_2_outputs(
        x0,
        f_sincos,
        &mpfr_x0.clone().sin(),
        &mpfr_x0.clone().cos(),
        "f_sincos".to_string(),
        0.5005,
    );
    test_method(
        x0,
        f_tan,
        &mpfr_x0.clone().tan(),
        "f_tan".to_string(),
        0.50097,
    );
    test_method(
        x0,
        f_acos,
        &mpfr_x0.clone().acos(),
        "f_acos".to_string(),
        0.5009,
    );
    test_method(
        x0,
        f_asin,
        &mpfr_x0.clone().asin(),
        "f_asin".to_string(),
        0.50097,
    );
    test_method_allow_not_normals(
        x0,
        f_atan,
        &mpfr_x0.clone().atan(),
        "f_atan".to_string(),
        0.5001,
    );
    // Powf currently not really bets handles extra large argument, ULP 10000 for extra large argument
    if x0.abs() < 1e12 && x1.abs() < 1e12 {
        test_method_2vals_ignore_nan(
            x0,
            x1,
            f_pow,
            &mpfr_x0.clone().pow(&mpfr_x1),
            "f_powf".to_string(),
            0.6,
        );
    }
});
