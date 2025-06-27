#![no_main]

use libfuzzer_sys::fuzz_target;
use pxfm::{
    f_acosf, f_asinf, f_cbrtf, f_cosf, f_coshf, f_cospif, f_exp2f, f_exp10f, f_expf, f_log2f,
    f_log10f, f_logf, f_powf, f_sinf, f_sinhf, f_sinpif, f_tanf, f_tanhf,
};
use rug::ops::Pow;
use rug::{Assign, Float};

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

fn test_method(value: f32, method: fn(f32) -> f32, mpfr_value: &Float, method_name: String) {
    let xr = method(value);
    let ulp = count_ulp(xr, mpfr_value);
    assert!(
        ulp <= 0.5,
        "ULP should be less than 0.5, but it was {}, on {} using {method_name} on {value} and MPFR {}",
        ulp,
        value,
        mpfr_value.to_f32(),
    );
}

fn test_method_allow_not_normals(
    value: f32,
    method: fn(f32) -> f32,
    mpfr_value: &Float,
    method_name: String,
) {
    let xr = method(value);
    let ulp = count_ulp(xr, mpfr_value);
    if !ulp.is_normal() {
        return;
    }
    assert!(
        ulp <= 0.5,
        "ULP should be less than 0.5, but it was {}, on {} using {method_name} on {value} and MPFR {}",
        value,
        ulp,
        mpfr_value.to_f32(),
    );
}

fn test_method_2vals_ignore_nan(
    value0: f32,
    value1: f32,
    method: fn(f32, f32) -> f32,
    mpfr_value: &Float,
    method_name: String,
) {
    let xr = method(value0, value1);
    let ulp = count_ulp(xr, mpfr_value);
    if !ulp.is_normal() {
        return;
    }
    assert!(
        ulp <= 0.5,
        "ULP should be less than 0.5, but it was {}, using {method_name} on x: {value0}, y: {value1}",
        ulp
    );
}

fuzz_target!(|data: (f32, f32)| {
    let x0 = data.0;
    let x1 = data.0;
    let mpfr_x0 = Float::with_val(70, x0);
    let mpfr_x1 = Float::with_val(70, x1);
    test_method(
        x0,
        f_cospif,
        &mpfr_x0.clone().cos_pi(),
        "f_cospif".to_string(),
    );
    test_method(
        x0,
        f_sinpif,
        &mpfr_x0.clone().sin_pi(),
        "f_sinpif".to_string(),
    );
    // test_method(x0, f_cbrtf, &mpfr_x0.clone().cbrt(), "f_cbrtf".to_string());
    // test_method(x0, f_logf, &mpfr_x0.clone().ln(), "f_logf".to_string());
    // test_method(x0, f_log2f, &mpfr_x0.clone().log2(), "f_log2f".to_string());
    // test_method(
    //     x0,
    //     f_log10f,
    //     &mpfr_x0.clone().log10(),
    //     "f_log10f".to_string(),
    // );
    // test_method(x0, f_expf, &mpfr_x0.clone().exp(), "f_expf".to_string());
    // test_method(x0, f_exp2f, &mpfr_x0.clone().exp2(), "f_exp2f".to_string());
    // test_method(
    //     x0,
    //     f_exp10f,
    //     &mpfr_x0.clone().exp10(),
    //     "f_exp10f".to_string(),
    // );
    // test_method(x0, f_sinf, &mpfr_x0.clone().sin(), "f_sinf".to_string());
    // test_method(x0, f_cosf, &mpfr_x0.clone().cos(), "f_cosf".to_string());
    // test_method(x0, f_coshf, &mpfr_x0.clone().cosh(), "f_coshf".to_string());
    // test_method(x0, f_sinhf, &mpfr_x0.clone().sinh(), "f_sinhf".to_string());
    // test_method(x0, f_tanhf, &mpfr_x0.clone().tanh(), "f_tanhf".to_string());
    // test_method_allow_not_normals(x0, f_tanf, &mpfr_x0.clone().tan(), "f_tanf".to_string());
    //
    // test_method(x0, f_acosf, &mpfr_x0.clone().acos(), "f_acosf".to_string());
    // test_method(x0, f_asinf, &mpfr_x0.clone().asin(), "f_asinf".to_string());
    //
    // test_method_2vals_ignore_nan(
    //     x0,
    //     x1,
    //     f_powf,
    //     &mpfr_x0.clone().pow(&mpfr_x1),
    //     "f_powf".to_string(),
    // );
});
