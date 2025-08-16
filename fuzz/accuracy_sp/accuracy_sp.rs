#![no_main]

use bessel::{bessel_i0, bessel_i1};
use libfuzzer_sys::fuzz_target;
use num_complex::Complex;
use pxfm::{
    f_acosf, f_acoshf, f_acospif, f_asinf, f_asinhf, f_asinpif, f_atan2f, f_atan2pif, f_atanhf,
    f_atanpif, f_cbrtf, f_cosf, f_coshf, f_cospif, f_cotf, f_cscf, f_erfcf, f_erff, f_exp2f,
    f_exp2m1f, f_exp10f, f_exp10m1f, f_expf, f_expm1f, f_hypotf, f_i0f, f_i1f, f_j0f, f_j1f, f_k0f,
    f_k1f, f_log1pf, f_log2f, f_log2p1f, f_log10f, f_log10p1f, f_logf, f_powf, f_rcbrtf, f_rerff,
    f_rsqrtf, f_secf, f_sincf, f_sinf, f_sinhf, f_sinpif, f_tanf, f_tanhf, f_tanpif, f_y0f, f_y1f,
};
use rug::ops::Pow;
use rug::{Assign, Float};
use std::ops::{Div, Mul};
use zbessel_rs::bessel_k;

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

#[track_caller]
fn test_method(value: f32, method: fn(f32) -> f32, mpfr_value: &Float, method_name: String) {
    let xr = method(value);
    let ulp = count_ulp(xr, mpfr_value);
    assert!(
        ulp <= 0.5,
        "ULP should be less than 0.5, but it was {}, on {}, result {xr} using {method_name} on {value} and MPFR {}",
        ulp,
        value,
        mpfr_value.to_f32(),
    );
}

fn test_method_max_ulp(
    value: f32,
    method: fn(f32) -> f32,
    mpfr_value: &Float,
    method_name: String,
    max_ulp: f32,
) {
    let xr = method(value);
    let ulp = count_ulp(xr, mpfr_value);
    assert!(
        ulp <= max_ulp,
        "ULP should be less than {max_ulp}, but it was {}, on {}, result {xr} using {method_name} on {value} and MPFR {}",
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
    if !xr.is_normal() {
        return;
    }
    let ulp = count_ulp(xr, mpfr_value);
    if !ulp.is_normal() {
        return;
    }
    assert!(
        ulp <= 0.5,
        "ULP should be less than 0.5, but it was {}, using {method_name} on x: {value0}, y: {value1}, MRFR {}",
        ulp,
        mpfr_value.to_f32()
    );
}

#[track_caller]
fn test_method_2vals(
    value0: f32,
    value1: f32,
    method: fn(f32, f32) -> f32,
    mpfr_value: &Float,
    method_name: String,
) {
    let xr = method(value0, value1);
    let ulp = count_ulp(xr, mpfr_value);
    assert!(
        ulp <= 0.5,
        "ULP should be less than 0.5, but it was {}, using {method_name} on x: {value0}, y: {value1}",
        ulp
    );
}

fn compound_m1_mpfr(x: f32, y: f32) -> Float {
    let mpfr_x0 = Float::with_val(70, x);
    let mpfr_x1 = Float::with_val(70, y);
    let log1p = mpfr_x0.log2_1p().mul(&mpfr_x1);
    let exp = log1p.exp2_m1();
    exp
}

fuzz_target!(|data: (f32, f32)| {
    let x0 = data.0;
    let x1 = data.1;
    let mpfr_x0 = Float::with_val(80, x0);
    let mpfr_x1 = Float::with_val(80, x1);

    // let compound_m1_mpfr = compound_m1_mpfr(x0, x1);
    //
    // test_method_2vals_ignore_nan(
    //     x0,
    //     x1,
    //     f_compound_m1f,
    //     &compound_m1_mpfr.clone(),
    //     "f_compoundm1f".to_string(),
    // );

    // let compound_mpfr = mpfr_x0.clone().add(&Float::with_val(70, 1.)).pow(&mpfr_x1);
    //
    // test_method_2vals_ignore_nan(
    //     x0,
    //     x1,
    //     f_compoundf,
    //     &compound_mpfr.clone(),
    //     "f_compoundf".to_string(),
    // );

    if x0 < 92. && x0.abs() > 0.01 {
        test_method(x0, f_i0f, &bessel_i0(x0 as f64, 100), "f_i0f".to_string());
    }

    if x0 < 91.9 && x0.abs() > 0.01 {
        test_method(x0, f_i1f, &bessel_i1(x0 as f64, 100), "f_i1f".to_string());
    }

    if x0 < 100. && x0.is_sign_positive() && x0.abs() > 0. {
        if let Ok(expected) = bessel_k(
            Complex {
                re: x0 as f64,
                im: 0.,
            },
            0.,
            1,
            1,
        ) {
            let e = expected.values[0].re;
            test_method(x0, f_k0f, &Float::with_val(53, e), "f_k0f".to_string());
        }
    }

    if x0 < 100. && x0.is_sign_positive() && x0.abs() > 0. {
        if let Ok(expected) = bessel_k(
            Complex {
                re: x0 as f64,
                im: 0.,
            },
            1.,
            1,
            1,
        ) {
            let e = expected.values[0].re;
            test_method(x0, f_k1f, &Float::with_val(53, e), "f_k1f".to_string());
        }
    }

    test_method(
        x0,
        f_rcbrtf,
        &mpfr_x0.clone().cbrt().recip(),
        "f_rcbrtf".to_string(),
    );

    test_method(
        x0,
        f_rsqrtf,
        &mpfr_x0.clone().recip_sqrt(),
        "f_rsqrtf".to_string(),
    );
    test_method(x0, f_y1f, &mpfr_x0.clone().y1(), "f_y1f".to_string());
    test_method(x0, f_y0f, &mpfr_x0.clone().y0(), "f_y0f".to_string());
    test_method(x0, f_cscf, &mpfr_x0.clone().csc(), "f_cscf".to_string());
    test_method(x0, f_secf, &mpfr_x0.clone().sec(), "f_secf".to_string());
    test_method(x0, f_cotf, &mpfr_x0.clone().cot(), "f_cotf".to_string());
    test_method(
        x0,
        f_rerff,
        &mpfr_x0.clone().erf().recip(),
        "f_rerff".to_string(),
    );

    if x0.abs() > 0.00000000000001 {
        let mpfr_x0 = Float::with_val(25, x0);
        test_method_max_ulp(x0, f_j1f, &mpfr_x0.clone().j1(), "f_j1f".to_string(), 0.5);
    } else {
        let mpfr_x0 = Float::with_val(23, x0);
        test_method_max_ulp(x0, f_j1f, &mpfr_x0.clone().j1(), "f_j1f".to_string(), 0.);
    }

    test_method(x0, f_j0f, &mpfr_x0.clone().j0(), "f_j0f".to_string());

    test_method(x0, f_erfcf, &mpfr_x0.clone().erfc(), "f_erfcf".to_string());
    test_method(x0, f_erff, &mpfr_x0.clone().erf(), "f_erff".to_string());
    let sinc_x0 = if x0 == 0. {
        Float::with_val(70, 1.)
    } else {
        mpfr_x0.clone().sin().div(&mpfr_x0)
    };
    test_method_max_ulp(x0, f_sincf, &sinc_x0, "f_sincf".to_string(), 0.5000);
    // TODO: fix subnormals for x86 without fma
    // TODO: Fix ULP should be less than 0.5, but it was 0.71068525, using f_hypotf on x: 0.000000000000000000000000000000000000000091771, y: 0.000000000000000000000000000000000000011754585, MRFR 0.000000000000000000000000000000000000011754944
    #[cfg(any(
        all(target_arch = "x86_64", target_feature = "fma"),
        target_arch = "aarch64"
    ))]
    test_method_2vals_ignore_nan(
        x0,
        x1,
        f_hypotf,
        &mpfr_x0.clone().hypot(&mpfr_x1),
        "f_hypotf".to_string(),
    );
    test_method(
        x0,
        f_atanhf,
        &mpfr_x0.clone().atanh(),
        "f_atanhf".to_string(),
    );
    test_method(
        x0,
        f_acoshf,
        &mpfr_x0.clone().acosh(),
        "f_acoshf".to_string(),
    );
    test_method(
        x0,
        f_asinhf,
        &mpfr_x0.clone().asinh(),
        "f_asinhf".to_string(),
    );
    test_method_2vals_ignore_nan(
        x0,
        x1,
        f_atan2pif,
        &mpfr_x0.clone().atan2_pi(&mpfr_x1),
        "f_atan2pif".to_string(),
    );
    test_method_2vals_ignore_nan(
        x0,
        x1,
        f_atan2f,
        &mpfr_x0.clone().atan2(&mpfr_x1),
        "f_atan2f".to_string(),
    );
    test_method(
        x0,
        f_atanpif,
        &mpfr_x0.clone().atan_pi(),
        "f_atanpif".to_string(),
    );
    test_method(
        x0,
        f_acospif,
        &mpfr_x0.clone().acos_pi(),
        "f_acospif".to_string(),
    );
    test_method(
        x0,
        f_asinpif,
        &mpfr_x0.clone().asin_pi(),
        "f_asinpif".to_string(),
    );
    test_method(
        x0,
        f_log10p1f,
        &mpfr_x0.clone().log10_1p(),
        "f_log10p1f".to_string(),
    );
    test_method(
        x0,
        f_log2p1f,
        &mpfr_x0.clone().log2_1p(),
        "f_log2p1f".to_string(),
    );
    test_method(
        x0,
        f_log1pf,
        &mpfr_x0.clone().ln_1p(),
        "f_log1pf".to_string(),
    );
    test_method(
        x0,
        f_exp10m1f,
        &mpfr_x0.clone().exp10_m1(),
        "f_exp10m1f".to_string(),
    );
    test_method(
        x0,
        f_expm1f,
        &mpfr_x0.clone().exp_m1(),
        "f_expm1f".to_string(),
    );
    test_method(
        x0,
        f_exp2m1f,
        &mpfr_x0.clone().exp2_m1(),
        "f_exp2m1f".to_string(),
    );
    test_method(
        x0,
        f_tanpif,
        &mpfr_x0.clone().tan_pi(),
        "f_tanpif".to_string(),
    );
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
    test_method(x0, f_cbrtf, &mpfr_x0.clone().cbrt(), "f_cbrtf".to_string());
    test_method(x0, f_logf, &mpfr_x0.clone().ln(), "f_logf".to_string());
    test_method(x0, f_log2f, &mpfr_x0.clone().log2(), "f_log2f".to_string());
    test_method(
        x0,
        f_log10f,
        &mpfr_x0.clone().log10(),
        "f_log10f".to_string(),
    );
    test_method(x0, f_expf, &mpfr_x0.clone().exp(), "f_expf".to_string());
    test_method(x0, f_exp2f, &mpfr_x0.clone().exp2(), "f_exp2f".to_string());
    test_method(
        x0,
        f_exp10f,
        &mpfr_x0.clone().exp10(),
        "f_exp10f".to_string(),
    );
    test_method(x0, f_sinf, &mpfr_x0.clone().sin(), "f_sinf".to_string());
    test_method(x0, f_cosf, &mpfr_x0.clone().cos(), "f_cosf".to_string());
    test_method(x0, f_coshf, &mpfr_x0.clone().cosh(), "f_coshf".to_string());
    test_method(x0, f_sinhf, &mpfr_x0.clone().sinh(), "f_sinhf".to_string());
    test_method(x0, f_tanhf, &mpfr_x0.clone().tanh(), "f_tanhf".to_string());
    test_method_allow_not_normals(x0, f_tanf, &mpfr_x0.clone().tan(), "f_tanf".to_string());

    test_method(x0, f_acosf, &mpfr_x0.clone().acos(), "f_acosf".to_string());
    test_method(x0, f_asinf, &mpfr_x0.clone().asin(), "f_asinf".to_string());

    test_method_2vals(
        x0,
        x1,
        f_powf,
        &mpfr_x0.clone().pow(&mpfr_x1),
        "f_powf".to_string(),
    );
});
