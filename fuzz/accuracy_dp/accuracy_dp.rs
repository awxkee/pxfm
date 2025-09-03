#![no_main]
#![allow(static_mut_refs)]

use bessel::{bessel_i0, bessel_i1};
use libfuzzer_sys::fuzz_target;
use pxfm::*;
use rug::ops::Pow;
use rug::{Assign, Float};
use std::ops::{Div, Mul, Sub};

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

#[track_caller]
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
        "ULP should be less than {max_ulp}, but it was {}, on {} result {}, using {method_name} and MPFR {}",
        ulp,
        value,
        xr,
        mpfr_value.to_f64(),
    );
}

#[track_caller]
fn test_method_ulp(
    value: f64,
    method: fn(f64) -> f64,
    mpfr_value: &Float,
    method_name: String,
    max_ulp: f64,
) -> Option<f64> {
    let xr = method(value);
    let ulp = count_ulp_f64(xr, mpfr_value);
    if ulp > max_ulp {
        return Some(ulp);
    } else {
        return None;
    }
}

#[derive(Copy, Clone, Debug)]
struct Outlier {
    ulp: f64,
    on: f64,
}

fn test_method_outlier(
    value: f64,
    method: fn(f64) -> f64,
    mpfr_value: &Float,
    method_name: String,
    max_ulp: f64,
) -> Option<Outlier> {
    let xr = method(value);
    let ulp = count_ulp_f64(xr, mpfr_value);
    if ulp > max_ulp {
        Some(Outlier { ulp, on: value })
    } else {
        None
    }
}

#[track_caller]
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
        "SIN ULP should be less than {max_ulp}, but it was {}, on {}, value {}, using {method_name} and MPFR {}",
        ulp,
        value,
        xr,
        mpfr_value0.to_f64(),
    );

    let ulp = count_ulp_f64(yr, mpfr_value1);
    assert!(
        ulp <= max_ulp,
        "COS ULP should be less than {max_ulp}, but it was {}, on {}, value {}, using {method_name} and MPFR {}",
        ulp,
        value,
        yr,
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

#[track_caller]
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
        "ULP should be less than {max_ulp}, but it was {}, using {method_name} on x: {value0}, y: {value1}, value {xr}, MPFR {}",
        ulp,
        mpfr_value.to_f64()
    );
}

#[track_caller]
fn test_method_2vals_ignore_nan1(
    value0: f64,
    value1: f64,
    method: fn(f64, f64) -> f64,
    mpfr_value: &Float,
    method_name: String,
    max_ulp: f64,
) {
    let xr = method(value0, value1);
    if !xr.is_normal() {
        return;
    }
    let ulp = count_ulp_f64(xr, mpfr_value);
    if !ulp.is_normal() {
        return;
    }
    assert!(
        ulp <= max_ulp,
        "ULP should be less than {max_ulp}, but it was {}, using {method_name} on x: {value0}, y: {value1}, value {xr}, MPFR {}",
        ulp,
        mpfr_value.to_f64()
    );
}

fn log1pmx(x: f64) -> Float {
    Float::with_val(250, x)
        .ln_1p()
        .sub(&Float::with_val(250, x))
}

#[track_caller]
fn track_ulp(
    value0: f64,
    value1: f64,
    method: fn(f64, f64) -> f64,
    mpfr_value: &Float,
    method_name: String,
    max_ulp: f64,
) {
    let xr = method(value0, value1);
    if !xr.is_normal() {
        return;
    }
    let ulp = count_ulp_f64(xr, mpfr_value);
    if !ulp.is_normal() {
        return;
    }
    if ulp > 0.5 {
        let r_upper = f64::mul_add(xr, f64::from_bits(0x3c91200000000000), xr);
        let r_lower = f64::mul_add(-xr, f64::from_bits(0x3c91200000000000), xr);
        if r_upper == r_lower {
            println!(
                "RU / RL should not match {r_upper}, {r_lower} be less than {max_ulp}, but it was {ulp}, using {method_name} on x: {value0}, y: {value1}"
            );
        }
    }
    // assert!(
    //     ulp <= max_ulp,
    //     "ULP should be less than {max_ulp}, but it was {}, using {method_name} on x: {value0}, y: {value1}",
    //     ulp
    // );
}

#[inline]
pub(crate) fn is_odd_integer(x: f64) -> bool {
    let x_u = x.to_bits();
    pub(crate) const EXP_MASK: u64 = (11 << 52) - 1;
    let x_e = (x_u & EXP_MASK) >> 52;
    let lsb = (x_u | EXP_MASK).trailing_zeros();
    const E_BIAS: u64 = (1u64 << (11 - 1u64)) - 1u64;

    const UNIT_EXPONENT: u64 = E_BIAS + 52;
    x_e + lsb as u64 == UNIT_EXPONENT
}

fn compound_m1_mpfr(x: f64, y: f64) -> Float {
    let mpfr_x0 = Float::with_val(150, x);
    let mpfr_x1 = Float::with_val(150, y);
    let log1p = mpfr_x0.log2_1p().mul(&mpfr_x1);
    let ln1pf = log1p.clone().to_f64();
    // if ln1pf < -745.0 {
    //     return Float::with_val(70, -1.0);
    // } else if ln1pf > 709.0 {
    //     return Float::with_val(69, f64::INFINITY);
    // }
    let exp = log1p.exp2_m1();
    exp
}

fn powm1(x: f64, y: f64) -> Float {
    Float::with_val(200, x)
        .pow(&Float::with_val(200, y))
        .sub(&Float::with_val(200, 1))
}

fn compound_mpfr(x: f64, y: f64) -> Float {
    let mpfr_x0 = Float::with_val(150, x);
    let mpfr_x1 = Float::with_val(150, y);
    let ln1p = mpfr_x0.ln_1p().mul(&mpfr_x1);
    if x == 1. {
        return Float::with_val(70, f64::INFINITY);
    }
    let ln1pf = ln1p.clone().to_f64();
    if ln1pf < -745.0 {
        Float::with_val(70, 1.0)
    } else if ln1pf > 709.0 {
        return Float::with_val(69, f64::INFINITY);
    } else {
        return ln1p.exp();
    }
}

fn mpfr_cosm1(x: f64) -> Float {
    let r = Float::with_val(100, (x) * 0.5).sin();
    let r = r.clone().mul(&r.clone());
    r.mul(Float::with_val(100, -2))
}

static mut MAX_ULP: f64 = 0.;

fuzz_target!(|data: (f64, f64)| {
    let x0 = data.0;
    let x1 = data.1;
    let mpfr_x0 = Float::with_val(100, x0);
    let mpfr_x1 = Float::with_val(100, x1);
    let sinc_x0 = if x0 == 0. {
        Float::with_val(100, 1.)
    } else {
        mpfr_x0.clone().sin().div(&mpfr_x0)
    };
    // let compound_m1_mpfr = compound_m1_mpfr(x0, x1);
    //
    // //TODO: MPFR computes wrong values on subnormals.
    // if x0.abs() > 0.000000000000000001 {
    //     test_method_2vals_ignore_nan1(
    //         x0,
    //         x1,
    //         f_compound_m1,
    //         &compound_m1_mpfr,
    //         "f_compound_m1".to_string(),
    //         0.56,
    //     );
    // }

    if x0.abs() > 2e-6 && x1.abs() > 2e-6 && x0.abs() < 20. && x1.abs() < 20. {
        test_method_2vals_ignore_nan(x0, x1, f_powm1, &powm1(x0, x1), "f_powm1".to_string(), 1.0);
    }

    if x0.abs() > 1e-55 {
        test_method(
            x0,
            f_log1pmx,
            &log1pmx(x0),
            "f_log1pmx".to_string(),
            0.500001,
        );
    }
    if x0.abs() > 1e-8 {
        test_method(x0, f_cosm1, &mpfr_cosm1(x0), "f_cosm1".to_string(), 0.5001);
    }
    if x0 < 1e+50 {
        test_method(
            x0,
            f_digamma,
            &mpfr_x0.clone().digamma(),
            "f_digamma".to_string(),
            1.2,
        );
    }
    if x0 < 1e+50 {
        test_method(
            x0,
            f_lgamma,
            &mpfr_x0.clone().ln_abs_gamma().0,
            "f_lgamma".to_string(),
            0.5,
        );
    }
    test_method(
        x0,
        f_tgamma,
        &mpfr_x0.clone().gamma(),
        "f_tgamma".to_string(),
        0.83,
    );
    test_method(
        x0,
        f_rerf,
        &mpfr_x0.clone().erf().recip(),
        "f_rerf".to_string(),
        0.50001,
    );
    test_method(
        x0,
        f_rsqrt,
        &mpfr_x0.clone().recip_sqrt(),
        "f_rsqrt".to_string(),
        0.5,
    );

    test_method(
        x0,
        f_rcbrt,
        &mpfr_x0.clone().cbrt().recip(),
        "f_rcbrt".to_string(),
        0.50000001,
    );
    // Only search for regression MPFR takes too long
    if x0.abs() < 15. && x0.abs() > 0.01 {
        test_method(x0, f_i0, &bessel_i0(x0, 70), "f_i0".to_string(), 0.5003);
    }

    // Only search for regression MPFR takes too long
    if x0.abs() < 15. && x0.abs() > 0.01 {
        test_method(x0, f_i1, &bessel_i1(x0, 70), "f_i1".to_string(), 0.5004);
    }

    test_method(x0, f_y1, &mpfr_x0.clone().y1(), "f_y1".to_string(), 0.5003);
    test_method(x0, f_y0, &mpfr_x0.clone().y0(), "f_y0".to_string(), 0.5);
    test_method(x0, f_csc, &mpfr_x0.clone().csc(), "f_csc".to_string(), 0.5);
    test_method(x0, f_sec, &mpfr_x0.clone().sec(), "f_sec".to_string(), 0.5);
    test_method(x0, f_cot, &mpfr_x0.clone().cot(), "f_cot".to_string(), 0.5);
    test_method(x0, f_j0, &mpfr_x0.clone().j0(), "f_j0".to_string(), 0.5);
    if x0 > 0.000000000000000001 {
        test_method(x0, f_j1, &mpfr_x0.clone().j1(), "f_j1".to_string(), 0.5);
    }
    test_method(x0, f_sinc, &sinc_x0, "f_sinc".to_string(), 0.5);
    test_method(
        x0,
        f_erfc,
        &mpfr_x0.clone().erfc(),
        "f_erfc".to_string(),
        0.5,
    );

    test_method(x0, f_erf, &mpfr_x0.clone().erf(), "f_erf".to_string(), 0.5);
    test_method(
        x0,
        f_atanh,
        &mpfr_x0.clone().atanh(),
        "f_atanh".to_string(),
        0.5,
    );
    test_method(
        x0,
        f_tanh,
        &mpfr_x0.clone().tanh(),
        "f_tanh".to_string(),
        0.5,
    );
    test_method(
        x0,
        f_cosh,
        &mpfr_x0.clone().cosh(),
        "f_cosh".to_string(),
        0.500,
    );
    test_method(
        x0,
        f_sinh,
        &mpfr_x0.clone().sinh(),
        "f_sinh".to_string(),
        0.5,
    );
    test_method(
        x0,
        f_asinh,
        &mpfr_x0.clone().asinh(),
        "f_asinh".to_string(),
        0.5,
    );
    test_method(
        x0,
        f_acosh,
        &mpfr_x0.clone().acosh(),
        "f_acosh".to_string(),
        0.5,
    );
    test_method_2vals_ignore_nan(
        x0,
        x1,
        f_atan2pi,
        &mpfr_x0.clone().atan2_pi(&mpfr_x1),
        "f_atan2pi".to_string(),
        0.5,
    );
    test_method_2vals_ignore_nan(
        x0,
        x1,
        f_atan2,
        &mpfr_x0.clone().atan2(&mpfr_x1),
        "f_atan2".to_string(),
        0.5,
    );
    test_method(
        x0,
        f_log10p1,
        &mpfr_x0.clone().log10_1p(),
        "f_log10p1".to_string(),
        0.5,
    );
    test_method(
        x0,
        f_log2p1,
        &mpfr_x0.clone().log2_1p(),
        "f_log2p1".to_string(),
        0.5,
    );
    test_method(
        x0,
        f_log1p,
        &mpfr_x0.clone().ln_1p(),
        "f_log1p".to_string(),
        0.5,
    );
    test_method(
        x0,
        f_asinpi,
        &mpfr_x0.clone().asin_pi(),
        "f_asinpi".to_string(),
        0.5,
    );
    test_method(
        x0,
        f_acospi,
        &mpfr_x0.clone().acos_pi(),
        "f_acospi".to_string(),
        0.5,
    );
    test_method(
        x0,
        f_atanpi,
        &mpfr_x0.clone().atan_pi(),
        "f_atanpi".to_string(),
        0.5,
    );
    test_method(
        x0,
        f_exp10m1,
        &mpfr_x0.clone().exp10_m1(),
        "f_exp10m1".to_string(),
        0.5,
    );
    test_method(
        x0,
        f_exp2m1,
        &mpfr_x0.clone().exp2_m1(),
        "f_exp2m1".to_string(),
        0.5,
    );
    test_method(
        x0,
        f_expm1,
        &mpfr_x0.clone().exp_m1(),
        "f_expm1".to_string(),
        0.5,
    );
    test_method(
        x0,
        f_cotpi,
        &mpfr_x0.clone().tan_pi().recip(),
        "f_cotpi".to_string(),
        0.5,
    );
    test_method(
        x0,
        f_tanpi,
        &mpfr_x0.clone().tan_pi(),
        "f_tanpi".to_string(),
        0.5,
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
        0.5,
    );
    test_method(
        x0,
        f_exp10,
        &mpfr_x0.clone().exp10(),
        "f_exp10".to_string(),
        0.5,
    );
    test_method(x0, f_sin, &mpfr_x0.clone().sin(), "f_sin".to_string(), 0.5);
    test_method(x0, f_cos, &mpfr_x0.clone().cos(), "f_cos".to_string(), 0.5);
    test_method_2_outputs(
        x0,
        f_sincos,
        &mpfr_x0.clone().sin(),
        &mpfr_x0.clone().cos(),
        "f_sincos".to_string(),
        0.5,
    );
    test_method_2_outputs(
        x0,
        f_sincospi,
        &mpfr_x0.clone().sin_pi(),
        &mpfr_x0.clone().cos_pi(),
        "f_sincospi".to_string(),
        0.5,
    );
    test_method(x0, f_tan, &mpfr_x0.clone().tan(), "f_tan".to_string(), 0.5);
    test_method(
        x0,
        f_acos,
        &mpfr_x0.clone().acos(),
        "f_acos".to_string(),
        0.5,
    );
    test_method(
        x0,
        f_asin,
        &mpfr_x0.clone().asin(),
        "f_asin".to_string(),
        0.5,
    );
    test_method(
        x0,
        f_atan,
        &mpfr_x0.clone().atan(),
        "f_atan".to_string(),
        0.5,
    );
    #[cfg(any(
        all(target_arch = "x86_64", target_feature = "fma"),
        target_arch = "aarch64"
    ))]
    test_method_2vals_ignore_nan(
        x0,
        x1,
        f_hypot,
        &mpfr_x0.clone().hypot(&mpfr_x1),
        "f_hypot".to_string(),
        0.5,
    );
    test_method_2vals_ignore_nan(
        x0,
        x1,
        f_pow,
        &mpfr_x0.clone().pow(&mpfr_x1),
        "f_pow".to_string(),
        0.5,
    );
    // let compound_mpfr = compound_mpfr(x0, x1);

    // //TODO: MPFR computes wrong values on subnormals.
    // if x0 > 0.000000000000000001 {
    //     test_method_2vals_ignore_nan1(
    //         x0,
    //         x1,
    //         f_compound,
    //         &compound_mpfr,
    //         "f_compound".to_string(),
    //         0.5,
    //     );
    // }
});
