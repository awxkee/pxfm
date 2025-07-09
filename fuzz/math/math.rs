#![no_main]

use libfuzzer_sys::fuzz_target;
use pxfm::{
    f_acos, f_acosf, f_acosh, f_acoshf, f_acospi, f_acospif, f_asin, f_asinf, f_asinhf, f_asinpi,
    f_asinpif, f_atan, f_atan2, f_atan2f, f_atan2pif, f_atanf, f_atanh, f_atanhf, f_atanpi,
    f_atanpif, f_cbrt, f_cbrtf, f_compound, f_compound_m1, f_compound_m1f, f_compoundf, f_cos,
    f_cosf, f_cosh, f_coshf, f_cospi, f_cospif, f_erff, f_exp, f_exp2, f_exp2f, f_exp2m1,
    f_exp2m1f, f_exp10, f_exp10f, f_exp10m1, f_exp10m1f, f_expf, f_expm1, f_expm1f, f_hypot, f_log,
    f_log1p, f_log1pf, f_log2, f_log2f, f_log2p1, f_log2p1f, f_log10, f_log10f, f_log10p1,
    f_log10p1f, f_logf, f_sin, f_sincos, f_sinf, f_sinh, f_sinhf, f_sinpi, f_sinpif, f_tanf,
    f_tanh, f_tanpi, f_tanpif,
};

fuzz_target!(|data: (f64, f32, f32, f64)| {
    let z_f32 = data.1;
    let z_f64 = data.0;
    let y_f64 = data.3;
    let y_f32 = data.2;

    _ = f_cbrtf(z_f32);
    _ = f_cbrt(z_f64);
    _ = f_atanf(z_f32);
    _ = f_cosf(z_f32);
    _ = f_exp(z_f64);
    _ = f_exp2(z_f64);
    _ = f_exp2f(z_f32);
    _ = f_exp10(z_f64);
    _ = f_exp10f(z_f32);
    _ = f_expf(z_f32);
    _ = f_log(z_f64);
    _ = f_log2(z_f64);
    _ = f_log10(z_f64);
    _ = f_logf(z_f32);
    _ = f_log2f(z_f32);
    _ = f_log10f(z_f32);
    _ = f_cosf(z_f32);
    _ = f_sinf(z_f32);
    _ = f_tanf(z_f32);
    _ = f_coshf(z_f32);
    _ = f_sinhf(z_f32);
    _ = f_acosf(z_f32);
    _ = f_asinf(z_f32);
    _ = f_sin(z_f64);
    _ = f_cos(z_f64);
    _ = f_sincos(z_f64);
    _ = f_atan(z_f64);
    _ = f_asin(z_f64);
    _ = f_acos(z_f64);
    _ = f_atan2f(z_f32, y_f32);
    _ = f_atan2pif(z_f32, y_f32);
    _ = f_atan2(z_f64, y_f64);
    _ = f_sinpif(z_f32);
    _ = f_exp10m1f(z_f32);
    _ = f_exp2m1f(z_f32);
    _ = f_expm1f(z_f32);
    _ = f_tanpif(z_f32);
    _ = f_sinpif(z_f32);
    _ = f_cospif(z_f32);
    _ = f_sinpi(z_f64);
    _ = f_cospi(z_f64);
    _ = f_tanpi(z_f64);
    _ = f_log1pf(z_f32);
    _ = f_log2p1f(z_f32);
    _ = f_log10p1f(z_f32);
    _ = f_asinpif(z_f32);
    _ = f_acospif(z_f32);
    _ = f_atanpif(z_f32);
    _ = f_expm1(z_f64);
    _ = f_exp2m1(z_f64);
    _ = f_exp10m1(z_f64);
    _ = f_atanpi(z_f64);
    _ = f_acospi(z_f64);
    _ = f_asinpi(z_f64);
    _ = f_log10p1(z_f64);
    _ = f_log2p1(z_f64);
    _ = f_log1p(z_f64);
    _ = f_cosh(z_f64);
    _ = f_tanh(z_f64);
    _ = f_sinh(z_f64);
    _ = f_acosh(z_f64);
    _ = f_acosh(z_f64);
    _ = f_asinhf(z_f32);
    _ = f_compoundf(z_f32, y_f32);
    _ = f_compound_m1f(z_f32, y_f32);
    _ = f_compound(z_f64, y_f64);
    _ = f_compound_m1(z_f64, y_f64);
    _ = f_acoshf(z_f32);
    _ = f_atanhf(z_f32);
    _ = f_atanh(z_f64);
    _ = f_tanh(z_f64);
    _ = f_hypot(z_f64, z_f64);
    _ = f_erff(z_f32);
});
