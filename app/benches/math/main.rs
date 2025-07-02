/*
 * // Copyright 2024 (c) the Radzivon Bartoshyk. All rights reserved.
 * //
 * // Use of this source code is governed by a BSD-style
 * // license that can be found in the LICENSE file.
 */
use criterion::{Criterion, criterion_group, criterion_main};
use pxfm::{
    exp, f_acos, f_acosf, f_acosh, f_acoshf, f_acospi, f_acospif, f_asin, f_asinf, f_asinh,
    f_asinhf, f_asinpi, f_asinpif, f_atan, f_atan2, f_atan2f, f_atan2pi, f_atan2pif, f_atanf,
    f_atanh, f_atanhf, f_atanpi, f_atanpif, f_cbrt, f_cbrtf, f_cos, f_cosf, f_cosh, f_coshf,
    f_cospi, f_cospif, f_exp, f_exp2, f_exp2f, f_exp2m1, f_exp2m1f, f_exp10, f_exp10f, f_exp10m1,
    f_exp10m1f, f_expf, f_expm1, f_expm1f, f_j1, f_log, f_log1p, f_log1pf, f_log2, f_log2f,
    f_log2p1, f_log2p1f, f_log10, f_log10f, f_log10p1, f_log10p1f, f_logf, f_pow, f_powf, f_sin,
    f_sincos, f_sincosf, f_sinf, f_sinh, f_sinhf, f_sinpi, f_sinpif, f_tan, f_tanf, f_tanh,
    f_tanhf, f_tanpi, f_tanpif, powf,
};
use std::hint::black_box;
use std::time::Duration;

pub fn criterion_benchmark(c: &mut Criterion) {
    let mut c = c.benchmark_group("Fast");
    c.warm_up_time(Duration::new(1, 100));
    c.sample_size(15);

    c.bench_function("libm: atanh", |b| {
        b.iter(|| {
            for i in 1..1000 {
                black_box(libm::atanh(i as f64 / 1000.0));
            }
        })
    });

    c.bench_function("system: atanh", |b| {
        b.iter(|| {
            for i in 1..1000 {
                black_box(f64::atanh(i as f64 / 1000.0));
            }
        })
    });

    c.bench_function("pxfm: f_atanh", |b| {
        b.iter(|| {
            for i in 1..1000 {
                black_box(f_atanh(i as f64 / 1000.0));
            }
        })
    });

    c.bench_function("libm: tanh", |b| {
        b.iter(|| {
            for i in 1..1000 {
                black_box(libm::tanh(i as f64 / 1000.0));
            }
        })
    });

    c.bench_function("system: tanh", |b| {
        b.iter(|| {
            for i in 1..1000 {
                black_box(f64::tanh(i as f64 / 1000.0));
            }
        })
    });

    c.bench_function("pxfm: f_tanh", |b| {
        b.iter(|| {
            for i in 1..1000 {
                black_box(f_tanh(i as f64 / 1000.0));
            }
        })
    });

    c.bench_function("libm: atanhf", |b| {
        b.iter(|| {
            for i in 1..1000 {
                black_box(libm::atanhf(i as f32 / 1000.0));
            }
        })
    });

    c.bench_function("system: atanhf", |b| {
        b.iter(|| {
            for i in 1..1000 {
                black_box(f32::atanh(i as f32 / 1000.0));
            }
        })
    });

    c.bench_function("pxfm: f_atanhf", |b| {
        b.iter(|| {
            for i in 1..1000 {
                black_box(f_atanhf(i as f32 / 1000.0));
            }
        })
    });

    c.bench_function("libm: acoshf", |b| {
        b.iter(|| {
            for i in 1..1000 {
                black_box(libm::acoshf(i as f32 / 1000.0));
            }
        })
    });

    c.bench_function("system: acoshf", |b| {
        b.iter(|| {
            for i in 1..1000 {
                black_box(f32::acosh(i as f32 / 1000.0));
            }
        })
    });

    c.bench_function("pxfm: f_acoshf", |b| {
        b.iter(|| {
            for i in 1..1000 {
                black_box(f_acoshf(i as f32 / 1000.0));
            }
        })
    });

    c.bench_function("libm: asinhf", |b| {
        b.iter(|| {
            for i in 1..1000 {
                black_box(libm::asinf(i as f32 / 1000.0));
            }
        })
    });

    c.bench_function("system: asinhf", |b| {
        b.iter(|| {
            for i in 1..1000 {
                black_box(f32::asinh(i as f32 / 1000.0));
            }
        })
    });

    c.bench_function("pxfm: f_asinhf", |b| {
        b.iter(|| {
            for i in 1..1000 {
                black_box(f_asinhf(i as f32 / 1000.0));
            }
        })
    });

    c.bench_function("libm: cosh", |b| {
        b.iter(|| {
            for i in 1..1000 {
                black_box(libm::cosh(i as f64 / 1000.0));
            }
        })
    });

    c.bench_function("system: cosh", |b| {
        b.iter(|| {
            for i in 1..1000 {
                black_box(f64::cosh(i as f64 / 1000.0));
            }
        })
    });

    c.bench_function("pxfm: f_cosh", |b| {
        b.iter(|| {
            for i in 1..1000 {
                black_box(f_cosh(i as f64 / 1000.0));
            }
        })
    });

    c.bench_function("libm: sinh", |b| {
        b.iter(|| {
            for i in 1..1000 {
                black_box(libm::sinh(i as f64 / 1000.0));
            }
        })
    });

    c.bench_function("system: sinh", |b| {
        b.iter(|| {
            for i in 1..1000 {
                black_box(f64::sinh(i as f64 / 1000.0));
            }
        })
    });

    c.bench_function("pxfm: f_sinh", |b| {
        b.iter(|| {
            for i in 1..1000 {
                black_box(f_sinh(i as f64 / 1000.0));
            }
        })
    });

    c.bench_function("libm: asinh", |b| {
        b.iter(|| {
            for i in 1..1000 {
                black_box(libm::asinh(i as f64 / 1000.0));
            }
        })
    });

    c.bench_function("system: asinh", |b| {
        b.iter(|| {
            for i in 1..1000 {
                black_box(f64::asinh(i as f64 / 1000.0));
            }
        })
    });

    c.bench_function("pxfm: f_asinh", |b| {
        b.iter(|| {
            for i in 1..1000 {
                black_box(f_asinh(i as f64 / 1000.0));
            }
        })
    });

    c.bench_function("libm: acosh", |b| {
        b.iter(|| {
            for i in 1..1000 {
                black_box(libm::acosh(i as f64 / 1000.0));
            }
        })
    });

    c.bench_function("system: acosh", |b| {
        b.iter(|| {
            for i in 1..1000 {
                black_box(f64::acosh(i as f64 / 1000.0));
            }
        })
    });

    c.bench_function("pxfm: f_acosh", |b| {
        b.iter(|| {
            for i in 1..1000 {
                black_box(f_acosh(i as f64 / 1000.0));
            }
        })
    });

    c.bench_function("pxfm: f_log10p1", |b| {
        b.iter(|| {
            for i in 1..1000 {
                black_box(f_log10p1(i as f64 / 1000.0));
            }
        })
    });

    c.bench_function("pxfm: f_log2p1", |b| {
        b.iter(|| {
            for i in 1..1000 {
                black_box(f_log2p1(i as f64 / 1000.0));
            }
        })
    });

    c.bench_function("system: ln_1p", |b| {
        b.iter(|| {
            for i in 1..1000 {
                black_box(f64::ln_1p(i as f64 / 1000.0));
            }
        })
    });

    c.bench_function("pxfm: f_log1p", |b| {
        b.iter(|| {
            for i in 1..1000 {
                black_box(f_log1p(i as f64 / 1000.0));
            }
        })
    });

    c.bench_function("pxfm: f_asinpi", |b| {
        b.iter(|| {
            for i in 1..1000 {
                black_box(f_asinpi(i as f64 / 1000.0));
            }
        })
    });

    c.bench_function("pxfm: f_acospi", |b| {
        b.iter(|| {
            for i in 1..1000 {
                black_box(f_acospi(i as f64 / 1000.0));
            }
        })
    });

    c.bench_function("pxfm: f_atanpi", |b| {
        b.iter(|| {
            for i in 1..1000 {
                black_box(f_atanpi(i as f64 / 1000.0));
            }
        })
    });

    c.bench_function("pxfm: f_exp10m1", |b| {
        b.iter(|| {
            for i in 1..1000 {
                black_box(f_exp10m1(i as f64 / 1000.0));
            }
        })
    });

    c.bench_function("pxfm: f_exp2m1", |b| {
        b.iter(|| {
            for i in 1..1000 {
                black_box(f_exp2m1(i as f64 / 1000.0));
            }
        })
    });

    c.bench_function("system: exp_m1", |b| {
        b.iter(|| {
            for i in 1..1000 {
                black_box(f64::exp_m1(i as f64 / 1000.0));
            }
        })
    });

    c.bench_function("pxfm: f_expm1", |b| {
        b.iter(|| {
            for i in 1..1000 {
                black_box(f_expm1(i as f64 / 1000.0));
            }
        })
    });

    c.bench_function("pxfm: f_atan2pif", |b| {
        b.iter(|| {
            for i in 1..1000 {
                black_box(f_atan2pif(i as f32 / 1000.0, i as f32 / 1000.0 + 0.5));
            }
        })
    });

    c.bench_function("pxfm: f_atanpif", |b| {
        b.iter(|| {
            for i in 1..1000 {
                black_box(f_atanpif(i as f32 / 1000.0));
            }
        })
    });

    c.bench_function("pxfm: f_acospif", |b| {
        b.iter(|| {
            for i in 1..1000 {
                black_box(f_acospif(i as f32 / 1000.0));
            }
        })
    });

    c.bench_function("pxfm: f_asinpif", |b| {
        b.iter(|| {
            for i in 1..1000 {
                black_box(f_asinpif(i as f32 / 1000.0));
            }
        })
    });

    c.bench_function("pxfm: f_log10p1f", |b| {
        b.iter(|| {
            for i in 1..1000 {
                black_box(f_log10p1f(i as f32 / 1000.0));
            }
        })
    });

    c.bench_function("pxfm: f_log2p1f", |b| {
        b.iter(|| {
            for i in 1..1000 {
                black_box(f_log2p1f(i as f32 / 1000.0));
            }
        })
    });

    c.bench_function("system: ln_1pf", |b| {
        b.iter(|| {
            for i in 1..1000 {
                black_box(f32::ln_1p(i as f32 / 1000.0));
            }
        })
    });

    c.bench_function("pxfm: f_log1pf", |b| {
        b.iter(|| {
            for i in 1..1000 {
                black_box(f_log1pf(i as f32 / 1000.0));
            }
        })
    });

    c.bench_function("pxfm: f_exp10m1f", |b| {
        b.iter(|| {
            for i in 1..1000 {
                black_box(f_exp10m1f(i as f32 / 1000.0));
            }
        })
    });

    c.bench_function("system: expm1f", |b| {
        b.iter(|| {
            for i in 1..1000 {
                black_box(f32::exp_m1(i as f32 / 1000.0));
            }
        })
    });

    c.bench_function("pxfm: f_expm1f", |b| {
        b.iter(|| {
            for i in 1..1000 {
                black_box(f_expm1f(i as f32 / 1000.0));
            }
        })
    });

    c.bench_function("pxfm: f_exp2m1f", |b| {
        b.iter(|| {
            for i in 1..1000 {
                black_box(f_exp2m1f(i as f32 / 1000.0));
            }
        })
    });

    c.bench_function("pxfm: f_tanpi", |b| {
        b.iter(|| {
            for i in 1..1000 {
                black_box(f_tanpi(i as f64 / 1000.0));
            }
        })
    });

    c.bench_function("pxfm: f_cospi", |b| {
        b.iter(|| {
            for i in 1..1000 {
                black_box(f_cospi(i as f64 / 1000.0));
            }
        })
    });

    c.bench_function("pxfm: f_sinpi", |b| {
        b.iter(|| {
            for i in 1..1000 {
                black_box(f_sinpi(i as f64 / 1000.0));
            }
        })
    });

    c.bench_function("pxfm: f_tanpif", |b| {
        b.iter(|| {
            for i in 1..1000 {
                black_box(f_tanpif(i as f32 / 1000.0));
            }
        })
    });

    c.bench_function("pxfm: f_sinpif", |b| {
        b.iter(|| {
            for i in 1..1000 {
                black_box(f_sinpif(i as f32 / 1000.0));
            }
        })
    });

    c.bench_function("pxfm: f_cospif", |b| {
        b.iter(|| {
            for i in 1..1000 {
                black_box(f_cospif(i as f32 / 1000.0));
            }
        })
    });

    c.bench_function("libm::j1", |b| {
        b.iter(|| {
            for i in 1..1000 {
                black_box(libm::j1(i as f64 / 1000.0));
            }
        })
    });

    c.bench_function("pxfm: j1", |b| {
        b.iter(|| {
            for i in 1..1000 {
                black_box(f_j1(i as f64 / 1000.0));
            }
        })
    });

    c.bench_function("libm::atan2", |b| {
        b.iter(|| {
            for i in 1..1000 {
                black_box(libm::atan2(i as f64 / 1000.0, i as f64 / 1000.0 + 0.5));
            }
        })
    });

    c.bench_function("system: atan2", |b| {
        b.iter(|| {
            for i in 1..1000 {
                black_box(f64::atan2(i as f64 / 1000.0, i as f64 / 1000.0 + 0.5));
            }
        })
    });

    c.bench_function("pxfm: f_atan2pi", |b| {
        b.iter(|| {
            for i in 1..1000 {
                black_box(f_atan2pi(i as f64 / 1000.0, i as f64 / 1000.0 + 0.5));
            }
        })
    });

    c.bench_function("pxfm: atan2", |b| {
        b.iter(|| {
            for i in 1..1000 {
                black_box(f_atan2(i as f64 / 1000.0, i as f64 / 1000.0 + 0.5));
            }
        })
    });

    c.bench_function("libm::atan2f", |b| {
        b.iter(|| {
            for i in 1..1000 {
                black_box(libm::atan2f(i as f32 / 1000.0, i as f32 / 1000.0 + 0.5));
            }
        })
    });

    c.bench_function("system: atan2f", |b| {
        b.iter(|| {
            for i in 1..1000 {
                black_box(f32::atan2(i as f32 / 1000.0, i as f32 / 1000.0 + 0.5));
            }
        })
    });

    c.bench_function("pxfm: atan2f", |b| {
        b.iter(|| {
            for i in 1..1000 {
                black_box(f_atan2f(i as f32 / 1000.0, i as f32 / 1000.0 + 0.5));
            }
        })
    });

    c.bench_function("libm::acos", |b| {
        b.iter(|| {
            for i in 1..1000 {
                black_box(libm::acos(i as f64 / 1000.0));
            }
        })
    });

    c.bench_function("system: acos", |b| {
        b.iter(|| {
            for i in 1..1000 {
                black_box(f64::acos(i as f64 / 1000.0));
            }
        })
    });

    c.bench_function("pxfm: acos", |b| {
        b.iter(|| {
            for i in 1..1000 {
                black_box(f_acos(i as f64 / 1000.0));
            }
        })
    });

    c.bench_function("libm::sincos", |b| {
        b.iter(|| {
            for i in 1..1000 {
                black_box(libm::sincos(i as f64 * 1000.0));
            }
        })
    });

    c.bench_function("system: sincos", |b| {
        b.iter(|| {
            for i in 1..1000 {
                black_box(f64::sin_cos(i as f64 * 1000.0));
            }
        })
    });

    c.bench_function("pxfm: FMA sincos", |b| {
        b.iter(|| {
            for i in 1..1000 {
                black_box(f_sincos(i as f64 * 1000.0));
            }
        })
    });

    c.bench_function("libm::tan", |b| {
        b.iter(|| {
            for i in 1..1000 {
                black_box(libm::tan(i as f64 * 1000.0));
            }
        })
    });

    c.bench_function("system: tan", |b| {
        b.iter(|| {
            for i in 1..1000 {
                black_box(f64::tan(i as f64 * 1000.0));
            }
        })
    });

    c.bench_function("pxfm: FMA tan", |b| {
        b.iter(|| {
            for i in 1..1000 {
                black_box(f_tan(i as f64 * 1000.0));
            }
        })
    });

    c.bench_function("libm::sin", |b| {
        b.iter(|| {
            for i in 1..1000 {
                black_box(libm::sin(i as f64 * 1000.0));
            }
        })
    });

    c.bench_function("system: sin", |b| {
        b.iter(|| {
            for i in 1..1000 {
                black_box(f64::sin(i as f64 * 1000.0));
            }
        })
    });

    c.bench_function("pxfm: FMA sin", |b| {
        b.iter(|| {
            for i in 1..1000 {
                black_box(f_sin(i as f64 * 1000.0));
            }
        })
    });

    c.bench_function("libm::cos", |b| {
        b.iter(|| {
            for i in 1..1000 {
                black_box(libm::cos(i as f64 * 1000.0));
            }
        })
    });

    c.bench_function("system: cos", |b| {
        b.iter(|| {
            for i in 1..1000 {
                black_box(f64::cos(i as f64 * 1000.0));
            }
        })
    });

    c.bench_function("pxfm: FMA cos", |b| {
        b.iter(|| {
            for i in 1..1000 {
                black_box(f_cos(i as f64 * 1000.0));
            }
        })
    });

    c.bench_function("libm::sincosf", |b| {
        b.iter(|| {
            for i in 1..1000 {
                black_box(libm::sincosf(i as f32));
            }
        })
    });

    c.bench_function("system: sincosf", |b| {
        b.iter(|| {
            for i in 1..1000 {
                black_box(f32::sin_cos(i as f32));
            }
        })
    });

    c.bench_function("pxfm: FMA sincosf", |b| {
        b.iter(|| {
            for i in 1..1000 {
                black_box(f_sincosf(i as f32));
            }
        })
    });

    c.bench_function("libm::tanf", |b| {
        b.iter(|| {
            for i in 1..1000 {
                black_box(libm::tanf(i as f32 / 10000.0 - 1.));
            }
        })
    });

    c.bench_function("system::tanf", |b| {
        b.iter(|| {
            for i in 1..1000 {
                black_box(f32::tan(i as f32 / 10000.0 - 1.));
            }
        })
    });

    c.bench_function("pxfm::tanf", |b| {
        b.iter(|| {
            for i in 1..1000 {
                black_box(f_tanf(i as f32 / 10000.0 - 1.));
            }
        })
    });

    c.bench_function("libm::cbrt", |b| {
        b.iter(|| {
            for i in 1..1000 {
                black_box(libm::cbrt(i as f64));
            }
        })
    });

    c.bench_function("system: cbrt", |b| {
        b.iter(|| {
            for i in 1..1000 {
                black_box(f64::cbrt(i as f64));
            }
        })
    });

    c.bench_function("pxfm: FMA cbrt", |b| {
        b.iter(|| {
            for i in 1..1000 {
                black_box(f_cbrt(i as f64));
            }
        })
    });

    c.bench_function("libm::log10", |b| {
        b.iter(|| {
            for i in 1..1000 {
                black_box(libm::log10(i as f64));
            }
        })
    });

    c.bench_function("system: log10", |b| {
        b.iter(|| {
            for i in 1..1000 {
                black_box(f64::log10(i as f64));
            }
        })
    });

    c.bench_function("pxfm: FMA log10", |b| {
        b.iter(|| {
            for i in 1..1000 {
                black_box(f_log10(i as f64));
            }
        })
    });

    c.bench_function("libm::log10f", |b| {
        b.iter(|| {
            for i in 1..1000 {
                black_box(libm::log10f(i as f32));
            }
        })
    });

    c.bench_function("system: log10f", |b| {
        b.iter(|| {
            for i in 1..1000 {
                black_box(f32::log10(i as f32));
            }
        })
    });

    c.bench_function("pxfm: FMA log10f", |b| {
        b.iter(|| {
            for i in 1..1000 {
                black_box(f_log10f(i as f32));
            }
        })
    });

    c.bench_function("libm::exp10", |b| {
        b.iter(|| {
            for i in 1..10000 {
                black_box(libm::exp10(i as f64 / 10000.0 - 1.));
            }
        })
    });

    c.bench_function("pxfm::exp10", |b| {
        b.iter(|| {
            for i in 1..1000 {
                black_box(f_exp10(i as f64 / 10000.0 - 1.));
            }
        })
    });

    c.bench_function("libm::exp10f", |b| {
        b.iter(|| {
            for i in 1..10000 {
                black_box(libm::exp10f(i as f32 / 10000.0 - 1.));
            }
        })
    });

    c.bench_function("pxfm::exp10f", |b| {
        b.iter(|| {
            for i in 1..1000 {
                black_box(f_exp10f(i as f32 / 10000.0 - 1.));
            }
        })
    });

    c.bench_function("libm::exp2f", |b| {
        b.iter(|| {
            for i in 1..10000 {
                black_box(libm::exp2f(i as f32 / 10000.0 - 1.));
            }
        })
    });

    c.bench_function("system::exp2f", |b| {
        b.iter(|| {
            for i in 1..1000 {
                black_box(f32::exp2(i as f32 / 10000.0 - 1.));
            }
        })
    });

    c.bench_function("pxfm::exp2f", |b| {
        b.iter(|| {
            for i in 1..1000 {
                black_box(f_exp2f(i as f32 / 10000.0 - 1.));
            }
        })
    });

    c.bench_function("libm::exp2", |b| {
        b.iter(|| {
            for i in 1..1000 {
                black_box(libm::exp2(i as f64 / 10000.0 - 1.));
            }
        })
    });

    c.bench_function("system::exp2", |b| {
        b.iter(|| {
            for i in 1..1000 {
                black_box(f64::exp2(i as f64 / 10000.0 - 1.));
            }
        })
    });

    c.bench_function("pxfm::exp2", |b| {
        b.iter(|| {
            for i in 1..1000 {
                black_box(f_exp2(i as f64 / 10000.0 - 1.));
            }
        })
    });

    c.bench_function("system::exp", |b| {
        b.iter(|| {
            for i in 1..1000 {
                black_box(f64::exp(i as f64 / 10000.0 - 1.));
            }
        })
    });

    c.bench_function("pxfm::exp", |b| {
        b.iter(|| {
            for i in 1..1000 {
                black_box(f_exp(i as f64 / 10000.0 - 1.));
            }
        })
    });

    c.bench_function("libm::cbrtf", |b| {
        b.iter(|| {
            for i in 1..1000 {
                black_box(libm::cbrtf(i as f32));
            }
        })
    });

    c.bench_function("system: cbrtf", |b| {
        b.iter(|| {
            for i in 1..1000 {
                black_box(f32::cbrt(i as f32));
            }
        })
    });

    c.bench_function("pxfm: FMA cbrtf", |b| {
        b.iter(|| {
            for i in 1..1000 {
                black_box(f_cbrtf(i as f32));
            }
        })
    });

    c.bench_function("libm::cosf", |b| {
        b.iter(|| {
            for i in 1..1000 {
                black_box(libm::cosf(i as f32));
            }
        })
    });

    c.bench_function("system: cosf", |b| {
        b.iter(|| {
            for i in 1..1000 {
                black_box(f32::cos(i as f32));
            }
        })
    });

    c.bench_function("pxfm: FMA cosf", |b| {
        b.iter(|| {
            for i in 1..1000 {
                black_box(f_cosf(i as f32));
            }
        })
    });

    c.bench_function("libm::sinf", |b| {
        b.iter(|| {
            for i in 1..1000 {
                black_box(libm::sinf(i as f32));
            }
        })
    });

    c.bench_function("system: sinf", |b| {
        b.iter(|| {
            for i in 1..1000 {
                black_box(f32::sin(i as f32));
            }
        })
    });

    c.bench_function("pxfm: FMA sinf", |b| {
        b.iter(|| {
            for i in 1..1000 {
                black_box(f_sinf(i as f32));
            }
        })
    });

    c.bench_function("libm::expf", |b| {
        b.iter(|| {
            for i in 1..1000 {
                black_box(libm::expf(i as f32));
            }
        })
    });

    c.bench_function("system: expf", |b| {
        b.iter(|| {
            for i in 1..1000 {
                black_box(f32::exp(i as f32));
            }
        })
    });

    c.bench_function("pxfm: FMA expf", |b| {
        b.iter(|| {
            for i in 1..1000 {
                black_box(f_expf(i as f32));
            }
        })
    });

    c.bench_function("libm::exp", |b| {
        b.iter(|| {
            for i in 1..1000 {
                black_box(libm::exp(i as f64));
            }
        })
    });

    c.bench_function("pxfm: FMA exp", |b| {
        b.iter(|| {
            for i in 1..1000 {
                black_box(f_exp(i as f64));
            }
        })
    });

    c.bench_function("pxfm: exp", |b| {
        b.iter(|| {
            for i in 1..1000 {
                black_box(exp(i as f64));
            }
        })
    });

    c.bench_function("libm::asinf", |b| {
        b.iter(|| {
            for i in 1..1000 {
                black_box(libm::asinf(i as f32 / 1000.0));
            }
        })
    });

    c.bench_function("system::asinf", |b| {
        b.iter(|| {
            for i in 1..1000 {
                black_box(f32::asin(i as f32 / 1000.0));
            }
        })
    });

    c.bench_function("pxfm: FMA asinf", |b| {
        b.iter(|| {
            for i in 1..1000 {
                black_box(f_asinf(i as f32 / 1000.0));
            }
        })
    });

    c.bench_function("libm::acosf", |b| {
        b.iter(|| {
            for i in 1..1000 {
                black_box(libm::acosf(i as f32 / 1000.0));
            }
        })
    });

    c.bench_function("system::acosf", |b| {
        b.iter(|| {
            for i in 1..1000 {
                black_box(f32::acos(i as f32 / 1000.0));
            }
        })
    });

    c.bench_function("pxfm: FMA acosf", |b| {
        b.iter(|| {
            for i in 1..1000 {
                black_box(f_acosf(i as f32 / 1000.0));
            }
        })
    });

    c.bench_function("libm::tanhf", |b| {
        b.iter(|| {
            for i in 1..1000 {
                black_box(libm::tanhf(i as f32));
            }
        })
    });

    c.bench_function("system::tanhf", |b| {
        b.iter(|| {
            for i in 1..1000 {
                black_box(f32::tanh(i as f32));
            }
        })
    });

    c.bench_function("pxfm: FMA tanhf", |b| {
        b.iter(|| {
            for i in 1..1000 {
                black_box(f_tanhf(i as f32));
            }
        })
    });

    c.bench_function("libm::sinhf", |b| {
        b.iter(|| {
            for i in 1..1000 {
                black_box(libm::sinhf(i as f32));
            }
        })
    });

    c.bench_function("system::sinhf", |b| {
        b.iter(|| {
            for i in 1..1000 {
                black_box(f32::sinh(i as f32));
            }
        })
    });

    c.bench_function("pxfm: FMA sinhf", |b| {
        b.iter(|| {
            for i in 1..1000 {
                black_box(f_sinhf(i as f32));
            }
        })
    });

    c.bench_function("libm::coshf", |b| {
        b.iter(|| {
            for i in 1..1000 {
                black_box(libm::coshf(i as f32));
            }
        })
    });

    c.bench_function("system::coshf", |b| {
        b.iter(|| {
            for i in 1..1000 {
                black_box(f32::cosh(i as f32));
            }
        })
    });

    c.bench_function("pxfm: FMA coshf", |b| {
        b.iter(|| {
            for i in 1..1000 {
                black_box(f_coshf(i as f32));
            }
        })
    });

    c.bench_function("libm::log2f", |b| {
        b.iter(|| {
            for i in 1..1000 {
                black_box(libm::log2f(i as f32));
            }
        })
    });

    c.bench_function("system::log2f", |b| {
        b.iter(|| {
            for i in 1..1000 {
                black_box(f32::log2(i as f32));
            }
        })
    });

    c.bench_function("pxfm: FMA log2f", |b| {
        b.iter(|| {
            for i in 1..1000 {
                black_box(f_log2f(i as f32));
            }
        })
    });

    c.bench_function("libm::log2", |b| {
        b.iter(|| {
            for i in 1..1000 {
                black_box(libm::log2(i as f64));
            }
        })
    });

    c.bench_function("system::log2", |b| {
        b.iter(|| {
            for i in 1..1000 {
                black_box(f64::log2(i as f64));
            }
        })
    });

    c.bench_function("pxfm: FMA log2", |b| {
        b.iter(|| {
            for i in 1..1000 {
                black_box(f_log2(i as f64));
            }
        })
    });

    c.bench_function("libm::log", |b| {
        b.iter(|| {
            for i in 1..1000 {
                black_box(libm::log(i as f64));
            }
        })
    });

    c.bench_function("system: log", |b| {
        b.iter(|| {
            for i in 1..1000 {
                black_box(f64::ln(i as f64));
            }
        })
    });

    c.bench_function("pxfm: FMA log", |b| {
        b.iter(|| {
            for i in 1..1000 {
                black_box(f_log(i as f64));
            }
        })
    });

    c.bench_function("libm::logf", |b| {
        b.iter(|| {
            for i in 1..1000 {
                black_box(libm::logf(i as f32));
            }
        })
    });

    c.bench_function("system::logf", |b| {
        b.iter(|| {
            for i in 1..1000 {
                black_box((i as f32).ln());
            }
        })
    });

    c.bench_function("pxfm: FMA logf", |b| {
        b.iter(|| {
            for i in 1..1000 {
                black_box(f_logf(i as f32));
            }
        })
    });

    c.bench_function("libm::powf", |b| {
        b.iter(|| {
            for i in 1..1000 {
                black_box(libm::powf(i as f32, 0.323221324312f32 * i as f32));
            }
        })
    });

    c.bench_function("pxfm: powf", |b| {
        b.iter(|| {
            for i in 1..1000 {
                black_box(powf(i as f32, 0.323221324312f32 * i as f32));
            }
        })
    });

    c.bench_function("system: powf", |b| {
        b.iter(|| {
            for i in 1..1000 {
                black_box(f32::powf(i as f32, 0.323221324312f32 * i as f32));
            }
        })
    });

    c.bench_function("pxfm: FMA powf", |b| {
        b.iter(|| {
            for i in 1..1000 {
                black_box(f_powf(i as f32, 0.323221324312f32 * i as f32));
            }
        })
    });

    c.bench_function("libm::asin", |b| {
        b.iter(|| {
            for i in 1..1000 {
                black_box(libm::asin(i as f64 / 100.0));
            }
        })
    });

    c.bench_function("system: asin", |b| {
        b.iter(|| {
            for i in 1..1000 {
                black_box(f64::asin(i as f64 / 100.0));
            }
        })
    });

    c.bench_function("pxfm: FMA asin", |b| {
        b.iter(|| {
            for i in 1..1000 {
                black_box(f_asin(i as f64 / 100.0));
            }
        })
    });

    c.bench_function("libm::atan", |b| {
        b.iter(|| {
            for i in 1..1000 {
                black_box(libm::atan(i as f64 / 100.0));
            }
        })
    });

    c.bench_function("system: atan", |b| {
        b.iter(|| {
            for i in 1..1000 {
                black_box(f64::atan(i as f64 / 100.0));
            }
        })
    });

    c.bench_function("pxfm: FMA atan", |b| {
        b.iter(|| {
            for i in 1..1000 {
                black_box(f_atan(i as f64 / 100.0));
            }
        })
    });

    c.bench_function("libm::pow", |b| {
        b.iter(|| {
            for i in 1..1000 {
                black_box(libm::pow(i as f64, 0.323221324312f64 * i as f64));
            }
        })
    });

    c.bench_function("system: pow", |b| {
        b.iter(|| {
            for i in 1..1000 {
                black_box(f64::powf(i as f64, 0.323221324312f64 * i as f64));
            }
        })
    });

    c.bench_function("pxfm: FMA pow", |b| {
        b.iter(|| {
            for i in 1..1000 {
                black_box(f_pow(i as f64, 0.323221324312f64 * i as f64));
            }
        })
    });

    c.bench_function("libm::atanf", |b| {
        b.iter(|| {
            for i in 1..1000 {
                black_box(libm::atanf(i as f32));
            }
        })
    });

    c.bench_function("system: atanf", |b| {
        b.iter(|| {
            for i in 1..1000 {
                black_box(f32::atan(i as f32));
            }
        })
    });

    c.bench_function("pxfm: FMA atanf", |b| {
        b.iter(|| {
            for i in 1..1000 {
                black_box(f_atanf(i as f32));
            }
        })
    });
    c.finish();
}

criterion_group!(benches, criterion_benchmark);
criterion_main!(benches);
