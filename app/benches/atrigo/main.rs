// #![feature(float_erf)]
/*
 * // Copyright 2024 (c) the Radzivon Bartoshyk. All rights reserved.
 * //
 * // Use of this source code is governed by a BSD-style
 * // license that can be found in the LICENSE file.
 */
use criterion::{Criterion, criterion_group, criterion_main};
use pxfm::{
    f_acos, f_acosf, f_acospi, f_acospif, f_asin, f_asinf, f_asinpi, f_asinpif, f_atan, f_atan2,
    f_atan2f, f_atan2pi, f_atan2pif, f_atanf, f_atanpi, f_atanpif, f_cos, f_cosf, f_cosm1, f_cospi,
    f_cospif, f_cot, f_cotf, f_cotpi, f_csc, f_cscf, f_secf, f_sin, f_sincf, f_sincos, f_sincosf,
    f_sincospi, f_sincospif, f_sinf, f_sinpi, f_sinpif, f_tan, f_tanf, f_tanpi, f_tanpif,
};
use std::hint::black_box;
use std::time::Duration;

pub fn criterion_benchmark(c: &mut Criterion) {
    let mut c = c.benchmark_group("Fast");
    c.warm_up_time(Duration::new(1, 100));
    c.sample_size(15);

    c.bench_function("pxfm: f_asinpi", |b| {
        b.iter(|| {
            for i in 1..1000 {
                black_box(f_asinpi(black_box(i as f64 / 1000.0)));
            }
        })
    });

    c.bench_function("pxfm: f_acospi", |b| {
        b.iter(|| {
            for i in 1..1000 {
                black_box(f_acospi(black_box(i as f64 / 1000.0)));
            }
        })
    });

    c.bench_function("pxfm: f_atanpi", |b| {
        b.iter(|| {
            for i in 1..1000 {
                black_box(f_atanpi(black_box(i as f64 / 1000.0)));
            }
        })
    });

    c.bench_function("pxfm: f_atan2pif", |b| {
        b.iter(|| {
            for i in 1..1000 {
                black_box(f_atan2pif(
                    black_box(i as f32 / 1000.0),
                    black_box(i as f32 / 1000.0 + 0.5),
                ));
            }
        })
    });

    c.bench_function("pxfm: f_atanpif", |b| {
        b.iter(|| {
            for i in 1..1000 {
                black_box(f_atanpif(black_box(i as f32 / 1000.0)));
            }
        })
    });

    c.bench_function("pxfm: f_acospif", |b| {
        b.iter(|| {
            for i in 1..1000 {
                black_box(f_acospif(black_box(i as f32 / 1000.0)));
            }
        })
    });

    c.bench_function("pxfm: f_asinpif", |b| {
        b.iter(|| {
            for i in 1..1000 {
                black_box(f_asinpif(black_box(i as f32 / 1000.0)));
            }
        })
    });

    c.bench_function("libm::atan2", |b| {
        b.iter(|| {
            for i in 1..1000 {
                black_box(libm::atan2(
                    black_box(i as f64 / 1000.0),
                    black_box(i as f64 / 1000.0 + 0.5),
                ));
            }
        })
    });

    c.bench_function("system: atan2", |b| {
        b.iter(|| {
            for i in 1..1000 {
                black_box(f64::atan2(
                    black_box(i as f64 / 1000.0),
                    black_box(i as f64 / 1000.0 + 0.5),
                ));
            }
        })
    });

    c.bench_function("pxfm: f_atan2pi", |b| {
        b.iter(|| {
            for i in 1..1000 {
                black_box(f_atan2pi(
                    black_box(i as f64 / 1000.0),
                    black_box(i as f64 / 1000.0 + 0.5),
                ));
            }
        })
    });

    c.bench_function("pxfm: atan2", |b| {
        b.iter(|| {
            for i in 1..1000 {
                black_box(f_atan2(
                    black_box(i as f64 / 1000.0),
                    black_box(i as f64 / 1000.0 + 0.5),
                ));
            }
        })
    });

    c.bench_function("libm::atan2f", |b| {
        b.iter(|| {
            for i in 1..1000 {
                black_box(libm::atan2f(
                    black_box(i as f32 / 1000.0),
                    black_box(i as f32 / 1000.0 + 0.5),
                ));
            }
        })
    });

    c.bench_function("system: atan2f", |b| {
        b.iter(|| {
            for i in 1..1000 {
                black_box(f32::atan2(
                    black_box(i as f32 / 1000.0),
                    black_box(i as f32 / 1000.0 + 0.5),
                ));
            }
        })
    });

    c.bench_function("pxfm: atan2f", |b| {
        b.iter(|| {
            for i in 1..1000 {
                black_box(f_atan2f(
                    black_box(i as f32 / 1000.0),
                    black_box(i as f32 / 1000.0 + 0.5),
                ));
            }
        })
    });

    c.bench_function("libm::acos", |b| {
        b.iter(|| {
            for i in 1..1000 {
                black_box(libm::acos(black_box(i as f64 / 1000.0)));
            }
        })
    });

    c.bench_function("system: acos", |b| {
        b.iter(|| {
            for i in 1..1000 {
                black_box(f64::acos(black_box(i as f64 / 1000.0)));
            }
        })
    });

    c.bench_function("pxfm: acos", |b| {
        b.iter(|| {
            for i in 1..1000 {
                black_box(f_acos(black_box(i as f64 / 1000.0)));
            }
        })
    });

    c.bench_function("libm::asinf", |b| {
        b.iter(|| {
            for i in 1..1000 {
                black_box(libm::asinf(black_box(i as f32 / 1000.0)));
            }
        })
    });

    c.bench_function("system::asinf", |b| {
        b.iter(|| {
            for i in 1..1000 {
                black_box(f32::asin(black_box(i as f32 / 1000.0)));
            }
        })
    });

    c.bench_function("pxfm: FMA asinf", |b| {
        b.iter(|| {
            for i in 1..1000 {
                black_box(f_asinf(black_box(i as f32 / 1000.0)));
            }
        })
    });

    c.bench_function("libm::acosf", |b| {
        b.iter(|| {
            for i in 1..1000 {
                black_box(libm::acosf(black_box(i as f32 / 1000.0)));
            }
        })
    });

    c.bench_function("system::acosf", |b| {
        b.iter(|| {
            for i in 1..1000 {
                black_box(f32::acos(black_box(i as f32 / 1000.0)));
            }
        })
    });

    c.bench_function("pxfm: FMA acosf", |b| {
        b.iter(|| {
            for i in 1..1000 {
                black_box(f_acosf(black_box(i as f32 / 1000.0)));
            }
        })
    });

    c.bench_function("libm::asin", |b| {
        b.iter(|| {
            for i in 1..1000 {
                black_box(libm::asin(black_box(i as f64 / 100.0)));
            }
        })
    });

    c.bench_function("system: asin", |b| {
        b.iter(|| {
            for i in 1..1000 {
                black_box(f64::asin(black_box(i as f64 / 100.0)));
            }
        })
    });

    c.bench_function("pxfm: FMA asin", |b| {
        b.iter(|| {
            for i in 1..1000 {
                black_box(f_asin(black_box(i as f64 / 100.0)));
            }
        })
    });

    c.bench_function("libm::atan", |b| {
        b.iter(|| {
            for i in 1..1000 {
                black_box(libm::atan(black_box(i as f64 / 100.0)));
            }
        })
    });

    c.bench_function("system: atan", |b| {
        b.iter(|| {
            for i in 1..1000 {
                black_box(f64::atan(black_box(i as f64 / 100.0)));
            }
        })
    });

    c.bench_function("pxfm: FMA atan", |b| {
        b.iter(|| {
            for i in 1..1000 {
                black_box(f_atan(black_box(i as f64 / 100.0)));
            }
        })
    });

    c.bench_function("libm::atanf", |b| {
        b.iter(|| {
            for i in 1..1000 {
                black_box(libm::atanf(black_box(i as f32)));
            }
        })
    });

    c.bench_function("system: atanf", |b| {
        b.iter(|| {
            for i in 1..1000 {
                black_box(f32::atan(black_box(i as f32)));
            }
        })
    });

    c.bench_function("pxfm: FMA atanf", |b| {
        b.iter(|| {
            for i in 1..1000 {
                black_box(f_atanf(black_box(i as f32)));
            }
        })
    });

    c.finish();
}

criterion_group!(benches, criterion_benchmark);
criterion_main!(benches);
