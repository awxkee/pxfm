// #![feature(float_erf)]
/*
 * // Copyright 2024 (c) the Radzivon Bartoshyk. All rights reserved.
 * //
 * // Use of this source code is governed by a BSD-style
 * // license that can be found in the LICENSE file.
 */
use criterion::{Criterion, criterion_group, criterion_main};
use pxfm::{
    f_cos, f_cosf, f_cosm1, f_cospi, f_cospif, f_cot, f_cotf, f_cotpi, f_csc, f_cscf, f_secf,
    f_sin, f_sincf, f_sincos, f_sincosf, f_sincospi, f_sincospif, f_sinf, f_sinpi, f_sinpif, f_tan,
    f_tanf, f_tanpi, f_tanpif,
};
use std::hint::black_box;
use std::time::Duration;

pub fn criterion_benchmark(c: &mut Criterion) {
    let mut c = c.benchmark_group("Fast");
    c.warm_up_time(Duration::new(1, 100));
    c.sample_size(15);

    c.bench_function("system: sincpi", |b| {
        b.iter(|| {
            for i in 1..1000 {
                black_box(black_box(
                    (i as f64 / 1000.0 * std::f64::consts::PI).sin()
                        / (i as f64 / 1000.0 * std::f64::consts::PI),
                ));
            }
        })
    });

    c.bench_function("pxfm: sincpi", |b| {
        b.iter(|| {
            for i in 1..1000 {
                black_box(pxfm::f_sincpi(i as f64 / 1000.0));
            }
        })
    });

    c.bench_function("system: sincpif", |b| {
        b.iter(|| {
            for i in 1..1000 {
                black_box(black_box(
                    (i as f32 / 1000.0 * std::f32::consts::PI).sin()
                        / (i as f32 / 1000.0 * std::f32::consts::PI),
                ));
            }
        })
    });

    c.bench_function("pxfm: sincpif", |b| {
        b.iter(|| {
            for i in 1..1000 {
                black_box(pxfm::f_sincpif(i as f32 / 1000.0));
            }
        })
    });

    c.bench_function("pxfm: sincf", |b| {
        b.iter(|| {
            for i in 1..1000 {
                black_box(f_sincf(i as f32 / 1000.0));
            }
        })
    });

    c.bench_function("system: sincf", |b| {
        b.iter(|| {
            for i in 1..1000 {
                black_box(black_box((i as f32 / 1000.0).sin() / (i as f32 / 1000.0)));
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

    c.bench_function("pxfm: f_cotpi", |b| {
        b.iter(|| {
            for i in 1..1000 {
                black_box(f_cotpi(i as f64 / 1000.0));
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

    c.bench_function("libm::sin_cos", |b| {
        b.iter(|| {
            for i in 1..1000 {
                black_box(libm::sincos(i as f64 / 1000.0));
            }
        })
    });

    c.bench_function("system: sin_cos", |b| {
        b.iter(|| {
            for i in 1..1000 {
                black_box(f64::sin_cos(i as f64 / 1000.0));
            }
        })
    });

    c.bench_function("pxfm: sin_cos", |b| {
        b.iter(|| {
            for i in 1..1000 {
                black_box(f_sincos(i as f64 / 1000.0));
            }
        })
    });

    c.bench_function("pxfm: sincospi", |b| {
        b.iter(|| {
            for i in 1..1000 {
                black_box(f_sincospi(i as f64 / 1000.0));
            }
        })
    });

    c.bench_function("pxfm: f_cot", |b| {
        b.iter(|| {
            for i in 1..1000 {
                black_box(f_cot(i as f64 * 1000.0));
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

    c.bench_function("pxfm: f_sin", |b| {
        b.iter(|| {
            for i in 1..1000 {
                black_box(f_sin(i as f64 * 1000.0));
            }
        })
    });

    c.bench_function("pxfm: f_csc", |b| {
        b.iter(|| {
            for i in 1..1000 {
                black_box(f_csc(i as f64 * 1000.0));
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

    c.bench_function("pxfm: cosm1", |b| {
        b.iter(|| {
            for i in 1..1000 {
                black_box(f_cosm1(i as f64 * 1000.0));
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

    c.bench_function("pxfm: sincosf", |b| {
        b.iter(|| {
            for i in 1..1000 {
                black_box(f_sincosf(i as f32));
            }
        })
    });

    c.bench_function("pxfm: sincospif", |b| {
        b.iter(|| {
            for i in 1..1000 {
                black_box(f_sincospif(i as f32));
            }
        })
    });

    c.bench_function("pxfm: f_cotf", |b| {
        b.iter(|| {
            for i in 1..1000 {
                black_box(f_cotf(i as f32 / 10000.0 - 1.));
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

    c.bench_function("pxfm: cosf", |b| {
        b.iter(|| {
            for i in 1..1000 {
                black_box(f_cosf(i as f32));
            }
        })
    });

    c.bench_function("pxfm: f_cosm1f", |b| {
        b.iter(|| {
            for i in 1..1000 {
                black_box(pxfm::f_cosm1f(i as f32));
            }
        })
    });

    c.bench_function("pxfm: secf", |b| {
        b.iter(|| {
            for i in 1..1000 {
                black_box(f_secf(i as f32));
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

    c.bench_function("pxfm: f_sinf", |b| {
        b.iter(|| {
            for i in 1..1000 {
                black_box(f_sinf(i as f32));
            }
        })
    });

    c.bench_function("pxfm: f_cscf", |b| {
        b.iter(|| {
            for i in 1..1000 {
                black_box(f_cscf(i as f32));
            }
        })
    });
    c.finish();
}

criterion_group!(benches, criterion_benchmark);
criterion_main!(benches);
