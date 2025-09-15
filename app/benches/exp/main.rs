// #![feature(float_erf)]
/*
 * // Copyright 2024 (c) the Radzivon Bartoshyk. All rights reserved.
 * //
 * // Use of this source code is governed by a BSD-style
 * // license that can be found in the LICENSE file.
 */
use criterion::{Criterion, criterion_group, criterion_main};
use pxfm::{
    exp, f_exp, f_exp2, f_exp2f, f_exp2m1, f_exp2m1f, f_exp10, f_exp10f, f_exp10m1, f_exp10m1f,
    f_expf, f_expm1, f_expm1f,
};
use std::hint::black_box;
use std::time::Duration;

pub fn criterion_benchmark(c: &mut Criterion) {
    let mut c = c.benchmark_group("Fast");
    c.warm_up_time(Duration::new(1, 100));
    c.sample_size(15);

    c.bench_function("pxfm: f_logisticf", |b| {
        b.iter(|| {
            for i in 1..1000 {
                black_box(pxfm::f_logisticf(i as f32 / 1000.0));
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

    c.bench_function("libm::expf", |b| {
        b.iter(|| {
            for i in 1..1000 {
                black_box(libm::expf(i as f32 / 100.0));
            }
        })
    });

    c.bench_function("system: expf", |b| {
        b.iter(|| {
            for i in 1..1000 {
                black_box(f32::exp(i as f32 / 100.0));
            }
        })
    });

    c.bench_function("pxfm: f_expf", |b| {
        b.iter(|| {
            for i in 1..1000 {
                black_box(f_expf(i as f32 / 100.0));
            }
        })
    });

    c.bench_function("libm::exp", |b| {
        b.iter(|| {
            for i in 1..1000 {
                black_box(libm::exp(i as f64 / 100.0));
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

    c.finish();
}

criterion_group!(benches, criterion_benchmark);
criterion_main!(benches);
