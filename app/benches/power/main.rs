// #![feature(float_erf)]
/*
 * // Copyright 2024 (c) the Radzivon Bartoshyk. All rights reserved.
 * //
 * // Use of this source code is governed by a BSD-style
 * // license that can be found in the LICENSE file.
 */
use criterion::{Criterion, criterion_group, criterion_main};
use pxfm::{f_compound, f_compound_m1, f_compound_m1f, f_compoundf, f_pow, f_powf, powf};
use std::hint::black_box;
use std::time::Duration;

pub fn criterion_benchmark(c: &mut Criterion) {
    let mut c = c.benchmark_group("Fast");
    c.warm_up_time(Duration::new(1, 100));
    c.sample_size(15);

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

    c.bench_function("pxfm: f_powm1f", |b| {
        b.iter(|| {
            for i in 1..1000 {
                black_box(pxfm::f_powm1f(i as f32, 0.323221324312f32 * i as f32));
            }
        })
    });

    c.bench_function("pxfm: f_compoundf", |b| {
        b.iter(|| {
            for i in 1..1000 {
                black_box(f_compoundf(i as f32, 0.323221324312f32 * i as f32));
            }
        })
    });

    c.bench_function("pxfm: f_compound_m1f", |b| {
        b.iter(|| {
            for i in 1..1000 {
                black_box(f_compound_m1f(i as f32, 0.323221324312f32 * i as f32));
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

    c.bench_function("pxfm: f_powm1", |b| {
        b.iter(|| {
            for i in 1..1000 {
                black_box(pxfm::f_powm1(i as f64, 0.0323221324312f64 * i as f64));
            }
        })
    });

    c.bench_function("pxfm: f_compound", |b| {
        b.iter(|| {
            for i in 1..1000 {
                black_box(f_compound(i as f64, 0.0323221324312f64 * i as f64));
            }
        })
    });

    c.bench_function("pxfm: f_compound_m1", |b| {
        b.iter(|| {
            for i in 1..1000 {
                black_box(f_compound_m1(i as f64, 0.0323221324312f64 * i as f64));
            }
        })
    });

    c.finish();
}

criterion_group!(benches, criterion_benchmark);
criterion_main!(benches);
