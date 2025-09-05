// #![feature(float_erf)]
/*
 * // Copyright 2024 (c) the Radzivon Bartoshyk. All rights reserved.
 * //
 * // Use of this source code is governed by a BSD-style
 * // license that can be found in the LICENSE file.
 */
use criterion::{Criterion, criterion_group, criterion_main};
use std::hint::black_box;
use std::time::Duration;

pub fn criterion_benchmark(c: &mut Criterion) {
    let mut c = c.benchmark_group("Fast");
    c.warm_up_time(Duration::new(1, 100));
    c.sample_size(15);

    c.bench_function("system: hypotf", |b| {
        b.iter(|| {
            for i in 1..1000 {
                black_box(f32::hypot(i as f32 / 1000., i as f32 / 324.));
            }
        })
    });

    c.bench_function("pxfm: f_hypotf", |b| {
        b.iter(|| {
            for i in 1..1000 {
                black_box(pxfm::f_hypotf(i as f32 / 1000.0, i as f32 / 324.));
            }
        })
    });

    c.bench_function("pxfm: f_cathetusf", |b| {
        b.iter(|| {
            for i in 1..1000 {
                black_box(pxfm::f_cathetusf(i as f32 / 1000.0, i as f32 / 324.));
            }
        })
    });

    c.bench_function("libm: hypot", |b| {
        b.iter(|| {
            for i in 1..1000 {
                black_box(libm::hypot(i as f64 / 1000.0, i as f64 / 324.));
            }
        })
    });

    c.bench_function("system: hypot", |b| {
        b.iter(|| {
            for i in 1..1000 {
                black_box(f64::hypot(i as f64 / 1000., i as f64 / 324.));
            }
        })
    });

    c.bench_function("pxfm: f_hypot", |b| {
        b.iter(|| {
            for i in 1..1000 {
                black_box(pxfm::f_hypot(i as f64 / 1000.0, i as f64 / 324.));
            }
        })
    });

    c.bench_function("pxfm: f_cathethus", |b| {
        b.iter(|| {
            for i in 1..1000 {
                black_box(pxfm::f_cathetus(i as f64 / 1000.0, i as f64 / 324.));
            }
        })
    });

    c.finish();
}

criterion_group!(benches, criterion_benchmark);
criterion_main!(benches);
