// #![feature(float_erf)]
/*
 * // Copyright 2024 (c) the Radzivon Bartoshyk. All rights reserved.
 * //
 * // Use of this source code is governed by a BSD-style
 * // license that can be found in the LICENSE file.
 */
use criterion::{Criterion, criterion_group, criterion_main};
use pxfm::{f_cbrt, f_cbrtf, f_rcbrt, f_rcbrtf};
use std::hint::black_box;
use std::time::Duration;

pub fn criterion_benchmark(c: &mut Criterion) {
    let mut c = c.benchmark_group("Fast");
    c.warm_up_time(Duration::new(1, 100));
    c.sample_size(15);

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

    c.bench_function("pxfm: cbrt", |b| {
        b.iter(|| {
            for i in 1..1000 {
                black_box(f_cbrt(i as f64));
            }
        })
    });

    c.bench_function("pxfm: rcbrt", |b| {
        b.iter(|| {
            for i in 1..1000 {
                black_box(f_rcbrt(i as f64));
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

    c.bench_function("pxfm: cbrtf", |b| {
        b.iter(|| {
            for i in 1..1000 {
                black_box(f_cbrtf(i as f32));
            }
        })
    });

    c.bench_function("pxfm: rcbrtf", |b| {
        b.iter(|| {
            for i in 1..1000 {
                black_box(f_rcbrtf(i as f32));
            }
        })
    });

    c.finish();
}

criterion_group!(benches, criterion_benchmark);
criterion_main!(benches);
