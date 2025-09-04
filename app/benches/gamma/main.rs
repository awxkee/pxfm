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

    c.bench_function("libm: lgamma", |b| {
        b.iter(|| {
            for i in 1..1000 {
                black_box(libm::lgamma(black_box(i as f64 / 1000.0 * 4.)));
            }
        })
    });

    c.bench_function("pxfm: f_lgamma", |b| {
        b.iter(|| {
            for i in 1..1000 {
                black_box(pxfm::f_lgamma(black_box(i as f64 / 1000.0 + 4.)));
            }
        })
    });

    c.bench_function("libm: tgamma", |b| {
        b.iter(|| {
            for i in 1..1000 {
                black_box(libm::tgamma(black_box(i as f64 / 1000.0 * 36.0)));
            }
        })
    });

    c.bench_function("pxfm: f_tgamma", |b| {
        b.iter(|| {
            for i in 1..1000 {
                black_box(pxfm::f_tgamma(black_box(i as f64 / 1000.0 * 36.0)));
            }
        })
    });

    c.bench_function("libm: lgammaf", |b| {
        b.iter(|| {
            for i in 1..1000 {
                black_box(libm::lgammaf(black_box(i as f32 / 1000.0 * 36.0)));
            }
        })
    });

    c.bench_function("pxfm: f_lgammaf", |b| {
        b.iter(|| {
            for i in 1..1000 {
                black_box(pxfm::f_lgammaf(black_box(i as f32 / 1000.0 * 36.0)));
            }
        })
    });

    c.bench_function("libm: tgammaf", |b| {
        b.iter(|| {
            for i in 1..1000 {
                black_box(libm::tgammaf(black_box(i as f32 / 1000.0 * 36.0)));
            }
        })
    });

    c.bench_function("pxfm: f_tgammaf", |b| {
        b.iter(|| {
            for i in 1..1000 {
                black_box(pxfm::f_tgammaf(black_box(i as f32 / 1000.0 * 36.0)));
            }
        })
    });

    c.bench_function("pxfm: f_digamma", |b| {
        b.iter(|| {
            for i in 1..1000 {
                black_box(pxfm::f_digamma(black_box(15. + i as f64 / 1000.0)));
            }
        })
    });

    c.bench_function("pxfm: f_trigamma", |b| {
        b.iter(|| {
            for i in 1..1000 {
                black_box(pxfm::f_trigamma(black_box(15. + i as f64 / 1000.0)));
            }
        })
    });

    c.bench_function("pxfm: f_digammaf", |b| {
        b.iter(|| {
            for i in 1..1000 {
                black_box(pxfm::f_digammaf(black_box(-i as f32 / 1000.0 * 100.0)));
            }
        })
    });

    c.bench_function("pxfm: f_trigammaf", |b| {
        b.iter(|| {
            for i in 1..1000 {
                black_box(pxfm::f_trigammaf(black_box(-i as f32 / 1000.0 * 100.0)));
            }
        })
    });

    c.finish();
}

criterion_group!(benches, criterion_benchmark);
criterion_main!(benches);
