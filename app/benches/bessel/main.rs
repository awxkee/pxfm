// #![feature(float_erf)]
/*
 * // Copyright 2024 (c) the Radzivon Bartoshyk. All rights reserved.
 * //
 * // Use of this source code is governed by a BSD-style
 * // license that can be found in the LICENSE file.
 */
use criterion::{Criterion, criterion_group, criterion_main};
use pxfm::{f_j0, f_j0f, f_j1, f_j1f, f_y0, f_y0f, f_y1, f_y1f};
use std::hint::black_box;
use std::time::Duration;

pub fn criterion_benchmark(c: &mut Criterion) {
    let mut c = c.benchmark_group("Fast");
    c.warm_up_time(Duration::new(1, 100));
    c.sample_size(15);

    c.bench_function("pxfm: i2f", |b| {
        b.iter(|| {
            for i in 1..1000 {
                black_box(pxfm::f_i2f(black_box(i as f32 / 50.0)));
            }
        })
    });

    c.bench_function("pxfm: i1ef", |b| {
        b.iter(|| {
            for i in 1..1000 {
                black_box(pxfm::f_i1ef(black_box(i as f32 / 50.0)));
            }
        })
    });

    c.bench_function("pxfm: i1f", |b| {
        b.iter(|| {
            for i in 1..1000 {
                black_box(pxfm::f_i1f(black_box(i as f32 / 50.0)));
            }
        })
    });

    c.bench_function("pxfm: i0ef", |b| {
        b.iter(|| {
            for i in 1..1000 {
                black_box(pxfm::f_i0ef(black_box(i as f32 / 50.0)));
            }
        })
    });

    c.bench_function("pxfm: i0f", |b| {
        b.iter(|| {
            for i in 1..1000 {
                black_box(pxfm::f_i0f(black_box(i as f32 / 50.0)));
            }
        })
    });

    c.bench_function("pxfm: i0e", |b| {
        b.iter(|| {
            for i in 1..1000 {
                black_box(pxfm::f_i0e(black_box(i as f64 / 50.0)));
            }
        })
    });

    c.bench_function("pxfm: i0", |b| {
        b.iter(|| {
            for i in 1..1000 {
                black_box(pxfm::f_i0(black_box(i as f64 / 50.0)));
            }
        })
    });

    c.bench_function("pxfm: i1", |b| {
        b.iter(|| {
            for i in 1..1000 {
                black_box(pxfm::f_i1(black_box(i as f64 / 50.0)));
            }
        })
    });

    c.bench_function("pxfm: i2", |b| {
        b.iter(|| {
            for i in 1..1000 {
                black_box(pxfm::f_i2(black_box(i as f64 / 50.0)));
            }
        })
    });

    c.bench_function("pxfm: k2f", |b| {
        b.iter(|| {
            for i in 1..1000 {
                black_box(pxfm::f_k2f(black_box(i as f32 / 50.0)));
            }
        })
    });

    c.bench_function("pxfm: k1f", |b| {
        b.iter(|| {
            for i in 1..1000 {
                black_box(pxfm::f_k1f(black_box(i as f32 / 50.0)));
            }
        })
    });

    c.bench_function("pxfm: k0ef", |b| {
        b.iter(|| {
            for i in 1..1000 {
                black_box(pxfm::f_k0ef(black_box(i as f32 / 50.0)));
            }
        })
    });

    c.bench_function("pxfm: k0f", |b| {
        b.iter(|| {
            for i in 1..1000 {
                black_box(pxfm::f_k0f(black_box(i as f32 / 50.0)));
            }
        })
    });

    c.bench_function("pxfm: k0e", |b| {
        b.iter(|| {
            for i in 1..1000 {
                black_box(pxfm::f_k0e(black_box(i as f64 / 50.0)));
            }
        })
    });

    c.bench_function("pxfm: k0", |b| {
        b.iter(|| {
            for i in 1..1000 {
                black_box(pxfm::f_k0(black_box(i as f64 / 50.0)));
            }
        })
    });

    c.bench_function("pxfm: k1e", |b| {
        b.iter(|| {
            for i in 1..1000 {
                black_box(pxfm::f_k1e(black_box(i as f64 / 50.0)));
            }
        })
    });

    c.bench_function("pxfm: k1", |b| {
        b.iter(|| {
            for i in 1..1000 {
                black_box(pxfm::f_k1(black_box(i as f64 / 50.0)));
            }
        })
    });

    c.bench_function("libm::j0f", |b| {
        b.iter(|| {
            for i in 1..1000 {
                black_box(libm::j0f(black_box(i as f32 / 1000.)));
            }
        })
    });

    c.bench_function("pxfm: f_j0f", |b| {
        b.iter(|| {
            for i in 1..1000 {
                black_box(f_j0f(black_box(i as f32 / 10.0)));
            }
        })
    });

    c.bench_function("libm::j1f", |b| {
        b.iter(|| {
            for i in 1..1000 {
                black_box(libm::j1f(black_box(i as f32 / 10.0)));
            }
        })
    });

    c.bench_function("libm::y0f", |b| {
        b.iter(|| {
            for i in 1..1000 {
                black_box(libm::y0f(black_box(i as f32 / 10.)));
            }
        })
    });

    c.bench_function("pxfm: f_y0f", |b| {
        b.iter(|| {
            for i in 1..1000 {
                black_box(f_y0f(black_box(i as f32 / 10.)));
            }
        })
    });

    c.bench_function("libm::y1f", |b| {
        b.iter(|| {
            for i in 1..1000 {
                black_box(libm::y1f(black_box(i as f32 / 10.)));
            }
        })
    });

    c.bench_function("pxfm: f_y1f", |b| {
        b.iter(|| {
            for i in 1..1000 {
                black_box(f_y1f(black_box(i as f32 / 10.)));
            }
        })
    });

    c.bench_function("pxfm: f_jincpi", |b| {
        b.iter(|| {
            for i in 1..1000 {
                black_box(pxfm::f_jincpi(black_box(i as f64 / 500.0)));
            }
        })
    });

    c.bench_function("pxfm: f_jincpif", |b| {
        b.iter(|| {
            for i in 1..1000 {
                black_box(pxfm::f_jincpif(black_box(i as f32 / 10.0)));
            }
        })
    });

    c.bench_function("pxfm: f_j1f", |b| {
        b.iter(|| {
            for i in 1..1000 {
                black_box(f_j1f(black_box(i as f32 / 10.0)));
            }
        })
    });

    c.bench_function("libm::j0", |b| {
        b.iter(|| {
            for i in 1..1000 {
                black_box(libm::j0(black_box(i as f64 / 100.)));
            }
        })
    });

    c.bench_function("pxfm: j0", |b| {
        b.iter(|| {
            for i in 1..1000 {
                black_box(f_j0(black_box(i as f64 / 100.)));
            }
        })
    });

    c.bench_function("libm::y0", |b| {
        b.iter(|| {
            for i in 1..1000 {
                black_box(libm::y0(black_box(i as f64 / 100.)));
            }
        })
    });

    c.bench_function("pxfm: y0", |b| {
        b.iter(|| {
            for i in 1..1000 {
                black_box(f_y0(black_box(i as f64 / 100.)));
            }
        })
    });

    c.bench_function("libm::j1", |b| {
        b.iter(|| {
            for i in 1..1000 {
                black_box(libm::j1(black_box(i as f64 / 100.)));
            }
        })
    });

    c.bench_function("pxfm: j1", |b| {
        b.iter(|| {
            for i in 1..1000 {
                black_box(f_j1(black_box(i as f64 / 100.)));
            }
        })
    });

    c.bench_function("libm::y1", |b| {
        b.iter(|| {
            for i in 1..1000 {
                black_box(libm::y1(black_box(i as f64 / 100.)));
            }
        })
    });

    c.bench_function("pxfm: f_y1", |b| {
        b.iter(|| {
            for i in 1..1000 {
                black_box(f_y1(black_box(i as f64 / 100.)));
            }
        })
    });

    c.finish();
}

criterion_group!(benches, criterion_benchmark);
criterion_main!(benches);
