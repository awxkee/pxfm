// #![feature(float_erf)]
/*
 * // Copyright 2024 (c) the Radzivon Bartoshyk. All rights reserved.
 * //
 * // Use of this source code is governed by a BSD-style
 * // license that can be found in the LICENSE file.
 */
use criterion::{Criterion, criterion_group, criterion_main};
use pxfm::{
    f_log, f_log1p, f_log1pf, f_log2, f_log2f, f_log2p1, f_log2p1f, f_log10, f_log10f, f_log10p1,
    f_log10p1f, f_logf,
};
use std::hint::black_box;
use std::time::Duration;

pub fn criterion_benchmark(c: &mut Criterion) {
    let mut c = c.benchmark_group("Fast");
    c.warm_up_time(Duration::new(1, 100));
    c.sample_size(15);

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

    c.bench_function("libm: log1p", |b| {
        b.iter(|| {
            for i in 1..1000 {
                black_box(libm::log1p(i as f64 / 1000.0));
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

    c.bench_function("pxfm: f_log1pmxf", |b| {
        b.iter(|| {
            for i in 1..1000 {
                black_box(pxfm::f_log1pmxf(i as f32 / 1000.0));
            }
        })
    });

    c.bench_function("libm::log10", |b| {
        b.iter(|| {
            for i in 1..1000 {
                black_box(libm::log10(i as f64 / 1000.0));
            }
        })
    });

    c.bench_function("system: log10", |b| {
        b.iter(|| {
            for i in 1..1000 {
                black_box(f64::log10(i as f64 / 1000.0));
            }
        })
    });

    c.bench_function("pxfm: log10", |b| {
        b.iter(|| {
            for i in 1..1000 {
                black_box(f_log10(i as f64 / 1000.0));
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

    c.bench_function("pxfm: log2", |b| {
        b.iter(|| {
            for i in 1..1000 {
                black_box(f_log2(i as f64));
            }
        })
    });

    c.bench_function("libm::log", |b| {
        b.iter(|| {
            for i in 1..1000 {
                black_box(libm::log(i as f64 / 1000.0));
            }
        })
    });

    c.bench_function("system: log", |b| {
        b.iter(|| {
            for i in 1..1000 {
                black_box(f64::ln(i as f64 / 1000.0));
            }
        })
    });

    c.bench_function("pxfm: log", |b| {
        b.iter(|| {
            for i in 1..1000 {
                black_box(f_log(i as f64 / 1000.0));
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

    c.finish();
}

criterion_group!(benches, criterion_benchmark);
criterion_main!(benches);
