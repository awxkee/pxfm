// #![feature(float_erf)]
/*
 * // Copyright 2024 (c) the Radzivon Bartoshyk. All rights reserved.
 * //
 * // Use of this source code is governed by a BSD-style
 * // license that can be found in the LICENSE file.
 */
use criterion::{Criterion, criterion_group, criterion_main};
use pxfm::{f_erf, f_erfc, f_erfcf, f_erff};
use std::hint::black_box;
use std::time::Duration;

pub fn criterion_benchmark(c: &mut Criterion) {
    let mut c = c.benchmark_group("Fast");
    c.warm_up_time(Duration::new(1, 100));
    c.sample_size(15);

    c.bench_function("pxfm: f_erfinvc", |b| {
        b.iter(|| {
            for i in 1..1000 {
                black_box(pxfm::f_erfcinv(black_box(i as f64 / 1000.0)));
            }
        })
    });

    c.bench_function("pxfm: f_erfinv", |b| {
        b.iter(|| {
            for i in 1..1000 {
                black_box(pxfm::f_erfinv(black_box(i as f64 / 1000.0)));
            }
        })
    });

    c.bench_function("pxfm: f_erfinvf", |b| {
        b.iter(|| {
            for i in 1..1000 {
                black_box(pxfm::f_erfinvf(black_box(i as f32 / 1000.0)));
            }
        })
    });

    c.bench_function("pxfm: f_erfcinvf", |b| {
        b.iter(|| {
            for i in 1..1000 {
                black_box(pxfm::f_erfcinvf(black_box(i as f32 / 1000.0)));
            }
        })
    });

    c.bench_function("pxfm: erfcx", |b| {
        b.iter(|| {
            for i in 1..1000 {
                black_box(pxfm::f_erfcx(black_box(black_box(i as f64 / 1000.0 + 8.))));
            }
        })
    });

    c.bench_function("pxfm: erfcxf", |b| {
        b.iter(|| {
            for i in 1..1000 {
                black_box(pxfm::f_erfcxf(black_box(black_box(i as f32 / 100.0))));
            }
        })
    });

    c.bench_function("libm: erfc", |b| {
        b.iter(|| {
            for i in 1..1000 {
                black_box(libm::erfc(black_box(i as f64 / 1000.0)));
            }
        })
    });

    // c.bench_function("system: erfcf", |b| {
    //     b.iter(|| {
    //         for i in 1..1000 {
    //             black_box(f64::erfc(i as f64 / 1000.0));
    //         }
    //     })
    // });

    c.bench_function("pxfm: f_erfc", |b| {
        b.iter(|| {
            for i in 1..1000 {
                black_box(f_erfc(black_box(i as f64 / 1000.0)));
            }
        })
    });

    c.bench_function("libm: erfcf", |b| {
        b.iter(|| {
            for i in 1..1000 {
                black_box(libm::erfcf(black_box(i as f32 / 1000.0)));
            }
        })
    });

    // c.bench_function("system: erfcf", |b| {
    //     b.iter(|| {
    //         for i in 1..1000 {
    //             black_box(f32::erfc(i as f32 / 1000.0));
    //         }
    //     })
    // });

    c.bench_function("pxfm: f_erfcf", |b| {
        b.iter(|| {
            for i in 1..1000 {
                black_box(f_erfcf(black_box(i as f32 / 1000.0)));
            }
        })
    });

    c.bench_function("libm: erf", |b| {
        b.iter(|| {
            for i in 1..1000 {
                black_box(libm::erf(black_box(i as f64 / 1000.0)));
            }
        })
    });

    // c.bench_function("system: erf", |b| {
    //     b.iter(|| {
    //         for i in 1..1000 {
    //             black_box(f64::erf(i as f64 / 1000.0));
    //         }
    //     })
    // });

    c.bench_function("pxfm: f_erf", |b| {
        b.iter(|| {
            for i in 1..1000 {
                black_box(f_erf(black_box(i as f64 / 1000.0)));
            }
        })
    });

    c.bench_function("libm: erff", |b| {
        b.iter(|| {
            for i in 1..1000 {
                black_box(libm::erff(black_box(i as f32 / 1000.0)));
            }
        })
    });

    c.bench_function("pxfm: rerf", |b| {
        b.iter(|| {
            for i in 1..1000 {
                black_box(pxfm::f_rerf(black_box(i as f64 / 1000.0)));
            }
        })
    });

    // c.bench_function("system: erff", |b| {
    //     b.iter(|| {
    //         for i in 1..1000 {
    //             black_box(f32::erf(i as f32 / 1000.0));
    //         }
    //     })
    // });

    c.bench_function("pxfm: f_erff", |b| {
        b.iter(|| {
            for i in 1..1000 {
                black_box(f_erff(black_box(i as f32 / 1000.0)));
            }
        })
    });

    c.finish();
}

criterion_group!(benches, criterion_benchmark);
criterion_main!(benches);
