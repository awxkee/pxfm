// #![feature(float_erf)]
/*
 * // Copyright 2024 (c) the Radzivon Bartoshyk. All rights reserved.
 * //
 * // Use of this source code is governed by a BSD-style
 * // license that can be found in the LICENSE file.
 */
use criterion::{Criterion, criterion_group, criterion_main};
use pxfm::{
    f_acosh, f_acoshf, f_asinh, f_asinhf, f_atanh, f_atanhf, f_cosh, f_coshf, f_sinh, f_sinhf,
    f_tanh, f_tanhf,
};
use std::hint::black_box;
use std::time::Duration;

pub fn criterion_benchmark(c: &mut Criterion) {
    let mut c = c.benchmark_group("Fast");
    c.warm_up_time(Duration::new(1, 100));
    c.sample_size(15);

    c.bench_function("libm: atanh", |b| {
        b.iter(|| {
            for i in 1..1000 {
                black_box(libm::atanh(black_box(i as f64 / 1000.0)));
            }
        })
    });

    c.bench_function("system: atanh", |b| {
        b.iter(|| {
            for i in 1..1000 {
                black_box(f64::atanh(black_box(i as f64 / 1000.0)));
            }
        })
    });

    c.bench_function("pxfm: f_atanh", |b| {
        b.iter(|| {
            for i in 1..1000 {
                black_box(f_atanh(black_box(i as f64 / 1000.0)));
            }
        })
    });

    c.bench_function("libm: tanh", |b| {
        b.iter(|| {
            for i in 1..1000 {
                black_box(libm::tanh(black_box(i as f64 / 1000.0)));
            }
        })
    });

    c.bench_function("system: tanh", |b| {
        b.iter(|| {
            for i in 1..1000 {
                black_box(f64::tanh(black_box(i as f64 / 1000.0)));
            }
        })
    });

    c.bench_function("pxfm: f_tanh", |b| {
        b.iter(|| {
            for i in 1..1000 {
                black_box(f_tanh(black_box(i as f64 / 1000.0)));
            }
        })
    });

    c.bench_function("libm: atanhf", |b| {
        b.iter(|| {
            for i in 1..1000 {
                black_box(libm::atanhf(black_box(i as f32 / 1000.0)));
            }
        })
    });

    c.bench_function("system: atanhf", |b| {
        b.iter(|| {
            for i in 1..1000 {
                black_box(f32::atanh(black_box(i as f32 / 1000.0)));
            }
        })
    });

    c.bench_function("pxfm: f_atanhf", |b| {
        b.iter(|| {
            for i in 1..1000 {
                black_box(f_atanhf(black_box(i as f32 / 1000.0)));
            }
        })
    });

    c.bench_function("libm: acoshf", |b| {
        b.iter(|| {
            for i in 1..1000 {
                black_box(libm::acoshf(black_box(i as f32 / 1000.0)));
            }
        })
    });

    c.bench_function("system: acoshf", |b| {
        b.iter(|| {
            for i in 1..1000 {
                black_box(f32::acosh(black_box(i as f32 / 1000.0)));
            }
        })
    });

    c.bench_function("pxfm: f_acoshf", |b| {
        b.iter(|| {
            for i in 1..1000 {
                black_box(f_acoshf(black_box(i as f32 / 1000.0)));
            }
        })
    });

    c.bench_function("libm: asinhf", |b| {
        b.iter(|| {
            for i in 1..1000 {
                black_box(libm::asinf(black_box(i as f32 / 1000.0)));
            }
        })
    });

    c.bench_function("system: asinhf", |b| {
        b.iter(|| {
            for i in 1..1000 {
                black_box(f32::asinh(black_box(i as f32 / 1000.0)));
            }
        })
    });

    c.bench_function("pxfm: f_asinhf", |b| {
        b.iter(|| {
            for i in 1..1000 {
                black_box(f_asinhf(black_box(i as f32 / 1000.0)));
            }
        })
    });

    c.bench_function("libm: cosh", |b| {
        b.iter(|| {
            for i in 1..1000 {
                black_box(libm::cosh(black_box(i as f64 / 1000.0)));
            }
        })
    });

    c.bench_function("system: cosh", |b| {
        b.iter(|| {
            for i in 1..1000 {
                black_box(f64::cosh(black_box(i as f64 / 1000.0)));
            }
        })
    });

    c.bench_function("pxfm: f_cosh", |b| {
        b.iter(|| {
            for i in 1..1000 {
                black_box(f_cosh(black_box(i as f64 / 1000.0)));
            }
        })
    });

    c.bench_function("libm: sinh", |b| {
        b.iter(|| {
            for i in 1..1000 {
                black_box(libm::sinh(black_box(i as f64 / 1000.0)));
            }
        })
    });

    c.bench_function("system: sinh", |b| {
        b.iter(|| {
            for i in 1..1000 {
                black_box(f64::sinh(black_box(i as f64 / 1000.0)));
            }
        })
    });

    c.bench_function("pxfm: f_sinh", |b| {
        b.iter(|| {
            for i in 1..1000 {
                black_box(f_sinh(black_box(i as f64 / 1000.0)));
            }
        })
    });

    c.bench_function("libm: asinh", |b| {
        b.iter(|| {
            for i in 1..1000 {
                black_box(libm::asinh(black_box(i as f64 / 1000.0)));
            }
        })
    });

    c.bench_function("system: asinh", |b| {
        b.iter(|| {
            for i in 1..1000 {
                black_box(f64::asinh(black_box(i as f64 / 1000.0)));
            }
        })
    });

    c.bench_function("pxfm: f_asinh", |b| {
        b.iter(|| {
            for i in 1..1000 {
                black_box(f_asinh(black_box(i as f64 / 1000.0)));
            }
        })
    });

    c.bench_function("libm: acosh", |b| {
        b.iter(|| {
            for i in 1..1000 {
                black_box(libm::acosh(black_box(i as f64 / 1000.0)));
            }
        })
    });

    c.bench_function("system: acosh", |b| {
        b.iter(|| {
            for i in 1..1000 {
                black_box(f64::acosh(black_box(i as f64 / 1000.0)));
            }
        })
    });

    c.bench_function("pxfm: f_acosh", |b| {
        b.iter(|| {
            for i in 1..1000 {
                black_box(f_acosh(black_box(i as f64 / 1000.0)));
            }
        })
    });

    c.bench_function("libm::tanhf", |b| {
        b.iter(|| {
            for i in 1..1000 {
                black_box(libm::tanhf(black_box(i as f32)));
            }
        })
    });

    c.bench_function("system::tanhf", |b| {
        b.iter(|| {
            for i in 1..1000 {
                black_box(f32::tanh(black_box(i as f32)));
            }
        })
    });

    c.bench_function("pxfm: FMA tanhf", |b| {
        b.iter(|| {
            for i in 1..1000 {
                black_box(f_tanhf(black_box(i as f32)));
            }
        })
    });

    c.bench_function("libm::sinhf", |b| {
        b.iter(|| {
            for i in 1..1000 {
                black_box(libm::sinhf(black_box(i as f32)));
            }
        })
    });

    c.bench_function("system::sinhf", |b| {
        b.iter(|| {
            for i in 1..1000 {
                black_box(f32::sinh(black_box(i as f32)));
            }
        })
    });

    c.bench_function("pxfm: FMA sinhf", |b| {
        b.iter(|| {
            for i in 1..1000 {
                black_box(f_sinhf(black_box(i as f32)));
            }
        })
    });

    c.bench_function("libm::coshf", |b| {
        b.iter(|| {
            for i in 1..1000 {
                black_box(libm::coshf(black_box(i as f32)));
            }
        })
    });

    c.bench_function("system::coshf", |b| {
        b.iter(|| {
            for i in 1..1000 {
                black_box(f32::cosh(black_box(i as f32)));
            }
        })
    });

    c.bench_function("pxfm: FMA coshf", |b| {
        b.iter(|| {
            for i in 1..1000 {
                black_box(f_coshf(black_box(i as f32)));
            }
        })
    });

    c.finish();
}

criterion_group!(benches, criterion_benchmark);
criterion_main!(benches);
