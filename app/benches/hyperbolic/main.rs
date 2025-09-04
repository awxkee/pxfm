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
                black_box(libm::atanh(i as f64 / 1000.0));
            }
        })
    });

    c.bench_function("system: atanh", |b| {
        b.iter(|| {
            for i in 1..1000 {
                black_box(f64::atanh(i as f64 / 1000.0));
            }
        })
    });

    c.bench_function("pxfm: f_atanh", |b| {
        b.iter(|| {
            for i in 1..1000 {
                black_box(f_atanh(i as f64 / 1000.0));
            }
        })
    });

    c.bench_function("libm: tanh", |b| {
        b.iter(|| {
            for i in 1..1000 {
                black_box(libm::tanh(i as f64 / 1000.0));
            }
        })
    });

    c.bench_function("system: tanh", |b| {
        b.iter(|| {
            for i in 1..1000 {
                black_box(f64::tanh(i as f64 / 1000.0));
            }
        })
    });

    c.bench_function("pxfm: f_tanh", |b| {
        b.iter(|| {
            for i in 1..1000 {
                black_box(f_tanh(i as f64 / 1000.0));
            }
        })
    });

    c.bench_function("libm: atanhf", |b| {
        b.iter(|| {
            for i in 1..1000 {
                black_box(libm::atanhf(i as f32 / 1000.0));
            }
        })
    });

    c.bench_function("system: atanhf", |b| {
        b.iter(|| {
            for i in 1..1000 {
                black_box(f32::atanh(i as f32 / 1000.0));
            }
        })
    });

    c.bench_function("pxfm: f_atanhf", |b| {
        b.iter(|| {
            for i in 1..1000 {
                black_box(f_atanhf(i as f32 / 1000.0));
            }
        })
    });

    c.bench_function("libm: acoshf", |b| {
        b.iter(|| {
            for i in 1..1000 {
                black_box(libm::acoshf(i as f32 / 1000.0));
            }
        })
    });

    c.bench_function("system: acoshf", |b| {
        b.iter(|| {
            for i in 1..1000 {
                black_box(f32::acosh(i as f32 / 1000.0));
            }
        })
    });

    c.bench_function("pxfm: f_acoshf", |b| {
        b.iter(|| {
            for i in 1..1000 {
                black_box(f_acoshf(i as f32 / 1000.0));
            }
        })
    });

    c.bench_function("libm: asinhf", |b| {
        b.iter(|| {
            for i in 1..1000 {
                black_box(libm::asinf(i as f32 / 1000.0));
            }
        })
    });

    c.bench_function("system: asinhf", |b| {
        b.iter(|| {
            for i in 1..1000 {
                black_box(f32::asinh(i as f32 / 1000.0));
            }
        })
    });

    c.bench_function("pxfm: f_asinhf", |b| {
        b.iter(|| {
            for i in 1..1000 {
                black_box(f_asinhf(i as f32 / 1000.0));
            }
        })
    });

    c.bench_function("libm: cosh", |b| {
        b.iter(|| {
            for i in 1..1000 {
                black_box(libm::cosh(i as f64 / 1000.0));
            }
        })
    });

    c.bench_function("system: cosh", |b| {
        b.iter(|| {
            for i in 1..1000 {
                black_box(f64::cosh(i as f64 / 1000.0));
            }
        })
    });

    c.bench_function("pxfm: f_cosh", |b| {
        b.iter(|| {
            for i in 1..1000 {
                black_box(f_cosh(i as f64 / 1000.0));
            }
        })
    });

    c.bench_function("libm: sinh", |b| {
        b.iter(|| {
            for i in 1..1000 {
                black_box(libm::sinh(i as f64 / 1000.0));
            }
        })
    });

    c.bench_function("system: sinh", |b| {
        b.iter(|| {
            for i in 1..1000 {
                black_box(f64::sinh(i as f64 / 1000.0));
            }
        })
    });

    c.bench_function("pxfm: f_sinh", |b| {
        b.iter(|| {
            for i in 1..1000 {
                black_box(f_sinh(i as f64 / 1000.0));
            }
        })
    });

    c.bench_function("libm: asinh", |b| {
        b.iter(|| {
            for i in 1..1000 {
                black_box(libm::asinh(i as f64 / 1000.0));
            }
        })
    });

    c.bench_function("system: asinh", |b| {
        b.iter(|| {
            for i in 1..1000 {
                black_box(f64::asinh(i as f64 / 1000.0));
            }
        })
    });

    c.bench_function("pxfm: f_asinh", |b| {
        b.iter(|| {
            for i in 1..1000 {
                black_box(f_asinh(i as f64 / 1000.0));
            }
        })
    });

    c.bench_function("libm: acosh", |b| {
        b.iter(|| {
            for i in 1..1000 {
                black_box(libm::acosh(i as f64 / 1000.0));
            }
        })
    });

    c.bench_function("system: acosh", |b| {
        b.iter(|| {
            for i in 1..1000 {
                black_box(f64::acosh(i as f64 / 1000.0));
            }
        })
    });

    c.bench_function("pxfm: f_acosh", |b| {
        b.iter(|| {
            for i in 1..1000 {
                black_box(f_acosh(i as f64 / 1000.0));
            }
        })
    });

    c.bench_function("libm::tanhf", |b| {
        b.iter(|| {
            for i in 1..1000 {
                black_box(libm::tanhf(i as f32));
            }
        })
    });

    c.bench_function("system::tanhf", |b| {
        b.iter(|| {
            for i in 1..1000 {
                black_box(f32::tanh(i as f32));
            }
        })
    });

    c.bench_function("pxfm: FMA tanhf", |b| {
        b.iter(|| {
            for i in 1..1000 {
                black_box(f_tanhf(i as f32));
            }
        })
    });

    c.bench_function("libm::sinhf", |b| {
        b.iter(|| {
            for i in 1..1000 {
                black_box(libm::sinhf(i as f32));
            }
        })
    });

    c.bench_function("system::sinhf", |b| {
        b.iter(|| {
            for i in 1..1000 {
                black_box(f32::sinh(i as f32));
            }
        })
    });

    c.bench_function("pxfm: FMA sinhf", |b| {
        b.iter(|| {
            for i in 1..1000 {
                black_box(f_sinhf(i as f32));
            }
        })
    });

    c.bench_function("libm::coshf", |b| {
        b.iter(|| {
            for i in 1..1000 {
                black_box(libm::coshf(i as f32));
            }
        })
    });

    c.bench_function("system::coshf", |b| {
        b.iter(|| {
            for i in 1..1000 {
                black_box(f32::cosh(i as f32));
            }
        })
    });

    c.bench_function("pxfm: FMA coshf", |b| {
        b.iter(|| {
            for i in 1..1000 {
                black_box(f_coshf(i as f32));
            }
        })
    });

    c.finish();
}

criterion_group!(benches, criterion_benchmark);
criterion_main!(benches);
