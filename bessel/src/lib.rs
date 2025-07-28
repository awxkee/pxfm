use rug::{Assign, Float};
use std::ops::Mul;

pub fn bessel_i0(x: f64, prec: u32) -> Float {
    let mut sum = Float::with_val(prec, 1); // First term: k = 0
    let mut term = Float::with_val(prec, 1);
    let mut k = 1u32;

    let half_x = Float::with_val(prec, x / 2.);

    let mut num = Float::with_val(prec, &half_x * &half_x); // (x/2)^2

    let mut k_fact_sq = Float::with_val(prec, 1); // (k!)^2

    let mut terms = 0usize;

    loop {
        // (x/2)^(2k) already stored in num
        // (k!)^2 stored in k_fact_sq
        term.assign(&num / &k_fact_sq);

        if term.clone().abs() < Float::with_val(prec, 1e-21) {
            break;
        }

        if terms > 1500 {
            break;
        }

        sum += &term;

        // Update num and den for next iteration
        k += 1;

        num *= &half_x;
        num *= &half_x;

        k_fact_sq *= k;
        k_fact_sq *= k;
        terms += 1;
    }

    sum
}

// Computes wrong values for small values
pub fn bessel_i1(x: f64, prec: u32) -> Float {
    let mut sum = Float::with_val(prec, 0);
    let mut k_fact = Float::with_val(prec, 1); // k!
    let mut kp1_fact = Float::with_val(prec, 1); // (k+1)!
    let x_half = Float::with_val(prec, x / 2.);
    let mut x_pow = Float::with_val(prec, &x_half); // (x/2)^(2k+1), starts at 1st power

    let max_terms = 1500;
    let epsilon = Float::with_val(prec, 1e-21);

    for k in 0..max_terms {
        if k > 0 {
            k_fact *= k;
            kp1_fact *= k + 1;
            x_pow *= &x_half;
            x_pow *= &x_half;
        }

        let denom = k_fact.clone().mul(kp1_fact.clone());
        let term = Float::with_val(prec, &x_pow / denom);

        if term.clone().abs() < epsilon {
            break;
        }

        sum += term;
    }

    sum
}
