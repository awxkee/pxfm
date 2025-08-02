use rug::ops::Pow;
use rug::{Assign, Float};
use std::ops::{Add, Div, Mul};

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

        // if term.clone().abs() < Float::with_val(prec, 1e-100) {
        //     break;
        // }

        if terms > 100 {
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
    use rug::Float;

    let mut sum = Float::with_val(prec, 0);
    let mut k_fact = Float::with_val(prec, 1); // k!
    let mut kp1_fact = Float::with_val(prec, 1); // (k+1)!

    let mut x_half = Float::with_val(prec, x);
    x_half /= 2; // accurate (x/2)

    let mut x_pow = x_half.clone(); // (x/2)^(2k+1), start with k=0: (x/2)^1

    let max_terms = 1500;
    let epsilon = Float::with_val(prec, 1e-70);

    for k in 0..max_terms {
        if k > 0 {
            k_fact *= k;
            kp1_fact *= k + 1;
            // multiply power by (x/2)^2
            x_pow *= &x_half;
            x_pow *= &x_half;
        }

        let denom = k_fact.clone().mul(kp1_fact.clone());
        let term = x_pow.clone().div(&denom);

        if term.clone().abs() < epsilon {
            break;
        }

        sum += term;
    }

    sum
}

/// Harmonic number H_n = 1 + 1/2 + ... + 1/n
fn harmonic(n: u32, prec: u32) -> Float {
    let mut sum = Float::with_val(prec, 0);
    for i in 1..=n {
        sum += Float::with_val(prec, 1) / i;
    }
    sum
}

/// Computes Bessel function K₀(x) using series for small x.
pub fn bessel_k0(x: f64, prec: u32) -> Float {
    let mut result = Float::with_val(prec, 0);
    let mut term = Float::with_val(prec, 1);
    let x2 = Float::with_val(prec, x * x).div(&Float::with_val(107, 4));

    let p = Float::parse(
        "0.57721566490153286060651209008240243104215933593992359880576723488486772677",
    )
    .unwrap();

    let gamma = Float::with_val(107, p);
    let i0 = bessel_i0(x, prec);

    // Series sum: ∑ [ (ψ(k+1) + γ) / (k!)² * (x²/4)^k ]
    let mut k_fact = Float::with_val(prec, 1);
    for k in 0..1000 {
        let psi_k1 = harmonic(k, prec) - &gamma;
        if k > 0 {
            k_fact *= k;
        }
        let denom = Float::with_val(prec, &k_fact * &k_fact);
        let num = Float::with_val(prec, &x2.clone().pow(k));
        term.assign(&psi_k1 * num / denom);
        result = result.add(&term.clone());

        if term.clone().abs() < Float::with_val(prec, 1e-70) {
            break;
        }
    }

    let log_term = Float::with_val(prec, x / 2.).ln() * i0;
    -log_term + result
}
