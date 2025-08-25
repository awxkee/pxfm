use rug::float::Constant;
use rug::ops::Pow;
use rug::{Assign, Float};
pub use std::ops::{Add, Div, Mul, Neg};

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

        if terms > 200 {
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
    let mut sum = Float::with_val(prec, 1);
    for i in 2..=n {
        sum += Float::with_val(prec, 1) / i;
    }
    sum
}

/// Computes Bessel function K₀(x) using series for small x.
pub fn bessel_k0(x: f64, prec: u32) -> Float {
    let mut term = Float::with_val(prec, 1);
    let x_over_two = Float::with_val(prec, x).div(&Float::with_val(prec, 2));
    let x2 = x_over_two.clone().mul(x_over_two.clone());

    let euler = Float::with_val(prec, Constant::Euler);
    let i0 = bessel_i0(x, prec);

    let mut result = x_over_two.clone().ln().add(euler.clone()).mul(i0).neg();

    let mut k_fact = Float::with_val(prec, 1);
    for k in 1..150 {
        let harmony = harmonic(k, prec);
        if k > 0 {
            k_fact *= k;
        }
        let denom = Float::with_val(prec, &k_fact * &k_fact);
        let num = x2.clone().pow(k);
        term.assign(&harmony * num / denom);
        result = result.add(&term.clone());
    }
    result
}

pub fn bessel_k1(x: f64, prec: u32) -> Float {
    assert!(
        x > 0.0,
        "K1(x) is singular at x = 0 and undefined for x <= 0 in this routine"
    );

    let dx = x;
    let x = Float::with_val(prec, x);
    let x2 = Float::with_val(prec, &x / 2); // x/2
    let euler = Float::with_val(prec, Constant::Euler); // γ

    // You need consistent I0/I1 implementations at `prec`.
    let i0 = bessel_i0(dx, prec);
    let i1 = bessel_i1(dx, prec);

    // K1(x) = I0(x)/x + (ln(x/2)+γ) I1(x) - Σ_{k>=0} H_{k+1} (x/2)^{2k+1} / (k!(k+1)!)
    let mut result = Float::with_val(prec, &i0 / &x);
    result += Float::with_val(prec, x2.clone().ln() + euler.clone()) * &i1;

    // Accumulate S = Σ H_{k+1} (x/2)^{2k+1} / (k!(k+1)!)
    let mut sum = Float::with_val(prec, 0);
    let mut k_fact = Float::with_val(prec, 1); // k!
    let mut power = Float::with_val(prec, x2.clone()); // (x/2)^{2k+1}, start at k=0 => (x/2)^1

    for k in 0..150 {
        // denom = k! * (k+1)!
        // (k+1)! = (k+1) * k!
        let kp1 = Float::with_val(prec, (k + 1) as i32);
        let denom = Float::with_val(prec, k_fact.clone().mul(k_fact.clone().mul(kp1)));

        // H_{k+1}
        let h_k1 = harmonic(k + 1, prec);

        // term
        let term = Float::with_val(prec, h_k1.mul(power.clone()).div(denom));
        sum += term;

        // advance k!, power
        let next = k + 1;
        if next > 0 {
            k_fact *= next;
        }
        power *= &x2; // (x/2)^{2k+1} -> (x/2)^{2k+2}
        power *= &x2; // -> (x/2)^{2(k+1)+1}
    }

    result -= sum;
    result
}

#[cfg(test)]
mod tests {
    use super::*;
    #[test]
    fn test_k0() {
        println!("{}", bessel_k0(5., 100));
    }

    #[test]
    fn test_k1() {
        println!("{}", bessel_k1(0.01, 100));
    }
}
