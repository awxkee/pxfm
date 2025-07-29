/*
 * // Copyright (c) Radzivon Bartoshyk 7/2025. All rights reserved.
 * //
 * // Redistribution and use in source and binary forms, with or without modification,
 * // are permitted provided that the following conditions are met:
 * //
 * // 1.  Redistributions of source code must retain the above copyright notice, this
 * // list of conditions and the following disclaimer.
 * //
 * // 2.  Redistributions in binary form must reproduce the above copyright notice,
 * // this list of conditions and the following disclaimer in the documentation
 * // and/or other materials provided with the distribution.
 * //
 * // 3.  Neither the name of the copyright holder nor the names of its
 * // contributors may be used to endorse or promote products derived from
 * // this software without specific prior written permission.
 * //
 * // THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
 * // AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
 * // IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
 * // DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE
 * // FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL
 * // DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR
 * // SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER
 * // CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY,
 * // OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
 * // OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
 */
use crate::bessel::j0::j0_maclaurin_series;
use crate::bessel::y0_coeffs::{LOG_NEG_DD, Y0_COEFFS};
use crate::bessel::y0f_coeffs::{Y0_ZEROS, Y0_ZEROS_VALUES};
use crate::double_double::DoubleDouble;
use crate::polyeval::{f_polyeval15, f_polyeval35};
use crate::pow_tables::POW_INVERSE;
use crate::sin_helper::sin_dd_small;
use crate::sincos_reduce::{AngleReduced, rem2pi_any};

/// Bessel of the second kind of order 0 (Y0)
///
/// Max ULP 0.5
pub fn f_y0(x: f64) -> f64 {
    if x < 0. {
        return f64::NAN;
    }

    if !x.is_normal() {
        if x == 0. {
            return f64::NEG_INFINITY;
        }

        if x.is_nan() {
            return x + x;
        }

        if x.is_infinite() {
            if x.is_sign_negative() {
                return f64::NAN;
            }
            return 0.;
        }
    }

    if x.abs() <= 1.35 {
        return y0_near_zero(x);
    }

    if x.abs() <= 77. {
        return y0_small_argument_path(x);
    }

    // Exceptions
    //TODO: This actually may be handled
    let bx = x.to_bits();
    if bx == 0x6e7c1d741dc52512u64 {
        return f64::from_bits(0x2696f860815bc669);
    }

    y0_asympt(x)
}

#[inline(always)]
pub(crate) fn log_poly(z: f64) -> DoubleDouble {
    /*
      See ./notes/y0_log.sollya
    */
    const P: [(u64, u64); 10] = [
        (0x3c755555556795ff, 0x3fd5555555555555),
        (0xba86f68980000000, 0xbfd0000000000000),
        (0xbc699b285263b391, 0x3fc999999999999a),
        (0xbc65526cf5c49dc3, 0xbfc5555555555555),
        (0xbc34fc756f340748, 0x3fc24924924924aa),
        (0xbc52e654a63b293e, 0xbfc0000000000023),
        (0xbc5c73c13a9c2171, 0x3fbc71c71c2042d5),
        (0xbc3d2af5e7ee68d8, 0xbfb999999934f78b),
        (0xbc590d76077808da, 0x3fb74612a55c3e99),
        (0xbc3447161ca8047c, 0xbfb5559a592aadc7),
    ];
    let x2 = DoubleDouble::from_exact_mult(z, z);
    let mut t = DoubleDouble::mul_f64_add(
        DoubleDouble::from_bit_pair(P[9]),
        z,
        DoubleDouble::from_bit_pair(P[8]),
    );
    t = DoubleDouble::mul_f64_add(t, z, DoubleDouble::from_bit_pair(P[7]));
    t = DoubleDouble::mul_f64_add(t, z, DoubleDouble::from_bit_pair(P[6]));
    t = DoubleDouble::mul_f64_add(t, z, DoubleDouble::from_bit_pair(P[5]));
    t = DoubleDouble::mul_f64_add(t, z, DoubleDouble::from_bit_pair(P[4]));
    t = DoubleDouble::mul_f64_add(t, z, DoubleDouble::from_bit_pair(P[3]));
    t = DoubleDouble::mul_f64_add(t, z, DoubleDouble::from_bit_pair(P[2]));
    t = DoubleDouble::mul_f64_add(t, z, DoubleDouble::from_bit_pair(P[1]));
    t = DoubleDouble::mul_f64_add(t, z, DoubleDouble::from_bit_pair(P[0]));
    t = DoubleDouble::quick_mult(t, x2);
    t = DoubleDouble::quick_mult_f64(t, z);
    DoubleDouble::mul_f64_add(x2, -0.5, t)
}

#[inline]
pub(crate) fn log_dd(x: f64) -> DoubleDouble {
    let x_u = x.to_bits();
    let mut m = x_u & 0xfffffffffffff;
    let mut e: i64 = ((x_u >> 52) & 0x7ff) as i64;

    let t;
    if e != 0 {
        t = m | (0x3ffu64 << 52);
        m = m.wrapping_add(1u64 << 52);
        e -= 0x3ff;
    } else {
        /* x is a subnormal double  */
        let k = m.leading_zeros() - 11;

        e = -0x3fei64 - k as i64;
        m = m.wrapping_shl(k);
        t = m | (0x3ffu64 << 52);
    }

    /* now |x| = 2^_e*_t = 2^(_e-52)*m with 1 <= _t < 2,
    and 2^52 <= _m < 2^53 */

    //   log(x) = log(t) + E · log(2)
    let mut t = f64::from_bits(t);

    // If m > sqrt(2) we divide it by 2 so ensure 1/sqrt(2) < t < sqrt(2)
    let c: usize = (m >= 0x16a09e667f3bcd) as usize;
    static CY: [f64; 2] = [1.0, 0.5];
    static CM: [u64; 2] = [44, 45];

    e = e.wrapping_add(c as i64);
    let be = e;
    let i = m >> CM[c];
    t *= CY[c];

    let r = f64::from_bits(POW_INVERSE[(i - 181) as usize]);
    let log_r = DoubleDouble::from_bit_pair(LOG_NEG_DD[(i - 181) as usize]);

    let z = f64::mul_add(r, t, -1.0);

    const LOG2_DD: DoubleDouble = DoubleDouble::new(
        f64::from_bits(0x3c7abc9e3b39803f),
        f64::from_bits(0x3fe62e42fefa39ef),
    );

    let tt = DoubleDouble::mul_f64_add(LOG2_DD, be as f64, log_r);

    let v = DoubleDouble::full_add_f64(tt, z);
    let p = log_poly(z);
    DoubleDouble::f64_add(v.hi, DoubleDouble::new(v.lo + p.lo, p.hi))
}

/**
Generated by SageMath:
Evaluates:
Y0(x) = 2/pi*(euler_gamma + log(x/2))*J0(x) - sum((-1)^m*(x/2)^(2*m)/(m!)^2*sum(1+1/2 + ... 1/m))
expressed as:
Y0(x)=log(x)*W0(x) - Z0(x)
```python
from sage.all import *

R = LaurentSeriesRing(RealField(300), 'x',default_prec=300)
x = R.gen()
N = 10  # Number of terms (adjust as needed)
gamma = RealField(300)(euler_gamma)
d2 = RealField(300)(2)
pi = RealField(300).pi()

# Define J0(x) Taylor expansion at x = 0
def j_series(n, x):
    return sum([(-1)**m * (x/2)**(ZZ(n) + ZZ(2)*ZZ(m)) / (ZZ(m).factorial() * (ZZ(m) + ZZ(n)).factorial()) for m in range(N)])

J0_series = j_series(0, x)

def z_series(x):
    return sum([(-1)**m * (x/2)**(ZZ(2)*ZZ(m)) / ZZ(m).factorial()**ZZ(2) * sum(RealField(300)(1)/RealField(300)(k) for k in range(1, m+1)) for m in range(1, N)])

W0 = (d2/pi) * J0_series
Z0 = -gamma * (d2/pi) * J0_series + RealField(300)(2).log() * (d2/pi) * J0_series + (d2/pi) * z_series(x)

# see the series
print(W0)
print(Z0)
```
**/
#[inline]
fn y0_near_zero(x: f64) -> f64 {
    const W: [(u64, u64); 15] = [
        (0xbc86b01ec5417056, 0x3fe45f306dc9c883),
        (0x3c66b01ec5417056, 0xbfc45f306dc9c883),
        (0xbc26b01ec5417056, 0x3f845f306dc9c883),
        (0xbbd67fe4a5feb897, 0xbf321bb945252402),
        (0x3b767fe4a5feb897, 0x3ed21bb945252402),
        (0xbaf5c2495706f745, 0xbe672db9f21b0f5f),
        (0x3a90c8209874dfad, 0x3df49a6c656d62ff),
        (0x3a12921e91b07dd0, 0xbd7ae90af76a4d0f),
        (0xb992921e91b07dd0, 0x3cfae90af76a4d0f),
        (0x39089b0d8a9228ca, 0xbc754331c053fdad),
        (0x3878d321ddfd3c6e, 0x3beb3749ebf0a0dd),
        (0x37e77548130d809b, 0xbb5cca5ae46eae67),
        (0xb73a848e7ca1c943, 0x3ac9976d3cd4293f),
        (0xb6c884706195a054, 0xba336206ff1ce731),
        (0xb6387a7d2389630d, 0x39995103e9f1818f),
    ];
    let x2 = DoubleDouble::from_exact_mult(x, x);
    let w = f_polyeval15(
        x2,
        DoubleDouble::from_bit_pair(W[0]),
        DoubleDouble::from_bit_pair(W[1]),
        DoubleDouble::from_bit_pair(W[2]),
        DoubleDouble::from_bit_pair(W[3]),
        DoubleDouble::from_bit_pair(W[4]),
        DoubleDouble::from_bit_pair(W[5]),
        DoubleDouble::from_bit_pair(W[6]),
        DoubleDouble::from_bit_pair(W[7]),
        DoubleDouble::from_bit_pair(W[8]),
        DoubleDouble::from_bit_pair(W[9]),
        DoubleDouble::from_bit_pair(W[10]),
        DoubleDouble::from_bit_pair(W[11]),
        DoubleDouble::from_bit_pair(W[12]),
        DoubleDouble::from_bit_pair(W[13]),
        DoubleDouble::from_bit_pair(W[14]),
    );

    const Z: [(u64, u64); 15] = [
        (0xbc5ddfd831a70821, 0x3fb2e4d699cbd01f),
        (0xbc6d93e63489aea6, 0xbfc6bbcb41034286),
        (0xbc1b88525c2e130b, 0x3f9075b1bbf41364),
        (0x3be097334e26e578, 0xbf41a6206b7b973d),
        (0x3b51c64a34c78cda, 0x3ee3e99794203bbd),
        (0xbb1c407b0f5b2805, 0xbe7bce4a600d3ea4),
        (0xbaa57d1e1e88c9ca, 0x3e0a6ee796b871b6),
        (0x3a3b6e7030a77899, 0xbd92393d82c6b2e4),
        (0x397fcfedacb03781, 0x3d131085da82054c),
        (0xb8e45f51f6118e46, 0xbc8f4ed4b492ebcc),
        (0xb89bd46046c3c8de, 0x3c04b7ac8a1b15d0),
        (0x37d1a206fb205e32, 0xbb769201941d0d49),
        (0x3782f38acbf23993, 0x3ae4987e587ab039),
        (0x36b691bdabf5672b, 0xba4ff1953e0a7c5b),
        (0x3636e1c8cd260e18, 0x39b55031dc5e1967),
    ];
    let z = f_polyeval15(
        x2,
        DoubleDouble::from_bit_pair(Z[0]),
        DoubleDouble::from_bit_pair(Z[1]),
        DoubleDouble::from_bit_pair(Z[2]),
        DoubleDouble::from_bit_pair(Z[3]),
        DoubleDouble::from_bit_pair(Z[4]),
        DoubleDouble::from_bit_pair(Z[5]),
        DoubleDouble::from_bit_pair(Z[6]),
        DoubleDouble::from_bit_pair(Z[7]),
        DoubleDouble::from_bit_pair(Z[8]),
        DoubleDouble::from_bit_pair(Z[9]),
        DoubleDouble::from_bit_pair(Z[10]),
        DoubleDouble::from_bit_pair(Z[11]),
        DoubleDouble::from_bit_pair(Z[12]),
        DoubleDouble::from_bit_pair(Z[13]),
        DoubleDouble::from_bit_pair(Z[14]),
    );
    let w_log = log_dd(x);
    DoubleDouble::mul_add(w, w_log, -z).to_f64()
}

/// This method on small range searches for nearest zero or extremum.
/// Then picks stored series expansion at the point end evaluates the poly at the point.
#[inline]
pub(crate) fn y0_small_argument_path(x: f64) -> f64 {
    let x_abs = x;

    // let avg_step = 74.607799 / 47.0;
    // let inv_step = 1.0 / avg_step;

    const INV_STEP: f64 = 0.6299609508652038;

    let fx = x_abs * INV_STEP;
    const Y0_ZEROS_COUNT: f64 = (Y0_ZEROS.len() - 1) as f64;
    let idx0 = fx.min(Y0_ZEROS_COUNT) as usize;
    let idx1 = fx.ceil().min(Y0_ZEROS_COUNT) as usize;

    let found_zero0 = DoubleDouble::from_bit_pair(Y0_ZEROS[idx0]);
    let found_zero1 = DoubleDouble::from_bit_pair(Y0_ZEROS[idx1]);

    let dist0 = (found_zero0.hi - x_abs).abs();
    let dist1 = (found_zero1.hi - x_abs).abs();

    let (found_zero, idx, dist) = if dist0 < dist1 {
        (found_zero0, idx0, dist0)
    } else {
        (found_zero1, idx1, dist1)
    };

    if idx == 0 {
        return j0_maclaurin_series(x);
    }

    let j1c = &Y0_COEFFS[idx - 1];
    let c0 = j1c;

    let r = DoubleDouble::full_add_f64(DoubleDouble::new(-found_zero.lo, -found_zero.hi), x_abs);

    // We hit exact zero, value, better to return it directly
    if dist == 0. {
        return f64::from_bits(Y0_ZEROS_VALUES[idx]);
    }

    let c = &c0[15..];

    let p0 = f_polyeval35(
        r.to_f64(),
        f64::from_bits(c[0].1),
        f64::from_bits(c[1].1),
        f64::from_bits(c[2].1),
        f64::from_bits(c[3].1),
        f64::from_bits(c[4].1),
        f64::from_bits(c[5].1),
        f64::from_bits(c[6].1),
        f64::from_bits(c[7].1),
        f64::from_bits(c[8].1),
        f64::from_bits(c[9].1),
        f64::from_bits(c[10].1),
        f64::from_bits(c[11].1),
        f64::from_bits(c[12].1),
        f64::from_bits(c[13].1),
        f64::from_bits(c[14].1),
        f64::from_bits(c[15].1),
        f64::from_bits(c[16].1),
        f64::from_bits(c[17].1),
        f64::from_bits(c[18].1),
        f64::from_bits(c[19].1),
        f64::from_bits(c[20].1),
        f64::from_bits(c[21].1),
        f64::from_bits(c[22].1),
        f64::from_bits(c[23].1),
        f64::from_bits(c[24].1),
        f64::from_bits(c[25].1),
        f64::from_bits(c[26].1),
        f64::from_bits(c[27].1),
        f64::from_bits(c[28].1),
        f64::from_bits(c[29].1),
        f64::from_bits(c[30].1),
        f64::from_bits(c[31].1),
        f64::from_bits(c[32].1),
        f64::from_bits(c[33].1),
        f64::from_bits(c[34].1),
    );

    let c = c0;

    let mut p_e = DoubleDouble::mul_f64_add(r, p0, DoubleDouble::from_bit_pair(c[14]));
    p_e = DoubleDouble::mul_add(p_e, r, DoubleDouble::from_bit_pair(c[13]));
    p_e = DoubleDouble::mul_add(p_e, r, DoubleDouble::from_bit_pair(c[12]));
    p_e = DoubleDouble::mul_add(p_e, r, DoubleDouble::from_bit_pair(c[11]));
    p_e = DoubleDouble::mul_add(p_e, r, DoubleDouble::from_bit_pair(c[10]));
    p_e = DoubleDouble::mul_add(p_e, r, DoubleDouble::from_bit_pair(c[9]));
    p_e = DoubleDouble::mul_add(p_e, r, DoubleDouble::from_bit_pair(c[8]));
    p_e = DoubleDouble::mul_add(p_e, r, DoubleDouble::from_bit_pair(c[7]));
    p_e = DoubleDouble::mul_add(p_e, r, DoubleDouble::from_bit_pair(c[6]));
    p_e = DoubleDouble::mul_add(p_e, r, DoubleDouble::from_bit_pair(c[5]));
    p_e = DoubleDouble::mul_add(p_e, r, DoubleDouble::from_bit_pair(c[4]));
    p_e = DoubleDouble::mul_add(p_e, r, DoubleDouble::from_bit_pair(c[3]));
    p_e = DoubleDouble::mul_add(p_e, r, DoubleDouble::from_bit_pair(c[2]));
    p_e = DoubleDouble::mul_add(p_e, r, DoubleDouble::from_bit_pair(c[1]));
    p_e = DoubleDouble::mul_add(p_e, r, DoubleDouble::from_bit_pair(c[0]));

    p_e.to_f64()
}

/*
   Evaluates:
   J1 = sqrt(2/(PI*x)) * beta(x) * sin(x - PI/4 - alpha(x))
*/
#[inline]
pub(crate) fn y0_asympt(x: f64) -> f64 {
    const SQRT_2_OVER_PI: DoubleDouble = DoubleDouble::new(
        f64::from_bits(0xbc8cbc0d30ebfd15),
        f64::from_bits(0x3fe9884533d43651),
    );
    const MPI_OVER_4: DoubleDouble = DoubleDouble::new(
        f64::from_bits(0xbc81a62633145c07),
        f64::from_bits(0xbfe921fb54442d18),
    );

    let recip = if x.to_bits() > 0x7fd000000000000u64 {
        DoubleDouble::quick_mult_f64(DoubleDouble::from_exact_safe_div(4.0, x), 0.25)
    } else {
        DoubleDouble::from_recip(x)
    };

    let alpha = crate::bessel::j0::j0_asympt_alpha(recip);
    let beta = crate::bessel::j0::j0_asympt_beta(recip);

    let AngleReduced { angle } = rem2pi_any(x);

    // Without full subtraction cancellation happens sometimes
    let x0pi34 = DoubleDouble::dd_sub(MPI_OVER_4, alpha);
    let r0 = DoubleDouble::dd_add(angle, x0pi34);

    let m_cos = sin_dd_small(r0);
    let z0 = DoubleDouble::quick_mult(beta, m_cos);
    let r_sqrt = DoubleDouble::from_rsqrt(x);
    let scale = DoubleDouble::quick_mult(SQRT_2_OVER_PI, r_sqrt);
    let p = DoubleDouble::quick_mult(scale, z0);
    let norm = DoubleDouble::from_full_exact_add(p.hi, p.lo);
    norm.to_f64()
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_y0_small_argument_path() {
        assert_eq!(f_y0(98.1760435789366), 0.0000000000000056889416242533015);
        assert_eq!(
            f_y0(91.8929453121571802176),
            -0.00000000000000007281665706677893
        );
        assert_eq!(f_y0(80.), -0.05562033908977);
        assert_eq!(f_y0(5.), -0.30851762524903376);
        assert_eq!(
            f_y0(f64::from_bits(0x3fec982eb8d417ea)),
            -0.000000000000000023389279284062102
        );
        assert_eq!(f_y0(f64::from_bits(0x3e04cdee58a47edd)), -13.58605001628649);
        assert_eq!(
            f_y0(0.89357696627916749),
            -0.000000000000000023389279284062102
        );
        assert!(f_y0(f64::NAN).is_nan());
        assert_eq!(f_y0(f64::INFINITY), 0.);
        assert!(f_y0(f64::NEG_INFINITY).is_nan());
    }
}
