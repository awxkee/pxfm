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
use crate::bessel::j0f::{j0f_asympt_alpha, j0f_asympt_beta, j1f_rsqrt};
use crate::bessel::y0f_coeffs::{Y0_ZEROS, Y0_ZEROS_VALUES, Y0F_COEFFS};
use crate::common::f_fmla;
use crate::double_double::DoubleDouble;
use crate::logs::{LOG_COEFFS, LOG_R_DD, LOG_RANGE_REDUCTION};
use crate::polyeval::{f_polyeval4, f_polyeval10, f_polyeval20};
use crate::sin_helper::sin_small;
use crate::sincos_reduce::rem2pif_any;

/// Bessel of the second kind of order 0 (Y0)
///
/// Max ULP 0.5
pub fn f_y0f(x: f32) -> f32 {
    if x < 0. {
        return f32::NAN;
    }

    if (x.to_bits() & 0x0007_ffff) == 0 {
        if x == 0. {
            return f32::NEG_INFINITY;
        }

        if x.is_nan() {
            return x + x;
        }

        if x.is_infinite() {
            if x.is_sign_negative() {
                return f32::NAN;
            }
            return 0.;
        }
    }

    let xb = x.to_bits();

    if xb <= 0x3faccccdu32 {
        // 1.35
        return y0f_near_zero(f32::from_bits(xb));
    }

    if xb <= 0x4296999au32 {
        // 75.3
        return y0f_small_argument_path(f32::from_bits(xb));
    }

    // Exceptions:
    let xb = x.to_bits();
    if xb == 0x5023e87f {
        return f32::from_bits(0x28085b2d);
    } else if xb == 0x48171521 {
        return f32::from_bits(0x2bd244ba);
    } else if xb == 0x4398c299 {
        return f32::from_bits(0x32c730db);
    } else if xb == 0x7f0e5a38 {
        return f32::from_bits(0x131f680b);
    } else if xb == 0x6ef9be45 {
        return f32::from_bits(0x987d8a8f);
    }

    y0f_asympt(x)
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
fn y0f_near_zero(x: f32) -> f32 {
    const W: [u64; 10] = [
        0x3fe45f306dc9c883,
        0xbfc45f306dc9c883,
        0x3f845f306dc9c883,
        0xbf321bb945252402,
        0x3ed21bb945252402,
        0xbe672db9f21b0f5f,
        0x3df49a6c656d62ff,
        0xbd7ae90af76a4d0f,
        0x3cfae90af76a4d0f,
        0xbc754331c053fdad,
    ];
    let dx = x as f64;
    let x2 = dx * dx;
    let w0 = f_polyeval10(
        x2,
        f64::from_bits(W[0]),
        f64::from_bits(W[1]),
        f64::from_bits(W[2]),
        f64::from_bits(W[3]),
        f64::from_bits(W[4]),
        f64::from_bits(W[5]),
        f64::from_bits(W[6]),
        f64::from_bits(W[7]),
        f64::from_bits(W[8]),
        f64::from_bits(W[9]),
    );
    const Z: [u64; 10] = [
        0x3fb2e4d699cbd01f,
        0xbfc6bbcb41034286,
        0x3f9075b1bbf41364,
        0xbf41a6206b7b973d,
        0x3ee3e99794203bbd,
        0xbe7bce4a600d3ea4,
        0x3e0a6ee796b871b6,
        0xbd92393d82c6b2e4,
        0x3d131085da82054c,
        0xbc8f4ed4b492ebcc,
    ];
    let z0 = f_polyeval10(
        x2,
        f64::from_bits(Z[0]),
        f64::from_bits(Z[1]),
        f64::from_bits(Z[2]),
        f64::from_bits(Z[3]),
        f64::from_bits(Z[4]),
        f64::from_bits(Z[5]),
        f64::from_bits(Z[6]),
        f64::from_bits(Z[7]),
        f64::from_bits(Z[8]),
        f64::from_bits(Z[9]),
    );
    let w_log = bessel_fast_log(dx);
    f_fmla(w0, w_log, -z0) as f32
}

#[inline]
pub(crate) fn bessel_fast_log(x: f64) -> f64 {
    let x_u = x.to_bits();

    const E_BIAS: u64 = (1u64 << (11 - 1u64)) - 1u64;

    let mut x_e: i32 = -(E_BIAS as i32);

    // log2(x) = log2(2^x_e * x_m)
    //         = x_e + log2(x_m)
    // Range reduction for log2(x_m):
    // For each x_m, we would like to find r such that:
    //   -2^-8 <= r * x_m - 1 < 2^-7
    let shifted = (x_u >> 45) as i32;
    let index = shifted & 0x7F;
    let r = f64::from_bits(LOG_RANGE_REDUCTION[index as usize]);

    // Add unbiased exponent. Add an extra 1 if the 8 leading fractional bits are
    // all 1's.
    x_e = x_e.wrapping_add(x_u.wrapping_add(1u64 << 45).wrapping_shr(52) as i32);
    let e_x = x_e as f64;

    const LOG_2_HI: f64 = f64::from_bits(0x3fe62e42fefa3800);
    const LOG_2_LO: f64 = f64::from_bits(0x3d2ef35793c76730);

    let log_r_dd = LOG_R_DD[index as usize];

    // hi is exact
    let hi = f_fmla(e_x, LOG_2_HI, f64::from_bits(log_r_dd.1));
    // lo errors ~ e_x * LSB(LOG_2_LO) + LSB(LOG_R[index].lo) + rounding err
    //           <= 2 * (e_x * LSB(LOG_2_LO) + LSB(LOG_R[index].lo))
    let lo = f_fmla(e_x, LOG_2_LO, f64::from_bits(log_r_dd.0));

    // Set m = 1.mantissa.
    let x_m = (x_u & 0x000F_FFFF_FFFF_FFFFu64) | 0x3FF0_0000_0000_0000u64;
    let m = f64::from_bits(x_m);

    let u;
    #[cfg(any(
        all(
            any(target_arch = "x86", target_arch = "x86_64"),
            target_feature = "fma"
        ),
        all(target_arch = "aarch64", target_feature = "neon")
    ))]
    {
        u = f_fmla(r, m, -1.0); // exact
    }
    #[cfg(not(any(
        all(
            any(target_arch = "x86", target_arch = "x86_64"),
            target_feature = "fma"
        ),
        all(target_arch = "aarch64", target_feature = "neon")
    )))]
    {
        use crate::logs::LOG_CD;
        let c_m = x_m & 0x3FFF_E000_0000_0000u64;
        let c = f64::from_bits(c_m);
        u = f_fmla(r, m - c, f64::from_bits(LOG_CD[index as usize])); // exact
    }

    let r1 = DoubleDouble::from_exact_add(hi, u);

    let u_sq = u * u;
    // Degree-7 minimax polynomial
    let p0 = f_fmla(
        u,
        f64::from_bits(LOG_COEFFS[1]),
        f64::from_bits(LOG_COEFFS[0]),
    );
    let p1 = f_fmla(
        u,
        f64::from_bits(LOG_COEFFS[3]),
        f64::from_bits(LOG_COEFFS[2]),
    );
    let p2 = f_fmla(
        u,
        f64::from_bits(LOG_COEFFS[5]),
        f64::from_bits(LOG_COEFFS[4]),
    );
    let p = f_polyeval4(u_sq, lo + r1.lo, p0, p1, p2);
    r1.hi + p
}

/// This method on small range searches for nearest zero or extremum.
/// Then picks stored series expansion at the point end evaluates the poly at the point.
#[inline]
fn y0f_small_argument_path(x: f32) -> f32 {
    let x_abs = x as f64;

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
        // Really should not happen here, but if it is then to log expansion
        return y0f_near_zero(x);
    }

    // We hit exact zero, value, better to return it directly
    if dist == 0. {
        return f64::from_bits(Y0_ZEROS_VALUES[idx]) as f32;
    }

    let c = &Y0F_COEFFS[idx - 1];

    let r = (x_abs - found_zero.hi) - found_zero.lo;

    let p = f_polyeval20(
        r,
        f64::from_bits(c[0]),
        f64::from_bits(c[1]),
        f64::from_bits(c[2]),
        f64::from_bits(c[3]),
        f64::from_bits(c[4]),
        f64::from_bits(c[5]),
        f64::from_bits(c[6]),
        f64::from_bits(c[7]),
        f64::from_bits(c[8]),
        f64::from_bits(c[9]),
        f64::from_bits(c[10]),
        f64::from_bits(c[11]),
        f64::from_bits(c[12]),
        f64::from_bits(c[13]),
        f64::from_bits(c[14]),
        f64::from_bits(c[15]),
        f64::from_bits(c[16]),
        f64::from_bits(c[17]),
        f64::from_bits(c[18]),
        f64::from_bits(c[19]),
    );

    p as f32
}

/*
   Evaluates:
   Y0 = sqrt(2/(PI*x)) * beta(x) * sin(x - PI/4 - alpha(x))
*/
#[inline]
fn y0f_asympt(x: f32) -> f32 {
    let dx = x as f64;

    let alpha = j0f_asympt_alpha(dx);
    let beta = j0f_asympt_beta(dx);

    let angle = rem2pif_any(x);

    const SQRT_2_OVER_PI: f64 = f64::from_bits(0x3fe9884533d43651);
    const MPI_OVER_4: f64 = f64::from_bits(0xbfe921fb54442d18);

    let x0pi34 = MPI_OVER_4 - alpha;
    let r0 = angle + x0pi34;

    let m_cos = sin_small(r0);

    let z0 = beta * m_cos;
    let scale = SQRT_2_OVER_PI * j1f_rsqrt(dx);

    (scale * z0) as f32
}

#[cfg(test)]
mod tests {
    use crate::f_y0f;

    #[test]
    fn test_y0f() {
        assert_eq!(f_y0f(90.5), 0.08254846);
        assert_eq!(f_y0f(77.5), 0.087678276);
        assert_eq!(f_y0f(1.5), 0.3824489);
        assert_eq!(f_y0f(0.5), -0.44451874);
        assert!(f_y0f(-1.).is_nan());
        assert_eq!(f_y0f(0.), f32::NEG_INFINITY);
        assert_eq!(f_y0f(f32::INFINITY), 0.);
        assert!(f_y0f(f32::NEG_INFINITY).is_nan());
    }
}
