/*
 * // Copyright (c) Radzivon Bartoshyk 6/2025. All rights reserved.
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
use crate::dyadic_float::{DyadicFloat128, DyadicSign};
use crate::logs::log_dyadic::{LOG_STEP_2, LOG_STEP_3, LOG_STEP_4};

// The tables keep only the 128-bit significands of -log(1 + s / D), D = 2^bits.
// The reciprocal bounds are -508 <= s <= 257 and 2^16 <= D <= 2^28.
//
// For s > 0, s - s^2/(2D) < D*log(1 + s/D) < s. The lower bound is
// above s - 1, so the binary exponent is floor(log2(s - 1)) - bits.
// For s = 1 the scaled logarithm lies between 1/2 and 1, giving -1 - bits.
//
// For s < 0, write n = -s. The log series gives
// n < -D*log(1 - n/D) < n + n^2/(2(D - n)).
// For n < 256 the increment is below 1. For 256 <= n <= 508 it is below 2,
// keeping the magnitude below 512. The binary exponent is floor(log2(n)) - bits.
//
// The distances from s and n are at least 1/(2(D + 257)) in scaled units,
// by integrating t/(1+t) or t/(1-t). The other bounds leave at least 1/4
// in scaled units at either binade boundary. The unscaled gaps exceed
// 2^-58, so rounding the significand to 128 bits preserves the exponent.
// Inline the reconstruction and addition so the running sum need not cross
// a call boundary at each reduction step.
#[inline(always)]
fn add_log_step(
    sum: DyadicFloat128,
    mantissa: u128,
    s: i32,
    bits: i16,
    zero_sign: DyadicSign,
) -> DyadicFloat128 {
    let (sign, exponent) = if s == 0 {
        (zero_sign, 0)
    } else {
        let (sign, n) = if s > 0 {
            (DyadicSign::Neg, s as u32 - 1)
        } else {
            (DyadicSign::Pos, s.unsigned_abs())
        };
        (sign, 31 - n.leading_zeros() as i16 - bits - 127)
    };
    sum + DyadicFloat128 {
        sign,
        exponent,
        mantissa,
    }
}

// The first reduction provides -2^-8 <= v*2^-60 < 2^-7.
// For bucket i > 0, let n = 128 + i and R = ceil(32640/n). The input
// mantissa m is in [n/128, (n+1)/128), and its reciprocal factor is R/256.
// Then R*n >= 32640 and R*(n+1) < 32641 + n + 32640/n <= 33024,
// since (n-128)*(n-255) <= 0. For i = 0, R = 256 and the bound is direct.
// The product of R <= 256 and the 53-bit significand M fits in i64.
//
// These steps apply r = 1 + s/2^bits to the current input u, yielding
// r*(1 + u) - 1 = u + s/2^bits + u*s/2^bits. Each expression below evaluates
// its three terms exactly in fixed-point arithmetic.
pub(super) fn log_range_reduction(v: i64, sum: DyadicFloat128) -> (DyadicFloat128, DyadicFloat128) {
    // For -2^-29 <= u < 2^-29, all three indices round to k = 0. The
    // reciprocal factors are all 1, and the last table contributes +0.
    if (-(1 << 31)..(1 << 31)).contains(&v) {
        let residual = DyadicFloat128::new(
            if v < 0 {
                DyadicSign::Neg
            } else {
                DyadicSign::Pos
            },
            -125,
            u128::from(v.unsigned_abs()) << 65,
        );
        let sum = if sum.mantissa == 0 {
            DyadicFloat128::zero()
        } else {
            sum
        };
        return (residual, sum);
    }

    // Step 2 rounds u*2^14 to k in [-64, 128], so |u - k*2^-14| <= 2^-15.
    // Rounding the reciprocal contributes at most (1 + k*2^-14)*2^-17
    // to the residual. Thus |vv2| in real units is at most
    // (1 + 257*2^-16)*2^-15 + (1 + 2^-7)*2^-17 < 80.5*2^-21.
    let k2 = ((v + (1 << 45)) >> 46) as i32;
    let idx2 = (k2 + 64) as usize;
    // The scaled reciprocal is 2^16 - 4*k + c, where c = 4*k^2/(2^14 + k).
    // Its derivative has the sign of k. Comparing 8*k^2 with
    // (2*j+1)*(2^14+k) at adjacent integers locates the half-integer
    // crossings of c: ±45/±46 for 1/2, 78/79 for 3/2, 101/102 for 5/2,
    // and 120/121 for 7/2. At -64, c < 3/2, and at 128, c < 9/2.
    let s2 = -4 * k2
        + (k2.unsigned_abs() > 45) as i32
        + (k2 > 78) as i32
        + (k2 > 101) as i32
        + (k2 > 120) as i32;
    let sum = add_log_step(sum, LOG_STEP_2[idx2], s2, 16, DyadicSign::Neg);
    // In integer units |vv2| + |s2*v| < 80.5*2^55 + 508*2^53 < 2^63,
    // so the product, the centered term, and their sum all fit in i64.
    let vv2 = (((s2 as i64) << 44) + v) * (1 << 16) + s2 as i64 * v; // 2^-76 units

    // Step 3 rounds to k in [-80, 80]. For D = 2^21, the exact scaled
    // reciprocal is D^2/(D + k) = D - k + k^2/(D + k).
    // The last term is below 1/2 even at |k| = 80, so s = -k.
    let k3 = ((vv2 + (1 << 54)) >> 55) as i32;
    let sum = add_log_step(
        sum,
        LOG_STEP_3[(k3 + 80) as usize],
        -k3,
        21,
        DyadicSign::Pos,
    );
    let s3 = -(k3 as i64);
    let vv3 = ((((s3 << 55) + vv2) as i128) << 21) + s3 as i128 * vv2 as i128; // 2^-97 units

    // In Step 3, u and k have the same sign unless k = 0. The centered
    // remainder is in [-2^-22, 2^-22), and 0 <= k*u*2^-21 < 6460*2^-42.
    // Thus the reduced argument is strictly above -64.5*2^-28 and below
    // 64.5*2^-28. Step 4 only needs k in [-64, 64], with s = -k by the
    // same reciprocal formula for D = 2^28.
    let k4 = (((vv3 >> 68) as i32) + 1) >> 1;
    let sum = add_log_step(
        sum,
        LOG_STEP_4[(k4 + 64) as usize],
        -k4,
        28,
        DyadicSign::Pos,
    );
    let s4 = -(k4 as i128);
    // The same centered-remainder bound gives
    // -2^-29 - 4128*2^-56 < vv4 in real units < 2^-29.
    let vv4 = (((s4 << 69) + vv3) << 28) + s4 * vv3; // 2^-125 units

    let residual = DyadicFloat128::new(
        if vv4 < 0 {
            DyadicSign::Neg
        } else {
            DyadicSign::Pos
        },
        -125,
        vv4.unsigned_abs(),
    );
    (residual, sum)
}
