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
use crate::common::*;
use crate::compoundf::{
    COMPOUNDF_EXP2_T, COMPOUNDF_EXP2_U, compoundf_exp2_poly2, compoundf_log2p1_accurate,
    compoundf_log2p1_fast,
};
use crate::double_double::DoubleDouble;
use crate::exponents::exp2m1_accurate_tiny;
use std::hint::black_box;

// INVLOG2 = 1/log(2) * (1 + eps1) with |eps1| < 2^-55.976
const INVLOG2: f64 = f64::from_bits(0x3ff71547652b82fe);

#[cold]
#[inline(never)]
fn as_compoundm1f_special(x: f32, y: f32) -> f32 {
    let nx = x.to_bits();
    let ny = y.to_bits();
    let ax: u32 = nx.wrapping_shl(1);
    let ay: u32 = ny.wrapping_shl(1);

    if ax == 0 || ay == 0 {
        // x or y is 0
        if ax == 0 {
            // compound(0,y) = 1 except for y = sNaN
            return if y.is_nan() { x + y } else { 0.0 };
        }

        if ay == 0 {
            // compound (x, 0)
            if x.is_nan() {
                return x + y;
            } // x = sNaN
            return if x < -1.0 {
                f32::NAN // rule (g)
            } else {
                0.0
            }; // rule (a)
        }
    }

    let mone = (-1.0f32).to_bits();
    if ay >= 0xffu32 << 24 {
        // y=Inf/NaN
        // the case x=0 was already checked above
        if ax > 0xffu32 << 24 {
            return x + y;
        } // x=NaN
        if ay == 0xffu32 << 24 {
            // y = +/-Inf
            if nx > mone {
                return f32::NAN;
            } // rule (g)
            let sy = ny >> 31; // sign bit of y
            if nx == mone {
                return if sy == 0 {
                    -1. // Rule (c)
                } else {
                    f32::INFINITY // Rule (b)
                };
            }
            if x < 0.0 {
                return if sy == 0 { -1. } else { f32::INFINITY };
            }
            if x > 0.0 {
                return if sy != 0 { -1. } else { f32::INFINITY };
            }
            return 0.0;
        }
        return x + y; // case y=NaN
    }

    if nx >= mone || nx >= 0xffu32 << 23 {
        // x is Inf, NaN or <= -1
        if ax == 0xffu32 << 24 {
            // x is +Inf or -Inf
            if (nx >> 31) != 0 {
                return f32::NAN;
            } // x = -Inf, rule (g)
            // (1 + Inf)^y = +Inf for y > 0, +0 for y < 0
            return (if (ny >> 31) != 0 { 1.0 / x } else { x }) - 1.;
        }
        if ax > 0xffu32 << 24 {
            return x + y;
        } // x is NaN
        if nx > mone {
            return f32::NAN; // x < -1.0: rule (g)
        }
        // now x = -1
        return if (ny >> 31) != 0 {
            // y < 0
            f32::INFINITY
        } else {
            // y > 0
            -1.0
        };
    }
    -1.
}

/* for |z| <= 2^-6, returns an approximation of 2^z
with absolute error < 2^-43.540  */
#[inline]
fn compoundf_expf_poly(z: f64) -> f64 {
    /* Q is a degree-4 polynomial generated by Sollya (cf q1.sollya)
    with absolute error < 2^-43.549 */
    const Q: [u64; 5] = [
        0x3fe62e42fefa39ef,
        0x3fcebfbdff8098eb,
        0x3fac6b08d7045dc3,
        0x3f83b2b276ce985d,
        0x3f55d8849c67ace4,
    ];
    let z2 = z * z;
    let c3 = dd_fmla(f64::from_bits(Q[4]), z, f64::from_bits(Q[3]));
    let c0 = dd_fmla(f64::from_bits(Q[1]), z, f64::from_bits(Q[0]));
    let c2 = dd_fmla(c3, z, f64::from_bits(Q[2]));
    dd_fmla(c2, z2, c0) * z
}

/* return the correct rounding of (1+x)^y, otherwise -1.0
where t is an approximation of y*log2(1+x) with absolute error < 2^-40.680,
assuming 0x1.7154759a0df53p-24 <= |t| <= 150
exact is non-zero iff (1+x)^y is exact or midpoint */
fn exp2m1_fast(t: f64) -> f64 {
    let k = t.round_ties_even(); // 0 <= |k| <= 150
    let mut r = t - k; // |r| <= 1/2, exact
    let mut v: u64 = (3.015625 + r).to_bits(); // 2.5 <= v <= 3.5015625
    // we add 2^-6 so that i is rounded to nearest
    let i: i32 = (v >> 46) as i32 - 0x10010; // 0 <= i <= 32
    r -= f64::from_bits(COMPOUNDF_EXP2_T[i as usize]); // exact
    // now |r| <= 2^-6
    // 2^t = 2^k * exp2_U[i][0] * 2^r
    let mut s = f64::from_bits(COMPOUNDF_EXP2_U[i as usize].1);
    let su = ((k as u64).wrapping_add(0x3ffu64)) << 52;
    s *= f64::from_bits(su);
    let q_poly = compoundf_expf_poly(r);
    v = q_poly.to_bits();
    /* the absolute error on exp2_U[i][0] is bounded by 2^-53.092, with
    exp2_U[i][0] < 2^0.5, and that on q1(r) is bounded by 2^-43.540,
    with |q1(r)| < 1.011, thus |v| < 1.43, and the absolute error on v is
    bounded by ulp(v) + 2^0.5s * 2^-43.540 + 2^-53.092 * 1.011 < 2^-43.035.
    Now t approximates u := y*log2(1+x) with |t-u| < 2^-40.680 thus
    2^u = 2^t * (1 + eps) with eps < 2^(2^-40.680)-1 < 2^-41.208.
    The total absolute error is thus bounded by 2^-43.035 + 2^-41.208
    < 2^-40.849. */
    let mut err: u64 = 0x3d61d00000000000; // 2^-40.849 < 0x1.1dp-41

    #[cfg(any(
        all(
            any(target_arch = "x86", target_arch = "x86_64"),
            target_feature = "fma"
        ),
        all(target_arch = "aarch64", target_feature = "neon")
    ))]
    {
        v = f_fmla(f64::from_bits(v), s, s - 1f64).to_bits();
    }
    #[cfg(not(any(
        all(
            any(target_arch = "x86", target_arch = "x86_64"),
            target_feature = "fma"
        ),
        all(target_arch = "aarch64", target_feature = "neon")
    )))]
    {
        let p0 = DoubleDouble::from_full_exact_add(s, -1.);
        let z = DoubleDouble::from_exact_mult(f64::from_bits(v), s);
        v = DoubleDouble::add(z, p0).to_f64().to_bits();
    }

    // in case of potential underflow, we defer to the accurate path
    if f64::from_bits(v) < f64::from_bits(0x3d61d00000000000) {
        return -1.0;
    }
    err = err.wrapping_add(((k as i64) << 52) as u64); // scale the error by 2^k too
    let lb = (f64::from_bits(v) - f64::from_bits(err)) as f32;
    let rb = (f64::from_bits(v) + f64::from_bits(err)) as f32;
    if lb != rb {
        return -1.0;
    } // rounding test failed

    f64::from_bits(v)
}
fn compoundf_exp2m1_accurate(x_dd: DoubleDouble, x: f32, y: f32) -> f32 {
    if y == 1.0 {
        let res = x;
        return res;
    }

    // check easy cases h+l is tiny thus 2^(h+l) rounds to 1, 1- or 1+
    // if x_dd.hi.abs() <= f64::from_bits(0x3fc0000000000000u64) {
    //     /* the relative error between h and y*log2(1+x) is bounded by
    //     (1 + 2^-48.445) * (1 + 2^-91.120) - 1 < 2^-48.444.
    //     2^h rounds to 1 to nearest for |h| <= H0 := 0x1.715476af0d4d9p-25.
    //     The above threshold is such that h*(1+2^-48.444) < H0. */
    //     return  exp2m1_accurate_tiny(x_dd.to_f64()) as f32;
    // }

    let k = x_dd.hi.round_ties_even(); // |k| <= 150

    // check easy cases h+l is tiny thus 2^(h+l) rounds to 1, 1- or 1+
    if k == 0. && x_dd.hi.abs() <= f64::from_bits(0x3e6715476af0d4c8) {
        /* the relative error between h and y*log2(1+x) is bounded by
        (1 + 2^-48.445) * (1 + 2^-91.120) - 1 < 2^-48.444.
        2^h rounds to 1 to nearest for |h| <= H0 := 0x1.715476af0d4d9p-25.
        The above threshold is such that h*(1+2^-48.444) < H0. */
        // let z0 = 1.0 + x_dd.hi * 0.5;
        // let k = Dekker::from_exact_sub(z0, 1.);
        // return k.to_f64() as f32;

        return exp2m1_accurate_tiny(x_dd.to_f64()) as f32;
    }

    let r = x_dd.hi - k; // |r| <= 1/2, exact
    // since r is an integer multiple of ulp(h), fast_two_sum() below is exact
    let mut v_dd = DoubleDouble::from_exact_add(r, x_dd.lo);
    let mut v = (3.015625 + v_dd.hi).to_bits(); // 2.5 <= v <= 3.5015625
    // we add 2^-6 so that i is rounded to nearest
    let i: i32 = ((v >> 46) as i32).wrapping_sub(0x10010); // 0 <= i <= 32
    // h is near (i-16)/2^5
    v_dd.hi -= f64::from_bits(COMPOUNDF_EXP2_T[i as usize]); // exact

    // now |h| <= 2^-6
    // 2^(h+l) = 2^k * exp2_U[i] * 2^(h+l)
    v_dd = DoubleDouble::from_exact_add(v_dd.hi, v_dd.lo);
    let q = compoundf_exp2_poly2(v_dd);

    /* we have 0.989 < qh < 1.011, |ql| < 2^-51.959, and
    |qh + ql - 2^(h+l)| < 2^-85.210 */
    let exp2u = DoubleDouble::from_bit_pair(COMPOUNDF_EXP2_U[i as usize]);
    let mut q = DoubleDouble::quick_mult(exp2u, q);

    q = DoubleDouble::from_exact_add(q.hi, q.lo);

    let mut du = (k as i64).wrapping_add(0x3ff).wrapping_shl(52) as u64;
    du = f64::from_bits(du).to_bits();
    let scale = f64::from_bits(du);

    q.hi *= scale;
    q.lo *= scale;

    let zf: DoubleDouble = if x >= 0. {
        // implies h >= 1 and the fast_two_sum pre-condition holds
        DoubleDouble::from_exact_add(q.hi, -1.0)
    } else {
        DoubleDouble::from_exact_add(-1.0, q.hi)
    };
    q.lo += zf.lo;
    q.hi = zf.hi;

    v = q.to_f64().to_bits();

    f64::from_bits(v) as f32
}

// at input, exact is non-zero iff (1+x)^y is exact
// x,y=0x1.0f6f1ap+1,0x1.c643bp+5: 49 identical bits after round bit
// x,y=0x1.ef272cp+15,-0x1.746ab2p+1: 55 identical bits after round bit
// x,y=0x1.07ffcp+0,-0x1.921a8ap+4: 47 identical bits after round bit
#[cold]
#[inline(never)]
fn compoundm1f_accurate(x: f32, y: f32) -> f32 {
    let mut v = compoundf_log2p1_accurate(x as f64);
    v = DoubleDouble::quick_mult_f64(v, y as f64);
    compoundf_exp2m1_accurate(v, x, y)
}

/// Computes compound (1.0 + x)^y - 1
///
/// Max ULP 0.5
#[inline]
pub fn f_compound_m1f(x: f32, y: f32) -> f32 {
    /* Rules from IEEE 754-2019 for compound (x, n) with n integer:
       (a) compound (x, 0) is 1 for x >= -1 or quiet NaN
       (b) compound (-1, n) is +Inf and signals the divideByZero exception for n < 0
       (c) compound (-1, n) is +0 for n > 0
       (d) compound (+/-0, n) is 1
       (e) compound (+Inf, n) is +Inf for n > 0
       (f) compound (+Inf, n) is +0 for n < 0
       (g) compound (x, n) is qNaN and signals the invalid exception for x < -1
       (h) compound (qNaN, n) is qNaN for n <> 0.
    */
    let mone = (-1.0f32).to_bits();
    let nx = x.to_bits();
    let ny = y.to_bits();
    if nx >= mone {
        return as_compoundm1f_special(x, y);
    } // x <= -1 
    // now x > -1

    let ax: u32 = nx.wrapping_shl(1);
    let ay: u32 = ny.wrapping_shl(1);

    if ax == 0 || ax >= 0xffu32 << 24 || ay == 0 || ay >= 0xffu32 << 24 {
        return as_compoundm1f_special(x, y);
    } // x=+-0 || x=+-inf/nan || y=+-0 || y=+-inf/nan

    // evaluate (1+x)^y explicitly for integer y in [-16,16] range and |x|<2^64
    if y.floor() == y && ay <= 0x83000000u32 && ax <= 0xbefffffeu32 {
        if ax <= 0x62000000u32 {
            return 1.0 + y * x;
        } // does it work for |x|<2^-29 and |y|<=16?
        let ky: i32 = (((ay & 0x00ffffff) | 1 << 24) >> (151 - (ay >> 24))) as i32;
        let s = 1.0 + x as f64;
        let mut p = 1.;
        let s2 = s * s;
        let s4 = s2 * s2;
        let s8 = s4 * s4;
        let s16 = s8 * s8;
        let sn: [f64; 6] = [1., s, s2, s4, s8, s16];
        p *= sn[(ky & 1) as usize];
        p *= sn[(ky & 2) as usize];
        p *= sn[(((ky >> 2) & 1) * 3) as usize];
        p *= sn[((ky >> 1) & 4) as usize];
        p *= sn[(((ky >> 4) & 1) * 5) as usize];
        let z = if (ny >> 31) != 0 { 1. / p } else { p };
        let k = DoubleDouble::from_full_exact_add(z, -1.).to_f64();
        return k as f32;
    }

    let xd = x as f64;
    let yd = y as f64;
    let tx = xd.to_bits();
    let ty = yd.to_bits();

    let l: f64 = if ax < 0x62000000u32 {
        // |x| < 2^-29
        /* |log2(1+x) - 1/log(2) * (x - x^2/2)| < 2^-59.584 * |log2(1+x)|
        (cf compoundf.sollya) */
        let t = xd - (xd * xd) * 0.5;
        /* since x is epresentable in binary32, x*x is exact, and so is (x * x) * 0.5.
           Thus the only error in the computation of t is the final rounding, which
           is bounded by ulp(t): t = (x - x^2/2) * (1 + eps2) with |eps2| < 2^-52
        */
        INVLOG2 * t
        /* since INVLOG2 = 1/log(2) * (1 + eps1) and
        and   t = (x - x^2/2) * (1 + eps2)
        let u = o(INVLOG2 * t) then u = INVLOG2 * t * (1 + eps3) with |eps3|<2^-53
        thus u = 1/log(2) * (x - x^2/2) * (1 + eps1)*(1 + eps2)*(1 + eps3)
        = 1/log(2) * (x - x^2/2) * (1 + eps4) with |eps4| < 2^-50.954
        Now Sollya says the relative error by approximating log2(1+x) by
        1/log(2) * (x - x^2/2) for |x| < 2^-29 is bounded by 2^-59.584
        (file compoundf.sollya), thus:
        u = log2(1+x) * (1+eps4)*(1+eps5) with |eps5| < 2^-59.584
        = log2(1+x) * (1+eps6) with |eps6| < 2^-50.950 */
    } else {
        compoundf_log2p1_fast(f64::from_bits(tx))
    };

    /* l approximates log2(1+x) with relative error < 2^-47.997,
    and 2^-149 <= |l| < 128 */

    let t: u64 = (l * f64::from_bits(ty)).to_bits();
    /* since 2^-149 <= |l| < 128 and 2^-149 <= |y| < 2^128, we have
    2^-298 <= |t| < 2^135, thus no underflow/overflow in double is possible.
    The relative error is bounded by (1+2^-47.997)*(1+2^-52)-1 < 2^-47.909 */

    // detect overflow/underflow
    if (t.wrapping_shl(1)) >= (0x406u64 << 53) {
        // |t| >= 128
        if t >= 0x3018bu64 << 46 {
            // t <= -150
            return black_box(f32::from_bits(0x00800000)) * black_box(f32::from_bits(0x00800000));
        } else if (t >> 63) == 0 {
            // t >= 128: overflow
            return black_box(f32::from_bits(0x7e800000)) * black_box(f32::from_bits(0x7e800000));
        }
    }

    /* since |t| < 150, the absolute error on t is bounded by
    150*2^-47.909 < 2^-40.680 */

    // 2^t rounds to 1 to nearest when |t| <= 0x1.715476ba97f14p-25
    if (t.wrapping_shl(1)) <= 0x3e6715476ba97f14u64 {
        return if (t >> 63) != 0 {
            black_box(1.0) - black_box(f32::from_bits(0x33000000))
        } else {
            black_box(1.0) + black_box(f32::from_bits(0x33000000))
        };
    }

    let res = exp2m1_fast(f64::from_bits(t));
    if res != -1.0 {
        return res as f32;
    }
    compoundm1f_accurate(x, y)
}

#[cfg(test)]
mod tests {
    use crate::compound_m1f::{compoundf_exp2m1_accurate, exp2m1_fast};
    use crate::double_double::DoubleDouble;
    use crate::f_compound_m1f;

    #[test]
    fn test_compoundf() {
        assert_eq!(
            f_compound_m1f(-0.000000000000001191123, -0.000000000000001191123),
            0.0000000000000000000000000000014187741
        );
        assert_eq!(f_compound_m1f(-0.000000000000001191123, 16.), 1.0);
        assert_eq!(f_compound_m1f(0.91123, 16.), 31695.21);
    }

    #[test]
    fn test_compoundf_expm1_fast() {
        assert_eq!(exp2m1_fast(3.764), 12.585539943149435);
    }

    #[test]
    fn test_compoundf_expm1_accurate() {
        assert_eq!(
            compoundf_exp2m1_accurate(DoubleDouble::new(0., 2.74), 12., 53.),
            5.680703,
        );
    }
}
