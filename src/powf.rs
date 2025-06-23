/*
 * // Copyright (c) Radzivon Bartoshyk 4/2025. All rights reserved.
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
use crate::expf::expf;
use crate::log2f::LOG2_R;
use crate::logf::{EXP_MASK_F32, f_polyeval3, logf};
use crate::pow::EXP2_MID1;

/// Power function for given value
#[inline]
pub const fn powf(d: f32, n: f32) -> f32 {
    let value = d.abs();
    let c = expf(n * logf(value));
    if n == 1. {
        return d;
    }
    if d < 0.0 {
        let y = n as i32;
        if y % 2 == 0 { c } else { -c }
    } else {
        c
    }
}

#[inline]
const fn is_integer(x: f32) -> bool {
    let x_u = x.to_bits();
    let x_e = (x_u & EXP_MASK_F32) >> 23;
    let lsb = (x_u | EXP_MASK_F32).trailing_zeros();
    const E_BIAS: u32 = (1u32 << (8 - 1u32)) - 1u32;
    const UNIT_EXPONENT: u32 = E_BIAS + 23;
    x_e + lsb >= UNIT_EXPONENT
}

#[inline]
fn is_odd_integer(x: f32) -> bool {
    let x_u = x.to_bits();
    let x_e = (x_u & EXP_MASK_F32) >> 23;
    let lsb = (x_u | EXP_MASK_F32).trailing_zeros();
    const E_BIAS: u32 = (1u32 << (8 - 1u32)) - 1u32;

    const UNIT_EXPONENT: u32 = E_BIAS + 23;
    x_e + lsb == UNIT_EXPONENT
}

static LOG2_R2_DD: [(u64, u64); 128] = [
    (0x0000000000000000, 0x0000000000000000),
    (0xbe6177c23362928c, 0x3f872c8000000000),
    (0xbe9179e0caa9c9ab, 0x3f97440000000000),
    (0xbe8c6cea541f5b70, 0x3fa184c000000000),
    (0xbe966c4d4e554434, 0x3fa773a000000000),
    (0xbe770700a00fdd55, 0x3fad6ec000000000),
    (0x3e853002a4e86631, 0x3fb1bb3000000000),
    (0x3e6fcd15f101c142, 0x3fb4c56000000000),
    (0x3e925b3eed319ced, 0x3fb7d60000000000),
    (0xbe94195120d8486f, 0x3fb960d000000000),
    (0x3e845b878e27d0d9, 0x3fbc7b5000000000),
    (0x3e9770744593a4cb, 0x3fbf9c9000000000),
    (0x3e9c673032495d24, 0x3fc097e000000000),
    (0xbe91eaa65b49696e, 0x3fc22db000000000),
    (0x3e9b2866f2850b22, 0x3fc3c6f800000000),
    (0x3e68ee37cd2ea9d3, 0x3fc494f800000000),
    (0x3e77e86f9c2154fb, 0x3fc633a800000000),
    (0x3e58e3cfc25f0ce6, 0x3fc7046000000000),
    (0x3e357f7a64ccd537, 0x3fc8a89800000000),
    (0xbe9a761c09fbd2ae, 0x3fc97c2000000000),
    (0x3e924bea9a2c66f3, 0x3fcb260000000000),
    (0xbe660002ccfe43f5, 0x3fcbfc6800000000),
    (0x3e969f220e97f22c, 0x3fcdac2000000000),
    (0xbe96164f64c210e0, 0x3fce858000000000),
    (0xbe70c1678ae89767, 0x3fd01d9c00000000),
    (0xbe9f26a05c813d57, 0x3fd08bd000000000),
    (0x3e74d8fc561c8d44, 0x3fd169c000000000),
    (0xbe9362ad8f7ca2d0, 0x3fd1d98400000000),
    (0x3e92b13cd6c4d042, 0x3fd249cc00000000),
    (0xbe91c8f11979a5db, 0x3fd32c0000000000),
    (0x3e8c2ab3edefe569, 0x3fd39de800000000),
    (0x3e57c3eca28e69ca, 0x3fd4106000000000),
    (0xbe734c4e99e1c6c6, 0x3fd4f6fc00000000),
    (0xbe9194a871b63619, 0x3fd56b2400000000),
    (0x3e8e3dd5c1c885ae, 0x3fd5dfdc00000000),
    (0xbe86ccf3b1129b7c, 0x3fd6552c00000000),
    (0xbe82f346e2bf924b, 0x3fd6cb1000000000),
    (0xbe8fa61aaa59c1d8, 0x3fd7b8a000000000),
    (0x3e990c11fd32a3ab, 0x3fd8304c00000000),
    (0x3e457f7a64ccd537, 0x3fd8a89800000000),
    (0x3e4249ba76fee235, 0x3fd9218000000000),
    (0xbe8aad2729b21ae5, 0x3fd99b0800000000),
    (0x3e971810a5e18180, 0x3fda8ff800000000),
    (0xbe46172fe015e13c, 0x3fdb0b6800000000),
    (0x3e75ec6c1bfbf89a, 0x3fdb877c00000000),
    (0x3e7678bf6cdedf51, 0x3fdc043800000000),
    (0x3e9c2d45fe43895e, 0x3fdc819c00000000),
    (0xbe99ee52ed49d71d, 0x3fdcffb000000000),
    (0x3e45786af187a96b, 0x3fdd7e6c00000000),
    (0x3e83ab0dc56138c9, 0x3fddfdd800000000),
    (0x3e9fe538ab34efb5, 0x3fde7df400000000),
    (0xbe9e4fee07aa4b68, 0x3fdefec800000000),
    (0xbe9172f32fe67287, 0x3fdf804c00000000),
    (0xbe99a83ff9ab9cc8, 0x3fe0014400000000),
    (0xbe968cb06cece193, 0x3fe042be00000000),
    (0x3e98cd71ddf82e20, 0x3fe0849400000000),
    (0x3e95e18ab2df3ae6, 0x3fe0c6ca00000000),
    (0x3e65dee4d9d8a273, 0x3fe1096000000000),
    (0x3e5fcd15f101c142, 0x3fe14c5600000000),
    (0xbe82474b0f992ba1, 0x3fe18fae00000000),
    (0x3e74b5a92a606047, 0x3fe1d36800000000),
    (0x3e916186fcf54bbd, 0x3fe2178600000000),
    (0x3e418efabeb7d722, 0x3fe25c0a00000000),
    (0xbe7e5fc7d238691d, 0x3fe2a0f400000000),
    (0x3e9f5809faf6283c, 0x3fe2e64400000000),
    (0x3e9f5809faf6283c, 0x3fe2e64400000000),
    (0x3e9c6e1dcd0cb449, 0x3fe32bfe00000000),
    (0x3e976e0e8f74b4d5, 0x3fe3722200000000),
    (0xbe7cb82c89692d99, 0x3fe3b8b200000000),
    (0xbe963161c5432aeb, 0x3fe3ffae00000000),
    (0x3e9458104c41b901, 0x3fe4471600000000),
    (0x3e9458104c41b901, 0x3fe4471600000000),
    (0xbe9cd9d0cde578d5, 0x3fe48ef000000000),
    (0x3e5b9884591add87, 0x3fe4d73800000000),
    (0x3e9c6042978605ff, 0x3fe51ff200000000),
    (0xbe9fc4c96b37dcf6, 0x3fe5692200000000),
    (0xbe72f346e2bf924b, 0x3fe5b2c400000000),
    (0xbe72f346e2bf924b, 0x3fe5b2c400000000),
    (0x3e9c4e4fbb68a4d1, 0x3fe5fcdc00000000),
    (0xbe89d499bd9b3226, 0x3fe6476e00000000),
    (0xbe8f89b355ede26f, 0x3fe6927800000000),
    (0xbe8f89b355ede26f, 0x3fe6927800000000),
    (0x3e753c7e319f6e92, 0x3fe6ddfc00000000),
    (0xbe9b291f070528c7, 0x3fe729fe00000000),
    (0x3e62967a451a7b48, 0x3fe7767c00000000),
    (0x3e62967a451a7b48, 0x3fe7767c00000000),
    (0x3e9244fcff690fce, 0x3fe7c37a00000000),
    (0x3e846fd97f5dc572, 0x3fe810fa00000000),
    (0x3e846fd97f5dc572, 0x3fe810fa00000000),
    (0xbe9f3a7352663e50, 0x3fe85efe00000000),
    (0x3e8b3cda690370b5, 0x3fe8ad8400000000),
    (0x3e8b3cda690370b5, 0x3fe8ad8400000000),
    (0x3e83226b211bf1d9, 0x3fe8fc9200000000),
    (0x3e8d24b136c101ee, 0x3fe94c2800000000),
    (0x3e8d24b136c101ee, 0x3fe94c2800000000),
    (0x3e97c40c7907e82a, 0x3fe99c4800000000),
    (0xbe9e81781d97ee91, 0x3fe9ecf600000000),
    (0xbe9e81781d97ee91, 0x3fe9ecf600000000),
    (0xbe96a77813f94e01, 0x3fea3e3000000000),
    (0xbe91cfdeb43cfd00, 0x3fea8ffa00000000),
    (0xbe91cfdeb43cfd00, 0x3fea8ffa00000000),
    (0xbe8f983f74d3138f, 0x3feae25600000000),
    (0xbe8e278ae1a1f51f, 0x3feb354600000000),
    (0xbe8e278ae1a1f51f, 0x3feb354600000000),
    (0xbe897552b7b5ea45, 0x3feb88cc00000000),
    (0xbe897552b7b5ea45, 0x3feb88cc00000000),
    (0xbe719b4f3c72c4f8, 0x3febdcea00000000),
    (0x3e8f7402d26f1a12, 0x3fec31a200000000),
    (0x3e8f7402d26f1a12, 0x3fec31a200000000),
    (0xbe82056d5dd31d96, 0x3fec86f800000000),
    (0xbe82056d5dd31d96, 0x3fec86f800000000),
    (0xbe76e46335aae723, 0x3fecdcec00000000),
    (0xbe9beb244c59f331, 0x3fed338200000000),
    (0xbe9beb244c59f331, 0x3fed338200000000),
    (0x3e416c071e93fd97, 0x3fed8aba00000000),
    (0x3e416c071e93fd97, 0x3fed8aba00000000),
    (0x3e9d8175819530c2, 0x3fede29800000000),
    (0x3e9d8175819530c2, 0x3fede29800000000),
    (0x3e851bd552842c1c, 0x3fee3b2000000000),
    (0x3e851bd552842c1c, 0x3fee3b2000000000),
    (0x3e9914e204f19d94, 0x3fee945200000000),
    (0x3e9914e204f19d94, 0x3fee945200000000),
    (0x3e9c55d997da24fd, 0x3feeee3200000000),
    (0x3e9c55d997da24fd, 0x3feeee3200000000),
    (0xbe9685c2d2298a6e, 0x3fef48c400000000),
    (0xbe9685c2d2298a6e, 0x3fef48c400000000),
    (0x3e97a4887bd74039, 0x3fefa40600000000),
    (0x0000000000000000, 0x3ff0000000000000),
];

/// Power function for given value using FMA
///
/// Max found ULP 0.5
#[inline]
pub fn f_powf(x: f32, y: f32) -> f32 {
    let mut x_u = x.to_bits();
    let x_abs = x_u & 0x7fff_ffff;
    let mut y = y;
    let y_u = y.to_bits();
    let y_abs = y_u & 0x7fff_ffff;
    let mut x = x;

    if (y_abs & 0x0007_ffff == 0) || (y_abs > 0x4f170000) {
        // y is signaling NaN
        if x.is_nan() || y.is_nan() {
            return f32::NAN;
        }

        // Exceptional exponents.
        if y == 0.0 {
            return 1.0;
        }

        match y_abs {
            0x7f80_0000 => {
                if x_abs > 0x7f80_0000 {
                    // pow(NaN, +-Inf) = NaN
                    return x;
                }
                if x_abs == 0x3f80_0000 {
                    // pow(+-1, +-Inf) = 1.0f
                    return 1.0;
                }
                if x == 0.0 && y_u == 0xff80_0000 {
                    // pow(+-0, -Inf) = +inf and raise FE_DIVBYZERO
                    return f32::INFINITY;
                }
                // pow (|x| < 1, -inf) = +inf
                // pow (|x| < 1, +inf) = 0.0f
                // pow (|x| > 1, -inf) = 0.0f
                // pow (|x| > 1, +inf) = +inf
                return if (x_abs < 0x3f80_0000) == (y_u == 0xff80_0000) {
                    f32::INFINITY
                } else {
                    0.
                };
            }
            _ => {
                match y_u {
                    0x3f00_0000 => {
                        // pow(x, 1/2) = sqrt(x)
                        if x == 0.0 || x_u == 0xff80_0000 {
                            // pow(-0, 1/2) = +0
                            // pow(-inf, 1/2) = +inf
                            // Make sure it is correct for FTZ/DAZ.
                            return x * x;
                        }
                        let r = x.sqrt();
                        return if r.to_bits() != 0x8000_0000 { r } else { 0.0 };
                    }
                    0x3f80_0000 => {
                        return x;
                    } // y = 1.0f
                    0x4000_0000 => return x * x, // y = 2.0f
                    _ => {
                        let is_int = is_integer(y);
                        if is_int && (y_u > 0x4000_0000) && (y_u <= 0x41c0_0000) {
                            // Check for exact cases when 2 < y < 25 and y is an integer.
                            let mut msb: i32 = if x_abs == 0 {
                                32 - 2
                            } else {
                                x_abs.leading_zeros() as i32
                            };
                            msb = if msb > 8 { msb } else { 8 };
                            let mut lsb: i32 = if x_abs == 0 {
                                0
                            } else {
                                x_abs.trailing_zeros() as i32
                            };
                            lsb = if lsb > 23 { 23 } else { lsb };
                            let extra_bits: i32 = 32 - 2 - lsb - msb;
                            let iter = y as i32;

                            if extra_bits * iter <= 23 + 2 {
                                // The result is either exact or exactly half-way.
                                // But it is exactly representable in double precision.
                                let x_d = x as f64;
                                let mut result = x_d;
                                for _ in 1..iter {
                                    result *= x_d;
                                }
                                return result as f32;
                            }
                        }

                        if y_abs > 0x4f17_0000 {
                            // if y is NaN
                            if y_abs > 0x7f80_0000 {
                                if x_u == 0x3f80_0000 {
                                    // x = 1.0f
                                    // pow(1, NaN) = 1
                                    return 1.0;
                                }
                                // pow(x, NaN) = NaN
                                return y;
                            }
                            // x^y will be overflow / underflow in single precision.  Set y to a
                            // large enough exponent but not too large, so that the computations
                            // won't be overflow in double precision.
                            y = f32::from_bits((y_u & 0x8000_0000).wrapping_add(0x4f800000u32));
                        }
                    }
                }
            }
        }
    }

    const E_BIAS: u32 = (1u32 << (8 - 1u32)) - 1u32;
    let mut ex = -(E_BIAS as i32);
    let mut sign: u64 = 0;

    if ((x_u & 0x801f_ffffu32) == 0) || x_u >= 0x7f80_0000u32 || x_u < 0x0080_0000u32 {
        if x.is_nan() {
            return f32::NAN;
        }

        if x_u == 0x3f80_0000 {
            return 1.;
        }

        let x_is_neg = x.to_bits() > 0x8000_0000;

        if x == 0.0 {
            let out_is_neg = x_is_neg && is_odd_integer(f32::from_bits(y_u));
            if y_u > 0x8000_0000u32 {
                // pow(0, negative number) = inf
                return if x_is_neg {
                    f32::NEG_INFINITY
                } else {
                    f32::INFINITY
                };
            }
            // pow(0, positive number) = 0
            return if out_is_neg { -0.0 } else { 0.0 };
        }

        if x_abs == 0x7f80_0000u32 {
            // x = +-Inf
            let out_is_neg = x_is_neg && is_odd_integer(f32::from_bits(y_u));
            if y_u >= 0x7fff_ffff {
                return if out_is_neg { -0.0 } else { 0.0 };
            }
            return if out_is_neg { -0.0 } else { 0.0 };
        }

        if x_abs > 0x7f80_0000 {
            // x is NaN.
            // pow (aNaN, 0) is already taken care above.
            return x;
        }

        // Normalize denormal inputs.
        if x_abs < 0x0080_0000u32 {
            ex = ex.wrapping_sub(64);
            x *= f32::from_bits(0x5f800000);
        }

        // x is finite and negative, and y is a finite integer.
        if x.is_sign_negative() {
            if is_integer(y) {
                x = -x;
                if is_odd_integer(y) {
                    sign = 0x8000_0000_0000_0000u64;
                }
            } else {
                // pow( negative, non-integer ) = NaN
                return f32::NAN;
            }
        }
    }

    // x^y = 2^( y * log2(x) )
    //     = 2^( y * ( e_x + log2(m_x) ) )
    // First we compute log2(x) = e_x + log2(m_x)
    x_u = x.to_bits();

    // Extract exponent field of x.
    ex = ex.wrapping_add((x_u >> 23) as i32);
    let e_x = ex as f64;
    // Use the highest 7 fractional bits of m_x as the index for look up tables.
    let x_mant = x_u & ((1u32 << 23) - 1);
    let idx_x = (x_mant >> (23 - 7)) as i32;
    // Add the hidden bit to the mantissa.
    // 1 <= m_x < 2
    let m_x = f32::from_bits(x_mant | 0x3f800000);

    // Reduced argument for log2(m_x):
    //   dx = r * m_x - 1.
    // The computation is exact, and -2^-8 <= dx < 2^-7.
    // Then m_x = (1 + dx) / r, and
    //   log2(m_x) = log2( (1 + dx) / r )
    //             = log2(1 + dx) - log2(r).

    let dx;
    #[cfg(any(
        all(
            any(target_arch = "x86", target_arch = "x86_64"),
            target_feature = "fma"
        ),
        all(target_arch = "aarch64", target_feature = "neon")
    ))]
    {
        use crate::logf::LOG_REDUCTION_F32;
        dx = f_fmlaf(
            m_x,
            f32::from_bits(LOG_REDUCTION_F32.0[idx_x as usize]),
            -1.0,
        ) as f64; // Exact.
    }
    #[cfg(not(any(
        all(
            any(target_arch = "x86", target_arch = "x86_64"),
            target_feature = "fma"
        ),
        all(target_arch = "aarch64", target_feature = "neon")
    )))]
    {
        use crate::log2::LOG_RANGE_REDUCTION;
        dx = f_fmla(
            m_x as f64,
            f64::from_bits(LOG_RANGE_REDUCTION[idx_x as usize]),
            -1.0,
        ); // Exact
    }

    // Degree-5 polynomial approximation:
    //   dx * P(dx) ~ log2(1 + dx)
    // Generated by Sollya with:
    // > P = fpminimax(log2(1 + x)/x, 5, [|D...|], [-2^-8, 2^-7]);
    // > dirtyinfnorm(log2(1 + x)/x - P, [-2^-8, 2^-7]);
    //   0x1.653...p-52
    const COEFFS: [u64; 6] = [
        0x3ff71547652b82fe,
        0xbfe71547652b7a07,
        0x3fdec709dc458db1,
        0xbfd715479c2266c9,
        0x3fd2776ae1ddf8f0,
        0xbfce7b2178870157,
    ];

    let dx2 = dx * dx; // Exact
    let c0 = f_fmla(dx, f64::from_bits(COEFFS[1]), f64::from_bits(COEFFS[0]));
    let c1 = f_fmla(dx, f64::from_bits(COEFFS[3]), f64::from_bits(COEFFS[2]));
    let c2 = f_fmla(dx, f64::from_bits(COEFFS[5]), f64::from_bits(COEFFS[4]));

    let p = f_polyeval3(dx2, c0, c1, c2);

    // s = e_x - log2(r) + dx * P(dx)
    // Approximation errors:
    //   |log2(x) - s| < ulp(e_x) + (bounds on dx) * (error bounds of P(dx))
    //                 = ulp(e_x) + 2^-7 * 2^-51
    //                 < 2^8 * 2^-52 + 2^-7 * 2^-43
    //                 ~ 2^-44 + 2^-50
    let s = f_fmla(dx, p, f64::from_bits(LOG2_R[idx_x as usize]) + e_x);

    // To compute 2^(y * log2(x)), we break the exponent into 3 parts:
    //   y * log(2) = hi + mid + lo, where
    //   hi is an integer
    //   mid * 2^6 is an integer
    //   |lo| <= 2^-7
    // Then:
    //   x^y = 2^(y * log2(x)) = 2^hi * 2^mid * 2^lo,
    // In which 2^mid is obtained from a look-up table of size 2^6 = 64 elements,
    // and 2^lo ~ 1 + lo * P(lo).
    // Thus, we have:
    //   hi + mid = 2^-6 * round( 2^6 * y * log2(x) )
    // If we restrict the output such that |hi| < 150, (hi + mid) uses (8 + 6)
    // bits, hence, if we use double precision to perform
    //   round( 2^6 * y * log2(x))
    // the lo part is bounded by 2^-7 + 2^(-(52 - 14)) = 2^-7 + 2^-38

    // In the following computations:
    //   y6  = 2^6 * y
    //   hm  = 2^6 * (hi + mid) = round(2^6 * y * log2(x)) ~ round(y6 * s)
    //   lo6 = 2^6 * lo = 2^6 * (y - (hi + mid)) = y6 * log2(x) - hm.
    let y6 = (y * f32::from_bits(0x42800000)) as f64; // Exact.
    let hm = (s * y6).round();

    let log2_rr = LOG2_R2_DD[idx_x as usize];

    // lo6 = 2^6 * lo.
    let lo6_hi = f_fmla(y6, e_x + f64::from_bits(log2_rr.1), -hm); // Exact
    // Error bounds:
    //   | (y*log2(x) - hm * 2^-6 - lo) / y| < err(dx * p) + err(LOG2_R_DD.lo)
    //                                       < 2^-51 + 2^-75
    let lo6 = f_fmla(y6, f_fmla(dx, p, f64::from_bits(log2_rr.0)), lo6_hi);

    // |2^(hi + mid) - exp2_hi_mid| <= ulp(exp2_hi_mid) / 2
    // Clamp the exponent part into smaller range that fits double precision.
    // For those exponents that are out of range, the final conversion will round
    // them correctly to inf/max float or 0/min float accordingly.
    let mut hm_i = hm as i64;
    hm_i = if hm_i > (1i64 << 15) {
        1 << 15
    } else if hm_i < (-(1i64 << 15)) {
        -(1 << 15)
    } else {
        hm_i
    };

    let idx_y = hm_i & 0x3f;

    // 2^hi
    let exp_hi_i = (hm_i >> 6).wrapping_shl(52);
    // 2^mid
    let exp_mid_i = EXP2_MID1[idx_y as usize].1;
    // (-1)^sign * 2^hi * 2^mid
    // Error <= 2^hi * 2^-53
    let exp2_hi_mid_i = (exp_hi_i.wrapping_add(exp_mid_i as i64) as u64).wrapping_add(sign);
    let exp2_hi_mid = f64::from_bits(exp2_hi_mid_i);

    // Degree-5 polynomial approximation P(lo6) ~ 2^(lo6 / 2^6) = 2^(lo).
    // Generated by Sollya with:
    // > P = fpminimax(2^(x/64), 5, [|1, D...|], [-2^-1, 2^-1]);
    // > dirtyinfnorm(2^(x/64) - P, [-0.5, 0.5]);
    // 0x1.a2b77e618f5c4c176fd11b7659016cde5de83cb72p-60
    const EXP2_COEFFS: [u64; 6] = [
        0x3ff0000000000000,
        0x3f862e42fefa39ef,
        0x3f0ebfbdff82a23a,
        0x3e8c6b08d7076268,
        0x3e03b2ad33f8b48b,
        0x3d75d870c4d84445,
    ];

    let lo6_sqr = lo6 * lo6;
    let d0 = f_fmla(
        lo6,
        f64::from_bits(EXP2_COEFFS[1]),
        f64::from_bits(EXP2_COEFFS[0]),
    );
    let d1 = f_fmla(
        lo6,
        f64::from_bits(EXP2_COEFFS[3]),
        f64::from_bits(EXP2_COEFFS[2]),
    );
    let d2 = f_fmla(
        lo6,
        f64::from_bits(EXP2_COEFFS[5]),
        f64::from_bits(EXP2_COEFFS[4]),
    );
    let pp = f_polyeval3(lo6_sqr, d0, d1, d2);

    let r = pp * exp2_hi_mid;
    r as f32
}

/// Power function for given value using FMA
#[inline]
pub fn dirty_powf(d: f32, n: f32) -> f32 {
    use crate::exp2f::dirty_exp2f;
    use crate::log2f::dirty_log2f;
    let value = d.abs();
    let lg = dirty_log2f(value);
    let c = dirty_exp2f(n * lg);
    if d < 0.0 {
        let y = n as i32;
        if y % 2 == 0 { c } else { -c }
    } else {
        c
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn powf_test() {
        assert!(
            (powf(2f32, 3f32) - 8f32).abs() < 1e-6,
            "Invalid result {}",
            powf(2f32, 3f32)
        );
        assert!(
            (powf(0.5f32, 2f32) - 0.25f32).abs() < 1e-6,
            "Invalid result {}",
            powf(0.5f32, 2f32)
        );
    }

    #[test]
    fn f_powf_test() {
        assert!(
            (f_powf(2f32, 3f32) - 8f32).abs() < 1e-6,
            "Invalid result {}",
            f_powf(2f32, 3f32)
        );
        assert!(
            (f_powf(0.5f32, 2f32) - 0.25f32).abs() < 1e-6,
            "Invalid result {}",
            f_powf(0.5f32, 2f32)
        );
    }

    #[test]
    fn dirty_powf_test() {
        println!("{}", dirty_powf(3., 3.));
        println!("{}", dirty_powf(27., 1. / 3.));
        assert!(
            (dirty_powf(2f32, 3f32) - 8f32).abs() < 1e-6,
            "Invalid result {}",
            dirty_powf(2f32, 3f32)
        );
        assert!(
            (dirty_powf(0.5f32, 2f32) - 0.25f32).abs() < 1e-6,
            "Invalid result {}",
            dirty_powf(0.5f32, 2f32)
        );
    }
}
