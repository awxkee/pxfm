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
use crate::dekker::Dekker;
use crate::dyadic_float::{DyadicFloat128, DyadicSign};
use crate::log1p::log1p_f64_dyadic;
use crate::log1p_dd::log1p_f64_dd;
use crate::polyeval::{f_polyeval8, f_polyeval18};
use crate::pow::{is_integer, is_odd_integer};
use crate::pow_exec::pow_expm1_1;
use crate::pow_tables::{EXP_T1_2_DYADIC, EXP_T2_2_DYADIC};

/// Computes (1+x)^y - 1
///
/// max found ULP 0.50013
//TODO: still many bad behaviours with ULP 0.501+- for some subnormals, too slow
#[inline]
pub fn f_compound_m1(x: f64, y: f64) -> f64 {
    /*
       Rules from IEEE 754-2019 for compound (x, n) with n integer:
           (a) compound (x, 0) is 1 for x >= -1 or quiet NaN
           (b) compound (-1, n) is +Inf and signals the divideByZero exception for n < 0
           (c) compound (-1, n) is +0 for n > 0
           (d) compound (+/-0, n) is 1
           (e) compound (+Inf, n) is +Inf for n > 0
           (f) compound (+Inf, n) is +0 for n < 0
           (g) compound (x, n) is qNaN and signals the invalid exception for x < -1
           (h) compound (qNaN, n) is qNaN for n <> 0.
    */

    let x_sign = x.is_sign_negative();
    let y_sign = y.is_sign_negative();

    let x_abs = x.to_bits() & 0x7fff_ffff_ffff_ffff;
    let y_abs = y.to_bits() & 0x7fff_ffff_ffff_ffff;

    const MANTISSA_MASK: u64 = (1u64 << 52) - 1;

    let y_mant = y.to_bits() & MANTISSA_MASK;
    let x_u = x.to_bits();
    let x_a = x_abs;
    let y_a = y_abs;

    // If x or y is signaling NaN
    if x.is_nan() || y.is_nan() {
        return f64::NAN;
    }

    let mut s = 1.0;

    // The double precision number that is closest to 1 is (1 - 2^-53), which has
    //   log2(1 - 2^-53) ~ -1.715...p-53.
    // So if |y| > |1075 / log2(1 - 2^-53)|, and x is finite:
    //   |y * log2(x)| = 0 or > 1075.
    // Hence, x^y will either overflow or underflow if x is not zero.
    if y_mant == 0
        || y_a > 0x43d7_4910_d52d_3052
        || x_u == 1f64.to_bits()
        || x_u >= f64::INFINITY.to_bits()
        || x_u < f64::MIN.to_bits()
    {
        // Exceptional exponents.
        if y == 0.0 {
            return 0.0;
        }

        // (h) compound(qNaN, n) is qNaN for n ≠ 0
        if x.is_nan() {
            if y != 0. {
                return x;
            } // propagate qNaN
            return 0.0;
        }

        // (d) compound(±0, n) is 1
        if x == 0.0 {
            return 0.0;
        }

        // (e, f) compound(+Inf, n)
        if x.is_infinite() && x > 0.0 {
            return if y > 0. { x } else { -1.0 };
        }

        // (g) compound(x, n) is qNaN and signals invalid for x < -1
        if x < -1.0 {
            // Optional: raise invalid explicitly
            return f64::NAN;
        }

        // (b, c) compound(-1, n)
        if x == -1.0 {
            return if y < 0. { f64::INFINITY } else { -1.0 };
        }

        match y_a {
            // 0x3fe0_0000_0000_0000 => {
            //     if x == 0.0 {
            //         return 0.0;
            //     }
            //     let z = Dekker::from_full_exact_add(x, 1.0).sqrt();
            //     if y_sign {
            //         const M_ONES: DyadicFloat128 = DyadicFloat128 {
            //             sign: DyadicSign::Neg,
            //             exponent: -127,
            //             mantissa: 0x80000000_00000000_00000000_00000000_u128,
            //         };
            //         let z = DyadicFloat128::new_from_f64(z.to_f64());
            //         (z.reciprocal() + M_ONES).fast_as_f64()
            //     } else {
            //         const M_ONES: DyadicFloat128 = DyadicFloat128 {
            //             sign: DyadicSign::Neg,
            //             exponent: -127,
            //             mantissa: 0x80000000_00000000_00000000_00000000_u128,
            //         };
            //         let z = DyadicFloat128::new_from_f64(z.to_f64());
            //         (z + M_ONES).fast_as_f64()
            //     };
            // }
            0x3ff0_0000_0000_0000 => {
                return if y_sign {
                    let z = DyadicFloat128::new_from_f64(x);
                    const ONES: DyadicFloat128 = DyadicFloat128 {
                        sign: DyadicSign::Pos,
                        exponent: -127,
                        mantissa: 0x80000000_00000000_00000000_00000000_u128,
                    };
                    const M_ONES: DyadicFloat128 = DyadicFloat128 {
                        sign: DyadicSign::Neg,
                        exponent: -127,
                        mantissa: 0x80000000_00000000_00000000_00000000_u128,
                    };
                    let p = (z + ONES).reciprocal() + M_ONES;
                    p.fast_as_f64()
                } else {
                    x
                };
            }
            0x4000_0000_0000_0000 => {
                const ONES: DyadicFloat128 = DyadicFloat128 {
                    sign: DyadicSign::Pos,
                    exponent: -127,
                    mantissa: 0x80000000_00000000_00000000_00000000_u128,
                };
                let z0 = DyadicFloat128::new_from_f64(x) + ONES;
                let z = z0 * z0;
                const M_ONES: DyadicFloat128 = DyadicFloat128 {
                    sign: DyadicSign::Neg,
                    exponent: -127,
                    mantissa: 0x80000000_00000000_00000000_00000000_u128,
                };
                return if y_sign {
                    (z.reciprocal() + M_ONES).fast_as_f64()
                } else {
                    f64::copysign((z + M_ONES).fast_as_f64(), x)
                };
            }
            _ => {}
        }

        // |y| > |1075 / log2(1 - 2^-53)|.
        if y_a >= 0x7ff0_0000_0000_0000 {
            // y is inf or nan
            if y_mant != 0 {
                // y is NaN
                // pow(1, NaN) = 1
                // pow(x, NaN) = NaN
                return if x_u == 1f64.to_bits() { 1.0 } else { y };
            }

            // Now y is +-Inf
            if f64::from_bits(x_abs).is_nan() {
                // pow(NaN, +-Inf) = NaN
                return x;
            }

            if x_a == 0x3ff0_0000_0000_0000 {
                // pow(+-1, +-Inf) = 1.0
                return 0.0;
            }

            if x == 0.0 && y_sign {
                // pow(+-0, -Inf) = +inf and raise FE_DIVBYZERO
                return f64::INFINITY;
            }
            // pow (|x| < 1, -inf) = +inf
            // pow (|x| < 1, +inf) = 0.0
            // pow (|x| > 1, -inf) = 0.0
            // pow (|x| > 1, +inf) = +inf
            return if (x_a < 1f64.to_bits()) == y_sign {
                f64::INFINITY
            } else {
                -1.0
            };
        }

        // y is finite and non-zero.

        if x_u == 1f64.to_bits() {
            // pow(1, y) = 1
            return 0.0;
        }

        if x == 0.0 {
            let out_is_neg = x_sign && is_odd_integer(y);
            if y_sign {
                // pow(0, negative number) = inf
                return if out_is_neg {
                    f64::NEG_INFINITY
                } else {
                    f64::INFINITY
                };
            }
            // pow(0, positive number) = 0
            return -1.0;
        }

        if x_a == f64::INFINITY.to_bits() {
            let out_is_neg = x_sign && is_odd_integer(y);
            if y_sign {
                return if out_is_neg { -1.0 } else { 1.0 };
            }
            return if out_is_neg {
                f64::NEG_INFINITY
            } else {
                f64::INFINITY
            };
        }

        if x_a > f64::INFINITY.to_bits() {
            // x is NaN.
            // pow (aNaN, 0) is already taken care above.
            return x;
        }

        // x is finite and negative, and y is a finite integer.
        if x_sign {
            if is_integer(y) {
                if is_odd_integer(y) {
                    // sign = -1.0;
                    static CS: [f64; 2] = [1.0, -1.0];

                    // set sign to 1 for y even, to -1 for y odd
                    let y_parity = if (y.abs()) >= f64::from_bits(0x4340000000000000) {
                        0usize
                    } else {
                        (y as i64 & 0x1) as usize
                    };
                    s = CS[y_parity];
                }
            } else {
                // pow( negative, non-integer ) = NaN
                return f64::NAN;
            }
        }

        // y is finite and non-zero.

        if x_u == 1f64.to_bits() {
            // pow(1, y) = 1
            return 0.0;
        }

        if x == 0.0 {
            let out_is_neg = x_sign && is_odd_integer(y);
            if y_sign {
                // pow(0, negative number) = inf
                return if out_is_neg {
                    f64::NEG_INFINITY
                } else {
                    f64::INFINITY
                };
            }
            // pow(0, positive number) = 0
            return if out_is_neg { -0.0 } else { 0.0 };
        }

        if x_a == f64::INFINITY.to_bits() {
            let out_is_neg = x_sign && is_odd_integer(y);
            if y_sign {
                return -1.;
            }
            return if out_is_neg {
                f64::NEG_INFINITY
            } else {
                f64::INFINITY
            };
        }

        if x_a > f64::INFINITY.to_bits() {
            // x is NaN.
            // pow (aNaN, 0) is already taken care above.
            return x;
        }
    }

    let ax = x.to_bits() & 0x7fff_ffff_ffff_ffff;
    let ay = y.to_bits() & 0x7fff_ffff_ffff_ffff;

    // evaluate (1+x)^y explicitly for integer y in [-16,16] range and |x|<2^64
    if y.floor() == y
        && ay <= 0x4030_0000_0000_0000u64
        && ax <= 0x43e0_0000_0000_0000u64
        && ax > 0x3cc0_0000_0000_0000
    {
        let iter_count = y.abs() as usize;
        const ONES: DyadicFloat128 = DyadicFloat128 {
            sign: DyadicSign::Pos,
            exponent: -127,
            mantissa: 0x80000000_00000000_00000000_00000000_u128,
        };
        const M_ONES: DyadicFloat128 = DyadicFloat128 {
            sign: DyadicSign::Neg,
            exponent: -127,
            mantissa: 0x80000000_00000000_00000000_00000000_u128,
        };
        let s = DyadicFloat128::new_from_f64(x) + ONES;
        let mut p = s;
        for _ in 0..iter_count - 1 {
            p = p * s;
        }
        return if y.is_sign_negative() {
            (p.reciprocal() + M_ONES).fast_as_f64()
        } else {
            (p + M_ONES).fast_as_f64()
        };
    }

    // approximate log(x)
    let l = log1p_f64_dd(x);

    let ey = ((y.to_bits() >> 52) & 0x7ff) as i32;
    if ey < 0x36 || ey >= 0x7f5 {
        return 0.;
    }

    let r = Dekker::quick_mult_f64(l, y);
    if r.hi.abs() > 1e-250 && r.hi.abs() < 70. && ey.abs() < 1050 {
        let res = pow_expm1_1(r, s);

        let res_min = res.hi + dd_fmla(f64::from_bits(0x3c9f066666666666), -res.hi, res.lo);
        let res_max = res.hi + dd_fmla(f64::from_bits(0x3c9f066666666666), res.hi, res.lo);
        if res_min == res_max {
            return res_max;
        }
    }

    compound_accurate(x, y)
}

fn expm1_dyadic_poly(x: DyadicFloat128) -> DyadicFloat128 {
    // compound_m1_expm1.sollya
    // Sollya:
    // pretty = proc(u) {
    //   return ~(floor(u*1000)/1000);
    // };
    //
    // // display = hexadecimal;
    //
    // n = 9;
    // P = 128;
    // N = 1;
    //
    // n = 7;
    // d = [-0.00016923,0.00016923];
    // f = 1;
    // w = 1/expm1(x);
    // p = remez(f, n, d, w);
    // Q = horner(fpminimax(expm1(x), [|1,2,3,4,5,6,7,8|], [|0,P...|], d, relative, floating));
    // e = -log2(dirtyinfnorm(Q * w - f, d));
    // print ("exp(x) :\n  Q(x) =", Q);
    //
    // for i from 1 to degree(Q) do print(coeff(Q, i));
    //
    // print ("  precision:", pretty(e));

    // Print in Sage:
    // from sage.all import *
    // # Sin coeffs
    // def format_hex(value):
    //     l = hex(value)[2:]
    //     n = 8
    //     x = [l[i:i + n] for i in range(0, len(l), n)]
    //     return "0x" + "_".join(x) + "_u128"
    //
    // def print_dyadic(value):
    //     (s, m, e) = RealField(128)(value).sign_mantissa_exponent();
    //     print("DyadicFloat128 {")
    //     print(f"    sign: DyadicSign::{'Pos' if s >= 0 else 'Neg'},")
    //     print(f"    exponent: {e},")
    //     print(f"    mantissa: {format_hex(m)},")
    //     print("},")
    //
    // arr = [1,
    // 0.50000000000000000000000000000000000005583598166406,
    // 0.166666666666666666666666666678803415539971886586616,
    // 4.1666666666666666666666666652567805725779758247843e-2,
    // 8.3333333333333333307906111754404186695749965214285e-3,
    // 1.38888888888888888969592305807840282305676282703644e-3,
    // 1.98412698564902881841609894942365238668102733754495e-4,
    // 2.4801587295717097909900971900592787228179147973606e-5]
    //
    // for num in arr:
    //     print_dyadic(num)
    const Q_2: [DyadicFloat128; 8] = [
        DyadicFloat128 {
            sign: DyadicSign::Pos,
            exponent: -127,
            mantissa: 0x80000000_00000000_00000000_00000000_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Pos,
            exponent: -128,
            mantissa: 0x80000000_00000000_00000000_00000012_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Pos,
            exponent: -130,
            mantissa: 0xaaaaaaaa_aaaaaaaa_aaaaaaae_8351142d_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Pos,
            exponent: -132,
            mantissa: 0xaaaaaaaa_aaaaaaaa_aaaaaa99_591f3615_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Pos,
            exponent: -134,
            mantissa: 0x88888888_88888885_880aef72_fd782d3b_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Pos,
            exponent: -137,
            mantissa: 0xb60b60b6_0b60b612_d0d4a4da_9e3f0bbd_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Pos,
            exponent: -140,
            mantissa: 0xd00d00d2_ba789e46_9bf951aa_87eb4ed7_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Pos,
            exponent: -143,
            mantissa: 0xd00d00cf_43459a92_1f8b060c_0584c599_u128,
        },
    ];
    f_polyeval8(
        x, Q_2[0], Q_2[1], Q_2[2], Q_2[3], Q_2[4], Q_2[5], Q_2[6], Q_2[7],
    ) * x
}

// |x| < 0.125
#[cold]
fn expm1_dyadic_tiny(x: DyadicFloat128) -> DyadicFloat128 {
    const Q_2: [DyadicFloat128; 18] = [
        DyadicFloat128 {
            sign: DyadicSign::Pos,
            exponent: -127,
            mantissa: 0x80000000_00000000_00000000_00000000_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Pos,
            exponent: -129,
            mantissa: 0xffffffff_ffffffff_ffffffff_ffffffec_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Pos,
            exponent: -130,
            mantissa: 0xaaaaaaaa_aaaaaaaa_aaaaaaaa_aaaa23cf_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Pos,
            exponent: -132,
            mantissa: 0xaaaaaaaa_aaaaaaaa_aaaaaaaa_aab21d4c_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Pos,
            exponent: -134,
            mantissa: 0x88888888_88888888_88888888_cbfc5fc3_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Pos,
            exponent: -137,
            mantissa: 0xb60b60b6_0b60b60b_60b60b5d_8078803c_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Pos,
            exponent: -140,
            mantissa: 0xd00d00d0_0d00d00d_00cfde9d_ccec5afa_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Pos,
            exponent: -143,
            mantissa: 0xd00d00d0_0d00d00d_00d14863_bae41ecd_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Pos,
            exponent: -146,
            mantissa: 0xb8ef1d2a_b6399c7d_65858919_0c26f5ee_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Pos,
            exponent: -149,
            mantissa: 0x93f27dbb_c4fae397_3894d69f_d6034859_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Pos,
            exponent: -153,
            mantissa: 0xd7322b3f_aa2716c8_b5e3ab07_194b120e_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Pos,
            exponent: -156,
            mantissa: 0x8f76c77f_c6c4cbda_eca94b44_23795040_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Pos,
            exponent: -160,
            mantissa: 0xb092309d_44a35a14_a638af0e_1067db7d_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Pos,
            exponent: -164,
            mantissa: 0xc9cba546_0070b9c7_ca2388ad_c9f3fc80_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Pos,
            exponent: -168,
            mantissa: 0xd73f9eea_d77e61bc_a59dd899_27a84c8b_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Pos,
            exponent: -172,
            mantissa: 0xd73f9f9c_fd658e0d_28f98933_ae7d0637_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Pos,
            exponent: -176,
            mantissa: 0xcaa0d334_a00ffa28_beb527a8_c05f021b_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Pos,
            exponent: -180,
            mantissa: 0xb412ae1e_c8fd4452_f35db8f0_4304ad4b_u128,
        },
    ];
    f_polyeval18(
        x, Q_2[0], Q_2[1], Q_2[2], Q_2[3], Q_2[4], Q_2[5], Q_2[6], Q_2[7], Q_2[8], Q_2[9], Q_2[10],
        Q_2[11], Q_2[12], Q_2[13], Q_2[14], Q_2[15], Q_2[16], Q_2[17],
    ) * x
}

// /* put in r an approximation of exp(x), for |x| < 744.45,
// with relative error < 2^-121.70 */
fn compound_expm1_dyadic(x: DyadicFloat128) -> DyadicFloat128 {
    // x < 0.125
    if x.exponent <= -130 {
        return expm1_dyadic_tiny(x);
    }

    const LOG2_INV: DyadicFloat128 = DyadicFloat128 {
        sign: DyadicSign::Pos,
        exponent: -115,
        mantissa: 0xb8aa_3b29_5c17_f0bc_0000_0000_0000_0000_u128,
    };

    const LOG2: DyadicFloat128 = DyadicFloat128 {
        sign: DyadicSign::Pos,
        exponent: -128,
        mantissa: 0xb172_17f7_d1cf_79ab_c9e3_b398_03f2_f6af_u128,
    };

    let mut bk = x * LOG2_INV;

    let unbiased = bk.biased_exponent();
    if unbiased >= 21 {
        return if x.sign == DyadicSign::Pos {
            DyadicFloat128 {
                sign: DyadicSign::Pos,
                exponent: 1270,
                mantissa: u128::MAX,
            }
        } else {
            DyadicFloat128 {
                sign: DyadicSign::Neg,
                exponent: -127,
                mantissa: 0x8000_0000_0000_0000_0000_0000_0000_0000_u128,
            }
        };
    }

    let k = bk.trunc_to_i64(); /* k = trunc(K) [rounded towards zero, exact] */
    /* The rounding error of mul_dint_int64() is bounded by 6 ulps, thus since
    |K| <= 4399162*log(2) < 3049267, the error on K is bounded by 2^-103.41.
    This error is divided by 2^12 below, thus yields < 2^-115.41. */
    bk = LOG2.mul_int64(k);
    bk.exponent -= 12;
    bk.sign = bk.sign.negate();
    let y = x + bk;

    let bm = k >> 12;
    let i2 = (k >> 6) & 0x3f;
    let i1 = k & 0x3f;

    let mut r = expm1_dyadic_poly(y);
    let di20 = EXP_T1_2_DYADIC[i2 as usize];
    let di21 = EXP_T2_2_DYADIC[i1 as usize];

    let m_ones = DyadicFloat128 {
        sign: DyadicSign::Neg,
        exponent: -127,
        mantissa: 0x80000000_00000000_00000000_00000000_u128,
    };
    // exp(x) = 2^k*(exp(r) - 1) + (2^k - 1) = 2^k*exp(r) - 2^k + 2^k - 1
    let mut pz = di20 * di21;
    pz.exponent += bm as i16;
    pz = pz + m_ones;

    r = di20 * r;
    r = di21 * r;

    r.exponent += bm as i16; /* exact */
    r = r + pz;
    r
}

#[cold]
fn compound_accurate(x: f64, y: f64) -> f64 {
    /* the idea of returning res_max instead of res_min is due to Laurent
    Théry: it is better in case of underflow since res_max = +0 always. */

    let f_y = DyadicFloat128::new_from_f64(y);

    let log_dyad = log1p_f64_dyadic(x);

    let r = log_dyad * f_y;

    let result = compound_expm1_dyadic(r);
    // 2^R.ex <= R < 2^(R.ex+1)

    result.fast_as_f64()
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_compound_exotic() {
        assert_eq!(f_compound_m1(
11944758478933760000000000000000000000000000000000000000000000000000000000000000000000000000000000000000.,
            -1242262631503757300000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000.,
        ), -1.);
    }

    #[test]
    fn test_compound_m1() {
        assert_eq!(f_compound_m1(3., 2.8927001953125), 54.154259038961406);
        assert_eq!(
            f_compound_m1(-0.43750000000000044, 19.),
            -0.9999821216263793
        );
        assert_eq!(
            f_compound_m1(127712., -2.0000000000143525),
            -0.9999999999386903
        );
        assert_eq!(
            f_compound_m1(-0.11718749767214207, 2893226081485815000000000000000.),
            -1.
        );
        assert_eq!(
            f_compound_m1(2418441935074801400000000., 512.),
            f64::INFINITY
        );
        assert_eq!(
            f_compound_m1(32.50198364245834, 128000.00000000093),
            f64::INFINITY
        );
        assert_eq!(
            f_compound_m1(1.584716796877785, 0.0000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000004168916810703412),
            0.0000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000003958869879428553
        );
        assert_eq!(
            f_compound_m1(
                -0.000000000000000000000000000000001997076793037533,
                366577337071337140000000000000000f64
            ),
            -0.5190938261758579
        );
        assert_eq!(f_compound_m1(2.1075630259863374, 0.5), 00.7628281328553664);
        assert_eq!(f_compound_m1(2.1078916412661783, 0.5), 0.7629213372315222);
        assert_eq!(f_compound_m1(3.0000000000001115, -0.5), -0.500000000000007);
        assert_eq!(
            f_compound_m1(0.0004873839215895903, 3.),
            0.0014628645098045245
        );

        assert_eq!(f_compound_m1(-0.483765364602732, 3.), -0.862424399516842);
        assert_eq!(f_compound_m1(3.0000001192092896, -2.), -0.9375000037252902);
        assert_eq!(f_compound_m1(29.38323424607434, -1.), -0.9670871115332561);

        assert_eq!(f_compound_m1(-0.4375, 4.), -0.8998870849609375);
        assert_eq!(
            f_compound_m1(-0.0039033182037826464, 3.),
            -0.011664306402886494
        );
        assert_eq!(
            f_compound_m1(0.000000000000000000000000000000000000007715336350455947,
                          -262034087537726030000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000.),
            -1.,
        );
        assert_eq!(f_compound_m1(10.000000059604645, 10.), 25937426005.44638);
        assert_eq!(f_compound_m1(10., -308.25471555814863), -1.0);
        assert_eq!(
            f_compound_m1(5.4172231599824623E-312, 9.4591068440831498E+164),
            -1.0
        );
        assert_eq!(
            f_compound_m1(5.8776567263633397E-39, 3.4223548116804511E-310),
            0.0
        );
        assert_eq!(
            f_compound_m1(5.8639503496997932E-148, -7.1936801558778956E+305),
            0.0
        );
        assert_eq!(
            f_compound_m1(0.9908447265624999,
                          -19032028850336152000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000.),
            -1.
        );
        assert_eq!(
            f_compound_m1(0.00000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000006952247559980936,
                          5069789834563405000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000.),
            -1.
        );
        assert_eq!(
            f_compound_m1(1.000000000000341,
                          -69261261804788370000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000.),
            -1.
        );
        assert_eq!(
            f_compound_m1(
                0.0000000000000001053438024827798,
                0.0000000000000001053438024827798
            ),
            0.000000000000000000000000000000011097316721530923
        );
        assert_eq!(
            f_compound_m1(
                0.00000000000000010755285551056508,
                0.00000000000000010755285551056508
            ),
            0.00000000000000000000000000000001156761672847649
        );

        assert_eq!(f_compound_m1(2.4324324, 1.4324324), 4.850778380908823);
        assert_eq!(f_compound_m1(2., 5.), 242.);
        assert_eq!(f_compound_m1(0.4324324, 126.4324324), 5.40545942023447e19);
        assert!(f_compound_m1(-0.4324324, 126.4324324).is_nan());
        assert_eq!(f_compound_m1(0.0, 0.0), 0.0);
        assert_eq!(f_compound_m1(0.0, -1. / 2.), 0.0);
        assert_eq!(f_compound_m1(-1., -1. / 2.), f64::INFINITY);
        assert_eq!(f_compound_m1(f64::INFINITY, -1. / 2.), -1.0);
        assert_eq!(f_compound_m1(f64::INFINITY, 1. / 2.), f64::INFINITY);
        assert_eq!(f_compound_m1(46.3828125, 46.3828125), 5.248159634773675e77);
    }

    #[test]
    fn test_expm1_dyadic() {
        let z = DyadicFloat128::new_from_f64(2.5);
        assert_eq!(compound_expm1_dyadic(z).fast_as_f64(), 11.182493960703473);
    }

    #[test]
    fn test_expm1_tiny_dyadic() {
        let z = DyadicFloat128::new_from_f64(0.1);
        assert_eq!(expm1_dyadic_tiny(z).fast_as_f64(), 0.10517091807564763);
    }
}
