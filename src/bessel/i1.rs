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
use crate::bessel::i0::i0_exp;
use crate::common::f_fmla;
use crate::double_double::DoubleDouble;
use crate::dyadic_float::{DyadicFloat128, DyadicSign};
use crate::horner::f_horner_polyeval23;
use crate::polyeval::{f_polyeval13, f_polyeval30};

/// Modified Bessel of the first kind order 1
///
/// Max found ULP 0.5003
pub fn f_i1(x: f64) -> f64 {
    if !x.is_normal() {
        if x == 0. {
            return 0.;
        }
        if x.is_infinite() {
            return if x.is_sign_positive() {
                f64::INFINITY
            } else {
                f64::NEG_INFINITY
            };
        }
        if x.is_nan() {
            return f64::NAN;
        }
    }

    let xb = x.to_bits() & 0x7fff_ffff_ffff_ffff;

    if xb >= 0x40864fe69ff9fec8u64 {
        return if x.is_sign_negative() {
            f64::NEG_INFINITY
        } else {
            f64::INFINITY
        };
    }

    static SIGN: [f64; 2] = [1., -1.];

    let sign_scale = SIGN[x.is_sign_negative() as usize];

    if xb < 0x401f000000000000u64 {
        // 7.75
        return i1_0_to_7p75(f64::from_bits(xb), sign_scale);
    }

    i1_asympt(f64::from_bits(xb), sign_scale)
}

/**
Computes
I1(x) = x/2 * (1 + 1 * (x/2)^2 + (x/2)^4 * P((x/2)^2))

Original series for Sollya generated in SageMath:

```python
from sage.all import *
from mpmath import mp, besseli, taylor

mp.prec = 500

def refined_approx(x):
    if x == 0:
        return 0
    p1 = besseli(1, x) / x * 2 - mp.mpf(1)
    x_over_two = mp.mpf(x) / 2
    x_over_two_sqr = x_over_two * x_over_two
    p2 = x_over_two_sqr * mp.mpf(0.5)
    p1 = p1 - p2
    x_over_two_p4 = x_over_two_sqr * x_over_two_sqr
    p1 = p1
    return p1

terms = 70

coeffs = taylor(lambda x: refined_approx(x), mp.mpf('0'), terms)

# Step 3: Build series in terms of y = (x/2)^2
R = PolynomialRing(RealField(450), 'y')
y = R.gen()
f = R(0)

for n in range(2, terms, 2):
    k = n // 2
    c = RealField(450)(coeffs[n])
    if n >= 1:
        f += R(c) * y**(k-1) * (4**k)
    else:
        f += R(c) * y**(k-1) * (4**k)

print(f)
```

See ./notes/bessel_sollya/bessel_i1_small.sollya for generation.

Poly relative err 2^(-111.212)
**/
#[inline]
fn i1_0_to_7p75(x: f64, sign_scale: f64) -> f64 {
    let dx = x;
    const ONE_OVER_4: f64 = 1. / 4.;
    let eval_x = DoubleDouble::quick_mult_f64(DoubleDouble::from_exact_mult(dx, dx), ONE_OVER_4);
    const C_HI: [u64; 13] = [
        0x3c18bce58901a2c3,
        0x3ba165e7c2d16482,
        0x3b253585cdc928c8,
        0x3aa69f7da8b356fb,
        0x3a254ad093b98735,
        0x39a1d02999d515bc,
        0x391aaac712cf2fce,
        0x3891f8e9b0232fad,
        0x3805c26d710171e9,
        0x377a7113ee429142,
        0x36a72afa13e6ad75,
        0x3685e04ea0ca6e89,
        0xb60be439ef13c652,
    ];

    let r = f_polyeval13(
        eval_x.to_f64(),
        f64::from_bits(C_HI[0]),
        f64::from_bits(C_HI[1]),
        f64::from_bits(C_HI[2]),
        f64::from_bits(C_HI[3]),
        f64::from_bits(C_HI[4]),
        f64::from_bits(C_HI[5]),
        f64::from_bits(C_HI[6]),
        f64::from_bits(C_HI[7]),
        f64::from_bits(C_HI[8]),
        f64::from_bits(C_HI[9]),
        f64::from_bits(C_HI[10]),
        f64::from_bits(C_HI[11]),
        f64::from_bits(C_HI[12]),
    );

    const C: [(u64, u64); 10] = [
        (0x3c55555555555555, 0x3fb5555555555555),
        (0x3c1c71c71c71c79b, 0x3f7c71c71c71c71c),
        (0xbbcf49f49f4a645e, 0x3f36c16c16c16c17),
        (0x3b85b66c77efc6ba, 0x3ee845c8a0ce5129),
        (0x3b3cbbc055c4d633, 0x3e927e4fb7789f5c),
        (0xbad604d9528a5187, 0x3e3522a43f65486a),
        (0xba7392e591b5e218, 0x3dd2c9758daf5cd0),
        (0x39f830be3775eba4, 0x3d6ab81ea75fcdf4),
        (0xb990dcd77d3ef549, 0x3cff17697cf1cf13),
        (0x392df0f806a62190, 0x3c8e2637bef9ff1e),
    ];

    let mut p = DoubleDouble::mul_f64_add(eval_x, r, DoubleDouble::from_bit_pair(C[9]));
    p = DoubleDouble::mul_add(eval_x, p, DoubleDouble::from_bit_pair(C[8]));
    p = DoubleDouble::mul_add(eval_x, p, DoubleDouble::from_bit_pair(C[7]));
    p = DoubleDouble::mul_add(eval_x, p, DoubleDouble::from_bit_pair(C[6]));
    p = DoubleDouble::mul_add(eval_x, p, DoubleDouble::from_bit_pair(C[5]));
    p = DoubleDouble::mul_add(eval_x, p, DoubleDouble::from_bit_pair(C[4]));
    p = DoubleDouble::mul_add(eval_x, p, DoubleDouble::from_bit_pair(C[3]));
    p = DoubleDouble::mul_add(eval_x, p, DoubleDouble::from_bit_pair(C[2]));
    p = DoubleDouble::mul_add(eval_x, p, DoubleDouble::from_bit_pair(C[1]));
    p = DoubleDouble::mul_add(eval_x, p, DoubleDouble::from_bit_pair(C[0]));

    let eval_sqr = DoubleDouble::quick_mult(eval_x, eval_x);

    let mut z = DoubleDouble::mul_f64_add_f64(eval_x, 0.5, 1.);
    z = DoubleDouble::mul_add(p, eval_sqr, z);

    let x_over_05 = DoubleDouble::from_exact_mult(x, 0.5);

    let r = DoubleDouble::quick_mult(z, x_over_05);

    let err = f_fmla(
        r.hi,
        f64::from_bits(0x3c20000000000000),
        f64::from_bits(0x3980000000000000),
    );
    let ub = r.hi + (r.lo + err);
    let lb = r.hi + (r.lo - err);
    if ub == lb {
        return r.to_f64() * sign_scale;
    }
    i1_0_to_7p5_hard(x, sign_scale)
}

#[cold]
#[inline(never)]
fn i1_0_to_7p5_hard(x: f64, sign_scale: f64) -> f64 {
    const P: [DyadicFloat128; 23] = [
        DyadicFloat128 {
            sign: DyadicSign::Pos,
            exponent: -131,
            mantissa: 0xaaaaaaaa_aaaaaaaa_aaaaaaaa_aaaaaaaa_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Pos,
            exponent: -135,
            mantissa: 0xe38e38e3_8e38e38e_38e38e38_e38e3ad0_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Pos,
            exponent: -139,
            mantissa: 0xb60b60b6_0b60b60b_60b60b60_b609b8cc_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Pos,
            exponent: -144,
            mantissa: 0xc22e4506_72894ab6_cd8efb11_d497871a_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Pos,
            exponent: -149,
            mantissa: 0x93f27dbb_c4fae397_780b69f4_926f09e4_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Pos,
            exponent: -155,
            mantissa: 0xa91521fb_2a434d3f_649f54e6_8d87bb1c_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Pos,
            exponent: -161,
            mantissa: 0x964bac6d_7ae67d8d_aec65bfd_1f6265fa_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Pos,
            exponent: -168,
            mantissa: 0xd5c0f53a_fe6fa17f_8c9309f3_dffc1e1a_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Pos,
            exponent: -175,
            mantissa: 0xf8bb4be7_8e7896b0_4e3d30b7_818b3816_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Pos,
            exponent: -182,
            mantissa: 0xf131bdf7_cff8d032_93e2aa28_3791a7ca_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Pos,
            exponent: -189,
            mantissa: 0xc5e72c48_0d1aeb77_d6bda068_c9b4e2b2_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Pos,
            exponent: -196,
            mantissa: 0x8b2f3e16_8a9fcdd2_bfee896d_351e9414_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Pos,
            exponent: -204,
            mantissa: 0xa9ac2e6e_5fc674fa_4d35755f_27d0abd8_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Pos,
            exponent: -212,
            mantissa: 0xb4fbed42_8d38f1f6_80b0cf5d_aae25ff6_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Pos,
            exponent: -220,
            mantissa: 0xaa5684f2_47c8645f_a8931821_f236d34a_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Pos,
            exponent: -228,
            mantissa: 0x8e814593_5ee16079_f0a5050a_66c3d544_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Pos,
            exponent: -237,
            mantissa: 0xd5573007_36641d12_8895577d_6cbd4444_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Pos,
            exponent: -245,
            mantissa: 0x8fbabd80_ed9f8d6d_1e301d22_b46e56a2_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Pos,
            exponent: -254,
            mantissa: 0xaefe68ca_093c3c26_2ce2eea8_9b6ec820_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Pos,
            exponent: -263,
            mantissa: 0xc8a8d97d_c60320a1_58888c31_cbd4ac2a_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Pos,
            exponent: -274,
            mantissa: 0xecc60a7a_17661b22_e1390e47_b31f6a46_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Pos,
            exponent: -278,
            mantissa: 0xf7964e72_ed6ab977_8bc70ede_d6163034_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Neg,
            exponent: -284,
            mantissa: 0xa9328051_1d8e8a54_31a21ed7_9e9491e8_u128,
        },
    ];

    let dx = DyadicFloat128::new_from_f64(x);
    let mut x_over_two = dx;
    x_over_two.exponent -= 1; // * 0.5
    let mut eval_x = dx * dx;
    eval_x.exponent -= 2; // * 0.25

    let p = f_horner_polyeval23(
        eval_x, P[0], P[1], P[2], P[3], P[4], P[5], P[6], P[7], P[8], P[9], P[10], P[11], P[12],
        P[13], P[14], P[15], P[16], P[17], P[18], P[19], P[20], P[21], P[22],
    );

    let mut eval_x_over_two = eval_x;
    eval_x_over_two.exponent -= 1; // * 0.5

    const ONES: DyadicFloat128 = DyadicFloat128 {
        sign: DyadicSign::Pos,
        exponent: -127,
        mantissa: 0x80000000_00000000_00000000_00000000_u128,
    };

    let eval_x_sqr = eval_x * eval_x;
    let r = (p * eval_x_sqr + eval_x_over_two + ONES) * x_over_two;
    r.fast_as_f64() * sign_scale
}

/**
Asymptotic expansion for I1.

Computes:
sqrt(x) * exp(-x) * I1(x) = Pn(1/x)/Qn(1/x)
hence:
I1(x) = Pn(1/x)/Qm(1/x)*exp(x)/sqrt(x)

Generated by Sollya:
```python
bessel_i1_approximant_big = library("./cmake-build-release/libbessel_sollya.dylib");

prec = 1000;

f = bessel_i1_approximant_big(1/x);
d = [1/713.98, 1/7.75];
w = 1;
pf = remez(f, 29, d, 1, 1e-36);

for i from 0 to degree(pf) do {
    print("'", coeff(pf, i), "',");
};
```
See ./notes/bessel_sollya/bessel_i1_big.sollya
**/
#[inline]
fn i1_asympt(x: f64, sign_scale: f64) -> f64 {
    let dx = x;
    let recip = DoubleDouble::from_recip(x);
    let z = f_polyeval30(
        recip,
        DoubleDouble::from_bit_pair((0xbc79892019caa6d3, 0x3fd9884533d43651)),
        DoubleDouble::from_bit_pair((0xbc697dbc64ad0e62, 0xbfc32633e6df29d3)),
        DoubleDouble::from_bit_pair((0xbc4868553039f2c6, 0xbfa7efc0e083df1c)),
        DoubleDouble::from_bit_pair((0x3c262ab46d250afc, 0xbfa4f1c8f1cb96bf)),
        DoubleDouble::from_bit_pair((0x3c3725e6f350b2d2, 0xbfad73bfa1ac95cd)),
        DoubleDouble::from_bit_pair((0x3c585d782fcc8a2e, 0xbfbc7a536a156d73)),
        DoubleDouble::from_bit_pair((0x3c281cb1ff2181c3, 0xbfc6e2ec43307f2c)),
        DoubleDouble::from_bit_pair((0xbcbb96ce81325b99, 0xc029f791fab12879)),
        DoubleDouble::from_bit_pair((0x3d35e6400c1a6d27, 0x4093677a9b1df97a)),
        DoubleDouble::from_bit_pair((0xbd81c2fa9c046cae, 0xc0f84034cc73eaa6)),
        DoubleDouble::from_bit_pair((0xbdfe2acaadb6e08f, 0x41580fb1ac1d97e9)),
        DoubleDouble::from_bit_pair((0x3e1f046a240bedfc, 0xc1b34679de6bff49)),
        DoubleDouble::from_bit_pair((0xbe64dccf496a0777, 0x42093665cfd5830f)),
        DoubleDouble::from_bit_pair((0x3eda0e215fd84ae7, 0xc25b28c4aa2c0e46)),
        DoubleDouble::from_bit_pair((0x3f4d45acd73886d4, 0x42a840d4046abf87)),
        DoubleDouble::from_bit_pair((0x3f9c97a93fea74c3, 0xc2f20a448bf6ec83)),
        DoubleDouble::from_bit_pair((0xbfd2a1309058bd4b, 0x43366c7ccc0c8a4e)),
        DoubleDouble::from_bit_pair((0xc017d40a12d5dbe6, 0xc3775353c42fe969)),
        DoubleDouble::from_bit_pair((0xc034bebc70f7dbb2, 0x43b44e5d4be7aa3a)),
        DoubleDouble::from_bit_pair((0xc089dd88abe28863, 0xc3ed8b0ca110b9f2)),
        DoubleDouble::from_bit_pair((0x40a3dd091a1e566e, 0x4421e612f7615754)),
        DoubleDouble::from_bit_pair((0xc0e2221903b12d9d, 0xc451f7074998df65)),
        DoubleDouble::from_bit_pair((0x410e90d2935f061a, 0x447da19db791c3ce)),
        DoubleDouble::from_bit_pair((0xc13c8a2749071af7, 0xc4a3d86387fa59e2)),
        DoubleDouble::from_bit_pair((0x4147b60d27af7e3d, 0x44c53a42bd3e2473)),
        DoubleDouble::from_bit_pair((0x4170c2d328f37b19, 0xc4e1b1024f0754ca)),
        DoubleDouble::from_bit_pair((0xc1891cb73bd261f5, 0x44f62166eab66554)),
        DoubleDouble::from_bit_pair((0x4170804b327e8716, 0xc503899c8a45877a)),
        DoubleDouble::from_bit_pair((0x419159386c2195c4, 0x4505b4153700642d)),
        DoubleDouble::from_bit_pair((0x419e7f4169901c41, 0xc4f6d0344a5382ba)),
    );

    let e = i0_exp(dx * 0.5);
    let r_sqrt = DoubleDouble::from_rsqrt(dx);
    DoubleDouble::quick_mult(z * r_sqrt * e, e).to_f64() * sign_scale
}

#[cfg(test)]
mod tests {
    use super::*;
    #[test]
    fn test_fi1() {
        assert_eq!(f_i1(7.482812501363189), 245.58002285881892);
        assert!(f_i1(f64::NAN).is_nan());
        assert_eq!(f_i1(f64::INFINITY), f64::INFINITY);
        assert_eq!(f_i1(f64::NEG_INFINITY), f64::NEG_INFINITY);
        assert_eq!(f_i1(0.01), 0.005000062500260418);
        assert_eq!(f_i1(-0.01), -0.005000062500260418);
        assert_eq!(f_i1(-9.01), -1040.752038018038);
        assert_eq!(f_i1(9.01), 1040.752038018038);
    }
}
