/*
 * // Copyright (c) Radzivon Bartoshyk 8/2025. All rights reserved.
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
use crate::bessel::i0::{bessel_rsqrt_hard, i0_exp};
use crate::common::f_fmla;
use crate::double_double::DoubleDouble;
use crate::dyadic_float::{DyadicFloat128, DyadicSign};
use crate::exponents::rational128_exp;
use crate::logs::{log_dd, log_dyadic};
use crate::polyeval::{f_horner_polyeval13, f_polyeval11, f_polyeval22};

/// Modified Bessel of the second kind order 1
///
/// Max ULP 0.5
pub fn f_k1(x: f64) -> f64 {
    if x < 0. {
        return f64::NAN;
    }

    if !x.is_normal() {
        if x == 0. {
            return f64::INFINITY;
        }
        if x.is_infinite() {
            return if x.is_sign_positive() { 0. } else { f64::NAN };
        }
        if x.is_nan() {
            return x + x;
        }
    }

    let xb = x.to_bits();

    if xb >= 0x4086140538aa7d38u64 {
        // 706.5025494880165
        return 0.;
    }

    if xb <= 0x3ff0000000000000 {
        return k1_small(x);
    }

    k1_asympt(x)
}

#[inline(always)]
fn i1_0_to_1_hard(x: f64) -> DyadicFloat128 {
    static P: [DyadicFloat128; 15] = [
        DyadicFloat128 {
            sign: DyadicSign::Pos,
            exponent: -131,
            mantissa: 0xaaaaaaaa_aaaaaaaa_aaaaaaaa_aaaaaaac_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Pos,
            exponent: -135,
            mantissa: 0xe38e38e3_8e38e38e_38e38e38_e38e0298_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Pos,
            exponent: -139,
            mantissa: 0xb60b60b6_0b60b60b_60b60b60_b7197728_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Pos,
            exponent: -144,
            mantissa: 0xc22e4506_72894ab6_cd8efb0d_d3d0066a_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Pos,
            exponent: -149,
            mantissa: 0x93f27dbb_c4fae397_780b71cc_6113c272_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Pos,
            exponent: -155,
            mantissa: 0xa91521fb_2a434d3f_648ce9bb_6361231e_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Pos,
            exponent: -161,
            mantissa: 0x964bac6d_7ae67d8d_caeaeaca_95245b64_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Pos,
            exponent: -168,
            mantissa: 0xd5c0f53a_fe6fa144_b6920c03_c832a616_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Pos,
            exponent: -175,
            mantissa: 0xf8bb4be7_8e78ed70_23d3c975_59c1dd1a_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Pos,
            exponent: -182,
            mantissa: 0xf131bdf7_cf9d2295_7b5b2f39_5fd63658_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Pos,
            exponent: -189,
            mantissa: 0xc5e72c48_52d48356_90683091_c1470906_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Pos,
            exponent: -196,
            mantissa: 0x8b2f3df0_b545dd95_8142e6bb_6b214694_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Pos,
            exponent: -204,
            mantissa: 0xa9ac4aff_88b809e3_1f38bca3_f3791d62_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Pos,
            exponent: -212,
            mantissa: 0xb4edad11_a4136582_20c52a39_15dab508_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Pos,
            exponent: -220,
            mantissa: 0xae8dd7b3_aae68b11_cb5f8a62_0a7ed726_u128,
        },
    ];

    let dx = DyadicFloat128::new_from_f64(x);
    let mut x_over_two = dx;
    x_over_two.exponent -= 1; // * 0.5
    let mut eval_x = dx * dx;
    eval_x.exponent -= 2; // * 0.25

    let mut p = P[14];
    for i in (0..14).rev() {
        p = eval_x * p + P[i];
    }

    let mut eval_x_over_two = eval_x;
    eval_x_over_two.exponent -= 1; // * 0.5

    const ONES: DyadicFloat128 = DyadicFloat128 {
        sign: DyadicSign::Pos,
        exponent: -127,
        mantissa: 0x80000000_00000000_00000000_00000000_u128,
    };

    let eval_x_sqr = eval_x * eval_x;
    (p * eval_x_sqr + eval_x_over_two + ONES) * x_over_two
}

#[inline]
fn i1_fast(x: f64) -> DoubleDouble {
    let dx = x;
    const ONE_OVER_4: f64 = 1. / 4.;
    let eval_x = DoubleDouble::quick_mult_f64(DoubleDouble::from_exact_mult(dx, dx), ONE_OVER_4);
    let r = f_polyeval11(
        eval_x.to_f64(),
        f64::from_bits(0x3e3522a43f65486a),
        f64::from_bits(0x3dd2c9758daf5c4f),
        f64::from_bits(0x3d6ab81ea760d0e0),
        f64::from_bits(0x3cff17697b809481),
        f64::from_bits(0x3c8e263939eab38f),
        f64::from_bits(0x3c18bbccb2c8332d),
        f64::from_bits(0x3ba1fad69e0b4040),
        f64::from_bits(0xbb463fe81b373085),
        f64::from_bits(0x3b3b0efc43abee12),
        f64::from_bits(0xbb1f88a199a5785b),
        f64::from_bits(0x3af095e96e675e05),
    );

    const C: [(u64, u64); 5] = [
        (0x3c55555555555552, 0x3fb5555555555555),
        (0x3c1c71c71c722160, 0x3f7c71c71c71c71c),
        (0xbbcf49f4a260e640, 0x3f36c16c16c16c17),
        (0x3b85b671d90e8cef, 0x3ee845c8a0ce5129),
        (0x3b3cb1dd9afdbed0, 0x3e927e4fb7789f5c),
    ];

    let mut p = DoubleDouble::mul_f64_add(eval_x, r, DoubleDouble::from_bit_pair(C[4]));
    p = DoubleDouble::mul_add(eval_x, p, DoubleDouble::from_bit_pair(C[3]));
    p = DoubleDouble::mul_add(eval_x, p, DoubleDouble::from_bit_pair(C[2]));
    p = DoubleDouble::mul_add(eval_x, p, DoubleDouble::from_bit_pair(C[1]));
    p = DoubleDouble::mul_add(eval_x, p, DoubleDouble::from_bit_pair(C[0]));

    let eval_sqr = DoubleDouble::quick_mult(eval_x, eval_x);

    let mut z = DoubleDouble::mul_f64_add_f64(eval_x, 0.5, 1.);
    z = DoubleDouble::mul_add(p, eval_sqr, z);

    let x_over_05 = DoubleDouble::from_exact_mult(x, 0.5);

    DoubleDouble::quick_mult(z, x_over_05)
}

/**
Series for
f(x) := BesselK(1, x) - Log(x)*BesselI(1, x) - 1/x

(z^17 (-6989 + 2520 EulerGamma - 2520 Log(2)))/4832746593583104000 + (
 z^15 (-1487 + 560 EulerGamma - 560 Log(2)))/3728971137024000 + (
 z^13 (-353 + 140 EulerGamma - 140 Log(2)))/4161798144000 + (
 z^9 (-131 + 60 EulerGamma - 60 Log(2)))/88473600 + (
 z^11 (-71 + 30 EulerGamma - 30 Log(2)))/5308416000 + (
 z^7 (-47 + 24 EulerGamma - 24 Log(2)))/442368 +
 1/64 z^3 (-5 + 4 EulerGamma - 4 Log(2)) + (
 z^5 (-5 + 3 EulerGamma - 3 Log(2)))/1152 +
 1/4 z (-1 + 2 EulerGamma - 2 Log(2))
**/
#[inline]
fn k1_small(x: f64) -> f64 {
    let rcp = DoubleDouble::from_recip(x);
    let x2 = DoubleDouble::from_exact_mult(x, x);

    const P: [(u64, u64); 13] = [
        (0xbc7037c12ba0b815, 0xbfd3b5b6028a83d6),
        (0xbc4037c12ba0b815, 0xbfb5dadb014541eb),
        (0x3c1e264dd50350dd, 0xbf7303ae729ff30f),
        (0x3bbeb7d012892972, 0xbf1d802af7a5dbc8),
        (0x3b188cf6afd16ea8, 0xbeba291822473f2f),
        (0xbaefd1fdc38805b2, 0xbe4e212a001aa46f),
        (0x3a74c549ad19d196, 0xbdd8630abd83ba61),
        (0x39fc11cc02c28333, 0xbd5d49398f1e78b6),
        (0x396ee0f0496b5608, 0xbcdb24176f948c55),
        (0xb8652b82b72eb179, 0xbbc80559d1876ef2),
        (0x37de784e3dddab40, 0xbb37f31cacac2b15),
        (0xb73338c4472cd7bf, 0xbaa4258454bb73f4),
        (0xb6a61d9412abaadf, 0xba0cfbb72084e258),
    ];

    let p = f_horner_polyeval13(
        x2,
        DoubleDouble::from_bit_pair(P[0]),
        DoubleDouble::from_bit_pair(P[1]),
        DoubleDouble::from_bit_pair(P[2]),
        DoubleDouble::from_bit_pair(P[3]),
        DoubleDouble::from_bit_pair(P[4]),
        DoubleDouble::from_bit_pair(P[5]),
        DoubleDouble::from_bit_pair(P[6]),
        DoubleDouble::from_bit_pair(P[7]),
        DoubleDouble::from_bit_pair(P[8]),
        DoubleDouble::from_bit_pair(P[9]),
        DoubleDouble::from_bit_pair(P[10]),
        DoubleDouble::from_bit_pair(P[11]),
        DoubleDouble::from_bit_pair(P[12]),
    );

    let lg = log_dd(x);
    let v_i = i1_fast(x);
    let z = DoubleDouble::mul_add(v_i, lg, rcp);
    let r = DoubleDouble::mul_f64_add(p, x, z);
    let err = f_fmla(
        r.hi,
        f64::from_bits(0x3c20000000000000), // 2^-61
        f64::from_bits(0x3a80000000000000), // 2^-87
    );
    let ub = r.hi + (r.lo + err);
    let lb = r.hi + (r.lo - err);
    if ub == lb {
        return r.to_f64();
    }
    k1_small_rational(x)
}

/**
Generated by SageMath:
```text
euler_gamma = R(euler)

lg = R(2).log()
expr = (z**17 * (R(-6989) + R(2520) * euler_gamma - R(2520) * lg))/R('4832746593583104000') \
        + (z**15 * (R(-1487)+R(560) * euler_gamma-R(560) *lg))/R('3728971137024000')\
        +(z**13 * (R(-353)+R(140) * euler_gamma-R(140) * lg))/R('4161798144000')\
        +(z**9 * (R(-131)+R(60) * euler_gamma-R(60) * lg))/R('88473600')\
        +(z**11 * (R(-71)+R(30) * euler_gamma - R(30) *lg))/R('5308416000')+(z**7 * (R(-47)+R(24) * euler_gamma - R(24) * lg))/R(442368)\
        +R(1/R(64)) * z**3 * (R(-5)+R(4) * euler_gamma-R(4) * lg)+ R(1/4) * z * (R(-1)+R(2) * euler_gamma-lg4)+(z**5 * (R(-5)+R(3) *euler_gamma-lg8))/R(1152)\
        + (z**21 * (R(-82451)+R(27720) * euler_gamma-R(27720) * lg))/R('8420577664659200409600000')\
        + (z**23 * (R(-42433)+R('13860') * euler_gamma-R('13860') * lg))/R('2223032503470028908134400000')\
        + (z**25 * (R('-1132133')+R('360360') * euler_gamma-R('360360')* lg))/R('36066479336297749005572505600000')\
        + (z**27 * (R('-1158863')+R('360360') * euler_gamma - R('360360') * lg))/R('26256396956824761276056784076800000')\
        + (z**29 * (R('-236749')+R('72072') * euler_gamma-R('72072') * lg))/R('4411074688746559894377539724902400000')\
        + (z**31 * (R('-4828073') + R('1441440') * euler_gamma - R('1441440') * lg))/R('84692634023933949972048762718126080000000')
```
**/
#[cold]
#[inline(never)]
fn k1_small_rational(x: f64) -> f64 {
    static P: [DyadicFloat128; 14] = [
        DyadicFloat128 {
            sign: DyadicSign::Neg,
            exponent: -129,
            mantissa: 0x9dadb014_541eb206_f8257417_02a02b59_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Neg,
            exponent: -131,
            mantissa: 0xaed6d80a_2a0f5903_7c12ba0b_815015ac_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Neg,
            exponent: -135,
            mantissa: 0x981d7394_ff98743b_36455f95_e46e4756_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Neg,
            exponent: -141,
            mantissa: 0xec0157bd_2ede3c29_05fdaeda_d1c27214_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Neg,
            exponent: -147,
            mantissa: 0xd148c112_39f977ce_e612a05d_22b00958_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Neg,
            exponent: -154,
            mantissa: 0xf1095000_d5237bfa_3fb87100_b63b2d98_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Neg,
            exponent: -161,
            mantissa: 0xc31855ec_1dd30567_56ca5cc5_cd489f65_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Neg,
            exponent: -169,
            mantissa: 0xea49cc78_f3c5ac7d_c67fa7af_99a21282_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Neg,
            exponent: -177,
            mantissa: 0xd920bb7c_a462a611_f0fb694a_9f804644_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Neg,
            exponent: -194,
            mantissa: 0xc02ace8c_3b7792a5_7056e5d6_2f1ec1ab_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Neg,
            exponent: -203,
            mantissa: 0xbf98e565_6158a430_f638444a_97fb71f9_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Neg,
            exponent: -212,
            mantissa: 0xa12c22a5_db9fa133_8c4472cd_7bedc36f_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Neg,
            exponent: -222,
            mantissa: 0xe7ddb904_2712c2c3_b2825575_5bd3e8be_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Neg,
            exponent: -231,
            mantissa: 0x9041bc8b_722802a8_cb714cd9_e56bd6bd_u128,
        },
    ];
    let dx = DyadicFloat128::new_from_f64(x);
    let x2 = dx * dx;

    let mut p = P[13];
    for i in (0..13).rev() {
        p = x2 * p + P[i];
    }

    let recip = DyadicFloat128::accurate_reciprocal(x);

    let lg = log_dyadic(x);
    let v_i = i1_0_to_1_hard(x);
    let z = v_i * lg + p * dx;
    let z0 = recip + z;
    z0.fast_as_f64()
}

/**
Generated by Wolfram Mathematica:
```text
<< FunctionApproximations`
ClearAll["Global`*"]
f[x_] := Sqrt[x] Exp[x] BesselK[1, x]
g[z_] := f[1/z]

fReg[s_] := Normal@Series[g[s], {s, 0, 70}]

r = MiniMaxApproximation[g[z], {z, {2^-23, 1}, 21, 2}, WorkingPrecision -> 53]
```
The series
```text
num_expr = R('1.2533141373155002546825431555831963990977483071423532')+R('11.581464036148641246581323440149428509260059307098616') * z+\
R('26.789405266295715810558337396548064861031194327368438')* z**2 + R('7.3649417349531436967188010018370139100180573487266502') * z**3 \
- R('1.7096604692508150493340167930493663713673704318166940') * z**4+R('1.0804255276094591411415638947116978390123141265771834') *z**5\
+ R('-1.0469297880831626083930679006829252269737874877936028') * z**6+R('1.3002369515375729922365615519495020332408157597613544') * z**7 \
+ R('-1.8716868006869984282123524738168938690121343397007868') * z**8+R('2.9031694730522518237131899238517881330616803775751557') * z**9\
+R('-4.5732446918809058964898044936311043412823198783007045') * z**10+R('6.9647249644245100963438249211100614638042076172726253') * z**11\
+R('-9.8478668366535772967541646404708100660109176278818808') * z**12+R('12.509168035842933091158686315206698671310735894796194')*  z ** 13\
+ R('-13.889426789762708072171011520013585808954445304261549') * z**14+R('13.160468639020739474423102361556011097677163142977284') * z**15\
+ R('-10.397245062168499847721119740813368446763822759405222') * z**16+R('6.6791470271098701770545337257192390539281302865779416')* z**17\
+ R('-3.3836343614489495447896172378322198542839494659295979') * z**18+R('1.2963434524119789098189498423127080263478530079657214') * z**19\
+ R('-0.35208923494132138684480968071673794226763563154951347') * z**20+R('0.060290645206928280777309624029614864219397241979017862') * z**21
+ R('-0.0048861904873423508267690490483453010423241530137713923') * z**22

den_expr = R('1') + R('8.8656713459366445640627099502526076307634900011454201') * z + R('18.167413600340225608199066162434028713917004023936959') * z**2

```
**/
#[inline]
fn k1_asympt(x: f64) -> f64 {
    let recip = DoubleDouble::from_recip(x);
    let e = i0_exp(x * 0.5);
    let r_sqrt = DoubleDouble::from_sqrt(x);

    const P: [(u64, u64); 22] = [
        (0xbc97b64b25244ad6, 0x3ff40d931ff62706),
        (0x3c9a1688cf334c34, 0x4026664d0189916c),
        (0xbcc7e17378873c58, 0x40390c1509d1e208),
        (0x3cbc4973f3c5d670, 0x401b3d731b85f2ac),
        (0x3c933564234461b2, 0xbff8fb80e04d69cf),
        (0x3c82498fbf50fe90, 0x3fef167068283d7e),
        (0x3c2f4dc704d7c59d, 0xbfed83372b6a7c7e),
        (0x3c68ea0b746e2ab0, 0x3ff1d3dae3d65ebf),
        (0xbc84fb7c307a8dcd, 0xbff8b994d39a4ae3),
        (0x3ca4e1e7f9af4249, 0x4002412f816c0a67),
        (0x3caa7719cafdd0de, 0xc00b00cc9475d532),
        (0x3cb818338b2aa469, 0x4013076e044ea75a),
        (0xbc92861cf8b4f4df, 0xc0188742789f0db3),
        (0xbcb865210d89e7f7, 0x401bf69171ba9a8f),
        (0x3cb3810304f23021, 0xc01b66ed0b78cfee),
        (0x3cb2fa6dcaaf5bca, 0x40167b5b0a43b01c),
        (0xbc9e7b2d015e8486, 0xc00e0e964b108bd6),
        (0xbc62ad01e6afddb9, 0x3fffb5bc4b1b46da),
        (0x3c13823e47c3c736, 0xbfe94b1ee6b990cb),
        (0xbc5091596e7d3991, 0x3fcc94fd84869a53),
        (0x3c3a10fdceafdcdc, 0xbfa456e7b9a102bb),
        (0xbbe97d1f7bfa68cf, 0x3f6b5d423e9b273a),
    ];

    let p0 = f_polyeval22(
        recip,
        DoubleDouble::from_bit_pair(P[0]),
        DoubleDouble::from_bit_pair(P[1]),
        DoubleDouble::from_bit_pair(P[2]),
        DoubleDouble::from_bit_pair(P[3]),
        DoubleDouble::from_bit_pair(P[4]),
        DoubleDouble::from_bit_pair(P[5]),
        DoubleDouble::from_bit_pair(P[6]),
        DoubleDouble::from_bit_pair(P[7]),
        DoubleDouble::from_bit_pair(P[8]),
        DoubleDouble::from_bit_pair(P[9]),
        DoubleDouble::from_bit_pair(P[10]),
        DoubleDouble::from_bit_pair(P[11]),
        DoubleDouble::from_bit_pair(P[12]),
        DoubleDouble::from_bit_pair(P[13]),
        DoubleDouble::from_bit_pair(P[14]),
        DoubleDouble::from_bit_pair(P[15]),
        DoubleDouble::from_bit_pair(P[16]),
        DoubleDouble::from_bit_pair(P[17]),
        DoubleDouble::from_bit_pair(P[18]),
        DoubleDouble::from_bit_pair(P[19]),
        DoubleDouble::from_bit_pair(P[20]),
        DoubleDouble::from_bit_pair(P[21]),
    );

    let mut q = DoubleDouble::from_bit_pair((0xbcd78a8ce1b4e974, 0x4030e43b1173fa56));
    q = DoubleDouble::mul_add(
        q,
        recip,
        DoubleDouble::from_bit_pair((0xbcc847253ca99018, 0x40211f4f6157a41d)),
    );
    q = DoubleDouble::mul_add_f64(q, recip, f64::from_bits(0x3ff0000000000000));

    let v = DoubleDouble::div(p0, q);
    let r = DoubleDouble::div(v, e * r_sqrt * e);

    let err = f_fmla(
        r.hi,
        f64::from_bits(0x3c60000000000000), // 2^-57
        f64::from_bits(0x3c10000000000000), // 2^-62
    );
    let ub = r.hi + (r.lo + err);
    let lb = r.hi + (r.lo - err);
    if ub != lb {
        return k1_asympt_hard(x);
    }
    r.to_f64()
}

/**
Generated by Wolfram:
```text
<< FunctionApproximations`
f[x_] := Sqrt[x] Exp[x] BesselK[1, x]
g[z_] := f[1/z]

fReg[s_] := Normal@Series[g[s], {s, 0, 90}]

r = MiniMaxApproximation[g[z], {z, {2^-23, 1}, 24, 3}, WorkingPrecision -> 60]
```

```python
numerator = R('1.25331413731550025121927282596719844509479202793343596081659') +\
R('18.1818090161526720433056569747515378514630065641474625981014') * z +\
R('84.0855308746947292990799286725009214888553349291444140660290') * z**2 +\
R('132.162365847816287381044941763470816630213463932078033240530') * z**3 +\
R('31.9226828676807644427339830797111472671372253510654795457144') * z**4 -\
R('6.55623468795551166985522946155773048406824737626044144871123') * z**5 +\
R('3.64813675168036216906333568031312648295984485156865732661545') * z**6 -\
R('3.09030769637907344287942414087144964418293628150871931378425') * z**7 +\
R('3.33386714181870765265684414717381480790589837093132038622747') * z**8 -\
R('4.15687509263144029996675634899519359876131330016912373219164') * z**9 +\
R('5.60168914866900308748492484302032849837652046070612403967421') * z**10 -\
R('7.74291668942803388088083086973577048045063529950220281515214') * z**11 +\
R('10.5182262725991203709271656214076983768897923475996777051390') * z**12 -\
R('13.5518704730239628983580812834313096171287792612392021785201') * z**13 +\
R('16.0756083151279354000726821130423204859625333785087750230562') * z**14 -\
R('17.1182809371435152577732819586801963660470433312228530729236') * z**15 +\
R('16.0011674817017251371454409588419263278595693541151351864594') * z**16 -\
R('12.8547201574813259441920569388886787474288736479192357179690') * z**17 +\
R('8.68486423302827636350417719675952558710622021096282923741362') * z**18 -\
R('4.81495629672548049328728197240253790905840957706216081678396') * z**19 +\
R('2.12449990606466172694756924925308136959633096070223397071382') * z**20 -\
R('0.715212551284256837854636761671218490310899339423567201757157') * z**21 +\
R('0.172124045885696320404925832426986371003557469307869009460024') * z**22 -\
R('0.0263211214093518652429028546952441524255064789438934805575092') * z**23 +\
R('0.00191881798965983631900602829845520749528121908949103217489718') * z**24

denominator = R('1') + R('14.1319847014545522182859614183872365899195687422959996226761') * z +\
R('61.9082401087861097277037880035610946116051551544124209371854') * z**2 +\
R('83.7882740830616503293819459356854835383521002198942068063204') * z**3
```
**/
#[cold]
#[inline(never)]
fn k1_asympt_hard(x: f64) -> f64 {
    static P: [DyadicFloat128; 25] = [
        DyadicFloat128 {
            sign: DyadicSign::Pos,
            exponent: -127,
            mantissa: 0xa06c98ff_b1382cb2_d9370160_5e6ea89f_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Pos,
            exponent: -123,
            mantissa: 0x91745849_13f30526_b9684bb8_b82db634_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Pos,
            exponent: -121,
            mantissa: 0xa82bcab3_eb3969ee_a0ffd314_025498c0_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Pos,
            exponent: -120,
            mantissa: 0x842990ce_e65bb768_2cdb5294_33e86bed_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Pos,
            exponent: -123,
            mantissa: 0xff61a78e_2a25a3b4_9f71f915_d7c1e4c8_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Neg,
            exponent: -125,
            mantissa: 0xd1ccacb0_356e14f3_a98c4b92_8b72f7b1_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Pos,
            exponent: -126,
            mantissa: 0xe97b1291_f3618b0c_0d3411fa_f169ba9e_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Neg,
            exponent: -126,
            mantissa: 0xc5c799ee_a19d2c31_6486d81a_bca8ac1d_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Pos,
            exponent: -126,
            mantissa: 0xd55e1449_d480814e_8b253c43_79c06b41_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Neg,
            exponent: -125,
            mantissa: 0x85051eea_0d169527_cb42b4f8_ffddb220_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Pos,
            exponent: -125,
            mantissa: 0xb3410999_fc86c28c_54c48d6b_0a9f0f62_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Neg,
            exponent: -125,
            mantissa: 0xf7c5f938_97df2fc8_36ff7aaa_8ce0f82f_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Pos,
            exponent: -124,
            mantissa: 0xa84aa7a1_cbdbf795_89732a19_627c8250_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Neg,
            exponent: -124,
            mantissa: 0xd8d47622_14416aac_8373673d_205486bf_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Pos,
            exponent: -123,
            mantissa: 0x809ad888_463dbc73_a198fa66_18e2f886_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Neg,
            exponent: -123,
            mantissa: 0x88f23d46_a62c4c2f_ee73c6f4_fe1270f9_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Pos,
            exponent: -123,
            mantissa: 0x80026418_bdd21999_a053bb54_80ad68b0_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Neg,
            exponent: -124,
            mantissa: 0xcdacef0b_39d4072e_62cb2e02_59904a60_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Pos,
            exponent: -124,
            mantissa: 0x8af53432_b0e7b5f8_52998d50_34116635_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Neg,
            exponent: -125,
            mantissa: 0x9a141f3a_435de2a3_5ded5d67_b49abf5b_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Pos,
            exponent: -126,
            mantissa: 0x87f7ce74_39c7613b_f70affe2_85a9ccb0_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Neg,
            exponent: -128,
            mantissa: 0xb7182b75_74610daf_73b44ee1_05507e76_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Pos,
            exponent: -130,
            mantissa: 0xb041492f_bcb106a3_45c37c92_ccce62f3_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Neg,
            exponent: -133,
            mantissa: 0xd79f6474_b5f54a25_809eddef_c039070e_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Pos,
            exponent: -137,
            mantissa: 0xfb80d906_6efafa9e_48980d5f_dcaa6be0_u128,
        },
    ];

    static Q: [DyadicFloat128; 4] = [
        DyadicFloat128 {
            sign: DyadicSign::Pos,
            exponent: -127,
            mantissa: 0x80000000_00000000_00000000_00000000_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Pos,
            exponent: -124,
            mantissa: 0xe21c9bfd_851d2f88_1bd6b241_3b20d294_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Pos,
            exponent: -122,
            mantissa: 0xf7a209b1_f09b77a1_449c5c22_148f5c3d_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Pos,
            exponent: -121,
            mantissa: 0xa79398a9_1e117f7a_c139750d_52237350_u128,
        },
    ];

    let recip = DyadicFloat128::accurate_reciprocal(x);
    let e = rational128_exp(x);
    let r_sqrt = bessel_rsqrt_hard(x, recip);

    let mut p0 = P[24];
    for i in (0..24).rev() {
        p0 = recip * p0 + P[i];
    }

    let mut q0 = Q[3];
    for i in (0..3).rev() {
        q0 = recip * q0 + Q[i];
    }

    let v = p0 * q0.reciprocal();
    let r = v * (e.reciprocal() * r_sqrt);
    r.fast_as_f64()
}

#[cfg(test)]
mod tests {
    use super::*;
    #[test]
    fn test_k1() {
        assert_eq!(f_k1(0.643), 1.184534109892725);
        assert_eq!(f_k1(0.964), 0.6402280656771248);
        assert_eq!(f_k1(2.964), 0.04192888446074039);
        assert_eq!(f_k1(8.43), 9.824733212831289e-5);
        assert_eq!(f_k1(16.43), 2.3142404075259965e-8);
        assert_eq!(f_k1(423.43), 7.793648638470207e-186);
        assert_eq!(f_k1(0.), f64::INFINITY);
        assert_eq!(f_k1(-0.), f64::INFINITY);
        assert!(f_k1(-0.5).is_nan());
        assert!(f_k1(f64::NEG_INFINITY).is_nan());
        assert_eq!(f_k1(f64::INFINITY), 0.);
    }
}
