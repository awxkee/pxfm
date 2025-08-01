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
use crate::bessel::i0::{bessel_rsqrt_hard, i0_0_to_3p5, i0_exp};
use crate::bessel::y0::log_dd;
use crate::common::f_fmla;
use crate::double_double::DoubleDouble;
use crate::dyadic_float::{DyadicFloat128, DyadicSign};
use crate::exponents::rational128_exp;
use crate::polyeval::{f_polyeval10, f_polyeval22, f_polyeval24};

/// Modified Bessel of the second kind order 0
///
/// Max ULP 0.5
pub fn f_k0(x: f64) -> f64 {
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

    if xb >= 0x40862e42fefa39f0u64 {
        // 709.7827128933841
        return 0.;
    }

    if xb <= 0x3ff0000000000000 {
        return k0_small_dd(x);
    }

    k0_asympt(x)
}

/**
K0(x) + log(x) * I0(x) = P(x^2)
hence
K0(x) = P(x^2) - log(x)*I0(x)

Series:
```python
euler_gamma = R(euler)

lg  = R(2).log()
lg4 = R(4).log()
lg64 = R(64).log()
lg4096 = R(4096).log()

expr2 = -euler_gamma + lg + R(1/4) *  z**2 *(1 - euler_gamma + lg) + (\
 z**12 *(R(49) - R(20) * euler_gamma + R(20) *lg))/R('42467328000') + (\
 z**10 *(R(137) - R(60) *euler_gamma + R(60) *lg))/R('884736000') + (\
 z**14 *(R(363) - R(140)* euler_gamma + R(140) *lg))/R('58265174016000') + (\
 z**16 *(R(761) - R(280) *euler_gamma + R(280) *lg))/R('29831769096192000') + (\
 z**18 *(R(7129) - R(2520) *euler_gamma + R(2520) *lg))/R('86989438684495872000') +\
 R(1/128)* z**4 *(R(3) - R(2) * euler_gamma + lg4) + (\
 z**6 *(R(11) - R(6) * euler_gamma + R(64).log()))/R(13824) + (\
 z**8 *(R(25) - R(12) * euler_gamma + R(4096).log()))/R(1769472)
```
**/
#[inline]
fn k0_small_dd(x: f64) -> f64 {
    let dx = DoubleDouble::from_exact_mult(x, x);
    const C: [(u64, u64); 10] = [
        (0x3c1be095d05c0a81, 0x3fbdadb014541eb2),
        (0x3c6037c12ba0b815, 0x3fd1dadb014541eb),
        (0x3c2037c12ba0b815, 0x3f99dadb014541eb),
        (0xbbe07eec845045e4, 0x3f4bb90e85debf56),
        (0x3b830f4c5f3df300, 0x3eef4747696cf839),
        (0x3b252fcaeee73fd1, 0x3e85d6b13b0d88ca),
        (0xbab3127e7a5114b0, 0x3e14c2b6e8177e1a),
        (0x3a309e5ad5685b51, 0x3d9ca0246d234e72),
        (0xb9a183e5b5dac36d, 0x3d1df24eb119a2f9),
        (0xb93ad95e64dfe5fc, 0x3c9896d55d330a18),
    ];
    let r = f_polyeval10(
        dx,
        DoubleDouble::from_bit_pair(C[0]),
        DoubleDouble::from_bit_pair(C[1]),
        DoubleDouble::from_bit_pair(C[2]),
        DoubleDouble::from_bit_pair(C[3]),
        DoubleDouble::from_bit_pair(C[4]),
        DoubleDouble::from_bit_pair(C[5]),
        DoubleDouble::from_bit_pair(C[6]),
        DoubleDouble::from_bit_pair(C[7]),
        DoubleDouble::from_bit_pair(C[8]),
        DoubleDouble::from_bit_pair(C[9]),
    );

    let vi_log = log_dd(x);
    let vi = i0_0_to_3p5(x);
    let r = DoubleDouble::mul_add(vi_log, -vi, r);
    r.to_f64()
}

/**
Generated in Wolfram

Computes sqrt(x)*exp(x)*K0(x)=Pn(1/x)/Qm(1/x)
hence
K0(x) = Pn(1/x)/Qm(1/x) / (sqrt(x) * exp(x))

```text
<< FunctionApproximations`
f[x_] := Sqrt[x] Exp[x] BesselK[0, x]
g[z_] := f[1/z]
r = MiniMaxApproximation[g[z], {z, {0.0000000000001, 1}, 21, 2}, WorkingPrecision -> 53]
```
**/
#[inline]
fn k0_asympt(x: f64) -> f64 {
    let recip = DoubleDouble::from_recip(x);
    let e = i0_exp(x * 0.5);
    let r_sqrt = DoubleDouble::from_sqrt(x);

    const P: [(u64, u64); 22] = [
        (0xbc9d19c3268f2679, 0x3ff40d931ff62706),
        (0xbcc32aceef776564, 0x4025455477b71ce2),
        (0xbcc1e5f3825a3b79, 0x40342caa45b3f6cc),
        (0x3c9970d1b035befa, 0xc0001963e879e9e9),
        (0xbc3317093756f755, 0x3feb6f010c6e9487),
        (0x3c8f65ccda8c1c67, 0xbfe49dbc440d6092),
        (0xbc69f65925c49b4f, 0x3fe57501e3f535b2),
        (0xbc8e9316972a8ad6, 0xbfeb5502cfaa3477),
        (0xbc8979206c11ead1, 0x3ff39530dd790d44),
        (0x3c88b4670722debb, 0xbffd8415c7d89a3a),
        (0xbcac97a762878f3f, 0x40061be0148be74d),
        (0xbca15d40d04f8aec, 0xc00f644bdae417da),
        (0x3caa7a75b76a92bc, 0x40144f79adfca031),
        (0xbcaf1e71de046ce3, 0xc017307e5d5be768),
        (0x3caae602742239d5, 0x4016b8ee8945c869),
        (0x3cb082d4076ff824, 0xc0129f1e45bd5d46),
        (0x3c815dd750f7ecf8, 0x4008d9db400c3166),
        (0x3c8572dfa95b4bd2, 0xbffa289e51c1dfb2),
        (0x3c74c6ed3f460bd8, 0x3fe4d010b9fe5e6e),
        (0xbc538cc8d46e2dc0, 0xbfc77495b50c2615),
        (0x3c24d4b496d939d8, 0x3fa0a516540369b0),
        (0x3bebe3b8583f8e21, 0xbf6654a059eb3ea1),
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

    let mut q = DoubleDouble::from_bit_pair((0xbcdd7c0915a96602, 0x40311a5a656052a4));
    q = DoubleDouble::mul_add(
        q,
        recip,
        DoubleDouble::from_bit_pair((0x3ca7816fede5dee4, 0x402138bea47588ef)),
    );
    q = DoubleDouble::mul_add_f64(q, recip, f64::from_bits(0x3ff0000000000000));

    let v = DoubleDouble::div(p0, q);
    let r = DoubleDouble::div(v, e * r_sqrt * e);

    let err = f_fmla(
        r.hi,
        f64::from_bits(0x3c62611186bae67f), // 2^-56.8
        f64::from_bits(0x3c094804b79e25f6), // 2^-62.34
    );
    let ub = r.hi + (r.lo + err);
    let lb = r.hi + (r.lo - err);
    if ub != lb {
        return k0_asympt_hard(x);
    }
    r.to_f64()
}

/**
Generated in Wolfram

Computes sqrt(x)*exp(x)*K0(x)=Pn(1/x)/Qm(1/x)
hence
K0(x) = Pn(1/x)/Qm(1/x) / (sqrt(x) * exp(x))

```text
<< FunctionApproximations`
f[x_] := Sqrt[x] Exp[x] BesselK[0, x]
g[z_] := f[1/z]
r = MiniMaxApproximation[g[z], {z, {2^-23, 1}, 23, 3}, WorkingPrecision -> 70]
```
**/
#[inline(never)]
#[cold]
fn k0_asympt_hard(x: f64) -> f64 {
    const P: [DyadicFloat128; 24] = [
        DyadicFloat128 {
            sign: DyadicSign::Pos,
            exponent: -127,
            mantissa: 0xa06c98ff_b1382cb2_7682afba_ae5bb2b7_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Pos,
            exponent: -123,
            mantissa: 0x88b1db65_eb09d512_35872892_7f54225d_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Pos,
            exponent: -121,
            mantissa: 0x8eb2357c_0039f396_d784c96a_b9538621_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Pos,
            exponent: -121,
            mantissa: 0xb0aac77f_ed145b45_75091ce7_75cda5e6_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Neg,
            exponent: -124,
            mantissa: 0x80254758_e25b2929_18edbb0d_24138310_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Pos,
            exponent: -126,
            mantissa: 0xc32c14bb_2673e613_52f0c7ea_83dadd24_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Neg,
            exponent: -126,
            mantissa: 0x81952ab0_1f7ae615_6d7ff168_de6d58cc_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Pos,
            exponent: -127,
            mantissa: 0xebfdd89a_c5fa59ad_c44d7abd_f29e685a_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Neg,
            exponent: -126,
            mantissa: 0x8299f9b6_2bfe99b6_41b82cf9_514c63e4_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Pos,
            exponent: -126,
            mantissa: 0xa2551b66_6b1b5db9_a6248623_7111f50e_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Neg,
            exponent: -126,
            mantissa: 0xd55c20dd_211041b7_20794f7f_ed21b1bb_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Pos,
            exponent: -125,
            mantissa: 0x8d2996e1_6f71a1ee_da21445a_06956e47_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Neg,
            exponent: -125,
            mantissa: 0xb47dbb70_949cea92_e708cda6_db36b180_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Pos,
            exponent: -125,
            mantissa: 0xd76e0b61_b4bfd9af_280e0080_5c2b1f86_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Neg,
            exponent: -125,
            mantissa: 0xe91bd545_21a1b801_7f738ebd_6f0ea106_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Pos,
            exponent: -125,
            mantissa: 0xdee8f756_023b4dc9_864d1d41_88bcad9d_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Neg,
            exponent: -125,
            mantissa: 0xb7ff9528_605693e7_175d68bd_803ec79e_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Pos,
            exponent: -125,
            mantissa: 0x8010bf6c_0fa7a2b2_035848ef_d59c6a75_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Neg,
            exponent: -126,
            mantissa: 0x927e4eba_7637c76c_d21ec7b2_2bccc0ea_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Pos,
            exponent: -127,
            mantissa: 0x856fbd4b_b8401151_7a8070e8_612ca817_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Neg,
            exponent: -129,
            mantissa: 0xb976216c_4f63c8e2_c4377163_4d9aba67_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Pos,
            exponent: -131,
            mantissa: 0xb830a5d3_619c8cd6_3223574d_9821436b_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Neg,
            exponent: -134,
            mantissa: 0xe84dd64e_ea661595_22dbfee3_1269fff2_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Pos,
            exponent: -137,
            mantissa: 0x8b8d2366_5712a95a_959dfaac_1ebc8619_u128,
        },
    ];

    const Q: [DyadicFloat128; 4] = [
        DyadicFloat128 {
            sign: DyadicSign::Pos,
            exponent: -127,
            mantissa: 0x80000000_00000000_00000000_00000000_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Pos,
            exponent: -124,
            mantissa: 0xdc221dd1_efa4ae26_809727db_87150bfe_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Pos,
            exponent: -122,
            mantissa: 0xea4ed875_b128a6a5_a4be0f75_b6b7f4d1_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Pos,
            exponent: -121,
            mantissa: 0x99d0e07e_7585770a_19ca394b_985f77c2_u128,
        },
    ];

    let recip = DyadicFloat128::accurate_reciprocal(x);
    let e = rational128_exp(x);
    let r_sqrt = bessel_rsqrt_hard(x, recip);

    let p0 = f_polyeval24(
        recip, P[0], P[1], P[2], P[3], P[4], P[5], P[6], P[7], P[8], P[9], P[10], P[11], P[12],
        P[13], P[14], P[15], P[16], P[17], P[18], P[19], P[20], P[21], P[22], P[23],
    );

    let mut q = Q[3];
    q = q * recip + Q[2];
    q = q * recip + Q[1];
    q = q * recip + Q[0];

    let v = p0 * q.reciprocal();
    let r = v * e.reciprocal() * r_sqrt;
    r.fast_as_f64()
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_k0() {
        assert_eq!(f_k0(0.11), 2.3332678776741127);
        assert_eq!(f_k0(0.643), 0.7241025575342853);
        assert_eq!(f_k0(0.964), 0.4433737413379138);
        assert_eq!(f_k0(2.964), 0.03621679838808167);
        assert_eq!(f_k0(423.43), 7.784461905543397e-186);
        assert_eq!(f_k0(0.), f64::INFINITY);
        assert_eq!(f_k0(-0.), f64::INFINITY);
        assert!(f_k0(-0.5).is_nan());
        assert!(f_k0(f64::NEG_INFINITY).is_nan());
        assert_eq!(f_k0(f64::INFINITY), 0.);
    }
}
