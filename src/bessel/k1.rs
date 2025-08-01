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
use crate::bessel::i0::i0_exp;
use crate::bessel::y0::log_dd;
use crate::double_double::DoubleDouble;
use crate::polyeval::{f_horner_polyeval13, f_polyeval11, f_polyeval22};

/// Modified Bessel of the second kind order 1
///
/// Max ULP 0.52
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
        //706.5025494880165
        return 0.;
    }

    if xb <= 0x3ff0000000000000 {
        return k1_small(x);
    }

    k1_asympt(x)
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
    let z0 = DoubleDouble::mul_f64_add(p, x, z);
    z0.to_f64()
}

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
    r.to_f64()
}

#[cfg(test)]
mod tests {
    use super::*;
    #[test]
    fn test_k1() {
        assert_eq!(f_k1(0.11), 8.935335323830637);
        assert_eq!(f_k1(0.643), 1.184534109892725);
        assert_eq!(f_k1(0.964), 0.6402280656771248);
        assert_eq!(f_k1(2.964), 0.04192888446074039);
        assert_eq!(f_k1(423.43), 7.793648638470207e-186);
        assert_eq!(f_k1(0.), f64::INFINITY);
        assert_eq!(f_k1(-0.), f64::INFINITY);
        assert!(f_k1(-0.5).is_nan());
        assert!(f_k1(f64::NEG_INFINITY).is_nan());
        assert_eq!(f_k1(f64::INFINITY), 0.);
    }
}
