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
use crate::common::f_fmla;
use crate::dekker::Dekker;
use crate::f_sin;

#[inline]
fn eval_sincf7(x: f32, c: &[u64; 7]) -> f32 {
    let x = x as f64;
    let x2 = x * x; // Exact

    let pw0 = f_fmla(x2, f64::from_bits(c[6]), f64::from_bits(c[5]));
    let pw1 = f_fmla(x2, pw0, f64::from_bits(c[4]));
    let pw2 = f_fmla(x2, pw1, f64::from_bits(c[3]));
    let pw3 = f_fmla(x2, pw2, f64::from_bits(c[2]));
    let pw4 = f_fmla(x2, pw3, f64::from_bits(c[1]));

    let p = f_fmla(x2, pw4, f64::from_bits(c[0]));
    p as f32
}

#[inline]
fn eval_sincf9(x: f32, c: &[u64; 9]) -> f32 {
    let x = x as f64;
    let x2 = x * x; // Exact

    let pwz0 = f_fmla(x2, f64::from_bits(c[8]), f64::from_bits(c[7]));
    let pw0 = f_fmla(x2, pwz0, f64::from_bits(c[6]));
    let pwz1 = f_fmla(x2, pw0, f64::from_bits(c[5]));
    let pw1 = f_fmla(x2, pwz1, f64::from_bits(c[4]));
    let pw2 = f_fmla(x2, pw1, f64::from_bits(c[3]));
    let pw3 = f_fmla(x2, pw2, f64::from_bits(c[2]));
    let pw4 = f_fmla(x2, pw3, f64::from_bits(c[1]));

    let p = f_fmla(x2, pw4, f64::from_bits(c[0]));
    p as f32
}

/// Computes sinc(x)
///
/// Max ULP 0.5
#[inline]
pub fn f_sincf(x: f32) -> f32 {
    if !x.is_finite() {
        return f32::NAN;
    }
    let x_abs = f32::from_bits(x.to_bits() & 0x7fff_ffff);
    if x_abs.to_bits() == 0 {
        return 1.0;
    }
    if x_abs < 0.25 {
        // |x| < 0.25
        const C: [u64; 4] = [
            0x3ff0000000000000,
            0xbfc55555550fde98,
            0x3f81110f70c67a6b,
            0xbf29f67b484037dc,
        ];
        let x = x as f64;
        let x2 = x * x; // Exact

        let p0 = f_fmla(x2, f64::from_bits(C[3]), f64::from_bits(C[2]));
        let p1 = f_fmla(x2, p0, f64::from_bits(C[1]));

        let p = f_fmla(x2, p1, f64::from_bits(C[0]));
        p as f32
    } else if x_abs < 2.85 {
        let coeffs = if x_abs < 0.75 {
            // |x| < 0.75
            static C: [u64; 7] = [
                0x3ff0000000000000,
                0xbfc5555555555549,
                0x3f8111111110fb74,
                0xbf2a01a019ca9bef,
                0x3ec71de362235231,
                0xbe5ae5f05f9152b3,
                0x3de5dc861a6a3c21,
            ];
            C
        } else if x_abs < 1.57 {
            // |x| < 1.57
            static C: [u64; 7] = [
                0x3ff0000000000000,
                0xbfc555555552065f,
                0x3f811111101242a8,
                0xbf2a019fa19b5214,
                0x3ec71dc6dd03cb48,
                0xbe5adedff8876f83,
                0x3de5184fd3c8edcd,
            ];
            C
        } else if x_abs < 2.1 {
            // |x| < 2.1
            static C: [u64; 7] = [
                0x3ff0000000000000,
                0xbfc5555552f41ecc,
                0x3f811110cc75224e,
                0xbf2a019352ab374e,
                0x3ec71ca0ef521a8f,
                0xbe5ac2a4d637f264,
                0x3de3f98370d40d85,
            ];
            C
        } else if x_abs < 2.45 {
            // |x| < 2.45
            static C: [u64; 7] = [
                0x3ff0000000000000,
                0xbfc5555538690083,
                0x3f81110ef6356d5f,
                0xbf2a015ecfd56289,
                0x3ec719aba6efd127,
                0xbe5a9796a6f71249,
                0x3de2fc7fd7fb04cb,
            ];
            C
        } else if x_abs < 2.61 {
            // |x| < 2.61
            static C: [u64; 7] = [
                0x3ff0000000000000,
                0xbfc55554f2b7b098,
                0x3f81110b443194fe,
                0xbf2a010e5d4a91bb,
                0x3ec7163df22ec949,
                0xbe5a72189ba0518c,
                0x3de2582235c59c7e,
            ];
            C
        } else {
            // |x| < 2.85
            static C: [u64; 7] = [
                0x3ff0000000000000,
                0xbfc5555462d7b375,
                0x3f811104d1bc6c26,
                0xbf2a0097e13d9374,
                0x3ec711fbc177f989,
                0xbe5a4adc1696833b,
                0x3de1c754f3c7bcfd,
            ];
            C
        };
        eval_sincf7(x, &coeffs)
    } else if x_abs < 6.45 {
        let coeffs = if x_abs < 3.2 {
            // |x| < 3.2
            static C: [u64; 9] = [
                0x3ff0000000000000,
                0xbfc5555555127248,
                0x3f8111110c5eb459,
                0xbf2a019f8ea35321,
                0x3ec71ddac14032fa,
                0xbe5ae5947e7233a3,
                0x3de609a0091b4cae,
                0xbd6a61f88bf0a150,
                0x3ce4ad7fcf8eec12,
            ];
            C
        } else if x_abs < 3.8 {
            // |x| < 3.8
            static C: [u64; 9] = [
                0x3ff0000000000000,
                0xbfc555554e8fd6d4,
                0x3f811110c91f046f,
                0xbf2a019ae1b235f3,
                0x3ec71dac4558dc3a,
                0xbe5ae366c2e87cbe,
                0x3de5f9d92c10a7d2,
                0xbd69e242b973cbc8,
                0x3ce2efbf486349c8,
            ];
            C
        } else if x_abs < 4.35 {
            // |x| < 4.35
            static C: [u64; 9] = [
                0x3ff0000000000000,
                0xbfc555550b573a0b,
                0x3f81110ed113dbce,
                0xbf2a018171ed5559,
                0x3ec71cf4ba66d7e7,
                0xbe5add288ca0d9c5,
                0x3de5d90a5ee880c8,
                0xbd69218ab5b64a6e,
                0x3ce1079b8a566e01,
            ];
            C
        } else if x_abs < 4.95 {
            // |x| < 4.95
            static C: [u64; 9] = [
                0x3ff0000000000000,
                0xbfc555532512b621,
                0x3f811103fc172ccd,
                0xbf2a01172538ca3f,
                0x3ec71aaed77d8fdc,
                0xbe5ace2a1549e291,
                0x3de59d6f695536c9,
                0xbd6819157f94659d,
                0x3cde1ce2e0fc7ef2,
            ];
            C
        } else if x_abs < 5.45 {
            // |x| < 5.45
            static C: [u64; 9] = [
                0x3ff0000000000000,
                0xbfc5554989dc90c6,
                0x3f8110d8b368b2ac,
                0xbf29ffc7d28c2f38,
                0x3ec715073f44efdc,
                0xbe5ab0c79a6a2d12,
                0x3de54186b79adaeb,
                0xbd66d8a684f63e33,
                0x3cda5c3db744ec2c,
            ];
            C
        } else if x_abs < 6.0 {
            // |x| < 6.0
            static C: [u64; 9] = [
                0x3ff0000000000000,
                0xbfc55523b744174b,
                0x3f81104d59f2c2fe,
                0xbf29fc55b96eff63,
                0x3ec708e1fa558332,
                0xbe5a7d48d3621025,
                0x3de4be3641d9988e,
                0xbd6563ae0d0ec9f1,
                0x3cd6cdfde196a174,
            ];
            C
        } else {
            // |x| < 6.45
            static C: [u64; 9] = [
                0x3ff0000000000000,
                0xbfc554adc87d6fe2,
                0x3f810ee1282ef843,
                0xbf29f4ca50711fbf,
                0x3ec6f29d7bec7cd6,
                0xbe5a2e46e410a5dd,
                0x3de415b67176a6a4,
                0xbd63d3abd4f19458,
                0x3cd39e93da3a418a,
            ];
            C
        };
        eval_sincf9(x, &coeffs)
    } else {
        let s = f_sin(x as f64);
        Dekker::from_exact_div(s, x as f64).to_f64() as f32
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    #[test]
    fn test_f_sincf() {
        assert_eq!(f_sincf(0.0), 1.0);
        assert_eq!(f_sincf(0.2), 0.99334663);
    }
}
