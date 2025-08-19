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
use crate::bessel::i0f::i0f_small;
use crate::common::f_fmla;
use crate::f_exp;
use crate::logs::simple_fast_log;
use crate::polyeval::{f_polyeval8, f_polyeval16};

/// Modified Bessel of the second kind order 0
///
/// Max ULP 0.5
///
/// This method have exactly one exception which is not correctly rounded with FMA.
pub fn f_k0f(x: f32) -> f32 {
    if x < 0. {
        return f32::NAN;
    }

    if (x.to_bits() & 0x0007_ffff) == 0 {
        if x == 0. {
            return f32::INFINITY;
        }
        if x.is_infinite() {
            return if x.is_sign_positive() { 0. } else { f32::NAN };
        }
        if x.is_nan() {
            return x + x;
        }
    }

    let xb = x.to_bits();

    if xb >= 0x42cbc4fbu32 {
        // 101.88473
        return 0.;
    }

    if xb <= 0x3f800000u32 {
        // 1.0
        return k0f_small(x);
    }

    k0f_asympt(x)
}

/**
K0(x) + log(x) * I0(x) = P(x^2)
hence
K0(x) = P(x^2) - log(x)*I0(x)

Series:
```python
euler_gamma = R(euler)

lg = R(2).log()
lg4 = R(4).log()
lg8 = R(8).log()

k0f_small = -euler_gamma + lg + 1/4 * z**2 * (R(1) - euler_gamma + lg) + ( \
 z**12 * (R(49) - 20 * euler_gamma + 20 * lg))/R(42467328000) + (\
 z**10 * (R(137) - 60 * euler_gamma + 60 * lg))/(884736000) + (\
 z**14 * (R(363) - 140 * euler_gamma + 140 * lg))/(58265174016000) + (\
 z**16 * (761 - 280  * euler_gamma  + 280 * lg))/R(29831769096192000) + \
 1/128 * z**4 * (3 - R(2) * euler_gamma + lg4) + (\
 z**6 * (11 - 6 * euler_gamma + R(64).log()))/R(13824) + (\
 z**8 * (25 - 12 * euler_gamma + R(4096).log()))/R(1769472)
```
**/
#[inline]
fn k0f_small(x: f32) -> f32 {
    let v_log = simple_fast_log(x as f64);
    let i0 = i0f_small(x);

    let dx = x as f64;

    let p = f_polyeval8(
        dx * dx,
        f64::from_bits(0x3fbdadb014541eb2),
        f64::from_bits(0x3fd1dadb014541eb),
        f64::from_bits(0x3f99dadb014541eb),
        f64::from_bits(0x3f4bb90e85debf56),
        f64::from_bits(0x3eef4747696cf839),
        f64::from_bits(0x3e85d6b13b0d88ca),
        f64::from_bits(0x3e14c2b6e8177e1a),
        f64::from_bits(0x3d9ca0246d234e72),
    );
    let c = f_fmla(-i0, v_log, p);
    c as f32
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
r = MiniMaxApproximation[g[z], {z, {0.0000000000001, 1}, 16, 2}, WorkingPrecision -> 53]
```
**/
#[inline]
fn k0f_asympt(x: f32) -> f32 {
    let dx = x as f64;
    let recip = 1. / dx;
    let e = f_exp(dx);
    let r_sqrt = dx.sqrt();

    let p0 = f_polyeval16(
        recip,
        f64::from_bits(0x3ff40d931ff626ed),
        f64::from_bits(0x402072fbcec226c2),
        f64::from_bits(0x40275186037f1723),
        f64::from_bits(0xbff146ffff2cc0e5),
        f64::from_bits(0x3fda7e3ca23841ac),
        f64::from_bits(0xbfd13462ad384f99),
        f64::from_bits(0x3fcd3539bf027d79),
        f64::from_bits(0xbfcc06b166c3fcd9),
        f64::from_bits(0x3fcb5b5a47741483),
        f64::from_bits(0xbfc8fa3c2bc30d81),
        f64::from_bits(0x3fc3e6b86fa6ba6c),
        f64::from_bits(0xbfba07af249e75b5),
        f64::from_bits(0x3faa4296b62b2eb2),
        f64::from_bits(0xbf92e8497eaaf439),
        f64::from_bits(0x3f71247eb3aab51a),
        f64::from_bits(0xbf3d3e0258904b11),
    );

    let mut q = f64::from_bits(0x402422f9b9281f24);
    q = f_fmla(q, recip, f64::from_bits(0x401abfc1f49353c1));
    q = f_fmla(q, recip, f64::from_bits(0x3ff0000000000000));
    let v = p0 / q;
    let pp = v / (e * r_sqrt);
    pp as f32
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_k0f() {
        assert_eq!(f_k0f(2.034804e-5), 10.918458);
        assert_eq!(f_k0f(0.010260499), 4.695535);
        assert_eq!(f_k0f(0.3260499), 1.2965646);
        assert_eq!(f_k0f(0.72341), 0.636511734);
        assert_eq!(f_k0f(0.), f32::INFINITY);
        assert_eq!(f_k0f(-0.), f32::INFINITY);
        assert!(f_k0f(-0.5).is_nan());
        assert!(f_k0f(f32::NEG_INFINITY).is_nan());
        assert_eq!(f_k0f(f32::INFINITY), 0.);
    }
}
