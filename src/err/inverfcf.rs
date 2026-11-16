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
use crate::err::{f_erfcinv, inverff::erfinv_core};

/// Complementary inverse error function
///
/// Max ulp 0.5
pub fn f_erfcinvf(x: f32) -> f32 {
    let ix = x.to_bits();
    let ux = ix.wrapping_shl(1);
    if ix >= 0x4000_0000u32 || ux == 0 {
        if x.is_infinite() {
            return f32::INFINITY;
        }
        if ux == 0 {
            return f32::INFINITY;
        }

        if ix == 0x3f80_0000u32 {
            return 0.;
        }
        // x > 2
        if ix == 0x4000_0000u32 {
            // x == 2.
            return f32::NEG_INFINITY;
        }
        return f32::NAN; // x == NaN, x < 0
    }

    // erfcinv(2 - x) = -erfcinv(x). Keep the complement for the logarithm
    // in the tail, where 1 - (1 - x) can lose significant bits.
    let (q, sign) = if x > 1. {
        (2. - x as f64, -1.)
    } else {
        (x as f64, 1.)
    };
    if q < f32::EPSILON as f64 * 0.5 {
        // Beyond the range of the single-precision erfinv approximation.
        return (f_erfcinv(q) as f32).copysign(sign);
    }
    let z = 1. - q;
    erfinv_core(z, (z as f32).to_bits(), q, sign)
}

#[cfg(test)]
mod tests {
    use super::f_erfcinvf;

    #[test]
    fn m_test() {
        assert_eq!(f_erfcinvf(2.), f32::NEG_INFINITY);
        assert!(f_erfcinvf(-1.).is_nan());
        assert_eq!(f_erfcinvf(0.), f32::INFINITY);
        assert!(f_erfcinvf(2.1).is_nan());
        assert_eq!(f_erfcinvf(0.5), 0.47693628);
        assert_eq!(f_erfcinvf(1.5), -0.47693628);
        assert_eq!(f_erfcinvf(1.95).to_bits(), 0xbfb1654e);
        assert_eq!(f_erfcinvf(f32::from_bits(0x3fffffff)).to_bits(), 0xc06f9c68);
        assert_eq!(f_erfcinvf(0.002), 2.1851242);
        assert_eq!(f_erfcinvf(1e-30).to_bits(), 0x410260bb);
        assert_eq!(f_erfcinvf(f32::from_bits(1)).to_bits(), 0x4120513e);
        assert_eq!(f_erfcinvf(1.002), -0.0017724329);
        assert!(f_erfcinvf(f32::NAN).is_nan());
    }
}
