/*
 * // Copyright (c) Radzivon Bartoshyk 6/2025. All rights reserved.
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
use crate::bits::EXP_MASK;
use crate::common::f_fmla;
use std::ops::{Add, Mul};

#[derive(Copy, Clone, Ord, PartialOrd, Eq, PartialEq)]
pub(crate) enum RationalSign {
    Pos,
    Neg,
}

const BITS: u32 = 128;

#[derive(Copy, Clone)]
pub(crate) struct RationalFloat128 {
    pub(crate) sign: RationalSign,
    pub(crate) exponent: i16,
    pub(crate) mantissa: u128,
}

#[inline]
fn f64_from_parts(sign: RationalSign, exp: u64, mantissa: u64) -> f64 {
    let r_sign = (if sign == RationalSign::Pos {
        0u64
    } else {
        1u64
    })
    .wrapping_shl(63);
    let r_exp = exp.wrapping_shl(52);
    f64::from_bits(r_sign | r_exp | mantissa)
}

#[inline]
fn mulhi_u128(a: u128, b: u128) -> u128 {
    let a_lo = a as u64;
    let a_hi = (a >> 64) as u64;
    let b_lo = b as u64;
    let b_hi = (b >> 64) as u64;

    let lo_lo = a_lo as u128 * b_lo as u128;
    let hi_lo = a_hi as u128 * b_lo as u128;
    let lo_hi = a_lo as u128 * b_hi as u128;
    let hi_hi = a_hi as u128 * b_hi as u128;

    // cross terms with carry propagation
    let mid1 = hi_lo + (lo_lo >> 64);
    let (mid2, carry) = mid1.overflowing_add(lo_hi);
    hi_hi + (mid2 >> 64) + if carry { 1 } else { 0 }
}

#[inline]
const fn explicit_exponent(x: f64) -> i16 {
    let exp = ((x.to_bits() >> 52) & ((1u64 << 11) - 1u64)) as i16 - 1023;
    if x == 0. {
        return 0;
    } else if x.is_subnormal() {
        const EXP_BIAS: u64 = (1u64 << (11 - 1u64)) - 1u64;
        return 1i16 - EXP_BIAS as i16;
    }
    exp
}

#[inline]
const fn explicit_mantissa(x: f64) -> u64 {
    const MASK: u64 = (1u64 << 52) - 1;
    let sig_bits = x.to_bits() & MASK;
    if x.is_subnormal() || x == 0. {
        return sig_bits;
    }
    (1u64 << 52) | sig_bits
}

impl RationalFloat128 {
    #[inline]
    fn zero() -> Self {
        Self {
            sign: RationalSign::Pos,
            exponent: 0,
            mantissa: 0,
        }
    }

    #[inline]
    pub(crate) const fn new_from_f64(x: f64) -> Self {
        let sign = if x.is_sign_negative() {
            RationalSign::Neg
        } else {
            RationalSign::Pos
        };
        let exponent = explicit_exponent(x) - 52;
        let mantissa = explicit_mantissa(x) as u128;
        let mut new_val = Self {
            sign,
            exponent,
            mantissa,
        };
        new_val.normalize();
        new_val
    }

    #[inline]
    fn shift_right(&mut self, amount: u32) {
        if amount < BITS {
            self.exponent += amount as i16;
            self.mantissa = self.mantissa.wrapping_shr(amount);
        } else {
            self.exponent = 0;
            self.mantissa = 0;
        }
    }

    #[inline]
    fn shift_left(&mut self, amount: u32) {
        if amount < BITS {
            self.exponent -= amount as i16;
            self.mantissa = self.mantissa.wrapping_shl(amount);
        } else {
            self.exponent = 0;
            self.mantissa = 0;
        }
    }

    #[inline]
    const fn normalize(&mut self) {
        if self.mantissa != 0 {
            let shift_length = self.mantissa.leading_zeros();
            self.exponent -= shift_length as i16;
            self.mantissa = self.mantissa.wrapping_shl(shift_length);
        }
    }

    #[inline]
    pub(crate) fn quick_add(&self, rhs: &Self) -> Self {
        if self.mantissa == 0 {
            return *rhs;
        }
        if rhs.mantissa == 0 {
            return *self;
        }
        let mut a = *self;
        let mut b = *rhs;

        // Align exponents
        if a.exponent > b.exponent {
            b.shift_right((a.exponent - b.exponent) as u32);
        } else if b.exponent > a.exponent {
            a.shift_right((b.exponent - a.exponent) as u32);
        }

        let mut result = RationalFloat128::zero();

        if a.sign == b.sign {
            // Addition
            result.sign = a.sign;
            result.exponent = a.exponent;
            result.mantissa = a.mantissa;
            let (sum, is_overflow) = result.mantissa.overflowing_add(b.mantissa);
            result.mantissa = sum;
            if is_overflow {
                // Mantissa addition overflow.
                result.shift_right(1);
                result.mantissa |= 1u128 << 127;
            }
            // Result is already normalized.
            return result;
        }

        // Subtraction
        if a.mantissa >= b.mantissa {
            result.sign = a.sign;
            result.exponent = a.exponent;
            result.mantissa = a.mantissa - b.mantissa;
        } else {
            result.sign = b.sign;
            result.exponent = b.exponent;
            result.mantissa = b.mantissa - a.mantissa;
        }

        result.normalize();
        result
    }

    #[inline]
    pub(crate) fn quick_mul(&self, rhs: &Self) -> Self {
        let mut result = RationalFloat128 {
            sign: if self.sign != rhs.sign {
                RationalSign::Neg
            } else {
                RationalSign::Pos
            },
            exponent: self.exponent + rhs.exponent + BITS as i16,
            mantissa: 0,
        };

        if !(self.mantissa == 0 || rhs.mantissa == 0) {
            result.mantissa = mulhi_u128(self.mantissa, rhs.mantissa);
            // Check the leading bit directly, should be faster than using clz in
            // normalize().
            if result.mantissa >> 127 == 0 {
                result.shift_left(1);
            }
        } else {
            result.mantissa = 0;
        }
        result
    }

    #[inline]
    pub(crate) fn fast_as_f64(&self) -> f64 {
        if self.mantissa == 0 {
            return if self.sign == RationalSign::Pos {
                0.
            } else {
                -0.0
            };
        }

        // Assume that it is normalized, and output is also normal.
        const PRECISION: u32 = 52 + 1;

        // SIG_MASK - FRACTION_MASK
        const SIG_MASK: u64 = (1u64 << 52) - 1;
        const FRACTION_MASK: u64 = (1u64 << 52) - 1;
        const IMPLICIT_MASK: u64 = SIG_MASK - FRACTION_MASK;
        const EXP_BIAS: u64 = (1u64 << (11 - 1u64)) - 1u64;

        let mut exp_hi = self.exponent as i64 + ((BITS - 1) as i64 + EXP_BIAS as i64);

        if exp_hi > 2 * EXP_BIAS as i64 {
            // Results overflow.
            let d_hi = f64_from_parts(self.sign, 2 * EXP_BIAS, IMPLICIT_MASK);
            // volatile prevents constant propagation that would result in infinity
            // always being returned no matter the current rounding mode.
            let two = 2.0f64;
            let r = two * d_hi;
            return r;
        }

        let mut denorm = false;
        let mut shift = BITS - PRECISION;
        if exp_hi <= 0 {
            // Output is denormal.
            denorm = true;
            shift = (BITS - PRECISION) + (1 - exp_hi) as u32;

            exp_hi = EXP_BIAS as i64;
        }

        let exp_lo = exp_hi - PRECISION as i64 - 1;
        //        let m_hi =
        //             shift >= MantissaType::BITS ? MantissaType(0) : mantissa >> shift;
        let m_hi = if shift >= BITS {
            0
        } else {
            self.mantissa >> shift
        };

        let d_hi = f64_from_parts(
            self.sign,
            exp_hi as u64,
            (m_hi as u64 & SIG_MASK) | IMPLICIT_MASK,
        );

        let round_mask = if shift > BITS {
            0
        } else {
            1u128 << (shift - 1)
        };
        let sticky_mask = round_mask - 1u128;

        let round_bit = (self.mantissa & round_mask) != 0;
        let sticky_bit = (self.mantissa & sticky_mask) != 0;
        let round_and_sticky = round_bit as i32 * 2 + sticky_bit as i32;

        let d_lo: f64;

        if exp_lo <= 0 {
            // d_lo is denormal, but the output is normal.
            let scale_up_exponent = 1 - exp_lo;
            let scale_up_factor = f64_from_parts(
                RationalSign::Pos,
                EXP_BIAS + scale_up_exponent as u64,
                IMPLICIT_MASK,
            );
            let scale_down_factor = f64_from_parts(
                RationalSign::Pos,
                EXP_BIAS - scale_up_exponent as u64,
                IMPLICIT_MASK,
            );

            d_lo = f64_from_parts(
                self.sign,
                (exp_lo + scale_up_exponent) as u64,
                IMPLICIT_MASK,
            );

            return f_fmla(d_lo, round_and_sticky as f64, d_hi * scale_up_factor)
                * scale_down_factor;
        }

        d_lo = f64_from_parts(self.sign, exp_lo as u64, IMPLICIT_MASK);

        const SIG_LEN: u64 = 53;
        // Still correct without FMA instructions if `d_lo` is not underflow.
        let r = f_fmla(d_lo, round_and_sticky as f64, d_hi);

        if denorm {
            // Exponent before rounding is in denormal range, simply clear the
            // exponent field.
            let clear_exp: u64 = (exp_hi as u64) << SIG_LEN;
            let mut r_bits: u64 = r.to_bits() - clear_exp;

            if r_bits & EXP_MASK == 0 {
                // Output is denormal after rounding, clear the implicit bit for 80-bit
                // long double.
                r_bits -= IMPLICIT_MASK;
            }

            return f64::from_bits(r_bits);
        }

        r
    }
}

impl Add<RationalFloat128> for RationalFloat128 {
    type Output = RationalFloat128;
    #[inline]
    fn add(self, rhs: RationalFloat128) -> Self::Output {
        self.quick_add(&rhs)
    }
}

impl Mul<RationalFloat128> for RationalFloat128 {
    type Output = RationalFloat128;
    #[inline]
    fn mul(self, rhs: RationalFloat128) -> Self::Output {
        self.quick_mul(&rhs)
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_rational_float() {
        let ones = RationalFloat128 {
            sign: RationalSign::Pos,
            exponent: -127,
            mantissa: 0x80000000_00000000_00000000_00000000_u128,
        };
        let cvt = ones.fast_as_f64();
        assert_eq!(cvt, 1.0);

        let minus_0_5 = RationalFloat128 {
            sign: RationalSign::Neg,
            exponent: -128,
            mantissa: 0x80000000_00000000_00000000_00000000_u128,
        };
        let cvt0 = minus_0_5.fast_as_f64();
        assert_eq!(cvt0, -1.0 / 2.0);

        let minus_1_f4 = RationalFloat128 {
            sign: RationalSign::Neg,
            exponent: -132,
            mantissa: 0xaaaaaaaa_aaaaaaaa_aaaaaaaa_aaaaaaab_u128,
        };
        let cvt0 = minus_1_f4.fast_as_f64();
        assert_eq!(cvt0, -1.0 / 24.0);

        let minus_1_f8 = RationalFloat128 {
            sign: RationalSign::Pos,
            exponent: -143,
            mantissa: 0xd00d00d0_0d00d00d_00d00d00_d00d00d0_u128,
        };
        let cvt0 = minus_1_f8.fast_as_f64();
        assert_eq!(cvt0, 1.0 / 40320.0);
    }

    #[test]
    fn rational_float_add() {
        let ones = RationalFloat128 {
            sign: RationalSign::Pos,
            exponent: -127,
            mantissa: 0x80000000_00000000_00000000_00000000_u128,
        };

        let cvt = ones.fast_as_f64();
        assert_eq!(cvt, 1.0);

        let minus_0_5 = RationalFloat128 {
            sign: RationalSign::Neg,
            exponent: -128,
            mantissa: 0x80000000_00000000_00000000_00000000_u128,
        };
        let cvt0 = ones.quick_add(&minus_0_5).fast_as_f64();
        assert_eq!(cvt0, 0.5);
    }

    #[test]
    fn rational_float_mul() {
        let ones = RationalFloat128 {
            sign: RationalSign::Pos,
            exponent: -127,
            mantissa: 0x80000000_00000000_00000000_00000000_u128,
        };

        let cvt = ones.fast_as_f64();
        assert_eq!(cvt, 1.0);

        let minus_0_5 = RationalFloat128 {
            sign: RationalSign::Neg,
            exponent: -128,
            mantissa: 0x80000000_00000000_00000000_00000000_u128,
        };
        let product = ones.quick_mul(&minus_0_5);
        let cvt0 = product.fast_as_f64();
        assert_eq!(cvt0, -0.5);
    }

    #[test]
    fn rational_round_trip() {
        let z00 = 0.0;
        let zvt00 = RationalFloat128::new_from_f64(z00);
        let b00 = zvt00.fast_as_f64();
        assert_eq!(b00, z00);

        let z0 = 1.0;
        let zvt0 = RationalFloat128::new_from_f64(z0);
        let b0 = zvt0.fast_as_f64();
        assert_eq!(b0, z0);

        let z1 = 0.5;
        let zvt1 = RationalFloat128::new_from_f64(z1);
        let b1 = zvt1.fast_as_f64();
        assert_eq!(b1, z1);

        let z2 = -0.5;
        let zvt2 = RationalFloat128::new_from_f64(z2);
        let b2 = zvt2.fast_as_f64();
        assert_eq!(b2, z2);

        let z3 = -532322.54324324232;
        let zvt3 = RationalFloat128::new_from_f64(z3);
        let b3 = zvt3.fast_as_f64();
        assert_eq!(b3, z3);
    }
}
