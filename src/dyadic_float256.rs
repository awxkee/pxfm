// /*
//  * // Copyright (c) Radzivon Bartoshyk 6/2025. All rights reserved.
//  * //
//  * // Redistribution and use in source and binary forms, with or without modification,
//  * // are permitted provided that the following conditions are met:
//  * //
//  * // 1.  Redistributions of source code must retain the above copyright notice, this
//  * // list of conditions and the following disclaimer.
//  * //
//  * // 2.  Redistributions in binary form must reproduce the above copyright notice,
//  * // this list of conditions and the following disclaimer in the documentation
//  * // and/or other materials provided with the distribution.
//  * //
//  * // 3.  Neither the name of the copyright holder nor the names of its
//  * // contributors may be used to endorse or promote products derived from
//  * // this software without specific prior written permission.
//  * //
//  * // THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
//  * // AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
//  * // IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
//  * // DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE
//  * // FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL
//  * // DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR
//  * // SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER
//  * // CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY,
//  * // OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
//  * // OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
//  */
// use crate::bits::EXP_MASK;
// use crate::common::f_fmla;
// use crate::dyadic_float::{DyadicSign, f64_from_parts};
// use crate::u256::{mulhi_u256, u256};
// use std::ops::{Add, Mul, Shr, Sub};
//
// const BITS: u32 = 256;
//
// #[derive(Copy, Clone, Debug)]
// pub(crate) struct DyadicFloat256 {
//     pub(crate) sign: DyadicSign,
//     pub(crate) exponent: i16,
//     pub(crate) mantissa: u256,
// }
//
// #[inline]
// const fn explicit_exponent(x: f64) -> i16 {
//     let exp = ((x.to_bits() >> 52) & ((1u64 << 11) - 1u64)) as i16 - 1023;
//     if x == 0. {
//         return 0;
//     } else if x.is_subnormal() {
//         const EXP_BIAS: u64 = (1u64 << (11 - 1u64)) - 1u64;
//         return 1i16 - EXP_BIAS as i16;
//     }
//     exp
// }
//
// #[inline]
// const fn explicit_mantissa(x: f64) -> u64 {
//     const MASK: u64 = (1u64 << 52) - 1;
//     let sig_bits = x.to_bits() & MASK;
//     if x.is_subnormal() || x == 0. {
//         return sig_bits;
//     }
//     (1u64 << 52) | sig_bits
// }
//
// impl DyadicFloat256 {
//     #[inline]
//     pub(crate) const fn zero() -> Self {
//         Self {
//             sign: DyadicSign::Pos,
//             exponent: 0,
//             mantissa: u256::ZERO,
//         }
//     }
//
//     #[inline]
//     pub(crate) fn new_from_f64(x: f64) -> Self {
//         let sign = if x.is_sign_negative() {
//             DyadicSign::Neg
//         } else {
//             DyadicSign::Pos
//         };
//         let exponent = explicit_exponent(x) - 52;
//         let mantissa = u256::from_u64(explicit_mantissa(x));
//         let mut new_val = Self {
//             sign,
//             exponent,
//             mantissa,
//         };
//         new_val.normalize();
//         new_val
//     }
//
//     #[inline]
//     pub(crate) fn new(sign: DyadicSign, exponent: i16, mantissa: u256) -> Self {
//         let mut new_item = DyadicFloat256 {
//             sign,
//             exponent,
//             mantissa,
//         };
//         new_item.normalize();
//         new_item
//     }
//
//     #[inline]
//     pub(crate) fn accurate_reciprocal(a: f64) -> Self {
//         // we convert 4/a and divide by 4 to avoid a spurious underflow
//         let mut r = DyadicFloat256::new_from_f64(4.0 / a); /* accurate to about 53 bits */
//         r.exponent -= 2;
//         /* we use Newton's iteration: r -> r + r*(1-a*r) */
//         let ba = DyadicFloat256::new_from_f64(-a);
//         let mut q = ba * r;
//         const F256_ONE: DyadicFloat256 = DyadicFloat256 {
//             sign: DyadicSign::Pos,
//             exponent: -(BITS as i16 - 1),
//             mantissa: u256 {
//                 hi: 0x8000_0000_0000_0000_0000_0000_0000_0000_u128,
//                 lo: 0x0u128,
//             },
//         };
//         q = F256_ONE + q;
//         q = r * q;
//
//         let mut q2 = ba * q;
//         q2 = F256_ONE + q2;
//         q2 = q * q2;
//
//         let mut q3 = ba * q2;
//         q3 = F256_ONE + q3;
//         q3 = q2 * q3;
//         r + q3
//     }
//
//     #[inline]
//     pub(crate) fn from_div_f64(a: f64, b: f64) -> Self {
//         let reciprocal = DyadicFloat256::accurate_reciprocal(b);
//         let da = DyadicFloat256::new_from_f64(a);
//         reciprocal * da
//     }
//
//     // /// Multiply self by integer scalar `b`.
//     // /// Returns a new normalized DyadicFloat256.
//     // #[inline]
//     // pub(crate) fn mul_int64(&self, b: i64) -> DyadicFloat256 {
//     //     if b == 0 {
//     //         return DyadicFloat256::zero();
//     //     }
//     //
//     //     let abs_b = b.unsigned_abs();
//     //     let sign = if (b < 0) ^ (self.sign == DyadicSign::Neg) {
//     //         DyadicSign::Neg
//     //     } else {
//     //         DyadicSign::Pos
//     //     };
//     //
//     //     let mut hi_prod = (self.mantissa >> 64).wrapping_mul(abs_b as u128);
//     //     let m = hi_prod.leading_zeros();
//     //     hi_prod <<= m;
//     //
//     //     let mut lo_prod = (self.mantissa & 0xffff_ffff_ffff_ffff).wrapping_mul(abs_b as u128);
//     //     lo_prod = (lo_prod << (m - 1)) >> 63;
//     //
//     //     let (mut product, overflow) = hi_prod.overflowing_add(lo_prod);
//     //
//     //     let mut result = DyadicFloat256 {
//     //         sign,
//     //         exponent: self.exponent + 64 - m as i16,
//     //         mantissa: product,
//     //     };
//     //
//     //     if overflow {
//     //         // Overflow means an implicit bit in the 129th place, which we shift down.
//     //         product += product & 0x1;
//     //         result.mantissa = (product >> 1) | (1u128 << 127);
//     //         result.shift_right(1);
//     //     }
//     //
//     //     result.normalize();
//     //     result
//     // }
//
//     #[inline]
//     fn shift_right(&mut self, amount: u32) {
//         if amount < BITS {
//             self.exponent += amount as i16;
//             self.mantissa = self.mantissa.shr(amount);
//         } else {
//             self.exponent = 0;
//             self.mantissa = u256::ZERO;
//         }
//     }
//
//     #[inline]
//     fn shift_left(&mut self, amount: u32) {
//         if amount < BITS {
//             self.exponent -= amount as i16;
//             self.mantissa = self.mantissa.wrapping_shl(amount);
//         } else {
//             self.exponent = 0;
//             self.mantissa = u256::ZERO;
//         }
//     }
//
//     // Don't forget to call if manually created
//     #[inline]
//     pub(crate) fn normalize(&mut self) {
//         if self.mantissa != u256::ZERO {
//             let shift_length = self.mantissa.leading_zeros();
//             self.exponent -= shift_length as i16;
//             self.mantissa = self.mantissa.wrapping_shl(shift_length);
//         }
//     }
//
//     #[inline]
//     pub(crate) fn negated(&self) -> Self {
//         Self {
//             sign: self.sign.negate(),
//             exponent: self.exponent,
//             mantissa: self.mantissa,
//         }
//     }
//
//     #[inline]
//     pub(crate) fn quick_sub(&self, rhs: &Self) -> Self {
//         self.quick_add(&rhs.negated())
//     }
//
//     #[inline]
//     pub(crate) fn quick_add(&self, rhs: &Self) -> Self {
//         if self.mantissa == u256::ZERO {
//             return *rhs;
//         }
//         if rhs.mantissa == u256::ZERO {
//             return *self;
//         }
//         let mut a = *self;
//         let mut b = *rhs;
//
//         let exp_diff = a.exponent.wrapping_sub(b.exponent);
//
//         // If exponent difference is too large, b is negligible
//         if exp_diff.abs() >= BITS as i16 {
//             return if a.sign == b.sign {
//                 // Adding very small number to large: return a
//                 return if a.exponent > b.exponent { a } else { b };
//             } else if a.exponent > b.exponent {
//                 a
//             } else {
//                 b
//             };
//         }
//
//         // Align exponents
//         if a.exponent > b.exponent {
//             b.shift_right((a.exponent - b.exponent) as u32);
//         } else if b.exponent > a.exponent {
//             a.shift_right((b.exponent - a.exponent) as u32);
//         }
//
//         let mut result = DyadicFloat256::zero();
//
//         if a.sign == b.sign {
//             // Addition
//             result.sign = a.sign;
//             result.exponent = a.exponent;
//             result.mantissa = a.mantissa;
//             let (sum, is_overflow) = result.mantissa.overflowing_add(b.mantissa);
//             result.mantissa = sum;
//             if is_overflow {
//                 // Mantissa addition overflow.
//                 result.shift_right(1);
//                 result.mantissa |= u256::ONE.wrapping_shl(255);
//             }
//             // Result is already normalized.
//             return result;
//         }
//
//         // Subtraction
//         if a.mantissa >= b.mantissa {
//             result.sign = a.sign;
//             result.exponent = a.exponent;
//             result.mantissa = a.mantissa.wrapping_sub(b.mantissa);
//         } else {
//             result.sign = b.sign;
//             result.exponent = b.exponent;
//             result.mantissa = b.mantissa.wrapping_sub(a.mantissa);
//         }
//
//         result.normalize();
//         result
//     }
//
//     #[inline]
//     pub(crate) fn quick_mul(&self, rhs: &Self) -> Self {
//         let mut result = DyadicFloat256 {
//             sign: if self.sign != rhs.sign {
//                 DyadicSign::Neg
//             } else {
//                 DyadicSign::Pos
//             },
//             exponent: self.exponent + rhs.exponent + BITS as i16,
//             mantissa: u256::ZERO,
//         };
//
//         if !(self.mantissa == u256::ZERO || rhs.mantissa == u256::ZERO) {
//             result.mantissa = mulhi_u256(self.mantissa, rhs.mantissa);
//             // Check the leading bit directly, should be faster than using clz in
//             // normalize().
//             if result.mantissa >> 255 == u256::ZERO {
//                 result.shift_left(1);
//             }
//         } else {
//             result.mantissa = u256::ZERO;
//         }
//         result
//     }
//
//     #[inline]
//     pub(crate) fn fast_as_f64(&self) -> f64 {
//         if self.mantissa == u256::ZERO {
//             return if self.sign == DyadicSign::Pos {
//                 0.
//             } else {
//                 -0.0
//             };
//         }
//
//         // Assume that it is normalized, and output is also normal.
//         const PRECISION: u32 = 52 + 1;
//
//         // SIG_MASK - FRACTION_MASK
//         const SIG_MASK: u64 = (1u64 << 52) - 1;
//         const FRACTION_MASK: u64 = (1u64 << 52) - 1;
//         const IMPLICIT_MASK: u64 = SIG_MASK - FRACTION_MASK;
//         const EXP_BIAS: u64 = (1u64 << (11 - 1u64)) - 1u64;
//
//         let mut exp_hi = self.exponent as i32 + ((BITS - 1) as i32 + EXP_BIAS as i32);
//
//         if exp_hi > 2 * EXP_BIAS as i32 {
//             // Results overflow.
//             let d_hi = f64_from_parts(self.sign, 2 * EXP_BIAS, IMPLICIT_MASK);
//             // volatile prevents constant propagation that would result in infinity
//             // always being returned no matter the current rounding mode.
//             let two = 2.0f64;
//             let r = two * d_hi;
//             return r;
//         }
//
//         let mut denorm = false;
//         let mut shift = BITS - PRECISION;
//         if exp_hi <= 0 {
//             // Output is denormal.
//             denorm = true;
//             shift = (BITS - PRECISION) + (1 - exp_hi) as u32;
//
//             exp_hi = EXP_BIAS as i32;
//         }
//
//         let exp_lo = exp_hi.wrapping_sub(PRECISION as i32).wrapping_sub(1);
//
//         let m_hi = if shift >= BITS {
//             0u64
//         } else {
//             (self.mantissa >> shift).to_u64()
//         };
//
//         let d_hi = f64_from_parts(self.sign, exp_hi as u64, (m_hi & SIG_MASK) | IMPLICIT_MASK);
//
//         let round_mask = if shift > BITS {
//             u256::ZERO
//         } else {
//             u256::ONE.wrapping_shl(shift.wrapping_sub(1))
//         };
//         let sticky_mask = round_mask.wrapping_sub(u256::ONE);
//
//         let round_bit = (self.mantissa & round_mask) != u256::ZERO;
//         let sticky_bit = (self.mantissa & sticky_mask) != u256::ZERO;
//         let round_and_sticky = round_bit as i32 * 2 + sticky_bit as i32;
//
//         let d_lo: f64;
//
//         if exp_lo <= 0 {
//             // d_lo is denormal, but the output is normal.
//             let scale_up_exponent = 1 - exp_lo;
//             let scale_up_factor = f64_from_parts(
//                 DyadicSign::Pos,
//                 EXP_BIAS + scale_up_exponent as u64,
//                 IMPLICIT_MASK,
//             );
//             let scale_down_factor = f64_from_parts(
//                 DyadicSign::Pos,
//                 EXP_BIAS - scale_up_exponent as u64,
//                 IMPLICIT_MASK,
//             );
//
//             d_lo = f64_from_parts(
//                 self.sign,
//                 (exp_lo + scale_up_exponent) as u64,
//                 IMPLICIT_MASK,
//             );
//
//             return f_fmla(d_lo, round_and_sticky as f64, d_hi * scale_up_factor)
//                 * scale_down_factor;
//         }
//
//         d_lo = f64_from_parts(self.sign, exp_lo as u64, IMPLICIT_MASK);
//
//         // Still correct without FMA instructions if `d_lo` is not underflow.
//         let r = f_fmla(d_lo, round_and_sticky as f64, d_hi);
//
//         if denorm {
//             const SIG_LEN: u64 = 52;
//             // Exponent before rounding is in denormal range, simply clear the
//             // exponent field.
//             let clear_exp: u64 = (exp_hi as u64) << SIG_LEN;
//             let mut r_bits: u64 = r.to_bits() - clear_exp;
//
//             if r_bits & EXP_MASK == 0 {
//                 // Output is denormal after rounding, clear the implicit bit for 80-bit
//                 // long double.
//                 r_bits -= IMPLICIT_MASK;
//             }
//
//             return f64::from_bits(r_bits);
//         }
//
//         r
//     }
//
//     // Approximate reciprocal - given a nonzero `a`, make a good approximation to 1/a.
//     // The method is Newton-Raphson iteration, based on quick_mul.
//     #[inline]
//     pub(crate) fn reciprocal(self) -> DyadicFloat256 {
//         // Computes the reciprocal using Newton-Raphson iteration:
//         // Given an approximation x ≈ 1/a, we refine via:
//         //     x' = x * (2 - a * x)
//         // This squares the error term: if ax ≈ 1 - e, then ax' ≈ 1 - e².
//
//         let guess = 1. / self.fast_as_f64();
//         let mut x = DyadicFloat256::new_from_f64(guess);
//
//         // The constant 2, which we'll need in every iteration
//         let twos = DyadicFloat256 {
//             sign: DyadicSign::Pos,
//             exponent: -(BITS as i16 - 2),
//             mantissa: u256 {
//                 hi: 0x80000000_00000000_00000000_00000000_u128,
//                 lo: 0x0u128,
//             },
//         };
//
//         x = x * (twos - (self * x));
//         x = x * (twos - (self * x));
//         x = x * (twos - (self * x));
//         x = x * (twos - (self * x));
//         x = x * (twos - (self * x));
//         x
//     }
//
//     // // Approximate reciprocal - given a nonzero `a`, make a good approximation to 1/a.
//     // // The method is Newton-Raphson iteration, based on quick_mul.
//     // // *This is very crude guess*
//     // #[inline]
//     // fn approximate_reciprocal(&self) -> DyadicFloat256 {
//     //     // Given an approximation x to 1/a, a better one is x' = x(2-ax).
//     //     //
//     //     // You can derive this by using the Newton-Raphson formula with the function
//     //     // f(x) = 1/x - a. But another way to see that it works is to say: suppose
//     //     // that ax = 1-e for some small error e. Then ax' = ax(2-ax) = (1-e)(1+e) =
//     //     // 1-e^2. So the error in x' is the square of the error in x, i.e. the number
//     //     // of correct bits in x' is double the number in x.
//     //
//     //     // An initial approximation to the reciprocal
//     //     let mut x = DyadicFloat256 {
//     //         sign: DyadicSign::Pos,
//     //         exponent: -32 - self.exponent - BITS as i16,
//     //         mantissa: self.mantissa >> (BITS - 32),
//     //     };
//     //     x.normalize();
//     //
//     //     // The constant 2, which we'll need in every iteration
//     //     let two = DyadicFloat256::new(DyadicSign::Pos, 1, 1);
//     //
//     //     // We expect at least 31 correct bits from our 32-bit starting approximation
//     //     let mut ok_bits = 31usize;
//     //
//     //     // The number of good bits doubles in each iteration, except that rounding
//     //     // errors introduce a little extra each time. Subtract a bit from our
//     //     // accuracy assessment to account for that.
//     //     while ok_bits < BITS as usize {
//     //         x = x * (two - (*self * x));
//     //         ok_bits = 2 * ok_bits - 1;
//     //     }
//     //
//     //     x
//     // }
// }
//
// impl Add<DyadicFloat256> for DyadicFloat256 {
//     type Output = DyadicFloat256;
//     #[inline]
//     fn add(self, rhs: DyadicFloat256) -> Self::Output {
//         self.quick_add(&rhs)
//     }
// }
//
// impl DyadicFloat256 {
//     #[inline]
//     pub(crate) fn biased_exponent(&self) -> i16 {
//         self.exponent + (BITS as i16 - 1)
//     }
//
//     // #[inline]
//     // pub(crate) fn trunc_to_i64(&self) -> i64 {
//     //     if self.exponent <= -(BITS as i16) {
//     //         // Absolute value of x is greater than equal to 0.5 but less than 1.
//     //         return 0;
//     //     }
//     //     let hi = self.mantissa >> 64;
//     //     let norm_exp = self.biased_exponent();
//     //     if norm_exp > 63 {
//     //         return if self.sign == DyadicSign::Neg {
//     //             i64::MIN
//     //         } else {
//     //             i64::MAX
//     //         };
//     //     }
//     //     let r: i64 = (hi >> (63 - norm_exp)) as i64;
//     //
//     //     if self.sign == DyadicSign::Neg { -r } else { r }
//     // }
//
//     #[inline]
//     pub(crate) fn round_to_nearest(&self) -> DyadicFloat256 {
//         if self.exponent == -(BITS as i16) {
//             // Absolute value of x is greater than equal to 0.5 but less than 1.
//             return DyadicFloat256 {
//                 sign: self.sign,
//                 exponent: -(BITS as i16 - 1),
//                 mantissa: u256 {
//                     hi: 0x80000000_00000000_00000000_00000000_u128,
//                     lo: 0x0u128,
//                 },
//             };
//         }
//         if self.exponent <= -((BITS + 1) as i16) {
//             // Absolute value of x is greater than equal to 0.5 but less than 1.
//             return DyadicFloat256 {
//                 sign: self.sign,
//                 exponent: 0,
//                 mantissa: u256::ZERO,
//             };
//         }
//         const FRACTION_LENGTH: u32 = BITS - 1;
//         let trim_size =
//             (FRACTION_LENGTH as i64).wrapping_sub(self.exponent as i64 + (BITS - 1) as i64) as u32;
//         let half_bit_set =
//             self.mantissa & (u256::ONE.wrapping_shl(trim_size.wrapping_sub(1))) != u256::ZERO;
//         let trunc_u: u256 = (self.mantissa >> trim_size).wrapping_shl(trim_size);
//         if trunc_u == self.mantissa {
//             return *self;
//         }
//
//         let truncated = DyadicFloat256::new(self.sign, self.exponent, trunc_u);
//
//         if !half_bit_set {
//             // Franctional part is less than 0.5 so round value is the
//             // same as the trunc value.
//             truncated
//         } else if self.sign == DyadicSign::Neg {
//             let ones = DyadicFloat256 {
//                 sign: DyadicSign::Pos,
//                 exponent: -(BITS as i16 - 1),
//                 mantissa: u256 {
//                     hi: 0x8000_0000_0000_0000_0000_0000_0000_0000_u128,
//                     lo: 0x0u128,
//                 },
//             };
//             truncated - ones
//         } else {
//             let ones = DyadicFloat256 {
//                 sign: DyadicSign::Pos,
//                 exponent: -(BITS as i16 - 1),
//                 mantissa: u256 {
//                     hi: 0x8000_0000_0000_0000_0000_0000_0000_0000_u128,
//                     lo: 0x0u128,
//                 },
//             };
//             truncated + ones
//         }
//     }
//
//     #[inline]
//     pub(crate) fn round_to_nearest_f64(&self) -> f64 {
//         self.round_to_nearest().fast_as_f64()
//     }
// }
//
// impl Sub<DyadicFloat256> for DyadicFloat256 {
//     type Output = DyadicFloat256;
//     #[inline]
//     fn sub(self, rhs: DyadicFloat256) -> Self::Output {
//         self.quick_sub(&rhs)
//     }
// }
//
// impl Mul<DyadicFloat256> for DyadicFloat256 {
//     type Output = DyadicFloat256;
//     #[inline]
//     fn mul(self, rhs: DyadicFloat256) -> Self::Output {
//         self.quick_mul(&rhs)
//     }
// }
//
// #[cfg(test)]
// mod tests {
//     use super::*;
//
//     #[test]
//     fn test_dyadic_float() {
//         let ones = DyadicFloat256 {
//             sign: DyadicSign::Pos,
//             exponent: -(BITS as i16 - 1),
//             mantissa: u256 {
//                 hi: 0x80000000_00000000_00000000_00000000_u128,
//                 lo: 0x0u128,
//             },
//         };
//         let cvt = ones.fast_as_f64();
//         assert_eq!(cvt, 1.0);
//
//         let minus_0_5 = DyadicFloat256 {
//             sign: DyadicSign::Neg,
//             exponent: -(BITS as i16),
//             mantissa: u256 {
//                 hi: 0x80000000_00000000_00000000_00000000_u128,
//                 lo: 0x0u128,
//             },
//         };
//         let cvt0 = minus_0_5.fast_as_f64();
//         assert_eq!(cvt0, -1.0 / 2.0);
//
//         let minus_1_f4 = DyadicFloat256 {
//             sign: DyadicSign::Neg,
//             exponent: -(BITS as i16 + 4),
//             mantissa: u256 {
//                 hi: 0xaaaaaaaa_aaaaaaaa_aaaaaaaa_aaaaaaab_u128,
//                 lo: 0x0u128,
//             },
//         };
//         let cvt0 = minus_1_f4.fast_as_f64();
//         assert_eq!(cvt0, -1.0 / 24.0);
//
//         let minus_1_f8 = DyadicFloat256 {
//             sign: DyadicSign::Pos,
//             exponent: -(BITS as i16 + 15),
//             mantissa: u256 {
//                 hi: 0xd00d00d0_0d00d00d_00d00d00_d00d00d0_u128,
//                 lo: 0x0u128,
//             },
//         };
//         let cvt0 = minus_1_f8.fast_as_f64();
//         assert_eq!(cvt0, 1.0 / 40320.0);
//     }
//
//     #[test]
//     fn dyadic_float_add() {
//         let ones = DyadicFloat256 {
//             sign: DyadicSign::Pos,
//             exponent: -(BITS as i16 - 1),
//             mantissa: u256 {
//                 hi: 0x80000000_00000000_00000000_00000000_u128,
//                 lo: 0x0u128,
//             },
//         };
//
//         let cvt = ones.fast_as_f64();
//         assert_eq!(cvt, 1.0);
//
//         let minus_0_5 = DyadicFloat256 {
//             sign: DyadicSign::Neg,
//             exponent: -(BITS as i16),
//             mantissa: u256 {
//                 hi: 0x80000000_00000000_00000000_00000000_u128,
//                 lo: 0x0u128,
//             },
//         };
//         let cvt0 = ones.quick_add(&minus_0_5).fast_as_f64();
//         assert_eq!(cvt0, 0.5);
//     }
//
//     #[test]
//     fn dyadic_float_mul() {
//         let ones = DyadicFloat256 {
//             sign: DyadicSign::Pos,
//             exponent: -(BITS as i16 - 1),
//             mantissa: u256 {
//                 hi: 0x80000000_00000000_00000000_00000000_u128,
//                 lo: 0x0u128,
//             },
//         };
//
//         let cvt = ones.fast_as_f64();
//         assert_eq!(cvt, 1.0);
//
//         let minus_0_5 = DyadicFloat256 {
//             sign: DyadicSign::Neg,
//             exponent: -(BITS as i16),
//             mantissa: u256 {
//                 hi: 0x80000000_00000000_00000000_00000000_u128,
//                 lo: 0x0u128,
//             },
//         };
//         let product = ones.quick_mul(&minus_0_5);
//         let cvt0 = product.fast_as_f64();
//         assert_eq!(cvt0, -0.5);
//
//         let twos = DyadicFloat256 {
//             sign: DyadicSign::Pos,
//             exponent: -(BITS as i16 - 2),
//             mantissa: u256 {
//                 hi: 0x80000000_00000000_00000000_00000000_u128,
//                 lo: 0x0u128,
//             },
//         };
//
//         let cvt = twos.fast_as_f64();
//         assert_eq!(cvt, 2.0);
//     }
//
//     #[test]
//     fn dyadic_round_trip() {
//         let z00 = 0.0;
//         let zvt00 = DyadicFloat256::new_from_f64(z00);
//         let b00 = zvt00.fast_as_f64();
//         assert_eq!(b00, z00);
//
//         let zvt000 = DyadicFloat256 {
//             sign: DyadicSign::Pos,
//             exponent: 0,
//             mantissa: u256::ZERO,
//         };
//         let b000 = zvt000.fast_as_f64();
//         assert_eq!(b000, z00);
//
//         let z0 = 1.0;
//         let zvt0 = DyadicFloat256::new_from_f64(z0);
//         let b0 = zvt0.fast_as_f64();
//         assert_eq!(b0, z0);
//
//         let z1 = 0.5;
//         let zvt1 = DyadicFloat256::new_from_f64(z1);
//         let b1 = zvt1.fast_as_f64();
//         assert_eq!(b1, z1);
//
//         let z2 = -0.5;
//         let zvt2 = DyadicFloat256::new_from_f64(z2);
//         let b2 = zvt2.fast_as_f64();
//         assert_eq!(b2, z2);
//
//         let z3 = -532322.54324324232;
//         let zvt3 = DyadicFloat256::new_from_f64(z3);
//         let b3 = zvt3.fast_as_f64();
//         assert_eq!(b3, z3);
//     }
//
//     #[test]
//     fn dyadic_float_reciprocal() {
//         let ones = DyadicFloat256 {
//             sign: DyadicSign::Pos,
//             exponent: -(BITS as i16 - 1),
//             mantissa: u256 {
//                 hi: 0x80000000_00000000_00000000_00000000_u128,
//                 lo: 0x0u128,
//             },
//         }
//         .reciprocal();
//
//         let cvt = ones.fast_as_f64();
//         assert_eq!(cvt, 1.0);
//
//         let minus_0_5 = DyadicFloat256::new_from_f64(4.).reciprocal();
//         let cvt0 = minus_0_5.fast_as_f64();
//         assert_eq!(cvt0, 0.25);
//     }
//
//     #[test]
//     fn dyadic_float_from_div() {
//         let from_div = DyadicFloat256::from_div_f64(1.0, 4.0);
//         let cvt = from_div.fast_as_f64();
//         assert_eq!(cvt, 0.25);
//     }
//
//     #[test]
//     fn dyadic_float_accurate_reciprocal() {
//         let from_div = DyadicFloat256::accurate_reciprocal(4.0);
//         let cvt = from_div.fast_as_f64();
//         assert_eq!(cvt, 0.25);
//     }
//     /*
//         #[test]
//         fn dyadic_float_mul_int() {
//             let from_div = DyadicFloat256::new_from_f64(4.0);
//             let m1 = from_div.mul_int64(-2);
//             assert_eq!(m1.fast_as_f64(), -8.0);
//
//             let from_div = DyadicFloat256::new_from_f64(-4.0);
//             let m1 = from_div.mul_int64(-2);
//             assert_eq!(m1.fast_as_f64(), 8.0);
//
//             let from_div = DyadicFloat256::new_from_f64(2.5);
//             let m1 = from_div.mul_int64(2);
//             assert_eq!(m1.fast_as_f64(), 5.0);
//         }
//     */
//     #[test]
//     fn dyadic_float_round() {
//         let from_div = DyadicFloat256::new_from_f64(2.5);
//         let m1 = from_div.round_to_nearest_f64();
//         assert_eq!(m1, 3.0);
//
//         let from_div = DyadicFloat256::new_from_f64(0.5);
//         let m1 = from_div.round_to_nearest_f64();
//         assert_eq!(m1, 1.0);
//
//         let from_div = DyadicFloat256::new_from_f64(-0.5);
//         let m1 = from_div.round_to_nearest_f64();
//         assert_eq!(m1, -1.0);
//
//         let from_div = DyadicFloat256::new_from_f64(-0.351);
//         let m1 = from_div.round_to_nearest_f64();
//         assert_eq!(m1, (-0.351f64).round());
//
//         let from_div = DyadicFloat256::new_from_f64(0.351);
//         let m1 = from_div.round_to_nearest_f64();
//         assert_eq!(m1, 0.351f64.round());
//
//         let z00 = 25.6;
//         let zvt00 = DyadicFloat256::new_from_f64(z00);
//         let b00 = zvt00.round_to_nearest_f64();
//         assert_eq!(b00, 26.);
//     }
//     /*
//     #[test]
//     fn dyadic_int_trunc() {
//         let from_div = DyadicFloat256::new_from_f64(-2.5);
//         let m1 = from_div.trunc_to_i64();
//         assert_eq!(m1, -2);
//
//         let from_div = DyadicFloat256::new_from_f64(2.5);
//         let m1 = from_div.trunc_to_i64();
//         assert_eq!(m1, 2);
//
//         let from_div = DyadicFloat256::new_from_f64(0.5);
//         let m1 = from_div.trunc_to_i64();
//         assert_eq!(m1, 0);
//
//         let from_div = DyadicFloat256::new_from_f64(-0.5);
//         let m1 = from_div.trunc_to_i64();
//         assert_eq!(m1, 0);
//
//         let from_div = DyadicFloat256::new_from_f64(-0.351);
//         let m1 = from_div.trunc_to_i64();
//         assert_eq!(m1, 0);
//
//         let from_div = DyadicFloat256::new_from_f64(0.351);
//         let m1 = from_div.trunc_to_i64();
//         assert_eq!(m1, 0);
//     }*/
// }
