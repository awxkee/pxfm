// /*
//  * // Copyright (c) Radzivon Bartoshyk 7/2025. All rights reserved.
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
// use crate::dekker::Dekker;
//
// #[derive(Clone, Copy, Debug)]
// pub(crate) struct TripleDoublePacked {
//     pub(crate) hi: u64,
//     pub(crate) mid: u64,
//     pub(crate) lo: u64,
// }
//
// #[derive(Clone, Copy, Debug)]
// pub(crate) struct TripleDouble {
//     pub(crate) hi: f64,
//     pub(crate) mid: f64,
//     pub(crate) lo: f64,
// }
//
// impl TripleDouble {
//     #[inline]
//     pub(crate) fn from_bit_pair(p0: (u64, u64, u64)) -> TripleDouble {
//         TripleDouble {
//             hi: f64::from_bits(p0.2),
//             mid: f64::from_bits(p0.1),
//             lo: f64::from_bits(p0.0),
//         }
//     }
// }
//
// impl TripleDoublePacked {
//     #[inline]
//     pub(crate) fn new(hi: u64, mid: u64, lo: u64) -> Self {
//         Self { hi, mid, lo }
//     }
//
//     #[inline]
//     pub(crate) const fn unpack(self) -> TripleDouble {
//         TripleDouble {
//             hi: f64::from_bits(self.hi),
//             mid: f64::from_bits(self.mid),
//             lo: f64::from_bits(self.lo),
//         }
//     }
// }
//
// impl TripleDouble {
//     pub(crate) const fn from_packed(u: TripleDoublePacked) -> TripleDouble {
//         TripleDouble {
//             hi: f64::from_bits(u.hi),
//             lo: f64::from_bits(u.lo),
//             mid: f64::from_bits(u.mid),
//         }
//     }
// }
//
// impl TripleDouble {
//     #[inline]
//     pub(crate) fn from_add_to_f64(x: f64, rhs: TripleDouble) -> TripleDouble {
//         let Dekker { hi: s0, lo: e0 } = Dekker::from_full_exact_add(x, rhs.hi);
//         let Dekker { hi: s1, lo: e1 } = Dekker::from_full_exact_add(rhs.mid, e0);
//         let Dekker { hi: s2, lo: e2 } = Dekker::from_full_exact_add(rhs.lo, e1);
//
//         let Dekker { hi: r1, lo: t1 } = Dekker::from_exact_add(s1, s2 + e2);
//         let Dekker { hi: r0, lo: r1 } = Dekker::from_exact_add(s0, r1);
//         TripleDouble::new(t1, r1, r0)
//     }
//
//     #[inline]
//     pub(crate) fn to_dekker(self) -> Dekker {
//         Dekker::new(self.lo + self.mid, self.hi)
//     }
//
//     #[inline]
//     pub(crate) const fn new(lo: f64, mid: f64, hi: f64) -> Self {
//         Self { hi, mid, lo }
//     }
//
//     #[inline]
//     pub(crate) const fn from_f64(x: f64) -> Self {
//         Self {
//             hi: x,
//             lo: 0.,
//             mid: 0.,
//         }
//     }
//
//     #[inline]
//     pub(crate) fn from_reciprocal(x: f64) -> TripleDouble {
//         let mut y = TripleDouble::from_f64(1. / x);
//         // r1 = 1 - x * y
//
//         let xy = y.mul_f64(x); // ~106 bits
//         let one = TripleDouble::from_f64(1.);
//         let r1 = one.sub(xy); // r1 ≈ residual
//         let correction1 = y.mul(r1); // y * (1 - x * y)
//         y = y.add(correction1); // improved y
//
//         let x0y = y.mul_f64(x); // ~106 bits
//         let r2 = one.sub(x0y); // r1 ≈ residual
//         let correction2 = y.mul(r2); // y * (1 - x * y)
//         y = y.add(correction2); // improved y
//
//         y
//     }
//
//     #[inline]
//     pub(crate) fn add(self, other: TripleDouble) -> TripleDouble {
//         let Dekker { hi: s0, lo: e0 } = Dekker::from_full_exact_add(self.hi, other.hi);
//         let Dekker { hi: s1, lo: e1 } = Dekker::from_full_exact_add(self.mid, other.mid);
//         let Dekker { hi: s2, lo: e2 } = Dekker::from_full_exact_add(self.lo, other.lo);
//
//         let Dekker { hi: t0, lo: e3 } = Dekker::from_full_exact_add(s1, e0);
//         let Dekker { hi: t1, lo: e4 } = Dekker::from_full_exact_add(s2, e1 + e3);
//         let t2 = e2 + e4;
//
//         let Dekker { hi, lo: mid } = Dekker::from_exact_add(s0, t0);
//         let Dekker { hi: mid, lo } = Dekker::from_exact_add(mid, t1 + t2);
//
//         TripleDouble { hi, mid, lo }
//     }
//
//     #[inline]
//     pub(crate) fn sub(self, other: TripleDouble) -> TripleDouble {
//         let Dekker { hi: s0, lo: e0 } = Dekker::from_full_exact_sub(self.hi, other.hi);
//         let Dekker { hi: s1, lo: e1 } = Dekker::from_full_exact_sub(self.mid, other.mid);
//         let Dekker { hi: s2, lo: e2 } = Dekker::from_full_exact_sub(self.lo, other.lo);
//
//         let Dekker { hi: t0, lo: e3 } = Dekker::from_full_exact_add(s1, e0);
//         let Dekker { hi: t1, lo: e4 } = Dekker::from_full_exact_add(s2, e1 + e3);
//         let t2 = e2 + e4;
//
//         let Dekker { hi, lo: mid } = Dekker::from_exact_add(s0, t0);
//         let Dekker { hi: mid, lo } = Dekker::from_exact_add(mid, t1 + t2);
//
//         TripleDouble { hi, mid, lo }
//     }
//
//     #[inline]
//     fn renorm(self) -> TripleDouble {
//         let Dekker { hi: s0, lo: e0 } = Dekker::from_full_exact_add(self.hi, self.mid);
//         let Dekker { hi: s1, lo: e1 } = Dekker::from_full_exact_add(e0, self.lo);
//         let Dekker { hi, lo: mid } = Dekker::from_full_exact_add(s0, s1);
//         let lo = e1;
//
//         TripleDouble { hi, mid, lo }
//     }
//
//     #[inline]
//     pub(crate) fn div_f64(self, b: f64) -> TripleDouble {
//         // Step 1: Approximate quotient
//         let q0 = self.hi / b;
//
//         // Step 2: Compute residual (hi - q0 * b), using FMA
//         let r0 = f64::mul_add(-q0, b, self.hi); // exact residual
//
//         // Step 3: Refine with mid-part
//         let r1 = r0 + self.mid;
//         let q1 = r1 / b;
//
//         // Step 4: Second residual (using remaining error and lo)
//         let r2 = f64::mul_add(-q1, b, r1 + self.lo);
//         let q2 = r2 / b;
//
//         // Step 5: Sum and renormalize
//         let Dekker { hi: s0, lo: e0 } = Dekker::from_full_exact_add(q0, q1);
//         let Dekker { hi: s1, lo: e1 } = Dekker::from_full_exact_add(e0, q2);
//         let Dekker { hi, lo: mid } = Dekker::from_full_exact_add(s0, s1);
//         let lo = e1;
//
//         TripleDouble { hi, mid, lo }
//     }
//
//     #[inline]
//     fn mul_f64(self, b: f64) -> TripleDouble {
//         let Dekker { hi: p0, lo: e0 } = Dekker::from_exact_mult(self.hi, b);
//         let Dekker { hi: p1, lo: e1 } = Dekker::from_exact_mult(self.mid, b);
//         let Dekker { hi: p2, lo: e2 } = Dekker::from_exact_mult(self.lo, b);
//
//         let Dekker { hi: s1, lo: t1 } = Dekker::from_full_exact_add(e0, p1);
//         let Dekker { hi: s2, lo: t2 } = Dekker::from_full_exact_add(t1, e1);
//         let Dekker { hi: s3, lo: t3 } = Dekker::from_full_exact_add(t2, p2);
//         let s4 = t3 + e2;
//
//         let Dekker { hi: r0, lo: r1_ } = Dekker::from_exact_add(p0, s1);
//         let Dekker { hi: r1, lo: r2_ } = Dekker::from_exact_add(r1_, s2);
//         let r2 = r2_ + s3 + s4;
//
//         TripleDouble {
//             hi: r0,
//             mid: r1,
//             lo: r2,
//         }
//     }
//
//     #[inline]
//     fn mul_short(self, rhs: TripleDouble) -> TripleDouble {
//         // Only partial product expansion, enough for one correction step
//         let Dekker { hi: p0, lo: e0 } =
//             Dekker::from_full_exact_add(self.hi * rhs.hi, self.hi * rhs.mid + self.mid * rhs.hi);
//         let Dekker { hi: p1, lo: e1 } = Dekker::from_full_exact_add(
//             e0,
//             self.hi * rhs.lo + self.mid * rhs.mid + self.lo * rhs.hi,
//         );
//         let Dekker { hi, lo: mid } = Dekker::from_full_exact_add(p0, p1);
//         let lo = e1; // error term
//         TripleDouble { hi, mid, lo }
//     }
//
//     #[inline]
//     pub(crate) fn mul(self, other: TripleDouble) -> TripleDouble {
//         // Step 1: compute leading product and error
//         let Dekker { hi: p0, lo: e0 } = Dekker::from_exact_mult(self.hi, other.hi);
//
//         // Step 2: compute cross terms
//         let Dekker { hi: p1a, lo: e1a } = Dekker::from_exact_mult(self.hi, other.mid);
//         let Dekker { hi: p1b, lo: e1b } = Dekker::from_exact_mult(self.mid, other.hi);
//
//         let Dekker { hi: s1, lo: e1 } = Dekker::from_full_exact_add(p1a, p1b);
//         let e1_total = e1 + e1a + e1b;
//
//         // Step 3: low-order terms
//         let Dekker { hi: p2a, lo: e2a } = Dekker::from_exact_mult(self.hi, other.lo);
//         let Dekker { hi: p2b, lo: e2b } = Dekker::from_exact_mult(self.mid, other.mid);
//         let Dekker { hi: p2c, lo: e2c } = Dekker::from_exact_mult(self.lo, other.hi);
//
//         let Dekker { hi: s2_1, lo: e2_1 } = Dekker::from_full_exact_add(p2a, p2b);
//         let Dekker { hi: s2, lo: e2_2 } = Dekker::from_full_exact_add(s2_1, p2c);
//         let e2_total = e2a + e2b + e2c + e2_1 + e2_2;
//
//         // Final accumulation
//         let Dekker { hi: r0, lo: t1 } = Dekker::from_exact_add(p0, e0);
//         let Dekker { hi: r1, lo: t2 } = Dekker::from_exact_add(s1, e1_total + t1);
//         let Dekker { hi: r2, lo: _ } = Dekker::from_exact_add(s2, e2_total + t2);
//
//         TripleDouble {
//             hi: r0,
//             mid: r1,
//             lo: r2,
//         }
//         .renorm()
//     }
//
//     #[inline]
//     pub(crate) fn div(self, rhs: TripleDouble) -> TripleDouble {
//         // Step 1: Initial rough reciprocal: 1 / rhs
//         let inv_hi = 1.0 / rhs.hi;
//
//         // Step 2: Multiply self by reciprocal estimate (TripleDouble × f64)
//         let q0 = self.mul_f64(inv_hi);
//
//         // Step 3: Compute rhs * q0 (TripleDouble × TripleDouble)
//         let prod = rhs.mul_short(q0);
//
//         // Step 4: Compute residual r = self - rhs * q0
//         let r = self.sub(prod);
//
//         // Step 5: Correction term = r.hi / rhs.hi
//         let delta = r.hi / rhs.hi;
//         let q1 = TripleDouble::from_add_to_f64(delta, q0);
//
//         q1.renorm()
//     }
//
//     #[inline]
//     pub(crate) fn recip(self) -> TripleDouble {
//         // Step 1: Initial approximation using f64 reciprocal
//         let approx = 1.0 / self.hi;
//         let mut y = TripleDouble::from_f64(approx);
//
//         // Step 2: Newton–Raphson: y = y * (2 - x * y)
//         // 1st refinement
//         let xy = self.mul(y);
//         let one = TripleDouble::from_f64(1.);
//         let diff = one.sub(xy);
//         let r = y.mul(diff);
//         y = y.add(r);
//
//         // Optional: 2nd refinement for full precision
//         let xy = self.mul(y);
//         let diff = one.sub(xy);
//         let r = y.mul(diff);
//         y = y.add(r);
//
//         y
//     }
//
//     // pub fn div(self, rhs: TripleDouble) -> TripleDouble {
//     //     // Step 1: Initial rough reciprocal: 1 / rhs
//     //     let inv_hi = 1.0 / rhs.hi;
//     //
//     //     // Step 2: Multiply self by reciprocal estimate (TripleDouble × f64)
//     //     let mut q0 = self.mul_f64(inv_hi);
//     //
//     //     // Step 3: Compute rhs * q0 (TripleDouble × TripleDouble)
//     //     let prod = rhs.mul_triple(q0);
//     //
//     //     // Step 4: Compute residual r = self - rhs * q0
//     //     let r = self.sub_triple(prod);
//     //
//     //     // Step 5: Correction term = r.hi / rhs.hi
//     //     let delta = r.hi / rhs.hi;
//     //     let q1 = q0.add_f64(delta);
//     //
//     //     q1.renormalize()
//     // }
// }
//
// /*
// impl TripleDouble {
//     #[inline(always)]
//     fn two_sum(a: f64, b: f64) -> (f64, f64) {
//         let s = a + b;
//         let bv = s - a;
//         let err = (a - (s - bv)) + (b - bv);
//         (s, err)
//     }
//
//     #[inline(always)]
//     fn quick_two_sum(a: f64, b: f64) -> (f64, f64) {
//         let s = a + b;
//         let err = b - (s - a);
//         (s, err)
//     }
//
//     #[inline(always)]
//     fn two_prod(a: f64, b: f64) -> (f64, f64) {
//         let p = Dekker::from_exact_mult(a, b);
//         (p.hi, p.lo)
//     }
//
//     /// Triple-double subtraction: self - other (approximate)
//     pub(crate) fn sub(self, other: TripleDouble) -> TripleDouble {
//         self.add(TripleDouble {
//             hi: -other.hi,
//             mid: -other.mid,
//             lo: -other.lo,
//         })
//     }
//
//     /// Multiply triple-double by f64 scalar
//    pub(crate) fn mul_f64(self, b: f64) -> TripleDouble {
//         let (p0, e0) = Self::two_prod(self.hi, b);
//         let (p1, e1) = Self::two_prod(self.mid, b);
//         let (p2, e2) = Self::two_prod(self.lo, b);
//
//         let (s0, c0) = Self::two_sum(p1, e0);
//         let (s1, c1) = Self::two_sum(p2, e1 + c0);
//         let lo = e2 + c1;
//
//         let (hi, mid) = Self::quick_two_sum(p0, s0);
//         let (mid, lo) = Self::quick_two_sum(mid, s1 + lo);
//
//         TripleDouble { hi, mid, lo }
//     }
//
//     fn mul_f64_scalar(a: f64, b: f64) -> TripleDouble {
//         let (p, e) = Self::two_prod(a, b);
//         TripleDouble {
//             hi: p,
//             mid: e,
//             lo: 0.0,
//         }
//     }
//
//     pub(crate) fn mul(a: TripleDouble, b: TripleDouble) -> TripleDouble {
//         let (p1, e1) = Self::two_prod(a.hi, b.hi);
//         let (p2, e2) = Self::two_prod(a.hi, b.mid);
//         let (p3, e3) = Self::two_prod(a.mid, b.hi);
//
//         let (s1, t1) = Self::two_sum(p1, p2);
//         let (s2, t2) = Self::two_sum(s1, p3);
//
//         let s3 = e1 + e2 + e3 + t1 + t2;
//
//         Self::renorm(s2, s3)
//     }
//
//     fn renorm(hi: f64, rest: f64) -> TripleDouble {
//         let (s1, e1) = Self::two_sum(hi, rest);
//         let (s2, e2) = Self::two_sum(s1, e1);
//         TripleDouble { hi: s2, mid: e2, lo: 0.0 }
//     }
//
//     pub(crate) fn recip_f64(x: f64) -> Self {
//         // Initial approximation in f64
//         let r0 = 1.0 / x;
//
//         // Promote to triple-double
//         let mut r = TripleDouble { hi: r0, mid: 0.0, lo: 0.0 };
//
//         // === Newton-Raphson iterations ===
//         // r ← r + r * (1 - x * r)
//         for _ in 0..2 {
//             // b * r (as TD)
//             let xr = Self::mul_f64(r, x);
//
//             // 1 - x * r
//             let e = TripleDouble {
//                 hi: 1.0 - xr.hi,
//                 mid: -xr.mid,
//                 lo: -xr.lo,
//             };
//
//             // r * e
//             let correction = Self::mul(r, e);
//
//             // r = r + correction
//             r = r.add(correction);
//         }
//
//         r
//     }
//
//     /// Division triple-double / f64
//     pub(crate) fn div_f64(self, b: f64) -> TripleDouble {
//         // Step 1: q0 = a / b
//         let q0 = self.hi / b;
//
//         // Step 2: compute r0 = a - q0 * b
//         let qb0 = Self::mul_f64_scalar(q0, b);
//         let r0 = self.sub(qb0);
//
//         // Step 3: q1 = r0.hi / b
//         let q1 = r0.hi / b;
//
//         // Step 4: r1 = r0 - q1 * b
//         let qb1 = Self::mul_f64_scalar(q1, b);
//         let r1 = r0.sub(qb1);
//
//         // Step 5: q2 = r1.hi / b
//         let q2 = r1.hi / b;
//
//         TripleDouble {
//             hi: q0,
//             mid: q1,
//             lo: q2,
//         }
//     }
//
//     pub(crate) fn div_f64_newton_rapshon(self, b: f64) -> TripleDouble {
//         let mut q = TripleDouble::div_f64(self, b);
//         // Compute e = 1 - b * q
//         let bq = Self::mul_f64(q, b);
//         let e = self.sub(bq);
//
//         // Refine: q = q + q * e
//         let qe = Self::div_f64(e, b);
//         q = q.add(qe);
//
//         q
//     }
//
//     #[inline]
//     pub(crate) fn to_f64(self) -> f64 {
//         (self.hi + self.mid) + self.lo
//     }
// }
// */
