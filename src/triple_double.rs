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

/*use crate::dekker::Dekker;
use std::ops::{Add, Sub};

#[derive(Clone, Copy, Debug)]
pub(crate) struct TripleDouble {
    pub(crate) hi: f64,
    pub(crate) mid: f64,
    pub(crate) lo: f64,
}

impl TripleDouble {
    #[inline(always)]
    fn two_sum(a: f64, b: f64) -> (f64, f64) {
        let s = a + b;
        let bv = s - a;
        let err = (a - (s - bv)) + (b - bv);
        (s, err)
    }

    #[inline(always)]
    fn quick_two_sum(a: f64, b: f64) -> (f64, f64) {
        let s = a + b;
        let err = b - (s - a);
        (s, err)
    }

    #[inline(always)]
    fn two_prod(a: f64, b: f64) -> (f64, f64) {
        let p = Dekker::from_exact_mult(a, b);
        (p.hi, p.lo)
    }

    pub(crate) fn add(self, other: TripleDouble) -> TripleDouble {
        let (s0, e0) = Self::two_sum(self.hi, other.hi);
        let (s1, e1) = Self::two_sum(self.mid, other.mid);
        let (s2, e2) = Self::two_sum(self.lo, other.lo);

        let (t0, e3) = Self::two_sum(s1, e0);
        let (t1, e4) = Self::two_sum(s2, e1 + e3);
        let t2 = e2 + e4;

        let (hi, mid) = Self::quick_two_sum(s0, t0);
        let (mid, lo) = Self::quick_two_sum(mid, t1 + t2);

        TripleDouble { hi, mid, lo }
    }

    /// Triple-double subtraction: self - other (approximate)
    pub(crate) fn sub(self, other: TripleDouble) -> TripleDouble {
        self.add(TripleDouble {
            hi: -other.hi,
            mid: -other.mid,
            lo: -other.lo,
        })
    }

    /// Multiply triple-double by f64 scalar
   pub(crate) fn mul_f64(self, b: f64) -> TripleDouble {
        let (p0, e0) = Self::two_prod(self.hi, b);
        let (p1, e1) = Self::two_prod(self.mid, b);
        let (p2, e2) = Self::two_prod(self.lo, b);

        let (s0, c0) = Self::two_sum(p1, e0);
        let (s1, c1) = Self::two_sum(p2, e1 + c0);
        let lo = e2 + c1;

        let (hi, mid) = Self::quick_two_sum(p0, s0);
        let (mid, lo) = Self::quick_two_sum(mid, s1 + lo);

        TripleDouble { hi, mid, lo }
    }

    fn mul_f64_scalar(a: f64, b: f64) -> TripleDouble {
        let (p, e) = Self::two_prod(a, b);
        TripleDouble {
            hi: p,
            mid: e,
            lo: 0.0,
        }
    }

    pub(crate) fn mul(a: TripleDouble, b: TripleDouble) -> TripleDouble {
        let (p1, e1) = Self::two_prod(a.hi, b.hi);
        let (p2, e2) = Self::two_prod(a.hi, b.mid);
        let (p3, e3) = Self::two_prod(a.mid, b.hi);

        let (s1, t1) = Self::two_sum(p1, p2);
        let (s2, t2) = Self::two_sum(s1, p3);

        let s3 = e1 + e2 + e3 + t1 + t2;

        Self::renorm(s2, s3)
    }

    fn renorm(hi: f64, rest: f64) -> TripleDouble {
        let (s1, e1) = Self::two_sum(hi, rest);
        let (s2, e2) = Self::two_sum(s1, e1);
        TripleDouble { hi: s2, mid: e2, lo: 0.0 }
    }

    pub(crate) fn recip_f64(x: f64) -> Self {
        // Initial approximation in f64
        let r0 = 1.0 / x;

        // Promote to triple-double
        let mut r = TripleDouble { hi: r0, mid: 0.0, lo: 0.0 };

        // === Newton-Raphson iterations ===
        // r ← r + r * (1 - x * r)
        for _ in 0..2 {
            // b * r (as TD)
            let xr = Self::mul_f64(r, x);

            // 1 - x * r
            let e = TripleDouble {
                hi: 1.0 - xr.hi,
                mid: -xr.mid,
                lo: -xr.lo,
            };

            // r * e
            let correction = Self::mul(r, e);

            // r = r + correction
            r = r.add(correction);
        }

        r
    }

    /// Division triple-double / f64
    pub(crate) fn div_f64(self, b: f64) -> TripleDouble {
        // Step 1: q0 = a / b
        let q0 = self.hi / b;

        // Step 2: compute r0 = a - q0 * b
        let qb0 = Self::mul_f64_scalar(q0, b);
        let r0 = self.sub(qb0);

        // Step 3: q1 = r0.hi / b
        let q1 = r0.hi / b;

        // Step 4: r1 = r0 - q1 * b
        let qb1 = Self::mul_f64_scalar(q1, b);
        let r1 = r0.sub(qb1);

        // Step 5: q2 = r1.hi / b
        let q2 = r1.hi / b;

        TripleDouble {
            hi: q0,
            mid: q1,
            lo: q2,
        }
    }

    pub(crate) fn div_f64_newton_rapshon(self, b: f64) -> TripleDouble {
        let mut q = TripleDouble::div_f64(self, b);
        // Compute e = 1 - b * q
        let bq = Self::mul_f64(q, b);
        let e = self.sub(bq);

        // Refine: q = q + q * e
        let qe = Self::div_f64(e, b);
        q = q.add(qe);

        q
    }

    #[inline]
    pub(crate) fn to_f64(self) -> f64 {
        (self.hi + self.mid) + self.lo
    }
}
*/
