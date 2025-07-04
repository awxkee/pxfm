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
#[allow(unused_imports)]
use crate::common::*;

#[derive(Copy, Clone, Default, Debug)]
pub(crate) struct Dekker {
    pub(crate) lo: f64,
    pub(crate) hi: f64,
}

impl Dekker {
    #[inline]
    pub(crate) const fn from_bit_pair(pair: (u64, u64)) -> Self {
        Self {
            lo: f64::from_bits(pair.0),
            hi: f64::from_bits(pair.1),
        }
    }

    #[inline]
    pub(crate) const fn new(lo: f64, hi: f64) -> Self {
        Dekker { lo, hi }
    }

    // Non FMA helper
    #[allow(dead_code)]
    #[inline]
    pub(crate) const fn split(a: f64) -> Dekker {
        // CN = 2^N.
        const CN: f64 = (1 << 27) as f64;
        const C: f64 = CN + 1.0;
        let t1 = C * a;
        let t2 = a - t1;
        let r_hi = t1 + t2;
        let r_lo = a - r_hi;
        Dekker::new(r_lo, r_hi)
    }

    // Non FMA helper
    #[allow(dead_code)]
    #[inline]
    fn from_exact_mult_impl_non_fma(asz: Dekker, a: f64, b: f64) -> Self {
        let bs = Dekker::split(b);

        let r_hi = a * b;
        let t1 = asz.hi * bs.hi - r_hi;
        let t2 = asz.hi * bs.lo + t1;
        let t3 = asz.lo * bs.hi + t2;
        let r_lo = asz.lo * bs.lo + t3;
        Dekker::new(r_lo, r_hi)
    }

    #[inline]
    pub(crate) const fn from_exact_add(a: f64, b: f64) -> Dekker {
        let r_hi = a + b;
        let t = r_hi - a;
        let r_lo = b - t;
        Dekker::new(r_lo, r_hi)
    }

    #[inline]
    pub(crate) const fn from_exact_sub(a: f64, b: f64) -> Dekker {
        let r_hi = a - b;
        let t = a - r_hi;
        let r_lo = t - b;
        Dekker::new(r_lo, r_hi)
    }

    #[inline]
    pub(crate) const fn from_full_exact_add(a: f64, b: f64) -> Dekker {
        let r_hi = a + b;
        let t1 = r_hi - a;
        let t2 = r_hi - t1;
        let t3 = b - t1;
        let t4 = a - t2;
        let r_lo = t3 + t4;
        Dekker::new(r_lo, r_hi)
    }

    #[allow(unused)]
    #[inline]
    pub(crate) fn dd_f64_mul_add(a: f64, b: f64, c: f64) -> f64 {
        let ddx2 = Dekker::from_exact_mult(a, b);
        let zv = Dekker::add_f64(ddx2, c);
        zv.to_f64()
    }

    // #[inline]
    // pub(crate) const fn from_full_exact_sub(a: f64, b: f64) -> Self {
    //     let r_hi = a - b;
    //     let t1 = r_hi - a;
    //     let t2 = r_hi - t1;
    //     let t3 = -b - t1;
    //     let t4 = a - t2;
    //     let r_lo = t3 + t4;
    //     Dekker::new(r_lo, r_hi)
    // }

    #[inline]
    pub(crate) fn add(a: Dekker, b: Dekker) -> Dekker {
        let s = a.hi + b.hi;
        let d = s - a.hi;
        let l = ((b.hi - d) + (a.hi + (d - s))) + (a.lo + b.lo);
        Dekker::new(l, s)
    }

    #[inline]
    pub(crate) fn sub(a: Dekker, b: Dekker) -> Dekker {
        let s = a.hi - b.hi;
        let d = s - a.hi;
        let l = ((-b.hi - d) + (a.hi + (d - s))) + (a.lo - b.lo);
        Dekker::new(l, s)
    }

    #[inline]
    pub(crate) fn from_exact_div(a: f64, b: f64) -> Self {
        #[cfg(any(
            all(
                any(target_arch = "x86", target_arch = "x86_64"),
                target_feature = "fma"
            ),
            all(target_arch = "aarch64", target_feature = "neon")
        ))]
        {
            let q_hi = a / b;
            let r = f_fmla(-q_hi, b, a);
            let q_lo = r / b;
            Self::new(q_lo, q_hi)
        }

        #[cfg(not(any(
            all(
                any(target_arch = "x86", target_arch = "x86_64"),
                target_feature = "fma"
            ),
            all(target_arch = "aarch64", target_feature = "neon")
        )))]
        {
            let q_hi = a / b;

            let p = Dekker::from_exact_mult(q_hi, b);
            let r = Dekker::from_exact_sub(a, p.hi);
            let r = r.hi + (r.lo - p.lo);
            let q_lo = r / b;

            Self::new(q_lo, q_hi)
        }
    }
    //
    // #[inline]
    // pub(crate) fn from_sqrt(x: f64) -> Self {
    //     let h = x.sqrt();
    //     /* h = sqrt(x) * (1 + e1) with |e1| < 2^-52
    //        thus h^2 = x * (1 + e2) with |e2| < 2^-50.999 */
    //     let e = -f_fmla (h, h, -x); // exact
    //     /* e = x - h^2 */
    //     let l = e / (h + h);
    //     Dekker::new(l, h)
    // }

    // #[inline]
    // pub(crate) fn div_dd_f64(a: Dekker, b: f64) -> Self {
    //     let q1 = a.hi / b;
    //     let r = f_fmla(-q1, b, a.hi);
    //     let r = r + a.lo;
    //     let q2 = r / b;
    //
    //     Dekker::new(q2, q1)
    // }
    //
    // #[inline]
    // pub(crate) fn neg(self) -> Self {
    //     Self {
    //         lo: -self.lo, hi: -self.hi,
    //     }
    // }

    // #[inline]
    // pub(crate) fn from_f64_div_dd(a: f64, b: Dekker) -> Self {
    //     let q1 = a / b.hi;
    //
    //     let prod = Dekker::from_exact_mult(q1, b.hi);
    //     let prod_lo = f_fmla(q1, b.lo, prod.lo);
    //     let rem = f_fmla(-1.0, prod.hi, a) - prod_lo;
    //
    //     let q2 = rem / b.hi;
    //
    //     Dekker::new(q2, q1)
    // }

    // #[inline]
    // pub(crate) fn mla_f64(a: Dekker, b: f64, c: f64) -> Self {
    //     let q = Dekker::mult_f64(a, b);
    //     Dekker::add_f64(q, c)
    // }
    //
    // #[inline]
    // pub(crate) fn mla_dd_f64(a: Dekker, b: Dekker, c: f64) -> Self {
    //     let q = Dekker::quick_mult(a, b);
    //     Dekker::add_f64(q, c)
    // }

    #[inline]
    pub(crate) fn div(a: Dekker, b: Dekker) -> Dekker {
        let q = 1.0 / b.hi;
        let r_hi = a.hi * q;
        #[cfg(any(
            all(
                any(target_arch = "x86", target_arch = "x86_64"),
                target_feature = "fma"
            ),
            all(target_arch = "aarch64", target_feature = "neon")
        ))]
        {
            let e_hi = f_fmla(b.hi, -r_hi, a.hi);
            let e_lo = f_fmla(b.lo, -r_hi, a.lo);
            let r_lo = q * (e_hi + e_lo);
            Dekker::new(r_lo, r_hi)
        }

        #[cfg(not(any(
            all(
                any(target_arch = "x86", target_arch = "x86_64"),
                target_feature = "fma"
            ),
            all(target_arch = "aarch64", target_feature = "neon")
        )))]
        {
            let b_hi_r_hi = Dekker::from_exact_mult(b.hi, -r_hi);
            let b_lo_r_hi = Dekker::from_exact_mult(b.lo, -r_hi);
            let e_hi = (a.hi + b_hi_r_hi.hi) + b_hi_r_hi.lo;
            let e_lo = (a.lo + b_lo_r_hi.hi) + b_lo_r_hi.lo;
            let r_lo = q * (e_hi + e_lo);
            Dekker::new(r_lo, r_hi)
        }
    }

    #[inline]
    pub(crate) fn from_exact_mult(a: f64, b: f64) -> Self {
        #[cfg(any(
            all(
                any(target_arch = "x86", target_arch = "x86_64"),
                target_feature = "fma"
            ),
            all(target_arch = "aarch64", target_feature = "neon")
        ))]
        {
            let r_hi = a * b;
            let r_lo = f_fmla(a, b, -r_hi);
            Dekker::new(r_lo, r_hi)
        }
        #[cfg(not(any(
            all(
                any(target_arch = "x86", target_arch = "x86_64"),
                target_feature = "fma"
            ),
            all(target_arch = "aarch64", target_feature = "neon")
        )))]
        {
            let splat = Dekker::split(a);
            Dekker::from_exact_mult_impl_non_fma(splat, a, b)
        }
    }

    // #[inline]
    // pub(crate) fn add_f64(&self, other: f64) -> DoubleDouble {
    //     let r = DoubleDouble::from_exact_add(self.hi, other);
    //     Dekker::from_exact_add(r.hi, r.lo + self.lo)
    // }

    #[inline]
    pub(crate) fn quick_mult(a: Dekker, b: Dekker) -> Self {
        #[cfg(any(
            all(
                any(target_arch = "x86", target_arch = "x86_64"),
                target_feature = "fma"
            ),
            all(target_arch = "aarch64", target_feature = "neon")
        ))]
        {
            let mut r = Dekker::from_exact_mult(a.hi, b.hi);
            let t1 = f_fmla(a.hi, b.lo, r.lo);
            let t2 = f_fmla(a.lo, b.hi, t1);
            r.lo = t2;
            r
        }
        #[cfg(not(any(
            all(
                any(target_arch = "x86", target_arch = "x86_64"),
                target_feature = "fma"
            ),
            all(target_arch = "aarch64", target_feature = "neon")
        )))]
        {
            Dekker::mult(a, b)
        }
    }

    #[inline]
    pub(crate) fn mult(a: Dekker, b: Dekker) -> Self {
        #[cfg(any(
            all(
                any(target_arch = "x86", target_arch = "x86_64"),
                target_feature = "fma"
            ),
            all(target_arch = "aarch64", target_feature = "neon")
        ))]
        {
            let ahlh = b.hi * a.lo;
            let alhh = b.lo * a.hi;
            let ahhh = b.hi * a.hi;
            let mut ahhl = f_fmla(b.hi, a.hi, -ahhh);
            ahhl += alhh + ahlh;
            let ch = ahhh + ahhl;
            let l = (ahhh - ch) + ahhl;
            Self { lo: l, hi: ch }
        }
        #[cfg(not(any(
            all(
                any(target_arch = "x86", target_arch = "x86_64"),
                target_feature = "fma"
            ),
            all(target_arch = "aarch64", target_feature = "neon")
        )))]
        {
            let zp = Dekker::from_exact_mult(a.hi, b.hi);

            let p2 = a.hi * b.lo;
            let p3 = a.lo * b.hi;

            let e1 = Dekker::from_exact_add(zp.lo, p2);
            let e2 = Dekker::from_exact_add(e1.hi, p3);

            let hi = zp.hi + e2.hi;
            let lo = (zp.hi - hi) + e2.hi + e1.lo + e2.lo;

            Self { lo, hi }
        }
    }

    #[inline]
    pub(crate) fn mult_f64(a: Dekker, b: f64) -> Self {
        #[cfg(any(
            all(
                any(target_arch = "x86", target_arch = "x86_64"),
                target_feature = "fma"
            ),
            all(target_arch = "aarch64", target_feature = "neon")
        ))]
        {
            let ahlh = b * a.lo;
            let ahhh = b * a.hi;
            let mut ahhl = f_fmla(b, a.hi, -ahhh);
            ahhl += ahlh;
            let ch = ahhh + ahhl;
            let l = (ahhh - ch) + ahhl;
            Dekker::new(l, ch)
        }
        #[cfg(not(any(
            all(
                any(target_arch = "x86", target_arch = "x86_64"),
                target_feature = "fma"
            ),
            all(target_arch = "aarch64", target_feature = "neon")
        )))]
        {
            let z = Dekker::from_exact_mult(a.hi, b);
            let lo_mul = a.lo * b;
            let s = z.lo + lo_mul;
            let r_hi = z.hi + s;
            let r_lo = s - (r_hi - z.hi);

            Dekker { hi: r_hi, lo: r_lo }
        }
    }

    #[inline]
    pub(crate) fn f64_mult(a: f64, b: Dekker) -> Dekker {
        #[cfg(any(
            all(
                any(target_arch = "x86", target_arch = "x86_64"),
                target_feature = "fma"
            ),
            all(target_arch = "aarch64", target_feature = "neon")
        ))]
        {
            let mut p = Dekker::from_exact_mult(a, b.hi);
            p.lo = f_fmla(a, b.lo, p.lo);
            p
        }
        #[cfg(not(any(
            all(
                any(target_arch = "x86", target_arch = "x86_64"),
                target_feature = "fma"
            ),
            all(target_arch = "aarch64", target_feature = "neon")
        )))]
        {
            let z = Dekker::from_exact_mult(a, b.hi);
            let lo_mul = a * b.lo;

            let s = z.lo + lo_mul;

            let r_hi = z.hi + s;
            let r_lo = s - (r_hi - z.hi);

            Dekker { hi: r_hi, lo: r_lo }
        }
    }

    #[inline]
    pub(crate) fn quick_mult_f64(a: Dekker, b: f64) -> Self {
        #[cfg(any(
            all(
                any(target_arch = "x86", target_arch = "x86_64"),
                target_feature = "fma"
            ),
            all(target_arch = "aarch64", target_feature = "neon")
        ))]
        {
            let h = b * a.hi;
            let l = b * a.lo + f_fmla(b, a.hi, -h);
            Self { lo: l, hi: h }
        }
        #[cfg(not(any(
            all(
                any(target_arch = "x86", target_arch = "x86_64"),
                target_feature = "fma"
            ),
            all(target_arch = "aarch64", target_feature = "neon")
        )))]
        {
            let h = a.hi * b;

            let a_split = Dekker::split(a.hi);
            let b_split = Dekker::split(b);

            let err = (((a_split.hi * b_split.hi - h) + a_split.hi * b_split.lo)
                + a_split.lo * b_split.hi)
                + a_split.lo * b_split.lo;

            let l = a.lo * b + err;

            Self { lo: l, hi: h }
        }
    }

    #[inline]
    pub(crate) fn add_f64(a: Dekker, b: f64) -> Self {
        let t = Dekker::from_exact_add(a.hi, b);
        let l = a.lo + t.lo;
        Self { lo: l, hi: t.hi }
    }

    #[inline]
    pub(crate) fn f64_add(b: f64, a: Dekker) -> Self {
        let t = Dekker::from_exact_add(b, a.hi);
        let l = a.lo + t.lo;
        Self { lo: l, hi: t.hi }
    }

    #[inline]
    pub(crate) const fn to_f64(self) -> f64 {
        self.lo + self.hi
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    #[test]
    fn test_f64_mult() {
        let d1 = 1.1231;
        let d2 = Dekker::new(1e-22, 3.2341);
        let p = Dekker::f64_mult(d1, d2);
        assert_eq!(p.hi, 3.6322177100000004);
        assert_eq!(p.lo, -1.971941841373783e-16);
    }

    #[test]
    fn test_mult_64() {
        let d1 = 1.1231;
        let d2 = Dekker::new(1e-22, 3.2341);
        let p = Dekker::mult_f64(d2, d1);
        assert_eq!(p.hi, 3.6322177100000004);
        assert_eq!(p.lo, -1.971941841373783e-16);
    }
}
