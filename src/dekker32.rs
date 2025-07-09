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
#![allow(unused, dead_code)]

#[derive(Copy, Clone, Default, Debug)]
pub(crate) struct Dekker32 {
    pub(crate) lo: f32,
    pub(crate) hi: f32,
}

impl Dekker32 {
    #[inline]
    pub(crate) const fn new(lo: f32, hi: f32) -> Self {
        Dekker32 { lo, hi }
    }

    // Non FMA helper
    #[allow(dead_code)]
    #[inline]
    pub(crate) const fn split(a: f32) -> Dekker32 {
        // CN = 2^N.
        const CN: f32 = (1 << 12) as f32;
        const C: f32 = CN + 1.0;
        let t1 = C * a;
        let t2 = a - t1;
        let r_hi = t1 + t2;
        let r_lo = a - r_hi;
        Dekker32::new(r_lo, r_hi)
    }

    // Non FMA helper
    #[allow(dead_code)]
    #[inline]
    fn from_exact_mult_impl_non_fma(asz: Dekker32, a: f32, b: f32) -> Self {
        let bs = Dekker32::split(b);

        let r_hi = a * b;
        let t1 = asz.hi * bs.hi - r_hi;
        let t2 = asz.hi * bs.lo + t1;
        let t3 = asz.lo * bs.hi + t2;
        let r_lo = asz.lo * bs.lo + t3;
        Dekker32::new(r_lo, r_hi)
    }

    #[inline]
    pub(crate) const fn from_full_exact_add(a: f32, b: f32) -> Dekker32 {
        let r_hi = a + b;
        let t1 = r_hi - a;
        let t2 = r_hi - t1;
        let t3 = b - t1;
        let t4 = a - t2;
        let r_lo = t3 + t4;
        Dekker32::new(r_lo, r_hi)
    }

    // // valid only for |a| > b
    // #[inline]
    // pub(crate) const fn from_exact_add(a: f32, b: f32) -> Dekker32 {
    //     let r_hi = a + b;
    //     let t = r_hi - a;
    //     let r_lo = b - t;
    //     Dekker32::new(r_lo, r_hi)
    // }

    #[inline]
    pub(crate) fn from_exact_mult(a: f32, b: f32) -> Self {
        #[cfg(any(
            all(
                any(target_arch = "x86", target_arch = "x86_64"),
                target_feature = "fma"
            ),
            all(target_arch = "aarch64", target_feature = "neon")
        ))]
        {
            use crate::common::f_fmlaf;
            let r_hi = a * b;
            let r_lo = f_fmlaf(a, b, -r_hi);
            Dekker32::new(r_lo, r_hi)
        }
        #[cfg(not(any(
            all(
                any(target_arch = "x86", target_arch = "x86_64"),
                target_feature = "fma"
            ),
            all(target_arch = "aarch64", target_feature = "neon")
        )))]
        {
            let splat = Dekker32::split(a);
            Dekker32::from_exact_mult_impl_non_fma(splat, a, b)
        }
    }

    #[allow(unused)]
    #[inline]
    pub(crate) fn dd_f32_mul_add(a: f32, b: f32, c: f32) -> f32 {
        let ddx2 = Dekker32::from_exact_mult(a, b);
        let zv = Dekker32::full_add_f32(ddx2, c);
        zv.to_f32()
    }

    #[inline]
    pub(crate) fn full_add_f32(a: Dekker32, b: f32) -> Self {
        let t = Dekker32::from_full_exact_add(a.hi, b);
        let l = a.lo + t.lo;
        Self { lo: l, hi: t.hi }
    }

    #[inline]
    pub(crate) fn to_f32(self) -> f32 {
        self.lo + self.hi
    }
}
