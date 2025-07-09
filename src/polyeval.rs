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
use crate::common::{f_fmla, f_fmlaf};
use crate::dekker::Dekker;
use crate::dyadic_float::DyadicFloat128;

pub(crate) trait PolyevalMla {
    fn polyeval_mla(a: Self, b: Self, c: Self) -> Self;
}

impl PolyevalMla for f64 {
    #[inline(always)]
    fn polyeval_mla(a: Self, b: Self, c: Self) -> Self {
        f_fmla(a, b, c)
    }
}

impl PolyevalMla for f32 {
    #[inline(always)]
    fn polyeval_mla(a: Self, b: Self, c: Self) -> Self {
        f_fmlaf(a, b, c)
    }
}

impl PolyevalMla for Dekker {
    #[inline(always)]
    fn polyeval_mla(a: Self, b: Self, c: Self) -> Self {
        Dekker::add(c, Dekker::mult(a, b))
    }
}

impl PolyevalMla for DyadicFloat128 {
    #[inline(always)]
    fn polyeval_mla(a: Self, b: Self, c: Self) -> Self {
        c.quick_add(&a.quick_mul(&b))
    }
}

#[inline(always)]
#[allow(clippy::too_many_arguments)]
pub(crate) fn f_polyeval6<T: PolyevalMla + Copy>(
    x: T,
    a0: T,
    a1: T,
    a2: T,
    a3: T,
    a4: T,
    a5: T,
) -> T {
    let t0 = T::polyeval_mla(x, a5, a4);
    let t2 = T::polyeval_mla(x, t0, a3);
    let t3 = T::polyeval_mla(x, t2, a2);
    let t4 = T::polyeval_mla(x, t3, a1);
    T::polyeval_mla(x, t4, a0)
}

// #[inline(always)]
// #[allow(clippy::too_many_arguments)]
// pub(crate) fn f_polyeval5<T: PolyevalMla + Copy>(
//     x: T,
//     a0: T,
//     a1: T,
//     a2: T,
//     a3: T,
//     a4: T,
// ) -> T {
//     let t2 = T::polyeval_mla(x, a4, a3);
//     let t3 = T::polyeval_mla(x, t2, a2);
//     let t4 = T::polyeval_mla(x, t3, a1);
//     T::polyeval_mla(x, t4, a0)
// }

#[inline(always)]
#[allow(clippy::too_many_arguments)]
pub(crate) fn f_polyeval9<T: PolyevalMla + Copy>(
    x: T,
    a0: T,
    a1: T,
    a2: T,
    a3: T,
    a4: T,
    a5: T,
    a6: T,
    a7: T,
    a8: T,
) -> T {
    let t0 = T::polyeval_mla(x, a8, a7);
    let t01 = T::polyeval_mla(x, t0, a6);
    let t1 = T::polyeval_mla(x, t01, a5);
    let t2 = T::polyeval_mla(x, t1, a4);
    let t3 = T::polyeval_mla(x, t2, a3);
    let t4 = T::polyeval_mla(x, t3, a2);
    let t5 = T::polyeval_mla(x, t4, a1);
    T::polyeval_mla(x, t5, a0)
}

#[inline(always)]
#[allow(clippy::too_many_arguments)]
pub(crate) fn f_polyeval10<T: PolyevalMla + Copy>(
    x: T,
    a0: T,
    a1: T,
    a2: T,
    a3: T,
    a4: T,
    a5: T,
    a6: T,
    a7: T,
    a8: T,
    a9: T,
) -> T {
    let z0 = T::polyeval_mla(x, a9, a8);
    let t0 = T::polyeval_mla(x, z0, a7);
    let t01 = T::polyeval_mla(x, t0, a6);
    let t1 = T::polyeval_mla(x, t01, a5);
    let t2 = T::polyeval_mla(x, t1, a4);
    let t3 = T::polyeval_mla(x, t2, a3);
    let t4 = T::polyeval_mla(x, t3, a2);
    let t5 = T::polyeval_mla(x, t4, a1);
    T::polyeval_mla(x, t5, a0)
}

// #[inline(always)]
// #[allow(clippy::too_many_arguments)]
// pub(crate) fn f_polyeval11<T: PolyevalMla + Copy>(
//     x: T,
//     a0: T,
//     a1: T,
//     a2: T,
//     a3: T,
//     a4: T,
//     a5: T,
//     a6: T,
//     a7: T,
//     a8: T,
//     a9: T,
//     a10: T,
// ) -> T {
//     let z00 = T::polyeval_mla(x, a10, a9);
//     let z0 = T::polyeval_mla(x, z00, a8);
//     let t0 = T::polyeval_mla(x, z0, a7);
//     let t01 = T::polyeval_mla(x, t0, a6);
//     let t1 = T::polyeval_mla(x, t01, a5);
//     let t2 = T::polyeval_mla(x, t1, a4);
//     let t3 = T::polyeval_mla(x, t2, a3);
//     let t4 = T::polyeval_mla(x, t3, a2);
//     let t5 = T::polyeval_mla(x, t4, a1);
//     T::polyeval_mla(x, t5, a0)
// }

#[inline(always)]
pub(crate) fn f_polyeval3<T: PolyevalMla + Copy>(x: T, a0: T, a1: T, a2: T) -> T {
    T::polyeval_mla(x, T::polyeval_mla(x, a2, a1), a0)
}

#[inline(always)]
#[allow(clippy::too_many_arguments)]
pub(crate) fn f_polyeval4<T: PolyevalMla + Copy>(x: T, a0: T, a1: T, a2: T, a3: T) -> T {
    let t2 = T::polyeval_mla(x, a3, a2);
    let t5 = T::polyeval_mla(x, t2, a1);
    T::polyeval_mla(x, t5, a0)
}

#[inline(always)]
#[allow(clippy::too_many_arguments)]
pub(crate) fn f_polyeval13<T: PolyevalMla + Copy>(
    x: T,
    a0: T,
    a1: T,
    a2: T,
    a3: T,
    a4: T,
    a5: T,
    a6: T,
    a7: T,
    a8: T,
    a9: T,
    a10: T,
    a11: T,
    a12: T,
) -> T {
    let z00 = T::polyeval_mla(x, a12, a11);
    let t0 = T::polyeval_mla(x, z00, a10);
    let k0 = T::polyeval_mla(x, t0, a9);
    let k1 = T::polyeval_mla(x, k0, a8);
    let z0 = T::polyeval_mla(x, k1, a7);
    let t0a = T::polyeval_mla(x, z0, a6);
    let t1 = T::polyeval_mla(x, t0a, a5);
    let t2 = T::polyeval_mla(x, t1, a4);
    let t3 = T::polyeval_mla(x, t2, a3);
    let t4 = T::polyeval_mla(x, t3, a2);
    let t5 = T::polyeval_mla(x, t4, a1);
    T::polyeval_mla(x, t5, a0)
}

// #[inline(always)]
// #[allow(clippy::too_many_arguments)]
// pub(crate) fn f_polyeval12<T: PolyevalMla + Copy>(
//     x: T,
//     a0: T,
//     a1: T,
//     a2: T,
//     a3: T,
//     a4: T,
//     a5: T,
//     a6: T,
//     a7: T,
//     a8: T,
//     a9: T,
//     a10: T,
//     a11: T,
// ) -> T {
//     let t0 = T::polyeval_mla(x, a11, a10);
//     let k0 = T::polyeval_mla(x, t0, a9);
//     let k1 = T::polyeval_mla(x, k0, a8);
//     let z0 = T::polyeval_mla(x, k1, a7);
//     let t0a = T::polyeval_mla(x, z0, a6);
//     let t1 = T::polyeval_mla(x, t0a, a5);
//     let t2 = T::polyeval_mla(x, t1, a4);
//     let t3 = T::polyeval_mla(x, t2, a3);
//     let t4 = T::polyeval_mla(x, t3, a2);
//     let t5 = T::polyeval_mla(x, t4, a1);
//     T::polyeval_mla(x, t5, a0)
// }

#[inline(always)]
#[allow(clippy::too_many_arguments)]
pub(crate) fn f_polyeval7<T: PolyevalMla + Copy>(
    x: T,
    a0: T,
    a1: T,
    a2: T,
    a3: T,
    a4: T,
    a5: T,
    a6: T,
) -> T {
    let t1 = T::polyeval_mla(x, a6, a5);
    let t2 = T::polyeval_mla(x, t1, a4);
    let t3 = T::polyeval_mla(x, t2, a3);
    let t4 = T::polyeval_mla(x, t3, a2);
    let t5 = T::polyeval_mla(x, t4, a1);
    T::polyeval_mla(x, t5, a0)
}

#[inline(always)]
#[allow(clippy::too_many_arguments)]
pub(crate) fn f_polyeval8<T: PolyevalMla + Copy>(
    x: T,
    a0: T,
    a1: T,
    a2: T,
    a3: T,
    a4: T,
    a5: T,
    a6: T,
    a7: T,
) -> T {
    let z0 = T::polyeval_mla(x, a7, a6);
    let t1 = T::polyeval_mla(x, z0, a5);
    let t2 = T::polyeval_mla(x, t1, a4);
    let t3 = T::polyeval_mla(x, t2, a3);
    let t4 = T::polyeval_mla(x, t3, a2);
    let t5 = T::polyeval_mla(x, t4, a1);
    T::polyeval_mla(x, t5, a0)
}

#[inline(always)]
#[allow(clippy::too_many_arguments)]
pub(crate) fn f_polyeval16<T: PolyevalMla + Copy>(
    x: T,
    a0: T,
    a1: T,
    a2: T,
    a3: T,
    a4: T,
    a5: T,
    a6: T,
    a7: T,
    a8: T,
    a9: T,
    a10: T,
    a11: T,
    a12: T,
    a13: T,
    a14: T,
    a15: T,
) -> T {
    let t1 = T::polyeval_mla(x, a15, a14);
    let t2 = T::polyeval_mla(x, t1, a13);
    let t3 = T::polyeval_mla(x, t2, a12);
    let t4 = T::polyeval_mla(x, t3, a11);
    let t5 = T::polyeval_mla(x, t4, a10);
    let t6 = T::polyeval_mla(x, t5, a9);
    let t7 = T::polyeval_mla(x, t6, a8);
    let t8 = T::polyeval_mla(x, t7, a7);
    let t9 = T::polyeval_mla(x, t8, a6);
    let t10 = T::polyeval_mla(x, t9, a5);
    let t11 = T::polyeval_mla(x, t10, a4);
    let t12 = T::polyeval_mla(x, t11, a3);
    let t13 = T::polyeval_mla(x, t12, a2);
    let t14 = T::polyeval_mla(x, t13, a1);
    T::polyeval_mla(x, t14, a0)
}

#[inline(always)]
#[allow(clippy::too_many_arguments)]
pub(crate) fn f_polyeval18<T: PolyevalMla + Copy>(
    x: T,
    a0: T,
    a1: T,
    a2: T,
    a3: T,
    a4: T,
    a5: T,
    a6: T,
    a7: T,
    a8: T,
    a9: T,
    a10: T,
    a11: T,
    a12: T,
    a13: T,
    a14: T,
    a15: T,
    a16: T,
    a17: T,
) -> T {
    let z00 = T::polyeval_mla(x, a17, a16);
    let z01 = T::polyeval_mla(x, z00, a15);
    let t1 = T::polyeval_mla(x, z01, a14);
    let t2 = T::polyeval_mla(x, t1, a13);
    let t3 = T::polyeval_mla(x, t2, a12);
    let t4 = T::polyeval_mla(x, t3, a11);
    let t5 = T::polyeval_mla(x, t4, a10);
    let t6 = T::polyeval_mla(x, t5, a9);
    let t7 = T::polyeval_mla(x, t6, a8);
    let t8 = T::polyeval_mla(x, t7, a7);
    let t9 = T::polyeval_mla(x, t8, a6);
    let t10 = T::polyeval_mla(x, t9, a5);
    let t11 = T::polyeval_mla(x, t10, a4);
    let t12 = T::polyeval_mla(x, t11, a3);
    let t13 = T::polyeval_mla(x, t12, a2);
    let t14 = T::polyeval_mla(x, t13, a1);
    T::polyeval_mla(x, t14, a0)
}
