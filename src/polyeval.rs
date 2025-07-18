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
use crate::double_double::DoubleDouble;
use crate::dyadic_float::DyadicFloat128;
use std::ops::Mul;

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

impl PolyevalMla for DoubleDouble {
    #[inline(always)]
    fn polyeval_mla(a: Self, b: Self, c: Self) -> Self {
        DoubleDouble::add(c, DoubleDouble::mult(a, b))
    }
}

impl PolyevalMla for DyadicFloat128 {
    #[inline(always)]
    fn polyeval_mla(a: Self, b: Self, c: Self) -> Self {
        c.quick_add(&a.quick_mul(&b))
    }
}

// impl PolyevalMla for DyadicFloat256 {
//     #[inline(always)]
//     fn polyeval_mla(a: Self, b: Self, c: Self) -> Self {
//         c.quick_add(&a.quick_mul(&b))
//     }
// }

#[inline(always)]
#[allow(clippy::too_many_arguments)]
pub(crate) fn f_polyeval6<T: PolyevalMla + Copy + Mul<T, Output = T>>(
    x: T,
    a0: T,
    a1: T,
    a2: T,
    a3: T,
    a4: T,
    a5: T,
) -> T {
    let x2 = x * x;

    let u0 = T::polyeval_mla(x, a5, a4);
    let u1 = T::polyeval_mla(x, a3, a2);
    let u2 = T::polyeval_mla(x, a1, a0);

    let v0 = T::polyeval_mla(x2, u0, u1);

    T::polyeval_mla(x2, v0, u2)
}

// #[inline(always)]
// #[allow(clippy::too_many_arguments)]
// pub(crate) fn f_polyeval5<T: PolyevalMla + Copy + Mul<T, Output = T>>(
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
pub(crate) fn f_polyeval10<T: PolyevalMla + Copy + Mul<T, Output = T>>(
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
    let x2 = x * x;
    let x4 = x2 * x2;
    let x8 = x4 * x4;

    let p0 = T::polyeval_mla(x, a1, a0);
    let p1 = T::polyeval_mla(x, a3, a2);
    let p2 = T::polyeval_mla(x, a5, a4);
    let p3 = T::polyeval_mla(x, a7, a6);
    let p4 = T::polyeval_mla(x, a9, a8);

    let q0 = T::polyeval_mla(x2, p1, p0);
    let q1 = T::polyeval_mla(x2, p3, p2);

    let r0 = T::polyeval_mla(x4, q1, q0);
    T::polyeval_mla(x8, p4, r0)
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
pub(crate) fn f_polyeval13<T: PolyevalMla + Copy + Mul<T, Output = T>>(
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
    let x2 = x * x;
    let x4 = x2 * x2;
    let x8 = x4 * x4;

    let t0 = T::polyeval_mla(x, a3, a2);
    let t1 = T::polyeval_mla(x, a1, a0);
    let t2 = T::polyeval_mla(x, a7, a6);
    let t3 = T::polyeval_mla(x, a5, a4);
    let t4 = T::polyeval_mla(x, a11, a10);
    let t5 = T::polyeval_mla(x, a9, a8);

    let q0 = T::polyeval_mla(x2, t0, t1);
    let q1 = T::polyeval_mla(x2, t2, t3);

    let q2 = T::polyeval_mla(x2, t4, t5);

    let q3 = a12;

    let r0 = T::polyeval_mla(x4, q1, q0);
    let r1 = T::polyeval_mla(x4, q3, q2);

    T::polyeval_mla(x8, r1, r0)
}

#[inline(always)]
#[allow(clippy::too_many_arguments)]
pub(crate) fn f_polyeval12<T: PolyevalMla + Copy>(
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
) -> T {
    let t0 = T::polyeval_mla(x, a11, a10);
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
//     let k0 = T::polyeval_mla(x, a10, a9);
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

// #[inline(always)]
// #[allow(clippy::too_many_arguments)]
// pub(crate) fn f_polyeval14<T: PolyevalMla + Copy + Mul<T, Output = T>>(
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
//     a12: T,
//     a13: T,
// ) -> T {
//     let x2 = x * x;
//     let x4 = x2 * x2;
//     let x8 = x4 * x4;
//
//     let g0 = T::polyeval_mla(x, a1, a0);
//     let g1 = T::polyeval_mla(x, a3, a2);
//     let g2 = T::polyeval_mla(x, a5, a4);
//     let g3 = T::polyeval_mla(x, a7, a6);
//     let g4 = T::polyeval_mla(x, a9, a8);
//     let g5 = T::polyeval_mla(x, a11, a10);
//     let g6 = T::polyeval_mla(x, a13, a12);
//
//     let h0 = T::polyeval_mla(x2, g1, g0);
//     let h1 = T::polyeval_mla(x2, g3, g2);
//     let h2 = T::polyeval_mla(x2, g5, g4);
//
//     let q0 = T::polyeval_mla(x4, h1, h0);
//     let q1 = T::polyeval_mla(x4, g6, h2);
//
//     let result = T::polyeval_mla(x8, q1, q0);
//     result
// }

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
pub(crate) fn f_polyeval16<T: PolyevalMla + Copy + Mul<T, Output = T>>(
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
    let x2 = x * x;
    let x4 = x2 * x2;
    let x8 = x4 * x4;

    let q0 = T::polyeval_mla(x, a1, a0);
    let q1 = T::polyeval_mla(x, a3, a2);
    let q2 = T::polyeval_mla(x, a5, a4);
    let q3 = T::polyeval_mla(x, a7, a6);
    let q4 = T::polyeval_mla(x, a9, a8);
    let q5 = T::polyeval_mla(x, a11, a10);
    let q6 = T::polyeval_mla(x, a13, a12);
    let q7 = T::polyeval_mla(x, a15, a14);

    let r0 = T::polyeval_mla(x2, q1, q0);
    let r1 = T::polyeval_mla(x2, q3, q2);
    let r2 = T::polyeval_mla(x2, q5, q4);
    let r3 = T::polyeval_mla(x2, q7, q6);

    let s0 = T::polyeval_mla(x4, r1, r0);
    let s1 = T::polyeval_mla(x4, r3, r2);

    T::polyeval_mla(x8, s1, s0)
}

// #[inline(always)]
// #[allow(clippy::too_many_arguments)]
// pub(crate) fn f_polyeval15<T: PolyevalMla + Copy + Mul<T, Output = T>>(
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
//     a12: T,
//     a13: T,
//     a14: T,
// ) -> T {
//     let x2 = x * x;
//     let x4 = x2 * x2;
//     let x8 = x4 * x4;
//
//     let e0 = T::polyeval_mla(x, a1, a0);
//     let e1 = T::polyeval_mla(x, a3, a2);
//     let e2 = T::polyeval_mla(x, a5, a4);
//     let e3 = T::polyeval_mla(x, a7, a6);
//     let e4 = T::polyeval_mla(x, a9, a8);
//     let e5 = T::polyeval_mla(x, a11, a10);
//     let e6 = T::polyeval_mla(x, a13, a12);
//
//     // Level 2
//     let f0 = T::polyeval_mla(x2, e1, e0);
//     let f1 = T::polyeval_mla(x2, e3, e2);
//     let f2 = T::polyeval_mla(x2, e5, e4);
//     let f3 = T::polyeval_mla(x2, a14, e6);
//
//     // Level 3
//     let g0 = T::polyeval_mla(x4, f1, f0);
//     let g1 = T::polyeval_mla(x4, f3, f2);
//
//     // Final
//     let result = T::polyeval_mla(x8, g1, g0);
//     result
// }

#[inline(always)]
#[allow(clippy::too_many_arguments)]
pub(crate) fn f_polyeval18<T: PolyevalMla + Copy + Mul<T, Output = T>>(
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
    let x2 = x * x;
    let x4 = x2 * x2;
    let x8 = x4 * x4;
    let x16 = x8 * x8;

    let q0 = T::polyeval_mla(x, a1, a0);
    let q1 = T::polyeval_mla(x, a3, a2);
    let q2 = T::polyeval_mla(x, a5, a4);
    let q3 = T::polyeval_mla(x, a7, a6);
    let q4 = T::polyeval_mla(x, a9, a8);
    let q5 = T::polyeval_mla(x, a11, a10);
    let q6 = T::polyeval_mla(x, a13, a12);
    let q7 = T::polyeval_mla(x, a15, a14);
    let q8 = T::polyeval_mla(x, a17, a16);

    let r0 = T::polyeval_mla(x2, q1, q0);
    let r1 = T::polyeval_mla(x2, q3, q2);
    let r2 = T::polyeval_mla(x2, q5, q4);
    let r3 = T::polyeval_mla(x2, q7, q6);

    let s0 = T::polyeval_mla(x4, r1, r0);
    let s1 = T::polyeval_mla(x4, r3, r2);

    let t0 = T::polyeval_mla(x8, s1, s0);

    T::polyeval_mla(x16, q8, t0)
}

// #[inline(always)]
// #[allow(clippy::too_many_arguments)]
// pub(crate) fn f_polyeval17<T: PolyevalMla + Copy>(
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
//     a12: T,
//     a13: T,
//     a14: T,
//     a15: T,
//     a16: T,
// ) -> T {
//     let z01 = T::polyeval_mla(x, a16, a15);
//     let t1 = T::polyeval_mla(x, z01, a14);
//     let t2 = T::polyeval_mla(x, t1, a13);
//     let t3 = T::polyeval_mla(x, t2, a12);
//     let t4 = T::polyeval_mla(x, t3, a11);
//     let t5 = T::polyeval_mla(x, t4, a10);
//     let t6 = T::polyeval_mla(x, t5, a9);
//     let t7 = T::polyeval_mla(x, t6, a8);
//     let t8 = T::polyeval_mla(x, t7, a7);
//     let t9 = T::polyeval_mla(x, t8, a6);
//     let t10 = T::polyeval_mla(x, t9, a5);
//     let t11 = T::polyeval_mla(x, t10, a4);
//     let t12 = T::polyeval_mla(x, t11, a3);
//     let t13 = T::polyeval_mla(x, t12, a2);
//     let t14 = T::polyeval_mla(x, t13, a1);
//     T::polyeval_mla(x, t14, a0)
// }

// #[inline(always)]
// #[allow(clippy::too_many_arguments)]
// pub(crate) fn f_polyeval19<T: PolyevalMla + Copy>(
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
//     a12: T,
//     a13: T,
//     a14: T,
//     a15: T,
//     a16: T,
//     a17: T,
//     a18: T,
// ) -> T {
//     let z000 = T::polyeval_mla(x, a18, a17);
//     let z00 = T::polyeval_mla(x, z000, a16);
//     let z01 = T::polyeval_mla(x, z00, a15);
//     let t1 = T::polyeval_mla(x, z01, a14);
//     let t2 = T::polyeval_mla(x, t1, a13);
//     let t3 = T::polyeval_mla(x, t2, a12);
//     let t4 = T::polyeval_mla(x, t3, a11);
//     let t5 = T::polyeval_mla(x, t4, a10);
//     let t6 = T::polyeval_mla(x, t5, a9);
//     let t7 = T::polyeval_mla(x, t6, a8);
//     let t8 = T::polyeval_mla(x, t7, a7);
//     let t9 = T::polyeval_mla(x, t8, a6);
//     let t10 = T::polyeval_mla(x, t9, a5);
//     let t11 = T::polyeval_mla(x, t10, a4);
//     let t12 = T::polyeval_mla(x, t11, a3);
//     let t13 = T::polyeval_mla(x, t12, a2);
//     let t14 = T::polyeval_mla(x, t13, a1);
//     T::polyeval_mla(x, t14, a0)
// }
