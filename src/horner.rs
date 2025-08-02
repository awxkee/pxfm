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
use crate::polyeval::PolyevalMla;
use std::ops::Mul;

#[inline(always)]
#[allow(clippy::too_many_arguments)]
pub(crate) fn f_horner_polyeval11<T: PolyevalMla + Copy + Mul<T, Output = T>>(
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
) -> T {
    let z00 = T::polyeval_mla(x, a10, a9);
    let z0 = T::polyeval_mla(x, z00, a8);
    let t0 = T::polyeval_mla(x, z0, a7);
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
pub(crate) fn f_horner_polyeval12<T: PolyevalMla + Copy + Mul<T, Output = T>>(
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
    let mut acc = a11;
    acc = T::polyeval_mla(x, acc, a10);
    acc = T::polyeval_mla(x, acc, a9);
    acc = T::polyeval_mla(x, acc, a8);
    acc = T::polyeval_mla(x, acc, a7);
    acc = T::polyeval_mla(x, acc, a6);
    acc = T::polyeval_mla(x, acc, a5);
    acc = T::polyeval_mla(x, acc, a4);
    acc = T::polyeval_mla(x, acc, a3);
    acc = T::polyeval_mla(x, acc, a2);
    acc = T::polyeval_mla(x, acc, a1);
    T::polyeval_mla(x, acc, a0)
}

// #[inline(always)]
// #[allow(clippy::too_many_arguments)]
// pub(crate) fn f_horner_polyeval9<T: PolyevalMla + Copy + Mul<T, Output = T>>(
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
// ) -> T {
//     let t0 = T::polyeval_mla(x, a8, a7);
//     let t01 = T::polyeval_mla(x, t0, a6);
//     let t1 = T::polyeval_mla(x, t01, a5);
//     let t2 = T::polyeval_mla(x, t1, a4);
//     let t3 = T::polyeval_mla(x, t2, a3);
//     let t4 = T::polyeval_mla(x, t3, a2);
//     let t5 = T::polyeval_mla(x, t4, a1);
//     T::polyeval_mla(x, t5, a0)
// }

// #[inline(always)]
// #[allow(clippy::too_many_arguments)]
// pub(crate) fn f_horner_polyeval17<T: PolyevalMla + Copy>(
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
//     let mut acc = a16;
//     acc = T::polyeval_mla(x, acc, a15);
//     acc = T::polyeval_mla(x, acc, a14);
//     acc = T::polyeval_mla(x, acc, a13);
//     acc = T::polyeval_mla(x, acc, a12);
//     acc = T::polyeval_mla(x, acc, a11);
//     acc = T::polyeval_mla(x, acc, a10);
//     acc = T::polyeval_mla(x, acc, a9);
//     acc = T::polyeval_mla(x, acc, a8);
//     acc = T::polyeval_mla(x, acc, a7);
//     acc = T::polyeval_mla(x, acc, a6);
//     acc = T::polyeval_mla(x, acc, a5);
//     acc = T::polyeval_mla(x, acc, a4);
//     acc = T::polyeval_mla(x, acc, a3);
//     acc = T::polyeval_mla(x, acc, a2);
//     acc = T::polyeval_mla(x, acc, a1);
//     T::polyeval_mla(x, acc, a0)
// }

#[allow(clippy::too_many_arguments)]
pub(crate) fn f_horner_polyeval18<T: PolyevalMla + Copy>(
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
    let mut acc = a17;
    acc = T::polyeval_mla(x, acc, a16);
    acc = T::polyeval_mla(x, acc, a15);
    acc = T::polyeval_mla(x, acc, a14);
    acc = T::polyeval_mla(x, acc, a13);
    acc = T::polyeval_mla(x, acc, a12);
    acc = T::polyeval_mla(x, acc, a11);
    acc = T::polyeval_mla(x, acc, a10);
    acc = T::polyeval_mla(x, acc, a9);
    acc = T::polyeval_mla(x, acc, a8);
    acc = T::polyeval_mla(x, acc, a7);
    acc = T::polyeval_mla(x, acc, a6);
    acc = T::polyeval_mla(x, acc, a5);
    acc = T::polyeval_mla(x, acc, a4);
    acc = T::polyeval_mla(x, acc, a3);
    acc = T::polyeval_mla(x, acc, a2);
    acc = T::polyeval_mla(x, acc, a1);
    T::polyeval_mla(x, acc, a0)
}

#[inline(always)]
#[allow(clippy::too_many_arguments)]
pub(crate) fn f_horner_polyeval23<T: PolyevalMla + Copy>(
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
    a18: T,
    a19: T,
    a20: T,
    a21: T,
    a22: T,
) -> T {
    let mut acc = a22;
    acc = T::polyeval_mla(x, acc, a21);
    acc = T::polyeval_mla(x, acc, a20);
    acc = T::polyeval_mla(x, acc, a19);
    acc = T::polyeval_mla(x, acc, a18);
    acc = T::polyeval_mla(x, acc, a17);
    acc = T::polyeval_mla(x, acc, a16);
    acc = T::polyeval_mla(x, acc, a15);
    acc = T::polyeval_mla(x, acc, a14);
    acc = T::polyeval_mla(x, acc, a13);
    acc = T::polyeval_mla(x, acc, a12);
    acc = T::polyeval_mla(x, acc, a11);
    acc = T::polyeval_mla(x, acc, a10);
    acc = T::polyeval_mla(x, acc, a9);
    acc = T::polyeval_mla(x, acc, a8);
    acc = T::polyeval_mla(x, acc, a7);
    acc = T::polyeval_mla(x, acc, a6);
    acc = T::polyeval_mla(x, acc, a5);
    acc = T::polyeval_mla(x, acc, a4);
    acc = T::polyeval_mla(x, acc, a3);
    acc = T::polyeval_mla(x, acc, a2);
    acc = T::polyeval_mla(x, acc, a1);
    T::polyeval_mla(x, acc, a0)
}
