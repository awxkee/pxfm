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
use crate::bits::set_exponent_f64;
use crate::common::{dd_fmla, f_fmla};
use crate::double_double::DoubleDouble;
use crate::dyadic_float::{DyadicFloat128, DyadicSign};
use crate::sincos_reduce_tables::ONE_TWENTY_EIGHT_OVER_PI;

#[derive(Debug)]
pub(crate) struct AngleReduced {
    pub(crate) angle: DoubleDouble,
}

#[derive(Default)]
pub(crate) struct LargeArgumentReduction {
    x_reduced: f64,
    idx: u64,
    y_hi: f64,
    y_lo: f64,
    // Low part of x * ONE_TWENTY_EIGHT_OVER_PI[idx][1].
    y_mid: DoubleDouble,
}

// For large range |x| >= 2^16, we perform the range reduction computations as:
//   u = x - k * pi/128 = (pi/128) * (x * (128/pi) - k).
// We use the exponent of x to find 4 double-chunks of 128/pi:
// c_hi, c_mid, c_lo, c_lo_2 such that:
//   1) ulp(round(x * c_hi, D, RN)) >= 2^8 = 256,
//   2) If x * c_hi = ph_hi + ph_lo and x * c_mid = pm_hi + pm_lo, then
//        min(ulp(ph_lo), ulp(pm_hi)) >= 2^-53.
// This will allow us to drop the high part ph_hi and the addition:
//   (ph_lo + pm_hi) mod 1
// can be exactly representable in a double precision.
// This will allow us to do split the computations as:
//   (x * 256/pi) ~ x * (c_hi + c_mid + c_lo + c_lo_2)    (mod 256)
//                ~ (ph_lo + pm_hi) + (pm_lo + x * c_lo) + x * c_lo_2.
// Then,
//   round(x * 128/pi) = round(ph_lo + pm_hi)    (mod 256)
// And the high part of fractional part of (x * 128/pi) can simply be:
//   {x * 128/pi}_hi = {ph_lo + pm_hi}.
// To prevent overflow when x is very large, we simply scale up
// (c_hi, c_mid, c_lo, c_lo_2) by a fixed power of 2 (based on the index) and
// scale down x by the same amount.
impl LargeArgumentReduction {
    #[cold]
    pub(crate) fn accurate(&self) -> DyadicFloat128 {
        // Sage math:
        // R = RealField(128)
        // π = R.pi()
        //
        // def format_hex(value):
        //     l = hex(value)[2:]
        //     n = 4
        //     x = [l[i:i + n] for i in range(0, len(l), n)]
        //     return "0x" + "_".join(x) + "_u128"
        //
        // def print_dyadic(value):
        //     (s, m, e) = RealField(128)(value).sign_mantissa_exponent();
        //     print("DyadicFloat128 {")
        //     print(f"    sign: DyadicSign::{'Pos' if s >= 0 else 'Neg'},")
        //     print(f"    exponent: {e},")
        //     print(f"    mantissa: {format_hex(m)},")
        //     print("};")
        //
        // print_dyadic(π/128)
        const PI_OVER_128_F128: DyadicFloat128 = DyadicFloat128 {
            sign: DyadicSign::Pos,
            exponent: -133,
            mantissa: 0xc90f_daa2_2168_c234_c4c6_628b_80dc_1cd1_u128,
        };

        // y_lo = x * c_lo + pm.lo
        let one_pi_rot = ONE_TWENTY_EIGHT_OVER_PI[self.idx as usize];
        let y_lo_0 = DyadicFloat128::new_from_f64(self.x_reduced * f64::from_bits(one_pi_rot.3));
        let y_lo_1 = DyadicFloat128::new_from_f64(self.y_lo) + y_lo_0;
        let y_mid_f128 = DyadicFloat128::new_from_f64(self.y_mid.lo) + y_lo_1;
        let y_hi_f128 =
            DyadicFloat128::new_from_f64(self.y_hi) + DyadicFloat128::new_from_f64(self.y_mid.hi);
        let y = y_hi_f128 + y_mid_f128;

        y * PI_OVER_128_F128
    }

    pub(crate) fn reduce(&mut self, x: f64) -> (u64, DoubleDouble) {
        const E_BIAS: u64 = (1u64 << (11 - 1u64)) - 1u64;
        let mut xbits = x.to_bits();
        let x_e = ((x.to_bits() >> 52) & 0x7ff) as i64;
        let x_e_m62 = x_e.wrapping_sub(E_BIAS as i64 + 62);
        self.idx = (x_e_m62 >> 4).wrapping_add(3) as u64;
        // Scale x down by 2^(-(16 * (idx - 3))
        xbits = set_exponent_f64(
            xbits,
            (x_e_m62 & 15)
                .wrapping_add(E_BIAS as i64)
                .wrapping_add(62i64) as u64,
        );
        // 2^62 <= |x_reduced| < 2^(62 + 16) = 2^78
        self.x_reduced = f64::from_bits(xbits);
        // x * c_hi = ph.hi + ph.lo exactly.
        let one_pi = ONE_TWENTY_EIGHT_OVER_PI[self.idx as usize];
        let ph = DoubleDouble::from_exact_mult(self.x_reduced, f64::from_bits(one_pi.0));
        // x * c_mid = pm.hi + pm.lo exactly.
        let pm = DoubleDouble::from_exact_mult(self.x_reduced, f64::from_bits(one_pi.1));
        // x * c_lo = pl.hi + pl.lo exactly.
        let pl = DoubleDouble::from_exact_mult(self.x_reduced, f64::from_bits(one_pi.2));
        // Extract integral parts and fractional parts of (ph.lo + pm.hi).
        let sum_hi = ph.lo + pm.hi;
        let kd = sum_hi.round();

        // x * 128/pi mod 1 ~ y_hi + y_mid + y_lo
        self.y_hi = (ph.lo - kd) + pm.hi; // Exact
        self.y_mid = DoubleDouble::from_exact_add(pm.lo, pl.hi);
        self.y_lo = pl.lo;

        // y_l = x * c_lo_2 + pl.lo
        let y_l = dd_fmla(self.x_reduced, f64::from_bits(one_pi.3), self.y_lo);
        let mut y = DoubleDouble::from_exact_add(self.y_hi, self.y_mid.hi);
        y.lo += self.y_mid.lo + y_l;

        // Digits of pi/128, generated by SageMath with:
        // import struct
        // from sage.all import *
        //
        // def double_to_hex(f):
        //     return "0x" + format(struct.unpack('<Q', struct.pack('<d', f))[0], '016x')
        //
        // R = RealField(128)
        // π = R.pi()
        //
        // RN = RealField(53)
        //
        // hi = RN(π/128)
        // lo = RN(π/128 - R(hi))
        //
        // print("lo: " + double_to_hex(lo))
        // print("hi: " + double_to_hex(hi))
        const PI_OVER_128_DD: DoubleDouble = DoubleDouble::new(
            f64::from_bits(0x3c31a62633145c07),
            f64::from_bits(0x3f9921fb54442d18),
        );

        // Error bound: with {a} denote the fractional part of a, i.e.:
        //   {a} = a - round(a)
        // Then,
        //   | {x * 128/pi} - (y_hi + y_lo) | <=  ulp(ulp(y_hi)) <= 2^-105
        //   | {x mod pi/128} - (u.hi + u.lo) | < 2 * 2^-6 * 2^-105 = 2^-110
        let u = DoubleDouble::quick_mult(y, PI_OVER_128_DD);

        ((kd as i64) as u64, u)
    }
}

static INVPI_2_62: [u64; 20] = [
    0x28be60db9391054a,
    0x7f09d5f47d4d3770,
    0x36d8a5664f10e410,
    0x7f9458eaf7aef158,
    0x6dc91b8e909374b8,
    0x1924bba82746487,
    0x3f877ac72c4a69cf,
    0xba208d7d4baed121,
    0x3a671c09ad17df90,
    0x4e64758e60d4ce7d,
    0x272117e2ef7e4a0e,
    0xc7fe25fff7816603,
    0xfbcbc462d6829b47,
    0xdb4d9fb3c9f2c26d,
    0xd3d18fd9a797fa8b,
    0x5d49eeb1faf97c5e,
    0xcf41ce7de294a4ba,
    0x9afed7ec47e35742,
    0x1580cc11bf1edaea,
    0xfc33ef0826bd0d87,
];

#[inline]
fn create_dd(c1: u64, c0: u64) -> DoubleDouble {
    let mut c1 = c1;
    let mut c0 = c0;
    if c1 != 0 {
        let e = c1.leading_zeros();
        if e != 0 {
            c1 = (c1 << e) | (c0 >> (64 - e));
            c0 = c0.wrapping_shl(e);
        }
        let f = 0x3fe - e;
        let t_u = ((f as u64) << 52) | ((c1 << 1) >> 12);
        let hi = f64::from_bits(t_u);
        c0 = (c1 << 53) | (c0 >> 11);
        let l = if (c0) != 0 {
            let g = c0.leading_zeros();
            if (g) != 0 {
                c0 = c0.wrapping_shl(g);
            }
            let t_u = (((f - 53 - g) as u64) << 52) | ((c0 << 1) >> 12);
            f64::from_bits(t_u)
        } else {
            0.
        };
        DoubleDouble::new(l, hi)
    } else if (c0) != 0 {
        let e = c0.leading_zeros();
        let f = 0x3fe - 64 - e;
        c0 = c0.wrapping_shl(e + 1); // most significant bit shifted out

        /* put the upper 52 bits of c0 into h */
        let t_u = ((f as u64) << 52) | (c0 >> 12);
        let hi = f64::from_bits(t_u);
        /* put the lower 12 bits of c0 into l */
        c0 = c0.wrapping_shl(52);
        let l = if (c0) != 0 {
            let g = c0.leading_zeros();
            c0 = c0.wrapping_shl(g + 1);
            let t_u = (((f - 64 - g) as u64) << 52) | (c0 >> 12);
            f64::from_bits(t_u)
        } else {
            0.
        };
        DoubleDouble::new(l, hi)
    } else {
        DoubleDouble::default()
    }
}

#[inline]
fn frac_2_pi(x: f64) -> DoubleDouble {
    if x <= f64::from_bits(0x401921fb54442d17)
    // x < 2*pi
    {
        /* | CH+CL - 1/(2pi) | < 2^-110.523 */
        const C: DoubleDouble = DoubleDouble::new(
            f64::from_bits(0xbc66b01ec5417056),
            f64::from_bits(0x3fc45f306dc9c883),
        );
        let mut z = DoubleDouble::quick_mult_f64(C, x);
        z.lo = f_fmla(C.lo, x, z.lo);
        z
    } else
    // x > 0x1.921fb54442d17p+2
    {
        let t = x.to_bits();
        let mut e = ((t >> 52) & 0x7ff) as i32; /* 1025 <= e <= 2046 */
        let m = (1u64 << 52) | (t & 0xfffffffffffffu64);
        let mut c0: u64;
        let mut c1: u64;
        let mut c2: u64;
        // x = m/2^53 * 2^(e-1022)
        if e <= 1074
        // 1025 <= e <= 1074: 2^2 <= x < 2^52
        {
            let mut u = m as u128 * INVPI_2_62[1] as u128;
            c0 = u as u64;
            c1 = (u >> 64) as u64;
            u = m as u128 * INVPI_2_62[0] as u128;
            c1 = c1.wrapping_add(u as u64);
            c2 = (u >> 64) as u64 + (c1 < (u as u64)) as u64;
            e = 1075 - e; // 1 <= e <= 50
        } else
        // 1075 <= e <= 2046, 2^52 <= x < 2^1024
        {
            let i = (e - 1138 + 63) / 64; // i = ceil((e-1138)/64), 0 <= i <= 15
            let mut u = m as u128 * INVPI_2_62[i as usize + 2] as u128;
            c0 = u as u64;
            c1 = (u >> 64) as u64;
            u = m as u128 * INVPI_2_62[i as usize + 1] as u128;
            c1 = c1.wrapping_add(u as u64);
            c2 = (u >> 64) as u64 + ((c1) < (u as u64)) as u64;
            u = m as u128 * INVPI_2_62[i as usize] as u128;
            c2 = c2.wrapping_add(u as u64);
            e = 1139 + (i << 6) - e; // 1 <= e <= 64
        }
        if e == 64 {
            c0 = c1;
            c1 = c2;
        } else {
            c0 = (c1 << (64 - e)) | c0 >> e;
            c1 = (c2 << (64 - e)) | c1 >> e;
        }
        create_dd(c1, c0)
    }
}

/// Returns x mod 2*PI
pub(crate) fn rem2pi_any(x: f64) -> AngleReduced {
    const TWO_PI: DoubleDouble = DoubleDouble::new(
        f64::from_bits(0x3cb1a62633145c07),
        f64::from_bits(0x401921fb54442d18),
    );
    let a = frac_2_pi(x);
    let z = DoubleDouble::quick_mult(a, TWO_PI);
    AngleReduced { angle: z }
}

#[inline]
pub(crate) fn rem2pif_any(x: f32) -> f64 {
    let x_abs = x.abs();
    let p = rem2pi_any(x_abs as f64);
    p.angle.to_f64()
}
