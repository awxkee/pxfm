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

#![allow(clippy::excessive_precision)]

use crate::bessel::i0::bessel_rsqrt_hard;
use crate::bessel::j1_coeffs::{
    J1_COEFFS, J1_COEFFS_RATIONAL128, J1_ZEROS, J1_ZEROS_RATIONAL, J1_ZEROS_VALUE,
};
use crate::double_double::DoubleDouble;
use crate::dyadic_float::{DyadicFloat128, DyadicSign};
use crate::horner::{f_horner_polyeval12, f_horner_polyeval18};
use crate::polyeval::{f_polyeval8, f_polyeval9, f_polyeval13, f_polyeval24};
use crate::sin_helper::{sin_dd_small, sin_f128_small};
use crate::sincos_reduce::{AngleReduced, rem2pi_any, rem2pi_f128};

/// Bessel J 1st order in f64
///
/// Max found ULP 0.5.
///
/// Note about accuracy:
/// - Close to zero Bessel have tiny values such that testing against MPFR must be done exactly
///   in the same precision, since any nearest representable number have ULP > 0.5,
///   for example `J1(0.000000000000000000000000000000000000023509886)` in single precision
///   have 0.7 ULP for any number with extended precision that would be represented in f32
///   Same applies to J1(4.4501477170144018E-309) in double precision and some others subnormal numbers
pub fn f_j1(x: f64) -> f64 {
    if !x.is_normal() {
        if x == 0. {
            return x;
        }
        if x.is_infinite() {
            return 0.;
        }
        if x.is_nan() {
            return x + x;
        }
    }

    let ax: u64 = x.to_bits() & 0x7fff_ffff_ffff_ffff;

    if ax < 0x4052a6784230fcf8u64 {
        // 74.60109
        if ax < 0x3fd0000000000000u64 {
            // 0.25
            return j1_maclaurin_series(x);
        }
        return j1_small_argument_path(x);
    }

    j1_asympt(x)
}

/*
   Evaluates:
   J1 = sqrt(2/(PI*x)) * beta(x) * cos(x - 3*PI/4 - alpha(x))
   discarding 1*PI/2 using identities gives:
   J1 = sqrt(2/(PI*x)) * beta(x) * sin(x - PI/4 - alpha(x))

   to avoid squashing small (-PI/4 - alpha(x)) into a large x actual expansion is:

   J1 = sqrt(2/(PI*x)) * beta(x) * sin((x mod 2*PI) - PI/4 - alpha(x))
*/
#[inline]
fn j1_asympt(x: f64) -> f64 {
    let origin_x = x;
    static SGN: [f64; 2] = [1., -1.];
    let sign_scale = SGN[x.is_sign_negative() as usize];
    let x = x.abs();

    const SQRT_2_OVER_PI: DoubleDouble = DoubleDouble::new(
        f64::from_bits(0xbc8cbc0d30ebfd15),
        f64::from_bits(0x3fe9884533d43651),
    );
    const MPI_OVER_4: DoubleDouble = DoubleDouble::new(
        f64::from_bits(0xbc81a62633145c07),
        f64::from_bits(0xbfe921fb54442d18),
    );

    let recip = if x.to_bits() > 0x7fd000000000000u64 {
        DoubleDouble::quick_mult_f64(DoubleDouble::from_exact_safe_div(4.0, x), 0.25)
    } else {
        DoubleDouble::from_recip(x)
    };

    let alpha = j1_asympt_alpha(recip);
    let beta = j1_asympt_beta(recip);

    let AngleReduced { angle } = rem2pi_any(x);

    // Without full subtraction cancellation happens sometimes
    let x0pi34 = DoubleDouble::dd_sub(MPI_OVER_4, alpha);
    let r0 = DoubleDouble::dd_add(angle, x0pi34);

    let m_sin = sin_dd_small(r0);
    let z0 = DoubleDouble::quick_mult(beta, m_sin);
    let r_sqrt = DoubleDouble::from_rsqrt(x);
    let scale = DoubleDouble::quick_mult(SQRT_2_OVER_PI, r_sqrt);
    let r = DoubleDouble::quick_mult(scale, z0);

    const ERR: f64 = f64::from_bits(0x39d0000000000000);

    let ub = r.hi + (r.lo + ERR);
    let lb = r.hi + (r.lo - ERR);

    if ub == lb {
        return r.to_f64() * sign_scale;
    }

    j1_asympt_hard(origin_x)
}

/**
Note expansion generation below: this is negative series expressed in Sage as positive,
so before any real evaluation `x=1/x` should be applied.

Generated by SageMath:
```python
def binomial_like(n, m):
    prod = QQ(1)
    z = QQ(4)*(n**2)
    for k in range(1,m + 1):
        prod *= (z - (2*k - 1)**2)
    return prod / (QQ(2)**(2*m) * (ZZ(m).factorial()))

R = LaurentSeriesRing(RealField(300), 'x',default_prec=300)
x = R.gen()

def Pn_asymptotic(n, y, terms=10):
    # now y = 1/x
    return sum( (-1)**m * binomial_like(n, 2*m) / (QQ(2)**(2*m)) * y**(QQ(2)*m) for m in range(terms) )

def Qn_asymptotic(n, y, terms=10):
    return sum( (-1)**m * binomial_like(n, 2*m + 1) / (QQ(2)**(2*m + 1)) * y**(QQ(2)*m + 1) for m in range(terms) )

P = Pn_asymptotic(1, x, 50)
Q = Qn_asymptotic(1, x, 50)

R_series = (-Q/P)

# alpha is atan(R_series) so we're doing Taylor series atan expansion on R_series

arctan_series_Z = sum([QQ(-1)**k * x**(QQ(2)*k+1) / RealField(700)(RealField(700)(2)*k+1) for k in range(25)])
alpha_series = arctan_series_Z(R_series)

# see the series
print(alpha_series)
```

See notes/bessel_asympt.ipynb for generation
**/
#[inline]
pub(crate) fn j1_asympt_alpha(recip: DoubleDouble) -> DoubleDouble {
    const C: [(u64, u64); 12] = [
        (0x0000000000000000, 0xbfd8000000000000),
        (0x0000000000000000, 0x3fc5000000000000),
        (0x3c6999999999999a, 0xbfd7bccccccccccd),
        (0x3cab6db6db6db6db, 0x4002f486db6db6db),
        (0x0000000000000000, 0xc03e9fbf40000000),
        (0x3d21745d1745d174, 0x4084997b55945d17),
        (0x3d789d89d89d89d9, 0xc0d4a914195269d9),
        (0xbdb999999999999a, 0x412cd1b53816aec1),
        (0xbdfe5a5a5a5a5a5a, 0xc18aa4095d419351),
        (0x3e7e0ca50d79435e, 0x41ef809305f11b9d),
        (0xbedff8b720000000, 0xc2572e6809ed618b),
        (0xbf64e5d8ae68b7a7, 0x42c4c5b6057839f9),
    ];

    // Doing (1/x)*(1/x) instead (1/(x*x)) to avoid spurious overflow/underflow
    let x2 = DoubleDouble::quick_mult(recip, recip);

    let mut p = DoubleDouble::mul_add(
        x2,
        DoubleDouble::from_bit_pair(C[11]),
        DoubleDouble::from_bit_pair(C[10]),
    );

    p = DoubleDouble::mul_add(x2, p, DoubleDouble::from_bit_pair(C[9]));
    p = DoubleDouble::mul_add(x2, p, DoubleDouble::from_bit_pair(C[8]));
    p = DoubleDouble::mul_add(x2, p, DoubleDouble::from_bit_pair(C[7]));
    p = DoubleDouble::mul_add(x2, p, DoubleDouble::from_bit_pair(C[6]));
    p = DoubleDouble::mul_add(x2, p, DoubleDouble::from_bit_pair(C[5]));
    p = DoubleDouble::mul_add(x2, p, DoubleDouble::from_bit_pair(C[4]));
    p = DoubleDouble::mul_add_f64(x2, p, f64::from_bits(C[3].1));
    p = DoubleDouble::mul_add(x2, p, DoubleDouble::from_bit_pair(C[2]));
    p = DoubleDouble::mul_add_f64(x2, p, f64::from_bits(C[1].1));
    p = DoubleDouble::mul_add_f64(x2, p, f64::from_bits(C[0].1));

    let z = DoubleDouble::quick_mult(p, recip);

    DoubleDouble::from_exact_add(z.hi, z.lo)
}

/**
Note expansion generation below: this is negative series expressed in Sage as positive,
so before any real evaluation `x=1/x` should be applied

Generated by SageMath:
```python
def binomial_like(n, m):
    prod = QQ(1)
    z = QQ(4)*(n**2)
    for k in range(1,m + 1):
        prod *= (z - (2*k - 1)**2)
    return prod / (QQ(2)**(2*m) * (ZZ(m).factorial()))

R = LaurentSeriesRing(RealField(300), 'x',default_prec=300)
x = R.gen()

def Pn_asymptotic(n, y, terms=10):
    # now y = 1/x
    return sum( (-1)**m * binomial_like(n, 2*m) / (QQ(2)**(2*m)) * y**(QQ(2)*m) for m in range(terms) )

def Qn_asymptotic(n, y, terms=10):
    return sum( (-1)**m * binomial_like(n, 2*m + 1) / (QQ(2)**(2*m + 1)) * y**(QQ(2)*m + 1) for m in range(terms) )

P = Pn_asymptotic(1, x, 50)
Q = Qn_asymptotic(1, x, 50)

def sqrt_series(s):
    val = S.valuation()
    lc = S[val]  # Leading coefficient
    b = lc.sqrt() * x**(val // 2)

    for _ in range(5):
        b = (b + S / b) / 2
        b = b
    return b

S = (P**2 + Q**2).truncate(50)

b_series = sqrt_series(S).truncate(30)
# see the beta series
print(b_series)
```

See notes/bessel_asympt.ipynb for generation
**/
#[inline]
pub(crate) fn j1_asympt_beta(recip: DoubleDouble) -> DoubleDouble {
    const C: [(u64, u64); 10] = [
        (0x0000000000000000, 0x3ff0000000000000), // 1
        (0x0000000000000000, 0x3fc8000000000000), // 2
        (0x0000000000000000, 0xbfc8c00000000000), // 3
        (0x0000000000000000, 0x3fe9c50000000000), // 4
        (0x0000000000000000, 0xc01ef5b680000000), // 5
        (0x0000000000000000, 0x40609860dd400000), // 6
        (0x0000000000000000, 0xc0abae9b7a06e000), // 7
        (0x0000000000000000, 0x41008711d41c1428), // 8
        (0xbdf7a00000000000, 0xc15ab70164c8be6e),
        (0xbe40e1f000000000, 0x41bc1055e24f297f),
    ];

    // Doing (1/x)*(1/x) instead (1/(x*x)) to avoid spurious overflow/underflow
    let x2 = DoubleDouble::quick_mult(recip, recip);

    let mut p = DoubleDouble::mul_add(
        x2,
        DoubleDouble::from_bit_pair(C[9]),
        DoubleDouble::from_bit_pair(C[8]),
    );

    p = DoubleDouble::mul_add_f64(x2, p, f64::from_bits(C[7].1)); // 8
    p = DoubleDouble::mul_add_f64(x2, p, f64::from_bits(C[6].1)); // 7
    p = DoubleDouble::mul_add_f64(x2, p, f64::from_bits(C[5].1)); // 6
    p = DoubleDouble::mul_add_f64(x2, p, f64::from_bits(C[4].1)); // 5
    p = DoubleDouble::mul_add_f64(x2, p, f64::from_bits(C[3].1)); // 4
    p = DoubleDouble::mul_add_f64(x2, p, f64::from_bits(C[2].1)); // 3
    p = DoubleDouble::mul_add_f64(x2, p, f64::from_bits(C[1].1)); // 2
    p = DoubleDouble::mul_add_f64(x2, p, f64::from_bits(C[0].1)); // 1
    p
}

/**
Generated in Sollya:
```text
pretty = proc(u) {
  return ~(floor(u*1000)/1000);
};

bessel_j1 = library("./cmake-build-release/libbessel_sollya.dylib");

f = bessel_j1(x)/x;
d = [0, 0.921];
w = 1;
pf = fpminimax(f, [|0,2,4,6,8,10,12,14,16,18,20,22,24|], [|107, 107, 107, 107, 107, D...|], d, absolute, floating);

w = 1;
or_f = bessel_j1(x);
pf1 = pf * x;
err_p = -log2(dirtyinfnorm(pf1*w-or_f, d));
print ("relative error:", pretty(err_p));

for i from 0 to degree(pf) by 2 do {
    print("'", coeff(pf, i), "',");
};
```
See ./notes/bessel_sollya/bessel_j1_at_zero.sollya
**/
#[inline]
pub(crate) fn j1_maclaurin_series(x: f64) -> f64 {
    const CL: [(u64, u64); 5] = [
        (0xb930000000000000, 0x3fe0000000000000),
        (0x39c8e80000000000, 0xbfb0000000000000),
        (0x3c05555554f3add7, 0x3f65555555555555),
        (0xbbac71c4eb0f8c94, 0xbf0c71c71c71c71c),
        (0xbb3f56b7a43206d4, 0x3ea6c16c16c16c17),
    ];

    let p = f_polyeval8(
        x * x,
        f64::from_bits(0xbe3845c8a0ce5129),
        f64::from_bits(0x3dc27e4fb7789ea2),
        f64::from_bits(0xbd4522a43f633af1),
        f64::from_bits(0x3cc2c97589d53f97),
        f64::from_bits(0xbc3ab8151dca7912),
        f64::from_bits(0x3baf08732286d1d4),
        f64::from_bits(0xbb10ac65637413f4),
        f64::from_bits(0xbae4d8336e4f779c),
    );

    let dx2 = DoubleDouble::from_exact_mult(x, x);

    let mut p_e = DoubleDouble::mul_f64_add(dx2, p, DoubleDouble::from_bit_pair(CL[4]));
    p_e = DoubleDouble::mul_add(dx2, p_e, DoubleDouble::from_bit_pair(CL[3]));
    p_e = DoubleDouble::mul_add(dx2, p_e, DoubleDouble::from_bit_pair(CL[2]));
    p_e = DoubleDouble::mul_add(dx2, p_e, DoubleDouble::from_bit_pair(CL[1]));
    p_e = DoubleDouble::mul_add(dx2, p_e, DoubleDouble::from_bit_pair(CL[0]));

    let px = DoubleDouble::quick_mult_f64(p_e, x);

    const ERR: f64 = f64::from_bits(0x395b1889a0146d8a);
    let ub = px.hi + (px.lo + ERR);
    let lb = px.hi + (px.lo - ERR);
    if ub != lb {
        return j1_maclaurin_series_hard(x);
    }

    px.to_f64()
}

/**
Taylor expansion at 0

Generated by SageMath:
```python
def print_expansion_at_0():
    print(f"static C: [DyadicFloat128; 13] = ")
    from mpmath import mp, j1, taylor, expm1
    poly = taylor(lambda val: j1(val), 0, 26)
    real_i = 0
    print("[")
    for i in range(1, len(poly), 2):
        print_dyadic(poly[i])
        real_i = real_i + 1
    print("],")

    print("];")

mp.prec = 180

print_expansion_at_0()
```
**/
#[cold]
#[inline(never)]
fn j1_maclaurin_series_hard(x: f64) -> f64 {
    static C: [DyadicFloat128; 13] = [
        DyadicFloat128 {
            sign: DyadicSign::Pos,
            exponent: -128,
            mantissa: 0x80000000_00000000_00000000_00000000_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Neg,
            exponent: -131,
            mantissa: 0x80000000_00000000_00000000_00000000_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Pos,
            exponent: -136,
            mantissa: 0xaaaaaaaa_aaaaaaaa_aaaaaaaa_aaaaaaab_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Neg,
            exponent: -142,
            mantissa: 0xe38e38e3_8e38e38e_38e38e38_e38e38e4_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Pos,
            exponent: -148,
            mantissa: 0xb60b60b6_0b60b60b_60b60b60_b60b60b6_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Neg,
            exponent: -155,
            mantissa: 0xc22e4506_72894ab6_cd8efb11_d33f5618_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Pos,
            exponent: -162,
            mantissa: 0x93f27dbb_c4fae397_780b69f5_333c725b_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Neg,
            exponent: -170,
            mantissa: 0xa91521fb_2a434d3f_649f5485_f169a743_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Pos,
            exponent: -178,
            mantissa: 0x964bac6d_7ae67d8d_aec68405_485dea03_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Neg,
            exponent: -187,
            mantissa: 0xd5c0f53a_fe6fa17f_8c7b0b68_39691f4e_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Pos,
            exponent: -196,
            mantissa: 0xf8bb4be7_8e7896b0_58fee362_01a4370c_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Neg,
            exponent: -205,
            mantissa: 0xf131bdf7_cff8d02e_e1ef6820_f9d58ab6_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Pos,
            exponent: -214,
            mantissa: 0xc5e72c48_0d1aec75_3caa2e0d_edd008ca_u128,
        },
    ];

    let rx = DyadicFloat128::new_from_f64(x);
    let dx = rx * rx;
    let p = f_polyeval13(
        dx, C[0], C[1], C[2], C[3], C[4], C[5], C[6], C[7], C[8], C[9], C[10], C[11], C[12],
    );
    (p * rx).fast_as_f64()
}

/// This method on small range searches for nearest zero or extremum.
/// Then picks stored series expansion at the point end evaluates the poly at the point.
#[inline]
pub(crate) fn j1_small_argument_path(x: f64) -> f64 {
    static SIGN: [f64; 2] = [1., -1.];
    let sign_scale = SIGN[x.is_sign_negative() as usize];
    let x_abs = f64::from_bits(x.to_bits() & 0x7fff_ffff_ffff_ffff);

    // let avg_step = 74.60109 / 47.0;
    // let inv_step = 1.0 / avg_step;

    const INV_STEP: f64 = 0.6300176043004198;

    let fx = x_abs * INV_STEP;
    const J1_ZEROS_COUNT: f64 = (J1_ZEROS.len() - 1) as f64;
    let idx0 = fx.min(J1_ZEROS_COUNT) as usize;
    let idx1 = fx.ceil().min(J1_ZEROS_COUNT) as usize;

    let found_zero0 = DoubleDouble::from_bit_pair(J1_ZEROS[idx0]);
    let found_zero1 = DoubleDouble::from_bit_pair(J1_ZEROS[idx1]);

    let dist0 = (found_zero0.hi - x_abs).abs();
    let dist1 = (found_zero1.hi - x_abs).abs();

    let (found_zero, idx, dist) = if dist0 < dist1 {
        (found_zero0, idx0, dist0)
    } else {
        (found_zero1, idx1, dist1)
    };

    if idx == 0 {
        return j1_maclaurin_series(x);
    }

    let r = DoubleDouble::full_add_f64(DoubleDouble::new(-found_zero.lo, -found_zero.hi), x_abs);

    // We hit exact zero, value, better to return it directly
    if dist == 0. {
        return f64::from_bits(J1_ZEROS_VALUE[idx]) * sign_scale;
    }

    let j1c = &J1_COEFFS[idx - 1];
    let c0 = j1c;

    let c = &c0[15..];

    let p0 = f_polyeval9(
        r.to_f64(),
        f64::from_bits(c[0].1),
        f64::from_bits(c[1].1),
        f64::from_bits(c[2].1),
        f64::from_bits(c[3].1),
        f64::from_bits(c[4].1),
        f64::from_bits(c[5].1),
        f64::from_bits(c[6].1),
        f64::from_bits(c[7].1),
        f64::from_bits(c[8].1),
    );

    let c = c0;

    let mut p_e = DoubleDouble::mul_f64_add(r, p0, DoubleDouble::from_bit_pair(c[14]));
    p_e = DoubleDouble::mul_add(p_e, r, DoubleDouble::from_bit_pair(c[13]));
    p_e = DoubleDouble::mul_add(p_e, r, DoubleDouble::from_bit_pair(c[12]));
    p_e = DoubleDouble::mul_add(p_e, r, DoubleDouble::from_bit_pair(c[11]));
    p_e = DoubleDouble::mul_add(p_e, r, DoubleDouble::from_bit_pair(c[10]));
    p_e = DoubleDouble::mul_add(p_e, r, DoubleDouble::from_bit_pair(c[9]));
    p_e = DoubleDouble::mul_add(p_e, r, DoubleDouble::from_bit_pair(c[8]));
    p_e = DoubleDouble::mul_add(p_e, r, DoubleDouble::from_bit_pair(c[7]));
    p_e = DoubleDouble::mul_add(p_e, r, DoubleDouble::from_bit_pair(c[6]));
    p_e = DoubleDouble::mul_add(p_e, r, DoubleDouble::from_bit_pair(c[5]));
    p_e = DoubleDouble::mul_add(p_e, r, DoubleDouble::from_bit_pair(c[4]));
    p_e = DoubleDouble::mul_add(p_e, r, DoubleDouble::from_bit_pair(c[3]));
    p_e = DoubleDouble::mul_add(p_e, r, DoubleDouble::from_bit_pair(c[2]));
    p_e = DoubleDouble::mul_add(p_e, r, DoubleDouble::from_bit_pair(c[1]));
    p_e = DoubleDouble::mul_add(p_e, r, DoubleDouble::from_bit_pair(c[0]));

    let p = DoubleDouble::from_full_exact_add(p_e.hi, p_e.lo);

    const ERR: f64 = f64::from_bits(0x3990000000000000);
    let ub = p.hi + (p.lo + ERR);
    let lb = p.hi + (p.lo - ERR);
    if ub != lb {
        return j1_small_argument_path_hard(x, idx, sign_scale);
    }
    p.to_f64() * sign_scale
}

#[cold]
#[inline(never)]
fn j1_small_argument_path_hard(x: f64, idx: usize, sign_scale: f64) -> f64 {
    let c = &J1_COEFFS_RATIONAL128[idx - 1];
    let zero = J1_ZEROS_RATIONAL[idx];
    let dx = DyadicFloat128::new_from_f64(x) - zero;

    let p = f_polyeval24(
        dx, c[0], c[1], c[2], c[3], c[4], c[5], c[6], c[7], c[8], c[9], c[10], c[11], c[12], c[13],
        c[14], c[15], c[16], c[17], c[18], c[19], c[20], c[21], c[22], c[23],
    );
    p.fast_as_f64() * sign_scale
}

/// see [j1_asympt_beta] for more info
pub(crate) fn j1_asympt_beta_hard(recip: DyadicFloat128) -> DyadicFloat128 {
    const C: [DyadicFloat128; 12] = [
        DyadicFloat128 {
            sign: DyadicSign::Pos,
            exponent: -127,
            mantissa: 0x80000000_00000000_00000000_00000000_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Pos,
            exponent: -130,
            mantissa: 0xc0000000_00000000_00000000_00000000_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Neg,
            exponent: -130,
            mantissa: 0xc6000000_00000000_00000000_00000000_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Pos,
            exponent: -128,
            mantissa: 0xce280000_00000000_00000000_00000000_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Neg,
            exponent: -125,
            mantissa: 0xf7adb400_00000000_00000000_00000000_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Pos,
            exponent: -120,
            mantissa: 0x84c306ea_00000000_00000000_00000000_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Neg,
            exponent: -116,
            mantissa: 0xdd74dbd0_37000000_00000000_00000000_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Pos,
            exponent: -110,
            mantissa: 0x84388ea0_e0a14000_00000000_00000000_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Neg,
            exponent: -105,
            mantissa: 0xd5b80b26_45f372f4_00000000_00000000_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Pos,
            exponent: -99,
            mantissa: 0xe082af12_794bf6f1_e1000000_00000000_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Neg,
            exponent: -92,
            mantissa: 0x94a06149_f30146bc_fe8ed000_00000000_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Pos,
            exponent: -86,
            mantissa: 0xf212edfc_42a62526_4fac2b0c_00000000_u128,
        },
    ];

    let x2 = recip * recip;

    f_horner_polyeval12(
        x2, C[0], C[1], C[2], C[3], C[4], C[5], C[6], C[7], C[8], C[9], C[10], C[11],
    )
}
//
/// See [j1_asympt_alpha] for the info
pub(crate) fn j1_asympt_alpha_hard(reciprocal: DyadicFloat128) -> DyadicFloat128 {
    const C: [DyadicFloat128; 18] = [
        DyadicFloat128 {
            sign: DyadicSign::Neg,
            exponent: -129,
            mantissa: 0xc0000000_00000000_00000000_00000000_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Pos,
            exponent: -130,
            mantissa: 0xa8000000_00000000_00000000_00000000_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Neg,
            exponent: -129,
            mantissa: 0xbde66666_66666666_66666666_66666666_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Pos,
            exponent: -126,
            mantissa: 0x97a436db_6db6db6d_b6db6db6_db6db6db_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Neg,
            exponent: -123,
            mantissa: 0xf4fdfa00_00000000_00000000_00000000_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Pos,
            exponent: -118,
            mantissa: 0xa4cbdaac_a2e8ba2e_8ba2e8ba_2e8ba2e9_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Neg,
            exponent: -113,
            mantissa: 0xa548a0ca_934ec4ec_4ec4ec4e_c4ec4ec5_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Pos,
            exponent: -108,
            mantissa: 0xe68da9c0_b5760666_66666666_66666666_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Neg,
            exponent: -102,
            mantissa: 0xd5204aea_0c9a8879_69696969_69696969_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Pos,
            exponent: -96,
            mantissa: 0xfc04982f_88dce9e0_ca50d794_35e50d79_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Neg,
            exponent: -89,
            mantissa: 0xb973404f_6b0c58ff_c5b90000_00000000_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Pos,
            exponent: -82,
            mantissa: 0xa62db02b_c1cfc563_44ea32e9_0b21642d_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Neg,
            exponent: -75,
            mantissa: 0xb220e7ff_443c1584_7e85f4e0_55eb851f_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Pos,
            exponent: -68,
            mantissa: 0xe10a255c_ca5e68cc_00c2d6c0_acdc8000_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Neg,
            exponent: -60,
            mantissa: 0xa573790c_5186f23b_5db502ea_d9fa5432_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Pos,
            exponent: -52,
            mantissa: 0x8c0ffedc_407a7015_453df84e_9c3f1d39_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Neg,
            exponent: -44,
            mantissa: 0x874226ed_c298a17a_d8c49a4e_dc9281a5_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Pos,
            exponent: -36,
            mantissa: 0x93cab36c_9ab9495c_310fa9cd_4b065359_u128,
        },
    ];

    let x2 = reciprocal * reciprocal;

    let p = f_horner_polyeval18(
        x2, C[0], C[1], C[2], C[3], C[4], C[5], C[6], C[7], C[8], C[9], C[10], C[11], C[12], C[13],
        C[14], C[15], C[16], C[17],
    );

    p * reciprocal
}

/*
   Evaluates:
   J1 = sqrt(2/(PI*x)) * beta(x) * cos(x - 3*PI/4 - alpha(x))
   discarding 1*PI/2 using identities gives:
   J1 = sqrt(2/(PI*x)) * beta(x) * sin(x - PI/4 - alpha(x))

   to avoid squashing small (-PI/4 - alpha(x)) into a large x actual expansion is:

   J1 = sqrt(2/(PI*x)) * beta(x) * sin((x mod 2*PI) - PI/4 - alpha(x))

   This method is required for situations where x*x or 1/(x*x) will overflow
*/
#[cold]
#[inline(never)]
fn j1_asympt_hard(x: f64) -> f64 {
    static SGN: [f64; 2] = [1., -1.];
    let sign_scale = SGN[x.is_sign_negative() as usize];
    let x = x.abs();

    const SQRT_2_OVER_PI: DyadicFloat128 = DyadicFloat128 {
        sign: DyadicSign::Pos,
        exponent: -128,
        mantissa: 0xcc42299e_a1b28468_7e59e280_5d5c7180_u128,
    };

    const MPI_OVER_4: DyadicFloat128 = DyadicFloat128 {
        sign: DyadicSign::Neg,
        exponent: -128,
        mantissa: 0xc90fdaa2_2168c234_c4c6628b_80dc1cd1_u128,
    };

    let x_dyadic = DyadicFloat128::new_from_f64(x);
    let recip = DyadicFloat128::accurate_reciprocal(x);

    let alpha = j1_asympt_alpha_hard(recip);
    let beta = j1_asympt_beta_hard(recip);

    let angle = rem2pi_f128(x_dyadic);

    let x0pi34 = MPI_OVER_4 - alpha;
    let r0 = angle + x0pi34;

    let m_sin = sin_f128_small(r0);

    let z0 = beta * m_sin;
    let r_sqrt = bessel_rsqrt_hard(x, recip);
    let scale = SQRT_2_OVER_PI * r_sqrt;
    let p = scale * z0;
    p.fast_as_f64() * sign_scale
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_j1() {
        assert_eq!(
            f_j1(162605674999778540000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000.),
            0.000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000008686943178258183
        );
        assert_eq!(f_j1(-6.1795701510782757E+307), 8.130935041593236e-155);
        assert_eq!(
            f_j1(0.000000000000000000000000000000000000008827127),
            0.0000000000000000000000000000000000000044135635
        );
        assert_eq!(
            f_j1(-0.000000000000000000000000000000000000008827127),
            -0.0000000000000000000000000000000000000044135635
        );
        assert_eq!(f_j1(5.4), -0.3453447907795863);
        assert_eq!(
            f_j1(77.743162408196766932633181568235159),
            0.09049267898021947
        );
        assert_eq!(
            f_j1(84.027189586293545175976760219782591),
            0.0870430264022591
        );
        assert_eq!(f_j1(f64::NEG_INFINITY), 0.0);
        assert_eq!(f_j1(f64::INFINITY), 0.0);
        assert!(f_j1(f64::NAN).is_nan());
    }
}
