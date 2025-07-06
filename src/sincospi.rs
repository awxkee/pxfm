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
use crate::atan::poly_dd_3;
use crate::common::{dd_fmla, f_fmla};
use crate::dekker::Dekker;
use crate::sincospi_tables::{
    SINCOS_PI_CM2, SINCOS_PI_SM2, SINCOS_PI_SN2, SINPI_CM1, SINPI_SM1, SINPI_SN1,
};

#[inline]
pub(crate) fn poly_dd_4(x: Dekker, poly: [(u64, u64); 4], l: f64) -> Dekker {
    let zch = poly[3];
    let ach = f64::from_bits(zch.0) + l;
    let acl = (f64::from_bits(zch.0) - ach) + l + f64::from_bits(zch.1);
    let mut ch = Dekker::new(acl, ach);

    let zch = poly[2];
    ch = Dekker::mult(ch, x);
    let th = ch.hi + f64::from_bits(zch.0);
    let tl = (f64::from_bits(zch.0) - th) + ch.hi;
    ch.hi = th;
    ch.lo += tl + f64::from_bits(zch.1);

    let zch = poly[1];
    ch = Dekker::mult(ch, x);
    let th = ch.hi + f64::from_bits(zch.0);
    let tl = (f64::from_bits(zch.0) - th) + ch.hi;
    ch.hi = th;
    ch.lo += tl + f64::from_bits(zch.1);

    let zch = poly[0];
    ch = Dekker::mult(ch, x);
    let th = ch.hi + f64::from_bits(zch.0);
    let tl = (f64::from_bits(zch.0) - th) + ch.hi;
    ch.hi = th;
    ch.lo += tl + f64::from_bits(zch.1);

    ch
}

#[inline]
fn sincosn(s: i32) -> (Dekker, Dekker) {
    let mut j: i32 = s & 0x3ff;
    let it: i32 = -((s >> 10) & 1);
    j = (!it & j).wrapping_sub(it << 10).wrapping_sub(it & j);
    let is = j >> 5;
    let ic = 0x20 - is;
    let jm = j & 0x1f;
    let ss = (s >> 11) & 1;
    let sc = ((s.wrapping_add(1024)) >> 11) & 1;

    let dsb = Dekker::from_bit_pair(SINPI_SN1[is as usize]);
    let dcb = Dekker::from_bit_pair(SINPI_SN1[ic as usize]);
    let dsl = Dekker::from_bit_pair(SINPI_SM1[jm as usize]);
    let dcl = Dekker::from_bit_pair(SINPI_CM1[jm as usize]);

    let sb = dsb.to_f64();
    let cb = dcb.to_f64();
    let mut ch = f_fmla(dcb.hi, dcl.hi, -dsb.hi * dsl.hi);
    let mut cl = f_fmla(
        dcl.hi,
        dcb.lo,
        f_fmla(-dsl.hi, dsb.lo, cb * dcl.lo - sb * dsl.lo),
    );
    let mut sh = f_fmla(dsb.hi, dcl.hi, dcb.hi * dsl.hi);
    let mut sl = f_fmla(
        dsl.hi,
        dcb.lo,
        f_fmla(dcl.hi, dsb.lo, cb * dsl.lo + sb * dcl.lo),
    );

    let tch = ch + cl;
    let tcl = (ch - tch) + cl;
    let tsh = sh + sl;
    let tsl = (sh - tsh) + sl;

    let sign_c = if sc == 1 { -0.0 } else { 0.0 };
    let sign_s = if ss == 1 { -0.0 } else { 0.0 };
    ch = f64::copysign(1.0, sign_c) * tch;
    cl = f64::copysign(1.0, sign_c) * tcl;

    sh = f64::copysign(1.0, sign_s) * tsh;
    sl = f64::copysign(1.0, sign_s) * tsl;
    let t_sin = Dekker::new(sl, sh);
    let t_cos = Dekker::new(cl, ch);
    (t_sin, t_cos)
}

#[inline]
fn sincosn2(s: i32) -> (Dekker, Dekker) {
    let mut j: i32 = s & 0x3ff;
    let it: i32 = -((s >> 10) & 1);
    j = (!it & j).wrapping_sub(it << 10).wrapping_sub(it & j);
    let is = j >> 5;
    let ic = 0x20 - is;
    let jm = j & 0x1f;
    let ass = (s >> 11) & 1;
    let asc = ((s.wrapping_add(1024)) >> 11) & 1;

    let sb = Dekker::from_bit_pair(SINCOS_PI_SN2[is as usize]);
    let cb = Dekker::from_bit_pair(SINCOS_PI_SN2[ic as usize]);
    let sl = Dekker::from_bit_pair(SINCOS_PI_SM2[jm as usize]);
    let cl = Dekker::from_bit_pair(SINCOS_PI_CM2[jm as usize]);

    let cc = Dekker::quick_mult(cl, cb);
    let ss = Dekker::quick_mult(sl, sb);
    let cs = Dekker::quick_mult(cl, sb);
    let sc = Dekker::quick_mult(sl, cb);

    let tc = Dekker::add(ss, cc);
    let ts = Dekker::add(Dekker::new(-sc.lo, -sc.hi), cs);
    let mut tc2 = Dekker::add(cb, Dekker::new(-tc.lo, -tc.hi));
    let mut ts2 = Dekker::add(sb, Dekker::new(-ts.lo, -ts.hi));

    let sgb_c = if asc == 1 { -0.0 } else { 0.0 };
    tc2.hi *= f64::copysign(1.0, sgb_c);
    tc2.lo *= f64::copysign(1.0, sgb_c);
    let sgb_s = if ass == 1 { -0.0 } else { 0.0 };
    ts2.hi *= f64::copysign(1.0, sgb_s);
    ts2.lo *= f64::copysign(1.0, sgb_s);
    (ts2, tc2)
}

#[inline]
pub(crate) fn poly_dd_2(x: Dekker, poly: [(u64, u64); 2], l: f64) -> Dekker {
    let zch = poly[1];
    let ach = f64::from_bits(zch.0) + l;
    let acl = (f64::from_bits(zch.0) - ach) + l + f64::from_bits(zch.1);
    let mut ch = Dekker::new(acl, ach);

    let zch = poly[0];
    ch = Dekker::quick_mult(ch, x);
    let th = ch.hi + f64::from_bits(zch.0);
    let tl = (f64::from_bits(zch.0) - th) + ch.hi;
    ch.hi = th;
    ch.lo += tl + f64::from_bits(zch.1);
    ch
}

#[cold]
fn as_cospi_zero(x: f64) -> f64 {
    let x2 = x * x;
    let dx2 = f_fmla(x, x, -x2);
    const CH: [(u64, u64); 2] = [
        (0xc013bd3cc9be45de, 0xbcb692b71366cc04),
        (0x40103c1f081b5ac4, 0xbcb32b33fda9113c),
    ];
    const CL: [u64; 2] = [0xbff55d3c7e3cbff9, 0x3fce1f50604fa0ff];
    let fl = x2 * f_fmla(x2, f64::from_bits(CL[1]), f64::from_bits(CL[0]));
    let mut f = poly_dd_2(Dekker::new(dx2, x2), CH, fl);
    f = Dekker::mult(Dekker::new(dx2, x2), f);
    let y = Dekker::from_exact_add(1., f.hi);
    let mut y1 = Dekker::from_exact_add(y.lo, fl);
    let mut t = y1.hi.to_bits();
    if (t & 0x000fffffffffffff) == 0 {
        let w = y1.lo.to_bits();
        if ((w ^ t) >> 63) != 0 {
            t = t.wrapping_sub(1);
        } else {
            t = t.wrapping_add(1);
        }
        y1.hi = f64::from_bits(t);
    }
    y.hi + y1.hi
}

#[cold]
fn as_sinpi_zero(x: f64) -> f64 {
    let x2 = x * x;
    let dx2 = dd_fmla(x, x, -x2);
    let x3 = x2 * x;
    let dx3 = dd_fmla(x2, x, -x3) + dx2 * x;
    const CH: [(u64, u64); 4] = [
        (0xc014abbce625be53, 0x3cb05511c68477be),
        (0x400466bc6775aae2, 0xbc96dc0cbefae1da),
        (0xbfe32d2cce62bd86, 0x3c8066bd54973829),
        (0x3fb50783487ee781, 0x3c5832989f39a743),
    ];
    const CL: [u64; 3] = [0xbf7e3074fde861f0, 0x3f3e8f4344534da6, 0xbef6f9cd7b8cb9db];

    let fl0 = f_fmla(x2, f64::from_bits(CL[2]), f64::from_bits(CL[1]));

    let fl = x2 * f_fmla(x2, fl0, f64::from_bits(CL[0]));
    let mut f = poly_dd_4(Dekker::new(dx2, x2), CH, fl);
    f = Dekker::mult(Dekker::new(dx3, x3), f);
    const PI0: f64 = f64::from_bits(0x4009200000000000);
    const PI1: f64 = f64::from_bits(0x3f4fb54442d1846a);
    const PI2: f64 = f64::from_bits(0xbbed9cceba3f91f2);
    let mut y0 = PI0 * x;
    let mut b = y0.to_bits() & 0x000fffffffffffff;
    b = b.wrapping_add(85u64 << 51);
    y0 = (y0 + f64::from_bits(b)) - f64::from_bits(b);
    let y0l = dd_fmla(PI0, x, -y0);
    let y1 = PI1 * x;
    let y2 = f_fmla(PI2, x, dd_fmla(PI1, x, -y1));
    let mut y = Dekker::add(Dekker::new(y2, y1), Dekker::new(0., y0l));
    y = Dekker::add(y, f);
    let k0 = Dekker::from_full_exact_add(y.hi, y0);
    let mut k1 = Dekker::from_full_exact_add(k0.lo, y.lo);

    let mut t = k1.hi.to_bits();
    if (t & 0x000fffffffffffff) == 0 {
        let w = k1.lo.to_bits();
        if ((w ^ t) >> 63) != 0 {
            t = t.wrapping_sub(1);
        } else {
            t = t.wrapping_add(1);
        }
        k1.hi = f64::from_bits(t);
    }
    k0.hi + k1.hi
}

#[cold]
fn as_sinpi_refine(iq: i32, z: f64) -> f64 {
    let x = z * f64::from_bits(0x3c00000000000000);
    let x2 = x * x;
    let dx2 = dd_fmla(x, x, -x2);
    const SH: [(u64, u64); 3] = [
        (0x400921fb54442d18, 0x3ca1a62633145c06),
        (0xbe94abbce625be53, 0x3b305511cbc65743),
        (0x3d0466bc6775aae1, 0xb8d9c3c168d990a0),
    ];
    const CH: [(u64, u64); 2] = [
        (0xbe93bd3cc9be45de, 0xbb3692b71366cc04),
        (0x3d103c1f081b5ac4, 0xb9b32b33fda9113c),
    ];
    let mut sl = poly_dd_3(
        Dekker::new(dx2, x2),
        SH,
        f64::from_bits(0xbb632d2cc920dcb4) * x2,
    );
    sl = Dekker::mult_f64(sl, x * f64::from_bits(0x3f30000000000000));
    let cll0 = x2
        * f_fmla(
            f64::from_bits(0x39ce1f50604fa0ff),
            x2,
            f64::from_bits(0xbb755d3c7e3cbff9),
        );
    let mut cl = poly_dd_2(Dekker::new(dx2, x2), CH, cll0);
    cl = Dekker::mult(cl, Dekker::new(dx2, x2));
    let (sb, cb) = sincosn2(iq);
    let cs = Dekker::mult(cl, sb);
    let sc = Dekker::mult(sl, cb);
    let mut ts = Dekker::add(sc, cs);
    ts = Dekker::add(sb, ts);
    ts.to_f64()
}

const SIN_COEFFS: [u64; 3] = [0x3b5921fb54442d18, 0xb204abbce625be51, 0x289466bc6044ba16];
const COS_COEFFS: [u64; 2] = [0xb6b3bd3cc9be45db, 0x2d503c1f00186416];

/// Computes sin(PI*x)
///
/// Max ULP 0.5
#[inline]
pub fn f_sinpi(x: f64) -> f64 {
    let ix = x.to_bits();
    let ax = ix & 0x7fff_ffff_ffff_ffff;
    if ax == 0 {
        return x;
    }
    let e: i32 = (ax >> 52) as i32;
    let m0 = (ix & 0x000fffffffffffff) | (1u64 << 52);
    let sgn: i64 = (ix as i64) >> 63;
    let m = ((m0 as i64) ^ sgn).wrapping_sub(sgn);
    let mut s: i32 = 1063i32.wrapping_sub(e);
    if s < 0 {
        if e == 0x7ff {
            if (ix << 12) == 0 {
                return f64::NAN;
            }
            return x + x; // case x=NaN
        }
        s = -s - 1;
        if s > 10 {
            return f64::copysign(0.0, x);
        }
        let iq: u64 = (m as u64).wrapping_shl(s as u32);
        if (iq & 2047) == 0 {
            return f64::copysign(0.0, x);
        }
        let (t_sin, _) = sincosn(m.wrapping_shl(s as u32) as i32);
        return t_sin.to_f64();
    }

    if ax <= 0x3fa2000000000000u64 {
        // |x| <= 0.03515625
        let p = Dekker::new(
            f64::from_bits(0x3ca1a62633145c07),
            f64::from_bits(0x400921fb54442d18),
        );

        if x.abs() < f64::from_bits(0x3c90000000000000) {
            if x.abs() < f64::from_bits(0x0350000000000000) {
                let t = x * f64::from_bits(0x4690000000000000);
                let zh = p.hi * t;
                let zl = f_fmla(p.lo, t, dd_fmla(p.hi, t, -zh));
                let r = zh + zl;
                let rs = r * f64::from_bits(0x3950000000000000);
                let rt = rs * f64::from_bits(0x4690000000000000);
                return dd_fmla((zh - rt) + zl, f64::from_bits(0x3950000000000000), rs);
            }
            let zh = p.hi * x;
            let zl = dd_fmla(p.lo, x, f_fmla(p.hi, x, -zh));
            return zh + zl;
        }
        let zh = p.hi * x;
        let mut zl = f_fmla(p.lo, x, dd_fmla(p.hi, x, -zh));
        let x2 = x * x;
        let x3 = x2 * x;
        let x4 = x2 * x2;
        let eps = x * f_fmla(
            x2,
            f64::from_bits(0x3d00000000000000),
            f64::from_bits(0x3990000000000000),
        );
        const C: [u64; 4] = [
            0xc014abbce625be51,
            0x400466bc67754b46,
            0xbfe32d2cc12a51f4,
            0x3fb5060540058476,
        ];

        let zl0 = f_fmla(x2, f64::from_bits(C[1]), f64::from_bits(C[0]));
        let zl1 = f_fmla(x2, f64::from_bits(C[3]), f64::from_bits(C[2]));

        zl += x3 * f_fmla(x4, zl1, zl0);
        let lb = (zl - eps) + zh;
        let ub = (zl + eps) + zh;
        if lb == ub {
            return lb;
        }
        return as_sinpi_zero(x);
    }

    let si = e - 1011;
    if si >= 0 && (m0 << (si.wrapping_add(1))) == 0 {
        // x is integer or half-integer
        if (m0 << si) == 0 {
            return f64::copysign(0.0, x); // x is integer
        }
        let t = (m0.wrapping_shl((si - 1) as u32)) >> 63;
        // t = 0 if |x| = 1/2 mod 2, t = 1 if |x| = 3/2 mod 2
        return if t == 0 {
            f64::copysign(1.0, x)
        } else {
            -f64::copysign(1.0, x)
        };
    }

    let mut iq: u64 = ((m >> s) & 8191) as u64;
    iq = (iq.wrapping_add(1)) >> 1;
    let k: i64 = (m as u64).wrapping_shl(e.wrapping_sub(1000) as u32) as i64;
    let z = k as f64;
    let z2 = z * z;

    let fs0 = f_fmla(
        z2,
        f64::from_bits(SIN_COEFFS[2]),
        f64::from_bits(SIN_COEFFS[1]),
    );

    let fc = f_fmla(
        z2,
        f64::from_bits(COS_COEFFS[1]),
        f64::from_bits(COS_COEFFS[0]),
    );

    let fs = f_fmla(z2, fs0, f64::from_bits(SIN_COEFFS[0]));
    let (t_sin, t_cos) = sincosn(iq as i32);
    const ERR: f64 = f64::from_bits(0x3c244a9a66cdb0d6);

    let r0 = f_fmla(t_sin.hi, z2 * fc, t_sin.lo);

    let r = f_fmla(t_cos.hi, z * fs, r0);
    let lb = (r - ERR) + t_sin.hi;
    let ub = (r + ERR) + t_sin.hi;
    if lb == ub {
        return lb;
    }
    as_sinpi_refine(iq as i32, z)
}

/// Computes cos(PI*x)
///
/// Max found ULP 0.5
#[inline]
pub fn f_cospi(x: f64) -> f64 {
    let ix = x.to_bits();
    let ax = ix & 0x7fff_ffff_ffff_ffff;
    if ax == 0 {
        return 1.0;
    }
    let e: i32 = (ax >> 52) as i32;
    // e is the unbiased exponent, we have 2^(e-1023) <= |x| < 2^(e-1022)
    let m: i64 = ((ix & 0x000fffffffffffff) | (1u64 << 52)) as i64;
    let mut s = 1063i32.wrapping_sub(e); // 2^(40-s) <= |x| < 2^(41-s)
    if s < 0 {
        // |x| >= 2^41
        if e == 0x7ff {
            // NaN or Inf
            if ix.wrapping_shl(12) == 0 {
                return f64::NAN;
            }
            return x + x; // NaN
        }
        s = -s - 1; // now 2^(41+s) <= |x| < 2^(42+s)
        if s > 11 {
            return 1.0;
        } // |x| >= 2^53
        let iq: u64 = (m as u64).wrapping_shl(s as u32).wrapping_add(1024);
        if (iq & 2047) == 0 {
            return 0.0;
        }
        let (sin, _) = sincosn(iq as i32);
        return sin.hi + sin.lo;
    }
    if ax <= 0x3f30000000000000u64 {
        // |x| <= 2^-12
        if ax <= 0x3e2ccf6429be6621u64 {
            return 1.0 - f64::from_bits(0x3c80000000000000);
        }
        let x2 = x * x;
        let x4 = x2 * x2;
        let eps = x2 * f64::from_bits(0x3cfa000000000000);

        const C: [u64; 4] = [
            0xc013bd3cc9be45dc,
            0x40103c1f081b0833,
            0xbff55d3c6fc9af15,
            0x3fce1d3ff2ae3f9a,
        ];

        let p0 = f_fmla(x2, f64::from_bits(C[3]), f64::from_bits(C[2]));
        let p1 = f_fmla(x2, f64::from_bits(C[1]), f64::from_bits(C[0]));

        let p = x2 * f_fmla(x4, p0, p1);
        let lb = (p - eps) + 1.;
        let ub = (p + eps) + 1.;
        if lb == ub {
            return lb;
        }
        return as_cospi_zero(x);
    }

    let si: i32 = e.wrapping_sub(1011);
    if si >= 0 && ((m as u64).wrapping_shl(si as u32) ^ 0x8000000000000000u64) == 0 {
        return 0.0;
    }

    let mut iq = ((m >> s).wrapping_add(2048)) & 8191;
    iq = (iq.wrapping_add(1)) >> 1;
    let k: i64 = (m as u64).wrapping_shl((e - 1000) as u32) as i64;
    let z = k as f64;
    let z2 = z * z;

    let fs0 = f_fmla(
        z2,
        f64::from_bits(SIN_COEFFS[2]),
        f64::from_bits(SIN_COEFFS[1]),
    );

    let fs = f_fmla(z2, fs0, f64::from_bits(SIN_COEFFS[0]));
    let fc = f_fmla(
        z2,
        f64::from_bits(COS_COEFFS[1]),
        f64::from_bits(COS_COEFFS[0]),
    );
    let (t_sin, t_cos) = sincosn(iq as i32);
    let err = z * f64::from_bits(0x3840000000000000);
    let r = f_fmla(t_cos.hi, z * fs, f_fmla(t_sin.hi, z2 * fc, t_sin.lo));
    let lb = (r - err) + t_sin.hi;
    let ub = (r + err) + t_sin.hi;
    if lb == ub {
        return lb;
    }
    as_sinpi_refine(iq as i32, z)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_sinpi() {
        assert_eq!(f_sinpi(0.0031909299901270445), 0.0100244343161398578);
        assert_eq!(f_sinpi(0.11909245901270445), 0.36547215190661003);
        assert_eq!(f_sinpi(0.99909245901270445), 0.0028511202357662186);
        assert!(f_sinpi(f64::INFINITY).is_nan());
    }

    #[test]
    fn test_cospi() {
        assert_eq!(0.9999497540959953, f_cospi(0.0031909299901270445));
        assert_eq!(0.9308216542079669, f_cospi(0.11909299901270445));
        assert_eq!(-0.1536194873288318, f_cospi(0.54909299901270445));
        assert!(f_cospi(f64::INFINITY).is_nan());
    }
}
