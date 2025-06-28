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
use crate::common::f_fmla;
use crate::dekker::Dekker;

static TANPI_REDUCE: [(u64, u64); 32] = [
    (0x0000000000000000, 0x0000000000000000),
    (0x3bfada13ceebab9d, 0x3fa927278a3b1162),
    (0x3c587d054f689d7a, 0x3fb936bb8c5b2da2),
    (0x3c52cfb5a746f62c, 0x3fc2fcac73a60640),
    (0x3c2ef5d367441946, 0x3fc975f5e0553158),
    (0x3c65a6d6c3c8b6a7, 0x3fd007fa758626ae),
    (0xbc6d704d1bfdb6e8, 0x3fd36a08355c63dc),
    (0x3c678e945dc3913c, 0x3fd6e649f7d78649),
    (0x3c708b2fb1366ea9, 0x3fda827999fcef32),
    (0x3c741522f15e53c5, 0x3fde450e0d273e7a),
    (0xbc8e564bcd1e635d, 0x3fe11ab7190834ec),
    (0xbc810b4421e6a4f8, 0x3fe32e1889047ffd),
    (0x3c87a8c52172b675, 0x3fe561b82ab7f990),
    (0xbc8aa7538e44e996, 0x3fe7bb99ed2990cf),
    (0xbc7a6db80fe796fe, 0x3fea43002ae42850),
    (0x3c78dcad85e60fbe, 0x3fed00cbc7384d2e),
    (0x0000000000000000, 0x3ff0000000000000),
    (0x3c9946cc0b66979f, 0x3ff1a73d55278c4b),
    (0xbc787e5ad9735569, 0x3ff37efd8d87607e),
    (0x3c86a085e3bc3af6, 0x3ff592d11142fa55),
    (0x3c9419fa6954928f, 0x3ff7f218e25a7461),
    (0xbc7b6fb77846d763, 0x3ffab1c35d8a74ea),
    (0x3c80fb3e75c7098e, 0x3ffdef13b73c1406),
    (0x3c97ce6cb463c972, 0x4000ea21d716fbf7),
    (0x3ca21165f626cdd5, 0x4003504f333f9de6),
    (0x3c7aca947bfb1dcc, 0x40065bc6cc825147),
    (0xbc9b7a14d0d691be, 0x400a5f59e90600dd),
    (0xbc889fcd637fbf3b, 0x400ff01305ecd8dc),
    (0x3c810706fed37f0e, 0x40141bfee2424771),
    (0xbcaae397239c5a0d, 0x401af73f4ca3310f),
    (0x3cc351daea79411d, 0x40244e6c595afdcc),
    (0xbc7b6e6b4de0cd24, 0x40345affed201b55),
];

#[inline]
pub(crate) fn poly_dd_5(x: Dekker, poly: [(u64, u64); 5], l: f64) -> Dekker {
    let zch = poly[4];
    let ach = f64::from_bits(zch.0) + l;
    let acl = (f64::from_bits(zch.0) - ach) + l + f64::from_bits(zch.1);
    let mut ch = Dekker::new(acl, ach);

    let zch = poly[3];
    ch = Dekker::mult(ch, x);
    let th = ch.hi + f64::from_bits(zch.0);
    let tl = (f64::from_bits(zch.0) - th) + ch.hi;
    ch.hi = th;
    ch.lo += tl + f64::from_bits(zch.1);

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

/// Computes tan(PI*x)
///
/// Max found ULP 5.0001
#[inline]
pub fn f_tanpi(x: f64) -> f64 {
    let ix = x.to_bits();
    let ax: u64 = ix & 0x7fff_ffff_ffff_ffff;
    let mut t: Dekker;
    let res: f64;
    if ax >= (0x3f3u64 << 52) {
        // |x| >= 0x1p-12
        if ax >= (0x42du64 << 52) {
            // |x| >= 0x1p+46
            if ax >= (0x7ffu64 << 52) {
                // NaN, Inf
                if ax > (0x7ffu64 << 52) {
                    return x + x;
                } // NaN
                return f64::NAN; // x=Inf
            }
            let e: i32 = (ax >> 52) as i32;
            let s = e.wrapping_sub(1069);
            if s > 6 {
                return f64::copysign(0., x);
            }
            let m: i64 = ax as i64;
            let sgn = (ix as i64) >> 63;
            let iq = ((((m ^ sgn) - sgn) as u64).wrapping_shl(s as u32)) & 127;
            if (iq & 31) == 0 {
                let jq: i64 = (iq >> 5) as i64;
                if (jq & 1) != 0 {
                    return if (jq & 2) != 0 {
                        f64::NEG_INFINITY
                    } else {
                        f64::INFINITY
                    };
                } else {
                    return if (jq ^ sgn) & 2 != 0 { -0.0 } else { 0.0 };
                }
            } else {
                let nh;
                let nl;
                if (iq & 32) != 0 {
                    let rr = TANPI_REDUCE[(32 - (iq & 31)) as usize];
                    nl = -f64::from_bits(rr.0);
                    nh = -f64::from_bits(rr.1);
                } else {
                    let rr = TANPI_REDUCE[(iq & 31) as usize];
                    nl = f64::from_bits(rr.0);
                    nh = f64::from_bits(rr.1);
                }
                return nh + nl;
            }
        }
        // now 0.000244140625 <= |x| < 70368744177664, we have 1011 <= e <= 1068
        let e: i32 = (ax >> 52) as i32;
        let s = 1068i32.wrapping_sub(e);
        let s1 = e.wrapping_sub(1011);
        let m = ((ax & 0x000fffffffffffff) | (1u64 << 52)) as i64;
        let mut ms = ((m as u64).wrapping_shl(s1 as u32) as i64) >> 63;
        let sgn = (ix as i64) >> 63;
        // m is the significand, 2^52 <= m < 2^53
        // 0 <= s1 <= 57 is a biased exponent, with s1=0 for 2^-12 <= |x| < 2^-11
        // 0 <= s <= 57 is another biased exponent, with s=0 for 2^45 <= |x| < 2^46
        // ms is the bit of weight 1/2 in x
        let mut iq: u64 = (((m ^ ms) >> s) & 63) as u64;
        iq = (iq.wrapping_add(1)) >> 1;
        ms ^= sgn;
        let sm: i64 = (m ^ sgn) - sgn; // sm = sign(x)*m
        let k: i64 = (sm as u64).wrapping_shl(e.wrapping_sub(1005) as u32) as i64; // 6 <= e-1005 <= 63
        // k contains the bits of m of weight <= 2^-7
        let mut z = k as f64;
        if ((k as u64).wrapping_shl(1)) == 0 {
            // x mod 2^-8 = 0
            if k == 0 {
                // x mod 2^-7 = 0
                if (iq & 31) == 0 {
                    let jq: i64 = sm >> (s.wrapping_add(6));
                    return if (jq & 1) != 0 {
                        if (jq & 2) != 0 {
                            f64::NEG_INFINITY
                        } else {
                            f64::INFINITY
                        }
                    } else if ((jq ^ sgn) & 2) != 0 {
                        -0.0
                    } else {
                        0.0
                    };
                }
                // avoid spurious inexact exception for x=1/4 mod 1/2
                let kq: u64 = ((m as u64).wrapping_shl(s1 as u32)) >> 58;
                if kq == 0x10 {
                    // |x| = 1/4 mod 1
                    return f64::copysign(1., x);
                }
                if kq == 0x30 {
                    // |x| = 3/4 mod 1
                    return -f64::copysign(1., x);
                }
            }
            z *= f64::copysign(1., x);
        }
        let mut z2 = z * z;
        let z4 = z2 * z2;
        let z3 = z * z2;
        const C: [u64; 4] = [
            0x3304abbce625be51,
            0x2a6466bc6776a9b1,
            0x21c45fff6eb26045,
            0x1924627663861052,
        ];

        let f0 = f_fmla(z2, f64::from_bits(C[3]), f64::from_bits(C[2]));
        let f1 = f_fmla(z2, f64::from_bits(C[1]), f64::from_bits(C[0]));

        let f = z3 * f_fmla(z4, f0, f1);
        let mut eps = z3 * f64::from_bits(0x2ff0000000000000)
            + f64::copysign(f64::from_bits(0x3970000000000000), z);

        const P: Dekker = Dekker::new(
            f64::from_bits(0x3841a62633145c07),
            f64::from_bits(0x3ba921fb54442d18),
        );

        t = Dekker::mult_f64(P, z);
        t = Dekker::add_f64(t, f);

        if iq == 32 {
            let ith = -1.0 / t.hi;
            t.lo = (f_fmla(ith, t.hi, 1.) + t.lo * ith) * ith;
            t.hi = ith;
        } else {
            let mut n = Dekker::from_bit_pair(TANPI_REDUCE[iq as usize]);
            static S2: [f64; 2] = [-1., 1.];
            n.hi *= S2[(ms + 1) as usize];
            n.lo *= S2[(ms + 1) as usize];
            let mut m = Dekker::mult(t, n);
            let z0 = Dekker::from_exact_sub(1.0, m.hi);
            m.hi = z0.hi;
            m.lo = z0.lo - m.lo;

            let z1 = Dekker::from_exact_add(n.hi, t.hi);
            n.hi = z1.hi;
            n.lo += z1.lo + t.lo;

            let imh = 1.0 / m.hi;
            t.hi = n.hi * imh;
            t.lo = f_fmla(n.hi, imh, -t.hi)
                + (n.lo + n.hi * (f_fmla(-m.hi, imh, 1.) - m.lo * imh)) * imh;
        }
        eps += eps * (t.hi * t.hi);
        let lb = t.hi + (t.lo - eps);
        let ub = t.hi + (t.lo + eps);
        if lb == ub {
            return lb;
        }
        z *= f64::from_bits(0x3c00000000000000);

        const CH: [(u64, u64); 5] = [
            (0x3f9921fb54442d18, 0x3c31a62633145c07),
            (0x3ed4abbce625be53, 0xbb705511c6842515),
            (0x3e1466bc6775aae2, 0xbaa6dc0d93fb2ece),
            (0x3d545fff9b48e95e, 0x39db226b250d2cc3),
            (0x3c945f472e3af011, 0x39190612a0755449),
        ];
        const CL: [u64; 3] = [0x3bd45f32f25dab7f, 0x3b145f3030c82af4, 0x3a5464490600a978];
        z2 = z * z;
        let dz2 = f_fmla(z, z, -z2);

        let tlo0 = f_fmla(z2, f64::from_bits(CL[2]), f64::from_bits(CL[1]));

        t.lo = z2 * f_fmla(z2, tlo0, f64::from_bits(CL[0]));
        t = poly_dd_5(Dekker::new(dz2, z2), CH, t.lo);
        t = Dekker::mult_f64(t, z);
        if iq == 32 {
            let ith = -1.0 / t.hi;
            t.lo = (f_fmla(ith, t.hi, 1.) + t.lo * ith) * ith;
            t.hi = ith;
        } else {
            let mut n = Dekker::from_bit_pair(TANPI_REDUCE[iq as usize]);
            static S2: [f64; 2] = [-1., 1.];
            n.hi *= S2[(ms + 1) as usize];
            n.lo *= S2[(ms + 1) as usize];
            let mut m = Dekker::mult(t, n);
            let z0 = Dekker::from_exact_sub(1.0, m.hi);
            m.hi = z0.hi;
            m.lo = z0.lo - m.lo;

            let z1 = Dekker::from_exact_add(n.hi, t.hi);
            n.hi = z1.hi;
            n.lo += z1.lo + t.lo;

            let imh = 1.0 / m.hi;
            t.hi = n.hi * imh;
            t.lo = f_fmla(n.hi, imh, -t.hi)
                + (n.lo + n.hi * (f_fmla(-m.hi, imh, 1.) - m.lo * imh)) * imh;
        }
        t = Dekker::from_exact_add(t.hi, t.lo);
        res = t.hi;
    } else {
        // |x| < 0.000244140625
        if ax == 0 {
            return x;
        }
        const P2: Dekker = Dekker::new(
            f64::from_bits(0x3ca1a62633145c07),
            f64::from_bits(0x400921fb54442d18),
        );
        if ax < (0x3cau64 << 52) {
            // |x| < 0x1p-53
            if ax < (0x36u64 << 52) {
                // |x| < 0x1p-969
                let e: i32 = (ax >> 52) as i32;
                let sc = (2045i64 - e as i64).wrapping_shl(52) as u64;
                let isc = 1i64.wrapping_add(e as i64).wrapping_shl(52) as u64;
                let z = x * f64::from_bits(sc);
                t = Dekker::mult_f64(P2, z);
                res = t.hi * f64::from_bits(isc);
                if res.abs() < f64::from_bits(0x0010000000000000) {
                    // we force underflow since the code below might be exact
                    let o = f64::copysign(f64::from_bits(0x0010000000000000), x);
                    let mut v0b = (o + res) * f64::from_bits(sc);
                    let v0h = res * f64::from_bits(sc);
                    t.lo += t.hi - v0h;
                    v0b += t.lo;
                    return v0b * f64::from_bits(isc) - o;
                }
            } else {
                // 0x1p-969 <= |x| < 0x1p-53
                t = Dekker::mult_f64(P2, x);
                res = t.hi;
            }
        } else {
            const C2: [u64; 3] = [0x4024abbce625be53, 0x404466bc6775aa2a, 0x40646000158496c2];
            let x2 = x * x;
            let x3 = x * x2;

            let f0 = f_fmla(x2, f64::from_bits(C2[2]), f64::from_bits(C2[1]));

            let f = x3 * f_fmla(x2, f0, f64::from_bits(C2[0]));
            let px = Dekker::mult_f64(P2, x);
            t = Dekker::add_f64(px, f);
            let eps =
                x * (x2 * f64::from_bits(0x3d01000000000000) + f64::from_bits(0x39a0000000000000));
            let lb = t.hi + (t.lo - eps);
            let ub = t.hi + (t.lo + eps);
            if lb == ub {
                return lb;
            }

            const CH: [(u64, u64); 3] = [
                (0x4024abbce625be53, 0xbcc05511c68476a8),
                (0x404466bc6775aae2, 0xbcd6dc0cbddc0e69),
                (0x40645fff9b48e95e, 0x3cea5047910ae0ef),
            ];

            const CL: [u64; 2] = [0x40845f472e3aed7d, 0x40a45f33be0e9598];

            let dx2 = f_fmla(x, x, -x2);
            let dx3 = f_fmla(x2, x, -x3) + dx2 * x;
            t.lo = x2 * f_fmla(x2, f64::from_bits(CL[1]), f64::from_bits(CL[0]));
            t = poly_dd_3(Dekker::new(dx2, x2), CH, t.lo);
            t = Dekker::mult(Dekker::new(dx3, x3), t);
            let z0 = Dekker::from_exact_add(px.hi, t.hi);
            t.hi = z0.hi;
            t.lo = t.lo + px.lo + z0.lo;
            t = Dekker::from_exact_add(t.hi, t.lo);
            res = t.hi;
        }
    }
    res
}

#[cfg(test)]
mod tests {

    use super::*;

    #[test]
    fn test_tanpi() {
        assert_eq!(-2867080569611329.5, f_tanpi(0.5000000000000001));
        #[cfg(any(
            all(
                any(target_arch = "x86", target_arch = "x86_64"),
                target_feature = "fma"
            ),
            all(target_arch = "aarch64", target_feature = "neon")
        ))]
        {
            assert_eq!(0.06704753721009375, f_tanpi(0.02131));
        }
        #[cfg(not(any(
            all(
                any(target_arch = "x86", target_arch = "x86_64"),
                target_feature = "fma"
            ),
            all(target_arch = "aarch64", target_feature = "neon")
        )))]
        {
            assert_eq!(0.06704753721009377, f_tanpi(0.02131));
        }
    }
}
