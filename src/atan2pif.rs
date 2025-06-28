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
use crate::atan2f::poly_dekker_generic;
use crate::common::f_fmla;
use crate::dekker::Dekker;

static ATANPIF_DD: [(u64, u64); 32] = [
    (0x3fd45f306dc9c883, 0xbc76b01ec5513324),
    (0xbfbb2995e7b7b604, 0x3c5e402b0c13eedc),
    (0x3fb04c26be3b06cf, 0xbc3571d178a53ef0),
    (0xbfa7483758e69c03, 0x3c0819a6ed7aaf38),
    (0x3fa21bb9452523ff, 0xbc3234d866fb9807),
    (0xbf9da1bace3cc54e, 0xbbfc84f6ada49294),
    (0x3f9912b1c23345dd, 0xbc3534890fbc1650),
    (0xbf95bade52f5f52a, 0x3c3f783bafc832f6),
    (0x3f932c69d084c5c0, 0x3c3042d155953025),
    (0xbf9127bcfb3e8c7d, 0xbc385aae199a7b6b),
    (0x3f8f0af43b11a731, 0x3c28f03563566630),
    (0xbf8c57e86801029e, 0x3c2dcdf3e3b38eb4),
    (0x3f8a136408617ea1, 0x3c0a71affb36c6c4),
    (0xbf8824ac7814ba37, 0x3c28928b295c0898),
    (0x3f86794e32ea5471, 0x3c20b4334fb41e63),
    (0xbf8501d57f643d97, 0x3c2516785bf1376e),
    (0x3f83adf02ff2400a, 0xbc1b0e30bb8c8076),
    (0xbf8267702f94faa0, 0xbc17a4d3a1850cc6),
    (0x3f810dce97099686, 0x3c2fcc208eee2571),
    (0xbf7eee49cdad8002, 0xbbf9109b3f1bab82),
    (0x3f7af93bc191a929, 0x3c1069fd3b47d7b0),
    (0xbf76240751b54675, 0xbc172dc8cfd03b6f),
    (0x3f70b61e84080884, 0x3c0825824c80941b),
    (0xbf66a72a8a74e3a5, 0x3c08786a82fd117e),
    (0x3f5aede3217d939d, 0xbbb93b626982e1fe),
    (0xbf4b66568f09ebee, 0xbbd704a39121d0a5),
    (0x3f373af3977fa973, 0xbbbaa050e2244ea3),
    (0xbf1fc69d85ed28c9, 0x3bb867f17b764ca0),
    (0x3f00c883a9270162, 0xbb96842833896dd9),
    (0xbed9a0b27b6dfe15, 0x3b6427fc2f4e1327),
    (0x3ea91e15e7ab5bdc, 0xbb2730dbc6279d0d),
    (0xbe67b1119c1ff867, 0x3b0145f9980759c4),
];

static OFF: [f32; 8] = [0.0, 0.5, 1.0, 0.5, -0.0, -0.5, -1.0, -0.5];
static OFF_F64: [f64; 8] = [0.0, 0.5, 1.0, 0.5, -0.0, -0.5, -1.0, -0.5];
static SGNF: [f32; 2] = [1., -1.];
static SGN: [f64; 2] = [1., -1.];

/// Computes atan(x/y * PI)
///
/// Max found ULP 0.5
#[inline]
pub fn f_atan2pif(y: f32, x: f32) -> f32 {
    let tx = x.to_bits();
    let ty: u32 = y.to_bits();
    let ux: u32 = tx;
    let uy: u32 = ty;
    let ax: u32 = ux & 0x7fff_ffff;
    let ay = uy & 0x7fff_ffff;
    if ay >= (0xff << 23) || ax >= (0xff << 23) {
        if ay > (0xff << 23) {
            return x + y;
        } // nan
        if ax > (0xff << 23) {
            return x + y;
        } // nan
        let yinf = ay == (0xff << 23);
        let xinf = ax == (0xff << 23);
        if yinf & xinf {
            return if (ux >> 31) != 0 {
                0.75 * SGNF[(uy >> 31) as usize]
            } else {
                0.25 * SGNF[(uy >> 31) as usize]
            };
        }
        if xinf {
            return if (ux >> 31) != 0 {
                SGNF[(uy >> 31) as usize]
            } else {
                0.0 * SGNF[(uy >> 31) as usize]
            };
        }
        if yinf {
            return 0.5 * SGNF[(uy >> 31) as usize];
        }
    }
    if ay == 0 {
        if (ay | ax) == 0 {
            let i: u32 = (uy >> 31) * 4 + (ux >> 31) * 2;
            return OFF[i as usize];
        }
        if (ux >> 31) == 0 {
            return 0.0 * SGNF[(uy >> 31) as usize];
        }
    }
    if ax == ay {
        static S: [f32; 4] = [0.25, 0.75, -0.25, -0.75];
        let i = (uy >> 31) * 2 + (ux >> 31);
        return S[i as usize];
    }
    let gt: usize = if ay > ax { 1 } else { 0 };
    let i: u32 = (uy >> 31) * 4 + (ux >> 31) * 2 + gt as u32;

    let zx = x as f64;
    let zy = y as f64;
    static M: [f64; 2] = [0., 1.];

    let mut z = f_fmla(M[gt], zx, M[1 - gt] * zy) / f_fmla(M[gt], zy, M[1 - gt] * zx);

    const CN: [u64; 7] = [
        0x3fd45f306dc9c883,
        0x3fe988d83a142ada,
        0x3fe747bebf492057,
        0x3fd2cc5645094ff3,
        0x3faa0521c711ab66,
        0x3f6881b8058b9a0d,
        0x3efb16ff514a0af0,
    ];

    let mut r = f64::from_bits(CN[0]);
    let z2 = z * z;
    z *= SGN[gt];
    // avoid spurious underflow in the polynomial evaluation excluding tiny arguments
    if z2 > f64::from_bits(0x3c90000000000000) {
        let z4 = z2 * z2;
        let z8 = z4 * z4;
        let mut cn0 = f_fmla(z2, f64::from_bits(CN[1]), r);
        let cn2 = f_fmla(z2, f64::from_bits(CN[3]), f64::from_bits(CN[2]));
        let mut cn4 = f_fmla(z2, f64::from_bits(CN[5]), f64::from_bits(CN[4]));
        let cn6 = f64::from_bits(CN[6]);
        cn0 += z4 * cn2;
        cn4 += z4 * cn6;
        cn0 += z8 * cn4;

        const CD: [u64; 7] = [
            0x3ff0000000000000,
            0x4006b8b143a3f6da,
            0x4008421201d18ed5,
            0x3ff8221d086914eb,
            0x3fd670657e3a07ba,
            0x3fa0f4951fd1e72d,
            0x3f4b3874b8798286,
        ];

        let mut cd0 = f_fmla(z2, f64::from_bits(CD[1]), f64::from_bits(CD[0]));
        let cd2 = f_fmla(z2, f64::from_bits(CD[3]), f64::from_bits(CD[2]));
        let mut cd4 = f_fmla(z2, f64::from_bits(CD[5]), f64::from_bits(CD[4]));
        let cd6 = f64::from_bits(CD[6]);
        cd0 += z4 * cd2;
        cd4 += z4 * cd6;
        cd0 += z8 * cd4;

        r = cn0 / cd0;
    }
    r = f_fmla(z, r, OFF[i as usize] as f64);
    let res = r.to_bits();
    if (res.wrapping_shl(1)) > 0x6d40000000000000 && ((res.wrapping_add(8)) & 0xfffffff) <= 16 {
        // |res| > 0x1p-149
        if ax == ay {
            static OFF2: [f64; 4] = [0.25, 0.75, -0.25, -0.75];
            r = OFF2[((uy >> 31) * 2 + (ux >> 31)) as usize];
        } else {
            let (zh, zl);
            if gt == 0 {
                zh = zy / zx;
                zl = f_fmla(zh, -zx, zy) / zx;
            } else {
                zh = zx / zy;
                zl = f_fmla(zh, -zy, zx) / zy;
            }
            let mut zd = Dekker::new(zl, zh);
            let zd2 = Dekker::quick_mult(zd, zd);
            let mut p = poly_dekker_generic(zd2, ATANPIF_DD);
            zd.hi *= SGN[gt];
            zd.lo *= SGN[gt];
            p = Dekker::quick_mult(zd, p);
            let sh = p.hi + OFF_F64[i as usize];
            let sl = ((OFF_F64[i as usize] - sh) + p.hi) + p.lo;
            let rf = sh as f32;
            let th = rf as f64;
            let dh = sh - th;
            let tm = dh + sl;
            r = th + tm;
            let d = (r - th).to_bits();
            if d.wrapping_shl(12) == 0 {
                let ad = f64::from_bits(d & 0x7fff_ffff_ffff_ffff);
                let am = tm.abs();
                if ad > am {
                    r -= f64::from_bits(d) * f64::from_bits(0x3f50000000000000);
                }
                if ad < am {
                    r += f64::from_bits(d) * f64::from_bits(0x3f50000000000000);
                }
            }
        }
    }
    r as f32
}
