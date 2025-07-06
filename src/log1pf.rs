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
use crate::common::{f_fmla, f_fmlaf};
use std::hint::black_box;

static X0: [u64; 32] = [
    0x3fef81f820000000,
    0x3fee9131ac000000,
    0x3fedae6077000000,
    0x3fecd85689000000,
    0x3fec0e0704000000,
    0x3feb4e81b5000000,
    0x3fea98ef60000000,
    0x3fe9ec8e95000000,
    0x3fe948b0fd000000,
    0x3fe8acb90f000000,
    0x3fe8181818000000,
    0x3fe78a4c81000000,
    0x3fe702e05c000000,
    0x3fe6816817000000,
    0x3fe6058160000000,
    0x3fe58ed231000000,
    0x3fe51d07eb000000,
    0x3fe4afd6a0000000,
    0x3fe446f865000000,
    0x3fe3e22cbd000000,
    0x3fe3813814000000,
    0x3fe323e34a000000,
    0x3fe2c9fb4e000000,
    0x3fe27350b9000000,
    0x3fe21fb781000000,
    0x3fe1cf06ae000000,
    0x3fe1811812000000,
    0x3fe135c811000000,
    0x3fe0ecf56c000000,
    0x3fe0a6810a000000,
    0x3fe0624dd3000000,
    0x3fe0204081000000,
];

static LIXB: [u64; 32] = [
    0x3f8fc0a8909b4218,
    0x3fa77458f51aac89,
    0x3fb341d793afb997,
    0x3fba926d3a5ebd2a,
    0x3fc0d77e7a8a823d,
    0x3fc44d2b6c557102,
    0x3fc7ab89040accec,
    0x3fcaf3c94ecab3d6,
    0x3fce27076d54e6c9,
    0x3fd0a324e3888ad5,
    0x3fd22941fc0c7357,
    0x3fd3a64c56ae3fdb,
    0x3fd51aad874af21f,
    0x3fd686c81d300ea0,
    0x3fd7eaf83c7fa9b5,
    0x3fd947941aa610ec,
    0x3fda9cec9a3f023b,
    0x3fdbeb4d9ea4156e,
    0x3fdd32fe7f35e5c7,
    0x3fde7442617b817a,
    0x3fdfaf588dd5ed10,
    0x3fe0723e5c635c39,
    0x3fe109f39d53c990,
    0x3fe19ee6b38a4668,
    0x3fe23130d7f93c3b,
    0x3fe2c0e9ec9b0b85,
    0x3fe34e289cb35ecc,
    0x3fe3d9026ad3d3f3,
    0x3fe4618bc1eadbbb,
    0x3fe4e7d8127dd8a9,
    0x3fe56bf9d5967092,
    0x3fe5ee02a926936e,
];

static LIX: [u64; 32] = [
    0x3f8fc0a890fc03e4,
    0x3fa77458f532dcfc,
    0x3fb341d793bbd1d1,
    0x3fba926d3a6ad563,
    0x3fc0d77e7a908e59,
    0x3fc44d2b6c5b7d1e,
    0x3fc7ab890410d909,
    0x3fcaf3c94ed0bff3,
    0x3fce27076d5af2e6,
    0x3fd0a324e38b90e3,
    0x3fd22941fc0f7966,
    0x3fd3a64c56b145ea,
    0x3fd51aad874df82d,
    0x3fd686c81d3314af,
    0x3fd7eaf83c82afc3,
    0x3fd947941aa916fb,
    0x3fda9cec9a42084a,
    0x3fdbeb4d9ea71b7c,
    0x3fdd32fe7f38ebd5,
    0x3fde7442617e8788,
    0x3fdfaf588dd8f31f,
    0x3fe0723e5c64df40,
    0x3fe109f39d554c97,
    0x3fe19ee6b38bc96f,
    0x3fe23130d7fabf43,
    0x3fe2c0e9ec9c8e8c,
    0x3fe34e289cb4e1d3,
    0x3fe3d9026ad556fb,
    0x3fe4618bc1ec5ec2,
    0x3fe4e7d8127f5bb1,
    0x3fe56bf9d597f399,
    0x3fe5ee02a9281675,
];

#[cold]
pub(crate) fn special_logf(x: f32) -> f32 {
    let t = x.to_bits();
    if t == 0xbf800000u32 {
        // +0.0
        return f32::NEG_INFINITY; // to raise FE_DIVBYZERO
    }
    if t == 0x7f800000u32 {
        return x;
    } // +inf
    let ax: u32 = t.wrapping_shl(1);
    if ax > 0xff000000u32 {
        return x + x;
    } // nan
    f32::NAN // to raise FE_INVALID
}

const B: [u64; 8] = [
    0x3ff0000000000000,
    0xbfe0000000000000,
    0x3fd5555555556f6b,
    0xbfd00000000029b9,
    0x3fc9999988d176e4,
    0xbfc55555418889a7,
    0x3fc24adeca50e2bc,
    0xbfc001ba33bf57cf,
];

#[cold]
fn log1pf_accurate(x: f32, z: f64, e: i32, j: usize) -> f32 {
    let z2 = z * z;
    let z4 = z2 * z2;

    let f0 = f_fmla(z, f64::from_bits(B[7]), f64::from_bits(B[6]));
    let f1 = f_fmla(z, f64::from_bits(B[5]), f64::from_bits(B[4]));
    let f2 = f_fmla(z, f64::from_bits(B[3]), f64::from_bits(B[2]));
    let f3 = f_fmla(z, f64::from_bits(B[1]), f64::from_bits(B[0]));

    let zf0 = f_fmla(z2, f0, f1);
    let zf1 = f_fmla(z2, f2, f3);
    let zf2 = f_fmla(z4, zf0, zf1);

    let f = z * zf2;
    const LN2L: f64 = f64::from_bits(0x3eb7f7d1cf79abca);
    const LN2H: f64 = f64::from_bits(0x3fe62e4000000000);
    let lh = LN2H * e as f64;
    let ll = LN2L * e as f64;
    let rl = f + ll + f64::from_bits(LIX[j]);
    let mut tr = (rl + lh).to_bits();
    if tr & 0xfffffffu64 == 0 {
        let x_bits = x.to_bits();
        if x_bits == 0xbc923d58u32 {
            return black_box(f32::from_bits(0xbc938f87)) - black_box(f32::from_bits(0x30000000));
        }
        if x_bits == 0xbd1d20afu32 {
            return black_box(f32::from_bits(0xbd203889)) + black_box(f32::from_bits(0x30800000));
        }
        if x_bits == 0x3efd81adu32 {
            return black_box(f32::from_bits(0x3ecdeee1)) + black_box(f32::from_bits(0x32000000));
        }
        tr = (f64::from_bits(tr) + 64. * (rl + (lh - f64::from_bits(tr)))).to_bits();
    } else if rl + (lh - f64::from_bits(tr)) == 0.0 {
        let x_bits = x.to_bits();
        if x_bits == 0x3ddbfec3u32 {
            return black_box(f32::from_bits(0x3dd0f671)) + black_box(f32::from_bits(0x31000000));
        }
        if x_bits == 0xbd1d20afu32 {
            return black_box(f32::from_bits(0xbd203889)) + black_box(f32::from_bits(0x30800000));
        }
        if x_bits == 0x3ca1e3f1u32 {
            return black_box(f32::from_bits(0x3ca04fc0)) + black_box(f32::from_bits(0x30000000));
        }
    }
    f64::from_bits(tr) as f32
}

#[cold]
fn log1pf_small(x: f32, z: f64, ax: u32) -> f32 {
    if ax < 0x33000000u32 {
        // |x| < 2.9802322387695312e-08
        if ax == 0 {
            return x;
        }
        let res = f_fmlaf(x, -x, x);
        return res;
    }
    let z2 = z * z;
    let z4 = z2 * z2;

    let f0 = f_fmla(z, f64::from_bits(B[6]), f64::from_bits(B[5]));
    let f1 = f_fmla(z, f64::from_bits(B[4]), f64::from_bits(B[3]));
    let f2 = f_fmla(z, f64::from_bits(B[2]), f64::from_bits(B[1]));

    let zf0 = f_fmla(z2, f64::from_bits(B[7]), f0);
    let zf1 = f_fmla(z2, f1, f2);

    let f = z2 * f_fmla(z4, zf0, zf1);
    let mut r = (z + f).to_bits();
    if (r & 0xfffffffu64) == 0 {
        r = r.wrapping_add(
            (f64::from_bits(0x40d0000000000000) * (f + (z - f64::from_bits(r)))).to_bits(),
        );
    }
    f64::from_bits(r) as f32
}

/// Computes log(x+1)
///
/// Max ULP 0.5
#[inline]
pub fn f_log1pf(x: f32) -> f32 {
    let mut z = x as f64;
    let t = x.to_bits();
    let ux: u32 = t;
    let ax = ux & 0x7fff_ffff;
    if ax < 0x3c880000u32 {
        // |x| < 0.0166015625
        log1pf_small(x, z, ax)
    } else {
        if ux >= 0xbf800000u32 || ax >= 0x7f800000u32 {
            return special_logf(x);
        }
        let tp = (z + 1.).to_bits();
        let mut e: i32 = (tp >> 52) as i32;
        let m52: u64 = tp & 0x000fffffffffffff;
        let j: usize = ((tp >> (52 - 5)) & 31) as usize;
        e -= 0x3ff;
        let xd = m52 | (0x3ffu64 << 52);
        z = f_fmla(f64::from_bits(xd), f64::from_bits(X0[j]), -1.);

        const C: [u64; 5] = [
            0xbd43902c33434e7f,
            0x3feffffffe1cbed5,
            0xbfdffffff7d1b014,
            0x3fd5564e0ed3613a,
            0xbfd0012232a00d4a,
        ];

        const LN2: f64 = f64::from_bits(0x3fe62e42fefa39ef);

        let z2 = z * z;
        let r0 = f_fmla(z, f64::from_bits(C[4]), f64::from_bits(C[3]));
        let r1 = f_fmla(z, f64::from_bits(C[2]), f64::from_bits(C[1]));

        let zr0 = f_fmla(z2, r0, r1);
        let zr1 = f_fmla(LN2, e as f64, f64::from_bits(LIXB[j]));

        let r = f_fmla(z, zr0, zr1);
        let ub = r as f32;
        let lb = (r + 2.2e-11) as f32;
        if ub != lb {
            return log1pf_accurate(x, z, e, j);
        }
        ub
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn log1pf_works() {
        assert_eq!(f_log1pf(0.0), 0.0);
        assert_eq!(f_log1pf(2.0), 1.0986123);
        assert_eq!(f_log1pf(-0.7), -1.2039728);
        assert_eq!(f_log1pf(-0.0000000000043243), -4.3243e-12);
    }
}
