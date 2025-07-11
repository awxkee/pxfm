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
use crate::common::{dd_fmla, f_fmla};
use crate::dekker::Dekker;
use crate::polyeval::f_polyeval9;
/*
  rf[64] and lf[64][2] are lookup tables for the fast path
  -ln(rf[][]) = lf[][1] + lf[][0]
  values are approximately from 1/sqrt(2) to sqrt(2)
*/
static RANGE_REDUCTION: [u64; 64] = [
    0x3ff6816818000000,
    0x3ff642c858000000,
    0x3ff6058160000000,
    0x3ff5c98828000000,
    0x3ff58ed230000000,
    0x3ff5555558000000,
    0x3ff51d07e8000000,
    0x3ff4e5e0a8000000,
    0x3ff4afd6a0000000,
    0x3ff47ae148000000,
    0x3ff446f868000000,
    0x3ff4141418000000,
    0x3ff3e22cc0000000,
    0x3ff3b13b10000000,
    0x3ff3813810000000,
    0x3ff3521cf8000000,
    0x3ff323e348000000,
    0x3ff2f684c0000000,
    0x3ff2c9fb50000000,
    0x3ff29e4128000000,
    0x3ff27350b8000000,
    0x3ff2492490000000,
    0x3ff21fb780000000,
    0x3ff1f70480000000,
    0x3ff1cf06b0000000,
    0x3ff1a7b960000000,
    0x3ff1811810000000,
    0x3ff15b1e60000000,
    0x3ff135c810000000,
    0x3ff1111110000000,
    0x3ff0ecf568000000,
    0x3ff0c97150000000,
    0x3ff0a68108000000,
    0x3ff0842108000000,
    0x3ff0624dd0000000,
    0x3ff0410410000000,
    0x3ff0204080000000,
    0x3ff0000000000000,
    0x3fef81f820000000,
    0x3fef07c1f0000000,
    0x3fee9131a8000000,
    0x3fee1e1e20000000,
    0x3fedae6078000000,
    0x3fed41d420000000,
    0x3fecd85688000000,
    0x3fec71c720000000,
    0x3fec0e0700000000,
    0x3febacf918000000,
    0x3feb4e81b8000000,
    0x3feaf286c0000000,
    0x3fea98ef60000000,
    0x3fea41a418000000,
    0x3fe9ec8e98000000,
    0x3fe9999998000000,
    0x3fe948b100000000,
    0x3fe8f9c190000000,
    0x3fe8acb910000000,
    0x3fe8618618000000,
    0x3fe8181818000000,
    0x3fe7d05f40000000,
    0x3fe78a4c80000000,
    0x3fe745d178000000,
    0x3fe702e060000000,
    0x3fe6c16c18000000,
];

static LOGP1_TABLE: [(u64, u64); 64] = [
    (0xbfd5d5bde3994000, 0xbd5f2f8281bade6a),
    (0xbfd522ae0438c000, 0x3d5c2843fdd367a4),
    (0xbfd4718dc171c000, 0xbd306c10c34c14b0),
    (0xbfd3c2526cb34000, 0x3d4cfa4e853f5890),
    (0xbfd314f1e0534000, 0xbd5ce3ac179bd856),
    (0xbfd269621934c000, 0xbd5b91f82deb8122),
    (0xbfd1bf995a9a8000, 0x3d546bbb83d7163e),
    (0xbfd1178e84a80000, 0x3d5b842e5a74bdb0),
    (0xbfd071385f4d4000, 0xbd5862715e5bb534),
    (0xbfcf991c6eb38000, 0xbd59bcbcbea0cdf8),
    (0xbfce530f10670000, 0xbd401101cb605958),
    (0xbfcd10380b658000, 0x3d50c38c81ad8f06),
    (0xbfcbd0874c3c0000, 0x3d53aa40992a6d82),
    (0xbfca93ed248b0000, 0x3d530f68780ae82e),
    (0xbfc95a5ac5f70000, 0xbd07d116989d0980),
    (0xbfc823c150518000, 0xbd51e0012ba619ca),
    (0xbfc6f0127cf58000, 0x3d454535d5671858),
    (0xbfc5bf407b540000, 0xbd5ed87db3498128),
    (0xbfc4913d94338000, 0xbd5aafde9c9fc39a),
    (0xbfc365fca3158000, 0xbd4015868c234000),
    (0xbfc23d7126ca0000, 0x3d5eff33f502c226),
    (0xbfc1178e72280000, 0x3d4b8521e874d358),
    (0xbfbfe89129dc0000, 0x3d454d75afe84568),
    (0xbfbda72783840000, 0xbd51a813f3fa7c1e),
    (0xbfbb6ac8afad0000, 0xbd56c6676f40963e),
    (0xbfb9335e4d590000, 0xbd52620b7957a7a6),
    (0xbfb700d2f4eb0000, 0x3d4f8ffee5598f38),
    (0xbfb4d31165200000, 0xbd5fab0f5bf42ca2),
    (0xbfb2aa0492470000, 0xbd37a3e970b1c3a8),
    (0xbfb08598a59e0000, 0xbd4d030435fecb50),
    (0xbfaccb7357de0000, 0x3d435084a4fb8ab8),
    (0xbfa894aa1ca00000, 0x3d432f36d60b44c4),
    (0xbfa466ae8a2e0000, 0x3d2c1bcce5be8110),
    (0xbfa0415d81e80000, 0x3d5777740b18714a),
    (0xbf98492470c80000, 0xbd4955c057693d94),
    (0xbf90205648940000, 0x3d44f71addb8be00),
    (0xbf801014f5880000, 0xbd3bcda4e198afb0),
    (0x0000000000000000, 0x3bcc000000000000),
    (0x3f8fc0a891000000, 0xbd5fe0df75092c5e),
    (0x3f9f829b1e780000, 0x3d298036ec7e0a10),
    (0x3fa7745938320000, 0x3d5ba010f49e5ff0),
    (0x3faf0a30a0120000, 0xbd53ab13c266d328),
    (0x3fb341d78b1c0000, 0xbd471798573e45d4),
    (0x3fb6f0d272e50000, 0x3d5ad32f072669fc),
    (0x3fba926d434b0000, 0xbd454e391e16ea38),
    (0x3fbe27074e2b0000, 0xbd2a302bbaf05590),
    (0x3fc0d77e8cd08000, 0x3d3cb4cd66e31f30),
    (0x3fc29552e9200000, 0xbd35b7a5bc474128),
    (0x3fc44d2b5e4b8000, 0xbd17062e8135f740),
    (0x3fc5ff3060a78000, 0x3d43d4c88fe1f4b0),
    (0x3fc7ab8904110000, 0xbd537b70004a6946),
    (0x3fc9525aa7f48000, 0xbd54a5885167c1ec),
    (0x3fcaf3c940008000, 0x3d5ff9d5953004ac),
    (0x3fcc8ff7cf9a8000, 0x3d4a21ec41d8219c),
    (0x3fce27075e2b0000, 0xbd3a322bf2f02ae8),
    (0x3fcfb9186b5e0000, 0x3d5f1548b8a33616),
    (0x3fd0a324e0f38000, 0x3d50e36401f7a006),
    (0x3fd1675cacabc000, 0xbd59f1fa55382a8a),
    (0x3fd22941fc0f8000, 0xbd3a69763deb0960),
    (0x3fd2e8e2bee10000, 0x3d5d30bc3ac91bda),
    (0x3fd3a64c59694000, 0x3d37a79cf4d73b28),
    (0x3fd4618bb81c4000, 0x3d5ec345197b22de),
    (0x3fd51aad7c2e0000, 0xbd3f4810a30aeba8),
    (0x3fd5d1bdbbd80000, 0x3d4394d2371c1d1c),
];

#[inline(always)]
pub(crate) fn log1p_tiny(z: Dekker) -> Dekker {
    /* The following is a degree-4 polynomial generated by Sollya for log1p(x) */
    const Q_1: [(u64, u64); 9] = [
        (0x0000000000000000, 0x3ff0000000000000),
        (0x0000000000000000, 0xbfe0000000000000),
        (0xbc85555555555555, 0x3fd5555555555556),
        (0x0000000000000000, 0xbfd0000000000000),
        (0xbc6999999999999a, 0x3fc999999999999a),
        (0x3c75555555555555, 0xbfc5555555555556),
        (0x3c62492492492492, 0x3fc2492492492492),
        (0x0000000000000000, 0xbfc0000000000000),
        (0x3c5c71c71c71c71c, 0x3fbc71c71c71c71c),
    ];
    Dekker::quick_mult(
        f_polyeval9(
            z,
            Dekker::from_bit_pair(Q_1[0]),
            Dekker::from_bit_pair(Q_1[1]),
            Dekker::from_bit_pair(Q_1[2]),
            Dekker::from_bit_pair(Q_1[3]),
            Dekker::from_bit_pair(Q_1[4]),
            Dekker::from_bit_pair(Q_1[5]),
            Dekker::from_bit_pair(Q_1[6]),
            Dekker::from_bit_pair(Q_1[7]),
            Dekker::from_bit_pair(Q_1[8]),
        ),
        z,
    )
}

#[inline]
pub(crate) fn log1p_f64_dd(x: f64) -> (Dekker, bool) {
    let ix = x.to_bits();
    let ax = ix.wrapping_shl(1);
    let mut ln1: f64;
    let mut ln0: f64;
    let cancel: bool;
    /* logp1 is expected to be used for x near 0, where it is more accurate than
    log(1+x), thus we expect x near 0 */
    if ax < 0x7f60000000000000u64 {
        // |x| < 0.0625
        let x2 = x * x;
        if ax < 0x7e60000000000000u64 {
            // |x| < 0x1p-12
            return (log1p_tiny(Dekker::new(0., x)), false);
        } else {
            const C: [u64; 12] = [
                0x3fd5555555555555,
                0xbfd0000000000000,
                0x3fc9999999999b41,
                0xbfc555555555583b,
                0x3fc24924923f39e0,
                0xbfbfffffffe42e43,
                0x3fbc71c75511d70b,
                0xbfb99999de10510f,
                0x3fb7457e81b175f6,
                0xbfb554fb43e54e0f,
                0x3fb3ed68744f3d18,
                0xbfb28558ad5a7ac4,
            ];
            let x3 = x2 * x;
            let x4 = x2 * x2;
            let hx = -0.5 * x;
            ln1 = dd_fmla(hx, x, x);
            ln0 = dd_fmla(hx, x, x - ln1);

            let f0 = f_fmla(x, f64::from_bits(C[11]), f64::from_bits(C[10]));
            let f1 = f_fmla(x, f64::from_bits(C[9]), f64::from_bits(C[8]));
            let f2 = f_fmla(x, f64::from_bits(C[7]), f64::from_bits(C[6]));
            let f3 = f_fmla(x, f64::from_bits(C[5]), f64::from_bits(C[4]));
            let f4 = f_fmla(x, f64::from_bits(C[3]), f64::from_bits(C[2]));
            let f5 = f_fmla(x, f64::from_bits(C[1]), f64::from_bits(C[0]));

            let f = (f5 + x2 * f4) + x4 * ((f3 + x2 * f2) + x4 * (f1 + x2 * f0));
            ln0 += x3 * f;
            let eps: f64 = f64::from_bits(0x3cb9400000000000);
            let lb = ln1 + (ln0 - eps);
            let ub = ln1 + (ln0 + eps);
            cancel = lb != ub;
        }
    } else {
        // |x| >= 0.0625
        const C: [u64; 6] = [
            0xbfe000000000003d,
            0x3fd5555555554cf5,
            0xbfcffffffeca2939,
            0x3fc99999a3661724,
            0xbfc555d345bfe6fd,
            0x3fc247b887a6e5ed,
        ];
        let t: u64;
        let dt: u64;
        if ix < 0x4340000000000000u64 || ix < 0xbff0000000000000u64 {
            /* 0.0625 < x < 0x1p+53 or -1 < x < -0.0625.*/
            let z = Dekker::from_exact_add(1.0, x);
            t = z.hi.to_bits();
            dt = z.lo.to_bits();
        } else if ix < 0x4690000000000000u64 {
            // x < 0x1p+106
            t = x.to_bits();
            dt = 1f64.to_bits();
        } else if ix < 0x7ff0000000000000u64 {
            // x < 0x1p+1024
            t = x.to_bits();
            dt = 0f64.to_bits();
        } else {
            if ax > 0xffe0000000000000u64 {
                return (Dekker::new(0., x + x), true);
            } // nan
            if ix == 0x7ff0000000000000u64 {
                return (Dekker::new(0., x), true);
            } // +inf
            if ix == 0xbff0000000000000u64 {
                // -1
                return (Dekker::new(0., -1. / 0.0), true);
            }
            return (Dekker::new(0., f64::NAN), true); // <-1
        }
        let j: i64 = (t as i64).wrapping_sub(0x3fe6a00000000000i64);
        let j1 = (j >> (52 - 6)) & 0x3f;
        let je = j >> 52;
        let eoff: u64 = (je as u64).wrapping_shl(52);
        let mut rs = RANGE_REDUCTION[j1 as usize];
        if je < 1022 {
            rs = rs.wrapping_sub(eoff);
        } else {
            rs = rs.wrapping_sub(1021 << 52);
            static SC: [u64; 3] = [0x3fe0000000000000, 0x3fd0000000000000, 0x3fc0000000000000];
            rs = (f64::from_bits(rs) * f64::from_bits(SC[(je - 1022) as usize])).to_bits();
        }
        let dh = f64::from_bits(rs) * f64::from_bits(t);
        let dl = dd_fmla(f64::from_bits(rs), f64::from_bits(t), -dh)
            + f64::from_bits(rs) * f64::from_bits(dt);
        let mut ddx = Dekker::from_exact_add(dh - 1.0, dl);
        let x2 = ddx.hi * ddx.hi;

        let f0 = f_fmla(ddx.hi, f64::from_bits(C[5]), f64::from_bits(C[4]));
        let f1 = f_fmla(ddx.hi, f64::from_bits(C[3]), f64::from_bits(C[2]));
        let f2 = f_fmla(ddx.hi, f64::from_bits(C[1]), f64::from_bits(C[0]));

        let f3 = f_fmla(x2, f0, f1);

        ddx.lo += x2 * f_fmla(x2, f3, f2);
        let bl1 = f64::from_bits(0x3fe62e42fefa4000) * je as f64;
        let bl0 = f64::from_bits(0xbd48432a1b0e2634) * je as f64;
        let logp1 = LOGP1_TABLE[j1 as usize];
        ln1 = f64::from_bits(logp1.0) + bl1;
        ln0 = f64::from_bits(logp1.1) + bl0;
        let zz = Dekker::add(Dekker::new(ln0, ln1), ddx);
        ln1 = zz.hi;
        ln0 = zz.lo;
        let eps: f64 = f64::from_bits(0x3bea000000000000);
        let lb = ln1 + (ln0 - eps);
        let ub = ln1 + (ln0 + eps);
        cancel = lb != ub;
    }
    (Dekker::new(ln0, ln1), cancel)
}

#[cfg(test)]
mod tests {
    use crate::log1p_dd::log1p_f64_dd;

    #[test]
    fn test_log1p_f64_dd() {
        println!("{:?}", log1p_f64_dd(1.225158611559834).0.to_f64());
    }
}
