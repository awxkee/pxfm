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
use crate::common::{dyad_fmla, f_fmla, min_normal_f64};
use crate::double_double::DoubleDouble;
use crate::dyadic_float::{DyadicFloat128, DyadicSign};
use crate::sincos_dyadic::{
    SIN_K_PI_OVER_128_F128, range_reduction_small_f128, sincos_eval_dyadic,
};
use crate::sincos_reduce::LargeArgumentReduction;

// For 2^-7 < |x| < 2^16, return k and u such that:
//   k = round(x * 128/pi)
//   x mod pi/128 = x - k * pi/128 ~ u.hi + u.lo
// Error bound:
//   |(x - k * pi/128) - (u_hi + u_lo)| <= max(ulp(ulp(u_hi)), 2^-119)
//                                      <= 2^-111.
#[inline]
pub(crate) fn range_reduction_small(x: f64) -> (DoubleDouble, u64) {
    const MPI_OVER_128: [u64; 3] = [0xbf9921fb54400000, 0xbd70b4611a600000, 0xbb43198a2e037073];
    const ONE_TWENTY_EIGHT_OVER_PI_D: f64 = f64::from_bits(0x40445f306dc9c883);
    let prod_hi = x * ONE_TWENTY_EIGHT_OVER_PI_D;
    let kd = prod_hi.round();

    // Let y = x - k * (pi/128)
    // Then |y| < pi / 256
    // With extra rounding errors, we can bound |y| < 1.6 * 2^-7.
    let y_hi = f_fmla(kd, f64::from_bits(MPI_OVER_128[0]), x); // Exact
    // |u.hi| < 1.6*2^-7
    let u_hi = f_fmla(kd, f64::from_bits(MPI_OVER_128[1]), y_hi);

    let u0 = y_hi - u_hi; // Exact
    // |u.lo| <= max(ulp(u.hi), |kd * MPI_OVER_128[2]|)
    let u1 = f_fmla(kd, f64::from_bits(MPI_OVER_128[1]), u0); // Exact
    let u_lo = f_fmla(kd, f64::from_bits(MPI_OVER_128[2]), u1);
    // Error bound:
    // |x - k * pi/128| - (u.hi + u.lo) <= ulp(u.lo)
    //                                  <= ulp(max(ulp(u.hi), kd*MPI_OVER_128[2]))
    //                                  <= 2^(-7 - 104) = 2^-111.
    (DoubleDouble::new(u_lo, u_hi), (kd as i64) as u64)
}

#[inline]
pub(crate) fn range_reduction_small_dd(x: DoubleDouble) -> (DoubleDouble, i64) {
    const MPI_OVER_128: [u64; 5] = [
        0xbf9921fb54400000,
        0xbd70b4611a600000,
        0xbb43198a2e000000,
        0xb91b839a25200000,
        0xb6b2704453400000,
    ];
    const ONE_TWENTY_EIGHT_OVER_PI_D: f64 = f64::from_bits(0x40445f306dc9c883);
    let prod_hi = DoubleDouble::quick_mult_f64(x, ONE_TWENTY_EIGHT_OVER_PI_D);
    let kd = prod_hi.to_f64().round();

    let p_hi = f64::from_bits(MPI_OVER_128[0]); // hi
    let p_mid = f64::from_bits(MPI_OVER_128[1]); // mid
    let p_lo = f64::from_bits(MPI_OVER_128[2]); // lo
    let p_lo_lo = f64::from_bits(MPI_OVER_128[3]); // lo_lo

    let q0 = DoubleDouble::from_exact_mult(kd, p_hi);
    let q1 = DoubleDouble::from_exact_mult(kd, p_mid);
    let q2 = DoubleDouble::from_exact_mult(kd, p_lo);
    let q3 = DoubleDouble::from_exact_mult(kd, p_lo_lo);

    let mut q = DoubleDouble::add(x, q0);
    q = DoubleDouble::add(q, q1);
    q = DoubleDouble::add(q, q2);
    q = DoubleDouble::add(q, q3);

    (q, kd as i64)
}

#[inline]
pub(crate) fn get_sin_k_rational(kk: u64) -> DyadicFloat128 {
    let idx = if (kk & 64) != 0 {
        64 - (kk & 63)
    } else {
        kk & 63
    };
    let mut ans = SIN_K_PI_OVER_128_F128[idx as usize];
    if (kk & 128) != 0 {
        ans.sign = DyadicSign::Neg;
    }
    ans
}

pub(crate) struct SinCos {
    pub(crate) v_sin: DoubleDouble,
    pub(crate) v_cos: DoubleDouble,
    pub(crate) err: f64,
}

#[inline]
pub(crate) fn sincos_eval(u: DoubleDouble) -> SinCos {
    // Evaluate sin(y) = sin(x - k * (pi/128))
    // We use the degree-7 Taylor approximation:
    //   sin(y) ~ y - y^3/3! + y^5/5! - y^7/7!
    // Then the error is bounded by:
    //   |sin(y) - (y - y^3/3! + y^5/5! - y^7/7!)| < |y|^9/9! < 2^-54/9! < 2^-72.
    // For y ~ u_hi + u_lo, fully expanding the polynomial and drop any terms
    // < ulp(u_hi^3) gives us:
    //   y - y^3/3! + y^5/5! - y^7/7! = ...
    // ~ u_hi + u_hi^3 * (-1/6 + u_hi^2 * (1/120 - u_hi^2 * 1/5040)) +
    //        + u_lo (1 + u_hi^2 * (-1/2 + u_hi^2 / 24))
    let u_hi_sq = u.hi * u.hi; // Error < ulp(u_hi^2) < 2^(-6 - 52) = 2^-58.
    // p1 ~ 1/120 + u_hi^2 / 5040.
    let p1 = f_fmla(
        u_hi_sq,
        f64::from_bits(0xbf2a01a01a01a01a),
        f64::from_bits(0x3f81111111111111),
    );
    // q1 ~ -1/2 + u_hi^2 / 24.
    let q1 = f_fmla(
        u_hi_sq,
        f64::from_bits(0x3fa5555555555555),
        f64::from_bits(0xbfe0000000000000),
    );
    let u_hi_3 = u_hi_sq * u.hi;
    // p2 ~ -1/6 + u_hi^2 (1/120 - u_hi^2 * 1/5040)
    let p2 = f_fmla(u_hi_sq, p1, f64::from_bits(0xbfc5555555555555));
    // q2 ~ 1 + u_hi^2 (-1/2 + u_hi^2 / 24)
    let q2 = f_fmla(u_hi_sq, q1, 1.0);
    let sin_lo = f_fmla(u_hi_3, p2, u.lo * q2);
    // Overall, |sin(y) - (u_hi + sin_lo)| < 2*ulp(u_hi^3) < 2^-69.

    // Evaluate cos(y) = cos(x - k * (pi/128))
    // We use the degree-8 Taylor approximation:
    //   cos(y) ~ 1 - y^2/2 + y^4/4! - y^6/6! + y^8/8!
    // Then the error is bounded by:
    //   |cos(y) - (...)| < |y|^10/10! < 2^-81
    // For y ~ u_hi + u_lo, fully expanding the polynomial and drop any terms
    // < ulp(u_hi^3) gives us:
    //   1 - y^2/2 + y^4/4! - y^6/6! + y^8/8! = ...
    // ~ 1 - u_hi^2/2 + u_hi^4(1/24 + u_hi^2 (-1/720 + u_hi^2/40320)) +
    //     + u_hi u_lo (-1 + u_hi^2/6)
    // We compute 1 - u_hi^2 accurately:
    //   v_hi + v_lo ~ 1 - u_hi^2/2
    // with error <= 2^-105.
    let u_hi_neg_half = (-0.5) * u.hi;

    let (mut v_lo, v_hi);

    #[cfg(any(
        all(
            any(target_arch = "x86", target_arch = "x86_64"),
            target_feature = "fma"
        ),
        all(target_arch = "aarch64", target_feature = "neon")
    ))]
    {
        v_hi = f_fmla(u.hi, u_hi_neg_half, 1.0);
        v_lo = 1.0 - v_hi; // Exact
        v_lo = f_fmla(u.hi, u_hi_neg_half, v_lo);
    }

    #[cfg(not(any(
        all(
            any(target_arch = "x86", target_arch = "x86_64"),
            target_feature = "fma"
        ),
        all(target_arch = "aarch64", target_feature = "neon")
    )))]
    {
        let u_hi_sq_neg_half = DoubleDouble::from_exact_mult(u.hi, u_hi_neg_half);
        let v = DoubleDouble::from_exact_add(1.0, u_hi_sq_neg_half.hi);
        v_lo = v.lo;
        v_lo += u_hi_sq_neg_half.lo;
        v_hi = v.hi;
    }

    // r1 ~ -1/720 + u_hi^2 / 40320
    let r1 = f_fmla(
        u_hi_sq,
        f64::from_bits(0x3efa01a01a01a01a),
        f64::from_bits(0xbf56c16c16c16c17),
    );
    // s1 ~ -1 + u_hi^2 / 6
    let s1 = f_fmla(u_hi_sq, f64::from_bits(0x3fc5555555555555), -1.0);
    let u_hi_4 = u_hi_sq * u_hi_sq;
    let u_hi_u_lo = u.hi * u.lo;
    // r2 ~ 1/24 + u_hi^2 (-1/720 + u_hi^2 / 40320)
    let r2 = f_fmla(u_hi_sq, r1, f64::from_bits(0x3fa5555555555555));
    // s2 ~ v_lo + u_hi * u_lo * (-1 + u_hi^2 / 6)
    let s2 = f_fmla(u_hi_u_lo, s1, v_lo);
    let cos_lo = f_fmla(u_hi_4, r2, s2);
    // Overall, |cos(y) - (v_hi + cos_lo)| < 2*ulp(u_hi^4) < 2^-75.

    let sin_u = DoubleDouble::from_exact_add(u.hi, sin_lo);
    let cos_u = DoubleDouble::from_exact_add(v_hi, cos_lo);

    let err = f_fmla(
        u_hi_3,
        f64::from_bits(0x3cc0000000000000),
        f64::from_bits(0x3960000000000000),
    );

    SinCos {
        v_sin: sin_u,
        v_cos: cos_u,
        err,
    }
}

pub(crate) static SIN_K_PI_OVER_128: [(u64, u64); 256] = [
    (0x0000000000000000, 0x0000000000000000),
    (0xbbfb1d63091a0130, 0x3f992155f7a3667e),
    (0xbc2912bd0d569a90, 0x3fa91f65f10dd814),
    (0xbc49a088a8bf6b2c, 0x3fb2d52092ce19f6),
    (0xbc3e2718d26ed688, 0x3fb917a6bc29b42c),
    (0x3c4a2704729ae56d, 0x3fbf564e56a9730e),
    (0x3c513000a89a11e0, 0x3fc2c8106e8e613a),
    (0x3c6531ff779ddac6, 0x3fc5e214448b3fc6),
    (0xbc626d19b9ff8d82, 0x3fc8f8b83c69a60b),
    (0xbc1af1439e521935, 0x3fcc0b826a7e4f63),
    (0xbc642deef11da2c4, 0x3fcf19f97b215f1b),
    (0x3c7824c20ab7aa9a, 0x3fd111d262b1f677),
    (0xbc75d28da2c4612d, 0x3fd294062ed59f06),
    (0x3c70c97c4afa2518, 0x3fd4135c94176601),
    (0xbc1efdc0d58cf620, 0x3fd58f9a75ab1fdd),
    (0xbc744b19e0864c5d, 0x3fd7088530fa459f),
    (0xbc672cedd3d5a610, 0x3fd87de2a6aea963),
    (0x3c66da81290bdbab, 0x3fd9ef7943a8ed8a),
    (0x3c65b362cb974183, 0x3fdb5d1009e15cc0),
    (0x3c56850e59c37f8f, 0x3fdcc66e9931c45e),
    (0x3c5e0d891d3c6841, 0x3fde2b5d3806f63b),
    (0xbc32ec1fc1b776b8, 0x3fdf8ba4dbf89aba),
    (0xbc8a5a014347406c, 0x3fe073879922ffee),
    (0xbc8ef23b69abe4f1, 0x3fe11eb3541b4b23),
    (0x3c8b25dd267f6600, 0x3fe1c73b39ae68c8),
    (0xbc85da743ef3770c, 0x3fe26d054cdd12df),
    (0xbc6efcc626f74a6f, 0x3fe30ff7fce17035),
    (0x3c7e3e25e3954964, 0x3fe3affa292050b9),
    (0x3c68076a2cfdc6b3, 0x3fe44cf325091dd6),
    (0x3c63c293edceb327, 0x3fe4e6cabbe3e5e9),
    (0xbc875720992bfbb2, 0x3fe57d69348ceca0),
    (0xbc7251b352ff2a37, 0x3fe610b7551d2cdf),
    (0xbc8bdd3413b26456, 0x3fe6a09e667f3bcd),
    (0x3c80d4ef0f1d915c, 0x3fe72d0837efff96),
    (0xbc70f537acdf0ad7, 0x3fe7b5df226aafaf),
    (0xbc76f420f8ea3475, 0x3fe83b0e0bff976e),
    (0xbc82c5e12ed1336d, 0x3fe8bc806b151741),
    (0x3c83d419a920df0b, 0x3fe93a22499263fb),
    (0xbc830ee286712474, 0x3fe9b3e047f38741),
    (0xbc7128bb015df175, 0x3fea29a7a0462782),
    (0x3c39f630e8b6dac8, 0x3fea9b66290ea1a3),
    (0xbc8926da300ffcce, 0x3feb090a58150200),
    (0xbc8bc69f324e6d61, 0x3feb728345196e3e),
    (0xbc8825a732ac700a, 0x3febd7c0ac6f952a),
    (0xbc76e0b1757c8d07, 0x3fec38b2f180bdb1),
    (0xbc52fb761e946603, 0x3fec954b213411f5),
    (0xbc5e7b6bb5ab58ae, 0x3feced7af43cc773),
    (0xbc84ef5295d25af2, 0x3fed4134d14dc93a),
    (0x3c7457e610231ac2, 0x3fed906bcf328d46),
    (0x3c883c37c6107db3, 0x3feddb13b6ccc23c),
    (0xbc8014c76c126527, 0x3fee212104f686e5),
    (0xbc616b56f2847754, 0x3fee6288ec48e112),
    (0x3c8760b1e2e3f81e, 0x3fee9f4156c62dda),
    (0x3c7e82c791f59cc2, 0x3feed740e7684963),
    (0x3c752c7adc6b4989, 0x3fef0a7efb9230d7),
    (0xbc7d7bafb51f72e6, 0x3fef38f3ac64e589),
    (0x3c7562172a361fd3, 0x3fef6297cff75cb0),
    (0x3c7ab256778ffcb6, 0x3fef8764fa714ba9),
    (0xbc87a0a8ca13571f, 0x3fefa7557f08a517),
    (0x3c81ec8668ecacee, 0x3fefc26470e19fd3),
    (0xbc887df6378811c7, 0x3fefd88da3d12526),
    (0x3c6521ecd0c67e35, 0x3fefe9cdad01883a),
    (0xbc6c57bc2e24aa15, 0x3feff621e3796d7e),
    (0xbc81354d4556e4cb, 0x3feffd886084cd0d),
    (0x0000000000000000, 0x3ff0000000000000),
    (0xbc81354d4556e4cb, 0x3feffd886084cd0d),
    (0xbc6c57bc2e24aa15, 0x3feff621e3796d7e),
    (0x3c6521ecd0c67e35, 0x3fefe9cdad01883a),
    (0xbc887df6378811c7, 0x3fefd88da3d12526),
    (0x3c81ec8668ecacee, 0x3fefc26470e19fd3),
    (0xbc87a0a8ca13571f, 0x3fefa7557f08a517),
    (0x3c7ab256778ffcb6, 0x3fef8764fa714ba9),
    (0x3c7562172a361fd3, 0x3fef6297cff75cb0),
    (0xbc7d7bafb51f72e6, 0x3fef38f3ac64e589),
    (0x3c752c7adc6b4989, 0x3fef0a7efb9230d7),
    (0x3c7e82c791f59cc2, 0x3feed740e7684963),
    (0x3c8760b1e2e3f81e, 0x3fee9f4156c62dda),
    (0xbc616b56f2847754, 0x3fee6288ec48e112),
    (0xbc8014c76c126527, 0x3fee212104f686e5),
    (0x3c883c37c6107db3, 0x3feddb13b6ccc23c),
    (0x3c7457e610231ac2, 0x3fed906bcf328d46),
    (0xbc84ef5295d25af2, 0x3fed4134d14dc93a),
    (0xbc5e7b6bb5ab58ae, 0x3feced7af43cc773),
    (0xbc52fb761e946603, 0x3fec954b213411f5),
    (0xbc76e0b1757c8d07, 0x3fec38b2f180bdb1),
    (0xbc8825a732ac700a, 0x3febd7c0ac6f952a),
    (0xbc8bc69f324e6d61, 0x3feb728345196e3e),
    (0xbc8926da300ffcce, 0x3feb090a58150200),
    (0x3c39f630e8b6dac8, 0x3fea9b66290ea1a3),
    (0xbc7128bb015df175, 0x3fea29a7a0462782),
    (0xbc830ee286712474, 0x3fe9b3e047f38741),
    (0x3c83d419a920df0b, 0x3fe93a22499263fb),
    (0xbc82c5e12ed1336d, 0x3fe8bc806b151741),
    (0xbc76f420f8ea3475, 0x3fe83b0e0bff976e),
    (0xbc70f537acdf0ad7, 0x3fe7b5df226aafaf),
    (0x3c80d4ef0f1d915c, 0x3fe72d0837efff96),
    (0xbc8bdd3413b26456, 0x3fe6a09e667f3bcd),
    (0xbc7251b352ff2a37, 0x3fe610b7551d2cdf),
    (0xbc875720992bfbb2, 0x3fe57d69348ceca0),
    (0x3c63c293edceb327, 0x3fe4e6cabbe3e5e9),
    (0x3c68076a2cfdc6b3, 0x3fe44cf325091dd6),
    (0x3c7e3e25e3954964, 0x3fe3affa292050b9),
    (0xbc6efcc626f74a6f, 0x3fe30ff7fce17035),
    (0xbc85da743ef3770c, 0x3fe26d054cdd12df),
    (0x3c8b25dd267f6600, 0x3fe1c73b39ae68c8),
    (0xbc8ef23b69abe4f1, 0x3fe11eb3541b4b23),
    (0xbc8a5a014347406c, 0x3fe073879922ffee),
    (0xbc32ec1fc1b776b8, 0x3fdf8ba4dbf89aba),
    (0x3c5e0d891d3c6841, 0x3fde2b5d3806f63b),
    (0x3c56850e59c37f8f, 0x3fdcc66e9931c45e),
    (0x3c65b362cb974183, 0x3fdb5d1009e15cc0),
    (0x3c66da81290bdbab, 0x3fd9ef7943a8ed8a),
    (0xbc672cedd3d5a610, 0x3fd87de2a6aea963),
    (0xbc744b19e0864c5d, 0x3fd7088530fa459f),
    (0xbc1efdc0d58cf620, 0x3fd58f9a75ab1fdd),
    (0x3c70c97c4afa2518, 0x3fd4135c94176601),
    (0xbc75d28da2c4612d, 0x3fd294062ed59f06),
    (0x3c7824c20ab7aa9a, 0x3fd111d262b1f677),
    (0xbc642deef11da2c4, 0x3fcf19f97b215f1b),
    (0xbc1af1439e521935, 0x3fcc0b826a7e4f63),
    (0xbc626d19b9ff8d82, 0x3fc8f8b83c69a60b),
    (0x3c6531ff779ddac6, 0x3fc5e214448b3fc6),
    (0x3c513000a89a11e0, 0x3fc2c8106e8e613a),
    (0x3c4a2704729ae56d, 0x3fbf564e56a9730e),
    (0xbc3e2718d26ed688, 0x3fb917a6bc29b42c),
    (0xbc49a088a8bf6b2c, 0x3fb2d52092ce19f6),
    (0xbc2912bd0d569a90, 0x3fa91f65f10dd814),
    (0xbbfb1d63091a0130, 0x3f992155f7a3667e),
    (0x0000000000000000, 0x0000000000000000),
    (0x3bfb1d63091a0130, 0xbf992155f7a3667e),
    (0x3c2912bd0d569a90, 0xbfa91f65f10dd814),
    (0x3c49a088a8bf6b2c, 0xbfb2d52092ce19f6),
    (0x3c3e2718d26ed688, 0xbfb917a6bc29b42c),
    (0xbc4a2704729ae56d, 0xbfbf564e56a9730e),
    (0xbc513000a89a11e0, 0xbfc2c8106e8e613a),
    (0xbc6531ff779ddac6, 0xbfc5e214448b3fc6),
    (0x3c626d19b9ff8d82, 0xbfc8f8b83c69a60b),
    (0x3c1af1439e521935, 0xbfcc0b826a7e4f63),
    (0x3c642deef11da2c4, 0xbfcf19f97b215f1b),
    (0xbc7824c20ab7aa9a, 0xbfd111d262b1f677),
    (0x3c75d28da2c4612d, 0xbfd294062ed59f06),
    (0xbc70c97c4afa2518, 0xbfd4135c94176601),
    (0x3c1efdc0d58cf620, 0xbfd58f9a75ab1fdd),
    (0x3c744b19e0864c5d, 0xbfd7088530fa459f),
    (0x3c672cedd3d5a610, 0xbfd87de2a6aea963),
    (0xbc66da81290bdbab, 0xbfd9ef7943a8ed8a),
    (0xbc65b362cb974183, 0xbfdb5d1009e15cc0),
    (0xbc56850e59c37f8f, 0xbfdcc66e9931c45e),
    (0xbc5e0d891d3c6841, 0xbfde2b5d3806f63b),
    (0x3c32ec1fc1b776b8, 0xbfdf8ba4dbf89aba),
    (0x3c8a5a014347406c, 0xbfe073879922ffee),
    (0x3c8ef23b69abe4f1, 0xbfe11eb3541b4b23),
    (0xbc8b25dd267f6600, 0xbfe1c73b39ae68c8),
    (0x3c85da743ef3770c, 0xbfe26d054cdd12df),
    (0x3c6efcc626f74a6f, 0xbfe30ff7fce17035),
    (0xbc7e3e25e3954964, 0xbfe3affa292050b9),
    (0xbc68076a2cfdc6b3, 0xbfe44cf325091dd6),
    (0xbc63c293edceb327, 0xbfe4e6cabbe3e5e9),
    (0x3c875720992bfbb2, 0xbfe57d69348ceca0),
    (0x3c7251b352ff2a37, 0xbfe610b7551d2cdf),
    (0x3c8bdd3413b26456, 0xbfe6a09e667f3bcd),
    (0xbc80d4ef0f1d915c, 0xbfe72d0837efff96),
    (0x3c70f537acdf0ad7, 0xbfe7b5df226aafaf),
    (0x3c76f420f8ea3475, 0xbfe83b0e0bff976e),
    (0x3c82c5e12ed1336d, 0xbfe8bc806b151741),
    (0xbc83d419a920df0b, 0xbfe93a22499263fb),
    (0x3c830ee286712474, 0xbfe9b3e047f38741),
    (0x3c7128bb015df175, 0xbfea29a7a0462782),
    (0xbc39f630e8b6dac8, 0xbfea9b66290ea1a3),
    (0x3c8926da300ffcce, 0xbfeb090a58150200),
    (0x3c8bc69f324e6d61, 0xbfeb728345196e3e),
    (0x3c8825a732ac700a, 0xbfebd7c0ac6f952a),
    (0x3c76e0b1757c8d07, 0xbfec38b2f180bdb1),
    (0x3c52fb761e946603, 0xbfec954b213411f5),
    (0x3c5e7b6bb5ab58ae, 0xbfeced7af43cc773),
    (0x3c84ef5295d25af2, 0xbfed4134d14dc93a),
    (0xbc7457e610231ac2, 0xbfed906bcf328d46),
    (0xbc883c37c6107db3, 0xbfeddb13b6ccc23c),
    (0x3c8014c76c126527, 0xbfee212104f686e5),
    (0x3c616b56f2847754, 0xbfee6288ec48e112),
    (0xbc8760b1e2e3f81e, 0xbfee9f4156c62dda),
    (0xbc7e82c791f59cc2, 0xbfeed740e7684963),
    (0xbc752c7adc6b4989, 0xbfef0a7efb9230d7),
    (0x3c7d7bafb51f72e6, 0xbfef38f3ac64e589),
    (0xbc7562172a361fd3, 0xbfef6297cff75cb0),
    (0xbc7ab256778ffcb6, 0xbfef8764fa714ba9),
    (0x3c87a0a8ca13571f, 0xbfefa7557f08a517),
    (0xbc81ec8668ecacee, 0xbfefc26470e19fd3),
    (0x3c887df6378811c7, 0xbfefd88da3d12526),
    (0xbc6521ecd0c67e35, 0xbfefe9cdad01883a),
    (0x3c6c57bc2e24aa15, 0xbfeff621e3796d7e),
    (0x3c81354d4556e4cb, 0xbfeffd886084cd0d),
    (0x0000000000000000, 0xbff0000000000000),
    (0x3c81354d4556e4cb, 0xbfeffd886084cd0d),
    (0x3c6c57bc2e24aa15, 0xbfeff621e3796d7e),
    (0xbc6521ecd0c67e35, 0xbfefe9cdad01883a),
    (0x3c887df6378811c7, 0xbfefd88da3d12526),
    (0xbc81ec8668ecacee, 0xbfefc26470e19fd3),
    (0x3c87a0a8ca13571f, 0xbfefa7557f08a517),
    (0xbc7ab256778ffcb6, 0xbfef8764fa714ba9),
    (0xbc7562172a361fd3, 0xbfef6297cff75cb0),
    (0x3c7d7bafb51f72e6, 0xbfef38f3ac64e589),
    (0xbc752c7adc6b4989, 0xbfef0a7efb9230d7),
    (0xbc7e82c791f59cc2, 0xbfeed740e7684963),
    (0xbc8760b1e2e3f81e, 0xbfee9f4156c62dda),
    (0x3c616b56f2847754, 0xbfee6288ec48e112),
    (0x3c8014c76c126527, 0xbfee212104f686e5),
    (0xbc883c37c6107db3, 0xbfeddb13b6ccc23c),
    (0xbc7457e610231ac2, 0xbfed906bcf328d46),
    (0x3c84ef5295d25af2, 0xbfed4134d14dc93a),
    (0x3c5e7b6bb5ab58ae, 0xbfeced7af43cc773),
    (0x3c52fb761e946603, 0xbfec954b213411f5),
    (0x3c76e0b1757c8d07, 0xbfec38b2f180bdb1),
    (0x3c8825a732ac700a, 0xbfebd7c0ac6f952a),
    (0x3c8bc69f324e6d61, 0xbfeb728345196e3e),
    (0x3c8926da300ffcce, 0xbfeb090a58150200),
    (0xbc39f630e8b6dac8, 0xbfea9b66290ea1a3),
    (0x3c7128bb015df175, 0xbfea29a7a0462782),
    (0x3c830ee286712474, 0xbfe9b3e047f38741),
    (0xbc83d419a920df0b, 0xbfe93a22499263fb),
    (0x3c82c5e12ed1336d, 0xbfe8bc806b151741),
    (0x3c76f420f8ea3475, 0xbfe83b0e0bff976e),
    (0x3c70f537acdf0ad7, 0xbfe7b5df226aafaf),
    (0xbc80d4ef0f1d915c, 0xbfe72d0837efff96),
    (0x3c8bdd3413b26456, 0xbfe6a09e667f3bcd),
    (0x3c7251b352ff2a37, 0xbfe610b7551d2cdf),
    (0x3c875720992bfbb2, 0xbfe57d69348ceca0),
    (0xbc63c293edceb327, 0xbfe4e6cabbe3e5e9),
    (0xbc68076a2cfdc6b3, 0xbfe44cf325091dd6),
    (0xbc7e3e25e3954964, 0xbfe3affa292050b9),
    (0x3c6efcc626f74a6f, 0xbfe30ff7fce17035),
    (0x3c85da743ef3770c, 0xbfe26d054cdd12df),
    (0xbc8b25dd267f6600, 0xbfe1c73b39ae68c8),
    (0x3c8ef23b69abe4f1, 0xbfe11eb3541b4b23),
    (0x3c8a5a014347406c, 0xbfe073879922ffee),
    (0x3c32ec1fc1b776b8, 0xbfdf8ba4dbf89aba),
    (0xbc5e0d891d3c6841, 0xbfde2b5d3806f63b),
    (0xbc56850e59c37f8f, 0xbfdcc66e9931c45e),
    (0xbc65b362cb974183, 0xbfdb5d1009e15cc0),
    (0xbc66da81290bdbab, 0xbfd9ef7943a8ed8a),
    (0x3c672cedd3d5a610, 0xbfd87de2a6aea963),
    (0x3c744b19e0864c5d, 0xbfd7088530fa459f),
    (0x3c1efdc0d58cf620, 0xbfd58f9a75ab1fdd),
    (0xbc70c97c4afa2518, 0xbfd4135c94176601),
    (0x3c75d28da2c4612d, 0xbfd294062ed59f06),
    (0xbc7824c20ab7aa9a, 0xbfd111d262b1f677),
    (0x3c642deef11da2c4, 0xbfcf19f97b215f1b),
    (0x3c1af1439e521935, 0xbfcc0b826a7e4f63),
    (0x3c626d19b9ff8d82, 0xbfc8f8b83c69a60b),
    (0xbc6531ff779ddac6, 0xbfc5e214448b3fc6),
    (0xbc513000a89a11e0, 0xbfc2c8106e8e613a),
    (0xbc4a2704729ae56d, 0xbfbf564e56a9730e),
    (0x3c3e2718d26ed688, 0xbfb917a6bc29b42c),
    (0x3c49a088a8bf6b2c, 0xbfb2d52092ce19f6),
    (0x3c2912bd0d569a90, 0xbfa91f65f10dd814),
    (0x3bfb1d63091a0130, 0xbf992155f7a3667e),
];

pub(crate) fn sin_dd_small(z: DoubleDouble) -> DoubleDouble {
    let x_e = (z.hi.to_bits() >> 52) & 0x7ff;
    const E_BIAS: u64 = (1u64 << (11 - 1u64)) - 1u64;

    if x_e < E_BIAS - 26 {
        let q = DoubleDouble::quick_mult_f64(z, f64::from_bits(0xbc90000000000000));
        return DoubleDouble::add(z, q);
    }

    let (u_f128, k) = range_reduction_small_dd(z);

    let sin_cos = sincos_eval_dd(u_f128);

    // Fast look up version, but needs 256-entry table.
    // cos(k * pi/128) = sin(k * pi/128 + pi/2) = sin((k + 64) * pi/128).
    let sk = SIN_K_PI_OVER_128[(k & 255) as usize];
    let ck = SIN_K_PI_OVER_128[((k.wrapping_add(64)) & 255) as usize];

    let sin_k = DoubleDouble::new(f64::from_bits(sk.0), f64::from_bits(sk.1));
    let cos_k = DoubleDouble::new(f64::from_bits(ck.0), f64::from_bits(ck.1));

    let sin_k_cos_y = DoubleDouble::quick_mult(sin_cos.v_cos, sin_k);
    let cos_k_sin_y = DoubleDouble::quick_mult(sin_cos.v_sin, cos_k);

    let mut rr = DoubleDouble::from_exact_add(sin_k_cos_y.hi, cos_k_sin_y.hi);
    rr.lo += sin_k_cos_y.lo + cos_k_sin_y.lo;
    rr
}

pub(crate) fn cos_dd_small(z: DoubleDouble) -> DoubleDouble {
    let x_e = (z.hi.to_bits() >> 52) & 0x7ff;
    const E_BIAS: u64 = (1u64 << (11 - 1u64)) - 1u64;

    if x_e < E_BIAS - 26 {
        let q = DoubleDouble::quick_mult_f64(z, f64::from_bits(0xbc90000000000000));
        return DoubleDouble::add(z, q);
    }

    let (u_f128, k) = range_reduction_small_dd(z);

    let sin_cos = sincos_eval_dd(u_f128);

    // cos(k * pi/128) = sin(k * pi/128 + pi/2) = sin((k + 64) * pi/128).
    let sk = SIN_K_PI_OVER_128[(k.wrapping_add(128) & 255) as usize];
    let ck = SIN_K_PI_OVER_128[((k.wrapping_add(64)) & 255) as usize];
    let msin_k = DoubleDouble::new(f64::from_bits(sk.0), f64::from_bits(sk.1));
    let cos_k = DoubleDouble::new(f64::from_bits(ck.0), f64::from_bits(ck.1));

    let sin_k_cos_y = DoubleDouble::quick_mult(sin_cos.v_cos, cos_k);
    let cos_k_sin_y = DoubleDouble::quick_mult(sin_cos.v_sin, msin_k);

    let mut rr = DoubleDouble::from_full_exact_add(sin_k_cos_y.hi, cos_k_sin_y.hi);
    rr.lo += sin_k_cos_y.lo + cos_k_sin_y.lo;
    rr
}

// pub(crate) fn sin_f128_small(z: DyadicFloat128) -> DyadicFloat128 {
//     let (u_f128, k) = range_reduction_small_f128_f128(z);
//
//     let sin_cos = sincos_eval_dyadic(&u_f128);
//     // cos(k * pi/128) = sin(k * pi/128 + pi/2) = sin((k + 64) * pi/128).
//     let sin_k_f128 = get_sin_k_rational(k);
//     let cos_k_f128 = get_sin_k_rational(k.wrapping_add(64));
//
//     // sin(x) = sin(k * pi/128 + u)
//     //        = sin(u) * cos(k*pi/128) + cos(u) * sin(k*pi/128)
//     (sin_k_f128 * sin_cos.v_cos) + (cos_k_f128 * sin_cos.v_sin)
// }

#[inline]
pub(crate) fn sin_small(z: f64) -> f64 {
    let x_e = (z.to_bits() >> 52) & 0x7ff;
    const E_BIAS: u64 = (1u64 << (11 - 1u64)) - 1u64;

    if x_e < E_BIAS - 26 {
        return dyad_fmla(z, f64::from_bits(0xbc90000000000000), z);
    }

    let (u_f128, k) = range_reduction_small(z);

    let sin_cos = sincos_eval(u_f128);

    // Fast look up version, but needs 256-entry table.
    // cos(k * pi/128) = sin(k * pi/128 + pi/2) = sin((k + 64) * pi/128).
    let sk = SIN_K_PI_OVER_128[(k & 255) as usize];
    let ck = SIN_K_PI_OVER_128[((k.wrapping_add(64)) & 255) as usize];

    let sin_k = DoubleDouble::new(f64::from_bits(sk.0), f64::from_bits(sk.1));
    let cos_k = DoubleDouble::new(f64::from_bits(ck.0), f64::from_bits(ck.1));

    let sin_k_cos_y = DoubleDouble::quick_mult(sin_cos.v_cos, sin_k);
    let cos_k_sin_y = DoubleDouble::quick_mult(sin_cos.v_sin, cos_k);

    let mut rr = DoubleDouble::from_full_exact_add(sin_k_cos_y.hi, cos_k_sin_y.hi);
    rr.lo += sin_k_cos_y.lo + cos_k_sin_y.lo;
    rr.to_f64()
}

#[inline]
pub(crate) fn cos_small(z: f64) -> f64 {
    let x_e = (z.to_bits() >> 52) & 0x7ff;
    const E_BIAS: u64 = (1u64 << (11 - 1u64)) - 1u64;

    if x_e < E_BIAS - 27 {
        // Signed zeros.
        if z == 0.0 {
            return 1.0;
        }
        // For |x| < 2^-26, |sin(x) - x| < ulp(x)/2.
        return 1.0 - min_normal_f64();
    }

    let (u_f128, k) = range_reduction_small(z);

    let sin_cos = sincos_eval(u_f128);

    // cos(k * pi/128) = sin(k * pi/128 + pi/2) = sin((k + 64) * pi/128).
    let sk = SIN_K_PI_OVER_128[(k.wrapping_add(128) & 255) as usize];
    let ck = SIN_K_PI_OVER_128[((k.wrapping_add(64)) & 255) as usize];
    let msin_k = DoubleDouble::new(f64::from_bits(sk.0), f64::from_bits(sk.1));
    let cos_k = DoubleDouble::new(f64::from_bits(ck.0), f64::from_bits(ck.1));

    let sin_k_cos_y = DoubleDouble::quick_mult(sin_cos.v_cos, cos_k);
    let cos_k_sin_y = DoubleDouble::quick_mult(sin_cos.v_sin, msin_k);

    let mut rr = DoubleDouble::from_full_exact_add(sin_k_cos_y.hi, cos_k_sin_y.hi);
    rr.lo += sin_k_cos_y.lo + cos_k_sin_y.lo;
    rr.to_f64()
}

#[inline]
pub(crate) fn sincos_eval_dd(u: DoubleDouble) -> SinCos {
    const SIN_C: [(u64, u64); 5] = [
        (0x0000000000000000, 0x3ff0000000000000),
        (0xbc65555555555555, 0xbfc5555555555555),
        (0x3c01111111111111, 0x3f81111111111111),
        (0xbb6a01a01a01a01a, 0xbf2a01a01a01a01a),
        (0xbb6c154f8ddc6c00, 0x3ec71de3a556c734),
    ];
    let x2 = DoubleDouble::quick_mult(u, u);
    let mut p = DoubleDouble::mul_add(
        x2,
        DoubleDouble::from_bit_pair(SIN_C[4]),
        DoubleDouble::from_bit_pair(SIN_C[3]),
    );

    p = DoubleDouble::mul_add(x2, p, DoubleDouble::from_bit_pair(SIN_C[2]));
    p = DoubleDouble::mul_add(x2, p, DoubleDouble::from_bit_pair(SIN_C[1]));
    p = DoubleDouble::mul_add(x2, p, DoubleDouble::from_bit_pair(SIN_C[0]));
    let sin_u = DoubleDouble::quick_mult(p, u);

    const COS_C: [(u64, u64); 5] = [
        (0x0000000000000000, 0x3ff0000000000000),
        (0x0000000000000000, 0xbfe0000000000000),
        (0x3c45555555555555, 0x3fa5555555555555),
        (0x3bef49f49f49f49f, 0xbf56c16c16c16c17),
        (0x3b3a01a01a01a01a, 0x3efa01a01a01a01a),
    ];

    let mut p = DoubleDouble::mul_add(
        x2,
        DoubleDouble::from_bit_pair(COS_C[4]),
        DoubleDouble::from_bit_pair(COS_C[3]),
    );

    p = DoubleDouble::mul_add(x2, p, DoubleDouble::from_bit_pair(COS_C[2]));
    p = DoubleDouble::mul_add(x2, p, DoubleDouble::from_bit_pair(COS_C[1]));
    p = DoubleDouble::mul_add(x2, p, DoubleDouble::from_bit_pair(COS_C[0]));

    let cos_u = p;
    SinCos {
        v_sin: sin_u,
        v_cos: cos_u,
        err: 0.,
    }
}

// pub(crate) fn cos_dd(z: Dekker) -> Dekker {
//     const EXP_BIAS: u64 = (1u64 << (11 - 1u64)) - 1u64;
//     let (u_f128, k) = range_reduction_small_dd(z);
//
//     let sin_cos = sincos_eval(u_f128);
//
//     let sk = SIN_K_PI_OVER_128[(k.wrapping_add(128) & 255) as usize];
//     let ck = SIN_K_PI_OVER_128[((k.wrapping_add(64)) & 255) as usize];
//     let msin_k = Dekker::new(f64::from_bits(sk.0), f64::from_bits(sk.1));
//     let cos_k = Dekker::new(f64::from_bits(ck.0), f64::from_bits(ck.1));
//
//     // sin(x) = sin(k * pi/128 + u)
//     //        = sin(u) * cos(k*pi/128) + cos(u) * sin(k*pi/128)
//     let sin_k_cos_y = Dekker::quick_mult(sin_cos.v_cos, cos_k);
//     let cos_k_sin_y = Dekker::quick_mult(sin_cos.v_sin, msin_k);
//
//     let mut rr = Dekker::from_exact_add(sin_k_cos_y.hi, cos_k_sin_y.hi);
//     rr.lo += sin_k_cos_y.lo + cos_k_sin_y.lo;
//     rr
// }

/// Sine for double precision
///
/// ULP 0.5
#[inline]
pub fn f_sin(x: f64) -> f64 {
    let x_e = (x.to_bits() >> 52) & 0x7ff;
    const E_BIAS: u64 = (1u64 << (11 - 1u64)) - 1u64;

    let y: DoubleDouble;
    let k;

    let mut argument_reduction = LargeArgumentReduction::default();

    // |x| < 2^32 (with FMA) or |x| < 2^23 (w/o FMA)
    if x_e < E_BIAS + 16 {
        // |x| < 2^-26
        if x_e < E_BIAS - 26 {
            // Signed zeros.
            if x == 0.0 {
                return x;
            }

            // For |x| < 2^-26, |sin(x) - x| < ulp(x)/2.
            return dyad_fmla(x, f64::from_bits(0xbc90000000000000), x);
        }

        // // Small range reduction.
        (y, k) = range_reduction_small(x);
    } else {
        // Inf or NaN
        if x_e > 2 * E_BIAS {
            // sin(+-Inf) = NaN
            return x + f64::NAN;
        }

        // Large range reduction.
        (k, y) = argument_reduction.reduce(x);
    }

    let r_sincos = sincos_eval(y);

    // Fast look up version, but needs 256-entry table.
    // cos(k * pi/128) = sin(k * pi/128 + pi/2) = sin((k + 64) * pi/128).
    let sk = SIN_K_PI_OVER_128[(k & 255) as usize];
    let ck = SIN_K_PI_OVER_128[((k.wrapping_add(64)) & 255) as usize];

    let sin_k = DoubleDouble::new(f64::from_bits(sk.0), f64::from_bits(sk.1));
    let cos_k = DoubleDouble::new(f64::from_bits(ck.0), f64::from_bits(ck.1));

    let sin_k_cos_y = DoubleDouble::quick_mult(r_sincos.v_cos, sin_k);
    let cos_k_sin_y = DoubleDouble::quick_mult(r_sincos.v_sin, cos_k);

    let mut rr = DoubleDouble::from_full_exact_add(sin_k_cos_y.hi, cos_k_sin_y.hi);
    rr.lo += sin_k_cos_y.lo + cos_k_sin_y.lo;

    let rlp = rr.lo + r_sincos.err;
    let rlm = rr.lo - r_sincos.err;

    let r_upper = rr.hi + rlp; // (rr.lo + ERR);
    let r_lower = rr.hi + rlm; // (rr.lo - ERR);

    // Ziv's accuracy test
    if r_upper == r_lower {
        return rr.to_f64();
    }

    const EXP_BIAS: u64 = (1u64 << (11 - 1u64)) - 1u64;
    let u_f128 = if x_e < EXP_BIAS + 16 {
        range_reduction_small_f128(x)
    } else {
        argument_reduction.accurate()
    };

    let sin_cos = sincos_eval_dyadic(&u_f128);

    // cos(k * pi/128) = sin(k * pi/128 + pi/2) = sin((k + 64) * pi/128).
    let sin_k_f128 = get_sin_k_rational(k);
    let cos_k_f128 = get_sin_k_rational(k.wrapping_add(64));

    // sin(x) = sin(k * pi/128 + u)
    //        = sin(u) * cos(k*pi/128) + cos(u) * sin(k*pi/128)
    let r = (sin_k_f128 * sin_cos.v_cos) + (cos_k_f128 * sin_cos.v_sin);
    r.fast_as_f64()
}

/// Cosine for double precision
///
/// ULP 0.5
#[inline]
pub fn f_cos(x: f64) -> f64 {
    let x_e = (x.to_bits() >> 52) & 0x7ff;
    const E_BIAS: u64 = (1u64 << (11 - 1u64)) - 1u64;

    let y: DoubleDouble;
    let k;

    let mut argument_reduction = LargeArgumentReduction::default();

    // |x| < 2^32 (with FMA) or |x| < 2^23 (w/o FMA)
    if x_e < E_BIAS + 16 {
        // |x| < 2^-26
        if x_e < E_BIAS - 7 {
            // |x| < 2^-26
            if x_e < E_BIAS - 27 {
                // Signed zeros.
                if x == 0.0 {
                    return 1.0;
                }
                // For |x| < 2^-26, |sin(x) - x| < ulp(x)/2.
                return 1.0 - min_normal_f64();
            }
            k = 0;
            y = DoubleDouble::new(0.0, x);
        } else {
            // // Small range reduction.
            (y, k) = range_reduction_small(x);
        }
    } else {
        // Inf or NaN
        if x_e > 2 * E_BIAS {
            // sin(+-Inf) = NaN
            return x + f64::NAN;
        }

        // Large range reduction.
        // k = argument_reduction.high_part(x);
        (k, y) = argument_reduction.reduce(x);
    }
    let r_sincos = sincos_eval(y);

    // Fast look up version, but needs 256-entry table.
    // cos(k * pi/128) = sin(k * pi/128 + pi/2) = sin((k + 64) * pi/128).
    let sk = SIN_K_PI_OVER_128[(k.wrapping_add(128) & 255) as usize];
    let ck = SIN_K_PI_OVER_128[((k.wrapping_add(64)) & 255) as usize];
    let msin_k = DoubleDouble::new(f64::from_bits(sk.0), f64::from_bits(sk.1));
    let cos_k = DoubleDouble::new(f64::from_bits(ck.0), f64::from_bits(ck.1));

    let sin_k_cos_y = DoubleDouble::quick_mult(r_sincos.v_cos, cos_k);
    let cos_k_sin_y = DoubleDouble::quick_mult(r_sincos.v_sin, msin_k);

    let mut rr = DoubleDouble::from_full_exact_add(sin_k_cos_y.hi, cos_k_sin_y.hi);
    rr.lo += sin_k_cos_y.lo + cos_k_sin_y.lo;
    let rlp = rr.lo + r_sincos.err;
    let rlm = rr.lo - r_sincos.err;

    let r_upper = rr.hi + rlp; // (rr.lo + ERR);
    let r_lower = rr.hi + rlm; // (rr.lo - ERR);

    // Ziv's accuracy test
    if r_upper == r_lower {
        return rr.to_f64();
    }

    const EXP_BIAS: u64 = (1u64 << (11 - 1u64)) - 1u64;
    let u_f128 = if x_e < EXP_BIAS + 16 {
        range_reduction_small_f128(x)
    } else {
        argument_reduction.accurate()
    };

    let sin_cos = sincos_eval_dyadic(&u_f128);

    // -sin(k * pi/128) = sin((k + 128) * pi/128)
    // cos(k * pi/128) = sin(k * pi/128 + pi/2) = sin((k + 64) * pi/128).
    let msin_k_f128 = get_sin_k_rational(k.wrapping_add(128));
    let cos_k_f128 = get_sin_k_rational(k.wrapping_add(64));

    // cos(x) = cos((k * pi/128 + u)
    //        = cos(u) * cos(k*pi/128) - sin(u) * sin(k*pi/128)
    let r = (cos_k_f128 * sin_cos.v_cos) + (msin_k_f128 * sin_cos.v_sin);
    r.fast_as_f64()
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn cos_test() {
        println!("{}", f_cos(4.750064341072175));
        assert_eq!(f_cos(0.0), 1.0);
        assert_eq!(f_cos(1.0), 0.5403023058681398);
        assert_eq!(f_cos(-0.5), 0.8775825618903728);
    }

    #[test]
    fn sin_test() {
        assert_eq!(f_sin(2.8477476437362989E-306), 2.8477476437362989E-306);
        assert_eq!(f_sin(0.0), 0.0);
        assert_eq!(f_sin(1.0), 0.8414709848078965);
        assert_eq!(f_sin(-0.5), -0.479425538604203);
    }
}
