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
use crate::rational_float::{RationalFloat128, RationalSign};
use crate::sin::ONE_TWENTY_EIGHT_OVER_PI;

#[inline]
pub(crate) fn range_reduction_small_f128(x: f64) -> RationalFloat128 {
    const PI_OVER_128_F128: RationalFloat128 = RationalFloat128 {
        sign: RationalSign::Pos,
        exponent: -133,
        mantissa: 0xc90f_daa2_2168_c234_c4c6_628b_80dc_1cd1_u128,
    };
    const ONE_TWENTY_EIGHT_OVER_PI_D: f64 = f64::from_bits(0x40445f306dc9c883);
    let prod_hi = x * ONE_TWENTY_EIGHT_OVER_PI_D;
    let kd = prod_hi.round();

    let mk_f128 = RationalFloat128::new_from_f64(-kd);
    let x_f128 = RationalFloat128::new_from_f64(x);
    let over_pi3 = ONE_TWENTY_EIGHT_OVER_PI[3];
    let p_hi = x_f128.quick_mul(&RationalFloat128::new_from_f64(f64::from_bits(over_pi3.0)));
    let p_mid = x_f128.quick_mul(&RationalFloat128::new_from_f64(f64::from_bits(over_pi3.1)));
    let p_lo = x_f128.quick_mul(&RationalFloat128::new_from_f64(f64::from_bits(over_pi3.2)));
    let s_hi = p_hi.quick_add(&mk_f128);
    let s_lo = p_mid.quick_add(&p_lo);
    let y = s_hi.quick_add(&s_lo);
    y.quick_mul(&PI_OVER_128_F128)
}

#[inline(always)]
pub(crate) fn r_fmla(
    a: &RationalFloat128,
    b: &RationalFloat128,
    c: &RationalFloat128,
) -> RationalFloat128 {
    a.quick_mul(b).quick_add(c)
}

#[inline(always)]
#[allow(clippy::too_many_arguments)]
pub(crate) fn r_polyeval7(
    x: &RationalFloat128,
    a0: &RationalFloat128,
    a1: &RationalFloat128,
    a2: &RationalFloat128,
    a3: &RationalFloat128,
    a4: &RationalFloat128,
    a5: &RationalFloat128,
    a6: &RationalFloat128,
) -> RationalFloat128 {
    let t1 = r_fmla(x, a6, a5); // a3 * x + a2
    let t2 = r_fmla(x, &t1, a4); // a3 * x + a2
    let t3 = r_fmla(x, &t2, a3); // (a3 * x + a2) * x + a1
    let t4 = r_fmla(x, &t3, a2); // (a3 * x + a2) * x + a1
    let t5 = r_fmla(x, &t4, a1); // (a3 * x + a2) * x + a1
    r_fmla(x, &t5, a0) // ((a3 * x + a2) * x + a1) * x + a0
}

#[inline(always)]
#[allow(clippy::too_many_arguments)]
pub(crate) fn r_polyeval9(
    x: &RationalFloat128,
    a0: &RationalFloat128,
    a1: &RationalFloat128,
    a2: &RationalFloat128,
    a3: &RationalFloat128,
    a4: &RationalFloat128,
    a5: &RationalFloat128,
    a6: &RationalFloat128,
    a7: &RationalFloat128,
    a8: &RationalFloat128,
) -> RationalFloat128 {
    let t0 = r_fmla(x, a8, a7); // a3 * x + a2
    let t0a = r_fmla(x, &t0, a6); // a3 * x + a2
    let t1 = r_fmla(x, &t0a, a5); // a3 * x + a2
    let t2 = r_fmla(x, &t1, a4); // a3 * x + a2
    let t3 = r_fmla(x, &t2, a3); // (a3 * x + a2) * x + a1
    let t4 = r_fmla(x, &t3, a2); // (a3 * x + a2) * x + a1
    let t5 = r_fmla(x, &t4, a1); // (a3 * x + a2) * x + a1
    r_fmla(x, &t5, a0) // ((a3 * x + a2) * x + a1) * x + a0
}

pub(crate) struct SinCosRational {
    pub(crate) v_sin: RationalFloat128,
    pub(crate) v_cos: RationalFloat128,
}

pub(crate) fn sincos_eval_rational(u: &RationalFloat128) -> SinCosRational {
    let u_sq = u.quick_mul(u);

    // sin(u) ~ x - x^3/3! + x^5/5! - x^7/7! + x^9/9! - x^11/11! + x^13/13!
    const SIN_COEFFS: [RationalFloat128; 7] = [
        RationalFloat128 {
            sign: RationalSign::Pos,
            exponent: -127,
            mantissa: 0x80000000_00000000_00000000_00000000_u128,
        }, // 1
        RationalFloat128 {
            sign: RationalSign::Neg,
            exponent: -130,
            mantissa: 0xaaaaaaaa_aaaaaaaa_aaaaaaaa_aaaaaaab_u128,
        }, // -1/3!
        RationalFloat128 {
            sign: RationalSign::Pos,
            exponent: -134,
            mantissa: 0x88888888_88888888_88888888_88888889_u128,
        }, // 1/5!
        RationalFloat128 {
            sign: RationalSign::Neg,
            exponent: -140,
            mantissa: 0xd00d00d0_0d00d00d_00d00d00_d00d00d0_u128,
        }, // -1/7!
        RationalFloat128 {
            sign: RationalSign::Pos,
            exponent: -146,
            mantissa: 0xb8ef1d2a_b6399c7d_560e4472_800b8ef2_u128,
        }, // 1/9!
        RationalFloat128 {
            sign: RationalSign::Neg,
            exponent: -153,
            mantissa: 0xd7322b3f_aa271c7f_3a3f25c1_bee38f10_u128,
        }, // -1/11!
        RationalFloat128 {
            sign: RationalSign::Pos,
            exponent: -160,
            mantissa: 0xb092309d_43684be5_1c198e91_d7b4269e_u128,
        }, // 1/13!
    ];

    // cos(u) ~ 1 - x^2/2 + x^4/4! - x^6/6! + x^8/8! - x^10/10! + x^12/12!
    const COS_COEFFS: [RationalFloat128; 7] = [
        RationalFloat128 {
            sign: RationalSign::Pos,
            exponent: -127,
            mantissa: 0x80000000_00000000_00000000_00000000_u128,
        }, // 1.0
        RationalFloat128 {
            sign: RationalSign::Neg,
            exponent: -128,
            mantissa: 0x80000000_00000000_00000000_00000000_u128,
        }, // 1/2
        RationalFloat128 {
            sign: RationalSign::Pos,
            exponent: -132,
            mantissa: 0xaaaaaaaa_aaaaaaaa_aaaaaaaa_aaaaaaab_u128,
        }, // 1/4!
        RationalFloat128 {
            sign: RationalSign::Neg,
            exponent: -137,
            mantissa: 0xb60b60b6_0b60b60b_60b60b60_b60b60b6_u128,
        }, // 1/6!
        RationalFloat128 {
            sign: RationalSign::Pos,
            exponent: -143,
            mantissa: 0xd00d00d0_0d00d00d_00d00d00_d00d00d0_u128,
        }, // 1/8!
        RationalFloat128 {
            sign: RationalSign::Neg,
            exponent: -149,
            mantissa: 0x93f27dbb_c4fae397_780b69f5_333c725b_u128,
        }, // 1/10!
        RationalFloat128 {
            sign: RationalSign::Pos,
            exponent: -156,
            mantissa: 0x8f76c77f_c6c4bdaa_26d4c3d6_7f425f60_u128,
        }, // 1/12!
    ];

    let sin_u = u.quick_mul(&r_polyeval7(
        &u_sq,
        &SIN_COEFFS[0],
        &SIN_COEFFS[1],
        &SIN_COEFFS[2],
        &SIN_COEFFS[3],
        &SIN_COEFFS[4],
        &SIN_COEFFS[5],
        &SIN_COEFFS[6],
    ));
    let cos_u = r_polyeval7(
        &u_sq,
        &COS_COEFFS[0],
        &COS_COEFFS[1],
        &COS_COEFFS[2],
        &COS_COEFFS[3],
        &COS_COEFFS[4],
        &COS_COEFFS[5],
        &COS_COEFFS[6],
    );
    SinCosRational {
        v_sin: sin_u,
        v_cos: cos_u,
    }
}

pub(crate) static SIN_K_PI_OVER_128_F128: [RationalFloat128; 65] = [
    RationalFloat128 {
        sign: RationalSign::Pos,
        exponent: 0,
        mantissa: 0,
    },
    RationalFloat128 {
        sign: RationalSign::Pos,
        exponent: -133,
        mantissa: 0xc90a_afbd_1b33_efc9_c539_edcb_fda0_cf2c_u128,
    },
    RationalFloat128 {
        sign: RationalSign::Pos,
        exponent: -132,
        mantissa: 0xc8fb_2f88_6ec0_9f37_6a17_954b_2b7c_5171_u128,
    },
    RationalFloat128 {
        sign: RationalSign::Pos,
        exponent: -131,
        mantissa: 0x96a9_0496_70cf_ae65_f775_7409_4d3c_35c4_u128,
    },
    RationalFloat128 {
        sign: RationalSign::Pos,
        exponent: -131,
        mantissa: 0xc8bd_35e1_4da1_5f0e_c739_6c89_4bbf_7389_u128,
    },
    RationalFloat128 {
        sign: RationalSign::Pos,
        exponent: -131,
        mantissa: 0xfab2_72b5_4b98_71a2_7047_29ae_56d7_8a37_u128,
    },
    RationalFloat128 {
        sign: RationalSign::Pos,
        exponent: -130,
        mantissa: 0x9640_8374_7309_d113_000a_89a1_1e07_c1fe_u128,
    },
    RationalFloat128 {
        sign: RationalSign::Pos,
        exponent: -130,
        mantissa: 0xaf10_a224_59fe_32a6_3fee_f3bb_58b1_f10d_u128,
    },
    RationalFloat128 {
        sign: RationalSign::Pos,
        exponent: -130,
        mantissa: 0xc7c5_c1e3_4d30_55b2_5cc8_c00e_4fcc_d850_u128,
    },
    RationalFloat128 {
        sign: RationalSign::Pos,
        exponent: -130,
        mantissa: 0xe05c_1353_f27b_17e5_0ebc_61ad_e6ca_83cd_u128,
    },
    RationalFloat128 {
        sign: RationalSign::Pos,
        exponent: -130,
        mantissa: 0xf8cf_cbd9_0af8_d57a_4221_dc4b_a772_598d_u128,
    },
    RationalFloat128 {
        sign: RationalSign::Pos,
        exponent: -129,
        mantissa: 0x888e_9315_8fb3_bb04_9841_56f5_5334_4306_u128,
    },
    RationalFloat128 {
        sign: RationalSign::Pos,
        exponent: -129,
        mantissa: 0x94a0_3176_acf8_2d45_ae4b_a773_da6b_f754_u128,
    },
    RationalFloat128 {
        sign: RationalSign::Pos,
        exponent: -129,
        mantissa: 0xa09a_e4a0_bb30_0a19_2f89_5f44_a303_cc0b_u128,
    },
    RationalFloat128 {
        sign: RationalSign::Pos,
        exponent: -129,
        mantissa: 0xac7c_d3ad_58fe_e7f0_811f_9539_84ef_f83e_u128,
    },
    RationalFloat128 {
        sign: RationalSign::Pos,
        exponent: -129,
        mantissa: 0xb844_2987_d22c_f576_9cc3_ef36_746d_e3b8_u128,
    },
    RationalFloat128 {
        sign: RationalSign::Pos,
        exponent: -129,
        mantissa: 0xc3ef_1535_754b_168d_3122_c2a5_9efd_dc37_u128,
    },
    RationalFloat128 {
        sign: RationalSign::Pos,
        exponent: -129,
        mantissa: 0xcf7b_ca1d_476c_516d_a812_90bd_baad_62e4_u128,
    },
    RationalFloat128 {
        sign: RationalSign::Pos,
        exponent: -129,
        mantissa: 0xdae8_804f_0ae6_015b_362c_b974_182e_3030_u128,
    },
    RationalFloat128 {
        sign: RationalSign::Pos,
        exponent: -129,
        mantissa: 0xe633_74c9_8e22_f0b4_2872_ce1b_fc7a_d1cd_u128,
    },
    RationalFloat128 {
        sign: RationalSign::Pos,
        exponent: -129,
        mantissa: 0xf15a_e9c0_37b1_d8f0_6c48_e9e3_420b_0f1e_u128,
    },
    RationalFloat128 {
        sign: RationalSign::Pos,
        exponent: -129,
        mantissa: 0xfc5d_26df_c4d5_cfda_27c0_7c91_1290_b8d1_u128,
    },
    RationalFloat128 {
        sign: RationalSign::Pos,
        exponent: -128,
        mantissa: 0x839c_3cc9_17ff_6cb4_bfd7_9717_f288_0abf_u128,
    },
    RationalFloat128 {
        sign: RationalSign::Pos,
        exponent: -128,
        mantissa: 0x88f5_9aa0_da59_1421_b892_ca83_61d8_c84c_u128,
    },
    RationalFloat128 {
        sign: RationalSign::Pos,
        exponent: -128,
        mantissa: 0x8e39_d9cd_7346_4364_bba4_cfec_bff5_4867_u128,
    },
    RationalFloat128 {
        sign: RationalSign::Pos,
        exponent: -128,
        mantissa: 0x9368_2a66_e896_f544_b178_2191_1e71_c16e_u128,
    },
    RationalFloat128 {
        sign: RationalSign::Pos,
        exponent: -128,
        mantissa: 0x987f_bfe7_0b81_a708_19ce_c845_ac87_a5c6_u128,
    },
    RationalFloat128 {
        sign: RationalSign::Pos,
        exponent: -128,
        mantissa: 0x9d7f_d149_0285_c9e3_e25e_3954_9638_ae68_u128,
    },
    RationalFloat128 {
        sign: RationalSign::Pos,
        exponent: -128,
        mantissa: 0xa267_9928_48ee_b0c0_3b51_67ee_359a_234e_u128,
    },
    RationalFloat128 {
        sign: RationalSign::Pos,
        exponent: -128,
        mantissa: 0xa736_55df_1f2f_489e_149f_6e75_9934_68a3_u128,
    },
    RationalFloat128 {
        sign: RationalSign::Pos,
        exponent: -128,
        mantissa: 0xabeb_49a4_6764_fd15_1bec_da80_89c1_a94c_u128,
    },
    RationalFloat128 {
        sign: RationalSign::Pos,
        exponent: -128,
        mantissa: 0xb085_baa8_e966_f6da_e4ca_d00d_5c94_bcd2_u128,
    },
    RationalFloat128 {
        sign: RationalSign::Pos,
        exponent: -128,
        mantissa: 0xb504_f333_f9de_6484_597d_89b3_754a_be9f_u128,
    },
    RationalFloat128 {
        sign: RationalSign::Pos,
        exponent: -128,
        mantissa: 0xb968_41bf_7ffc_b21a_9de1_e3b2_2b8b_f4db_u128,
    },
    RationalFloat128 {
        sign: RationalSign::Pos,
        exponent: -128,
        mantissa: 0xbdae_f913_557d_76f0_ac85_320f_528d_6d5d_u128,
    },
    RationalFloat128 {
        sign: RationalSign::Pos,
        exponent: -128,
        mantissa: 0xc1d8_705f_fcbb_6e90_bdf0_715c_b8b2_0bd7_u128,
    },
    RationalFloat128 {
        sign: RationalSign::Pos,
        exponent: -128,
        mantissa: 0xc5e4_0358_a8ba_05a7_43da_25d9_9267_326b_u128,
    },
    RationalFloat128 {
        sign: RationalSign::Pos,
        exponent: -128,
        mantissa: 0xc9d1_124c_931f_da7a_8335_241b_e169_3225_u128,
    },
    RationalFloat128 {
        sign: RationalSign::Pos,
        exponent: -128,
        mantissa: 0xcd9f_023f_9c3a_059e_23af_31db_7179_a4aa_u128,
    },
    RationalFloat128 {
        sign: RationalSign::Pos,
        exponent: -128,
        mantissa: 0xd14d_3d02_313c_0eed_744f_ea20_e8ab_ef92_u128,
    },
    RationalFloat128 {
        sign: RationalSign::Pos,
        exponent: -128,
        mantissa: 0xd4db_3148_750d_1819_f630_e8b6_dac8_3e69_u128,
    },
    RationalFloat128 {
        sign: RationalSign::Pos,
        exponent: -128,
        mantissa: 0xd848_52c0_a80f_fcdb_24b9_fe00_6635_74a4_u128,
    },
    RationalFloat128 {
        sign: RationalSign::Pos,
        exponent: -128,
        mantissa: 0xdb94_1a28_cb71_ec87_2c19_b632_53da_43fc_u128,
    },
    RationalFloat128 {
        sign: RationalSign::Pos,
        exponent: -128,
        mantissa: 0xdebe_0563_7ca9_4cfb_4b19_aa71_fec3_ae6d_u128,
    },
    RationalFloat128 {
        sign: RationalSign::Pos,
        exponent: -128,
        mantissa: 0xe1c5_978c_05ed_8691_f4e8_a837_2f8c_5810_u128,
    },
    RationalFloat128 {
        sign: RationalSign::Pos,
        exponent: -128,
        mantissa: 0xe4aa_5909_a08f_a7b4_1227_85ae_67f5_515d_u128,
    },
    RationalFloat128 {
        sign: RationalSign::Pos,
        exponent: -128,
        mantissa: 0xe76b_d7a1_e63b_9786_1251_2952_9d48_a92f_u128,
    },
    RationalFloat128 {
        sign: RationalSign::Pos,
        exponent: -128,
        mantissa: 0xea09_a68a_6e49_cd62_15ad_45b4_a1b5_e823_u128,
    },
    RationalFloat128 {
        sign: RationalSign::Pos,
        exponent: -128,
        mantissa: 0xec83_5e79_946a_3145_7e61_0231_ac1d_6181_u128,
    },
    RationalFloat128 {
        sign: RationalSign::Pos,
        exponent: -128,
        mantissa: 0xeed8_9db6_6611_e307_86f8_c20f_b664_b01b_u128,
    },
    RationalFloat128 {
        sign: RationalSign::Pos,
        exponent: -128,
        mantissa: 0xf109_0827_b437_25fd_6712_7db3_5b28_7316_u128,
    },
    RationalFloat128 {
        sign: RationalSign::Pos,
        exponent: -128,
        mantissa: 0xf314_4762_4708_8f74_a548_6bdc_455d_56a2_u128,
    },
    RationalFloat128 {
        sign: RationalSign::Pos,
        exponent: -128,
        mantissa: 0xf4fa_0ab6_316e_d2ec_163c_5c7f_03b7_18c5_u128,
    },
    RationalFloat128 {
        sign: RationalSign::Pos,
        exponent: -128,
        mantissa: 0xf6ba_073b_424b_19e8_2c79_1f59_cc1f_fc23_u128,
    },
    RationalFloat128 {
        sign: RationalSign::Pos,
        exponent: -128,
        mantissa: 0xf853_f7dc_9186_b952_c7ad_c6b4_9888_91bb_u128,
    },
    RationalFloat128 {
        sign: RationalSign::Pos,
        exponent: -128,
        mantissa: 0xf9c7_9d63_272c_4628_4504_ae08_d19b_2980_u128,
    },
    RationalFloat128 {
        sign: RationalSign::Pos,
        exponent: -128,
        mantissa: 0xfb14_be7f_bae5_8156_2172_a361_fd2a_722f_u128,
    },
    RationalFloat128 {
        sign: RationalSign::Pos,
        exponent: -128,
        mantissa: 0xfc3b_27d3_8a5d_49ab_2567_78ff_cb5c_1769_u128,
    },
    RationalFloat128 {
        sign: RationalSign::Pos,
        exponent: -128,
        mantissa: 0xfd3a_abf8_4528_b50b_eae6_bd95_1c1d_abbe_u128,
    },
    RationalFloat128 {
        sign: RationalSign::Pos,
        exponent: -128,
        mantissa: 0xfe13_2387_0cfe_9a3d_90cd_1d95_9db6_74ef_u128,
    },
    RationalFloat128 {
        sign: RationalSign::Pos,
        exponent: -128,
        mantissa: 0xfec4_6d1e_8929_2cf0_4139_0efd_c726_e9ef_u128,
    },
    RationalFloat128 {
        sign: RationalSign::Pos,
        exponent: -128,
        mantissa: 0xff4e_6d68_0c41_d0a9_0f66_8633_f1ab_858a_u128,
    },
    RationalFloat128 {
        sign: RationalSign::Pos,
        exponent: -128,
        mantissa: 0xffb1_0f1b_cb6b_ef1d_421e_8eda_af59_453e_u128,
    },
    RationalFloat128 {
        sign: RationalSign::Pos,
        exponent: -128,
        mantissa: 0xffec_4304_2668_65d9_5657_5523_6696_1732_u128,
    },
    RationalFloat128 {
        sign: RationalSign::Pos,
        exponent: -127,
        mantissa: 0x8000_0000_0000_0000_0000_0000_0000_0000_u128,
    },
];
