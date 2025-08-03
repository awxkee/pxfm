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
use crate::dyadic_float::{DyadicFloat128, DyadicSign};

/**
Remez poly at zero and extremums for Y0

```python
def compute_intervals(zeros):
    intervals = []
    for i in range(0, len(zeros)):
        if i == 0:
            a = 2 - zeros[i]
            b = (zeros[i] + zeros[i + 1]) / 2 + 0.03 - zeros[i]
            intervals.append((RealField(18)(a), RealField(18)(b), RealField(110)(zeros[i])))
        elif i + 1 > len(zeros) - 1:
            a = (zeros[i - 1] + zeros[i]) / 2 - 0.03 - zeros[i]
            b = (zeros[i]) + 0.83 + 0.03 - zeros[i]
            intervals.append((RealField(18)(a), RealField(18)(b), RealField(110)(zeros[i])))
        else:
            a = (zeros[i - 1] + zeros[i]) / 2 - zeros[i] - 0.03
            b = (zeros[i] + zeros[i + 1]) / 2 + 0.03  - zeros[i]
            intervals.append((RealField(18)(a), RealField(18)(b), RealField(110)(zeros[i])))
    return intervals

intervals = compute_intervals(y0_zeros)
# print(intervals)

def build_sollya_script(a, b, zero, deg):
    return f"""
prec = 200;
bessel_y0 = library("/Users/radzivon/RustroverProjects/pxfm/notes/bessel_sollya/cmake-build-release/libbessel_sollya.dylib");
f = bessel_y0(x + {zero});
d = [{a}, {b}];
pf = remez(f, {deg}, d);
for i from 0 to degree(pf) do {{
    write(coeff(pf, i)) >> "coefficients.txt";
    write("\\n") >> "coefficients.txt";
}};
"""

def load_coefficients(filename):
    with open(filename, "r") as f:
        return [RealField(500)(line.strip()) for line in f if line.strip()]

def call_sollya_on_interval(a, b, zero, degree=12):
    sollya_script = build_sollya_script(a, b, zero, degree)
    with open("tmp_interval.sollya", "w") as f:
        f.write(sollya_script)
    import subprocess
    if os.path.exists("coefficients.txt"):
        os.remove("coefficients.txt")
    try:
        result = subprocess.run(
            ["sollya", "tmp_interval.sollya"],
            check=True,
            capture_output=True,
            text=True
        )
    except subprocess.CalledProcessError as e:
        return

def print_remez_coeffs(poly):
    print("[")
    for i in range(len(poly)):
        coeff = poly[i]
        print_dyadic(coeff)
    print("],")

degree = 27

print(f"pub(crate) static Y0_COEFFS_RATIONAL128_REMEZ: [[DyadicFloat128; {degree + 1}]; {len(intervals)}] = [")
for i in range(0, len(intervals)):
    interval = intervals[i]
    call_sollya_on_interval(interval[0], interval[1], interval[2], degree)
    coeffs = load_coefficients(f"coefficients.txt")
    print_remez_coeffs(coeffs)
print("];")
```
**/
pub(crate) static Y0_COEFFS_RATIONAL128_REMEZ: [[DyadicFloat128; 28]; 47] = [
    [
        DyadicFloat128 {
            sign: DyadicSign::Pos,
            exponent: -128,
            mantissa: 0x85524221_780a56b6_d3e528f8_032b4a4b_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Pos,
            exponent: -192,
            mantissa: 0xe6e66301_cc8c4339_7cfec8a9_fbebccff_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Neg,
            exponent: -129,
            mantissa: 0x85524221_780a56df_e6dc5a69_44ff6a66_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Pos,
            exponent: -132,
            mantissa: 0xa1cfd60b_292fa8a5_82712abe_f91207eb_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Pos,
            exponent: -134,
            mantissa: 0x86957a74_991ddf53_24e1e44d_76746ba7_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Pos,
            exponent: -138,
            mantissa: 0xfb8b235f_547c5631_38e7fab0_a7b6c7b8_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Neg,
            exponent: -137,
            mantissa: 0xa225ed45_b4ea6cf5_39366260_e9078af6_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Pos,
            exponent: -139,
            mantissa: 0xc2611033_80dfeb43_9717caf9_91d4b2e7_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Neg,
            exponent: -140,
            mantissa: 0x8bd5a57c_685f8035_c83a6db9_9d1a8e40_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Pos,
            exponent: -142,
            mantissa: 0xed7dcdcc_d3878cc6_9aa7d55b_7b515587_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Neg,
            exponent: -143,
            mantissa: 0xc5dcb1ce_b0103aa3_df9b2e03_7aaed252_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Pos,
            exponent: -144,
            mantissa: 0xa45b3a8c_13c0cb42_1b0e2525_4ee581dc_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Neg,
            exponent: -145,
            mantissa: 0x89b2662c_43fc6a12_ad376312_c2ca56ac_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Pos,
            exponent: -147,
            mantissa: 0xe8357e27_4bb44852_2f66ae9f_4cc9884d_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Neg,
            exponent: -148,
            mantissa: 0xc4cdbaa9_8ea2d56f_1618297a_b0544812_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Pos,
            exponent: -149,
            mantissa: 0xa78aba33_444d9fb2_48db3b13_d2376e94_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Neg,
            exponent: -150,
            mantissa: 0x8f29beab_f2c3915f_b0fa7789_258a1a6d_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Pos,
            exponent: -152,
            mantissa: 0xf5e501f3_a7c38f6e_c223010f_4abcd554_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Neg,
            exponent: -153,
            mantissa: 0xd4adfe38_bbef8edc_448397ed_0c23d22c_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Pos,
            exponent: -154,
            mantissa: 0xb5cdc5d5_153f3324_ad252ca0_ab1f4f42_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Neg,
            exponent: -155,
            mantissa: 0x96a66cef_866e16fd_1bef2282_1d4f5469_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Pos,
            exponent: -156,
            mantissa: 0x8b4e8271_e3c3328f_e9b24d67_2cdad1e8_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Neg,
            exponent: -157,
            mantissa: 0x96d623ff_ee589521_27695d85_f76d15e4_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Pos,
            exponent: -159,
            mantissa: 0xd0b5ece9_22a7f87e_a333265a_2e99211d_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Neg,
            exponent: -164,
            mantissa: 0x8c7abe80_6686405c_59937860_e3b9a46e_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Pos,
            exponent: -162,
            mantissa: 0xd6189a58_a7aa626c_007c5b32_e78df314_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Neg,
            exponent: -161,
            mantissa: 0xc45fe49b_69d22f1a_6794ae0c_6146f4eb_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Pos,
            exponent: -162,
            mantissa: 0x83a515fe_c9ad3a01_870c6a65_48ea07eb_u128,
        },
    ],
    [
        DyadicFloat128 {
            sign: DyadicSign::Neg,
            exponent: -222,
            mantissa: 0xe7235dcb_f1d9b985_ea12b68d_bc657cb6_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Neg,
            exponent: -129,
            mantissa: 0xce1a12b5_095060d2_01372240_3ad5be14_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Pos,
            exponent: -132,
            mantissa: 0xd04e4948_1b3f74c4_25237a00_9ffc60e0_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Pos,
            exponent: -132,
            mantissa: 0xefb6acdf_a875e756_27ec9c5d_4664e0e5_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Neg,
            exponent: -135,
            mantissa: 0xe08b7ee2_cc04b31b_7f8f7289_cbfe61d4_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Neg,
            exponent: -136,
            mantissa: 0x8f195e27_7c5206bf_3066a0b0_ef0526be_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Pos,
            exponent: -140,
            mantissa: 0xccc13b28_b3129c8f_a4c39e08_8ebd01ad_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Pos,
            exponent: -142,
            mantissa: 0xd5960ff6_6e7e1f97_cfed28a4_340b8347_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Neg,
            exponent: -145,
            mantissa: 0xa431b8f3_5331fc64_be7d40a8_7250c4f8_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Neg,
            exponent: -149,
            mantissa: 0x9d90b420_44dcc05e_5496c261_e936d4f8_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Neg,
            exponent: -153,
            mantissa: 0xbbb82c3a_0d580277_4e61e388_1d1e1c9e_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Pos,
            exponent: -153,
            mantissa: 0xa2975004_a447ed4f_84281042_e6e40f71_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Neg,
            exponent: -156,
            mantissa: 0xfb49c414_d7566667_1d671993_14a3093e_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Pos,
            exponent: -158,
            mantissa: 0xd2d0f9d0_1e9aeb3a_db8d2ede_792a373a_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Neg,
            exponent: -160,
            mantissa: 0xcca739ac_676fae66_b5d2f0a9_5fcee799_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Pos,
            exponent: -162,
            mantissa: 0xc4572e98_d7ae09a3_54bba4fc_370f782a_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Neg,
            exponent: -164,
            mantissa: 0xba994bc3_f9ec6157_cd4a8a94_07347a6e_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Pos,
            exponent: -166,
            mantissa: 0xb21e2bc7_8426745e_f8c75c2b_cf4c8af2_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Neg,
            exponent: -168,
            mantissa: 0xaaa3c069_d5d0d2a4_d44d2d80_2b7c327f_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Pos,
            exponent: -170,
            mantissa: 0xa3e34d1b_493f446c_926be7b9_2e303f51_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Neg,
            exponent: -172,
            mantissa: 0x9dc47053_1379a117_5a5c0e46_061ac5d8_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Pos,
            exponent: -174,
            mantissa: 0x97ef7416_653d414f_94bb00b0_bd01fd47_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Neg,
            exponent: -176,
            mantissa: 0x92cc1f58_9cd75bd5_8c355097_5f9fdde7_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Pos,
            exponent: -178,
            mantissa: 0x922066c4_245c34f4_52404830_80a68731_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Neg,
            exponent: -180,
            mantissa: 0x8b26343b_e2b86768_2eb95b41_1b78c4dc_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Pos,
            exponent: -183,
            mantissa: 0xc504153e_e1cc8d12_fa0954d7_41228b06_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Neg,
            exponent: -184,
            mantissa: 0x90f50328_e90cc3e8_0f4a0327_d643f5fb_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Pos,
            exponent: -185,
            mantissa: 0x91acb5cc_f59ff7aa_476f7bae_5e830119_u128,
        },
    ],
    [
        DyadicFloat128 {
            sign: DyadicSign::Neg,
            exponent: -129,
            mantissa: 0xae3e2ab7_860ccc8e_5abc76d2_ff094304_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Pos,
            exponent: -232,
            mantissa: 0xf22a25c3_8470cad1_718741d7_2e3d5d98_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Pos,
            exponent: -130,
            mantissa: 0xae3e2ab7_860ccc8e_5abc76d2_f78280e4_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Neg,
            exponent: -134,
            mantissa: 0xab26a58f_686bf369_a89555e1_01a9cb75_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Neg,
            exponent: -134,
            mantissa: 0xd0aec96f_f1f13437_aec78a36_76413c65_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Pos,
            exponent: -138,
            mantissa: 0xda1c2182_34ffb448_2adba1bb_d8a343c9_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Pos,
            exponent: -139,
            mantissa: 0xbd449266_17c8a382_ec090e29_821dfc90_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Neg,
            exponent: -143,
            mantissa: 0xafb4da54_51e02881_0edfa10c_d91bf259_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Neg,
            exponent: -145,
            mantissa: 0xc2dca075_b03fcae5_2c8f4f9e_365e90a9_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Pos,
            exponent: -149,
            mantissa: 0xaad607ae_18b163b3_3b111e91_e02acff6_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Pos,
            exponent: -152,
            mantissa: 0xeab1f0c5_a3f5c3e4_b450e891_f971c5e5_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Neg,
            exponent: -156,
            mantissa: 0x9f561d22_f0379d75_acf800f8_0369be6e_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Neg,
            exponent: -158,
            mantissa: 0x803b38ea_ad457aab_497d164a_bdddfc98_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Pos,
            exponent: -162,
            mantissa: 0x82709f46_96b1dad0_382e3b16_f570d3c6_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Neg,
            exponent: -170,
            mantissa: 0xb7159836_8c7f10d9_8200cee4_ad451fc1_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Pos,
            exponent: -169,
            mantissa: 0xa68c7be4_047acea9_66782c39_e4466d43_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Neg,
            exponent: -171,
            mantissa: 0xa4c044ac_003fa2cb_7e4c95e9_954d8c62_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Pos,
            exponent: -174,
            mantissa: 0xd2da1d18_b40113e2_cb1d2b23_f1647e43_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Neg,
            exponent: -176,
            mantissa: 0x8e4e26ae_3f173b72_9e326d19_5dd5301b_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Pos,
            exponent: -179,
            mantissa: 0xc99e565d_bb1777c4_d63adeec_930255d2_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Neg,
            exponent: -181,
            mantissa: 0x8e335af2_ce087235_4179f3c7_0bc11381_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Pos,
            exponent: -184,
            mantissa: 0xc8298310_ffdb6505_75a60345_dc484012_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Neg,
            exponent: -186,
            mantissa: 0x8d38a158_c4766da0_6d62c3a3_53fd90d2_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Pos,
            exponent: -189,
            mantissa: 0xc7534da4_fd0760ff_41d6e4f6_04be1853_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Neg,
            exponent: -191,
            mantissa: 0x8c35b13b_9e5c5791_425823bc_3e63bb79_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Pos,
            exponent: -194,
            mantissa: 0xcbf71803_30eed7bd_51f00563_9cced75e_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Neg,
            exponent: -196,
            mantissa: 0xa101e963_8ef71a3d_3953de7e_4e50410f_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Pos,
            exponent: -199,
            mantissa: 0xbbce5b34_7a798290_fdaf698b_82393e9e_u128,
        },
    ],
    [
        DyadicFloat128 {
            sign: DyadicSign::Pos,
            exponent: -233,
            mantissa: 0x9420ccaf_22e24227_de4b69db_60f43030_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Pos,
            exponent: -129,
            mantissa: 0x99a66503_4bd2d4e5_76f69065_ab201a4a_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Neg,
            exponent: -133,
            mantissa: 0xad77b08f_e26aba83_8cd90da7_dc3833b8_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Neg,
            exponent: -132,
            mantissa: 0xc4b4e326_5fa29292_53687556_fc6f0aca_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Pos,
            exponent: -136,
            mantissa: 0xd978a54a_a93e59ab_353487eb_8be29cce_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Pos,
            exponent: -136,
            mantissa: 0x8e9af42f_ef151b54_cf4c5f51_0810aaa4_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Neg,
            exponent: -140,
            mantissa: 0x936eb8f1_c8e45927_0ba54528_c9445780_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Neg,
            exponent: -142,
            mantissa: 0xc0bbf27f_29192023_0c75cfb4_f42c88f7_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Pos,
            exponent: -146,
            mantissa: 0xb549139a_97cb50ca_4aa96c9a_a6e122eb_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Pos,
            exponent: -148,
            mantissa: 0x9a553ab9_f78a6d6b_48cc60d5_74d57343_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Neg,
            exponent: -152,
            mantissa: 0x85140a6b_d62242d3_16f7086c_55cfb10c_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Neg,
            exponent: -155,
            mantissa: 0xa2151e99_acd6cb20_177298ad_721008fe_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Pos,
            exponent: -160,
            mantissa: 0xfd36284a_b3639358_d9504254_90ef69f9_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Pos,
            exponent: -163,
            mantissa: 0xf8a72ac2_411ebad4_d95c9a2f_44ab1632_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Neg,
            exponent: -167,
            mantissa: 0xb8a61495_f1cb10ba_0801a442_ab1d2466_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Neg,
            exponent: -170,
            mantissa: 0x833ce6ac_2a4f0ee4_3396a0b0_dceda0ca_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Pos,
            exponent: -175,
            mantissa: 0x99b2ead3_c9089ae4_282d2213_aa45d444_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Pos,
            exponent: -178,
            mantissa: 0x9c678359_89a696b6_832c99f1_6f47a0a0_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Neg,
            exponent: -182,
            mantissa: 0x95ead425_f9b22358_c6ecd86a_c309ba06_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Pos,
            exponent: -187,
            mantissa: 0x9374996a_26cabba2_1660d31a_cd4b93a7_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Neg,
            exponent: -189,
            mantissa: 0x9388a823_29896be5_01eb0c7b_9b86861d_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Pos,
            exponent: -192,
            mantissa: 0xc6f9f30e_d14f2670_90b40cca_2f2cf773_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Neg,
            exponent: -195,
            mantissa: 0xcd8931e4_30ab9e5d_ab1c5e56_68fd4582_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Pos,
            exponent: -198,
            mantissa: 0xd9f7a582_977349db_06e04e97_e2767870_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Neg,
            exponent: -201,
            mantissa: 0xed3d75b9_568d9258_74788d2e_95989014_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Pos,
            exponent: -204,
            mantissa: 0xfedde63c_3d7fcd6f_dea0671a_2e517fcf_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Neg,
            exponent: -206,
            mantissa: 0x98943818_fcf66d92_41028889_626f0c73_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Pos,
            exponent: -209,
            mantissa: 0xc2e729a0_bd1dff12_45f87c17_3ccbd690_u128,
        },
    ],
    [
        DyadicFloat128 {
            sign: DyadicSign::Pos,
            exponent: -129,
            mantissa: 0x8afcc9fe_755ae11d_c672a53c_58fddb22_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Neg,
            exponent: -234,
            mantissa: 0xb540ec57_60a8a4af_0e8c0154_14871adf_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Neg,
            exponent: -130,
            mantissa: 0x8afcc9fe_755ae11d_c672a53c_58f356b6_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Pos,
            exponent: -135,
            mantissa: 0xac77b77d_f6bcb62e_3c6fe763_509796fa_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Pos,
            exponent: -134,
            mantissa: 0xb1caeff2_4fe69d78_2433690f_a4e87e38_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Neg,
            exponent: -139,
            mantissa: 0xfd8a8825_20e000da_2422c161_0bde7679_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Neg,
            exponent: -139,
            mantissa: 0xafc4508e_aab27a88_21cea57d_6f0fcdd4_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Pos,
            exponent: -144,
            mantissa: 0xf9be9135_4125b49a_4c8f8b49_760264ec_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Pos,
            exponent: -145,
            mantissa: 0xb7bdd588_279aa6e0_4e87dc86_8cba7d30_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Neg,
            exponent: -150,
            mantissa: 0xf8622d1c_126b7c4b_a68f188a_d44f1eed_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Neg,
            exponent: -152,
            mantissa: 0xeff08ef8_9638a822_26710037_d82f66a4_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Pos,
            exponent: -156,
            mantissa: 0x988d6e17_3a9f6d5a_d8faf2a8_9258f66a_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Pos,
            exponent: -159,
            mantissa: 0xd69a70c2_8270adbf_f3da8668_0fb0b50f_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Neg,
            exponent: -164,
            mantissa: 0xfee4030d_7458ab8d_d6b73f5b_c5aab0a9_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Neg,
            exponent: -166,
            mantissa: 0x8ca4c653_677cdaa2_e8e8c92c_5ef5e50d_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Pos,
            exponent: -171,
            mantissa: 0x9d148eac_c8fb5b3f_3df98a9e_9bc2552f_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Pos,
            exponent: -174,
            mantissa: 0x8b7bd740_656ffbe7_f26cf332_de625e8c_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Neg,
            exponent: -179,
            mantissa: 0x8fc6c06f_65b2c922_e176854f_c2b5b6c2_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Neg,
            exponent: -183,
            mantissa: 0xe13da9a8_6be7387d_997cec23_e390ca24_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Pos,
            exponent: -188,
            mantissa: 0xe8401a70_706a60fc_0d6e34e4_56fe776b_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Pos,
            exponent: -192,
            mantissa: 0xfd610091_06b3bc25_9c895ad2_189c2898_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Neg,
            exponent: -197,
            mantissa: 0x9f492bdd_df648fd2_46a5113e_d958296e_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Neg,
            exponent: -200,
            mantissa: 0xee70435e_caa75ed4_529079f6_052b3d90_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Pos,
            exponent: -204,
            mantissa: 0xe22fa277_56348a78_8cd48c55_2c64ab93_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Neg,
            exponent: -208,
            mantissa: 0xdfc16a29_fac36e53_5e361a15_2dd9c344_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Pos,
            exponent: -211,
            mantissa: 0xf780d229_7312183a_b7510e2c_ea21c3b0_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Neg,
            exponent: -214,
            mantissa: 0xfd5d0d64_7e0a0823_50317c2c_710a4086_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Pos,
            exponent: -217,
            mantissa: 0xcddd7cc8_b7f6bd15_716a1178_044a6209_u128,
        },
    ],
    [
        DyadicFloat128 {
            sign: DyadicSign::Pos,
            exponent: -232,
            mantissa: 0xf29cf725_f4cc5bed_f21d0af0_76fc147a_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Neg,
            exponent: -130,
            mantissa: 0xffb1ae63_95cf85c2_e152e0e7_29373036_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Pos,
            exponent: -134,
            mantissa: 0xc81b228f_fabe28fd_375f45ba_a54fbf21_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Pos,
            exponent: -132,
            mantissa: 0xa7333d38_aab57ac6_8b79e78b_70156497_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Neg,
            exponent: -136,
            mantissa: 0x8192f720_f4885d34_38f2fac1_eef9be7c_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Neg,
            exponent: -137,
            mantissa: 0xff11c8a7_dc893c43_292f724a_c72c3ea6_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Pos,
            exponent: -141,
            mantissa: 0xbfc26be2_87261cc7_af5ba756_96db52ec_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Pos,
            exponent: -142,
            mantissa: 0xb57eeabd_f0f09c2a_9fade6cd_01dfd1a1_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Neg,
            exponent: -146,
            mantissa: 0x82029d5f_a1c35037_7bf6d881_9526863f_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Neg,
            exponent: -148,
            mantissa: 0x95754f62_46c7e687_83b49641_54594014_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Pos,
            exponent: -153,
            mantissa: 0xc9f5cf8e_eda745d0_b939100c_73b13ead_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Pos,
            exponent: -155,
            mantissa: 0xa1451d1f_18509b3f_0643e79c_66e616b1_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Neg,
            exponent: -160,
            mantissa: 0xccec6217_bd010033_93826bad_7f47cdf3_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Neg,
            exponent: -163,
            mantissa: 0xf64c0282_2d52c324_07946f99_8fd9d789_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Pos,
            exponent: -167,
            mantissa: 0x9309b2ce_74ca7673_154d5522_c7dbd59e_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Pos,
            exponent: -170,
            mantissa: 0x8c79815e_dd4bfcff_584c9cb0_9282911d_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Neg,
            exponent: -175,
            mantissa: 0x9de4e65b_8212807a_9c1a000e_a76e9094_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Neg,
            exponent: -179,
            mantissa: 0xf86e2033_df47f0f1_9728b7de_47e3ad6b_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Pos,
            exponent: -183,
            mantissa: 0x837248fa_2c0bb4af_ccb29a60_49eb5c4b_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Pos,
            exponent: -187,
            mantissa: 0xb009dd63_b97ea768_48f8ff41_62e3177d_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Neg,
            exponent: -192,
            mantissa: 0xb0b5a6b0_64ceeb0c_1c760717_777c30a2_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Neg,
            exponent: -196,
            mantissa: 0xc9fdb7e0_8b53257f_134b98fc_253bc101_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Pos,
            exponent: -201,
            mantissa: 0xbcb1d1dd_4efb420e_1bba6db7_826e8888_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Pos,
            exponent: -205,
            mantissa: 0xc9bf942c_b1a95261_e36fb9b8_5d312eed_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Neg,
            exponent: -210,
            mantissa: 0xc1bfc972_1ae385e1_fcc4bc56_d7196500_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Neg,
            exponent: -214,
            mantissa: 0x8525f47d_ef52100f_e83e85cc_d391b0f4_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Pos,
            exponent: -221,
            mantissa: 0xce067125_d6286579_ac0a77a8_1b75de8d_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Pos,
            exponent: -223,
            mantissa: 0xfa37f349_8c515cbd_3df4a0f0_65548099_u128,
        },
    ],
    [
        DyadicFloat128 {
            sign: DyadicSign::Neg,
            exponent: -130,
            mantissa: 0xee0a750a_744fc6b2_eb38843e_96f24072_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Neg,
            exponent: -231,
            mantissa: 0xd3a86d82_4d1de399_fad0aef6_bd06d983_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Pos,
            exponent: -131,
            mantissa: 0xee0a750a_744fc6b2_eb38843e_97825e94_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Neg,
            exponent: -136,
            mantissa: 0xd81bff4e_79528de0_1f80d14f_9c929c8e_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Neg,
            exponent: -134,
            mantissa: 0x9b3ebeb0_4725d3e3_351eb626_f46ad108_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Pos,
            exponent: -139,
            mantissa: 0xa55f7ab1_b727bbb1_2462fe52_1ea44624_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Pos,
            exponent: -139,
            mantissa: 0x9ec53309_1492cab7_4b6b6c12_d0b6b697_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Neg,
            exponent: -144,
            mantissa: 0xae7f4950_d1622900_a3f0097a_3595b0d8_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Neg,
            exponent: -145,
            mantissa: 0xab8c0a50_d5180908_7e84c579_90c08b93_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Pos,
            exponent: -150,
            mantissa: 0xba0a719e_4dd66a56_ce81b61e_0e8a2e9d_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Pos,
            exponent: -152,
            mantissa: 0xe53834f3_9e8ec9a4_836a20e8_27ff6b99_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Neg,
            exponent: -157,
            mantissa: 0xf059d382_ff923433_4058c9fc_eeb6b81c_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Neg,
            exponent: -159,
            mantissa: 0xd0aeeb16_7cf87d51_16120bfc_e6904e2a_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Pos,
            exponent: -164,
            mantissa: 0xd19ca475_6d095e66_17948226_0eda3446_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Pos,
            exponent: -166,
            mantissa: 0x8a0cd8fe_a1477a47_31651fe6_97427176_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Neg,
            exponent: -171,
            mantissa: 0x84417b05_c53ba665_1e1593ff_ba6371b5_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Neg,
            exponent: -174,
            mantissa: 0x8af2e057_ba4e9d9c_356882cb_9831bdf7_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Pos,
            exponent: -180,
            mantissa: 0xfd8f08af_044d92a9_a041c0d8_3ed4e2cd_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Pos,
            exponent: -183,
            mantissa: 0xdc1c131b_84307a14_9f3dd52b_446b4663_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Neg,
            exponent: -188,
            mantissa: 0xbf2ea9cc_3a6362af_82be5fe1_d70e2be7_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Neg,
            exponent: -191,
            mantissa: 0x8ced71c8_ffabe8ee_0a30c787_21d525e4_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Pos,
            exponent: -197,
            mantissa: 0xe95b93c1_c86a2e01_d1db8d6e_6dca0f2f_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Pos,
            exponent: -200,
            mantissa: 0x94c5db96_77b07938_25b1f8b7_7ea85d10_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Neg,
            exponent: -206,
            mantissa: 0xea85e111_6a87a50d_de541446_cd0e4ebf_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Neg,
            exponent: -209,
            mantissa: 0x84561069_9f1b7d2a_69343d93_040528d4_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Pos,
            exponent: -215,
            mantissa: 0xc8aa5055_958e2962_a8b170aa_42b2989c_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Pos,
            exponent: -219,
            mantissa: 0xc425952e_e80dbdf2_5dd5fca6_40f1b5a2_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Neg,
            exponent: -224,
            mantissa: 0x8e1f8cf6_98ef859a_154e15dc_68298480_u128,
        },
    ],
    [
        DyadicFloat128 {
            sign: DyadicSign::Pos,
            exponent: -232,
            mantissa: 0xf2253ffa_3c29d7b8_2eff9292_1c952294_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Pos,
            exponent: -130,
            mantissa: 0xdf99513a_ca003889_cc656551_965add8e_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Neg,
            exponent: -134,
            mantissa: 0x85e16c27_32909f86_8c6d4236_f193274b_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Neg,
            exponent: -132,
            mantissa: 0x936559c5_459b3ca6_55d6e290_2bf2fd6e_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Pos,
            exponent: -137,
            mantissa: 0xaf81f238_b2eb9357_342207c4_135f5747_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Pos,
            exponent: -137,
            mantissa: 0xe5553b71_a4c967ae_99faf970_1e06132d_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Neg,
            exponent: -141,
            mantissa: 0x862f8c62_3681836f_56abd10c_cc5a3357_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Neg,
            exponent: -142,
            mantissa: 0xa7857bea_366a419f_58bad5e7_b5f05dac_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Pos,
            exponent: -147,
            mantissa: 0xbe8f1404_a710cf7d_f87fa811_39914178_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Pos,
            exponent: -148,
            mantissa: 0x8d683989_1423c60e_d2125602_f798fa57_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Neg,
            exponent: -153,
            mantissa: 0x9b010e46_59e1a19f_e69a18b2_b745b7be_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Neg,
            exponent: -155,
            mantissa: 0x9b8d7363_ac8c04da_8182a88f_a7ae49bf_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Pos,
            exponent: -160,
            mantissa: 0xa36703bd_45f6b898_65c946f5_d6775570_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Pos,
            exponent: -163,
            mantissa: 0xf11ae9ab_25b1c402_97934dc7_e90f63b2_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Neg,
            exponent: -168,
            mantissa: 0xf204712c_e397d85c_4776be5e_2fce0e31_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Neg,
            exponent: -170,
            mantissa: 0x8af83d64_3f4190ce_992b1011_a6e6d585_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Pos,
            exponent: -175,
            mantissa: 0x85227bbf_b547736c_da622ab6_d85aceaa_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Pos,
            exponent: -179,
            mantissa: 0xf7ec55a9_7c7607a0_877b7e76_22c7c871_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Neg,
            exponent: -184,
            mantissa: 0xe2a75dad_69ab3794_51155cba_fd3131ff_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Neg,
            exponent: -187,
            mantissa: 0xb05a1e5d_b02d0a39_1071218d_02d018a5_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Pos,
            exponent: -192,
            mantissa: 0x99eba556_2e234aa3_ee26d9ef_c916d508_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Pos,
            exponent: -196,
            mantissa: 0xccdfdda1_4526846e_a2d76ded_d7d29d54_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Neg,
            exponent: -201,
            mantissa: 0xaae04adc_e9ef12a1_efad248e_90045a50_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Neg,
            exponent: -205,
            mantissa: 0xc614eb5a_cb209074_656b40ae_682c9721_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Pos,
            exponent: -210,
            mantissa: 0x9e02ea2d_06761188_e53fd261_c0a4fdb6_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Pos,
            exponent: -214,
            mantissa: 0xa1ff5797_f6a885a7_c7d8e51c_6c41ad7c_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Neg,
            exponent: -220,
            mantissa: 0xf62e24f0_eaf4424c_76ac9e69_7cc93d87_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Neg,
            exponent: -224,
            mantissa: 0xe4150d8f_ba0446ce_08327ff7_c7fa1292_u128,
        },
    ],
    [
        DyadicFloat128 {
            sign: DyadicSign::Pos,
            exponent: -130,
            mantissa: 0xd38115f0_426cc97b_a12cd0fc_91ee0039_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Neg,
            exponent: -232,
            mantissa: 0x8a79b78a_a85f606e_048b5818_1cf093f4_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Neg,
            exponent: -131,
            mantissa: 0xd38115f0_426cc97b_a12cd0fc_91c8d1d6_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Pos,
            exponent: -136,
            mantissa: 0x97703ce8_10588c99_9b1cd06c_62689e70_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Pos,
            exponent: -134,
            mantissa: 0x8b18c8e1_855313a1_f1e6ce26_23d20e03_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Neg,
            exponent: -140,
            mantissa: 0xebc03753_97e3b38f_e5b76706_525c9272_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Neg,
            exponent: -139,
            mantissa: 0x907cab81_58ea6e06_13a7ebcb_2b7a0aa4_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Pos,
            exponent: -144,
            mantissa: 0x80308441_14a55333_c341e781_57c19968_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Pos,
            exponent: -145,
            mantissa: 0x9f1cc65e_23975115_05ce3b25_d403e6a2_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Neg,
            exponent: -150,
            mantissa: 0x8ddf0c0e_328b1020_a18ef1d2_e7caafde_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Neg,
            exponent: -152,
            mantissa: 0xd87c4dbe_30f9d387_f7065af6_a99ab042_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Pos,
            exponent: -157,
            mantissa: 0xbe1fc2c4_1024db23_6f96a2e0_67f72cdb_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Pos,
            exponent: -159,
            mantissa: 0xc811b30a_6426dd12_b013e1de_aa9213d6_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Neg,
            exponent: -164,
            mantissa: 0xab249068_9c3ee275_766933af_17f5fbd0_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Neg,
            exponent: -166,
            mantissa: 0x85f0c089_2ae7ec9c_1fe5c3b1_6ce2ae51_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Pos,
            exponent: -172,
            mantissa: 0xddcb892c_35fd6173_96b68be0_fd5733ee_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Pos,
            exponent: -174,
            mantissa: 0x880e18f0_ef9112a0_62d90544_5a270cfc_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Neg,
            exponent: -180,
            mantissa: 0xd94c66d3_b73c91cd_112332c0_06ef6b07_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Neg,
            exponent: -183,
            mantissa: 0xd9156622_0eb09cfc_7dabf85f_35d90e52_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Pos,
            exponent: -188,
            mantissa: 0xa6e51bd1_98c30c19_228d8513_8d7eaced_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Pos,
            exponent: -191,
            mantissa: 0x8bb74a76_ecd1b899_c624a1fe_7663f05a_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Neg,
            exponent: -197,
            mantissa: 0xcea786be_7918701d_97a80ea3_375c630b_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Neg,
            exponent: -200,
            mantissa: 0x943ac27f_0e6d07cf_485bd995_32fc79e6_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Pos,
            exponent: -206,
            mantissa: 0xd2e01211_168606a8_099fc6d5_e17bc887_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Pos,
            exponent: -209,
            mantissa: 0x83e56a89_713fa4cf_adc33ba6_1a2d6d63_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Neg,
            exponent: -215,
            mantissa: 0xb48c2c46_861caa5e_4f79cffb_a32825a8_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Neg,
            exponent: -219,
            mantissa: 0xc69cf73b_a4d85aab_62591315_0e81d6ae_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Pos,
            exponent: -224,
            mantissa: 0x8690e689_5fbb05e7_d345aed3_ea7dfded_u128,
        },
    ],
    [
        DyadicFloat128 {
            sign: DyadicSign::Neg,
            exponent: -232,
            mantissa: 0xe66cecd5_4d8ab415_9e42ebba_ce7ffc33_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Neg,
            exponent: -130,
            mantissa: 0xc92e1acc_47714ab7_3873c84e_337da5d1_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Pos,
            exponent: -135,
            mantissa: 0xc312a49b_3ac8efc4_ac7f4ee8_a9969dd0_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Pos,
            exponent: -132,
            mantissa: 0x85228901_ceb51342_85dad4aa_438d750c_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Neg,
            exponent: -137,
            mantissa: 0x809d9c67_dc948e6b_31acc39a_47e8cb58_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Neg,
            exponent: -137,
            mantissa: 0xd125090a_fb3421a2_a14cf0f6_1a43e635_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Pos,
            exponent: -142,
            mantissa: 0xc7c8a109_bbfd6430_96776e97_b2482e70_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Pos,
            exponent: -142,
            mantissa: 0x9ae8be76_00b97bee_5fa1ea41_3ccd5d33_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Neg,
            exponent: -147,
            mantissa: 0x9141d498_8a04af19_6de8d848_1ddd0a39_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Neg,
            exponent: -148,
            mantissa: 0x84cf389c_97aa6f53_112e033f_d1f8e3c0_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Pos,
            exponent: -154,
            mantissa: 0xf2ef00f1_732b2e92_31e3cbca_c1e4d07c_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Pos,
            exponent: -155,
            mantissa: 0x9442c2a7_54582f56_889ce5dd_6313d3f0_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Neg,
            exponent: -160,
            mantissa: 0x8398614c_2fe90f43_513b942d_eba7558d_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Neg,
            exponent: -163,
            mantissa: 0xe8ba1e5d_9708ce1f_963d9920_d727a700_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Pos,
            exponent: -168,
            mantissa: 0xc7c3943a_2f0f3224_3253e8a1_71bca443_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Pos,
            exponent: -170,
            mantissa: 0x878b3659_4d3e89e1_722c3846_07aa0d6e_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Neg,
            exponent: -176,
            mantissa: 0xe08733e0_0271cb93_6fea354f_ce2963ab_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Neg,
            exponent: -179,
            mantissa: 0xf3d760dd_19f4477c_b7628ece_e35ff65e_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Pos,
            exponent: -184,
            mantissa: 0xc2a1f096_4de60063_1afc3b2f_d083385c_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Pos,
            exponent: -187,
            mantissa: 0xae9c9cfb_694f8b10_1e81a207_297db93e_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Neg,
            exponent: -192,
            mantissa: 0x863dd08c_96d0e2bc_fe17b15e_98f3c83e_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Neg,
            exponent: -196,
            mantissa: 0xcbec8bb9_fd706b57_768052ab_b26f6e3c_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Pos,
            exponent: -201,
            mantissa: 0x96fdca8f_f255afb1_02dc74e8_880e1d26_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Pos,
            exponent: -205,
            mantissa: 0xc6027b45_45d5763d_83ec88a4_d056ca40_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Neg,
            exponent: -210,
            mantissa: 0x8d3c1d94_82f8691c_b733bb21_7c2ab56d_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Neg,
            exponent: -214,
            mantissa: 0xa26fafd6_415803c0_3cbc8699_64c752a4_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Pos,
            exponent: -220,
            mantissa: 0xddfaee2c_761342d4_de0a0c00_43f8e3b3_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Pos,
            exponent: -224,
            mantissa: 0xe4cf9a5e_198716d3_d13fba69_96805e19_u128,
        },
    ],
    [
        DyadicFloat128 {
            sign: DyadicSign::Neg,
            exponent: -130,
            mantissa: 0xc03c0e19_21173a7b_6d18ace5_075f8f56_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Pos,
            exponent: -231,
            mantissa: 0xef58fddd_7c2b2e26_bbedcf57_887f935b_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Pos,
            exponent: -131,
            mantissa: 0xc03c0e19_21173a7b_6d18ace5_06f57059_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Neg,
            exponent: -137,
            mantissa: 0xe34919ba_6ab0fb56_a1298073_28941ffc_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Neg,
            exponent: -135,
            mantissa: 0xfdf36fc2_0423f932_419a9d40_5c9c93cb_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Pos,
            exponent: -140,
            mantissa: 0xb27a388d_42f8d21e_07785ceb_eba07b52_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Pos,
            exponent: -139,
            mantissa: 0x84ffa158_6bd340ec_cb97360c_be47ab19_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Neg,
            exponent: -145,
            mantissa: 0xc530342d_2cbec91a_ba5a7ce6_e94d5ac8_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Neg,
            exponent: -145,
            mantissa: 0x94169353_a61c7003_b62a7623_dae9c058_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Pos,
            exponent: -151,
            mantissa: 0xdedabd44_4d26045b_c56ad848_315dee82_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Pos,
            exponent: -152,
            mantissa: 0xcbece095_306774cc_bf7e9a66_14ae3c0f_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Neg,
            exponent: -157,
            mantissa: 0x98d54348_b9c7d7fd_43e16f75_993701b2_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Neg,
            exponent: -159,
            mantissa: 0xbea43777_1cd42586_cc78d494_6cbe430d_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Pos,
            exponent: -164,
            mantissa: 0x8cbf92c1_0a59e891_ed6a3574_a003b243_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Pos,
            exponent: -166,
            mantissa: 0x80eebfc9_0caff69a_b53fe3e9_3ed591f2_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Neg,
            exponent: -172,
            mantissa: 0xba405cf7_d543037b_c6bd6e1d_7a809d60_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Neg,
            exponent: -174,
            mantissa: 0x8419f50f_c1a30508_863cf454_7ee91836_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Pos,
            exponent: -180,
            mantissa: 0xb9dfbcb2_76e1c932_6656e27a_ad95cffd_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Pos,
            exponent: -183,
            mantissa: 0xd448eb82_7f694bd1_d0bba774_60ac4970_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Neg,
            exponent: -188,
            mantissa: 0x91116e78_4b929531_91b38310_3505f8e8_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Neg,
            exponent: -191,
            mantissa: 0x896dce90_cdb11ffd_47768198_79f59801_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Pos,
            exponent: -197,
            mantissa: 0xb6231c75_b68554b2_0040246c_beb013f6_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Pos,
            exponent: -200,
            mantissa: 0x927f94e6_1a52f678_f0c02499_16349af8_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Neg,
            exponent: -206,
            mantissa: 0xbc17f051_76aa3f58_67888ee7_2a7c8909_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Neg,
            exponent: -209,
            mantissa: 0x82dcc516_9c5d662c_5e0d581e_706e3e1e_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Pos,
            exponent: -215,
            mantissa: 0xa2b8ced7_8019cf22_07e4ea4d_a4d175fe_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Pos,
            exponent: -219,
            mantissa: 0xc5a8b9af_ca568751_f4a92d91_d03047bc_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Neg,
            exponent: -225,
            mantissa: 0xf413f77b_34d4f181_56ec3098_31efc988_u128,
        },
    ],
    [
        DyadicFloat128 {
            sign: DyadicSign::Pos,
            exponent: -235,
            mantissa: 0x95667cea_0277e51c_4cbe3bb5_e107fadc_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Pos,
            exponent: -130,
            mantissa: 0xb8627b36_55a3fa9f_dea7e912_7caa2a37_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Neg,
            exponent: -135,
            mantissa: 0x96339883_8f49b3c0_64ed0631_ce543689_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Neg,
            exponent: -133,
            mantissa: 0xf4925c2d_0b9b0607_f94de436_a8868a8b_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Pos,
            exponent: -138,
            mantissa: 0xc6b61b26_c9040dd7_9682f073_96ae6461_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Pos,
            exponent: -137,
            mantissa: 0xc123d816_b587b2b9_7c0c71e9_9dcbfa4e_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Neg,
            exponent: -142,
            mantissa: 0x9bb092c9_c2b0eb5b_d166cc1f_a7db7b3a_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Neg,
            exponent: -142,
            mantissa: 0x9036d119_515b4c2f_c370458f_8cb0fa0d_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Pos,
            exponent: -148,
            mantissa: 0xe57c26dc_38a88249_67862e6a_7f262394_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Pos,
            exponent: -149,
            mantissa: 0xf9b388bb_98364fb7_47114e6e_ec6cc4e3_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Neg,
            exponent: -154,
            mantissa: 0xc32a5800_803667e0_74e81a0e_8a7bd61c_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Neg,
            exponent: -155,
            mantissa: 0x8ccf6445_8c3b2db1_a8dc38c3_f62ebfda_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Pos,
            exponent: -161,
            mantissa: 0xd758d424_982fd9b9_a2ed7b26_26d4a869_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Pos,
            exponent: -163,
            mantissa: 0xdf33e0cb_e24602f9_5e20842d_9e30fecc_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Neg,
            exponent: -168,
            mantissa: 0xa66fd05e_b4c4c02d_6aab2b71_a09b3652_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Neg,
            exponent: -170,
            mantissa: 0x83249b7c_bc9179ba_b7dec103_0e4d8afa_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Pos,
            exponent: -176,
            mantissa: 0xbe41da8c_31116e74_5218acb5_689a8cb9_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Pos,
            exponent: -179,
            mantissa: 0xedb9d6bf_d45e6470_6dadc7df_5dd7ad7d_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Neg,
            exponent: -184,
            mantissa: 0xa7765f85_c8d557d8_0397a22e_ee35b729_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Neg,
            exponent: -187,
            mantissa: 0xab55b4be_059b7bb8_b7749544_cadbbb7e_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Pos,
            exponent: -193,
            mantissa: 0xea2837f1_cd7d9f7e_5aaa3bc1_61dfc189_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Pos,
            exponent: -196,
            mantissa: 0xc92ec0d5_1563c9d4_9e7992e0_ac4953f1_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Neg,
            exponent: -201,
            mantissa: 0x8544ec53_331996e3_d243e341_fd6c5016_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Neg,
            exponent: -205,
            mantissa: 0xc43aea35_74489a63_5c8f54f2_230dad00_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Pos,
            exponent: -211,
            mantissa: 0xfbf16033_f0548403_b5452d18_53143a7e_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Pos,
            exponent: -214,
            mantissa: 0xa1940948_253c30c3_103a41d8_b9167610_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Neg,
            exponent: -220,
            mantissa: 0xc7ce7786_be110732_cb0e59c4_765bc52d_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Neg,
            exponent: -224,
            mantissa: 0xe3d73222_0142f70f_46376b1c_93980d5e_u128,
        },
    ],
    [
        DyadicFloat128 {
            sign: DyadicSign::Pos,
            exponent: -130,
            mantissa: 0xb16ca6cb_f42cdc5a_1fdf4186_30888fac_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Pos,
            exponent: -231,
            mantissa: 0xbe5d41b5_d99f6cd2_a2d14d4f_0f5ef45e_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Neg,
            exponent: -131,
            mantissa: 0xb16ca6cb_f42cdc5a_1fdf4186_30d06fe1_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Pos,
            exponent: -137,
            mantissa: 0xb2a40daa_f577d527_d8842f90_55e6d617_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Pos,
            exponent: -135,
            mantissa: 0xeafc2bd1_53083ae0_260eb40c_dfcfce96_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Neg,
            exponent: -140,
            mantissa: 0x8d00ad16_65bc7541_f34ab398_8fb02879_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Neg,
            exponent: -140,
            mantissa: 0xf75a757e_c30a5a25_e6922e9e_0192adb1_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Pos,
            exponent: -145,
            mantissa: 0x9d3dbcf2_6c07df65_6d14f5f2_4e283bf0_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Pos,
            exponent: -145,
            mantissa: 0x8aa76a2c_c6977e9d_27660f71_6d02014d_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Neg,
            exponent: -151,
            mantissa: 0xb3fbb175_149d0944_4753ed61_eed92886_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Neg,
            exponent: -152,
            mantissa: 0xc0760f02_177311e1_1bb6d8cf_181fd8f2_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Pos,
            exponent: -158,
            mantissa: 0xfa94b295_cef07ae4_b4d3d11e_8349fb51_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Pos,
            exponent: -159,
            mantissa: 0xb56bde91_23faa509_f11bc1cf_c0193df7_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Neg,
            exponent: -165,
            mantissa: 0xea72c17d_a1aa97ec_85472d62_912baa01_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Neg,
            exponent: -167,
            mantissa: 0xf75b77e5_db113c11_fb048a35_08d1a57e_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Pos,
            exponent: -172,
            mantissa: 0x9d8f55da_e0e16313_ee79ea8e_081c1157_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Pos,
            exponent: -175,
            mantissa: 0xff497844_15a1fac1_4d9bad84_51a4a057_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Neg,
            exponent: -180,
            mantissa: 0x9f90338a_b45ed341_510c1345_df361342_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Neg,
            exponent: -183,
            mantissa: 0xce718560_b361dafc_e06deb22_1391d38f_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Pos,
            exponent: -189,
            mantissa: 0xfc6db74a_7ae2bd1d_e14bf7f4_691b6ced_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Pos,
            exponent: -191,
            mantissa: 0x86648602_df90b904_2b39eb41_4ddd5f4c_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Neg,
            exponent: -197,
            mantissa: 0xa0635dd2_9a1b3b7e_72f76619_6d8aceaa_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Neg,
            exponent: -200,
            mantissa: 0x8ff1c101_01913170_d1d6ed27_86034f80_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Pos,
            exponent: -206,
            mantissa: 0xa76c88c2_59a4df91_815c0698_283ca74c_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Pos,
            exponent: -209,
            mantissa: 0x8119b8e7_bd52c3c1_62556b39_245a3b73_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Neg,
            exponent: -215,
            mantissa: 0x923a92c6_e5cca605_275fd917_30074692_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Neg,
            exponent: -219,
            mantissa: 0xc3a7a40c_4bea255f_3dcbb93f_cbc342b1_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Pos,
            exponent: -225,
            mantissa: 0xdcd1ae58_13452557_2a8d482f_31636ba1_u128,
        },
    ],
    [
        DyadicFloat128 {
            sign: DyadicSign::Neg,
            exponent: -231,
            mantissa: 0xd1c91ee4_e0397ac1_619373ff_b3ef0849_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Neg,
            exponent: -130,
            mantissa: 0xab3269be_1bebd601_defc08eb_73ffb70a_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Pos,
            exponent: -136,
            mantissa: 0xf0774764_232cccbe_0882e861_9c3c4ff5_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Pos,
            exponent: -133,
            mantissa: 0xe3620ae4_b8da54f8_b105f168_5df0e726_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Neg,
            exponent: -138,
            mantissa: 0x9f624e5c_a0faa841_bf948cfc_609c8947_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Neg,
            exponent: -137,
            mantissa: 0xb429d968_23c4243d_7d7a3d1a_ff515a9b_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Pos,
            exponent: -143,
            mantissa: 0xfb19549c_b7ea6ce6_2b73d4da_b9d1a1af_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Pos,
            exponent: -142,
            mantissa: 0x87357d00_d7ab05f6_cef554c9_af74ef88_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Neg,
            exponent: -148,
            mantissa: 0xbaa0923d_66fac55a_af56cce6_3c1ce256_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Neg,
            exponent: -149,
            mantissa: 0xeb9af45f_5c7e5cb1_b7a5c4e1_a0194b89_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Pos,
            exponent: -154,
            mantissa: 0xa07353e5_7e04acd8_5a295b7d_afbb520e_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Pos,
            exponent: -155,
            mantissa: 0x85ce7d71_5599b9ae_709b26c7_4f2c6305_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Neg,
            exponent: -161,
            mantissa: 0xb34043f3_393cd5a5_236c00bc_20ac65f9_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Neg,
            exponent: -163,
            mantissa: 0xd5a4c737_b4e3942a_90117204_d24954ae_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Pos,
            exponent: -168,
            mantissa: 0x8c5ab3a4_4a57db8f_d7504ef1_5bdf514f_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Pos,
            exponent: -171,
            mantissa: 0xfcced011_e3918ca8_73946212_90d4db94_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Neg,
            exponent: -176,
            mantissa: 0xa285c0b2_e1325773_39c8d53f_631ea91e_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Neg,
            exponent: -179,
            mantissa: 0xe69c3b44_628490bc_d50a50a7_a660ffc2_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Pos,
            exponent: -184,
            mantissa: 0x90d081de_84f537bf_d03359f6_70356a85_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Pos,
            exponent: -187,
            mantissa: 0xa72a07d7_7528812b_31b04643_f74ffcfd_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Neg,
            exponent: -193,
            mantissa: 0xcccd9e7c_5cbc9c08_06b2c7e7_ffc6f230_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Neg,
            exponent: -196,
            mantissa: 0xc5478e0c_d5ae42af_192cdf38_86d20b41_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Pos,
            exponent: -202,
            mantissa: 0xeb8e18aa_aeffcee5_0c6bcbc7_049966be_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Pos,
            exponent: -205,
            mantissa: 0xc144f030_c884ff42_48a89372_57f1215c_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Neg,
            exponent: -211,
            mantissa: 0xe0c0c2ee_07ac47c0_db1b299e_92a316d9_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Neg,
            exponent: -214,
            mantissa: 0x9fbdf876_69988c45_39d669cb_3ed77bc3_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Pos,
            exponent: -220,
            mantissa: 0xb3bda05b_3c5006dd_52248419_721cc043_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Pos,
            exponent: -224,
            mantissa: 0xe1aa01f9_ff66555a_5829680e_6e594ab3_u128,
        },
    ],
    [
        DyadicFloat128 {
            sign: DyadicSign::Neg,
            exponent: -130,
            mantissa: 0xa5951c78_d5cda1e9_088e9ff2_518894fc_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Neg,
            exponent: -233,
            mantissa: 0xd9b61dcf_23f88156_cdffbbb5_7992d787_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Pos,
            exponent: -131,
            mantissa: 0xa5951c78_d5cda1e9_088e9ff2_519a7a14_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Neg,
            exponent: -137,
            mantissa: 0x912d2e39_f9119b99_01064ed4_0e41f2c6_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Neg,
            exponent: -135,
            mantissa: 0xdba86c4d_4d9af55b_4030d03b_db6f56c6_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Pos,
            exponent: -141,
            mantissa: 0xe5eda78e_aedde5ee_494d6cc6_710ca7c7_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Pos,
            exponent: -140,
            mantissa: 0xe7f499fe_3685c26b_84727190_206e84fe_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Neg,
            exponent: -145,
            mantissa: 0x80f7886c_1cdddf69_5d62e0af_15ebdf8a_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Neg,
            exponent: -145,
            mantissa: 0x829bad2c_453955c6_00946bfa_5b011350_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Pos,
            exponent: -151,
            mantissa: 0x94d7f3f5_3dd216c0_8e4886a0_4e755e94_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Pos,
            exponent: -152,
            mantissa: 0xb647f00a_f3308133_6f1328c0_96c3ea8e_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Neg,
            exponent: -158,
            mantissa: 0xd1500f3d_781204d4_b26a3f4a_0a58ae4d_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Neg,
            exponent: -159,
            mantissa: 0xacd99828_690490c2_8c55a49b_e6464f64_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Pos,
            exponent: -165,
            mantissa: 0xc6052873_c29696b4_f901f041_b00cdef2_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Pos,
            exponent: -167,
            mantissa: 0xed1a1379_aab482db_c2664fe9_5acca325_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Neg,
            exponent: -172,
            mantissa: 0x869da421_4e7423e6_1e4d93fb_0e3c7a76_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Neg,
            exponent: -175,
            mantissa: 0xf62370a6_6aece49e_64d4f951_625ebd96_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Pos,
            exponent: -180,
            mantissa: 0x89e1ed3a_02efb4eb_3e629eb6_72174c3a_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Pos,
            exponent: -183,
            mantissa: 0xc821da50_780b2155_845c8403_58d09ef2_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Neg,
            exponent: -189,
            mantissa: 0xdc817a59_391cf919_06c502f6_087aaa6f_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Neg,
            exponent: -191,
            mantissa: 0x82edf3a8_167be6c4_e03700a3_04cba67b_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Pos,
            exponent: -197,
            mantissa: 0x8d87aee4_f585594e_300158ee_2c93d006_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Pos,
            exponent: -200,
            mantissa: 0x8cdaa850_1c074278_6035478a_bfe6cf93_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Neg,
            exponent: -206,
            mantissa: 0x951f096e_9d10ed6f_4f47c0c0_210692e4_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Neg,
            exponent: -210,
            mantissa: 0xfda3de36_4c1cb95a_96cd7834_eb9b9ad7_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Pos,
            exponent: -215,
            mantissa: 0x835a8e39_f234c95d_5f51aa77_0438de9f_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Pos,
            exponent: -219,
            mantissa: 0xc0d8599e_f82629ba_3931d96b_5c3c494e_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Neg,
            exponent: -225,
            mantissa: 0xc7a662f9_66bc6c72_723199d5_a4b80373_u128,
        },
    ],
    [
        DyadicFloat128 {
            sign: DyadicSign::Neg,
            exponent: -231,
            mantissa: 0x97e54b6b_3d59e125_bf8fd89f_cf88675a_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Pos,
            exponent: -130,
            mantissa: 0xa07c7fef_84d2f9f3_474ffad3_fda56cae_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Neg,
            exponent: -136,
            mantissa: 0xc61bd14e_22c37570_fde8d0da_062ee212_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Neg,
            exponent: -133,
            mantissa: 0xd5584c98_a5904a90_8e398efb_fccd0bc6_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Pos,
            exponent: -138,
            mantissa: 0x837b9ae6_07b532d5_501f9a08_89d6eab3_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Pos,
            exponent: -137,
            mantissa: 0xa9694d03_401701be_02c39288_5769b24f_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Neg,
            exponent: -143,
            mantissa: 0xcfe02633_ae076660_20efdac2_4f463358_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Neg,
            exponent: -143,
            mantissa: 0xff24412f_69e4704b_e4c1b6ad_a9f5918b_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Pos,
            exponent: -148,
            mantissa: 0x9b5e96ac_7fa2da64_a2cee4dc_fa94f40b_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Pos,
            exponent: -149,
            mantissa: 0xdf43f165_69673650_f749261c_af0fd5c0_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Neg,
            exponent: -154,
            mantissa: 0x868fef6e_3658bb90_0f8d54a0_c1158174_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Neg,
            exponent: -156,
            mantissa: 0xfedf5883_eee67810_64ddbed1_5f1b94ba_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Pos,
            exponent: -161,
            mantissa: 0x97a26809_90b422b1_26948927_a51e5c71_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Pos,
            exponent: -163,
            mantissa: 0xcc909b9f_50792edc_d1e3a93c_a5b68f98_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Neg,
            exponent: -169,
            mantissa: 0xefb5c1d6_2b4c08ce_216f72f5_b45d9e5a_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Neg,
            exponent: -171,
            mantissa: 0xf36002db_15ffcc9f_db463dc2_f10f2fd6_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Pos,
            exponent: -176,
            mantissa: 0x8c24648e_f38b082b_92ad620e_7a918aea_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Pos,
            exponent: -179,
            mantissa: 0xdf2a9868_7e79043f_9ce3eec4_77e64c39_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Neg,
            exponent: -185,
            mantissa: 0xfc2a05c5_82dea3c7_eeb28da6_a2540684_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Neg,
            exponent: -187,
            mantissa: 0xa28f5cdf_82293f79_48db6efc_9cf82638_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Pos,
            exponent: -193,
            mantissa: 0xb3f8f6e9_e35e8d40_2605415a_eb339610_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Pos,
            exponent: -196,
            mantissa: 0xc0b49095_7433ecb6_acb69f17_269e3222_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Neg,
            exponent: -202,
            mantissa: 0xd0d06f6a_534cf850_62ec2eab_8043b68a_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Neg,
            exponent: -205,
            mantissa: 0xbd8ccec2_356fd795_cea0b3ea_a39ef464_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Pos,
            exponent: -211,
            mantissa: 0xc8dd8505_b0885f8e_4fdf32ee_4426b693_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Pos,
            exponent: -214,
            mantissa: 0x9d3afc5b_efa06771_104c583d_b89090be_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Neg,
            exponent: -220,
            mantissa: 0xa1d63467_ddc1e9c2_5ec9cd09_376d35e8_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Neg,
            exponent: -224,
            mantissa: 0xde984a66_6ea3688c_50058425_5978569b_u128,
        },
    ],
    [
        DyadicFloat128 {
            sign: DyadicSign::Pos,
            exponent: -130,
            mantissa: 0x9bd56756_4c3dc8cc_bbc142cd_0066df30_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Neg,
            exponent: -231,
            mantissa: 0x83133c8c_bfced639_37c05d5c_9bb32ede_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Neg,
            exponent: -131,
            mantissa: 0x9bd56756_4c3dc8cc_bbc142cd_0040b4d4_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Pos,
            exponent: -138,
            mantissa: 0xf1fedfeb_280a4efa_5914828d_b9a060be_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Pos,
            exponent: -135,
            mantissa: 0xcef3d19d_e1d4bbb2_3d0140a5_368fe0bd_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Neg,
            exponent: -141,
            mantissa: 0xc00ec88f_de8362ad_de2183d7_9aed2de7_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Neg,
            exponent: -140,
            mantissa: 0xdb02d2d6_f1d3fc62_da18ee1b_e0821b20_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Pos,
            exponent: -146,
            mantissa: 0xd8513040_a251e6c6_4327e236_76d5a032_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Pos,
            exponent: -146,
            mantissa: 0xf76759a0_d6c19499_31511fd7_33b2e2c5_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Neg,
            exponent: -152,
            mantissa: 0xfb11fe3d_f4fcf96c_b1030890_b1a1d63b_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Neg,
            exponent: -152,
            mantissa: 0xad4e54a5_06beea82_5757e6ab_0554e509_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Pos,
            exponent: -158,
            mantissa: 0xb1c5d17a_faf35902_5c74264c_26e6ce15_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Pos,
            exponent: -159,
            mantissa: 0xa50da10b_f7c31f5d_20c4c950_a318e546_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Neg,
            exponent: -165,
            mantissa: 0xa9863897_39c6a256_56c76aa5_ecca11e6_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Neg,
            exponent: -167,
            mantissa: 0xe3729047_448a31b3_a7b1b943_5f06af33_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Pos,
            exponent: -173,
            mantissa: 0xe875a981_d40251ef_5cdd9301_3e2b5072_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Pos,
            exponent: -175,
            mantissa: 0xed35701e_f7509a01_0268da66_2dfbc698_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Neg,
            exponent: -181,
            mantissa: 0xf02f7c55_75643d56_7db6bd70_e5db8652_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Neg,
            exponent: -183,
            mantissa: 0xc1bcb611_d886bd96_fd996b5b_6a4f9d8b_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Pos,
            exponent: -189,
            mantissa: 0xc1b81df0_f1bdef39_c333e7b3_6d368c70_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Pos,
            exponent: -192,
            mantissa: 0xfe907785_857708d9_dcca2271_6c7d0094_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Neg,
            exponent: -198,
            mantissa: 0xfac082e7_d2057e06_c6830bf9_79655eae_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Neg,
            exponent: -200,
            mantissa: 0x89773a8c_90d7e293_514e00eb_081a4e4c_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Pos,
            exponent: -206,
            mantissa: 0x852571f2_8260b745_0d394feb_e317176d_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Pos,
            exponent: -210,
            mantissa: 0xf86b5fcf_1084a292_2467df20_fd6164b5_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Neg,
            exponent: -216,
            mantissa: 0xec4cbf12_87cc5002_e18ac460_5ee0ee57_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Neg,
            exponent: -219,
            mantissa: 0xbd79c291_fa6c57aa_1d74d3d3_fc53ff6d_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Pos,
            exponent: -225,
            mantissa: 0xb4a838c0_53fe1b3e_cd349a6c_9de96e49_u128,
        },
    ],
    [
        DyadicFloat128 {
            sign: DyadicSign::Pos,
            exponent: -231,
            mantissa: 0x8456143e_06a2a686_054014ee_329e6c2c_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Neg,
            exponent: -130,
            mantissa: 0x97903724_c84e3969_47964468_6c2c55c8_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Pos,
            exponent: -136,
            mantissa: 0xa6dfb54f_dc06bbf5_63707a8b_86e442fd_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Pos,
            exponent: -133,
            mantissa: 0xc99b2219_8c7685a4_dd7d0d63_24eca052_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Neg,
            exponent: -139,
            mantissa: 0xddb551ea_74f3c21a_a877fbe0_8d9de826_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Neg,
            exponent: -137,
            mantissa: 0xa05552ec_a5ec6641_c0978fbc_6d92f1f0_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Pos,
            exponent: -143,
            mantissa: 0xafb0db33_094d3011_9d1c00b3_33672286_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Pos,
            exponent: -143,
            mantissa: 0xf20ac1c8_f96158ec_442366f7_5e3f7e44_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Neg,
            exponent: -148,
            mantissa: 0x83d3d03a_2fba4d39_2c614a3f_f79873aa_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Neg,
            exponent: -149,
            mantissa: 0xd4754bdb_3802b83f_36680a4e_ba24d350_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Pos,
            exponent: -155,
            mantissa: 0xe58859d7_ff57f694_c0d7d588_f4c8e295_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Pos,
            exponent: -156,
            mantissa: 0xf36bc43c_fcc50825_745218ed_74d19180_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Neg,
            exponent: -161,
            mantissa: 0x822223c8_202a714a_2f934330_d6fafd4d_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Neg,
            exponent: -163,
            mantissa: 0xc42864cd_bb6a817c_3a2994a0_1a14dec4_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Pos,
            exponent: -169,
            mantissa: 0xcf27d47f_a97c8214_2972e3cb_25a3bc50_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Pos,
            exponent: -171,
            mantissa: 0xea5960d7_ecb27b54_9d72380c_b0a3ee00_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Neg,
            exponent: -177,
            mantissa: 0xf403f158_013cf7cc_a335448b_d4dc9dec_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Neg,
            exponent: -179,
            mantissa: 0xd7ca152c_f0eab346_aec8e90b_b199f80e_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Pos,
            exponent: -185,
            mantissa: 0xdd326208_458e789b_41a3a543_37d34bb8_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Pos,
            exponent: -187,
            mantissa: 0x9dd30a6a_52396298_7d44afd5_eb09e664_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Neg,
            exponent: -193,
            mantissa: 0x9f0e72c2_c4b0a387_e4fdcf90_12730585_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Neg,
            exponent: -196,
            mantissa: 0xbbcf736e_46767893_d4aa7058_a3541b73_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Pos,
            exponent: -202,
            mantissa: 0xb9e3dd9d_624c9311_5e30cc48_b4afc664_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Pos,
            exponent: -205,
            mantissa: 0xb965237b_e3de23a0_51f1eec8_2fd1964f_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Neg,
            exponent: -211,
            mantissa: 0xb40e458a_fd402fe8_ceeae475_aafcba2f_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Neg,
            exponent: -214,
            mantissa: 0x9a497b5b_0aeb5a2f_44f4a375_d2e638ee_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Pos,
            exponent: -220,
            mantissa: 0x92036a61_85354c18_0583c31e_22fa80b5_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Pos,
            exponent: -224,
            mantissa: 0xdae93a93_3dcda701_0c10e509_ac4d39d5_u128,
        },
    ],
    [
        DyadicFloat128 {
            sign: DyadicSign::Neg,
            exponent: -130,
            mantissa: 0x93a040cf_8e554c91_9b23808a_b20d7391_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Neg,
            exponent: -231,
            mantissa: 0xa7cabc1d_952731fc_b36d157f_3a97a9a1_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Pos,
            exponent: -131,
            mantissa: 0x93a040cf_8e554c91_9b23808a_b2394acd_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Neg,
            exponent: -138,
            mantissa: 0xcdb79be8_5007ca96_f41ef0f1_0eeaa724_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Neg,
            exponent: -135,
            mantissa: 0xc4346ba0_0df97111_249d3675_d9ecafe4_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Pos,
            exponent: -141,
            mantissa: 0xa3853e5e_5a1b52da_ffbf843b_6326e5e7_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Pos,
            exponent: -140,
            mantissa: 0xcff4a671_e9f4154b_29365930_75acaa9b_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Neg,
            exponent: -146,
            mantissa: 0xb8b30730_f8786135_2a68b135_b7ac61dd_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Neg,
            exponent: -146,
            mantissa: 0xeb6fe6d8_1301431a_68d0d6dc_5540e0b7_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Pos,
            exponent: -152,
            mantissa: 0xd73ce1ae_80de5100_6603a221_bb0c9fb3_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Pos,
            exponent: -152,
            mantissa: 0xa56285f3_65016c21_cb83eb29_8942481f_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Neg,
            exponent: -158,
            mantissa: 0x992c2049_087a2710_8687eba8_ac1f8a1c_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Neg,
            exponent: -159,
            mantissa: 0x9e04b42e_02f125f7_85f2bfa8_3be694e9_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Pos,
            exponent: -165,
            mantissa: 0x92ec2cc9_338710a0_8c1c4f6b_941ae3ad_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Pos,
            exponent: -167,
            mantissa: 0xda849447_a295c307_dca33e29_750d0ac2_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Neg,
            exponent: -173,
            mantissa: 0xcac36864_3486bda5_14264802_7bec139b_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Neg,
            exponent: -175,
            mantissa: 0xe4b9d0ba_ffdda7ca_9b1fbd20_43085a94_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Pos,
            exponent: -181,
            mantissa: 0xd2eb0b8b_d38242f0_9128d861_f6d62a88_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Pos,
            exponent: -183,
            mantissa: 0xbb7d4f70_d2017e40_2aad82bc_fc761ab7_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Neg,
            exponent: -189,
            mantissa: 0xab483e7f_5caac56a_d99d9ead_304ad858_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Neg,
            exponent: -192,
            mantissa: 0xf73a6f2f_9ff4ea1c_386dc363_52098ad9_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Pos,
            exponent: -198,
            mantissa: 0xdf37a506_e32c4566_1abd939d_322fb030_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Pos,
            exponent: -200,
            mantissa: 0x85f462ad_74b7e8b1_64609c2f_e972e2de_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Neg,
            exponent: -207,
            mantissa: 0xee9e0ecd_9bf897bb_49660d08_aac1a6bd_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Neg,
            exponent: -210,
            mantissa: 0xf2d6bd45_bb57915b_91662c8f_bb161e3f_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Pos,
            exponent: -216,
            mantissa: 0xd514684b_8fec5a45_419e868b_6ff78812_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Pos,
            exponent: -219,
            mantissa: 0xb9c15da0_8243054d_166d37d0_9ed30ece_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Neg,
            exponent: -225,
            mantissa: 0xa3c622a1_132843b2_9cf629f7_d3504832_u128,
        },
    ],
    [
        DyadicFloat128 {
            sign: DyadicSign::Neg,
            exponent: -232,
            mantissa: 0x8db42b68_9f9327f1_6976dc7f_eca46c71_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Pos,
            exponent: -130,
            mantissa: 0x8ffaf5ee_e9e1cc83_bfacb15a_99f2e924_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Neg,
            exponent: -136,
            mantissa: 0x8f101a99_2321e297_74a5c0bd_a4185772_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Neg,
            exponent: -133,
            mantissa: 0xbf9a836a_50918b48_5ab795a4_0a05f064_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Pos,
            exponent: -139,
            mantissa: 0xbe32e498_1629d4e4_40488269_61f7390f_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Pos,
            exponent: -137,
            mantissa: 0x988bc0b1_9a8a8410_c2f905bb_e40c1fed_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Neg,
            exponent: -143,
            mantissa: 0x96fcd7d2_90949d00_fdb8e881_d6155ef7_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Neg,
            exponent: -143,
            mantissa: 0xe6aea54e_bc50a24e_1fd8ab43_3eac0376_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Pos,
            exponent: -149,
            mantissa: 0xe3396a6b_7bbb2cb8_64e95340_dfe524ed_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Pos,
            exponent: -149,
            mantissa: 0xcaf25881_0ca77e8c_9dbc5e22_0a789611_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Neg,
            exponent: -155,
            mantissa: 0xc69131b0_a07092e9_b0ebd5c8_5139ef51_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Neg,
            exponent: -156,
            mantissa: 0xe927e4a2_c5fbf410_4cc669b9_8f69dc00_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Pos,
            exponent: -162,
            mantissa: 0xe2316ced_31cb85ec_8fd92ea3_40e78cc4_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Pos,
            exponent: -163,
            mantissa: 0xbc757d8a_cb26fbd8_22a2d4b4_66b08ad1_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Neg,
            exponent: -169,
            mantissa: 0xb4faea4c_5b3058cb_f5cd5a61_c013399d_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Neg,
            exponent: -171,
            mantissa: 0xe1e272d8_7d5a4699_7ca9f260_53348933_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Pos,
            exponent: -177,
            mantissa: 0xd665788c_1a7c5a41_8ad04fb7_5a1a0843_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Pos,
            exponent: -179,
            mantissa: 0xd0b0e31d_5d8c18d1_a6e2e33a_c203b840_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Neg,
            exponent: -185,
            mantissa: 0xc3811257_b3d8172f_10670560_0d364cf3_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Neg,
            exponent: -187,
            mantissa: 0x99243860_5d0a635a_16980eb6_cb092c8a_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Pos,
            exponent: -193,
            mantissa: 0x8d6e520c_f62b1cd1_3527fe5e_66161d10_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Pos,
            exponent: -196,
            mantissa: 0xb6d41e23_a733a9d5_1df22546_bddf44a7_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Neg,
            exponent: -202,
            mantissa: 0xa648321e_30604d72_525e3284_06d04053_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Neg,
            exponent: -205,
            mantissa: 0xb5092b5d_f6f91d0c_06f4b73f_e8c643b9_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Pos,
            exponent: -211,
            mantissa: 0xa2011f8e_fe2417dd_8ca07470_5d9ea207_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Pos,
            exponent: -214,
            mantissa: 0x9718d97d_f8f86350_ec36abd5_b4898893_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Neg,
            exponent: -220,
            mantissa: 0x841baf73_8d2ed064_62e805d6_cb34c0b6_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Neg,
            exponent: -224,
            mantissa: 0xd6d6ee23_3ced6aae_7a3a5ca8_079e5c44_u128,
        },
    ],
    [
        DyadicFloat128 {
            sign: DyadicSign::Pos,
            exponent: -130,
            mantissa: 0x8c979313_d3a7177b_f1ead7a9_6c2609a2_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Neg,
            exponent: -232,
            mantissa: 0xf8cd416a_0ac9a259_76051a17_cbbd9c2f_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Neg,
            exponent: -131,
            mantissa: 0x8c979313_d3a7177b_f1ead7a9_6c088f36_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Pos,
            exponent: -138,
            mantissa: 0xb1ad9c57_ff34bc4e_e8467c31_a2447177_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Pos,
            exponent: -135,
            mantissa: 0xbaf67555_fbfc30df_21a7186b_d2cdad3a_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Neg,
            exponent: -141,
            mantissa: 0x8d64f00b_27be071a_bcc72ad5_5249b08a_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Neg,
            exponent: -140,
            mantissa: 0xc6623251_ad857528_ad88f90f_89edaba7_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Pos,
            exponent: -146,
            mantissa: 0xa00a6cdf_9c4e14d0_cac58619_b9debc97_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Pos,
            exponent: -146,
            mantissa: 0xe0f82d16_c28b2571_7e630aa2_92969d6a_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Neg,
            exponent: -152,
            mantissa: 0xbb0e8c75_d8257559_ba4117c3_c7496fd0_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Neg,
            exponent: -152,
            mantissa: 0x9e5cd788_1f07abd6_8b23791c_13ccb751_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Pos,
            exponent: -158,
            mantissa: 0x859e5aad_df2e9212_2c4dd29c_8ad3f897_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Pos,
            exponent: -159,
            mantissa: 0x97ae3b1d_99b33438_350b401c_8ad5732a_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Neg,
            exponent: -165,
            mantissa: 0x80bbb8ed_ad2f45d4_9e711ff0_e367a6e1_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Neg,
            exponent: -167,
            mantissa: 0xd252f3cf_9cfdca95_df01ef69_91a7858f_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Pos,
            exponent: -173,
            mantissa: 0xb2894662_6cce60f9_08cecb5b_c9a251f8_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Pos,
            exponent: -175,
            mantissa: 0xdcc8abad_4ed7edf7_613eb3d0_bf7f47e0_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Neg,
            exponent: -181,
            mantissa: 0xbab239d2_be4a296a_ff2cdbba_53403181_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Neg,
            exponent: -183,
            mantissa: 0xb58387f6_73fe1877_7b314206_b1515a82_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Pos,
            exponent: -189,
            mantissa: 0x9870fcf8_ac4ad9a9_bf0f1b7e_51f2f0a5_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Pos,
            exponent: -192,
            mantissa: 0xf00d94df_564ed9bf_52f2b971_db434d89_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Neg,
            exponent: -198,
            mantissa: 0xc7c340fb_658a7b6b_a2ea0614_ba405971_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Neg,
            exponent: -200,
            mantissa: 0x8270b576_7130625e_6a785757_605f0313_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Pos,
            exponent: -207,
            mantissa: 0xd6b73c08_e552c2c0_c5eef9bc_e717ccd7_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Pos,
            exponent: -210,
            mantissa: 0xed1e4abd_c27ab4c4_50aac7ed_cf262b59_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Neg,
            exponent: -216,
            mantissa: 0xc0c30e71_8a252807_afdf0a42_758dbcde_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Neg,
            exponent: -219,
            mantissa: 0xb5d8f465_c207b539_587624ed_dec4455a_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Pos,
            exponent: -225,
            mantissa: 0x94d970cc_4e3cddd8_66487315_12c3580c_u128,
        },
    ],
    [
        DyadicFloat128 {
            sign: DyadicSign::Pos,
            exponent: -232,
            mantissa: 0x95314f91_7e49fcfb_17f1b102_790f91fc_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Neg,
            exponent: -130,
            mantissa: 0x896eaaea_5f1598f3_d526d953_c45eadba_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Pos,
            exponent: -137,
            mantissa: 0xf8d7718c_0c68c524_b40b06d2_54ceadca_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Pos,
            exponent: -133,
            mantissa: 0xb6f32121_54187dcc_378541e1_c595e201_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Neg,
            exponent: -139,
            mantissa: 0xa57efc4f_e530d2bc_acce6d32_de62b59c_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Neg,
            exponent: -137,
            mantissa: 0x91c67e09_d63b8735_70d33f97_08df3262_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Pos,
            exponent: -143,
            mantissa: 0x838ce89f_007292c6_57a57a5b_44d88214_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Pos,
            exponent: -143,
            mantissa: 0xdcba3c0d_2935a05c_ac7ca419_d265c395_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Neg,
            exponent: -149,
            mantissa: 0xc6641538_6ba978ea_55150a77_63b896e1_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Neg,
            exponent: -149,
            mantissa: 0xc285743c_612ddb7a_3189297d_a96d5524_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Pos,
            exponent: -155,
            mantissa: 0xaddd39f1_4176f336_9e533f19_7bed023a_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Pos,
            exponent: -156,
            mantissa: 0xdff09cb7_41b22250_88edcd3d_01d5f661_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Neg,
            exponent: -162,
            mantissa: 0xc6bedf25_120d5756_af520ace_881636eb_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Neg,
            exponent: -163,
            mantissa: 0xb56fdebb_9ba46ac4_72dfe12c_cff7d7b5_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Pos,
            exponent: -169,
            mantissa: 0x9fa8ee59_cc5693fe_6d9326ac_716a091f_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Pos,
            exponent: -171,
            mantissa: 0xda06985f_42f9e5e7_ebcf114c_ecc55db0_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Neg,
            exponent: -177,
            mantissa: 0xbdfaa243_96ce618a_32d9f779_534fdc7f_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Neg,
            exponent: -179,
            mantissa: 0xc9f8b1d0_49bcf080_516d9e6c_e25deaea_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Pos,
            exponent: -185,
            mantissa: 0xae0e53e9_afd71543_4e598b2d_45ae6800_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Pos,
            exponent: -187,
            mantissa: 0x949da1d3_c365e621_e9cd7add_36cbedc5_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Neg,
            exponent: -194,
            mantissa: 0xfd0ef810_dce0ece8_c5b2d1f8_5b2366f5_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Neg,
            exponent: -196,
            mantissa: 0xb1e8599f_7898a87e_d9627a4e_5df8dfe3_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Pos,
            exponent: -202,
            mantissa: 0x957edd56_18cf7488_1e35a177_c5ef84e0_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Pos,
            exponent: -205,
            mantissa: 0xb0a14008_8e169813_682a897a_7957661e_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Neg,
            exponent: -211,
            mantissa: 0x925ce9c8_8dbf06f7_71e78a33_e7662804_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Neg,
            exponent: -214,
            mantissa: 0x93cb6157_f6a0d207_6c5d5988_14a73cbf_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Pos,
            exponent: -221,
            mantissa: 0xefd7a384_fd394413_ec979973_80944ed2_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Pos,
            exponent: -224,
            mantissa: 0xd28e2547_37a4467b_77367ffd_a2c91397_u128,
        },
    ],
    [
        DyadicFloat128 {
            sign: DyadicSign::Neg,
            exponent: -130,
            mantissa: 0x8679f74c_7bb4d5e7_f23af3e2_f2fde2e6_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Neg,
            exponent: -232,
            mantissa: 0xcf0ed7bb_7527199f_eaabed5f_910b795c_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Pos,
            exponent: -131,
            mantissa: 0x8679f74c_7bb4d5e7_f23af3e2_f314544a_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Neg,
            exponent: -138,
            mantissa: 0x9b7a28fb_70a7d851_a5ea6757_61af759c_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Neg,
            exponent: -135,
            mantissa: 0xb2e82ca4_54a34cd7_bb425eec_cd9776fe_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Pos,
            exponent: -142,
            mantissa: 0xf7aae2d0_68b6691e_d6777d1a_61a7cf65_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Pos,
            exponent: -140,
            mantissa: 0xbdff64e8_ae81c3e0_b1454751_606285f6_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Neg,
            exponent: -146,
            mantissa: 0x8c62a4f9_47293e9c_7a9033c3_b0d4241a_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Neg,
            exponent: -146,
            mantissa: 0xd7baa277_563b2e6e_a74babca_c506a038_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Pos,
            exponent: -152,
            mantissa: 0xa4740f6f_9c58e951_45c85d21_cb02d12b_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Pos,
            exponent: -152,
            mantissa: 0x981a0f3c_6880a1bd_8241e556_e26dd01d_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Neg,
            exponent: -159,
            mantissa: 0xeb9e964a_174b8023_83f8e77f_c05739d8_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Neg,
            exponent: -159,
            mantissa: 0x91f61394_e9d369a0_49a511a8_26ba28ea_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Pos,
            exponent: -166,
            mantissa: 0xe3c72266_f6459a52_e35243a1_5f27eac1_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Pos,
            exponent: -167,
            mantissa: 0xcad37b08_3da73033_8add369b_9ccc9f9e_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Neg,
            exponent: -173,
            mantissa: 0x9e8e3bef_4501eb71_cad0e84c_fb4d49f5_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Neg,
            exponent: -175,
            mantissa: 0xd5674a20_a0e2fceb_495c5207_5e42708d_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Pos,
            exponent: -181,
            mantissa: 0xa67e5af4_24b27c60_16b9e79f_492dfb8f_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Pos,
            exponent: -183,
            mantissa: 0xafde1d87_6499e79b_e14d63cf_2ac67a46_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Neg,
            exponent: -189,
            mantissa: 0x888ae246_0459a10d_1de4a1de_3926b789_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Neg,
            exponent: -192,
            mantissa: 0xe927597d_7d16dee6_cc7427dd_d9a00663_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Pos,
            exponent: -198,
            mantissa: 0xb3bc84f9_02d98fec_9cb16bec_1d8eaad7_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Pos,
            exponent: -201,
            mantissa: 0xfdff678a_f8f66b6a_ff6b7622_ffb0cd0f_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Neg,
            exponent: -207,
            mantissa: 0xc211deb5_0da68922_bf0c08e0_88be5036_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Neg,
            exponent: -210,
            mantissa: 0xe768f78b_199efb40_122faafd_b4e8e4f5_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Pos,
            exponent: -216,
            mantissa: 0xaf0375df_f170a0ee_df72f200_8b657416_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Pos,
            exponent: -219,
            mantissa: 0xb1df4bfc_aa33be0a_b73d8a4c_1902e911_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Neg,
            exponent: -225,
            mantissa: 0x87b0cac2_f0d77fa4_0fa9a191_bd3f0283_u128,
        },
    ],
    [
        DyadicFloat128 {
            sign: DyadicSign::Pos,
            exponent: -231,
            mantissa: 0x9232f6a0_03dee785_6bf8d1e0_419c9ba9_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Pos,
            exponent: -130,
            mantissa: 0x83b412be_d6ab4d79_5805277a_1c33952e_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Neg,
            exponent: -137,
            mantissa: 0xdb0153df_55523d8b_7cf24601_9ad7d2bd_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Neg,
            exponent: -133,
            mantissa: 0xaf5e1177_e8495f4a_a886493e_c8f098e0_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Pos,
            exponent: -139,
            mantissa: 0x91b53020_a17307f4_9bf9b27c_d90e2caa_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Pos,
            exponent: -137,
            mantissa: 0x8bd2417d_546c29dd_cf66c597_8c7efc9f_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Neg,
            exponent: -144,
            mantissa: 0xe7e00448_3c3d6f4a_94db2836_a464b1fb_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Neg,
            exponent: -143,
            mantissa: 0xd3ebd894_5629bdae_124644a9_54ec1454_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Pos,
            exponent: -149,
            mantissa: 0xaf20cf8d_ab82195a_6c98550b_6980534f_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Pos,
            exponent: -149,
            mantissa: 0xbb01e751_1a921fd9_a92739eb_fa1afe61_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Neg,
            exponent: -155,
            mantissa: 0x99d05a1b_64412339_dc606a8d_1d6baa01_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Neg,
            exponent: -156,
            mantissa: 0xd7a320bc_22f9fe18_76f4e1a7_514d2dd3_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Pos,
            exponent: -162,
            mantissa: 0xb04d6466_c0bbd0fb_9babdd28_2dcfbb46_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Pos,
            exponent: -163,
            mantissa: 0xaf08dc3d_69c814c6_5499684a_dddcd502_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Neg,
            exponent: -169,
            mantissa: 0x8e13ab2c_81c8cd15_911e72bf_e86a6a3a_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Neg,
            exponent: -171,
            mantissa: 0xd2c3241e_ec48e7ce_e19a2cdb_b30e5a6d_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Pos,
            exponent: -177,
            mantissa: 0xa9a69163_004b12e2_63e3f224_dcb5bc7b_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Pos,
            exponent: -179,
            mantissa: 0xc3aad699_16283200_cee8719b_64de85b2_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Neg,
            exponent: -185,
            mantissa: 0x9c04f294_62f3828a_aeea2345_ce0ac934_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Neg,
            exponent: -187,
            mantissa: 0x904d053c_6b496d0f_efc0ec1e_96253bdc_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Pos,
            exponent: -194,
            mantissa: 0xe3bc2168_a12f6ecd_bbe6dc27_74c884ff_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Pos,
            exponent: -196,
            mantissa: 0xad229c04_a51c7020_5290001c_171195da_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Neg,
            exponent: -202,
            mantissa: 0x87154302_0bf3b6d1_3fd69742_d484b454_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Neg,
            exponent: -205,
            mantissa: 0xac479f7b_94611db3_3127ea61_24d0001a_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Pos,
            exponent: -211,
            mantissa: 0x84cb807f_1a5f9495_88704ee5_fec0f92d_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Pos,
            exponent: -214,
            mantissa: 0x9078e389_339fa826_14088708_92196ae4_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Neg,
            exponent: -221,
            mantissa: 0xda7d3ae1_30e3c1b8_f46c26e1_b97f3259_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Neg,
            exponent: -224,
            mantissa: 0xce2f8b1d_0ea86f08_02a2e4eb_f8ec5ac5_u128,
        },
    ],
    [
        DyadicFloat128 {
            sign: DyadicSign::Pos,
            exponent: -130,
            mantissa: 0x81185d48_79438a22_cd6f32b2_cf38288c_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Neg,
            exponent: -234,
            mantissa: 0xa71feca9_827ad8c7_069fc49d_28779e09_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Neg,
            exponent: -131,
            mantissa: 0x81185d48_79438a22_cd6f32b2_cf33fc1f_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Pos,
            exponent: -138,
            mantissa: 0x898b2504_3266f30b_a1be11ff_263d25d0_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Pos,
            exponent: -135,
            mantissa: 0xabce0ded_e7e4caeb_733f95c1_f47ff748_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Neg,
            exponent: -142,
            mantissa: 0xdb3f0c89_e3345b50_8fa64c35_f877a173_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Neg,
            exponent: -140,
            mantissa: 0xb6936f27_ec62e9dc_3b0bfb76_e557bf14_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Pos,
            exponent: -147,
            mantissa: 0xf8da9058_31c29a13_9ff954e0_81f67153_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Pos,
            exponent: -146,
            mantissa: 0xcf80f3e0_c84d016e_93d7c50e_03429554_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Neg,
            exponent: -152,
            mantissa: 0x92038d8b_cb013725_e108bbc1_df1048ed_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Neg,
            exponent: -152,
            mantissa: 0x927c749c_e721dff9_92cfb8da_113721ae_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Pos,
            exponent: -159,
            mantissa: 0xd1ab31ac_aa31f4a7_e896f8f4_5f62327c_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Pos,
            exponent: -159,
            mantissa: 0x8cc8906a_4ea837b6_a61b7cfb_3ae8e695_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Neg,
            exponent: -166,
            mantissa: 0xcb3b8f0d_8af21e42_34a4c4e7_f12e22a4_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Neg,
            exponent: -167,
            mantissa: 0xc3f6f09a_76bcfca4_a322ac46_7b073c1d_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Pos,
            exponent: -173,
            mantissa: 0x8de66fdd_7d0fe9c9_2871627b_22c35c7d_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Pos,
            exponent: -175,
            mantissa: 0xce91536e_2d2bd63a_70448bd2_f4092225_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Neg,
            exponent: -181,
            mantissa: 0x958071b3_980e70ee_93fceda1_9cfb560c_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Neg,
            exponent: -183,
            mantissa: 0xaa91b391_596e34a0_904dc906_6078c11b_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Pos,
            exponent: -190,
            mantissa: 0xf616724d_4c272b64_0e654bc3_2644b222_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Pos,
            exponent: -192,
            mantissa: 0xe29697da_024a7675_b0988765_337d33a7_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Neg,
            exponent: -198,
            mantissa: 0xa29171cf_a2e11bdf_627393d5_7f4989d3_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Neg,
            exponent: -201,
            mantissa: 0xf759f69d_f24e4918_df511f3a_ea96e1e0_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Pos,
            exponent: -207,
            mantissa: 0xb0329433_a7127a15_13829799_f3eb6d28_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Pos,
            exponent: -210,
            mantissa: 0xe1d04dcb_7a80ca5d_22ed4ff1_2812f379_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Neg,
            exponent: -216,
            mantissa: 0x9f806aef_63c40846_d8e98680_62dd2cc4_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Neg,
            exponent: -219,
            mantissa: 0xade9ff9b_9a5c5b13_84edf214_7eecb003_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Pos,
            exponent: -226,
            mantissa: 0xf830ad5c_ebe59782_461344b9_0b2b1107_u128,
        },
    ],
    [
        DyadicFloat128 {
            sign: DyadicSign::Neg,
            exponent: -234,
            mantissa: 0x85fc0f3d_b5fbf935_e9581fc0_0c2917be_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Neg,
            exponent: -131,
            mantissa: 0xfd459fcd_721ba7df_7d0f7113_47c7cfc0_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Pos,
            exponent: -137,
            mantissa: 0xc2b039db_fd166783_32d3279a_c7c1e8f1_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Pos,
            exponent: -133,
            mantissa: 0xa8a73297_5974b2cc_98d306d1_62fe1416_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Neg,
            exponent: -139,
            mantissa: 0x81914c38_c750d0c9_da9f6413_b647b6e1_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Neg,
            exponent: -137,
            mantissa: 0x8688a9fe_f218d440_81df8246_e38a28c7_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Pos,
            exponent: -144,
            mantissa: 0xce5a2ae1_4eb539c3_65adad42_aa7411d3_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Pos,
            exponent: -143,
            mantissa: 0xcc10b3d8_e296140a_c9f679de_7c76a2dc_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Neg,
            exponent: -149,
            mantissa: 0x9c0c8f7b_47559573_291f0fa6_90126b20_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Neg,
            exponent: -149,
            mantissa: 0xb44329a2_d9897dc1_631f3d87_45465a4d_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Pos,
            exponent: -155,
            mantissa: 0x894b0008_6d1c3b5e_a193c007_67f2e10f_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Pos,
            exponent: -156,
            mantissa: 0xd0201b3c_fb440c39_2ddc83e3_fb3b1d9a_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Neg,
            exponent: -162,
            mantissa: 0x9db3e605_cb602f74_430a5ff7_a04fdd3a_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Neg,
            exponent: -163,
            mantissa: 0xa9301c5c_ab77eacb_e637edce_1f884510_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Pos,
            exponent: -170,
            mantissa: 0xfed03742_07d19ea1_26a970b6_05d6b6b2_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Pos,
            exponent: -171,
            mantissa: 0xcc0ef320_6870023c_2ebfba5b_5ae5e366_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Neg,
            exponent: -177,
            mantissa: 0x988fe2b6_25c443ac_914b1202_25833f2e_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Neg,
            exponent: -179,
            mantissa: 0xbdc7b144_cfe418ed_66cbf11e_1d627cba_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Pos,
            exponent: -185,
            mantissa: 0x8cbbe593_a06e4a2b_1c3cfaa8_4db39256_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Pos,
            exponent: -187,
            mantissa: 0x8c383759_8eeeeb53_89a7af1f_45ca36e4_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Neg,
            exponent: -194,
            mantissa: 0xce17679a_89de0cdf_152ae9a8_c99785e7_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Neg,
            exponent: -196,
            mantissa: 0xa88f39a0_f8cfaa7a_1752fb8f_673f8740_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Pos,
            exponent: -203,
            mantissa: 0xf550c785_0de4e68f_506ff8eb_051809ce_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Pos,
            exponent: -205,
            mantissa: 0xa80c87fe_ca8847b7_8941253d_7f43c30c_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Neg,
            exponent: -212,
            mantissa: 0xf1fdb91c_4f993714_1b268e92_0e29153b_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Neg,
            exponent: -214,
            mantissa: 0x8d314b20_d680945b_2b4273fa_804e8006_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Pos,
            exponent: -221,
            mantissa: 0xc7c2b13e_a359b5ed_5b78134c_e4e7dfaa_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Pos,
            exponent: -224,
            mantissa: 0xc9d25fd7_91f1be65_f5cbd4c7_c0715dbd_u128,
        },
    ],
    [
        DyadicFloat128 {
            sign: DyadicSign::Neg,
            exponent: -131,
            mantissa: 0xf89fd860_737e65e4_9a75f3b8_28094fe9_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Pos,
            exponent: -232,
            mantissa: 0xfde534d0_aef7b9ba_9ea35ce8_8dc8a06f_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Pos,
            exponent: -132,
            mantissa: 0xf89fd860_737e65e4_9a75f3b8_27da4801_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Neg,
            exponent: -139,
            mantissa: 0xf59f37e6_23e003b5_4e2d264a_9d34c96d_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Neg,
            exponent: -135,
            mantissa: 0xa57ba65e_bbdf7b4f_923df34f_b3ea4c71_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Pos,
            exponent: -142,
            mantissa: 0xc3dd8f5f_5763a71a_d8e49f8f_37b969a7_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Pos,
            exponent: -140,
            mantissa: 0xaff31490_18a8721b_83c220ed_8a52093b_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Neg,
            exponent: -147,
            mantissa: 0xde879532_aaf35b32_16b66912_022f7ac4_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Neg,
            exponent: -146,
            mantissa: 0xc8206f41_8324eca0_c7906fbe_352cdc4c_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Pos,
            exponent: -152,
            mantissa: 0x82bfbd3b_4c9eb8ba_ce259532_8259cf55_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Pos,
            exponent: -152,
            mantissa: 0x8d6b6ea8_82e1eee2_9af8d079_4dc574ef_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Neg,
            exponent: -159,
            mantissa: 0xbc1537ab_9c66d457_e299aff7_e118da29_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Neg,
            exponent: -159,
            mantissa: 0x8813ee03_6a229640_a1ad8834_4eb0f24f_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Pos,
            exponent: -166,
            mantissa: 0xb6b2d89d_8b7b24e3_8f31d184_3ab91ffa_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Pos,
            exponent: -167,
            mantissa: 0xbdacf4f6_b9b3e71a_6c3c77b1_03c34c38_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Neg,
            exponent: -174,
            mantissa: 0xffc00795_24d0043c_972e77bd_4555663a_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Neg,
            exponent: -175,
            mantissa: 0xc83df218_2a575b3f_ecdc4d47_b5e5134e_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Pos,
            exponent: -181,
            mantissa: 0x87170ec7_ad6ba695_993e25e7_17dc60ca_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Pos,
            exponent: -183,
            mantissa: 0xa59d4587_da809575_f5f0078b_ba114835_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Neg,
            exponent: -190,
            mantissa: 0xdf036cb5_b8b2e69f_88d0f518_6630e4c5_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Neg,
            exponent: -192,
            mantissa: 0xdc611399_dd254a1b_f276702f_02ed47f3_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Pos,
            exponent: -198,
            mantissa: 0x93c6e674_92011f7e_c32050f8_5a18fae6_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Pos,
            exponent: -201,
            mantissa: 0xf0fda8b5_abe20ac8_0b91e6a8_cee029be_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Neg,
            exponent: -207,
            mantissa: 0xa0ac8a85_849b742c_6a1438be_53051d5a_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Neg,
            exponent: -210,
            mantissa: 0xdc6420d7_41c80154_adeee798_af27e3db_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Pos,
            exponent: -216,
            mantissa: 0x91eaf88f_9910689c_5f7f944f_95e741e0_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Pos,
            exponent: -219,
            mantissa: 0xaa079893_8e50b94b_b8b7542c_331e94ef_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Neg,
            exponent: -226,
            mantissa: 0xe3ba7aa4_23b70a83_dc705529_bb1d2f36_u128,
        },
    ],
    [
        DyadicFloat128 {
            sign: DyadicSign::Pos,
            exponent: -237,
            mantissa: 0x9c53a772_8a68de73_fa9a0c92_d20ebd03_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Pos,
            exponent: -131,
            mantissa: 0xf4393e2b_951612d4_f85db1f5_08846904_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Neg,
            exponent: -137,
            mantissa: 0xae8f7849_559ca442_3863f698_5dcd3324_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Neg,
            exponent: -133,
            mantissa: 0xa2a73d39_cab1aef9_225df145_785a988d_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Pos,
            exponent: -140,
            mantissa: 0xe8661d3f_d369e999_eb39ffae_52774c86_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Pos,
            exponent: -137,
            mantissa: 0x81cc6965_e816f798_b8692ac1_9423929f_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Neg,
            exponent: -144,
            mantissa: 0xb92d8484_fce32cab_c1681981_ff38caa5_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Neg,
            exponent: -143,
            mantissa: 0xc50117af_f7223b64_1b0f4ec0_0e6d0e43_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Pos,
            exponent: -149,
            mantissa: 0x8c2d8f54_bd2a5701_a79dea66_c19cb00c_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Pos,
            exponent: -149,
            mantissa: 0xae2b88b3_a00d73d4_1543bd70_2c9a6841_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Neg,
            exponent: -156,
            mantissa: 0xf6ffb28a_70ed445b_14c59a13_da5c7785_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Neg,
            exponent: -156,
            mantissa: 0xc94c6d1e_27d386da_6d09630b_84b2afbe_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Pos,
            exponent: -162,
            mantissa: 0x8e1a382f_56fd85bd_4d8b83ac_3ff35104_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Pos,
            exponent: -163,
            mantissa: 0xa3d5c762_94e286db_1259c1b9_bdc6a538_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Neg,
            exponent: -170,
            mantissa: 0xe6125178_03e1f2e3_7120bdfc_3bd089bc_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Neg,
            exponent: -171,
            mantissa: 0xc5de6227_3611292a_661b676a_2291ab83_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Pos,
            exponent: -177,
            mantissa: 0x8a106fcd_da1ac1f6_4643e54c_719b2bf9_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Pos,
            exponent: -179,
            mantissa: 0xb84afc59_9335eda7_0c272cc1_96d745e9_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Neg,
            exponent: -186,
            mantissa: 0xff5c3f27_32846cc8_cc0bc656_64af9bcf_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Neg,
            exponent: -187,
            mantissa: 0x8860657a_5d3f6a97_03b3a3b2_86c565cd_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Pos,
            exponent: -194,
            mantissa: 0xbb7a6510_c1ef9594_9c8aed8c_cf057a66_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Pos,
            exponent: -196,
            mantissa: 0xa4340bbf_354a551e_0d6f5114_e41f8ecd_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Neg,
            exponent: -203,
            mantissa: 0xdfc9ed64_eeaba2a3_d1ba5bb7_4a37730c_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Neg,
            exponent: -205,
            mantissa: 0xa3f961f5_27a0f532_43a62ed2_f516d297_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Pos,
            exponent: -212,
            mantissa: 0xdd65ca34_0d60da45_63582bde_5b1cc21e_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Pos,
            exponent: -214,
            mantissa: 0x89fecb7f_e8b17b6a_2e5da9ce_c1b38208_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Neg,
            exponent: -221,
            mantissa: 0xb74bff84_45bf43f8_21956ce4_4ff62232_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Neg,
            exponent: -224,
            mantissa: 0xc586931b_8927b421_9a389013_4839d17f_u128,
        },
    ],
    [
        DyadicFloat128 {
            sign: DyadicSign::Pos,
            exponent: -131,
            mantissa: 0xf00c6d60_e0bf1723_071767c3_4284558a_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Neg,
            exponent: -232,
            mantissa: 0xe63ff7c7_61f0688a_0d77ef64_dcbef4c7_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Neg,
            exponent: -132,
            mantissa: 0xf00c6d60_e0bf1723_071767c3_425c936b_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Pos,
            exponent: -139,
            mantissa: 0xdd10de8a_e8151ea1_209dccc5_46ebec58_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Pos,
            exponent: -135,
            mantissa: 0x9fcf06d8_3f3f7a71_69e1a748_21da8b9c_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Neg,
            exponent: -142,
            mantissa: 0xb05bbe2f_13b112de_47bb95c6_30a37139_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Neg,
            exponent: -140,
            mantissa: 0xa9fcd8d2_e9145f3c_958fd82a_ed98831f_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Pos,
            exponent: -147,
            mantissa: 0xc8859c40_96129404_47ff4390_a97bb13b_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Pos,
            exponent: -146,
            mantissa: 0xc17736fe_ff6f577b_6529564c_ac07cc32_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Neg,
            exponent: -153,
            mantissa: 0xebe6154d_30176cae_819f133b_de194bac_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Neg,
            exponent: -152,
            mantissa: 0x88d2bec3_229abc95_f0f5ead1_75236bc1_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Pos,
            exponent: -159,
            mantissa: 0xa9ea00a8_ca216967_3cba0bc5_6e79b9c5_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Pos,
            exponent: -159,
            mantissa: 0x83c8af7c_3a04663f_c9a5208a_ed91bff6_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Neg,
            exponent: -166,
            mantissa: 0xa555e5b0_a0e3adfb_bc388296_63859a3c_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Neg,
            exponent: -167,
            mantissa: 0xb7e5c717_56d0aba5_05fdc4f2_586dde88_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Pos,
            exponent: -174,
            mantissa: 0xe7e7ed97_d7580d1c_7c32f1c6_bc92dbb4_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Pos,
            exponent: -175,
            mantissa: 0xc2629b4e_9719e452_8e17c27a_bac9cfce_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Neg,
            exponent: -182,
            mantissa: 0xf588ebf3_2c0e85f7_dd628c7a_12da7054_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Neg,
            exponent: -183,
            mantissa: 0xa0fcd5a3_195ffb22_683c4355_6c1b1d1f_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Pos,
            exponent: -190,
            mantissa: 0xcb28aca0_3f4acb32_d2318d0a_1d03eae1_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Pos,
            exponent: -192,
            mantissa: 0xd6871e67_f1cf34b4_1e7dca29_50853c53_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Neg,
            exponent: -198,
            mantissa: 0x86f7074c_99aa6415_6333fd70_21a53ec9_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Neg,
            exponent: -201,
            mantissa: 0xeaf01a3b_5c0d8238_beb5b8ee_8aa32c5f_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Pos,
            exponent: -207,
            mantissa: 0x93229afa_429fd198_d6415c2c_8f7ba52c_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Pos,
            exponent: -210,
            mantissa: 0xd72d7f52_2d0e8264_79005fb0_b5fa10df_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Neg,
            exponent: -216,
            mantissa: 0x85fd006e_b94d4ac9_cb554d52_896a7451_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Neg,
            exponent: -219,
            mantissa: 0xa64168fc_4c95864f_b39bf9b4_0fc063c0_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Pos,
            exponent: -226,
            mantissa: 0xd1a69d24_f8126f0b_771750e1_161ab239_u128,
        },
    ],
    [
        DyadicFloat128 {
            sign: DyadicSign::Pos,
            exponent: -237,
            mantissa: 0xb1e6007a_ede058f6_fd60c598_2f98ed58_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Neg,
            exponent: -131,
            mantissa: 0xec149cd5_b119c554_5b03cea4_566f70c6_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Pos,
            exponent: -137,
            mantissa: 0x9dad2a42_2b3876a0_dfaa598e_c8f4bfaf_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Pos,
            exponent: -133,
            mantissa: 0x9d3ff8b1_16f7bd9c_23ea176a_8dec13d1_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Neg,
            exponent: -140,
            mantissa: 0xd1f5e23b_53036920_01fd8691_5ca46429_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Neg,
            exponent: -138,
            mantissa: 0xfb0d6ef1_d463085e_971ca14d_29f7e2fa_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Pos,
            exponent: -144,
            mantissa: 0xa7622ed0_23c57f46_f27ec347_991f1dba_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Pos,
            exponent: -143,
            mantissa: 0xbe9d940a_cf91e040_a5654f41_39853aaf_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Neg,
            exponent: -150,
            mantissa: 0xfd9f90de_6bfe1dc4_56f8fdf2_e67c54d0_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Neg,
            exponent: -149,
            mantissa: 0xa8a2d43b_691395f7_0b70e2a4_c2bcce24_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Pos,
            exponent: -156,
            mantissa: 0xdfb266fb_e02abd51_467a45d2_0e732d32_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Pos,
            exponent: -156,
            mantissa: 0xc3110256_213e9071_fa8da863_e4809198_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Neg,
            exponent: -162,
            mantissa: 0x80e06a44_444087d5_5aca1949_4462a5ce_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Neg,
            exponent: -163,
            mantissa: 0x9eeb658e_5e6c7744_c8e4cae0_1c2ddbeb_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Pos,
            exponent: -170,
            mantissa: 0xd1010645_9e106a1f_4e0f2715_3d6ee293_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Pos,
            exponent: -171,
            mantissa: 0xc0254eb8_894970dd_e2c7ce93_76300b49_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Neg,
            exponent: -178,
            mantissa: 0xfb50f7c8_c17153d0_bed9928d_f69f6838_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Neg,
            exponent: -179,
            mantissa: 0xb32e416f_e9888cdc_3fe20df2_5a138191_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Pos,
            exponent: -186,
            mantissa: 0xe8e5e23a_570ccc2b_517830ee_024f6857_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Pos,
            exponent: -187,
            mantissa: 0x84c41b94_7682216e_72f6eaa6_ee398fa7_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Neg,
            exponent: -194,
            mantissa: 0xab5ed0cd_dec839fe_1b5b6914_46f95a6e_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Neg,
            exponent: -196,
            mantissa: 0xa012e3f4_367e41f2_fa5959c5_2ff98858_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Pos,
            exponent: -203,
            mantissa: 0xcd0bd128_2308c5dd_46bdc3e3_f438fc53_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Pos,
            exponent: -205,
            mantissa: 0xa0130908_28a03eae_b4998a5a_1a13c882_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Neg,
            exponent: -212,
            mantissa: 0xcb5ae840_3a25b1ae_1c02372f_341d92a3_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Neg,
            exponent: -214,
            mantissa: 0x86e79259_4612484d_09eae383_f2645245_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Pos,
            exponent: -221,
            mantissa: 0xa8c81959_006e4220_3bb2cf56_200d2558_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Pos,
            exponent: -224,
            mantissa: 0xc15695e0_1d35495b_89a8e986_27b36005_u128,
        },
    ],
    [
        DyadicFloat128 {
            sign: DyadicSign::Neg,
            exponent: -131,
            mantissa: 0xe84d90f1_b605eadd_bb013e9d_8856fa09_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Neg,
            exponent: -232,
            mantissa: 0x9d69aaa6_6dee1751_6f1dd378_ee32e780_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Pos,
            exponent: -132,
            mantissa: 0xe84d90f1_b605eadd_bb013e9d_88706e4f_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Neg,
            exponent: -139,
            mantissa: 0xc8587352_3af225c4_bf2087ce_53f4ce8e_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Neg,
            exponent: -135,
            mantissa: 0x9aadc827_dfbf7281_2ebaa309_5edf396e_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Pos,
            exponent: -142,
            mantissa: 0x9fe22ce8_f12d21e3_b983f117_7c46c449_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Pos,
            exponent: -140,
            mantissa: 0xa496630e_8cefeeb6_b3885e88_77839c89_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Neg,
            exponent: -147,
            mantissa: 0xb5e78885_0156869c_982b2d66_a3a9917a_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Neg,
            exponent: -146,
            mantissa: 0xbb6a27b5_41a91492_2b725241_95119da1_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Pos,
            exponent: -153,
            mantissa: 0xd630f7e5_9e1b6004_0c015d0a_c052fa1b_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Pos,
            exponent: -152,
            mantissa: 0x84a1b7e5_51cb29eb_84629add_c9852efa_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Neg,
            exponent: -159,
            mantissa: 0x9a75b04a_fcb4db3c_529ae4c9_7fbdb157_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Neg,
            exponent: -160,
            mantissa: 0xffb314a2_f59b431b_6f823871_9f1bb200_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Pos,
            exponent: -166,
            mantissa: 0x9682dfc7_fa8c0590_8c68114a_e7f5d5f6_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Pos,
            exponent: -167,
            mantissa: 0xb292fdda_b1feffd6_6b833ff9_3f28f989_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Neg,
            exponent: -174,
            mantissa: 0xd375c9be_3fb60327_c79f8310_4e085aac_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Neg,
            exponent: -175,
            mantissa: 0xbcf48234_2c4d5102_6c04e791_e1c42ce5_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Pos,
            exponent: -182,
            mantissa: 0xe04cc5a3_615422a5_b1bd4151_ce9f5ab7_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Pos,
            exponent: -183,
            mantissa: 0x9caaff2f_d886d69a_ae60bb50_315099ee_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Neg,
            exponent: -190,
            mantissa: 0xb9f644ba_cdd21c72_bbe5ecaf_9eda69e1_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Neg,
            exponent: -192,
            mantissa: 0xd105ebda_2aa49738_ab7f997a_9890da3e_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Pos,
            exponent: -199,
            mantissa: 0xf79cb165_ed01b717_af08b30b_d054e84d_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Pos,
            exponent: -201,
            mantissa: 0xe5328312_2c1279a2_7d0e101e_8a58acf9_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Neg,
            exponent: -207,
            mantissa: 0x87465718_e13fa89e_ec7a3b50_57b7140d_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Neg,
            exponent: -210,
            mantissa: 0xd230e3b5_5d1989f5_c680a514_473a081f_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Pos,
            exponent: -217,
            mantissa: 0xf6f36d05_82e7a851_81f38b32_6957ad79_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Pos,
            exponent: -219,
            mantissa: 0xa29d0634_e18100a9_56fbd2aa_570b0321_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Neg,
            exponent: -226,
            mantissa: 0xc1a0fc7f_c4a803a4_3ad8bbe9_583f1eff_u128,
        },
    ],
    [
        DyadicFloat128 {
            sign: DyadicSign::Pos,
            exponent: -232,
            mantissa: 0x870e7c8f_65e67d63_df27cb45_2e9e3147_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Pos,
            exponent: -131,
            mantissa: 0xe4b38001_8fb009e1_f19f1029_59fd9988_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Neg,
            exponent: -137,
            mantissa: 0x8f59503e_8784ef3e_63105adb_3bb1da75_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Neg,
            exponent: -133,
            mantissa: 0x9859b76e_d511a337_fca7d825_9121dbc2_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Pos,
            exponent: -140,
            mantissa: 0xbee96ef4_25b9dedf_950b4573_289e5dd6_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Pos,
            exponent: -138,
            mantissa: 0xf34b2a9f_05c526e0_c590f22d_1a630fd7_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Neg,
            exponent: -144,
            mantissa: 0x9842d51a_d0b236f4_143398d2_218f8f8f_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Neg,
            exponent: -143,
            mantissa: 0xb8ccf925_e6570536_d5433ab5_25e3b821_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Pos,
            exponent: -150,
            mantissa: 0xe6dd5b53_b30bb335_97a0c50d_7371371c_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Pos,
            exponent: -149,
            mantissa: 0xa3953e7f_dfc8f26d_6c6bdbab_d14160af_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Neg,
            exponent: -156,
            mantissa: 0xcbcfa908_e26d70c8_cffd7347_80f8bccf_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Neg,
            exponent: -156,
            mantissa: 0xbd5a580d_4816db67_507ad48a_0b44f91f_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Pos,
            exponent: -163,
            mantissa: 0xeb1c8272_6243df16_f882b0a3_24f02391_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Pos,
            exponent: -163,
            mantissa: 0x9a641fd8_e113c8d1_1336aeb5_47752b56_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Neg,
            exponent: -170,
            mantissa: 0xbee7ef89_f2e65eef_2abdfa95_88095f6f_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Neg,
            exponent: -171,
            mantissa: 0xbad80eed_67cc3f63_d316412c_77adbb9b_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Pos,
            exponent: -178,
            mantissa: 0xe5ea0f53_57e9369a_4c8e2811_9cb9e4a8_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Pos,
            exponent: -179,
            mantissa: 0xae6a3755_d0f4bcfd_bda090cf_9b8bd721_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Neg,
            exponent: -186,
            mantissa: 0xd56f6916_8e136501_806403f0_ce2de341_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Neg,
            exponent: -187,
            mantissa: 0x81607f59_664e6eb9_59c16aaa_df1f93b9_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Pos,
            exponent: -194,
            mantissa: 0x9d5877d8_3949f1cf_68aedada_eff31fc8_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Pos,
            exponent: -196,
            mantissa: 0x9c2b2573_922a73f1_ba64b4d9_a1100b9f_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Neg,
            exponent: -203,
            mantissa: 0xbca4a6cd_513cda13_be8f44e7_ac0fcc49_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Neg,
            exponent: -205,
            mantissa: 0x9c5b6330_d9e5cf58_0b9f6079_43e224b4_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Pos,
            exponent: -212,
            mantissa: 0xbb7afea5_52395e6a_d6b5151f_5e72496c_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Pos,
            exponent: -214,
            mantissa: 0x83ef08b6_66327b36_2d4cec7c_6db4fc40_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Neg,
            exponent: -221,
            mantissa: 0x9bf0f115_32c10238_d224759c_b26a4744_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Neg,
            exponent: -224,
            mantissa: 0xbd48ccc9_5d2a6cef_dd0c5777_eb4911b0_u128,
        },
    ],
    [
        DyadicFloat128 {
            sign: DyadicSign::Pos,
            exponent: -131,
            mantissa: 0xe14309a3_d8d9c46b_96213062_6b5e46a3_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Pos,
            exponent: -232,
            mantissa: 0xe675a86f_f6cf4f9a_a649ee45_c5ce54f3_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Neg,
            exponent: -132,
            mantissa: 0xe14309a3_d8d9c46b_96213062_6b815115_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Pos,
            exponent: -139,
            mantissa: 0xb6abdffd_9bfea8f9_00d7117f_f813be8a_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Pos,
            exponent: -135,
            mantissa: 0x9603077a_a9f8bd29_626a7a0b_c7de5368_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Neg,
            exponent: -142,
            mantissa: 0x91d203db_916971df_91b90d34_53276618_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Neg,
            exponent: -140,
            mantissa: 0x9faaac0f_6341eb46_865a401a_ba26f646_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Pos,
            exponent: -147,
            mantissa: 0xa5fda6d5_17fb0c3f_68a08b44_9fbbd69b_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Pos,
            exponent: -146,
            mantissa: 0xb5e3484c_d797b76a_eb0dcdff_9039d89e_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Neg,
            exponent: -153,
            mantissa: 0xc399a04c_3e93b36a_4e2fc63c_f858c0bc_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Neg,
            exponent: -152,
            mantissa: 0x80ca8f81_6959b86e_5ce0cf60_93ec45bf_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Pos,
            exponent: -159,
            mantissa: 0x8d30ba7f_2afdf1c2_da9623b3_9af68d65_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Pos,
            exponent: -160,
            mantissa: 0xf87652de_01424ed9_a8001c3f_2de41b6d_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Neg,
            exponent: -166,
            mantissa: 0x89be58fe_62ad19b8_e7e494f2_44ac466c_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Neg,
            exponent: -167,
            mantissa: 0xada7bf95_d8cf67ec_83e56db0_b540a4be_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Pos,
            exponent: -174,
            mantissa: 0xc1c9a76b_cbd5e19d_a79a9d54_cd9872dd_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Pos,
            exponent: -175,
            mantissa: 0xb7e94fe5_e79c4e61_54d61abf_8977f4f6_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Neg,
            exponent: -182,
            mantissa: 0xcddf251c_35908de7_2715ad2f_75e2f40b_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Neg,
            exponent: -183,
            mantissa: 0x98a1dc7b_7d3c7b42_c1ac6b19_3dfda31c_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Pos,
            exponent: -190,
            mantissa: 0xaaf8f7c6_3298656f_eb0fa14e_e7aff52c_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Pos,
            exponent: -192,
            mantissa: 0xcbd904d3_cfcf950b_8ffb5b41_71781854_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Neg,
            exponent: -199,
            mantissa: 0xe41151ec_c34e3f5c_5aeba01d_68094872_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Neg,
            exponent: -201,
            mantissa: 0xdfc36cae_64788d40_a88d2ac2_c906458c_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Pos,
            exponent: -208,
            mantissa: 0xf9ac8be8_59fbacbf_5d220fac_15791266_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Pos,
            exponent: -210,
            mantissa: 0xcd6fc387_e853da44_08be1c48_6c1aef1d_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Neg,
            exponent: -217,
            mantissa: 0xe45a5994_c8db290b_db6d5527_46a20ad4_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Neg,
            exponent: -219,
            mantissa: 0x9f1d6b57_d59161a5_3a583f64_1f64629c_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Pos,
            exponent: -226,
            mantissa: 0xb36453a7_954808a7_68bd4abc_a98c08a8_u128,
        },
    ],
    [
        DyadicFloat128 {
            sign: DyadicSign::Pos,
            exponent: -233,
            mantissa: 0x8bf58935_5d7c3898_ff2d5099_c4edd107_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Neg,
            exponent: -131,
            mantissa: 0xddf92300_ce06a066_1be2d50f_c9e71bfb_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Pos,
            exponent: -137,
            mantissa: 0x831120cc_8a06bcfb_78472e09_d96d198c_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Pos,
            exponent: -133,
            mantissa: 0x93e1a0b6_04c4bdc6_ea961f1e_42f661a7_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Neg,
            exponent: -140,
            mantissa: 0xae93cee4_3e7b45d7_fecde400_21f10330_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Neg,
            exponent: -138,
            mantissa: 0xec352fad_6ded6029_3ebb76e5_40c9fe1f_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Pos,
            exponent: -144,
            mantissa: 0x8b4820b1_a6124383_4269a8f1_8434f7e1_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Pos,
            exponent: -143,
            mantissa: 0xb37ade8d_f6d70a11_7a9bed16_8b3ac96e_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Neg,
            exponent: -150,
            mantissa: 0xd34d9638_beaa4d97_e876da0b_dfb1db97_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Neg,
            exponent: -149,
            mantissa: 0x9ef2755a_2b65e64c_d89a7f7c_47ef8100_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Pos,
            exponent: -156,
            mantissa: 0xbaaf55de_f348a7cb_5b5c92a4_ae40205f_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Pos,
            exponent: -156,
            mantissa: 0xb817f3dd_4fc61422_2c0100aa_825c360b_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Neg,
            exponent: -163,
            mantissa: 0xd79003f7_5712b0c0_69df67bb_a40cf475_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Neg,
            exponent: -163,
            mantissa: 0x9634ba34_a472a406_06a7d4c9_503b9bbf_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Pos,
            exponent: -170,
            mantissa: 0xaf3bd6cd_4ed78331_c05e6d1a_fe7c3468_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Pos,
            exponent: -171,
            mantissa: 0xb5ebde30_22e67cc0_5dbd042d_d51cabf5_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Neg,
            exponent: -178,
            mantissa: 0xd3515100_51960306_5d5a670f_30b9b9b9_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Neg,
            exponent: -179,
            mantissa: 0xa9f78000_d47e22ba_1b00a5e5_14af6046_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Pos,
            exponent: -186,
            mantissa: 0xc4759bfe_0e7e9793_833b856c_d9303018_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Pos,
            exponent: -188,
            mantissa: 0xfc64197a_8da553f3_512dbe10_16144d17_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Neg,
            exponent: -194,
            mantissa: 0x91100501_b02297fe_b4e73df6_a47dc3f3_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Neg,
            exponent: -196,
            mantissa: 0x987acb6d_feed70cb_cbfe5527_02367260_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Pos,
            exponent: -203,
            mantissa: 0xae36e80b_0186684b_48e9d787_15814118_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Pos,
            exponent: -205,
            mantissa: 0x98d2753b_4e3e97fd_63ef454d_7eda6d09_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Neg,
            exponent: -212,
            mantissa: 0xad73de91_dc840643_1caef383_494f10cb_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Neg,
            exponent: -214,
            mantissa: 0x8116b94c_069e3e04_80429da9_abea48fe_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Pos,
            exponent: -221,
            mantissa: 0x908a9ab1_ad0ab774_ae56464a_75bff858_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Pos,
            exponent: -224,
            mantissa: 0xb960db4f_6d3cf5cd_41733ab5_1943ab37_u128,
        },
    ],
    [
        DyadicFloat128 {
            sign: DyadicSign::Neg,
            exponent: -131,
            mantissa: 0xdad3110c_c5395e68_2d46d7e5_437df41f_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Neg,
            exponent: -234,
            mantissa: 0xa89f4f88_e42252f2_7e38d053_34f92f4b_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Pos,
            exponent: -132,
            mantissa: 0xdad3110c_c5395e68_2d46d7e5_43840072_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Neg,
            exponent: -139,
            mantissa: 0xa7738e7b_3e1215ac_051a12dc_cb463305_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Neg,
            exponent: -135,
            mantissa: 0x91be015a_317b7782_58705d87_22727559_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Pos,
            exponent: -142,
            mantissa: 0x85b3db7e_ff630cb9_a0f441de_6a7b67ae_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Pos,
            exponent: -140,
            mantissa: 0x9b28b3ec_5e197c95_8bbe569a_12862227_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Neg,
            exponent: -147,
            mantissa: 0x98435bd1_4b734bbe_4f4f9c1d_8b6c9ca6_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Neg,
            exponent: -146,
            mantissa: 0xb0d0a0a1_2e153299_c6189117_3900843c_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Pos,
            exponent: -153,
            mantissa: 0xb389591c_c8f3fd4b_4f37fe61_950d9bfa_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Pos,
            exponent: -153,
            mantissa: 0xfa839b41_cf240b81_0735a89d_3b53e995_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Neg,
            exponent: -159,
            mantissa: 0x81b36815_c8dc1e6b_2d77b369_0841bc94_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Neg,
            exponent: -160,
            mantissa: 0xf1c7cded_2fb6f81b_7d3083bf_f34dc3f0_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Pos,
            exponent: -167,
            mantissa: 0xfd51e8e5_80f329f4_7f210ca9_10fd2938_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Pos,
            exponent: -167,
            mantissa: 0xa918bd15_2dcd0c35_853ca255_195fd6a9_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Neg,
            exponent: -174,
            mantissa: 0xb2668342_3de8e53d_791b92a2_ce4bac0a_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Neg,
            exponent: -175,
            mantissa: 0xb33775ff_5ab9e1a8_8cfd6fbb_36cd1c27_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Pos,
            exponent: -182,
            mantissa: 0xbdc5aa78_bc2e7dd5_4c75d1dd_fee83aa4_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Pos,
            exponent: -183,
            mantissa: 0x94db8651_9769e690_ceab0654_97ba7cf0_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Neg,
            exponent: -190,
            mantissa: 0x9dd4247b_c9d37bcd_dff464ba_6dbd0af9_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Neg,
            exponent: -192,
            mantissa: 0xc6fb2a91_1fc022c5_df8fb4ac_08dfebc3_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Pos,
            exponent: -199,
            mantissa: 0xd2dc803f_a7a3775f_7e996204_aa05f135_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Pos,
            exponent: -201,
            mantissa: 0xda9fd0d4_e8575bcc_0b4d690e_7f0354d7_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Neg,
            exponent: -208,
            mantissa: 0xe737e73d_462942b1_dbfa9f4b_76dc7d1f_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Neg,
            exponent: -210,
            mantissa: 0xc8e9a429_d5db3a79_2ec42a98_18a18722_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Pos,
            exponent: -217,
            mantissa: 0xd3d6b7a6_efa88dd9_bc22e8c5_617bd6a3_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Pos,
            exponent: -219,
            mantissa: 0x9bc3c9e1_900f097e_2ebffed7_e35ff774_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Neg,
            exponent: -226,
            mantissa: 0xa6b2bae0_5e1ce794_087572eb_3823a0eb_u128,
        },
    ],
    [
        DyadicFloat128 {
            sign: DyadicSign::Pos,
            exponent: -232,
            mantissa: 0xe6517958_6a11a45f_9e6dbd98_a13caa56_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Pos,
            exponent: -131,
            mantissa: 0xd7ce5a16_68453adc_b67aca94_7ba74333_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Neg,
            exponent: -138,
            mantissa: 0xf0e336bb_0b71b6b2_0dd33ee9_8a2fcd4a_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Neg,
            exponent: -133,
            mantissa: 0x8fc87e88_09ac45bf_2e66443a_2668af75_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Pos,
            exponent: -140,
            mantissa: 0xa071f584_dbeb666b_d95a04a4_e9fd0f54_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Pos,
            exponent: -138,
            mantissa: 0xe5b417fa_3893a0fa_098f6d81_f20b5b84_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Neg,
            exponent: -144,
            mantissa: 0x800b40cd_30eebc6c_656611e9_da93fbb0_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Neg,
            exponent: -143,
            mantissa: 0xae968e21_032366cd_a2964ea6_7f8f8502_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Pos,
            exponent: -150,
            mantissa: 0xc2587f45_429aa6a8_575eaa34_555d0382_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Pos,
            exponent: -149,
            mantissa: 0x9aacecb5_f7acc0ea_ddf136cc_05b869a3_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Neg,
            exponent: -156,
            mantissa: 0xabd05f3a_4cbd6584_b54652a4_f84f8cf3_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Neg,
            exponent: -156,
            mantissa: 0xb33be47a_9db530e4_0bbbb54b_6ff3c21a_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Pos,
            exponent: -163,
            mantissa: 0xc68d710d_bab50675_c3e9f375_5de6f54f_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Pos,
            exponent: -163,
            mantissa: 0x925370d6_32f43543_7c7f4999_88765112_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Neg,
            exponent: -170,
            mantissa: 0xa1900b81_a5ca5e74_d7204c92_ac472465_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Neg,
            exponent: -171,
            mantissa: 0xb1570223_e3e44086_05c6b920_7d856fc6_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Pos,
            exponent: -178,
            mantissa: 0xc30c8735_4384143f_7b1d8557_18eff201_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Pos,
            exponent: -179,
            mantissa: 0xa5cf0b11_7cac2ef1_059a3baf_8993dda3_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Neg,
            exponent: -186,
            mantissa: 0xb58f9f28_f3998f23_58b2ff5c_a328d5c7_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Neg,
            exponent: -188,
            mantissa: 0xf66a098a_a14fcccf_ebe96e2b_a6e5776e_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Pos,
            exponent: -194,
            mantissa: 0x863eb76e_83d8bb4c_2462238a_bfc7febe_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Pos,
            exponent: -196,
            mantissa: 0x94ff0efd_6f040a79_c9e3fd65_d78ed1d6_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Neg,
            exponent: -203,
            mantissa: 0xa175ab8a_87cf20b6_c230c5af_28b01250_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Neg,
            exponent: -205,
            mantissa: 0x95771bb8_fa846964_74c8a6dc_fb6b7ec9_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Pos,
            exponent: -212,
            mantissa: 0xa100da06_419ef6cc_9b85936b_b50dda28_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Pos,
            exponent: -215,
            mantissa: 0xfcbde450_a218fc29_4e41d297_193a08c6_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Neg,
            exponent: -221,
            mantissa: 0x86622859_a5f7d9ee_477a8ff3_708ef92e_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Neg,
            exponent: -224,
            mantissa: 0xb5a0815d_3690fc92_8e249b1e_aa2945c5_u128,
        },
    ],
    [
        DyadicFloat128 {
            sign: DyadicSign::Pos,
            exponent: -131,
            mantissa: 0xd4e8c1de_022a299a_6abbfd83_9f7b8c3a_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Pos,
            exponent: -233,
            mantissa: 0xd3b939d1_7be585d4_bd784c91_4ad6fbee_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Neg,
            exponent: -132,
            mantissa: 0xd4e8c1de_022a299a_6abbfd83_9f89ed20_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Pos,
            exponent: -139,
            mantissa: 0x9a3bcc3a_0086e02e_04bdc699_68d0a692_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Pos,
            exponent: -135,
            mantissa: 0x8dd114cd_5c453ee0_10403b1b_164c5d5d_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Neg,
            exponent: -143,
            mantissa: 0xf658cfad_ee32482b_81c54324_bba4ab53_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Neg,
            exponent: -140,
            mantissa: 0x97029383_e642cb9c_eba864b0_498b98a5_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Pos,
            exponent: -147,
            mantissa: 0x8c52ed3b_f6ec14b3_88b2fa53_b698e05a_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Pos,
            exponent: -146,
            mantissa: 0xac235b11_7ac95445_fe557d86_39196589_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Neg,
            exponent: -153,
            mantissa: 0xa58bf4bd_ac42515e_872e9afa_840e4c65_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Neg,
            exponent: -153,
            mantissa: 0xf3fbb33a_6fd557a5_7779541b_d83817ed_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Pos,
            exponent: -160,
            mantissa: 0xef5a8276_9a9e2d69_97afdccc_efd81e57_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Pos,
            exponent: -160,
            mantissa: 0xeb9712f8_c8204ace_eb04cd07_b18955b2_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Neg,
            exponent: -167,
            mantissa: 0xe9f03b58_100eb10d_222c75ec_8a54ecb0_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Neg,
            exponent: -167,
            mantissa: 0xa4dc1298_157c3fe9_db7b1d55_08707143_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Pos,
            exponent: -174,
            mantissa: 0xa4e9aa3a_fe2e23b0_bf3d2ebc_b7535062_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Pos,
            exponent: -175,
            mantissa: 0xaed64b11_fbd9cc6b_f66f10a1_233c1de2_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Neg,
            exponent: -182,
            mantissa: 0xaf9ea898_6e2134fe_d2dcde87_a27f3d7f_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Neg,
            exponent: -183,
            mantissa: 0x9152575b_4baed63c_fc75970a_c083ad3a_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Pos,
            exponent: -190,
            mantissa: 0x923cf1d6_84827930_f7db32e8_83144243_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Pos,
            exponent: -192,
            mantissa: 0xc266dea4_8eda10ea_9aede360_f9d48beb_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Neg,
            exponent: -199,
            mantissa: 0xc3a2f6b4_75d84c10_bae52a30_e8465e42_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Neg,
            exponent: -201,
            mantissa: 0xd5c3d258_799081b9_7e422517_d32cb4f6_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Pos,
            exponent: -208,
            mantissa: 0xd6d43a9d_a797d3c8_2a095247_82c25598_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Pos,
            exponent: -210,
            mantissa: 0xc49cd768_75cab1de_d58316b4_b3b534ca_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Neg,
            exponent: -217,
            mantissa: 0xc51e2838_05db1f6d_9ab63708_3f7a7072_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Neg,
            exponent: -219,
            mantissa: 0x98902525_147a198c_ff6f49b7_cc7190a5_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Pos,
            exponent: -226,
            mantissa: 0x9b55eca6_6d4f528a_d9c385df_d6e53221_u128,
        },
    ],
    [
        DyadicFloat128 {
            sign: DyadicSign::Neg,
            exponent: -234,
            mantissa: 0xff798a37_c5e4616a_eec01447_ebe85e3b_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Neg,
            exponent: -131,
            mantissa: 0xd2203ed6_3914b51a_c1791b73_da45fabc_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Pos,
            exponent: -138,
            mantissa: 0xde5d26f6_0ed222ae_47b28b9e_f05b3664_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Pos,
            exponent: -133,
            mantissa: 0x8c01e32e_57efd403_1154de50_4a678fb3_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Neg,
            exponent: -140,
            mantissa: 0x941ef959_f2536cde_60b50489_65838d68_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Neg,
            exponent: -138,
            mantissa: 0xdfb4af44_92ce6760_e0fe3707_0f1c8063_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Pos,
            exponent: -145,
            mantissa: 0xec79bac8_916a8930_f86ed6a9_31539cd5_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Pos,
            exponent: -143,
            mantissa: 0xaa123647_0210c206_1f6bfa1a_29d613f9_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Neg,
            exponent: -150,
            mantissa: 0xb3888de4_733d21b2_983b9c2c_401cc402_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Neg,
            exponent: -149,
            mantissa: 0x96b950fb_1760c347_9de56eba_de3716d2_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Pos,
            exponent: -156,
            mantissa: 0x9ece51f5_8f117b7a_2cae04c4_9ec88fd7_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Pos,
            exponent: -156,
            mantissa: 0xaeba5371_12a8a47f_ca2661c8_a8aea46a_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Neg,
            exponent: -163,
            mantissa: 0xb7a5a105_97cabb3c_16c8ad06_eff54900_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Neg,
            exponent: -163,
            mantissa: 0x8eb7cc2e_e1f64b9b_f8486207_a0f5842b_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Pos,
            exponent: -170,
            mantissa: 0x958ec088_0ef542f7_e62f629e_6a40a667_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Pos,
            exponent: -171,
            mantissa: 0xad10cb84_4fdfd4f9_d387ff50_6b2b050e_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Neg,
            exponent: -178,
            mantissa: 0xb4ba57f8_be3c6e38_944c6331_ecc96a83_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Neg,
            exponent: -179,
            mantissa: 0xa1ea4497_50de74ce_29ebfc20_828d725c_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Pos,
            exponent: -186,
            mantissa: 0xa8690da2_ffc0645d_43fc7fb8_e774b2b4_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Pos,
            exponent: -188,
            mantissa: 0xf0cb5a04_395987fd_b7c1c476_12198ccc_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Neg,
            exponent: -195,
            mantissa: 0xf955ea6d_aecdec7f_6bfbafea_86b29b06_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Neg,
            exponent: -196,
            mantissa: 0x91b4ceb3_35ff632b_3846ac35_49092372_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Pos,
            exponent: -203,
            mantissa: 0x962189f8_cec85dbb_008a31e2_4dd0e173_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Pos,
            exponent: -205,
            mantissa: 0x92478549_4259bade_c7321e0b_2961da49_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Neg,
            exponent: -212,
            mantissa: 0x95e88dfb_46ad8700_0a551237_3f5ace2c_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Neg,
            exponent: -215,
            mantissa: 0xf78e6965_b9cdbc91_4ef8b483_c93b94c8_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Pos,
            exponent: -222,
            mantissa: 0xfa98e54c_345074b1_03af10b5_f56b3b80_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Pos,
            exponent: -224,
            mantissa: 0xb20836eb_2b4b0b9d_4713dce5_a4035f22_u128,
        },
    ],
    [
        DyadicFloat128 {
            sign: DyadicSign::Neg,
            exponent: -131,
            mantissa: 0xcf72f770_d4be2e0a_256ad93e_a559b1ed_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Neg,
            exponent: -232,
            mantissa: 0xa017b3d0_5ac1a11d_3382c127_6a42faef_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Pos,
            exponent: -132,
            mantissa: 0xcf72f770_d4be2e0a_256ad93e_a56e569b_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Neg,
            exponent: -139,
            mantissa: 0x8eaaec9f_2cdeb933_10091641_45e66fce_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Neg,
            exponent: -135,
            mantissa: 0x8a310c9c_a54ce00b_443f4a26_6802b3d3_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Pos,
            exponent: -143,
            mantissa: 0xe3e995db_237151ca_8e50931e_16c464c8_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Pos,
            exponent: -140,
            mantissa: 0x932ccd76_a78e5002_b2ac3089_82f666b7_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Neg,
            exponent: -147,
            mantissa: 0x81dd2c76_42e08f83_e6730dd6_cd1649a0_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Neg,
            exponent: -146,
            mantissa: 0xa7cf1cb0_ddceea04_d0ddb3db_a7b1f74a_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Pos,
            exponent: -153,
            mantissa: 0x9946a6e1_9a03c9ea_9ee2c6a2_042d094e_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Pos,
            exponent: -153,
            mantissa: 0xeded3d63_87fc960c_5ac532a3_92b9148b_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Neg,
            exponent: -160,
            mantissa: 0xddbe178c_38a26e85_7b622da0_1f114832_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Neg,
            exponent: -160,
            mantissa: 0xe5d61263_b02b38a6_8035475c_e054a207_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Pos,
            exponent: -167,
            mantissa: 0xd8e1adb2_f65272dd_eedcc349_06166244_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Pos,
            exponent: -167,
            mantissa: 0xa0e91fc9_ef3cfc4a_755362f0_efa2b443_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Neg,
            exponent: -174,
            mantissa: 0x99046c03_34c60cca_85dc9725_032e0e06_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Neg,
            exponent: -175,
            mantissa: 0xaabe0b31_286cdd1c_ea955cd1_29b75961_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Pos,
            exponent: -182,
            mantissa: 0xa31b8985_4ebfbc53_ad577889_827fa263_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Pos,
            exponent: -183,
            mantissa: 0x8e010c9a_52cf8ac3_298eb52f_a7c956b9_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Neg,
            exponent: -190,
            mantissa: 0x87f688ce_4b9509fd_63bcb1b7_5d72d5d8_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Neg,
            exponent: -192,
            mantissa: 0xbe16b228_eb60838a_83628b71_604fa611_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Pos,
            exponent: -199,
            mantissa: 0xb619abe4_eb82bbb7_81286fd2_4bcdae97_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Pos,
            exponent: -201,
            mantissa: 0xd12b3250_56348b6a_e1da2db5_5856a775_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Neg,
            exponent: -208,
            mantissa: 0xc8362dcb_30d4a82e_c3f69033_99411021_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Neg,
            exponent: -210,
            mantissa: 0xc086fa2f_a40093a3_19aae864_9f62c454_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Pos,
            exponent: -217,
            mantissa: 0xb7f1fb33_947e59f5_b739f649_5c9fa9f7_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Pos,
            exponent: -219,
            mantissa: 0x9581bba9_910c6d4d_1617cabc_14d90336_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Neg,
            exponent: -226,
            mantissa: 0x91249f9c_f4a1fa1e_3f52351f_67b74f0f_u128,
        },
    ],
    [
        DyadicFloat128 {
            sign: DyadicSign::Neg,
            exponent: -234,
            mantissa: 0x80b6dced_8f5be326_a784b826_23129675_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Pos,
            exponent: -131,
            mantissa: 0xccdf39fd_1f7e6122_d28ddb06_63ee5237_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Neg,
            exponent: -138,
            mantissa: 0xce19246d_3abba77c_cee7448b_a694959a_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Neg,
            exponent: -133,
            mantissa: 0x88838a3f_6d3ffeb1_982a290f_e08a978f_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Pos,
            exponent: -140,
            mantissa: 0x894c05f4_ec28fe60_7a2932a5_da7036e7_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Pos,
            exponent: -138,
            mantissa: 0xda270910_f3028749_405ed567_d2595bd2_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Neg,
            exponent: -145,
            mantissa: 0xdb3dbc60_eece2abf_f9d68123_f3017efd_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Neg,
            exponent: -143,
            mantissa: 0xa5e24f4d_ade7eaa4_c8e8e9b7_8c674d56_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Pos,
            exponent: -150,
            mantissa: 0xa6817fdc_5e7e7f41_3486d723_5922aac2_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Pos,
            exponent: -149,
            mantissa: 0x930e18a2_b721e3c6_5a4ec8d6_9a4fbba9_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Neg,
            exponent: -156,
            mantissa: 0x935a003a_26f78aec_2333db47_99615b20_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Neg,
            exponent: -156,
            mantissa: 0xaa8928d1_5a111455_10e22487_bd6b5f45_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Pos,
            exponent: -163,
            mantissa: 0xaa807ca9_8889696f_e837cede_ee2df6f7_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Pos,
            exponent: -163,
            mantissa: 0x8b5a749b_7c54b018_3f723e08_3d02639a_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Neg,
            exponent: -170,
            mantissa: 0x8af38594_0c171e6a_40728355_04b91243_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Neg,
            exponent: -171,
            mantissa: 0xa9118601_e594fa61_b5ffe70e_e0541b5b_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Pos,
            exponent: -178,
            mantissa: 0xa80c84ac_3c9498af_3b93c831_72618dff_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Pos,
            exponent: -179,
            mantissa: 0x9e4326ec_5ec3ba92_4029f393_33a63072_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Neg,
            exponent: -186,
            mantissa: 0x9cbd7baf_190c608f_e8a1b338_ac4b1223_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Neg,
            exponent: -188,
            mantissa: 0xeb80e4c9_be4ad6ca_e1e977e8_32661b7f_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Pos,
            exponent: -195,
            mantissa: 0xe84b2ec9_302a6542_d866fce4_e3c73ed7_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Pos,
            exponent: -196,
            mantissa: 0x8e98ce27_d0a17ced_068cc466_a228fa46_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Neg,
            exponent: -203,
            mantissa: 0x8c060f04_1fc26b81_12b762ce_f531c1f9_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Neg,
            exponent: -205,
            mantissa: 0x8f41829b_b0632c31_081e2454_d21f2b36_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Pos,
            exponent: -212,
            mantissa: 0x8bfaf48d_7cb2cb69_ae9e179c_0b374718_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Pos,
            exponent: -215,
            mantissa: 0xf29d049a_d0fa1b07_1f7d245d_7dc92d82_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Neg,
            exponent: -222,
            mantissa: 0xea49ccc1_ef00ff59_ed3b3d78_ea936fbf_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Neg,
            exponent: -224,
            mantissa: 0xae976cca_82f23560_4560952e_741b1674_u128,
        },
    ],
    [
        DyadicFloat128 {
            sign: DyadicSign::Pos,
            exponent: -131,
            mantissa: 0xca637ac4_c3845548_87ace5f3_a666195e_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Pos,
            exponent: -232,
            mantissa: 0xfb905625_02c85699_948d7103_bf250cea_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Neg,
            exponent: -132,
            mantissa: 0xca637ac4_c3845548_87ace5f3_a684f94e_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Pos,
            exponent: -139,
            mantissa: 0x847a78dc_7d86de6e_04860dd9_d21f571b_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Pos,
            exponent: -135,
            mantissa: 0x86d498e3_bb638a24_e9ab8515_8c3c4cd9_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Neg,
            exponent: -143,
            mantissa: 0xd3aa8694_a0794fe9_710c3011_1c3b5fe3_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Neg,
            exponent: -140,
            mantissa: 0x8f9dcbc5_23e79e08_bb9d2029_bec9c1b3_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Pos,
            exponent: -148,
            mantissa: 0xf1475012_396f9d2b_57bd9228_4b44c9aa_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Pos,
            exponent: -146,
            mantissa: 0xa3c98576_3c65bc6f_e5d186a7_9c8f424c_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Neg,
            exponent: -153,
            mantissa: 0x8e719627_9a728707_d4290eee_6eb3518f_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Neg,
            exponent: -153,
            mantissa: 0xe84a874c_e459c497_4a46ed18_561df4e6_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Pos,
            exponent: -160,
            mantissa: 0xce2d0b01_f3cf3e59_9e82a8d2_3e5483c6_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Pos,
            exponent: -160,
            mantissa: 0xe078bbc6_7cc4343a_200f3e1d_7dd6f94c_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Neg,
            exponent: -167,
            mantissa: 0xc9c7f098_4bb482d1_e127bfa9_9ac64fcb_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Neg,
            exponent: -167,
            mantissa: 0x9d385ed7_6d5f7bd9_df539347_504d12b9_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Pos,
            exponent: -174,
            mantissa: 0x8e777694_54dde3db_02c84807_a3b73987_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Pos,
            exponent: -175,
            mantissa: 0xa6e7cae2_e54a4afb_495a6c7d_f920fc55_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Neg,
            exponent: -182,
            mantissa: 0x97fc938b_b2be4147_7913a3a5_5c021ec2_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Neg,
            exponent: -183,
            mantissa: 0x8ae2d1f0_122175e8_a8c96ed9_b73ad454_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Pos,
            exponent: -191,
            mantissa: 0xfd9e4b46_8b9233b0_2f4599f5_1fa5b398_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Pos,
            exponent: -192,
            mantissa: 0xba0571d4_69ff1632_fef5dfbe_79bff398_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Neg,
            exponent: -199,
            mantissa: 0xaa02a7d3_062a684e_e61436ec_7f402d94_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Neg,
            exponent: -201,
            mantissa: 0xccd19917_f60b9c8d_b40ea505_8097b58e_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Pos,
            exponent: -208,
            mantissa: 0xbb1ee000_adee529f_aad16d5c_95f85854_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Pos,
            exponent: -210,
            mantissa: 0xbca5489c_25c60496_7103b3db_5c87adc6_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Neg,
            exponent: -217,
            mantissa: 0xac1d3afe_9b849e72_a7f9c00a_44eb73ae_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Neg,
            exponent: -219,
            mantissa: 0x92974c40_a70d1034_58c47fb9_a1dd6194_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Pos,
            exponent: -226,
            mantissa: 0x87f58751_7565b27d_0dadad9e_8ff25822_u128,
        },
    ],
    [
        DyadicFloat128 {
            sign: DyadicSign::Neg,
            exponent: -234,
            mantissa: 0xeb17bb29_488f7a3a_b5572360_e69f2bca_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Neg,
            exponent: -131,
            mantissa: 0xc7fe4dcd_098fb000_8eb84d9f_354f4472_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Pos,
            exponent: -138,
            mantissa: 0xbfb9218a_5dcce9b4_dd57b364_9034d8ae_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Pos,
            exponent: -133,
            mantissa: 0x8544e2e8_83a5cdef_26878d26_8ae364e3_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Neg,
            exponent: -141,
            mantissa: 0xff757599_a32f3ce4_b90c4d65_65b5c54c_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Neg,
            exponent: -138,
            mantissa: 0xd4fdd1b5_3d0473e8_67289834_eb21de6f_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Pos,
            exponent: -145,
            mantissa: 0xcc001a3b_879d8c1b_52848ce4_634661d2_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Pos,
            exponent: -143,
            mantissa: 0xa1fd25dc_65d2ba36_991f5d55_c1bd453b_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Neg,
            exponent: -150,
            mantissa: 0x9afa079d_1c815819_48481030_ddf353f3_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Neg,
            exponent: -149,
            mantissa: 0x8fa32dd1_4960643a_2c4d7eb3_12214f4a_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Pos,
            exponent: -156,
            mantissa: 0x893447ba_7c082f9c_7bbf4262_34a573b9_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Pos,
            exponent: -156,
            mantissa: 0xa69fbf1f_eef4aa2e_0600e272_ba7878ca_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Neg,
            exponent: -163,
            mantissa: 0x9ed774e1_c63a52cf_ce9a5440_69ad9513_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Neg,
            exponent: -163,
            mantissa: 0x88350997_734d35ca_2a1f2ee1_fa774a7a_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Pos,
            exponent: -170,
            mantissa: 0x81873bb2_eb63ab7c_1e979461_b110b37b_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Pos,
            exponent: -171,
            mantissa: 0xa5526120_38562114_fc8399d6_3e9a694a_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Neg,
            exponent: -178,
            mantissa: 0x9cc3a013_895af7b1_df3ab96d_a3575751_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Neg,
            exponent: -179,
            mantissa: 0x9ad43cca_7cf9112d_ee4f35f6_15a8d734_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Pos,
            exponent: -186,
            mantissa: 0x925509c7_b921b971_bf8bdfcb_f63a88cb_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Pos,
            exponent: -188,
            mantissa: 0xe683f5bd_25ae465a_949792a5_b42077fa_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Neg,
            exponent: -195,
            mantissa: 0xd90f96ad_ddfbc2b1_837fba69_1c62ba0c_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Neg,
            exponent: -196,
            mantissa: 0x8ba7dc09_06bfd155_0a5585a2_6880ff8b_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Pos,
            exponent: -203,
            mantissa: 0x82f7a636_f0a040a3_d86fe47e_0416d177_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Pos,
            exponent: -205,
            mantissa: 0x8c62ba2a_d6174c8e_fc91282f_8c0950c3_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Neg,
            exponent: -212,
            mantissa: 0x830fc37a_d0fe0268_cb923bac_e673c5d9_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Neg,
            exponent: -215,
            mantissa: 0xede721aa_27269e71_8be6f574_ef452a82_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Pos,
            exponent: -222,
            mantissa: 0xdb990c19_e945d72e_f0a98aa8_e060c20a_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Pos,
            exponent: -224,
            mantissa: 0xab4d1f16_82e9ba12_909e283a_58de1237_u128,
        },
    ],
    [
        DyadicFloat128 {
            sign: DyadicSign::Neg,
            exponent: -131,
            mantissa: 0xc5ae6581_ea2cd4fb_17d79504_95a3fb2e_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Neg,
            exponent: -232,
            mantissa: 0xf8819c1f_53b0f6dd_4749f308_d6fd7b19_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Pos,
            exponent: -132,
            mantissa: 0xc5ae6581_ea2cd4fb_17d79504_95c11402_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Neg,
            exponent: -140,
            mantissa: 0xf6e4c5c9_7f6f94fa_c8162a0c_602625bc_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Neg,
            exponent: -135,
            mantissa: 0x83b3ea92_67543e57_6365540d_7c1524f5_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Pos,
            exponent: -143,
            mantissa: 0xc542d423_e44e1492_a2b7860f_6a94441f_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Pos,
            exponent: -140,
            mantissa: 0x8c4d7d51_5ffe6044_9b849530_6fe2c3c0_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Neg,
            exponent: -148,
            mantissa: 0xe0e918d3_4dd8109c_e2e39715_a65ac880_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Neg,
            exponent: -146,
            mantissa: 0xa009cded_e4e9ac2e_c3200594_186023e7_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Pos,
            exponent: -153,
            mantissa: 0x84d342ce_471b9fa6_999ff138_40afa156_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Pos,
            exponent: -153,
            mantissa: 0xe307e2fc_da82f8f6_a2305823_db7555bf_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Neg,
            exponent: -160,
            mantissa: 0xc0569c8f_b912621e_9748c93a_7e074f0e_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Neg,
            exponent: -160,
            mantissa: 0xdb74aa6e_b649122e_ba397aa1_f4587f46_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Pos,
            exponent: -167,
            mantissa: 0xbc56e3b6_50e8d58b_60b85656_af461d18_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Pos,
            exponent: -167,
            mantissa: 0x99c33e41_279be042_ce13520f_7f072dad_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Neg,
            exponent: -174,
            mantissa: 0x850f62fa_6e54813f_2721a70f_beae02a1_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Neg,
            exponent: -175,
            mantissa: 0xa34d6439_5b47e3d1_206f240d_63d4d4d9_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Pos,
            exponent: -182,
            mantissa: 0x8e0db49b_f88fa5f0_72913b4b_c7a1d997_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Pos,
            exponent: -183,
            mantissa: 0x87f34381_9e6e26e8_8a5f081f_09c895ea_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Neg,
            exponent: -191,
            mantissa: 0xed3ba2cc_76e4bbff_16616ff7_46c77c56_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Neg,
            exponent: -192,
            mantissa: 0xb62e3c85_9b699ea3_7fc606c4_a23327ad_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Pos,
            exponent: -199,
            mantissa: 0x9f2a7c67_d95a31ec_e04cbe1e_4a880525_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Pos,
            exponent: -201,
            mantissa: 0xc8b2c293_04fb7f3e_216c218c_d8966163_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Neg,
            exponent: -208,
            mantissa: 0xaf59a9e8_6b8aee1a_b8d1cb0a_d709fd6a_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Neg,
            exponent: -210,
            mantissa: 0xb8f4d521_cf5c6421_86dfdbb7_fb9aa406_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Pos,
            exponent: -217,
            mantissa: 0xa1732700_3c4d595b_d48680a8_1b6ccfd4_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Pos,
            exponent: -219,
            mantissa: 0x8fcf506f_756c4df8_a987fad2_6e657d0e_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Neg,
            exponent: -227,
            mantissa: 0xff511190_93360e5e_4f88ddcf_630c91fb_u128,
        },
    ],
    [
        DyadicFloat128 {
            sign: DyadicSign::Neg,
            exponent: -232,
            mantissa: 0xc93dfd46_da958f03_7c85f4d9_1cfd732a_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Pos,
            exponent: -131,
            mantissa: 0xc3728dd9_77126822_c808c6b6_c2fd2a29_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Neg,
            exponent: -138,
            mantissa: 0xb2f1dd69_0a75bd7e_10483f42_11c5254b_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Neg,
            exponent: -133,
            mantissa: 0x823eb768_ace18aa0_57bbd86e_a32066b2_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Pos,
            exponent: -141,
            mantissa: 0xee7251c0_e7126d2e_b7d5ec4b_d23836d2_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Pos,
            exponent: -138,
            mantissa: 0xd02dc9c2_16d4d214_5ef4c306_97a68a66_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Neg,
            exponent: -145,
            mantissa: 0xbe71efb0_587507cf_6bea4375_786d6ff7_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Neg,
            exponent: -143,
            mantissa: 0x9e5a809c_6672567a_7d983c94_f6338653_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Pos,
            exponent: -150,
            mantissa: 0x90b74252_ef63a396_c1f7aa0c_a723ba19_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Pos,
            exponent: -149,
            mantissa: 0x8c71aa60_0a99d9d7_337b289b_6376092b_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Neg,
            exponent: -156,
            mantissa: 0x802a4981_cf2f82a7_8620563d_b6fd0318_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Neg,
            exponent: -156,
            mantissa: 0xa2f6a51b_99f4e5c2_1a67e6b7_8ca07b89_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Pos,
            exponent: -163,
            mantissa: 0x9471769c_478ba754_40e74bec_f3381648_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Pos,
            exponent: -163,
            mantissa: 0x8541fdc5_572e5dd9_3cd33f30_3eb767ec_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Neg,
            exponent: -171,
            mantissa: 0xf23a2d8f_20858db8_56a27690_82f408cd_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Neg,
            exponent: -171,
            mantissa: 0xa1cd5743_96d719c2_b2efa2e5_8672c1cb_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Pos,
            exponent: -178,
            mantissa: 0x92abd4ec_d28cf6af_f6c12c7e_8c233033_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Pos,
            exponent: -179,
            mantissa: 0x97989b23_e5a91c7b_b51bc702_5bc3fbb8_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Neg,
            exponent: -186,
            mantissa: 0x8901c4f9_cd36419d_d312bb7e_a5ae9c91_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Neg,
            exponent: -188,
            mantissa: 0xe1ce598a_f751d378_8754b9a4_12495d6e_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Pos,
            exponent: -195,
            mantissa: 0xcb62d3ba_45105ee7_a5833bf9_d6d5648c_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Pos,
            exponent: -196,
            mantissa: 0x88dee7f8_559f432d_d5a87d50_e3131d8d_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Neg,
            exponent: -204,
            mantissa: 0xf5a3df9e_869c9061_2bfe8fc1_0226953e_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Neg,
            exponent: -205,
            mantissa: 0x89a8c925_10af85b8_8f4c9746_d970a764_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Pos,
            exponent: -213,
            mantissa: 0xf60a2259_1b3b4ee1_173f2d2c_0ca541df_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Pos,
            exponent: -215,
            mantissa: 0xe969db43_cddabb64_6bf1c58f_e92d3684_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Neg,
            exponent: -222,
            mantissa: 0xce51fe13_51a9f4bb_675021c8_743a3cf8_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Neg,
            exponent: -224,
            mantissa: 0xa827dfc6_9a7fe3db_4971facf_8a0ed79e_u128,
        },
    ],
    [
        DyadicFloat128 {
            sign: DyadicSign::Pos,
            exponent: -131,
            mantissa: 0xc149ab61_7db3dcab_4e4fb606_81f88e4a_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Pos,
            exponent: -232,
            mantissa: 0x8765f1aa_60515e28_8bced7c5_5b64281c_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Neg,
            exponent: -132,
            mantissa: 0xc149ab61_7db3dcab_4e4fb606_8207b64e_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Pos,
            exponent: -140,
            mantissa: 0xe6cb272f_26550207_dc3866d5_18272d72_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Pos,
            exponent: -135,
            mantissa: 0x80c86733_5cbe19ad_ae94b604_802cc371_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Neg,
            exponent: -143,
            mantissa: 0xb86afec5_d6bda842_4eea1f3c_641e4bef_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Neg,
            exponent: -140,
            mantissa: 0x89350a82_34d828f1_2f0bace2_1340f7e7_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Pos,
            exponent: -148,
            mantissa: 0xd24f542c_0a8023f2_7ca9fc03_84a3e155_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Pos,
            exponent: -146,
            mantissa: 0x9c887ae0_acfa9f2c_33843de8_c5ed5639_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Neg,
            exponent: -154,
            mantissa: 0xf87a5edf_10647fae_a92546e6_e19f8691_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Neg,
            exponent: -153,
            mantissa: 0xde1b4bc1_4d9382e9_82f3750e_d5c6b8b2_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Pos,
            exponent: -160,
            mantissa: 0xb3f94ace_6a3e1306_a9431f30_ded0ffa4_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Pos,
            exponent: -160,
            mantissa: 0xd6c0e042_d5f8f0ee_37b22fae_f0bc00d2_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Neg,
            exponent: -167,
            mantissa: 0xb0508558_44c56085_28006915_4f0e6aad_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Neg,
            exponent: -167,
            mantissa: 0x9683fecc_c248fe29_a6e38b53_399931df_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Pos,
            exponent: -175,
            mantissa: 0xf9443d59_a5d25c1f_fbce6c6b_21d1e09d_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Pos,
            exponent: -175,
            mantissa: 0x9fe96237_a1dfc674_1e3b504c_104c595c_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Neg,
            exponent: -182,
            mantissa: 0x852410a5_96c14835_bfdbee1f_4ffbfba4_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Neg,
            exponent: -183,
            mantissa: 0x852e6927_3d704e60_9920f433_174d069b_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Pos,
            exponent: -191,
            mantissa: 0xde813a7f_23aec4c1_3d0ce932_4d09272e_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Pos,
            exponent: -192,
            mantissa: 0xb28c8cc9_a1f3777d_ca81c595_b40795e7_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Neg,
            exponent: -199,
            mantissa: 0x95664265_cc8ff240_018c6c2b_281ed815_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Neg,
            exponent: -201,
            mantissa: 0xc4ca97e4_53505137_8e992bff_236bd7ab_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Pos,
            exponent: -208,
            mantissa: 0xa4ba47f1_c2b94b67_351785d9_7f482df1_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Pos,
            exponent: -210,
            mantissa: 0xb572abf2_433a2497_79f0e2b2_e1688e63_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Neg,
            exponent: -217,
            mantissa: 0x97cdd31d_db1b6ce9_2aee86c3_cca8d23b_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Neg,
            exponent: -219,
            mantissa: 0x8d2818c3_bab62062_927c45ff_cb350f47_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Pos,
            exponent: -227,
            mantissa: 0xf047b1d0_14025af1_06e86a10_df1960d7_u128,
        },
    ],
    [
        DyadicFloat128 {
            sign: DyadicSign::Pos,
            exponent: -232,
            mantissa: 0xa72d4b69_e6323ebc_1001db4f_f23052c6_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Neg,
            exponent: -131,
            mantissa: 0xbf32b76a_bd0689ea_89f3a432_08089cd4_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Pos,
            exponent: -138,
            mantissa: 0xa7867e31_75b0e3a6_b0c67353_d604760f_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Pos,
            exponent: -134,
            mantissa: 0xfed5d31f_4e32a973_6dd17549_bf293598_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Neg,
            exponent: -141,
            mantissa: 0xdf3dd68e_282a34eb_a0fa7724_27416763_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Neg,
            exponent: -138,
            mantissa: 0xcbad603b_864c64c1_3f85366a_fd262f22_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Pos,
            exponent: -145,
            mantissa: 0xb25357c1_dd42fdd2_4642bbd0_06c99c71_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Pos,
            exponent: -143,
            mantissa: 0x9af35a0d_b6bbdab2_6f70e416_fde5b682_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Neg,
            exponent: -150,
            mantissa: 0x87896828_6f90c602_f26e5db9_b5efc695_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Neg,
            exponent: -149,
            mantissa: 0x8973a211_1762df91_ffa19918_3c01735a_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Pos,
            exponent: -157,
            mantissa: 0xf0254605_44df3060_af134d2f_d776c224_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Pos,
            exponent: -156,
            mantissa: 0x9f876af3_5069d97f_6d69f6f7_174b8267_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Neg,
            exponent: -163,
            mantissa: 0x8b1feb8a_bac9d6fd_de22fd80_096aba77_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Neg,
            exponent: -163,
            mantissa: 0x827c77f2_8cb29eef_1903fe84_eaf6da91_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Pos,
            exponent: -171,
            mantissa: 0xe320bfbd_08ea053c_8f3cedb7_b42a694b_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Pos,
            exponent: -171,
            mantissa: 0x9e7d157f_d7cf1c13_e8f10763_3c475c62_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Neg,
            exponent: -178,
            mantissa: 0x899a742b_d6e4907f_284b88b3_51714890_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Neg,
            exponent: -179,
            mantissa: 0x948bd70c_58cca90f_caa7e7a1_f846fece_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Pos,
            exponent: -186,
            mantissa: 0x809da1c4_ffddc107_5a83e768_5a6339ce_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Pos,
            exponent: -188,
            mantissa: 0xdd5a608a_0f4db93e_c31956d7_66e06ef4_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Neg,
            exponent: -195,
            mantissa: 0xbf0f540b_f0f89858_4b2df842_51d26211_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Neg,
            exponent: -196,
            mantissa: 0x863b0e34_9155bdba_d0d97114_0dd87b43_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Pos,
            exponent: -204,
            mantissa: 0xe6ecd4de_8865e73b_1fc1a81c_60258fc6_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Pos,
            exponent: -205,
            mantissa: 0x871157dc_bf0034e9_53429932_1c9c8cac_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Neg,
            exponent: -213,
            mantissa: 0xe77c6ed6_34b88b48_ec467fba_85eed694_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Neg,
            exponent: -215,
            mantissa: 0xe5222a30_c7c842dc_f2e70396_f32ec990_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Pos,
            exponent: -222,
            mantissa: 0xc247cc67_06e45555_bad79ad3_3b080235_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Pos,
            exponent: -224,
            mantissa: 0xa5261cef_79e86851_cac1224a_64c2834b_u128,
        },
    ],
    [
        DyadicFloat128 {
            sign: DyadicSign::Neg,
            exponent: -131,
            mantissa: 0xbd2cbf5b_b52f1a5a_5fe083fe_5ecd327e_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Pos,
            exponent: -232,
            mantissa: 0x975b8f38_044047a6_6370d0d8_8630818a_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Pos,
            exponent: -132,
            mantissa: 0xbd2cbf5b_b52f1a5a_5fe083fe_5ebcf793_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Neg,
            exponent: -140,
            mantissa: 0xd85ecf5b_0a98a662_3a9b0102_8a926257_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Neg,
            exponent: -136,
            mantissa: 0xfc18dcb1_4d669279_7d3d3758_35466f11_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Pos,
            exponent: -143,
            mantissa: 0xace8e206_4add3546_7c268b8f_4fe15f12_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Pos,
            exponent: -140,
            mantissa: 0x864e9aea_5f0c88eb_dcd68d1b_49e81a92_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Neg,
            exponent: -148,
            mantissa: 0xc5388df6_34905995_33f4678c_a4de9437_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Neg,
            exponent: -146,
            mantissa: 0x993f2179_c833708c_07952b05_a417d512_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Pos,
            exponent: -154,
            mantissa: 0xe912d546_f6b669b2_3dbc7dc1_491973dd_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Pos,
            exponent: -153,
            mantissa: 0xd97c1d4a_3576a65e_183a001a_095be5b2_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Neg,
            exponent: -160,
            mantissa: 0xa8df7837_5aecf6df_8f7c0bbc_c6904e87_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Neg,
            exponent: -160,
            mantissa: 0xd2558cc2_ce6739f4_50be91eb_6e8c384f_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Pos,
            exponent: -167,
            mantissa: 0xa581e719_cec1ce6f_819c4f87_7b729203_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Pos,
            exponent: -167,
            mantissa: 0x937595dd_531e3906_17d77666_24306324_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Neg,
            exponent: -175,
            mantissa: 0xea19eb50_d47e1cfb_feab1e25_d95114c8_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Neg,
            exponent: -175,
            mantissa: 0x9cb6eca1_45beb838_6c8e2223_e8f5a091_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Pos,
            exponent: -183,
            mantissa: 0xfa38429f_aecba39f_5bc8234b_0bce0c61_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Pos,
            exponent: -183,
            mantissa: 0x8290aeed_f37c5df7_53ff9271_cbc64808_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Neg,
            exponent: -191,
            mantissa: 0xd1368294_2a0e5bd3_bf9918f6_00901484_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Neg,
            exponent: -192,
            mantissa: 0xaf1c3aa7_90bed686_f807066b_4b876a9a_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Pos,
            exponent: -199,
            mantissa: 0x8c920001_4625a7fb_fc346bd4_e29fd561_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Pos,
            exponent: -201,
            mantissa: 0xc1153d3f_177d3085_acabd3e6_75a15b94_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Neg,
            exponent: -208,
            mantissa: 0x9b1b7a79_0f80fc7e_c8bafbb7_9a8fbac0_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Neg,
            exponent: -210,
            mantissa: 0xb21be630_317775a6_22516393_2c5c08e3_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Pos,
            exponent: -217,
            mantissa: 0x8f2b1b06_117cde99_731bc988_2a36ba3e_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Pos,
            exponent: -219,
            mantissa: 0x8a9ea82d_1ce04847_bd53bc9a_c06713e8_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Neg,
            exponent: -227,
            mantissa: 0xfd650d4b_ae2936bc_3ccfc8fe_9065a75e_u128,
        },
    ],
];
