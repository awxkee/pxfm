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
use crate::common::f_fmla;
use crate::dyadic_float::{DyadicFloat128, DyadicSign};
use crate::polyeval::f_polyeval8;

#[inline]
fn poly_exp_f128(x: DyadicFloat128) -> DyadicFloat128 {
    static COEFFS_128: [DyadicFloat128; 8] = [
        DyadicFloat128 {
            sign: DyadicSign::Pos,
            exponent: -127,
            mantissa: 0x80000000_00000000_00000000_00000000_u128,
        }, // 1.0
        DyadicFloat128 {
            sign: DyadicSign::Pos,
            exponent: -127,
            mantissa: 0x80000000_00000000_00000000_00000000_u128,
        }, // 1.0
        DyadicFloat128 {
            sign: DyadicSign::Pos,
            exponent: -128,
            mantissa: 0x80000000_00000000_00000000_00000000_u128,
        }, // 0.5
        DyadicFloat128 {
            sign: DyadicSign::Pos,
            exponent: -130,
            mantissa: 0xaaaaaaaa_aaaaaaaa_aaaaaaaa_aaaaaaab_u128,
        }, // 1/6
        DyadicFloat128 {
            sign: DyadicSign::Pos,
            exponent: -132,
            mantissa: 0xaaaaaaaa_aaaaaaaa_aaaaaaaa_aaaaaaab_u128,
        }, // 1/24
        DyadicFloat128 {
            sign: DyadicSign::Pos,
            exponent: -134,
            mantissa: 0x88888888_88888888_88888888_88888889_u128,
        }, // 1/120
        DyadicFloat128 {
            sign: DyadicSign::Pos,
            exponent: -137,
            mantissa: 0xb60b60b6_0b60b60b_60b60b60_b60b60b6_u128,
        }, // 1/720
        DyadicFloat128 {
            sign: DyadicSign::Pos,
            exponent: -140,
            mantissa: 0xd00d00d0_0d00d00d_00d00d00_d00d00d0_u128,
        }, // 1/5040
    ];

    f_polyeval8(
        x,
        COEFFS_128[0],
        COEFFS_128[1],
        COEFFS_128[2],
        COEFFS_128[3],
        COEFFS_128[4],
        COEFFS_128[5],
        COEFFS_128[6],
        COEFFS_128[7],
    )
}

#[cold]
#[inline(never)]
pub(crate) fn rational128_exp(x: f64) -> DyadicFloat128 {
    const LOG2_E: f64 = f64::from_bits(0x3ff71547652b82fe);
    let tmp = f_fmla(x, LOG2_E, f64::from_bits(0x4148000000040000));
    let k = (tmp.to_bits() >> 19) as i32;
    let kd = k as f64;

    let idx1: usize = ((k >> 6) & 0x3f) as usize;
    let idx2 = (k & 0x3f) as usize;

    const MLOG_2_EXP2_M12_HI: f64 = f64::from_bits(0xbf262e42ff000000);
    const MLOG_2_EXP2_M12_MID_30: f64 = f64::from_bits(0x3d0718432a000000);
    const MLOG_2_EXP2_M12_LO: f64 = f64::from_bits(0x3b0b0e2633fe0685);

    let t1 = f_fmla(kd, MLOG_2_EXP2_M12_HI, x); // exact
    let t2 = kd * MLOG_2_EXP2_M12_MID_30; // exact
    let t3 = kd * MLOG_2_EXP2_M12_LO; // Error < 2^-133

    let dx = DyadicFloat128::new_from_f64(t1)
        + (DyadicFloat128::new_from_f64(t2) + DyadicFloat128::new_from_f64(t3));

    let exp_mid1 = DyadicFloat128::new_from_f64(f64::from_bits(TD_EXP2_MID1[idx1].2))
        + (DyadicFloat128::new_from_f64(f64::from_bits(TD_EXP2_MID1[idx1].1))
            + DyadicFloat128::new_from_f64(f64::from_bits(TD_EXP2_MID1[idx1].0)));

    let exp_mid2 = DyadicFloat128::new_from_f64(f64::from_bits(TD_EXP2_MID2[idx2].2))
        + (DyadicFloat128::new_from_f64(f64::from_bits(TD_EXP2_MID2[idx2].1))
            + DyadicFloat128::new_from_f64(f64::from_bits(TD_EXP2_MID2[idx2].0)));

    let exp_mid = exp_mid1 * exp_mid2;

    let p = poly_exp_f128(dx);

    let mut r = exp_mid * p;

    r.exponent += (kd as i64 >> 12) as i16;

    r
}

// Lookup table for 2^(k * 2^-6) with k = 0..63.
// Generated by Sollya with:
// > display=hexadecimal;
// > prec = 500;
// > for i from 0 to 63 do {
//     a = 2^(i * 2^-6);
//     b = round(a, D, RN);
//     c = round(a - b, D, RN);
//     d = round(a - b - c, D, RN);
//     print("{", d, ",", c, ",", b, "},");
//   };
static TD_EXP2_MID1: [(u64, u64, u64); 64] = [
    (0x0000000000000000, 0x0000000000000000, 0x3ff0000000000000),
    (0xb919085b0a3d74d5, 0xbc719083535b085d, 0x3ff02c9a3e778061),
    (0x39105ff94f8d257e, 0x3c8d73e2a475b465, 0x3ff059b0d3158574),
    (0x39015820d96b414f, 0x3c6186be4bb284ff, 0x3ff0874518759bc8),
    (0xb9367c9bd6ebf74c, 0x3c98a62e4adc610b, 0x3ff0b5586cf9890f),
    (0xb8e5aa76994e9ddb, 0x3c403a1727c57b53, 0x3ff0e3ec32d3d1a2),
    (0x3929d58b988f562d, 0xbc96c51039449b3a, 0x3ff11301d0125b51),
    (0xb932fe7bb4c76416, 0xbc932fbf9af1369e, 0x3ff1429aaea92de0),
    (0x3924f2406aa13ff0, 0xbc819041b9d78a76, 0x3ff172b83c7d517b),
    (0x390ad36183926ae8, 0x3c8e5b4c7b4968e4, 0x3ff1a35beb6fcb75),
    (0x391ea62d0881b918, 0x3c9e016e00a2643c, 0x3ff1d4873168b9aa),
    (0xb90781dbc16f1ea4, 0x3c8dc775814a8495, 0x3ff2063b88628cd6),
    (0xb924d89f9af532e0, 0x3c99b07eb6c70573, 0x3ff2387a6e756238),
    (0x391277393a461b77, 0x3c82bd339940e9d9, 0x3ff26b4565e27cdd),
    (0x390de54485604690, 0x3c8612e8afad1255, 0x3ff29e9df51fdee1),
    (0xb91ee9d8f8cb9307, 0x3c90024754db41d5, 0x3ff2d285a6e4030b),
    (0x3917b7b2f09cd0d9, 0x3c86f46ad23182e4, 0x3ff306fe0a31b715),
    (0xb93406a2ea6cfc6b, 0x3c932721843659a6, 0x3ff33c08b26416ff),
    (0x39387e3e12516bfa, 0xbc963aeabf42eae2, 0x3ff371a7373aa9cb),
    (0x3909b0b1ff17c296, 0xbc75e436d661f5e3, 0x3ff3a7db34e59ff7),
    (0xb92808ba68fa8fb7, 0x3c8ada0911f09ebc, 0x3ff3dea64c123422),
    (0xb8d32b43eafc6518, 0xbc5ef3691c309278, 0x3ff4160a21f72e2a),
    (0xb8d0ac312de3d922, 0x3c489b7a04ef80d0, 0x3ff44e086061892d),
    (0x390e1eebae743ac0, 0x3c73c1a3b69062f0, 0x3ff486a2b5c13cd0),
    (0x38ec06c7745c2b39, 0x3c7d4397afec42e2, 0x3ff4bfdad5362a27),
    (0xb8f1aa1fd7b685cd, 0xbc94b309d25957e3, 0x3ff4f9b2769d2ca7),
    (0x390fa733951f214c, 0xbc807abe1db13cad, 0x3ff5342b569d4f82),
    (0xb90ff86852a613ff, 0x3c99bb2c011d93ad, 0x3ff56f4736b527da),
    (0xb92744ee506fdafe, 0x3c96324c054647ad, 0x3ff5ab07dd485429),
    (0xb9395f9ab75fa7d6, 0x3c9ba6f93080e65e, 0x3ff5e76f15ad2148),
    (0x3905d8e757cfb991, 0xbc9383c17e40b497, 0x3ff6247eb03a5585),
    (0x3934a337f4dc0a3b, 0xbc9bb60987591c34, 0x3ff6623882552225),
    (0x39357d3e3adec175, 0xbc9bdd3413b26456, 0x3ff6a09e667f3bcd),
    (0x38ca59f88abbe778, 0xbc6bbe3a683c88ab, 0x3ff6dfb23c651a2f),
    (0xb92269796953a4c3, 0xbc816e4786887a99, 0x3ff71f75e8ec5f74),
    (0xb938f8e7fa19e5e8, 0xbc90245957316dd3, 0x3ff75feb564267c9),
    (0xb8e4217a932d10d4, 0xbc841577ee04992f, 0x3ff7a11473eb0187),
    (0x38f70a1427f8fcdf, 0x3c705d02ba15797e, 0x3ff7e2f336cf4e62),
    (0x38f0f6ad65cbbac1, 0xbc9d4c1dd41532d8, 0x3ff82589994cce13),
    (0xb92f16f65181d921, 0xbc9fc6f89bd4f6ba, 0x3ff868d99b4492ed),
    (0xb9130644a7836333, 0x3c96e9f156864b27, 0x3ff8ace5422aa0db),
    (0x38d3bf26d2b85163, 0x3c85cc13a2e3976c, 0x3ff8f1ae99157736),
    (0x390697e257ac0db2, 0xbc675fc781b57ebc, 0x3ff93737b0cdc5e5),
    (0x3937edb9d7144b6f, 0xbc9d185b7c1b85d1, 0x3ff97d829fde4e50),
    (0x3916376b7943085c, 0x3c7c7c46b071f2be, 0x3ff9c49182a3f090),
    (0x392354084551b4fb, 0xbc9359495d1cd533, 0x3ffa0c667b5de565),
    (0xb90bfd7adfd63f48, 0xbc9d2f6edb8d41e1, 0x3ffa5503b23e255d),
    (0x3928b16ae39e8cb9, 0x3c90fac90ef7fd31, 0x3ffa9e6b5579fdbf),
    (0x393a7fbc3ae675ea, 0x3c97a1cd345dcc81, 0x3ffae89f995ad3ad),
    (0x3902babc0edda4d9, 0xbc62805e3084d708, 0x3ffb33a2b84f15fb),
    (0x390aa64481e1ab72, 0xbc75584f7e54ac3b, 0x3ffb7f76f2fb5e47),
    (0x3929a164050e1258, 0x3c823dd07a2d9e84, 0x3ffbcc1e904bc1d2),
    (0x39199e51125928da, 0x3c811065895048dd, 0x3ffc199bdd85529c),
    (0xb92fc44c329d5cb2, 0x3c92884dff483cad, 0x3ffc67f12e57d14b),
    (0x391d8765566b032e, 0x3c7503cbd1e949db, 0x3ffcb720dcef9069),
    (0xb93e7044039da0f6, 0xbc9cbc3743797a9c, 0x3ffd072d4a07897c),
    (0xb90ab053b05531fc, 0x3c82ed02d75b3707, 0x3ffd5818dcfba487),
    (0x3937f6246f0ec615, 0x3c9c2300696db532, 0x3ffda9e603db3285),
    (0x393b7225a944efd6, 0xbc91a5cd4f184b5c, 0x3ffdfc97337b9b5f),
    (0x3921e92cb3c2d278, 0x3c839e8980a9cc8f, 0x3ffe502ee78b3ff6),
    (0xb92fc0f242bbf3de, 0xbc9e9c23179c2893, 0x3ffea4afa2a490da),
    (0x393f6dd5d229ff69, 0x3c9dc7f486a4b6b0, 0x3ffefa1bee615a27),
    (0xb914019bffc80ef3, 0x3c99d3e12dd8a18b, 0x3fff50765b6e4540),
    (0x38fdc060c36f7651, 0x3c874853f3a5931e, 0x3fffa7c1819e90d8),
];

// Lookup table for 2^(k * 2^-12) with k = 0..63.
// Generated by Sollya with:
// > display=hexadecimal;
// > prec = 500;
// > for i from 0 to 63 do {
//     a = 2^(i * 2^-12);
//     b = round(a, D, RN);
//     c = round(a - b, D, RN);
//     d = round(a - b - c, D, RN);
//     print("{", d, ",", c, ",", b, "},");
//   };
static TD_EXP2_MID2: [(u64, u64, u64); 64] = [
    (0x0000000000000000, 0x0000000000000000, 0x3ff0000000000000),
    (0x39339726694630e3, 0x3c9ae8e38c59c72a, 0x3ff000b175effdc7),
    (0x38fe5e06ddd31156, 0xbc57b5d0d58ea8f4, 0x3ff00162f3904052),
    (0x3905a0768b51f609, 0x3c94115cb6b16a8e, 0x3ff0021478e11ce6),
    (0x390d008403605217, 0xbc8d7c96f201bb2f, 0x3ff002c605e2e8cf),
    (0x39289bc16f765708, 0x3c984711d4c35e9f, 0x3ff003779a95f959),
    (0xb924535b7f8c1e2d, 0xbc80484245243777, 0x3ff0042936faa3d8),
    (0xb938ba92f6b25456, 0xbc94b237da2025f9, 0x3ff004dadb113da0),
    (0xb8e30c72e81f4294, 0xbc75e00e62d6b30d, 0x3ff0058c86da1c0a),
    (0xb9134a5384e6f0b9, 0x3c9a1d6cedbb9481, 0x3ff0063e3a559473),
    (0x393f8d0580865d2e, 0xbc94acf197a00142, 0x3ff006eff583fc3d),
    (0xb90002bcb3ae9a99, 0xbc6eaf2ea42391a5, 0x3ff007a1b865a8ca),
    (0x390c3c5aedee9851, 0x3c7da93f90835f75, 0x3ff0085382faef83),
    (0x3927217851d1ec6e, 0xbc86a79084ab093c, 0x3ff00905554425d4),
    (0xb9180cbca335a7c3, 0x3c986364f8fbe8f8, 0x3ff009b72f41a12b),
    (0xb91706bd4eb22595, 0xbc882e8e14e3110e, 0x3ff00a6910f3b6fd),
    (0xb90b55dd523f3c08, 0xbc84f6b2a7609f71, 0x3ff00b1afa5abcbf),
    (0x39190a1e207cced1, 0xbc7e1a258ea8f71b, 0x3ff00bcceb7707ec),
    (0x39178d0472db37c5, 0x3c74362ca5bc26f1, 0x3ff00c7ee448ee02),
    (0xb92bcd4db3cb52fe, 0x3c9095a56c919d02, 0x3ff00d30e4d0c483),
    (0xb8fcf1b131575ec2, 0xbc6406ac4e81a645, 0x3ff00de2ed0ee0f5),
    (0xb8f6aaa1fa7ff913, 0x3c9b5a6902767e09, 0x3ff00e94fd0398e0),
    (0x39168f236dff3218, 0xbc991b2060859321, 0x3ff00f4714af41d3),
    (0xb92e8bb58067e60a, 0x3c8427068ab22306, 0x3ff00ff93412315c),
    (0x393d4cd5e1d71fdf, 0x3c9c1d0660524e08, 0x3ff010ab5b2cbd11),
    (0x393e4ecf350ebe88, 0xbc9e7bdfb3204be8, 0x3ff0115d89ff3a8b),
    (0x3926a2aa2c89c4f8, 0x3c8843aa8b9cbbc6, 0x3ff0120fc089ff63),
    (0x3911ca368a20ed05, 0xbc734104ee7edae9, 0x3ff012c1fecd613b),
    (0x38dedb1095d925cf, 0xbc72b6aeb6176892, 0x3ff0137444c9b5b5),
    (0xb90488c78eded75f, 0x3c7a8cd33b8a1bb3, 0x3ff01426927f5278),
    (0xb8e7480f5ea1b3c9, 0x3c72edc08e5da99a, 0x3ff014d8e7ee8d2f),
    (0xb90ae45989a04dd5, 0x3c857ba2dc7e0c73, 0x3ff0158b4517bb88),
    (0x392bf48007d80987, 0x3c9b61299ab8cdb7, 0x3ff0163da9fb3335),
    (0x3921aa91a059292c, 0xbc990565902c5f44, 0x3ff016f0169949ed),
    (0x391b6663292855f5, 0x3c870fc41c5c2d53, 0x3ff017a28af25567),
    (0x393e7fbca6793d94, 0x3c94b9a6e145d76c, 0x3ff018550706ab62),
    (0xb915b9f5c7de3b93, 0xbc7008eff5142bf9, 0x3ff019078ad6a19f),
    (0x3914638bf2f6acab, 0xbc977669f033c7de, 0x3ff019ba16628de2),
    (0xb92ab237b9a069c5, 0xbc909bb78eeead0a, 0x3ff01a6ca9aac5f3),
    (0x3933ab358be97cef, 0x3c9371231477ece5, 0x3ff01b1f44af9f9e),
    (0xb914027b2294bb64, 0x3c75e7626621eb5b, 0x3ff01bd1e77170b4),
    (0x390656394426c990, 0xbc9bc72b100828a5, 0x3ff01c8491f08f08),
    (0x390bf9785189bdd8, 0xbc6ce39cbbab8bbe, 0x3ff01d37442d5070),
    (0x3927c12f86114fe3, 0x3c816996709da2e2, 0x3ff01de9fe280ac8),
    (0xb92653d5d24b5d28, 0xbc8c11f5239bf535, 0x3ff01e9cbfe113ef),
    (0x39204a0cdc1d86d7, 0x3c8e1d4eb5edc6b3, 0x3ff01f4f8958c1c6),
    (0x392c678c46149782, 0xbc9afb99946ee3f0, 0x3ff020025a8f6a35),
    (0x39348524e1e9df70, 0xbc98f06d8a148a32, 0x3ff020b533856324),
    (0x3929953ea727ff0b, 0xbc82bf310fc54eb6, 0x3ff02168143b0281),
    (0xb93ccfbbec22d28e, 0xbc9c95a035eb4175, 0x3ff0221afcb09e3e),
    (0x3939e2bb6e181de1, 0xbc9491793e46834d, 0x3ff022cdece68c4f),
    (0x391f17609ae29308, 0xbc73e8d0d9c49091, 0x3ff02380e4dd22ad),
    (0xb91c7dc2c476bfb8, 0xbc9314aa16278aa3, 0x3ff02433e494b755),
    (0xb92fab994971d4a3, 0x3c848daf888e9651, 0x3ff024e6ec0da046),
    (0x392848b62cbdd0af, 0x3c856dc8046821f4, 0x3ff02599fb483385),
    (0xb92bf603ba715d0c, 0x3c945b42356b9d47, 0x3ff0264d1244c719),
    (0x39189434e751e1aa, 0xbc7082ef51b61d7e, 0x3ff027003103b10e),
    (0xb9103b54fd64e8ac, 0x3c72106ed0920a34, 0x3ff027b357854772),
    (0x3927785ea0acc486, 0xbc9fd4cf26ea5d0f, 0x3ff0286685c9e059),
    (0xb92ce447fdb35ff9, 0xbc909f8775e78084, 0x3ff02919bbd1d1d8),
    (0x38f5b884aab5642a, 0x3c564cbba902ca27, 0x3ff029ccf99d720a),
    (0xb93cfb3e46d7c1c0, 0x3c94383ef231d207, 0x3ff02a803f2d170d),
    (0xb8f0d40cee4b81af, 0x3c94a47a505b3a47, 0x3ff02b338c811703),
    (0x3926ae7d36d7c1f7, 0x3c9e47120223467f, 0x3ff02be6e199c811),
];

#[cfg(test)]
mod tests {
    use crate::exponents::exp_f128::rational128_exp;

    #[test]
    fn test_exp() {
        assert_eq!(rational128_exp(2.).fast_as_f64(), 7.38905609893065);
    }
}
