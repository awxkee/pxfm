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
use crate::dekker::Dekker;
use crate::dyadic_float::{DyadicFloat128, DyadicSign};
use crate::f_sin;
use crate::sincos_dyadic::r_fmla;

/*
   Sage math:
   def format_hex(value):
   l = hex(value)[2:]
   n = 8
   x = [l[i:i + n] for i in range(0, len(l), n)]
   return "0x" + "_".join(x) + "_u128"

   def print_dyadic(value):
       (s, m, e) = RealField(128)(value).sign_mantissa_exponent();
       print("DyadicFloat128 {")
       print(f"    sign: DyadicSign::{'Pos' if s >= 0 else 'Neg'},")
       print(f"    exponent: {e},")
       print(f"    mantissa: {format_hex(m)},")
       print("},")

   arr = [1,
   -0.16666666666666666666666663759455503219383933132492,
   8.333333333333333333321622304375846892928056249074e-3,
   -1.9841269841269841109350457617809654358034836364015e-4,
   2.7557319223984845104409221628075830059151085168279e-6,
   -2.5052108381793520761125795078702331987183520443341e-8,
   1.6059036834818088311380777393878900204395001831901e-10,
   -7.6401944360594586348474147669857868784676013162697e-13]

   for num in arr:
       print_dyadic(num)
*/
static DYADIC_0P0_TO_0_25: [DyadicFloat128; 8] = [
    DyadicFloat128 {
        sign: DyadicSign::Pos,
        exponent: -127,
        mantissa: 0x80000000_00000000_00000000_00000000_u128,
    },
    DyadicFloat128 {
        sign: DyadicSign::Neg,
        exponent: -130,
        mantissa: 0xaaaaaaaa_aaaaa800_00000000_00000000_u128,
    },
    DyadicFloat128 {
        sign: DyadicSign::Pos,
        exponent: -134,
        mantissa: 0x88888888_88888800_00000000_00000000_u128,
    },
    DyadicFloat128 {
        sign: DyadicSign::Neg,
        exponent: -140,
        mantissa: 0xd00d00d0_0d00d000_00000000_00000000_u128,
    },
    DyadicFloat128 {
        sign: DyadicSign::Pos,
        exponent: -146,
        mantissa: 0xb8ef1d2a_b631e800_00000000_00000000_u128,
    },
    DyadicFloat128 {
        sign: DyadicSign::Neg,
        exponent: -153,
        mantissa: 0xd7322b3f_238ed800_00000000_00000000_u128,
    },
    DyadicFloat128 {
        sign: DyadicSign::Pos,
        exponent: -160,
        mantissa: 0xb0922b91_9f23c800_00000000_00000000_u128,
    },
    DyadicFloat128 {
        sign: DyadicSign::Neg,
        exponent: -168,
        mantissa: 0xd70d6724_bed4c000_00000000_00000000_u128,
    },
];

/*
   Sage math:
   def format_hex(value):
   l = hex(value)[2:]
   n = 8
   x = [l[i:i + n] for i in range(0, len(l), n)]
   return "0x" + "_".join(x) + "_u128"

   def print_dyadic(value):
       (s, m, e) = RealField(128)(value).sign_mantissa_exponent();
       print("DyadicFloat128 {")
       print(f"    sign: DyadicSign::{'Pos' if s >= 0 else 'Neg'},")
       print(f"    exponent: {e},")
       print(f"    mantissa: {format_hex(m)},")
       print("},")

   arr = [1,
   -0.16666666666666666666652467408447354578599623840058,
   8.333333333333333322066227569718928911050320613513e-3,
   -1.9841269841269803319127574978448925981835256668121e-4,
   2.7557319223915544151944427792496719633161442895723e-6,
   -2.5052108307927863349060866846635968331974050926464e-8,
   1.6058993051954926021431109116837627701217886346218e-10,
   -7.6288391097349136337505644035050873914638104845285e-13]

   for num in arr:
       print_dyadic(num)
*/
static DYADIC_0_25_TO_0_35: [DyadicFloat128; 8] = [
    DyadicFloat128 {
        sign: DyadicSign::Pos,
        exponent: -127,
        mantissa: 0x80000000_00000000_00000000_00000000_u128,
    },
    DyadicFloat128 {
        sign: DyadicSign::Neg,
        exponent: -130,
        mantissa: 0xaaaaaaaa_aaaaa800_00000000_00000000_u128,
    },
    DyadicFloat128 {
        sign: DyadicSign::Pos,
        exponent: -134,
        mantissa: 0x88888888_88888800_00000000_00000000_u128,
    },
    DyadicFloat128 {
        sign: DyadicSign::Neg,
        exponent: -140,
        mantissa: 0xd00d00d0_0d006000_00000000_00000000_u128,
    },
    DyadicFloat128 {
        sign: DyadicSign::Pos,
        exponent: -146,
        mantissa: 0xb8ef1d2a_b4328800_00000000_00000000_u128,
    },
    DyadicFloat128 {
        sign: DyadicSign::Neg,
        exponent: -153,
        mantissa: 0xd7322b34_7e657000_00000000_00000000_u128,
    },
    DyadicFloat128 {
        sign: DyadicSign::Pos,
        exponent: -160,
        mantissa: 0xb0920c05_1be27000_00000000_00000000_u128,
    },
    DyadicFloat128 {
        sign: DyadicSign::Neg,
        exponent: -168,
        mantissa: 0xd6bb9443_80712000_00000000_00000000_u128,
    },
];

/*
   Sage math:
   def format_hex(value):
   l = hex(value)[2:]
   n = 8
   x = [l[i:i + n] for i in range(0, len(l), n)]
   return "0x" + "_".join(x) + "_u128"

   def print_dyadic(value):
       (s, m, e) = RealField(128)(value).sign_mantissa_exponent();
       print("DyadicFloat128 {")
       print(f"    sign: DyadicSign::{'Pos' if s >= 0 else 'Neg'},")
       print(f"    exponent: {e},")
       print(f"    mantissa: {format_hex(m)},")
       print("},")

   arr = [1,
   -0.16666666666666666666666665431284625027515646277734,
   8.3333333333333333333327623938421076994176503878756e-3,
   -1.98412698412698412686838285684373080553452431753257e-4,
   2.755731922398588930143680288443694184320533479124e-6,
   -2.5052108385440717687220969400875175283842309493642e-8,
   1.6059043836333236616166829243490112374308232224007e-10,
   -7.6471635748909637567300641342863708348011679270333e-13,
   2.8114252130721001205097328473556375427344650805421e-15,
   -8.1828914745930529226989291705301970732437343257723e-18]

   for num in arr:
       print_dyadic(num)
*/
static DYADIC_0_35_TO_0_5: [DyadicFloat128; 10] = [
    DyadicFloat128 {
        sign: DyadicSign::Pos,
        exponent: -127,
        mantissa: 0x80000000_00000000_00000000_00000000_u128,
    },
    DyadicFloat128 {
        sign: DyadicSign::Neg,
        exponent: -130,
        mantissa: 0xaaaaaaaa_aaaaa800_00000000_00000000_u128,
    },
    DyadicFloat128 {
        sign: DyadicSign::Pos,
        exponent: -134,
        mantissa: 0x88888888_88888800_00000000_00000000_u128,
    },
    DyadicFloat128 {
        sign: DyadicSign::Neg,
        exponent: -140,
        mantissa: 0xd00d00d0_0d00d000_00000000_00000000_u128,
    },
    DyadicFloat128 {
        sign: DyadicSign::Pos,
        exponent: -146,
        mantissa: 0xb8ef1d2a_b6399800_00000000_00000000_u128,
    },
    DyadicFloat128 {
        sign: DyadicSign::Neg,
        exponent: -153,
        mantissa: 0xd7322b3f_aa1da800_00000000_00000000_u128,
    },
    DyadicFloat128 {
        sign: DyadicSign::Pos,
        exponent: -160,
        mantissa: 0xb092309d_2c582800_00000000_00000000_u128,
    },
    DyadicFloat128 {
        sign: DyadicSign::Neg,
        exponent: -168,
        mantissa: 0xd73f9eef_82383800_00000000_00000000_u128,
    },
    DyadicFloat128 {
        sign: DyadicSign::Pos,
        exponent: -176,
        mantissa: 0xca95a431_ff989800_00000000_00000000_u128,
    },
    DyadicFloat128 {
        sign: DyadicSign::Neg,
        exponent: -184,
        mantissa: 0x96f29cc8_61a58000_00000000_00000000_u128,
    },
];

/*
   Sage math:
   def format_hex(value):
   l = hex(value)[2:]
   n = 8
   x = [l[i:i + n] for i in range(0, len(l), n)]
   return "0x" + "_".join(x) + "_u128"

   def print_dyadic(value):
       (s, m, e) = RealField(128)(value).sign_mantissa_exponent();
       print("DyadicFloat128 {")
       print(f"    sign: DyadicSign::{'Pos' if s >= 0 else 'Neg'},")
       print(f"    exponent: {e},")
       print(f"    mantissa: {format_hex(m)},")
       print("},")

   arr = [1,
   -0.16666666666666666666666457755319003657011473920948,
   8.3333333333333333332799500865042489779193054365186e-3,
   -1.9841269841269841209647234639869759609278513863403e-4,
   2.7557319223985851340171180878330048974263404613262e-6,
   -2.50521083854253292969531023672633482348464287071005e-8,
   1.6059043832297815024612439555735437671335294789624e-10,
   -7.6471629049979143065955605064518111859983827437937e-13,
   2.8113607377069428950749441559509686261993662551303e-15,
   -8.15530791168233478835577453776451246995095753727e-18]

   for num in arr:
       print_dyadic(num)
*/
static DYADIC_0_5_TO_0_7: [DyadicFloat128; 10] = [
    DyadicFloat128 {
        sign: DyadicSign::Pos,
        exponent: -127,
        mantissa: 0x80000000_00000000_00000000_00000000_u128,
    },
    DyadicFloat128 {
        sign: DyadicSign::Neg,
        exponent: -130,
        mantissa: 0xaaaaaaaa_aaaaa800_00000000_00000000_u128,
    },
    DyadicFloat128 {
        sign: DyadicSign::Pos,
        exponent: -134,
        mantissa: 0x88888888_88888800_00000000_00000000_u128,
    },
    DyadicFloat128 {
        sign: DyadicSign::Neg,
        exponent: -140,
        mantissa: 0xd00d00d0_0d00d000_00000000_00000000_u128,
    },
    DyadicFloat128 {
        sign: DyadicSign::Pos,
        exponent: -146,
        mantissa: 0xb8ef1d2a_b6395000_00000000_00000000_u128,
    },
    DyadicFloat128 {
        sign: DyadicSign::Neg,
        exponent: -153,
        mantissa: 0xd7322b3f_a98c5000_00000000_00000000_u128,
    },
    DyadicFloat128 {
        sign: DyadicSign::Pos,
        exponent: -160,
        mantissa: 0xb092309c_6dc6e800_00000000_00000000_u128,
    },
    DyadicFloat128 {
        sign: DyadicSign::Neg,
        exponent: -168,
        mantissa: 0xd73f9db3_291e8800_00000000_00000000_u128,
    },
    DyadicFloat128 {
        sign: DyadicSign::Pos,
        exponent: -176,
        mantissa: 0xca9473b8_1093c000_00000000_00000000_u128,
    },
    DyadicFloat128 {
        sign: DyadicSign::Neg,
        exponent: -184,
        mantissa: 0x96705a4d_24a24000_00000000_00000000_u128,
    },
];

/*
   Sage math:
   def format_hex(value):
   l = hex(value)[2:]
   n = 8
   x = [l[i:i + n] for i in range(0, len(l), n)]
   return "0x" + "_".join(x) + "_u128"

   def print_dyadic(value):
       (s, m, e) = RealField(128)(value).sign_mantissa_exponent();
       print("DyadicFloat128 {")
       print(f"    sign: DyadicSign::{'Pos' if s >= 0 else 'Neg'},")
       print(f"    exponent: {e},")
       print(f"    mantissa: {format_hex(m)},")
       print("},")

   arr = [1,
   -0.16666666666666666666646596082488799120297575619405,
   8.3333333333333333303034230814156232744345975825785e-3,
   -1.98412698412698392416773188737223235113727427914005e-4,
   2.755731922398510054277513740669281275414831803543e-6,
   -2.505210838524429747541859402731849124539145861329e-8,
   1.6059043804009079988086180387031993498811476720109e-10,
   -7.6471601038206065567009808851818755693151221391883e-13,
   2.8111998741297096040595618349158665759292753445849e-15,
   -8.114261525007097674139405605791126614246779754623e-18]

   for num in arr:
       print_dyadic(num)
*/
static DYADIC_0_7_TO_0_85: [DyadicFloat128; 10] = [
    DyadicFloat128 {
        sign: DyadicSign::Pos,
        exponent: -127,
        mantissa: 0x80000000_00000000_00000000_00000000_u128,
    },
    DyadicFloat128 {
        sign: DyadicSign::Neg,
        exponent: -130,
        mantissa: 0xaaaaaaaa_aaaaa800_00000000_00000000_u128,
    },
    DyadicFloat128 {
        sign: DyadicSign::Pos,
        exponent: -134,
        mantissa: 0x88888888_88888800_00000000_00000000_u128,
    },
    DyadicFloat128 {
        sign: DyadicSign::Neg,
        exponent: -140,
        mantissa: 0xd00d00d0_0d00c800_00000000_00000000_u128,
    },
    DyadicFloat128 {
        sign: DyadicSign::Pos,
        exponent: -146,
        mantissa: 0xb8ef1d2a_b633c800_00000000_00000000_u128,
    },
    DyadicFloat128 {
        sign: DyadicSign::Neg,
        exponent: -153,
        mantissa: 0xd7322b3f_a2de8800_00000000_00000000_u128,
    },
    DyadicFloat128 {
        sign: DyadicSign::Pos,
        exponent: -160,
        mantissa: 0xb0923097_35e11800_00000000_00000000_u128,
    },
    DyadicFloat128 {
        sign: DyadicSign::Neg,
        exponent: -168,
        mantissa: 0xd73f9888_578fe000_00000000_00000000_u128,
    },
    DyadicFloat128 {
        sign: DyadicSign::Pos,
        exponent: -176,
        mantissa: 0xca917c0f_eecc1800_00000000_00000000_u128,
    },
    DyadicFloat128 {
        sign: DyadicSign::Neg,
        exponent: -184,
        mantissa: 0x95ae8443_c2813800_00000000_00000000_u128,
    },
];

/*
   Sage math:
    def format_hex(value):
        l = hex(value)[2:]
        n = 8
        x = [l[i:i + n] for i in range(0, len(l), n)]
        return "0x" + "_".join(x) + "_u128"

    def print_dyadic(value):
        (s, m, e) = RealField(128)(value).sign_mantissa_exponent();
        print("DyadicFloat128 {")
        print(f"    sign: DyadicSign::{'Pos' if s >= 0 else 'Neg'},")
        print(f"    exponent: {e},")
        print(f"    mantissa: {format_hex(m)},")
        print("},")

    arr = [1,
    -0.166666666666666666666666655046098333674241424421013,
    8.3333333333333333333331831421004348543883293530357e-3,
    -1.9841269841269841269753150083731137729917315898448e-4,
    2.7557319223985890621576620783905775088236557741698e-6,
    -2.5052108385441711522991020660308014545508855347942e-8,
    1.60590438368204277871370189786500351551966543276647e-10,
    -7.6471637316812555218399795763893819204414937261132e-13,
    2.8114572428040154437944346095476833593632491361069e-15,
    -8.2206285242834984003520751985570987913105938287195e-18,
    1.9570332899853812841548996540992839174401821857085e-20,
    -3.807458690891922170447691589270102843057880400501e-23]

    for num in arr:
        print_dyadic(num)
*/
static DYADIC_0_85_TO_1_0: [DyadicFloat128; 12] = [
    DyadicFloat128 {
        sign: DyadicSign::Pos,
        exponent: -127,
        mantissa: 0x80000000_00000000_00000000_00000000_u128,
    },
    DyadicFloat128 {
        sign: DyadicSign::Neg,
        exponent: -130,
        mantissa: 0xaaaaaaaa_aaaaa800_00000000_00000000_u128,
    },
    DyadicFloat128 {
        sign: DyadicSign::Pos,
        exponent: -134,
        mantissa: 0x88888888_88888800_00000000_00000000_u128,
    },
    DyadicFloat128 {
        sign: DyadicSign::Neg,
        exponent: -140,
        mantissa: 0xd00d00d0_0d00d000_00000000_00000000_u128,
    },
    DyadicFloat128 {
        sign: DyadicSign::Pos,
        exponent: -146,
        mantissa: 0xb8ef1d2a_b639a000_00000000_00000000_u128,
    },
    DyadicFloat128 {
        sign: DyadicSign::Neg,
        exponent: -153,
        mantissa: 0xd7322b3f_aa270800_00000000_00000000_u128,
    },
    DyadicFloat128 {
        sign: DyadicSign::Pos,
        exponent: -160,
        mantissa: 0xb092309d_4359f000_00000000_00000000_u128,
    },
    DyadicFloat128 {
        sign: DyadicSign::Neg,
        exponent: -168,
        mantissa: 0xd73f9f39_8d00b800_00000000_00000000_u128,
    },
    DyadicFloat128 {
        sign: DyadicSign::Pos,
        exponent: -176,
        mantissa: 0xca963b73_917e1000_00000000_00000000_u128,
    },
    DyadicFloat128 {
        sign: DyadicSign::Neg,
        exponent: -184,
        mantissa: 0x97a4d213_93548000_00000000_00000000_u128,
    },
    DyadicFloat128 {
        sign: DyadicSign::Pos,
        exponent: -193,
        mantissa: 0xb8d62956_e54f7800_00000000_00000000_u128,
    },
    DyadicFloat128 {
        sign: DyadicSign::Neg,
        exponent: -202,
        mantissa: 0xb81e0e3c_b4c52800_00000000_00000000_u128,
    },
];

/*
   Sage math:
   def format_hex(value):
   l = hex(value)[2:]
   n = 8
   x = [l[i:i + n] for i in range(0, len(l), n)]
   return "0x" + "_".join(x) + "_u128"

   def print_dyadic(value):
       (s, m, e) = RealField(128)(value).sign_mantissa_exponent();
       print("DyadicFloat128 {")
       print(f"    sign: DyadicSign::{'Pos' if s >= 0 else 'Neg'},")
       print(f"    exponent: {e},")
       print(f"    mantissa: {format_hex(m)},")
       print("},")

   arr = [1,
   -0.16666666666666666666666643664074995503870539433371,
   8.3333333333333333333310754191591155962798591128976e-3,
   -1.98412698412698412688344569559121358875731681453906e-4,
   2.75573192239858903833548257322620497320002006087195e-6,
   -2.5052108385441670816902962774374525805905780365301e-8,
   1.6059043836815637392008153483480068323823016116603e-10,
   -7.6471637312879804931918019928763953970133515753401e-13,
   2.8114572205599720535422372191081920816278443975254e-15,
   -8.2206202274994648164166381371359674475851646872964e-18,
   1.95684899378927996422263752178740040059682095522686e-20,
   -3.7889439405371494879185903453568935723791975549856e-23]

   for num in arr:
       print_dyadic(num)
*/
static DYADIC_1_0_TO_1_2: [DyadicFloat128; 12] = [
    DyadicFloat128 {
        sign: DyadicSign::Pos,
        exponent: -127,
        mantissa: 0x80000000_00000000_00000000_00000000_u128,
    },
    DyadicFloat128 {
        sign: DyadicSign::Neg,
        exponent: -130,
        mantissa: 0xaaaaaaaa_aaaaa800_00000000_00000000_u128,
    },
    DyadicFloat128 {
        sign: DyadicSign::Pos,
        exponent: -134,
        mantissa: 0x88888888_88888800_00000000_00000000_u128,
    },
    DyadicFloat128 {
        sign: DyadicSign::Neg,
        exponent: -140,
        mantissa: 0xd00d00d0_0d00d000_00000000_00000000_u128,
    },
    DyadicFloat128 {
        sign: DyadicSign::Pos,
        exponent: -146,
        mantissa: 0xb8ef1d2a_b6399800_00000000_00000000_u128,
    },
    DyadicFloat128 {
        sign: DyadicSign::Neg,
        exponent: -153,
        mantissa: 0xd7322b3f_aa26a800_00000000_00000000_u128,
    },
    DyadicFloat128 {
        sign: DyadicSign::Pos,
        exponent: -160,
        mantissa: 0xb092309d_43200800_00000000_00000000_u128,
    },
    DyadicFloat128 {
        sign: DyadicSign::Neg,
        exponent: -168,
        mantissa: 0xd73f9f39_5d757000_00000000_00000000_u128,
    },
    DyadicFloat128 {
        sign: DyadicSign::Pos,
        exponent: -176,
        mantissa: 0xca963b58_ad4b6000_00000000_00000000_u128,
    },
    DyadicFloat128 {
        sign: DyadicSign::Neg,
        exponent: -184,
        mantissa: 0x97a4c80b_d8605000_00000000_00000000_u128,
    },
    DyadicFloat128 {
        sign: DyadicSign::Pos,
        exponent: -193,
        mantissa: 0xb8d1b499_ea8e5800_00000000_00000000_u128,
    },
    DyadicFloat128 {
        sign: DyadicSign::Neg,
        exponent: -202,
        mantissa: 0xb738daa6_b09fe000_00000000_00000000_u128,
    },
];

#[inline(always)]
#[allow(clippy::too_many_arguments)]
pub(crate) fn r_polyeval8(
    x: &DyadicFloat128,
    a0: &DyadicFloat128,
    a1: &DyadicFloat128,
    a2: &DyadicFloat128,
    a3: &DyadicFloat128,
    a4: &DyadicFloat128,
    a5: &DyadicFloat128,
    a6: &DyadicFloat128,
    a7: &DyadicFloat128,
) -> DyadicFloat128 {
    let z0 = r_fmla(x, a7, a6); // a3 * x + a2
    let t1 = r_fmla(x, &z0, a5); // a3 * x + a2
    let t2 = r_fmla(x, &t1, a4); // a3 * x + a2
    let t3 = r_fmla(x, &t2, a3); // (a3 * x + a2) * x + a1
    let t4 = r_fmla(x, &t3, a2); // (a3 * x + a2) * x + a1
    let t5 = r_fmla(x, &t4, a1); // (a3 * x + a2) * x + a1
    r_fmla(x, &t5, a0) // ((a3 * x + a2) * x + a1) * x + a0
}

#[inline(always)]
#[allow(clippy::too_many_arguments)]
pub(crate) fn r_polyeval10(
    x: &DyadicFloat128,
    a0: &DyadicFloat128,
    a1: &DyadicFloat128,
    a2: &DyadicFloat128,
    a3: &DyadicFloat128,
    a4: &DyadicFloat128,
    a5: &DyadicFloat128,
    a6: &DyadicFloat128,
    a7: &DyadicFloat128,
    a8: &DyadicFloat128,
    a9: &DyadicFloat128,
) -> DyadicFloat128 {
    let t0 = r_fmla(x, a9, a8); // a3 * x + a2
    let z0 = r_fmla(x, &t0, a7); // a3 * x + a2
    let t0a = r_fmla(x, &z0, a6); // a3 * x + a2
    let t1 = r_fmla(x, &t0a, a5); // a3 * x + a2
    let t2 = r_fmla(x, &t1, a4); // a3 * x + a2
    let t3 = r_fmla(x, &t2, a3); // (a3 * x + a2) * x + a1
    let t4 = r_fmla(x, &t3, a2); // (a3 * x + a2) * x + a1
    let t5 = r_fmla(x, &t4, a1); // (a3 * x + a2) * x + a1
    r_fmla(x, &t5, a0) // ((a3 * x + a2) * x + a1) * x + a0
}

#[inline(always)]
#[allow(clippy::too_many_arguments)]
pub(crate) fn r_polyeval12(
    x: &DyadicFloat128,
    a0: &DyadicFloat128,
    a1: &DyadicFloat128,
    a2: &DyadicFloat128,
    a3: &DyadicFloat128,
    a4: &DyadicFloat128,
    a5: &DyadicFloat128,
    a6: &DyadicFloat128,
    a7: &DyadicFloat128,
    a8: &DyadicFloat128,
    a9: &DyadicFloat128,
    a10: &DyadicFloat128,
    a11: &DyadicFloat128,
) -> DyadicFloat128 {
    let t0 = r_fmla(x, a11, a10);
    let k0 = r_fmla(x, &t0, a9);
    let k1 = r_fmla(x, &k0, a8);
    let z0 = r_fmla(x, &k1, a7);
    let t0a = r_fmla(x, &z0, a6);
    let t1 = r_fmla(x, &t0a, a5);
    let t2 = r_fmla(x, &t1, a4);
    let t3 = r_fmla(x, &t2, a3);
    let t4 = r_fmla(x, &t3, a2);
    let t5 = r_fmla(x, &t4, a1);
    r_fmla(x, &t5, a0)
}

fn eval_sinc8_dyadic(x: f64, dd: &[DyadicFloat128; 8]) -> f64 {
    let x = DyadicFloat128::new_from_f64(x);
    let x2 = x.quick_mul(&x);
    let p = r_polyeval8(
        &x2, &dd[0], &dd[1], &dd[2], &dd[3], &dd[4], &dd[5], &dd[6], &dd[7],
    );

    p.fast_as_f64()
}

fn eval_sinc10_dyadic(x: f64, dd: &[DyadicFloat128; 10]) -> f64 {
    let x = DyadicFloat128::new_from_f64(x);
    let x2 = x.quick_mul(&x);
    let p = r_polyeval10(
        &x2, &dd[0], &dd[1], &dd[2], &dd[3], &dd[4], &dd[5], &dd[6], &dd[7], &dd[8], &dd[9],
    );

    p.fast_as_f64()
}

fn eval_sinc12_dyadic(x: f64, dd: &[DyadicFloat128; 12]) -> f64 {
    let x = DyadicFloat128::new_from_f64(x);
    let x2 = x.quick_mul(&x);
    let p = r_polyeval12(
        &x2, &dd[0], &dd[1], &dd[2], &dd[3], &dd[4], &dd[5], &dd[6], &dd[7], &dd[8], &dd[9],
        &dd[10], &dd[11],
    );

    p.fast_as_f64()
}

#[inline(always)]
#[allow(clippy::too_many_arguments)]
fn f_polyeval5_skip0(x: f64, a1: f64, a2: f64, a3: f64, a4: f64, a5: f64) -> f64 {
    let t0 = f_fmla(x, a5, a4); // a3 * x + a2
    let t2 = f_fmla(x, t0, a3); // a3 * x + a2
    let t3 = f_fmla(x, t2, a2); // (a3 * x + a2) * x + a1
    let t4 = f_fmla(x, t3, a1); // (a3 * x + a2) * x + a1
    x * t4
}

#[inline(always)]
#[allow(clippy::too_many_arguments)]
fn f_polyeval7_skip0(x: f64, a1: f64, a2: f64, a3: f64, a4: f64, a5: f64, a6: f64, a7: f64) -> f64 {
    let z0 = f_fmla(x, a7, a6);
    let z1 = f_fmla(x, z0, a5);
    let t0 = f_fmla(x, z1, a4);
    let t2 = f_fmla(x, t0, a3);
    let t3 = f_fmla(x, t2, a2);
    let t4 = f_fmla(x, t3, a1);
    x * t4
}

#[inline(always)]
#[allow(clippy::too_many_arguments)]
fn f_polyeval9_skip0(
    x: f64,
    a1: f64,
    a2: f64,
    a3: f64,
    a4: f64,
    a5: f64,
    a6: f64,
    a7: f64,
    a8: f64,
    a9: f64,
) -> f64 {
    let k0 = f_fmla(x, a9, a8);
    let k1 = f_fmla(x, k0, a7);
    let k2 = f_fmla(x, k1, a6);
    let z1 = f_fmla(x, k2, a5);
    let t0 = f_fmla(x, z1, a4);
    let t2 = f_fmla(x, t0, a3);
    let t3 = f_fmla(x, t2, a2);
    let t4 = f_fmla(x, t3, a1);
    x * t4
}

fn eval_sinc(x: f64, c: &[u64; 6], c_dyadic: &[DyadicFloat128; 8], err: f64) -> f64 {
    let z2 = Dekker::from_exact_mult(x, x);
    let x2 = z2.to_f64();
    let p = f_polyeval5_skip0(
        x2,
        f64::from_bits(c[1]),
        f64::from_bits(c[2]),
        f64::from_bits(c[3]),
        f64::from_bits(c[4]),
        f64::from_bits(c[5]),
    );
    let p0 = Dekker::from_exact_add(f64::from_bits(c[0]), p);

    let err = err;
    let lb = p0.hi + (p0.lo - err);
    let ub = p0.hi + (p0.lo + err);
    // Ziv's accuracy test
    if lb != ub {
        return eval_sinc8_dyadic(x, c_dyadic);
    }
    lb
}

fn eval_sinc2(x: f64, c: &[u64; 8], c_dyadic: &[DyadicFloat128; 10], err: f64) -> f64 {
    let z2 = Dekker::from_exact_mult(x, x);
    let x2 = z2.to_f64();
    let p = f_polyeval7_skip0(
        x2,
        f64::from_bits(c[1]),
        f64::from_bits(c[2]),
        f64::from_bits(c[3]),
        f64::from_bits(c[4]),
        f64::from_bits(c[5]),
        f64::from_bits(c[6]),
        f64::from_bits(c[7]),
    );
    let p0 = Dekker::from_exact_add(f64::from_bits(c[0]), p);

    let err = err;
    let lb = p0.hi + (p0.lo - err);
    let ub = p0.hi + (p0.lo + err);
    // Ziv's accuracy test
    if lb != ub {
        return eval_sinc10_dyadic(x, c_dyadic);
    }
    lb
}

fn eval_sinc3(x: f64, c: &[u64; 10], c_dyadic: &[DyadicFloat128; 12], err: f64) -> f64 {
    let z2 = Dekker::from_exact_mult(x, x);
    let x2 = z2.to_f64();
    let p = f_polyeval9_skip0(
        x2,
        f64::from_bits(c[1]),
        f64::from_bits(c[2]),
        f64::from_bits(c[3]),
        f64::from_bits(c[4]),
        f64::from_bits(c[5]),
        f64::from_bits(c[6]),
        f64::from_bits(c[7]),
        f64::from_bits(c[8]),
        f64::from_bits(c[9]),
    );
    let p0 = Dekker::from_exact_add(f64::from_bits(c[0]), p);

    let err = err;
    let lb = p0.hi + (p0.lo - err);
    let ub = p0.hi + (p0.lo + err);
    // Ziv's accuracy test
    if lb != ub {
        return eval_sinc12_dyadic(x, c_dyadic);
    }
    lb
}

/// Computes sinc(x)
///
/// Max ULP 0.6 on [-1;1], after 1.0
pub fn f_sinc(x: f64) -> f64 {
    if !x.is_finite() {
        return f64::NAN;
    }
    let x_abs = f64::from_bits(x.to_bits() & 0x7fff_ffff_ffff_ffff);
    if x_abs.to_bits() == 0 {
        return 1.0;
    }
    if x_abs < 0.25 {
        static C: [u64; 6] = [
            0x3ff0000000000000,
            0xbfc5555555555555,
            0x3f81111111110d8f,
            0xbf2a01a019cd434e,
            0x3ec71de23e7a591f,
            0xbe5add41c334be4f,
        ];
        return eval_sinc(
            x,
            &C,
            &DYADIC_0P0_TO_0_25,
            f64::from_bits(0x3bd4bfe5b10e6f11),
        );
    } else if x_abs < 0.35 {
        static C: [u64; 6] = [
            0x3ff0000000000000,
            0xbfc5555555555531,
            0x3f811111111090f9,
            0xbf2a01a017396460,
            0x3ec71ddc04a28b5e,
            0xbe5ad1a21a7a1409,
        ];
        return eval_sinc(
            x,
            &C,
            &DYADIC_0_25_TO_0_35,
            f64::from_bits(0x3bdb20923105c2ea),
        );
    } else if x_abs < 0.5 {
        static C: [u64; 8] = [
            0x3ff0000000000000,
            0xbfc5555555555555,
            0x3f8111111111105d,
            0xbf2a01a019ff2ee7,
            0x3ec71de3a0df5f5b,
            0xbe5ae63c52fed023,
            0x3de6088734c4c504,
            0xbd623ecc53ecab44,
        ];
        return eval_sinc2(
            x,
            &C,
            &DYADIC_0_35_TO_0_5,
            f64::from_bits(0x3b60945ac3361eb1),
        );
    } else if x_abs < 0.7 {
        static C: [u64; 8] = [
            0x3ff0000000000000,
            0xbfc5555555555555,
            0x3f811111111110b3,
            0xbf2a01a01a00f2b8,
            0x3ec71de3a4aa37f1,
            0xbe5ae644a0b6772a,
            0x3de611c1c0f10ee5,
            0xbd6a8bb78de6424c,
        ];
        return eval_sinc2(
            x,
            &C,
            &DYADIC_0_5_TO_0_7,
            f64::from_bits(0x3b6270b0494c849a),
        );
    } else if x_abs < 0.85 {
        static C: [u64; 8] = [
            0x3ff0000000000000,
            0xbfc5555555555552,
            0x3f81111111110eac,
            0xbf2a01a019feaf85,
            0x3ec71de3a3572930,
            0xbe5ae643c7da0c34,
            0x3de6117b92ba06a0,
            0xbd6a7aa99c6ef724,
        ];
        return eval_sinc2(
            x,
            &C,
            &DYADIC_0_7_TO_0_85,
            f64::from_bits(0x3b3e13f14f87cd1b),
        );
    } else if x_abs < 1.0 {
        static C: [u64; 10] = [
            0x3ff0000000000000,
            0xbfc5555555555555,
            0x3f811111111110df,
            0xbf2a01a01a016ccb,
            0x3ec71de3a538c9d0,
            0xbe5ae64552167295,
            0x3de6123be34bb9f9,
            0xbd6ae207905e7e2a,
            0x3ce75bcec8aae2f0,
            0x3c7b0201cc269b71,
        ];
        return eval_sinc3(
            x,
            &C,
            &DYADIC_0_85_TO_1_0,
            f64::from_bits(0x3aa462950ecf078a),
        );
    } else if x_abs < 1.12 {
        static C: [u64; 10] = [
            0x3ff0000000000000,
            0xbfc5555555555555,
            0x3f811111111110eb,
            0xbf2a01a01a018268,
            0x3ec71de3a5498ca2,
            0xbe5ae645609843d2,
            0x3de6124374313d8e,
            0xbd6ae6c8aaaa907d,
            0x3ce9064ed1fbd96c,
            0xbc54863d7851eb24,
        ];
        return eval_sinc3(
            x,
            &C,
            &DYADIC_1_0_TO_1_2,
            f64::from_bits(0x3a6bb06f9c34b85a),
        );
    }

    // Refine with Newton-raphson
    let one = Dekker::new(0.0, 1.0);
    let mut recip = Dekker::from_exact_div(1.0, x);

    // First Newton-Raphson refinement
    let prod1 = Dekker::mult_f64(recip, x); // x * y₀
    let err1 = Dekker::sub(one, prod1); // 1 - x*y₀
    recip = Dekker::add(recip, Dekker::mult(recip, err1)); // y₁ = y₀ + y₀*(1 - x*y₀)

    // Second Newton-Raphson refinement
    let prod2 = Dekker::mult_f64(recip, x);
    let err2 = Dekker::sub(one, prod2);
    recip = Dekker::add(recip, Dekker::mult(recip, err2));

    // Now compute sin(x)/x ≈ sin(x) * (1/x refined)
    let sinx = f_sin(x);
    let result = Dekker::f64_mult(sinx, recip);
    result.to_f64()
}

#[cfg(test)]
mod tests {
    use super::*;
    #[test]
    fn test_sinc() {
        assert_eq!(f_sinc(0.1), 0.9983341664682815);
        assert_eq!(f_sinc(0.9), 0.870363232919426);
        assert_eq!(f_sinc(-0.1), 0.9983341664682815);
        assert_eq!(f_sinc(-0.9), 0.870363232919426);
    }
}
