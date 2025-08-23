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
use crate::common::f_fmla;
use crate::double_double::DoubleDouble;
use crate::logs::fast_log_dd;
use crate::polyeval::{f_polyeval4, f_polyeval6, f_polyeval7, f_polyeval8, f_polyeval10};

#[cold]
fn inverf_0p06_to_0p75(x: f64) -> f64 {
    // First step rational approximant is generated, but it's ill-conditioned, thus
    // we're using taylor expansion to create Newton form at the point.
    //
    // <<FunctionApproximations`
    // ClearAll["Global`*"]
    // f[x_]:=InverseErf[x]/x
    // g[x_] =f[Sqrt[x]];
    // {err0,approx}=MiniMaxApproximation[g[z],{z,{0.06,0.75},9,9},WorkingPrecision->75, MaxIterations->100]
    // num=Numerator[approx][[1]];
    // den=Denominator[approx][[1]];
    // poly=den;
    // coeffs=CoefficientList[poly,z];
    // TableForm[Table[Row[{"'",NumberForm[coeffs[[i+1]],{50,50}, ExponentFunction->(Null&)],"',"}],{i,0,Length[coeffs]-1}]]
    // x0=SetPrecision[0.5625,75];
    // NumberForm[Series[num[x],{x,x0,50}], ExponentFunction->(Null&)]
    // coeffs=Table[SeriesCoefficient[num[x],{x,x0,k}],{k,0,9}];
    // TableForm[Table[Row[{"'",NumberForm[coeffs[[i+1]],{50,50}, ExponentFunction->(Null&)],"',"}],{i,0,Length[coeffs]-1}]];
    const P: [(u64, u64); 10] = [
        (0xbc3e06eda42202a0, 0x3f93c2fc5d00e0c8),
        (0xbc6eb374406b33b4, 0xbfc76fcfd022e3ff),
        (0xbc857822d7ffd282, 0x3fe6f8443546010a),
        (0x3c68269c66dfb28a, 0xbff80996754ceb79),
        (0x3c543dce8990a9f9, 0x3ffcf778d5ef0504),
        (0xbc72fc55f73765f6, 0xbff433be821423d0),
        (0xbc66d05fb37c8592, 0x3fdf15f19e9d8da4),
        (0x3c56dfb85e83a2c5, 0xbfb770b6827e0829),
        (0x3bff1472ecdfa403, 0x3f7a98a2980282bb),
        (0x3baffb33d69d6276, 0xbf142a246fd2c07c),
    ];
    let x2 = DoubleDouble::from_exact_mult(x, x);
    let z = DoubleDouble::full_add_f64(x2, -0.5625);
    let mut num = DoubleDouble::mul_add(
        DoubleDouble::from_bit_pair(P[9]),
        z,
        DoubleDouble::from_bit_pair(P[8]),
    );
    num = DoubleDouble::mul_add(z, num, DoubleDouble::from_bit_pair(P[7]));
    num = DoubleDouble::mul_add(z, num, DoubleDouble::from_bit_pair(P[6]));
    num = DoubleDouble::mul_add(z, num, DoubleDouble::from_bit_pair(P[5]));
    num = DoubleDouble::mul_add(z, num, DoubleDouble::from_bit_pair(P[4]));
    num = DoubleDouble::mul_add(z, num, DoubleDouble::from_bit_pair(P[3]));
    num = DoubleDouble::mul_add(z, num, DoubleDouble::from_bit_pair(P[2]));
    num = DoubleDouble::mul_add(z, num, DoubleDouble::from_bit_pair(P[1]));
    num = DoubleDouble::mul_add(z, num, DoubleDouble::from_bit_pair(P[0]));
    // <<FunctionApproximations`
    // ClearAll["Global`*"]
    // f[x_]:=InverseErf[x]/x
    // g[x_] =f[Sqrt[x]];
    // {err0,approx}=MiniMaxApproximation[g[z],{z,{0.06,0.75},9,9},WorkingPrecision->75, MaxIterations->100]
    // num=Numerator[approx][[1]];
    // den=Denominator[approx][[1]];
    // coeffs=CoefficientList[poly,z];
    // TableForm[Table[Row[{"'",NumberForm[coeffs[[i+1]],{50,50}, ExponentFunction->(Null&)],"',"}],{i,0,Length[coeffs]-1}]]
    // x0=SetPrecision[0.5625,75];
    // NumberForm[Series[den[x],{x,x0,50}], ExponentFunction->(Null&)]
    // coeffs=Table[SeriesCoefficient[den[x],{x,x0,k}],{k,0,9}];
    // TableForm[Table[Row[{"'",NumberForm[coeffs[[i+1]],{50,50}, ExponentFunction->(Null&)],"',"}],{i,0,Length[coeffs]-1}]];
    const Q: [(u64, u64); 10] = [
        (0xbc36337f24e57cb9, 0x3f92388d5d757e3a),
        (0xbc63dfae43d60e0b, 0xbfc6ca7da581358c),
        (0xbc77656389bd0e62, 0x3fe7c82ce417b4e0),
        (0xbc93679667bef2f0, 0xbffad58651fd1a51),
        (0x3ca2c6cb9eb17fb4, 0x4001bdb67e93a242),
        (0xbc9b58961ba253bc, 0xbffbdaeff6fbb81c),
        (0x3c7861f549c6aa61, 0x3fe91b12cf47da3a),
        (0xbc696dfd665b2f5e, 0xbfc7c5d0ffb7f1da),
        (0x3c1552b0ec0ba7b3, 0x3f939ada247f7609),
        (0xbbcaa226fb7b30a8, 0xbf41be65038ccfe6),
    ];

    let mut den = DoubleDouble::mul_add(
        DoubleDouble::from_bit_pair(Q[9]),
        z,
        DoubleDouble::from_bit_pair(Q[8]),
    );
    den = DoubleDouble::mul_add(z, den, DoubleDouble::from_bit_pair(Q[7]));
    den = DoubleDouble::mul_add(z, den, DoubleDouble::from_bit_pair(Q[6]));
    den = DoubleDouble::mul_add(z, den, DoubleDouble::from_bit_pair(Q[5]));
    den = DoubleDouble::mul_add(z, den, DoubleDouble::from_bit_pair(Q[4]));
    den = DoubleDouble::mul_add(z, den, DoubleDouble::from_bit_pair(Q[3]));
    den = DoubleDouble::mul_add(z, den, DoubleDouble::from_bit_pair(Q[2]));
    den = DoubleDouble::mul_add(z, den, DoubleDouble::from_bit_pair(Q[1]));
    den = DoubleDouble::mul_add(z, den, DoubleDouble::from_bit_pair(Q[0]));
    let r = DoubleDouble::div(num, den);
    let k = DoubleDouble::quick_mult_f64(r, x);
    f64::copysign(k.to_f64(), x)
}

/// Inverse error function
///
/// Max ulp 0.7, ulp 0.5 on [0;0.75]
pub fn f_erfinv(x: f64) -> f64 {
    let ax = x.to_bits() & 0x7fff_ffff_ffff_ffff;
    if ax >= 0x3ff0000000000000u64 {
        // |x| > 1
        if ax == 0x3ff0000000000000u64 {
            return if x.is_sign_negative() {
                f64::NEG_INFINITY
            } else {
                f64::INFINITY
            };
        }
        return f64::NAN;
    }
    if ax == 0 {
        return 0.;
    }

    let z = f64::from_bits(ax);

    if ax <= 0x3f8374bc6a7ef9db {
        // 0.0095
        // for small |x| using taylor series first 3 terms
        // Generated by SageMath:
        // from mpmath import mp, erf
        //
        // mp.prec = 100
        //
        // def inverf_series(n_terms):
        //     from mpmath import taylor
        //     series_erf = taylor(mp.erfinv, 0, n_terms)
        //     return series_erf
        //
        // ser = inverf_series(10)
        // for i in range(1, len(ser), 2):
        //     k = ser[i]
        //     print("f64::from_bits(" + double_to_hex(RealField(100)(k)) + "),")
        let z2 = DoubleDouble::from_exact_mult(z, z);
        let p = f_fmla(
            z2.hi,
            f64::from_bits(0x3fb62847c47dda48),
            f64::from_bits(0x3fc053c2c0ab91c5),
        );
        let mut r = DoubleDouble::mul_f64_add(
            z2,
            p,
            DoubleDouble::from_bit_pair((0xbc33ea2ef8dde075, 0x3fcdb29fb2fee5e4)),
        );
        r = DoubleDouble::mul_add(
            z2,
            r,
            DoubleDouble::from_bit_pair((0xbc8618f13eb7ca89, 0x3fec5bf891b4ef6b)),
        );
        // (rh + rl) * z = rh * z + rl*z
        let v = DoubleDouble::quick_mult_f64(r, z);
        return f64::copysign(v.to_f64(), x);
    } else if ax <= 0x3faeb851eb851eb8 {
        // 0.06
        // for |x| < 0.06 using taylor series first 5 terms
        // Generated by SageMath:
        // from mpmath import mp, erf
        //
        // mp.prec = 100
        //
        // def inverf_series(n_terms):
        //     from mpmath import taylor
        //     series_erf = taylor(mp.erfinv, 0, n_terms)
        //     return series_erf
        //
        // ser = inverf_series(10)
        // for i in range(1, len(ser), 2):
        //     k = ser[i]
        //     print("f64::from_bits(" + double_to_hex(RealField(100)(k)) + "),")
        let z2 = DoubleDouble::from_exact_mult(z, z);
        let p = f_polyeval4(
            z2.hi,
            f64::from_bits(0x3fb62847c47dda48),
            f64::from_bits(0x3fb0a13189c6ef7a),
            f64::from_bits(0x3faa7c85c89bb08b),
            f64::from_bits(0x3fa5eeb1d488e312),
        );
        let mut r = DoubleDouble::mul_f64_add(
            z2,
            p,
            DoubleDouble::from_bit_pair((0x3c2cec68daff0d80, 0x3fc053c2c0ab91c5)),
        );
        r = DoubleDouble::mul_add(
            z2,
            r,
            DoubleDouble::from_bit_pair((0xbc33ea2ef8dde075, 0x3fcdb29fb2fee5e4)),
        );
        r = DoubleDouble::mul_add(
            z2,
            r,
            DoubleDouble::from_bit_pair((0xbc8618f13eb7ca89, 0x3fec5bf891b4ef6b)),
        );
        // (rh + rl) * z = rh * z + rl*z
        let v = DoubleDouble::quick_mult_f64(r, z);
        return f64::copysign(v.to_f64(), x);
    }

    if ax <= 0x3fe8000000000000u64 {
        // |x| < 0.75

        // First step rational approximant is generated, but it's ill-conditioned, thus
        // we're using taylor expansion to create Newton form at the point.
        //
        // <<FunctionApproximations`
        // ClearAll["Global`*"]
        // f[x_]:=InverseErf[x]/x
        // g[x_] =f[Sqrt[x]];
        // {err0,approx}=MiniMaxApproximation[g[z],{z,{0.06,0.75},9,9},WorkingPrecision->75, MaxIterations->100]
        // num=Numerator[approx][[1]];
        // den=Denominator[approx][[1]];
        // poly=den;
        // coeffs=CoefficientList[poly,z];
        // TableForm[Table[Row[{"'",NumberForm[coeffs[[i+1]],{50,50}, ExponentFunction->(Null&)],"',"}],{i,0,Length[coeffs]-1}]]
        // x0=SetPrecision[0.5625,75];
        // NumberForm[Series[num[x],{x,x0,50}], ExponentFunction->(Null&)]
        // coeffs=Table[SeriesCoefficient[num[x],{x,x0,k}],{k,0,9}];
        // TableForm[Table[Row[{"'",NumberForm[coeffs[[i+1]],{50,50}, ExponentFunction->(Null&)],"',"}],{i,0,Length[coeffs]-1}]];
        const P: [(u64, u64); 3] = [
            (0xbc3e06eda42202a0, 0x3f93c2fc5d00e0c8),
            (0xbc6eb374406b33b4, 0xbfc76fcfd022e3ff),
            (0xbc857822d7ffd282, 0x3fe6f8443546010a),
        ];
        let x2 = DoubleDouble::from_exact_mult(x, x);
        let z = DoubleDouble::full_add_f64(x2, -0.5625);
        let ps_num = f_polyeval7(
            z.hi,
            f64::from_bits(0xbff80996754ceb79),
            f64::from_bits(0x3ffcf778d5ef0504),
            f64::from_bits(0xbff433be821423d0),
            f64::from_bits(0x3fdf15f19e9d8da4),
            f64::from_bits(0xbfb770b6827e0829),
            f64::from_bits(0x3f7a98a2980282bb),
            f64::from_bits(0xbf142a246fd2c07c),
        );
        let mut num = DoubleDouble::mul_f64_add(z, ps_num, DoubleDouble::from_bit_pair(P[2]));
        num = DoubleDouble::mul_add(z, num, DoubleDouble::from_bit_pair(P[1]));
        num = DoubleDouble::mul_add(z, num, DoubleDouble::from_bit_pair(P[0]));

        // <<FunctionApproximations`
        // ClearAll["Global`*"]
        // f[x_]:=InverseErf[x]/x
        // g[x_] =f[Sqrt[x]];
        // {err0,approx}=MiniMaxApproximation[g[z],{z,{0.06,0.75},9,9},WorkingPrecision->75, MaxIterations->100]
        // num=Numerator[approx][[1]];
        // den=Denominator[approx][[1]];
        // coeffs=CoefficientList[poly,z];
        // TableForm[Table[Row[{"'",NumberForm[coeffs[[i+1]],{50,50}, ExponentFunction->(Null&)],"',"}],{i,0,Length[coeffs]-1}]]
        // x0=SetPrecision[0.5625,75];
        // NumberForm[Series[den[x],{x,x0,50}], ExponentFunction->(Null&)]
        // coeffs=Table[SeriesCoefficient[den[x],{x,x0,k}],{k,0,9}];
        // TableForm[Table[Row[{"'",NumberForm[coeffs[[i+1]],{50,50}, ExponentFunction->(Null&)],"',"}],{i,0,Length[coeffs]-1}]];
        const Q: [(u64, u64); 3] = [
            (0xbc36337f24e57cb9, 0x3f92388d5d757e3a),
            (0xbc63dfae43d60e0b, 0xbfc6ca7da581358c),
            (0xbc77656389bd0e62, 0x3fe7c82ce417b4e0),
        ];

        let ps_den = f_polyeval7(
            z.hi,
            f64::from_bits(0xbffad58651fd1a51),
            f64::from_bits(0x4001bdb67e93a242),
            f64::from_bits(0xbffbdaeff6fbb81c),
            f64::from_bits(0x3fe91b12cf47da3a),
            f64::from_bits(0xbfc7c5d0ffb7f1da),
            f64::from_bits(0x3f939ada247f7609),
            f64::from_bits(0xbf41be65038ccfe6),
        );

        let mut den = DoubleDouble::mul_f64_add(z, ps_den, DoubleDouble::from_bit_pair(Q[2]));
        den = DoubleDouble::mul_add(z, den, DoubleDouble::from_bit_pair(Q[1]));
        den = DoubleDouble::mul_add(z, den, DoubleDouble::from_bit_pair(Q[0]));
        let r = DoubleDouble::div(num, den);
        let k = DoubleDouble::quick_mult_f64(r, x);
        let err = f_fmla(
            k.hi,
            f64::from_bits(0x3c10000000000000), // 2^-62
            f64::from_bits(0x3b50000000000000), // 2^-74
        );
        let ub = k.hi + (k.lo + err);
        let lb = k.hi + (k.lo - err);
        if ub == lb {
            return f64::copysign(k.to_f64(), x);
        }
        return inverf_0p06_to_0p75(x);
    }

    let q = DoubleDouble::from_full_exact_add(1.0, -z);

    let mut zeta = fast_log_dd(q);
    zeta = DoubleDouble::from_exact_add(zeta.hi, zeta.lo);
    zeta = -zeta;
    let zeta_sqrt = zeta.fast_sqrt();

    if zeta_sqrt.hi < 3.0 {
        const Y: f64 = 0.807220458984375;
        let xs = DoubleDouble::full_add_f64(zeta_sqrt, -1.125);

        const P: [f64; 11] = [
            -0.131102781679951906451,
            -0.163794047193317060787,
            0.117030156341995252019,
            0.387079738972604337464,
            0.337785538912035898924,
            0.142869534408157156766,
            0.0290157910005329060432,
            0.00214558995388805277169,
            -0.679465575181126350155e-6,
            0.285225331782217055858e-7,
            -0.681149956853776992068e-9,
        ];

        let p_num_s = f_polyeval10(
            xs.hi, P[1], P[2], P[3], P[4], P[5], P[6], P[7], P[8], P[9], P[10],
        );
        let p_num = DoubleDouble::mul_f64_add_f64(xs, p_num_s, P[0]);

        const Q: [f64; 8] = [
            1.0,
            3.46625407242567245975,
            5.38168345707006855425,
            4.77846592945843778382,
            2.59301921623620271374,
            0.848854343457902036425,
            0.152264338295331783612,
            0.01105924229346489121,
        ];

        let p_den_s = f_polyeval7(xs.hi, Q[1], Q[2], Q[3], Q[4], Q[5], Q[6], Q[7]);
        let p_den = DoubleDouble::mul_f64_add_f64(xs, p_den_s, Q[0]);
        let r = DoubleDouble::div(p_num, p_den);

        let r0 = DoubleDouble::quick_mult(r, zeta_sqrt);
        let p = DoubleDouble::mul_f64_add(zeta_sqrt, Y, r0);
        f64::copysign(p.to_f64(), x)
    } else if zeta_sqrt.hi < 6.0 {
        const Y: f64 = 0.93995571136474609375;
        let xs = DoubleDouble::full_add_f64(zeta_sqrt, -3.0);

        const P: [f64; 9] = [
            -0.0350353787183177984712,
            -0.00222426529213447927281,
            0.0185573306514231072324,
            0.00950804701325919603619,
            0.00187123492819559223345,
            0.000157544617424960554631,
            0.460469890584317994083e-5,
            -0.230404776911882601748e-9,
            0.266339227425782031962e-11,
        ];

        let p_num_s = f_polyeval8(xs.hi, P[1], P[2], P[3], P[4], P[5], P[6], P[7], P[8]);
        let p_num = DoubleDouble::mul_f64_add_f64(xs, p_num_s, P[0]);

        const Q: [f64; 7] = [
            1.0,
            1.3653349817554063097,
            0.762059164553623404043,
            0.220091105764131249824,
            0.0341589143670947727934,
            0.00263861676657015992959,
            0.764675292302794483503e-4,
        ];

        let p_den_s = f_polyeval6(xs.hi, Q[1], Q[2], Q[3], Q[4], Q[5], Q[6]);
        let p_den = DoubleDouble::mul_f64_add_f64(xs, p_den_s, Q[0]);
        let r = DoubleDouble::div(p_num, p_den);

        let r0 = DoubleDouble::quick_mult(r, zeta_sqrt);
        let p = DoubleDouble::mul_f64_add(zeta_sqrt, Y, r0);
        f64::copysign(p.to_f64(), x)
    } else if zeta_sqrt.hi < 18.0 {
        const Y: f64 = 0.98362827301025390625;
        let xs = DoubleDouble::full_add_f64(zeta_sqrt, -6.0);

        const P: [f64; 9] = [
            -0.0167431005076633737133,
            -0.00112951438745580278863,
            0.00105628862152492910091,
            0.000209386317487588078668,
            0.149624783758342370182e-4,
            0.449696789927706453732e-6,
            0.462596163522878599135e-8,
            -0.281128735628831791805e-13,
            0.99055709973310326855e-16,
        ];

        let p_num_s = f_polyeval8(xs.hi, P[1], P[2], P[3], P[4], P[5], P[6], P[7], P[8]);
        let p_num = DoubleDouble::mul_f64_add_f64(xs, p_num_s, P[0]);

        const Q: [f64; 7] = [
            1.0,
            0.591429344886417493481,
            0.138151865749083321638,
            0.0160746087093676504695,
            0.000964011807005165528527,
            0.275335474764726041141e-4,
            0.282243172016108031869e-6,
        ];

        let p_den_s = f_polyeval6(xs.hi, Q[1], Q[2], Q[3], Q[4], Q[5], Q[6]);
        let p_den = DoubleDouble::mul_f64_add_f64(xs, p_den_s, Q[0]);
        let r = DoubleDouble::div(p_num, p_den);

        let r0 = DoubleDouble::quick_mult(r, zeta_sqrt);
        let p = DoubleDouble::mul_f64_add(zeta_sqrt, Y, r0);
        f64::copysign(p.to_f64(), x)
    } else if zeta_sqrt.hi < 44.0 {
        const Y: f64 = 0.99714565277099609375;
        let xs = DoubleDouble::full_add_f64(zeta_sqrt, -18.0);

        const P: [f64; 8] = [
            -0.0024978212791898131227,
            -0.779190719229053954292e-5,
            0.254723037413027451751e-4,
            0.162397777342510920873e-5,
            0.396341011304801168516e-7,
            0.411632831190944208473e-9,
            0.145596286718675035587e-11,
            -0.116765012397184275695e-17,
        ];

        let p_num_s = f_polyeval7(xs.hi, P[1], P[2], P[3], P[4], P[5], P[6], P[7]);
        let p_num = DoubleDouble::mul_f64_add_f64(xs, p_num_s, P[0]);

        const Q: [f64; 7] = [
            1.0,
            0.207123112214422517181,
            0.0169410838120975906478,
            0.000690538265622684595676,
            0.145007359818232637924e-4,
            0.144437756628144157666e-6,
            0.509761276599778486139e-9,
        ];

        let p_den_s = f_polyeval6(xs.hi, Q[1], Q[2], Q[3], Q[4], Q[5], Q[6]);
        let p_den = DoubleDouble::mul_f64_add_f64(xs, p_den_s, Q[0]);
        let r = DoubleDouble::div(p_num, p_den);

        let r0 = DoubleDouble::quick_mult(r, zeta_sqrt);
        let p = DoubleDouble::mul_f64_add(zeta_sqrt, Y, r0);
        f64::copysign(p.to_f64(), x)
    } else {
        const Y: f64 = 0.99941349029541015625;
        let xs = DoubleDouble::full_add_f64(zeta_sqrt, -44.0);

        const P: [f64; 8] = [
            -0.000539042911019078575891,
            -0.28398759004727721098e-6,
            0.899465114892291446442e-6,
            0.229345859265920864296e-7,
            0.225561444863500149219e-9,
            0.947846627503022684216e-12,
            0.135880130108924861008e-14,
            -0.348890393399948882918e-21,
        ];

        let p_num_s = f_polyeval7(xs.hi, P[1], P[2], P[3], P[4], P[5], P[6], P[7]);
        let p_num = DoubleDouble::mul_f64_add_f64(xs, p_num_s, P[0]);

        const Q: [f64; 7] = [
            1.0,
            0.0845746234001899436914,
            0.00282092984726264681981,
            0.468292921940894236786e-4,
            0.399968812193862100054e-6,
            0.161809290887904476097e-8,
            0.231558608310259605225e-11,
        ];

        let p_den_s = f_polyeval6(xs.hi, Q[1], Q[2], Q[3], Q[4], Q[5], Q[6]);
        let p_den = DoubleDouble::mul_f64_add_f64(xs, p_den_s, Q[0]);
        let r = DoubleDouble::div(p_num, p_den);

        let r0 = DoubleDouble::quick_mult(r, zeta_sqrt);
        let p = DoubleDouble::mul_f64_add(zeta_sqrt, Y, r0);
        f64::copysign(p.to_f64(), x)
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_erfinv() {
        assert_eq!(f_erfinv(0.001000000000084706), 0.0008862271575416209);
        assert_eq!(f_erfinv(-0.001000000000084706), -0.0008862271575416209);
        assert_eq!(f_erfinv(0.71), 0.7482049711849852);
        assert_eq!(f_erfinv(-0.71), -0.7482049711849852);
        assert_eq!(f_erfinv(0.41), 0.381014610957532);
        assert_eq!(f_erfinv(-0.41), -0.381014610957532);
        assert_eq!(f_erfinv(0.32), 0.29165547581744206);
        assert_eq!(f_erfinv(-0.32), -0.29165547581744206);
        assert_eq!(f_erfinv(0.82), 0.9480569762323499);
        assert_eq!(f_erfinv(-0.82), -0.9480569762323499);
        assert_eq!(f_erfinv(0.05), 0.044340387910005497);
        assert_eq!(f_erfinv(-0.05), -0.044340387910005497);
        assert_eq!(f_erfinv(0.99), 1.8213863677184494);
        assert_eq!(f_erfinv(-0.99), -1.8213863677184494);
        assert_eq!(f_erfinv(0.9900000000867389), 1.821386369839293);
        assert_eq!(f_erfinv(0.99999), 3.123413274341571);
        assert!(f_erfinv(f64::NEG_INFINITY).is_nan());
        assert!(f_erfinv(f64::INFINITY).is_nan());
        assert!(f_erfinv(f64::NAN).is_nan());
    }
}
