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

#[derive(Clone)]
pub(crate) struct J1TaylorExtendedSeries {
    pub(crate) a0: (u64, u64),
    pub(crate) a1: (u64, u64),
    pub(crate) a2: (u64, u64),
    pub(crate) a3: (u64, u64),
    pub(crate) a4: (u64, u64),
    pub(crate) a5: (u64, u64),
    pub(crate) a6: (u64, u64),
    pub(crate) a7: (u64, u64),
    pub(crate) a8: (u64, u64),
    pub(crate) a9: (u64, u64),
    pub(crate) a10: (u64, u64),
    pub(crate) a11: (u64, u64),
    pub(crate) a12: (u64, u64),
    pub(crate) a13: (u64, u64),
    pub(crate) c: [u64; 8],
}

/**
J1 zeros and extremums on [-76;76]

Generated in SageMath:
```python
from mpmath import mp, mpf, findroot, j1
from sage.all import *
import struct

DR = RealField(52)
DD = RealField(190)

def double_to_hex(f):
    packed = struct.pack('>d', float(f))
    return '0x' + packed.hex()

def split_double_double(x):
    x_hi = DR(x)  # convert to f64
    x_lo = x - DD(x_hi)
    return (x_lo,x_hi)

def print_double_double(mark, x):
    splat = split_double_double(x)
    print(f"{mark}({double_to_hex(splat[0])}, {double_to_hex(splat[1])}),")

zeros = []

# Step size to detect sign changes
mp.prec = 150

step = mpf("0.001")
epsilon = mpf("1e-35")
x = mpf("0.0")

previous_zero = R120(0)
j1_zeros = []

while x < mpf("76.0"):
    f1 = besselj(1, x)
    f2 = besselj(1, x + step)
    if f1 * f2 < 0:
        zero = findroot(lambda t: j1(t), (x, x + step), solver='bisect', tol=mp.mpf("1e-41"))
        previous_zero = zero
        j1_zeros.append(zero)
    if previous_zero is not None and abs(x - mpf(f'{round(x)}')) < epsilon:
        zeros.append(previous_zero)
    x += step

j1_extrema = []

x = mpf("0.0")
while x < mpf("76.0"):
    d1 = mp.diff(lambda t: j1(t), x)
    d2 = mp.diff(lambda t: j1(t), x + step)
    if d1 * d2 < 0:
        extremum = findroot(lambda t: mp.diff(lambda u: j1(u), t), (x, x + step), solver='bisect', tol=mp.mpf("1e-41"))
        j1_extrema.append(extremum)
    x += step

j1_zeros.extend(j1_extrema)
j1_zeros = sorted(j1_zeros)

print("static J1_ZEROS: [(u64, u64); 46] = [")
for z in j1_zeros:
    k = split_double_double(DD(z))
    hi = double_to_hex(k[1])
    lo = double_to_hex(k[0])
    print(f"({lo}, {hi}),")

print("];")
```
See notes/bessel_j1_taylor.ipynb
**/
pub(crate) static J1_ZEROS: [(u64, u64); 48] = [
    (0x0, 0x0),
    (0x3c5616d820cfdaeb, 0x3ffd757d1fec8a3a),
    (0xbca60155a9d1b257, 0x400ea75575af6f09),
    (0x3ca5c646a75d7539, 0x40155365bc032467),
    (0xbc9b226d9d243828, 0x401c0ff5f3b47250),
    (0xbca63e17ec20a31d, 0x402112980f0b88a1),
    (0x3cc02610a51562b6, 0x402458d0d0bdfc29),
    (0x3cc9a84d3a5fedc1, 0x40276979797ee5ac),
    (0x3cb2bce7fd18e693, 0x402aa5baf310e5a2),
    (0xbcc6932b987094f1, 0x402dba284a17ac59),
    (0xbcdd2a68e88ab318, 0x4030787b360508c5),
    (0xbca022f6b2b54db9, 0x403203f9a24e6527),
    (0xbcd21830197e9e86, 0x40339da8e7416ca4),
    (0x3cdeaafeaf8ec1af, 0x40352a1424a1a9fa),
    (0xbcc1bf33afef88f2, 0x4036c294e3d4d8ac),
    (0xbcb2d773b50cf8b9, 0x40384fb31dee1635),
    (0x3cc1a2686480d882, 0x4039e7570dcea106),
    (0x3cd0bdee27293d79, 0x403b75014427514d),
    (0xbcb42ce39ec976fb, 0x403d0bfcf471fccc),
    (0x3cbda49c2c143483, 0x403e9a179fba532a),
    (0xbcdbe3a1cd066b67, 0x404018476e6b2bf0),
    (0xbce6b00c1279ef0a, 0x4040df82ebd54e32),
    (0xbced5fbbff045068, 0x4041aa890dc5e97c),
    (0x3cd7d864bbf17a30, 0x404271eb1b80430e),
    (0x3cc9eafeca0ca4fc, 0x40433cc523d5cb69),
    (0xbce5cecac300a9a1, 0x40440447e50db184),
    (0x3cc489bd556e5109, 0x4044cefcf1734b62),
    (0x3cdd0fd96f29c211, 0x4045969bc7271083),
    (0x3ce4f716f3179d90, 0x404661315d6b133f),
    (0xbce158b763edd0e8, 0x404728e892a88fc9),
    (0xbcef3950a842db79, 0x4047f36312028ad6),
    (0x3ce97656bbc2396e, 0x4048bb2fa2037de3),
    (0x3ce85d7bdb30baf1, 0x404985928f96d51e),
    (0xbce71f8560ac9f18, 0x404a4d71fcb56f8c),
    (0x3ce3d41e041caa68, 0x404b17c038c2018c),
    (0xbcde6d04716d8d21, 0x404bdfb06eb790aa),
    (0x3cda139ce2cd08ac, 0x404ca9ec5a82324b),
    (0x3cc8b5cc7b4501c1, 0x404d71eb98682f07),
    (0xbcb12e6ef2e594e2, 0x404e3c1731d64f1e),
    (0x3cb399bfca430021, 0x404f0423f99b4b53),
    (0x3cdfd1ee8286358a, 0x404fce40efb1156e),
    (0x3c800660b51502f0, 0x40504b2cfcbb084d),
    (0x3ced3cacfc720418, 0x4050b034dde75b42),
    (0x3cc4b877d4f6d900, 0x40511446f60f1458),
    (0xbcee669304bfe748, 0x40517948db63675c),
    (0x3cfad20ca758a714, 0x4051dd600b743a9b),
    (0x3cf8eb4a94b63936, 0x4052425c7dcacdf6),
    (0xbcfa196892f68386, 0x4052a67859bc641e),
];

/**
Precomputed values in exact Bessel J1 zero.

Generated by MPFR:
```text
let mut arr = vec![];
for zeros in J1_ZEROS.iter() {
    let mpfr = Float::with_val(107, f64::from_bits(zeros.1)).j1();
    arr.push(mpfr.to_f64().to_bits());
}
println!(
    "arr: [{}]",
    arr.iter()
        .map(|x| format!("0x{:016X}", x))
        .collect::<Vec<_>>()
        .join(", ")
);
```
**/
pub(crate) static J1_ZEROS_VALUE: [u64; 48] = [
    0x0000000000000000,
    0x3FE29EA3D19F035F,
    0xBC91B9C1C3FB286F,
    0xBFD626EE83500BF2,
    0x3C8049770CE74C2E,
    0x3FD17DBF09D40D25,
    0x3CA0212F4E592523,
    0xBFCDDCEB4CE1BF4A,
    0xBC905DCC62D0D222,
    0x3FCA7F63FEA81F26,
    0xBCB6EB905BA2ABFA,
    0xBFC810F50225B04B,
    0x3CAA10B2F7B4E69D,
    0x3FC633E7F7F05301,
    0xBC97BC6D5A660382,
    0xBFC4B71D4CA2CC69,
    0xBC961C29FAC28FDF,
    0x3FC37DFA8F5A550A,
    0xBC87E3B01386785F,
    0xBFC2768D29C69936,
    0x3CAF5EFD41F756B6,
    0x3FC194EBA75B32F9,
    0xBCBF89DCEDB3EA9B,
    0xBFC0D0D36473E98C,
    0xBC9AAAF726A29E97,
    0x3FC02455675AB6D2,
    0x3C9451B6225ACBFB,
    0xBFBF161D0C28B48C,
    0xBCB40032091A4E00,
    0x3FBE0357C158B119,
    0xBCBCCB5A05A6E4AA,
    0xBFBD0B36E5737458,
    0xBCB5C457E4A6A2F1,
    0x3FBC29AE8400A320,
    0x3CB13169F65EFC7C,
    0xBFBB5B8273B75055,
    0xBCA5FB7DBD93E256,
    0x3FBA9E13A0DB6429,
    0xBC7C3482175F80D7,
    0xBFB9EF3BB2213B0B,
    0xBCA977092852774B,
    0x3FB94D3276914E51,
    0x3CB6D73591BFEB5D,
    0xBFB8B67A2481077D,
    0x3CB735BC851F7831,
    0x3FB829D06FEE9266,
    0x3CC29C7C75EEB12F,
    0xBFB7A62320798175,
];

/**
Following search for J1 (see [J1_ZEROS]) zeros and extremums:
at each zero and extremum we're doing Taylor series expansion
one that should be enough to cover whole interval between zero or peak
which is PI/4

Generated in SageMath and Sollya:
```python
def compute_intervals(zeros):
    intervals = []
    for i in range(0, len(zeros)):
        if i == 0:
            a = (zeros[i]) / 2 - 0.05 - zeros[i]
            b = (zeros[i] + zeros[i + 1]) / 2 + 0.05 - zeros[i]
            intervals.append((RealField(18)(a), RealField(18)(b), RealField(110)(zeros[i])))
        elif i + 1 > len(zeros) - 1:
            a = (zeros[i - 1] + zeros[i]) / 2 - 0.05 - zeros[i]
            b = (zeros[i]) + 0.83 + 0.05 - zeros[i]
            intervals.append((RealField(18)(a), RealField(18)(b), RealField(110)(zeros[i])))
        else:
            a = (zeros[i - 1] + zeros[i]) / 2 - zeros[i] - 0.05
            b = (zeros[i] + zeros[i + 1]) / 2 + 0.05  - zeros[i]
            intervals.append((RealField(18)(a), RealField(18)(b), RealField(110)(zeros[i])))
    return intervals

intervals = compute_intervals(j1_zeros)
# print(intervals)

def build_sollya_script(a, b, zero, deg):
    return f"""
prec = 500;
bessel_j1 = library("./notes/bessel_sollya/cmake-build-release/libbessel_sollya.dylib");
f = bessel_j1(x + {zero});
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
    print("J1TaylorExtendedSeries {")
    print_double_double("a0: ", poly[0])
    print_double_double("a1: ", poly[1])
    print_double_double("a2: ", poly[2])
    print_double_double("a3: ", poly[3])
    print_double_double("a4: ", poly[4])
    print_double_double("a5: ", poly[5])
    print_double_double("a6: ", poly[6])
    print_double_double("a7: ", poly[7])
    print_double_double("a8: ", poly[8])
    print_double_double("a9: ", poly[9])
    print_double_double("a10: ", poly[10])
    print_double_double("a11: ", poly[11])
    print_double_double("a12: ", poly[12])
    print_double_double("a13: ", poly[13])
    print("c: [")
    for i in range(14, len(poly)):
        coeff = poly[i]
        print(f"{double_to_hex(coeff)},")
    print("],")
    print("},")


degree = 21

print(f"pub(crate) static J1_COEFFS: [J1TaylorExtendedSeries; {len(intervals)}] = [")
for i in range(0, len(intervals)):
    interval = intervals[i]
    call_sollya_on_interval(interval[0], interval[1], interval[2], degree)
    coeffs = load_coefficients(f"coefficients.txt")
    print_remez_coeffs(coeffs)
print("];")
```
**/
pub(crate) static J1_COEFFS: [J1TaylorExtendedSeries; 47] = [
    J1TaylorExtendedSeries {
        a0: (0x3c61f1c32444b011, 0x3fe29ea3d19f035f),
        a1: (0xb6f8500beb501283, 0xba5893b073dc3c7d),
        a2: (0xbc6e3631ad953349, 0xbfca41115c5df243),
        a3: (0xbc18acc496ba163f, 0x3f78d1448e6fed48),
        a4: (0x3c0e85a2c4c2a951, 0x3f8c441a2f9de22b),
        a5: (0x3bd94bf2390ffca8, 0xbf386671c18b088a),
        a6: (0x3bd953433969cbcd, 0xbf39e2504ddc7608),
        a7: (0xbb5c2ac0a4b9cdf1, 0x3ee34ccbca0c75d1),
        a8: (0x3b75c19030dc1713, 0x3eda4973784d1087),
        a9: (0x3ad4bc4c86170bfe, 0xbe81045322aaab4a),
        a10: (0xbb13786263c018b0, 0xbe70fae0da6cdcd7),
        a11: (0x3abaa5d08cb5c25a, 0x3e13546cef5ed920),
        a12: (0xba977aca5508657c, 0x3dfe5ee82e66311d),
        a13: (0xba08171190487e0e, 0xbd9ec80cc8c9ee06),
        c: [
            0xbd83eb2e991f8d89,
            0x3d222bfce889a6f5,
            0x3d03fb3323985824,
            0xbca0901dbb7bdc76,
            0xbc7fa628b1caa3aa,
            0x3c1811d52a833825,
            0x3bf41324eaf0514d,
            0xbb9036caa886ca4f,
        ],
    },
    J1TaylorExtendedSeries {
        a0: (0xb64cb9ce523201ef, 0xb9daf540f249f94a),
        a1: (0x3c62de114375fe8e, 0xbfd9c6cf582cbf7f),
        a2: (0xbc46b8d732eceb5f, 0x3faae8a39f51ad04),
        a3: (0xbc47767d9604e787, 0x3fab589d1da13905),
        a4: (0x3c0e65e36b50ffb5, 0xbf7537544c331da7),
        a5: (0x3bc117be8d01b8c3, 0xbf624b3409959064),
        a6: (0x3bb8c9b78bf08f74, 0x3f26e4c2d5354224),
        a7: (0x3b98631fccc89483, 0x3f083a06e30c4109),
        a8: (0x3b4c8bd49aa9c434, 0xbec9799d4c9f2549),
        a9: (0x3b45cfcd887e5a3c, 0xbea33825cd2e2c17),
        a10: (0xbaf38ba1c0748318, 0x3e617069233e9176),
        a11: (0xbab14fe958eb5a01, 0x3e34569b22afc436),
        a12: (0x3a9b97cd69610928, 0xbdf03b9e9651149b),
        a13: (0x3a3c281f22098b97, 0xbdbec62310b0d355),
        c: [
            0x3d75ec84e4973333,
            0x3d417a40cacddff7,
            0xbcf67cb1ea6a9236,
            0xbcbee80391cd2a6e,
            0x3c721fbe7ad45711,
            0x3c35e32d70a81990,
            0xbbe77df9940f194c,
            0xbbab09312e642465,
        ],
    },
    J1TaylorExtendedSeries {
        a0: (0xbc782d6271607093, 0xbfd626ee83500bf2),
        a1: (0xb69e3afb79dcb7aa, 0x3a06536578ebe315),
        a2: (0x3c6ae8952e6d3654, 0x3fc55f6bec9ef962),
        a3: (0x3c0d30f199360dee, 0xbf83d23336fd10e4),
        a4: (0x3c2695c35f7ccef3, 0xbf88c77a983a0814),
        a5: (0x3bea120604a003f7, 0x3f45cdc98db1cbe2),
        a6: (0x3bd34310c08689ff, 0x3f373576ff46ee3b),
        a7: (0x3b8a79af191763da, 0xbef2461447d7b423),
        a8: (0x3b6b01ceb89294ad, 0xbed7b853456b6eaa),
        a9: (0x3b31e9100f5ed49b, 0x3e90abfc68274a98),
        a10: (0xbb0b5dfcb96c18d2, 0x3e6ea7a1ee26124a),
        a11: (0x3acfcc4a3734436f, 0xbe235c0413e014df),
        a12: (0xba93659fab408697, 0xbdfb5c5d512fb37f),
        a13: (0xba4e307557ca6cc6, 0x3daf4c5e26ffd93c),
        c: [
            0x3d81e4c4338cc63d,
            0xbd32addefeae5a4c,
            0xbd01e4fac74c45c8,
            0x3cb12a0e65886769,
            0x3c7c4215974fc917,
            0xbc2912b2f128cffc,
            0xbbf1f0562e313f57,
            0x3b9fe6f4dbac33bc,
        ],
    },
    J1TaylorExtendedSeries {
        a0: (0xb6085a5c1dff2f55, 0xb9a6d878548bb9a8),
        a1: (0x3c7af22d033efcfb, 0x3fd33518b3874e8a),
        a2: (0xbc23b4d62966d2ca, 0xbf95e70dc60362bf),
        a3: (0xbc476d871641a91c, 0xbfa80c83bdeee5b0),
        a4: (0x3c03dac1f09bb62a, 0x3f69a4b292e3de42),
        a5: (0x3bd887f514c0eb20, 0x3f613fbc7d698217),
        a6: (0xbbac5217d3c34969, 0xbf207358bbdbff91),
        a7: (0xbb64a49a8dfea15e, 0xbf0796a751f89051),
        a8: (0xbb6e5f7a01ae7c88, 0x3ec4255b015aded4),
        a9: (0xbb2e3a67fb431c4e, 0x3ea3026e0ce97ab9),
        a10: (0x3af15b6c41ccc646, 0xbe5d48dcdae92f27),
        a11: (0xbad2c3d805d69130, 0xbe344639d7eeb0df),
        a12: (0xba75c8097b0b2363, 0x3dec62ccb4a322a7),
        a13: (0xba51abc92d8cb0a2, 0x3dbecae92e85f35a),
        c: [
            0xbd73bb6898c67420,
            0xbd4183edbc361f20,
            0x3cf4ae3e600e66f2,
            0x3cbefbb36a1305dc,
            0xbc70f26b98964b53,
            0xbc35ed2502282928,
            0x3be629d035de2b74,
            0x3ba9d5ac93cea8a8,
        ],
    },
    J1TaylorExtendedSeries {
        a0: (0xbc7d48dbfa0ea2ca, 0x3fd17dbf09d40d25),
        a1: (0x368535961cf8b870, 0xb9fb8b6c63a84946),
        a2: (0xbc61eb914d30abf9, 0xbfc1404bf647c28f),
        a3: (0x3c098a23a9b1bfaa, 0x3f74f4df2769f830),
        a4: (0x3c24ae93d1ab24f6, 0x3f85c6285429b66d),
        a5: (0x3bdddb1f3dfe719d, 0xbf3d68ab722881bd),
        a6: (0xbbd8cd540e6e76de, 0xbf356acb6452d860),
        a7: (0x3b61032e6a718639, 0x3eec10b47cf7ef69),
        a8: (0xbb4133e4863c23a4, 0x3ed67eaae97bbc86),
        a9: (0xbb24a9fa4de576f3, 0xbe8bb6530c63f2df),
        a10: (0xbb0a55e515929277, 0xbe6d87201e450ed8),
        a11: (0xbacceecc330fa4cc, 0x3e20f47c83ec5589),
        a12: (0x3a78c8b06ae9dd0b, 0x3dfa98331f6e9d6d),
        a13: (0xba236b7512cab83c, 0xbdac70414a24f17c),
        c: [
            0xbd817c057a51d300,
            0x3d316fea16507525,
            0x3d0189bd25fa371e,
            0xbcb05af6121bc01b,
            0xbc7bbd7da965743c,
            0x3c2845b45394bf53,
            0x3bf19d6ea53a06b9,
            0xbb9e7196ef9d4fe1,
        ],
    },
    J1TaylorExtendedSeries {
        a0: (0x363c5a7903ea3096, 0x39a79d763979ea5a),
        a1: (0xbc5052a3a2545760, 0xbfcff654544ebcd1),
        a2: (0xbc01b402d4b54951, 0x3f89223ff2c0785b),
        a3: (0xbc323a2755287955, 0x3fa4b0c5d5da6789),
        a4: (0x3bb3ea5867e2d7b6, 0xbf5f91a9ee0d2897),
        a5: (0x3bfc41f71f2b192d, 0xbf5f51c2489b9e6f),
        a6: (0x3bad7a53f6c2e7b5, 0x3f16b4c9ca0f770d),
        a7: (0xbbadd1ad776490db, 0x3f063c5475439cb2),
        a8: (0x3b5ae463dd2c77cc, 0xbebe3725daf69867),
        a9: (0x3b03fa0cc84a1b77, 0xbea25c1238b32e59),
        a10: (0xbafc1af0857f7c4b, 0x3e57486f6b9aa94b),
        a11: (0x3ac249934b9f25b4, 0x3e33e3bf2482780e),
        a12: (0x3a81d66372ad9477, 0xbde78a38a73e6ef8),
        a13: (0xba5c7d95254adc05, 0xbdbe844eb6b27265),
        c: [
            0x3d70e24abb2ebf61,
            0x3d41797e5f08b80f,
            0xbcf21fc0d245d977,
            0xbcbf0c13b54cfe38,
            0x3c6e43e02bd7ff6a,
            0x3c360881b9c90c65,
            0xbbe414c97f5b516b,
            0xbba9be9b1e319fd5,
        ],
    },
    J1TaylorExtendedSeries {
        a0: (0x3c50f8942d3f953d, 0xbfcddceb4ce1bf4a),
        a1: (0x3695c238ef734832, 0x39f40a43a4ae7bd5),
        a2: (0x3c433d5334a67ab6, 0x3fbda52116c0a640),
        a3: (0xbbe72468b449b811, 0xbf6a9da4603b67ea),
        a4: (0xbc0c5e08007404f7, 0xbf8331e74ea59ab8),
        a5: (0xbbdac452f9b0a025, 0x3f33e5cb6eba6eaa),
        a6: (0x3bd0ad74e95bb6eb, 0x3f33885fe9afa541),
        a7: (0xbb74e8d3ae497e2b, 0xbee494c0f4b0680b),
        a8: (0x3b78a1543d0c4570, 0xbed512b9d37762d7),
        a9: (0x3b2c1a679da9625e, 0x3e85a861082bfb7f),
        a10: (0x3b0d25a3fa6413c8, 0x3e6c323ea0a042bd),
        a11: (0xbaa7df9979973410, 0xbe1bcc962f7b92a2),
        a12: (0xba9a2714861162d6, 0xbdf9bc94e2f28f5a),
        a13: (0xba3b4de01c79f449, 0x3da82bc6fcfbaaac),
        c: [
            0x3d81141ce7a8a4c9,
            0xbd2e79ccb3afe225,
            0xbd013e1fa40eefa4,
            0x3cad36d3baa90935,
            0x3c7b65cea883b0da,
            0xbc260e62ac37dc69,
            0xbbf1734ea328bb81,
            0x3b9bc59d25df91b9,
        ],
    },
    J1TaylorExtendedSeries {
        a0: (0x3640392b055bd990, 0xb9a66f9fb1c1a259),
        a1: (0x3c6c8c66d2e432ad, 0x3fcbf3337873a7d8),
        a2: (0x3c25e81c4bca2558, 0xbf80c83a2d7add33),
        a3: (0x3c44192692b75256, 0xbfa251858011816b),
        a4: (0x3bd0475b67bb2a6e, 0x3f559eb160bf72d8),
        a5: (0x3b9e04f5bc4fbc86, 0x3f5c5bce33af2d77),
        a6: (0x3b7eb9780cb68cb7, 0xbf10413e306e0039),
        a7: (0xbb9e3c3d2cf898c6, 0xbf04a6704d05ad0b),
        a8: (0xbb5b2fd12c90bba5, 0x3eb6c43eedfed6c9),
        a9: (0x3b32a7ff57b2d995, 0x3ea16abd7815de74),
        a10: (0xbaf2132c8541e2be, 0xbe5257f16f5d4340),
        a11: (0xbac1952928900082, 0xbe332db1b4b2ff9f),
        a12: (0x3a84bdb0e9bbeb8e, 0x3de33acccf7bf184),
        a13: (0x3a57941f45c1c665, 0x3dbdc8f56825a2db),
        c: [
            0xbd6c6513384c9791,
            0xbd413585a9ef9985,
            0x3cef322e698ffdee,
            0x3cbec74a36a389b8,
            0xbc6a8a7ce28eac54,
            0xbc35f2f0a1f5699d,
            0x3be1e13f584dd63e,
            0x3ba99d63c4adf4a6,
        ],
    },
    J1TaylorExtendedSeries {
        a0: (0xbc2639709548a5a8, 0x3fca7f63fea81f26),
        a1: (0xb68508422eae4c65, 0xb9eefb28a0d6e2bf),
        a2: (0x3c4341d92ec91ce7, 0xbfba60afb06640cf),
        a3: (0x3c0aa0cf82550cc4, 0x3f62c1e930935d3c),
        a4: (0x3c2253174ead6b36, 0x3f814506466d7f1f),
        a5: (0x3bcf5948329a3551, 0xbf2cca8c0c0eaa3f),
        a6: (0xbbc232bfb72e29eb, 0xbf31df821cc1377e),
        a7: (0xbb6f26ec631c8982, 0x3edee8814ed0ac45),
        a8: (0x3b2da177dbef9bc7, 0x3ed3a365a4199dd1),
        a9: (0x3b2c72232430884f, 0xbe80ed2f9c3e458f),
        a10: (0x3b00f1cbad7d2e11, 0xbe6ab3b37c5271ae),
        a11: (0x3aa97f8af75bdc52, 0x3e1684d6e62b5cf4),
        a12: (0x3a92145bef2c1f7f, 0x3df8b105a5120394),
        a13: (0xba44c1b2e55b0437, 0xbda42dc5991c7a4d),
        c: [
            0xbd808d6405f0b94c,
            0x3d2a1520351409b8,
            0x3d00d7c0fe660214,
            0xbca98501cdfe00c0,
            0xbc7aec79ef103b93,
            0x3c239845fef853b1,
            0x3bf139ee9e89aab9,
            0xbb98e4f87cea7e8e,
        ],
    },
    J1TaylorExtendedSeries {
        a0: (0x36407fe2e6ff3bd1, 0x39a1b705f31b85d4),
        a1: (0x3c6e9557ccd1641f, 0xbfc925c6fca08f55),
        a2: (0x3c091bef68b0d406, 0x3f786dd32e059b0e),
        a3: (0x3c3dac1b11b6e257, 0x3fa09463bbd0367f),
        a4: (0x3be231ff7f9940b9, 0xbf4fda0298c8768b),
        a5: (0xbbf01855815db634, 0xbf59f4be60758fb1),
        a6: (0xbbac47cf40eb0c65, 0x3f0877991af9d1bb),
        a7: (0x3b8565cf445c5de2, 0x3f032cb00ee8c1f3),
        a8: (0x3b483966d1b02c7d, 0xbeb19d8ce8c35f58),
        a9: (0x3b42e138a1c0f159, 0xbea06a042fbba455),
        a10: (0xbacd5b9cb3e34577, 0x3e4d3a689e677726),
        a11: (0xbadc7e13a6da0f68, 0x3e325108c4ce2b71),
        a12: (0xba7e70d4a53137b4, 0xbddf7b8e9ab51af9),
        a13: (0x3a5202625e77df0b, 0xbdbcc40d05654d4d),
        c: [
            0x3d67cd76e2b9de9a,
            0x3d40c5877046dc59,
            0xbceaafec18f899e6,
            0xbcbe36dda2616796,
            0x3c671966bed5b9e7,
            0x3c35abafd03c5927,
            0xbbdf935b0d8ba815,
            0xbba9558d84dce902,
        ],
    },
    J1TaylorExtendedSeries {
        a0: (0x3c2a5f193800699c, 0xbfc810f50225b04b),
        a1: (0x368c535c1975e7c8, 0x39e8e57896ccd60a),
        a2: (0x3c5462bc86bdea70, 0x3fb7fdf97ac36b1f),
        a3: (0x3beaefb0c894ddd0, 0xbf5c3c256a8caa05),
        a4: (0xbc191dbda4b6912e, 0xbf7f98feb7286b47),
        a5: (0xbbcd0c657061ba1d, 0x3f25f6559e5686e2),
        a6: (0xbbd9de6e2d04eb59, 0x3f3080f57ac215af),
        a7: (0xbb6a036554dd764e, 0xbed80c51397e5eba),
        a8: (0x3b1bbb536eeaaa48, 0xbed256db543cd140),
        a9: (0xbb0b8971b2541027, 0x3e7af7598a219825),
        a10: (0xbad16e4e0b2a05c5, 0x3e69398226ca2300),
        a11: (0x3aad46c1b4813a0b, 0xbe1260985d9258f1),
        a12: (0x3a8e99ad58c40ef0, 0xbdf792bb3eea6423),
        a13: (0xb9f022e39678d430, 0x3da0d862695ece04),
        c: [
            0x3d7fe52adc0d4f43,
            0xbd26380138d03249,
            0xbd005a1bd594c094,
            0x3ca62051db57812a,
            0x3c7a4ddc723b36e9,
            0xbc213ff7b677d55a,
            0xbbf0ea55aa3f0798,
            0x3b96238007b3a2b9,
        ],
    },
    J1TaylorExtendedSeries {
        a0: (0xb6452b3310c1aaa3, 0xb9a6caa60b80c40a),
        a1: (0x3c62da0057f8558c, 0x3fc70c511227d5aa),
        a2: (0xbbfb574e4f9ef56d, 0xbf72ccb0e97558da),
        a3: (0x3c2e61277db2e337, 0xbf9e7dc08e70e99a),
        a4: (0x3ba7794d19fb3cb9, 0x3f48acdc5b058c0e),
        a5: (0xbbf340f4a2767c76, 0x3f580503724ad30a),
        a6: (0xbba1157631734811, 0xbf032ee4ca1fcafb),
        a7: (0x3ba6a5fb4a5950b0, 0xbf01e5d2836c8d99),
        a8: (0xbb3072dc0d2e739a, 0x3eac129f077bb163),
        a9: (0x3b3553c308d30f16, 0x3e9ef161591181a2),
        a10: (0xbac0ca7749962717, 0xbe47b9bb07f19f78),
        a11: (0x3ada4f160f80a6c8, 0xbe316f3937595d9f),
        a12: (0x3a7f51a981b7fab7, 0x3dda0bc8665b5467),
        a13: (0x3a16de59bc74bac4, 0x3dbba135f99ab7b7),
        c: [
            0xbd640d543d11cd05,
            0xbd403d05922f0354,
            0x3ce6db5df45427c8,
            0x3cbd745a4b8a5d8e,
            0xbc6413ea67ced0af,
            0xbc353f5df48a0964,
            0x3bdbcd2a6635d3a5,
            0x3ba8ea0871d7c2f0,
        ],
    },
    J1TaylorExtendedSeries {
        a0: (0xbc6b166d180d5a3d, 0x3fc633e7f7f05301),
        a1: (0x368603e86ddafd73, 0xb9e4e3a80f1408c5),
        a2: (0xbc100659a0043462, 0xbfb6273784c1c06e),
        a3: (0x3bfcb74bd5360255, 0x3f563ae94ade18d4),
        a4: (0x3c02a45d3e4654d8, 0x3f7d4666536c88b9),
        a5: (0x3b825e4c72e71577, 0xbf216d528345ca11),
        a6: (0x3bbb949409711381, 0xbf2ec0dcdbb7c5fe),
        a7: (0xbb7cad8c39ded90b, 0x3ed34e966b0b09f8),
        a8: (0xbb7ae00200e1ef00, 0x3ed135c64dc2d8d0),
        a9: (0xbb05ace0f7565a79, 0xbe75f7bc78b5fc2b),
        a10: (0x3af619c7e1ee8836, 0xbe67dc35b0764091),
        a11: (0xbaa0c306d5ba51c1, 0x3e0e6d697361eb15),
        a12: (0xba7a1421e480124b, 0x3df679e370497c85),
        a13: (0x3a261dcb6f5951b6, 0xbd9c595f278b788b),
        c: [
            0xbd7ea36aef38de5c,
            0x3d22fd66b6082dce,
            0x3cffa04c6bb5e53b,
            0xbca32f49448ccf69,
            0xbc7996207844d9df,
            0x3c1e4dd7ef2fcd8d,
            0x3bf0892eac74695c,
            0xbb93a2d1dcde613d,
        ],
    },
    J1TaylorExtendedSeries {
        a0: (0x3631a9a688e4b0cf, 0x39a2c4683733b8b9),
        a1: (0x3c6a47ab4241a433, 0xbfc5664e13b70622),
        a2: (0x3c04d78c24f34117, 0x3f6e16555e108dc6),
        a3: (0x3c1fe75afdbea881, 0x3f9c5e1ad9fb2f40),
        a4: (0x3be099fea32a0a37, 0xbf43d369f958e56a),
        a5: (0xbbf3ed71086574b3, 0xbf566f4ec27a96e9),
        a6: (0xbb9b31d62a43eb34, 0x3eff0de0532652d5),
        a7: (0xbb747d1a2d398962, 0x3f00cf264341409e),
        a8: (0xbb4f86de410445e2, 0xbea6f46d51e5766e),
        a9: (0xbb3844e9610225ff, 0xbe9d407f7c248d45),
        a10: (0xbac45fdf7fbe109e, 0x3e43a33cd9df668d),
        a11: (0x3ac0094085e3705a, 0x3e309901b0a816eb),
        a12: (0x3a7cc1b6217e8b33, 0xbdd5d856a5843200),
        a13: (0xba4862738bdcbef4, 0xbdba7cbcd8fc17ef),
        c: [
            0x3d610b62c2e53186,
            0x3d3f56a09dbd4c11,
            0xbce3ae682028b016,
            0xbcbc97753bc859be,
            0x3c617f9d54bc257c,
            0x3c34bbe5e3b83ba4,
            0xbbd87f62ce0c30d9,
            0xbba864df2d69a3a6,
        ],
    },
    J1TaylorExtendedSeries {
        a0: (0xbc4f6f3391278ee3, 0xbfc4b71d4ca2cc69),
        a1: (0x3678fc8afa380163, 0x39e1ef19ea49504a),
        a2: (0x3c5422c1a1a78d76, 0x3fb4ae245697fba6),
        a3: (0xbbf4ff57301398c7, 0xbf5215e4e1a5f1d6),
        a4: (0xbbf259ec364eef6d, 0xbf7b633ed6d9cf61),
        a5: (0x3bab33baee217f66, 0x3f1c7f17b4b7dbbd),
        a6: (0xbb793c006f2c175a, 0x3f2ce01b8b6aa34c),
        a7: (0x3b5f2ea4843f882a, 0xbecfced71b11e35b),
        a8: (0x3b5985e78b19c147, 0xbed03c9d5823261d),
        a9: (0xbaf52e9947953967, 0x3e724508091063b2),
        a10: (0x3ac6a4d2747f11ac, 0x3e66a2d20111e2fe),
        a11: (0x3aa3637869a41174, 0xbe0995a18f8e692c),
        a12: (0x3a9f5c05ad68d218, 0xbdf572d1a074eb44),
        a13: (0xba1ae65621b43636, 0x3d981df03c1a100b),
        c: [
            0x3d7d6895e471b380,
            0xbd205887f1aa0d32,
            0xbcfe86700acf9733,
            0x3ca0b3c82d0e6bed,
            0x3c78d254dc038740,
            0xbc1aaa95daaebae5,
            0xbbf01d49c56e1f4d,
            0x3b916dba679a9496,
        ],
    },
    J1TaylorExtendedSeries {
        a0: (0x363087e4187831b1, 0xb99cc9ce3377dcb9),
        a1: (0x3c6316f8ffd298b5, 0x3fc40f90793605bb),
        a2: (0xbbd411ad327ce8bd, 0xbf68c833077fbeae),
        a3: (0x3c051eb6f02cd242, 0xbf9aa0ce0421d1a8),
        a4: (0xbbec0fe7d45a82f8, 0x3f405fa598ef5d1d),
        a5: (0x3bff2085766194f4, 0x3f551d30d78ab526),
        a6: (0x3b92b3e8a5676a50, 0xbef9c5807675c5f6),
        a7: (0xbb7cb75ba818ed24, 0xbeffc1bbf57e3ae2),
        a8: (0x3b3a75d361f67e5c, 0x3ea32dfea2518ce6),
        a9: (0x3b1c16d7a1d01e7c, 0x3e9bc212085dcbc6),
        a10: (0x3ade9ccc9ffbcb83, 0xbe408b946d64c5bb),
        a11: (0xbac4dadc757e02ae, 0xbe2fa8f9d8da7371),
        a12: (0x3a496baee8e66ed4, 0x3dd293fe14af0cff),
        a13: (0xba5caa6057f9d34a, 0x3db96544cb75afd9),
        c: [
            0xbd5d4750746319dc,
            0xbd3e34181241fc7f,
            0x3ce112aa24315e2a,
            0x3cbbb1657668df47,
            0xbc5ea76a0640025d,
            0xbc342ca5ca3f84e1,
            0x3bd5a7656d96923a,
            0x3ba7d0424c130652,
        ],
    },
    J1TaylorExtendedSeries {
        a0: (0x3c4f5ffd01952bba, 0x3fc37dfa8f5a550a),
        a1: (0xb67e790c9f358b8b, 0xb9debd8a96fce99b),
        a2: (0xbc5c4cd2161806c6, 0xbfb3775c1a04f09c),
        a3: (0x3bd3d5629ef162bd, 0x3f4e2b4810a46c60),
        a4: (0xbc1b976f4bba7a04, 0x3f79d151a72b83a8),
        a5: (0x3bbcef84225c3faf, 0xbf17d8e5a090e4e6),
        a6: (0xbbcb7eac8c11f8f2, 0xbf2b49a6427386a0),
        a7: (0x3b593fe488c3fe75, 0x3ecac10957ddd2eb),
        a8: (0x3b3e5a3702e79695, 0x3ececa620745d3d3),
        a9: (0x3b09dfe9ec66f633, 0xbe6eefc7e795dcde),
        a10: (0x3affa4f697ae8ac6, 0xbe658c5d2a0da418),
        a11: (0xba84b7e207352397, 0x3e05d4721f44a987),
        a12: (0xba921ea422260dc7, 0x3df481ce2314a48d),
        a13: (0xba36fa5ad6a01172, 0xbd94c0d3279f6d36),
        c: [
            0xbd7c3ea7073609b4,
            0x3d1c61d5b30e5585,
            0x3cfd72bee69c1d9e,
            0xbc9d427c87f60ddb,
            0xbc780c60bc4216be,
            0x3c178f71a8ed9a2e,
            0x3bef592e2a63e3c3,
            0xbb8f065cc0074401,
        ],
    },
    J1TaylorExtendedSeries {
        a0: (0xb637fb11e489955b, 0x399ab9d07dd97ec9),
        a1: (0x3c689d1f481859c8, 0xbfc2f2072e638cf4),
        a2: (0x3c0f482572ea117b, 0x3f64df208bbd44f1),
        a3: (0xbc282c4cefffbadf, 0x3f992bb5e1e159fc),
        a4: (0x3bb5967522dbe710, 0xbf3ba181c06897cd),
        a5: (0xbbea6566f6bc4dfb, 0xbf53fe9d5baa4a3d),
        a6: (0xbb98b5b735102bc2, 0x3ef5d17602b01cac),
        a7: (0xbb84f1499db7d259, 0x3efe26d3747fe829),
        a8: (0xbb4681b52f6e9c5d, 0xbea0509768ab6ecb),
        a9: (0xbb3176c34f3b4b02, 0xbe9a70f232d9d06c),
        a10: (0xbad3b98f2d2d4767, 0x3e3c509252de33ec),
        a11: (0x3a98b1b6bf4c3e95, 0x3e2e454fee071173),
        a12: (0xba727b711bbf9f79, 0xbdd0015b062b92c0),
        a13: (0x3a5469a0cc109c78, 0xbdb860e95adf89e0),
        c: [
            0x3d59691e90d140cf,
            0x3d3d1ce79981f7ec,
            0xbcdddcd182c99d6d,
            0xbcbacd109ea1fcfc,
            0x3c5b03f342badf5c,
            0x3c3399abee987a20,
            0xbbd339ed067c0b80,
            0xbba7346da1bbfbee,
        ],
    },
    J1TaylorExtendedSeries {
        a0: (0xbc6b9fbd8965378e, 0xbfc2768d29c69936),
        a1: (0x36745d65ce201adf, 0x39daca2a5c741c91),
        a2: (0xbc592c5350db9236, 0x3fb271811730b0ef),
        a3: (0xbbcd42069434710b, 0xbf49a8df96a1225e),
        a4: (0xbbde8194d7b0f8be, 0xbf787c81cf1c6fc4),
        a5: (0x3ba09874f66ede09, 0x3f14549cdbb77978),
        a6: (0xbbcc0cd6e11aca49, 0x3f29ed2568116e19),
        a7: (0xbb639613521ed3d4, 0xbec6e4136f033ace),
        a8: (0xbb686a2390c87124, 0xbecd53330316cde7),
        a9: (0xbb0b07baaace1d5a, 0x3e6a983b5782dfcb),
        a10: (0x3b09b5373b692edb, 0x3e64952ba7c5a1d7),
        a11: (0x3aa8d322732c3776, 0xbe02df3ad6f82e89),
        a12: (0x3a9026e5be51484d, 0xbdf3a70f9a89c833),
        a13: (0xba2b43ec59c4a5a4, 0x3d920e086c18b52c),
        c: [
            0x3d7b29a554a4c0f2,
            0xbd18dbe090bcdce1,
            0xbcfc6bd9cb6d1a64,
            0x3c99ce5c24289978,
            0x3c774aa9bb8cbc36,
            0xbc14ec9a97b72924,
            0xbbee76d57c368758,
            0x3b8bb9462bafdc9d,
        ],
    },
    J1TaylorExtendedSeries {
        a0: (0x3634a6608477254f, 0xb99ca14830120445),
        a1: (0x3c51f9b16832f6f8, 0x3fc1ff5eec6a01cd),
        a2: (0xbc0f89ce0cda9130, 0xbf61e438b722c3b5),
        a3: (0x3c39a4b7b3e710ee, 0xbf97ed5fffc1c774),
        a4: (0xbbdc35da1f4ffca2, 0x3f37b7997babd9ca),
        a5: (0xbbe39da3ed21d513, 0x3f53081def9612c5),
        a6: (0x3b855304acc82cc7, 0xbef2c5f5edafc4e9),
        a7: (0xbb9dbb13656c0a49, 0xbefcc11a59e13739),
        a8: (0x3b333aefd6d59ef2, 0x3e9c2c3a1b8014a3),
        a9: (0xbb32e77a577c67ee, 0x3e9946d1dab7bd01),
        a10: (0xbaca12c850b9b161, 0xbe388db61946be57),
        a11: (0xbac0a5c8a264caef, 0xbe2d04d33be580ea),
        a12: (0xba68b88b50a6fe92, 0x3dcbe64386d2abed),
        a13: (0x3a5d3ef6acec5887, 0x3db77142e0e44c0b),
        c: [
            0xbd5645847644e34d,
            0xbd3c15e96b26243f,
            0x3cda545e32fe14cf,
            0x3cb9f0a9ae6d2b9f,
            0xbc57f711ec08fb63,
            0xbc33082924cf233e,
            0x3bd128b8eacdfbb5,
            0x3ba6973ec797e108,
        ],
    },
    J1TaylorExtendedSeries {
        a0: (0x3c54fa3fb220bfbc, 0x3fc194eba75b32f9),
        a1: (0x364c6dde47765f5a, 0xb9d7f6ae3097f5de),
        a2: (0x3c59f5fdd1333d04, 0xbfb190f7dc27362b),
        a3: (0xbbd9624469d54380, 0x3f462bb47a5c5f7f),
        a4: (0x3c106edbc92d5fb5, 0x3f7756ef20f5d2e2),
        a5: (0x3bb123163a50d25a, 0xbf1198b0ba97ecfb),
        a6: (0xbbafc9ed2aedfdbb, 0xbf28be8cf9358d55),
        a7: (0x3b346f45e74bea1c, 0x3ec3dd6f7c8cc3c0),
        a8: (0xbb2e47950796bc38, 0x3ecc09c80ee7f9af),
        a9: (0x3b060ae4dcd08be6, 0xbe6728e46a451e33),
        a10: (0xbb0940ec7c915ba4, 0xbe63b9111350861d),
        a11: (0xba9379f0dba83fae, 0x3e0080fddad62c65),
        a12: (0xba9e17ac0ccc5cb8, 0x3df2e111e88da3d0),
        a13: (0xba200bc229aea8a0, 0xbd8fbae88bdc038c),
        c: [
            0xbd7a2a4fbe8cc149,
            0x3d15f540eabd06ed,
            0x3cfb74c082bf1541,
            0xbc96eb77252542f4,
            0xbc7690d91531a3cd,
            0x3c12b027080cddb2,
            0x3bed99324cc7ae9a,
            0xbb88e343e7f8a1ca,
        ],
    },
    J1TaylorExtendedSeries {
        a0: (0xb5f27eb8425661cf, 0x39962cf43b966f44),
        a1: (0x3c6e71c482be66a3, 0xbfc12dd57bf18ada),
        a2: (0xbbe9a8a828b34c4a, 0x3f5f1e1e7f393e83),
        a3: (0x3c3286f932c278f9, 0x3f96d9afe88301fa),
        a4: (0xbbd36032a30c0f63, 0xbf34a538a482979b),
        a5: (0x3bef838dce4d1f81, 0xbf52316250b4ae37),
        a6: (0xbb94c1d517dee0e9, 0x3ef05f11577b4627),
        a7: (0x3b6772ba4ac4040b, 0x3efb86bad42fc220),
        a8: (0x3b208e81ea43285b, 0xbe98a1b3a9e92749),
        a9: (0x3b32ad3b2765f435, 0xbe983dcaf3f8fcc5),
        a10: (0xbab609d5d4d6c24f, 0x3e3589a7ca5fdce6),
        a11: (0x3acb8552cd51a5a1, 0x3e2be3ee3298bb9a),
        a12: (0xba6bd2445a7b23d2, 0xbdc8913f1d0fd9cb),
        a13: (0xba5c2559c5798254, 0xbdb695c38666083c),
        c: [
            0x3d53b25d3628107e,
            0x3d3b20c42e600774,
            0xbcd76504fcbff5dd,
            0xbcb91f513d4121da,
            0x3c5565e0edbf306d,
            0x3c327b32c2b3d24d,
            0xbbcecb287af7cca0,
            0xbba5fc9195b2cc7d,
        ],
    },
    J1TaylorExtendedSeries {
        a0: (0xbc27736b1f56b116, 0xbfc0d0d36473e98c),
        a1: (0x36737d96fd1440df, 0x39d501585ea2f23e),
        a2: (0xbc517d3bbb94e24e, 0x3fb0cda9974abe2b),
        a3: (0xbbe2a43b5d51fd62, 0xbf4367f38f201c25),
        a4: (0xbc17d13219074276, 0xbf7656b75e3c242e),
        a5: (0x3b94ba096117af23, 0x3f0ed82abf7489f1),
        a6: (0xbbcb0c2e3a6f8c80, 0x3f27b4e5b83eeb36),
        a7: (0x3b03c7eb52e929c2, 0xbec171fd0fb670e7),
        a8: (0xbb62d31f3950453a, 0xbecae62b4ad017fb),
        a9: (0xbb00f0b64258dae6, 0x3e64648495a7b49f),
        a10: (0xbb0e990274a8174b, 0x3e62f42a577135a9),
        a11: (0xba95ce40a4198307, 0xbdfd286e7fa32718),
        a12: (0xba941a4beba1849c, 0xbdf22dbcbf7697b7),
        a13: (0x3a187150250a0206, 0x3d8c222accc1ff20),
        c: [
            0x3d793fc7f6add916,
            0xbd138c7e2088ec77,
            0xbcfa8e4ecdf267e7,
            0x3c947e81e1b1564c,
            0x3c75e0beb4e50892,
            0xbc10c93c7bae0756,
            0xbbecc3a35357e730,
            0x3b8673177cbfe872,
        ],
    },
    J1TaylorExtendedSeries {
        a0: (0x3624134a8d41599c, 0xb98fcc67d9b3fa2e),
        a1: (0x3c61a13e2fee571c, 0x3fc076826cc2c191),
        a2: (0xbbcb789ff805c707, 0xbf5b62885e0070c6),
        a3: (0x3c35dbe9d71f0f81, 0xbf95e7f53001e4b1),
        a4: (0xbbddb8eb7dfede07, 0x3f322ebeb8dc2202),
        a5: (0xbbf72618e4c99ed1, 0x3f517444a7a04cd0),
        a6: (0x3b7b0405bc4e7e1e, 0xbeece06f1f1fcd7e),
        a7: (0x3b80166c21642e01, 0xbefa7006e6ad9cfe),
        a8: (0x3b0b35e09b163489, 0x3e95c42f02cf15ca),
        a9: (0x3b09161546627386, 0x3e9750ca5e1366b4),
        a10: (0x3a97a1855ba1a75a, 0xbe3314982df7ea98),
        a11: (0x3ac854e72235a416, 0xbe2aded75306b3b3),
        a12: (0xba35698aadc318e5, 0x3dc5d47847d8d6c2),
        a13: (0xba20d11b486bc8b1, 0x3db5ccf44d287881),
        c: [
            0xbd518fce3de79e02,
            0xbd3a3d6bcac93767,
            0x3cd4efbd458a8709,
            0x3cb85a4b04332e92,
            0xbc5339a997da2fdd,
            0xbc31f472e8277652,
            0x3bcbc765504a9b9d,
            0x3ba566c8c1ba2c6e,
        ],
    },
    J1TaylorExtendedSeries {
        a0: (0x3c61a6e0255395bf, 0x3fc02455675ab6d2),
        a1: (0x365367486cde7557, 0xb9d2de73ada15fb9),
        a2: (0xbc5a58b308383b05, 0xbfb021c155a720df),
        a3: (0xbbea19c0ff637b61, 0x3f412be56fc1449a),
        a4: (0x3c01bad3bc189838, 0x3f75749d556ad61c),
        a5: (0xbb75625a283720ea, 0xbf0b51f1f9bea93e),
        a6: (0x3bc81d09e4bd7902, 0xbf26c96a07e236bd),
        a7: (0xbb20f9ab32c5170e, 0x3ebef3a7abd5ac6b),
        a8: (0x3b6fdd4845ee401f, 0x3ec9e207c2574339),
        a9: (0xbb09dc04c0220da2, 0xbe6220b96eef8058),
        a10: (0xbb05927018ffdf07, 0xbe624317cb296732),
        a11: (0xba8811d8e1cb3b44, 0x3df9fc2f2cd3922d),
        a12: (0xba999a1bf13d99fe, 0x3df18ae8347e7886),
        a13: (0x3a16d3a20dc483a9, 0xbd89254042403af2),
        c: [
            0xbd78687dcda89378,
            0x3d1187909207f2bd,
            0x3cf9b835fb1ec657,
            0xbc92711ffd2e6ed8,
            0xbc753b03a8916a20,
            0x3c0e526a007563d4,
            0x3bebf7f332b975f4,
            0xbb8458e403ac3253,
        ],
    },
    J1TaylorExtendedSeries {
        a0: (0x363c1447b69d2534, 0x39950aaadeb6d336),
        a1: (0x3c5d7cc417171540, 0xbfbfa8b41711c83a),
        a2: (0x3bf6219a483ff763, 0x3f5857d3969997d1),
        a3: (0xbc395ccf34fbf4f3, 0x3f9511c6dadaaa12),
        a4: (0x3bd13cc5b308d757, 0xbf302c289dbdbd4f),
        a5: (0xbbde8aacfc00175d, 0xbf50cc2238d229f9),
        a6: (0x3b831e2617d2326d, 0x3ee9b64d5c63668f),
        a7: (0xbb967f7b56783af8, 0x3ef976fb023f0f79),
        a8: (0xbb2380291bc0d03f, 0xbe93693ba0b5ba70),
        a9: (0xbb30b121d8bfc480, 0xbe967b952987350c),
        a10: (0x3ad0db7f7cf0dd1d, 0x3e310cb79a2adda2),
        a11: (0x3ab8913f7b6e4b72, 0x3e29f2079f8e397e),
        a12: (0x3a6106519236841e, 0xbdc38d957eaa4065),
        a13: (0xba57f4822cc5b330, 0xbdb51511e93ba451),
        c: [
            0x3d4f8bb4d99f4248,
            0x3d396afe820b40ac,
            0xbcd2dc3c2dbcacd9,
            0xbcb7a1c8c9e23adb,
            0x3c515fdc9bf2a2c1,
            0x3c3174ac4fd3369d,
            0xbbc9303dcfe48c28,
            0xbba4d735960aa6a9,
        ],
    },
    J1TaylorExtendedSeries {
        a0: (0x3c50e4250a158ea6, 0xbfbf161d0c28b48c),
        a1: (0x366dc7a84e6855a3, 0x39d15e98bc4f5aba),
        a2: (0x3c421e360c403723, 0x3faf11d837aa6f64),
        a3: (0x3bd18f48bbb1671b, 0xbf3eab76da4d07a0),
        a4: (0x3c1ebe3bae7abdca, 0xbf74ab329f067aea),
        a5: (0xbb5e7fc40f37d1c8, 0x3f086ada57bc1c51),
        a6: (0x3bc03923b5382583, 0x3f25f6e78f11ab9a),
        a7: (0x3b5e67e026bc41b0, 0xbebbb271f54c8965),
        a8: (0x3b5e5c46ba1bbecf, 0xbec8f85328c26cb7),
        a9: (0xbb0a85c43724284f, 0x3e603f82aebdeac2),
        a10: (0xbb0536ca2f424d19, 0x3e61a3010279a191),
        a11: (0xba7e6dc2924bd831, 0xbdf75660809cdf7c),
        a12: (0xba952044a4ea0700, 0xbdf0f693177496cb),
        a13: (0x3a15c87ba753a5b8, 0x3d86a2e61267fa26),
        c: [
            0x3d77a2a8e0177051,
            0xbd0fa4c8e5a3e4ba,
            0xbcf8f19243174f40,
            0x3c90b166c469cad9,
            0x3c749fa3dbc72ed5,
            0xbc0b87a328df903d,
            0xbbeb36e448efb603,
            0x3b828735647bb88e,
        ],
    },
    J1TaylorExtendedSeries {
        a0: (0xb62adda9c544c77b, 0xb991f83094712e8d),
        a1: (0x3c0020b401658d22, 0x3fbe8727daa3daed),
        a2: (0xbbea4d873564af91, 0xbf55d353e2854a37),
        a3: (0x3c3361836c53a869, 0xbf94524d4813cc25),
        a4: (0xbbb707373cd38b67, 0x3f2d037574e28370),
        a5: (0x3bf7f1d7bf66c5e1, 0x3f50356bb747a763),
        a6: (0x3b73031e95498842, 0xbee7156bfccef376),
        a7: (0xbb9952364e168345, 0xbef896d7dc819faf),
        a8: (0x3b2e77ec7c88cfc7, 0x3e9172c6dadf4149),
        a9: (0xbb336a1da23e2ce5, 0x3e95baae8efc2e31),
        a10: (0x3aa978284691a008, 0xbe2eb347eb4d6930),
        a11: (0x3ab17bb519c6ec40, 0xbe291a60a72a20de),
        a12: (0x3a2cbd46db08a6d2, 0x3dc1a345a9a69456),
        a13: (0x3a51c4f5fe32a3b5, 0x3db46c56b01902c0),
        c: [
            0xbd4c84bb373836f5,
            0xbd38a83e6e4018ce,
            0x3cd11796a765f757,
            0x3cb6f56745fe6fc1,
            0xbc4f9315f38a72c8,
            0xbc30fc146b30875b,
            0x3bc6f36c610d3146,
            0x3ba44e7ef0b748c9,
        ],
    },
    J1TaylorExtendedSeries {
        a0: (0xbc5b4c98f0d3c92b, 0x3fbe0357c158b119),
        a1: (0xb65f88e59081b616, 0xb9cf55f7cb2f5d4e),
        a2: (0xbc410072ccac970b, 0xbfadffc2fc1a91f5),
        a3: (0x3bd769e0c7cacb9c, 0x3f3b9b82ae07da44),
        a4: (0xbc02d8dc24dc5d87, 0x3f73f64e05320ac6),
        a5: (0x3bafe4482dc51019, 0xbf05fe4b66cf19d9),
        a6: (0xbbc540acbf833d6d, 0xbf2539518e1b00f5),
        a7: (0xbb39451dd0346b5e, 0x3eb8f8d01c487905),
        a8: (0xbb5f6e2e9e6accde, 0x3ec825045b97e2dc),
        a9: (0x3afe7000c228e315, 0xbe5d565f3bb61deb),
        a10: (0xbb05c7db0f3b8c37, 0xbe611186586f4f6f),
        a11: (0x3a8fd2e1933f5344, 0x3df51a66915819aa),
        a12: (0xba9dc692d0babf24, 0x3df06ef52f670c53),
        a13: (0xba0719fca0438ecb, 0xbd848215e95d871f),
        c: [
            0xbd76ec8421e88d69,
            0x3d0cbaf5fdbdfeb4,
            0x3cf8393fcf2a7d1f,
            0xbc8e6234ecb09bae,
            0xbc740e3b95d09b98,
            0x3c091e0a8d886ba9,
            0x3bea8095dc7b9bd8,
            0xbb80f22969163527,
        ],
    },
    J1TaylorExtendedSeries {
        a0: (0xb5f66a067a7fde82, 0x3993f21a518bd848),
        a1: (0xbc5cb1f28997c986, 0xbfbd8293aa55d18f),
        a2: (0xbbd0e0b7130ad0e2, 0x3f53b6beb83f2596),
        a3: (0x3c36c091c5e1640e, 0x3f93a5ccbc12a67b),
        a4: (0x3bc80bb59a3454bb, 0xbf2a3765d26aa42b),
        a5: (0x3bc464656594e1e5, 0xbf4f5ab33748c215),
        a6: (0xbb86c47bc7c01607, 0x3ee4df6f1c257a5c),
        a7: (0xbb937260a9186c82, 0x3ef7cbd49c315be0),
        a8: (0x3b214f25495206a2, 0xbe8f96098cf07175),
        a9: (0x3b0c3e7ff1df25d4, 0xbe950b37dd43531f),
        a10: (0xbac9b0ccc7504385, 0x3e2bd2e6405c604e),
        a11: (0xbacaa3283c597a1e, 0x3e285530df0d4b6e),
        a12: (0xba3d58d2e3ee4d00, 0xbdc0029e2192fef7),
        a13: (0x3a5a284ba2f65dba, 0xbdb3d11aeba72cdd),
        c: [
            0x3d49ef077ddaff49,
            0x3d37f3d211caa03c,
            0xbccf2617841f5b4b,
            0xbcb6547847159407,
            0x3c4cd584a885473d,
            0x3c308a8db8f08190,
            0xbbc501eaa69f2c15,
            0xbba3ccd1bbe8fa54,
        ],
    },
    J1TaylorExtendedSeries {
        a0: (0xbc57ac02118cdbe7, 0xbfbd0b36e5737458),
        a1: (0xb6490c045585ba6e, 0x39cc5fde9ec2e219),
        a2: (0xbc494ec699a4242b, 0x3fad082ce3c6b59b),
        a3: (0x3bdbdb6b77183080, 0xbf3905d00c5e6800),
        a4: (0xbc14d8fdc887065d, 0xbf7352b073fdac7b),
        a5: (0x3ba8a3978ba6dfbb, 0x3f03f1ccfec2fc88),
        a6: (0xbb8f7838c3e2f1fd, 0x3f248d74583834bc),
        a7: (0x3b535fab9bd4f960, 0xbeb6a9ef0d896bae),
        a8: (0xbb67c5ab2b76db14, 0xbec764d9798d6a80),
        a9: (0x3ad78503bf911dad, 0x3e5aa785d6736f5c),
        a10: (0xbae2a50177f4ad94, 0x3e608cae36118cd7),
        a11: (0x3a4a618903cd3483, 0xbdf332ddfb39cd83),
        a12: (0xba7c2b7b5187846b, 0xbdefe502ff1a7ccb),
        a13: (0xba0f240552dea02a, 0x3d82afc83349b844),
        c: [
            0x3d764468c08a9e99,
            0xbd0a399e8621ba47,
            0xbcf78e094d0b7989,
            0x3c8bc9e06c7da22c,
            0x3c738636d3be3177,
            0xbc07054780e2c29c,
            0xbbe9d4c50e7ec1b1,
            0x3b7f20af6ac6f947,
        ],
    },
    J1TaylorExtendedSeries {
        a0: (0xb62d3cac6a7f5de5, 0xb98d1a03e12e50c5),
        a1: (0xbc49df1f0f8d232b, 0x3fbc96700bf039e2),
        a2: (0x3bd298b7ee5cb852, 0xbf51ec0b5de4befe),
        a3: (0xbc263977044e1b4b, 0xbf93095734a24496),
        a4: (0xbbbe43e8581af8c7, 0x3f27d74e12285cb2),
        a5: (0x3bbff97a03e9ec76, 0x3f4e636fe259352c),
        a6: (0x3b7416b81c6c6e69, 0xbee2fe11972bc0c6),
        a7: (0x3b833db120d1e2d4, 0xbef712e4d44c4a74),
        a8: (0xbb13ce2e09cc1aca, 0x3e8cc3adabae0452),
        a9: (0x3b2b4912b88de4cb, 0x3e946ad2d9cbeb5c),
        a10: (0x3acf383229e33953, 0xbe295d81ae83f613),
        a11: (0x3acd4de8f8985893, 0xbe27a02aefea3d5e),
        a12: (0xba021ff850e40131, 0x3dbd3a949a7205ce),
        a13: (0xba5f1a8512b7a21f, 0x3db341e0bb1935eb),
        c: [
            0xbd47b550a4cf80e2,
            0xbd374c654a17d4b1,
            0x3ccc8616f1d35810,
            0x3cb5be2d9ad06cf0,
            0xbc4a73e964fea49c,
            0xbc301fccb16cd771,
            0x3bc34f6eef22a6b1,
            0x3ba3521a2fc6e158,
        ],
    },
    J1TaylorExtendedSeries {
        a0: (0x3c2c279ff462a21b, 0x3fbc29ae8400a320),
        a1: (0x366610f47bb63bcf, 0xb9c98eecfc2045ce),
        a2: (0xbc47ac9bcf28dedf, 0xbfac27138da31c2b),
        a3: (0xbbbe5ab17b6e56fc, 0x3f36d141fcbea853),
        a4: (0x3c12b724e90d0d8d, 0x3f72bdc71062acd6),
        a5: (0x3baf5f2449e22519, 0xbf0231cf643ffc17),
        a6: (0x3b6e62387b32be61, 0xbf23f0bf3b3fe8be),
        a7: (0xbb4ce6fed1ca1602, 0x3eb4b05e955de175),
        a8: (0x3b6c765f7160873c, 0x3ec6b52b868fa5e2),
        a9: (0x3ad53fc9ec6d9129, 0xbe585a7aa3e84cc3),
        a10: (0xbb0fcbb1cefb5ca9, 0xbe6012d3384c915f),
        a11: (0xba8db199b09a335d, 0x3df18f8c487254c7),
        a12: (0xba7b78971366e9e1, 0x3deeffc4029f1d2f),
        a13: (0xba212be0ccbd1765, 0xbd811d5a3db029fc),
        c: [
            0xbd75a8d84cc94337,
            0x3d080dfa2957835c,
            0x3cf6eec055f75e5a,
            0xbc8987d894466509,
            0xbc7306ed66634428,
            0x3c053014a5223127,
            0x3be932f89d1b6dba,
            0xbb7cb347dfedd35d,
        ],
    },
    J1TaylorExtendedSeries {
        a0: (0x362e91e6c398724b, 0x3987eba7e0e8186c),
        a1: (0x3c58fff451519210, 0xbfbbbf246914235f),
        a2: (0xbbb59d6393d49e4b, 0x3f5062daee35411a),
        a3: (0x3c1bef9e89608542, 0x3f927a96f174b6d1),
        a4: (0x3bc79c690c414ee0, 0xbf25cdb5dea9c121),
        a5: (0x3bd1f78099dcc721, 0xbf4d818348f98a0f),
        a6: (0xbb8cf9ee9a454625, 0x3ee160aab829409d),
        a7: (0xbb9dfe77a45d0635, 0x3ef6698d6ee99eb9),
        a8: (0xbb05cf85193f4fa1, 0xbe8a5633d8f0b3bf),
        a9: (0x3b3185e284b4d378, 0xbe93d788d61154a7),
        a10: (0xbacb14c451b5ff80, 0x3e273ec2ae0084ac),
        a11: (0x3aa3b7a88d9756f7, 0x3e26f958f6235de8),
        a12: (0x3a1b1d5116d3ecf0, 0xbdbad0939c438e6f),
        a13: (0x3a4156e749c27deb, 0xbdb2bd56309cebc2),
        c: [
            0x3d45c709d6f30288,
            0x3d36b0b8fe6422b3,
            0xbcca3ceca74bee19,
            0xbcb531b13d69d602,
            0x3c485ef0860cddbb,
            0x3c2f76dcf05959c9,
            0xbbc1d1f147956fc7,
            0xbba2de1b72e3b7b7,
        ],
    },
    J1TaylorExtendedSeries {
        a0: (0xbc4b1bd5a08c3e5e, 0xbfbb5b8273b75055),
        a1: (0xb63105d01841394f, 0x39c7a5b471ac802f),
        a2: (0xbc1906639445a0d3, 0x3fab59418c36a684),
        a3: (0x3bd4c6e49f5b24b4, 0xbf34eafeaa92aa79),
        a4: (0xbc1e8c245002eff0, 0xbf7235801af9be44),
        a5: (0x3b9bcf38558a0827, 0x3f00af9747d0be92),
        a6: (0xbbaf4147bb1d33c8, 0x3f23611db0e1566f),
        a7: (0xbb5fad412b158d31, 0xbeb2fbe414da1250),
        a8: (0xbb64d9cb0030c89d, 0xbec613ccbb9cbe59),
        a9: (0xbad5be5553b9ee6d, 0x3e565cf274e84d31),
        a10: (0x3afc3f796bc69dfa, 0x3e5f452996e3dc22),
        a11: (0x3a512effea5b7a70, 0xbdf023f5382da41b),
        a12: (0xba81075b7054928b, 0xbdee2be24fbac4cd),
        a13: (0x39e1ea217de9fe5c, 0x3d7f7ed37410e1bd),
        c: [
            0x3d75187e9972bb08,
            0xbd06293c9fc14195,
            0xbcf65a495030c123,
            0x3c878dc081b595bb,
            0x3c728fb3269eb1a7,
            0xbc0393a9c5af363c,
            0xbbe89a9c757bfc4c,
            0x3b7a8ed0787e145f,
        ],
    },
    J1TaylorExtendedSeries {
        a0: (0x3610cdbe8213145c, 0xb97fb5a853b0a6f5),
        a1: (0xbc5024304247af3a, 0x3fbaf9cb49c4f935),
        a2: (0x3bd43675d91e2701, 0xbf4e1d930b513228),
        a3: (0xbc26b1ae5fffa9b6, 0xbf91f7a8fec6eba8),
        a4: (0xbbccaafbda8d2c0d, 0x3f240a55310866fc),
        a5: (0xbbe93eb326275203, 0x3f4cb20c812fd3aa),
        a6: (0x3b7bc265bcc2c947, 0xbedff51953c6b6cc),
        a7: (0x3b7ca74622ccb1a5, 0xbef5cdc48f5d75eb),
        a8: (0xbb2102780b474378, 0x3e883b091952c721),
        a9: (0x3b38eb9445d28445, 0x3e934fb685e58ab7),
        a10: (0x3ace22d168fa72ab, 0xbe2566fc4369ab65),
        a11: (0xbaccdcd0c54dcefc, 0xbe265f0f7de2971c),
        a12: (0x3a4cf08c29789ad6, 0x3db8b61e5f9b6070),
        a13: (0xba5a9924fae47775, 0x3db24253069b727c),
        c: [
            0xbd441732d700702b,
            0xbd361fa985917b61,
            0x3cc83c171b6e658d,
            0x3cb4ae327c5c1800,
            0xbc468a4050d84c91,
            0xbc2eba0d3577ba8a,
            0x3bc0813dc8754581,
            0x3ba2708032d85bb6,
        ],
    },
    J1TaylorExtendedSeries {
        a0: (0x3c5ffacf3e2414f1, 0x3fba9e13a0db6429),
        a1: (0x3652a6dc76a56bb3, 0xb9c68da9fcc7e65e),
        a2: (0xbc305f5aef1ca77c, 0xbfaa9c1ca2161b9b),
        a3: (0xbbd5ab98a7a3bffb, 0x3f3344a09efdc635),
        a4: (0x3be8bf0db932c378, 0x3f71b82c430a2381),
        a5: (0xbb93a99eb26711dd, 0xbefebfb97bca01f2),
        a6: (0xbb9911b56b5a60fb, 0xbf22dcdb1bc1d038),
        a7: (0x3b5eacee2324f656, 0x3eb180047f0b79ae),
        a8: (0x3b443bf48034c40b, 0x3ec57eeeee84d0d0),
        a9: (0xbafb4f7732f6cb9b, 0xbe54a0c699c8318b),
        a10: (0x3af2afbf1df3e52e, 0xbe5e7594e8a2c758),
        a11: (0xba8ca658aa96d4e8, 0x3dedccbbb4c0bb4b),
        a12: (0xba8d5ef4334d02b0, 0x3ded6766337b4b85),
        a13: (0x3a1db5b1ef72a3f2, 0xbd7d1a05cea2c532),
        c: [
            0xbd74922fb4f34d01,
            0x3d047fa35d268da5,
            0x3cf5cfa05fb66f4f,
            0xbc85d00af6170d9b,
            0xbc721fe0ec577779,
            0x3c02273d92bf1f78,
            0x3be80b119507d3bb,
            0xbb78a89c0f6fb53e,
        ],
    },
    J1TaylorExtendedSeries {
        a0: (0xb62e693f0d9702c7, 0x3988c5158d2d6de4),
        a1: (0xbc55d35a88f1dee0, 0xbfba4407e04298d1),
        a2: (0x3bdf6dddbf890972, 0x3f4bcc9df0cf00b2),
        a3: (0xbc3c3cb8ccc6dc49, 0x3f917f0266db2149),
        a4: (0x3ba3d9765edce892, 0xbf2280a052234a05),
        a5: (0x3bb902b6dd06e589, 0xbf4bf2ada1f44071),
        a6: (0xbb64af9eabd78133, 0x3edd83d58032b48d),
        a7: (0x3b3d675c5fcb3900, 0x3ef53dd972d8f232),
        a8: (0xbb21587e747ea18c, 0xbe8663c1fe202028),
        a9: (0xbb264a61c81c5b01, 0xbe92d1fbf2203ff6),
        a10: (0xba7aed74bf35206a, 0x3e23c9f0b759c5ee),
        a11: (0xba9176f04b051697, 0x3e25cfe1b012696a),
        a12: (0xba17413203174503, 0xbdb6ddc079576a20),
        a13: (0x3a49d8bd4248ab6c, 0xbdb1cfd495041ff9),
        c: [
            0x3d429b7af5016d82,
            0x3d3598302aa34d44,
            0xbcc677f66f19a253,
            0xbcb432eca5d26b74,
            0x3c44ebd68b53e09d,
            0x3c2e084e7437a4f8,
            0xbbbead32feaa163b,
            0xbba208ea174a9ed2,
        ],
    },
    J1TaylorExtendedSeries {
        a0: (0x3c423404089af1e2, 0xbfb9ef3bb2213b0b),
        a1: (0x3660b480f5ff98b6, 0x39c43e4962dc90d4),
        a2: (0x3c342c5cb509d404, 0x3fa9ed82007a9a45),
        a3: (0xbbad41b78055063f, 0xbf31d2fdeeb29f8a),
        a4: (0x3c189945a84e1573, 0xbf71446866ff1b83),
        a5: (0x3b94bf1996b6da38, 0x3efc73b684f93259),
        a6: (0xbbc79f1528c73490, 0x3f22628de594b6c9),
        a7: (0xbb2c38e6aa0b7492, 0xbeb03303c1427449),
        a8: (0x3b69b28959e7c8d6, 0xbec4f51007c51087),
        a9: (0xbaf7b22a2057fb3d, 0x3e531adfa36f2214),
        a10: (0xbafc9193bec7644e, 0x3e5db4f306b1908d),
        a11: (0xba42de55730b77fa, 0xbdeb9e33598f8a5b),
        a12: (0xba87292f522499c5, 0xbdecb09ed6eca7e1),
        a13: (0x3a029caec2960f0e, 0x3d7afe118315af33),
        c: [
            0x3d7414e442942bc4,
            0xbd0307c12a5c157f,
            0xbcf54dda9620c820,
            0x3c84455dab5a0431,
            0x3c71b6d9333c6b88,
            0xbc00e39e2d8ff4ee,
            0xbbe783b87245b473,
            0x3b76f7a0097228c4,
        ],
    },
    J1TaylorExtendedSeries {
        a0: (0xb605ac43a84fd975, 0xb97d63069997c23d),
        a1: (0x3c5728ab934a24b7, 0x3fb99be744018c90),
        a2: (0xbbdb8852606e8468, 0xbf49c3f52a2af724),
        a3: (0xbc2f281d89c30e7b, 0xbf910f5ca51f98b0),
        a4: (0xbbc3c6687fed5cd4, 0x3f2126c8e8ca2766),
        a5: (0x3babb470b9714b43, 0x3f4b416f7d4fc313),
        a6: (0xbb6778fd98e14a6c, 0xbedb5e2e5580e1ce),
        a7: (0xbb988d25da9203b1, 0xbef4b862279de756),
        a8: (0xbb2cda2c362563da, 0x3e84c5071b39dc13),
        a9: (0xbb3e32e398c638e0, 0x3e925d2fc3b19021),
        a10: (0x3a7e51c8ffb9063f, 0xbe225df322279967),
        a11: (0x3ab8a93743ec8b1c, 0xbe254a971eb6fe38),
        a12: (0x3a500f899a42d03b, 0x3db53cc6c99218f9),
        a13: (0x3a2836dcd9460b2d, 0x3db164f95180d21a),
        c: [
            0xbd414b9ef3e72c2e,
            0xbd3519623b3d9bab,
            0x3cc4e726e730eb97,
            0x3cb3bf29fff61de1,
            0xbc437b8a036af53a,
            0xbc2d60cb4d55e23e,
            0x3bbc98f1afaed941,
            0x3ba1a6f3bd6ffbbc,
        ],
    },
    J1TaylorExtendedSeries {
        a0: (0xbc54096ec86381e0, 0x3fb94d3276914e51),
        a1: (0x366f27241668bff0, 0xb9c2e635ca6a2453),
        a2: (0x3c4baee1d6cb47ff, 0xbfa94bac1950e319),
        a3: (0x3bdbfe3ad953a76a, 0x3f308d4ff8f2059e),
        a4: (0xbc1cca1c3a2beb82, 0x3f70d90d29bfeecd),
        a5: (0xbb987ff7587f29eb, 0xbefa6d56162f7fb4),
        a6: (0xbbc7e6eed5339f10, 0xbf21f107da23807d),
        a7: (0x3b42dcf0b0aa1a89, 0x3eae1a626277437b),
        a8: (0x3b65c98cc8e5e4c4, 0x3ec474eafd0cc642),
        a9: (0xbaec9f27a52ab9df, 0xbe51c27144f42d78),
        a10: (0x3aeaabae89ee7a4f, 0xbe5d01999b1a4fe1),
        a11: (0x3a8aa3070f036eff, 0x3de9b014801da92a),
        a12: (0x3a82d261d1072e03, 0x3dec061740bc79e7),
        a13: (0xba118577cdd472f5, 0xbd791f8c297364fe),
        c: [
            0xbd739fb558a1bd55,
            0x3d01b9f74e4ed656,
            0x3cf4d4252fc06ec2,
            0xbc82e618e09b3cd5,
            0xbc71540a2fbf4464,
            0x3bff85bfd57cd339,
            0x3be703f6136e8b65,
            0xbb757450d93398fc,
        ],
    },
    J1TaylorExtendedSeries {
        a0: (0xb61cba1ef04e7fed, 0x3989b1835a2335f1),
        a1: (0x3c5e213a1a4b3876, 0xbfb8ffc9bd24fe08),
        a2: (0x3be5e8a5d6a378df, 0x3f47f7d46ab33721),
        a3: (0x3c2f532ddb1c6c84, 0x3f90a7a725d3fbc4),
        a4: (0xbbb3f088ed4f7fe5, 0xbf1fea1728f216b4),
        a5: (0xbbdf63ef123c231d, 0xbf4a9cac69f0ed64),
        a6: (0x3b78f66ceaed25e0, 0x3ed977f48ff1056b),
        a7: (0x3b8a770044bcebfe, 0x3ef43c2d8e698c10),
        a8: (0x3b2abb994b44965f, 0xbe8355d1a6765ea6),
        a9: (0x3b18b37f4d1cc4d3, 0xbe91f0553501d121),
        a10: (0x3a91463b99acc51c, 0x3e211b47f6a4481f),
        a11: (0x3ac2e0b6bc45d4ab, 0x3e24ce23303889a5),
        a12: (0xba50a9232c28a1e3, 0xbdb3ca98df620d5a),
        a13: (0xba53ba791a5ae954, 0xbdb100fc74652309),
        c: [
            0x3d4020f11e11f1c1,
            0x3d34a26ef210ac8b,
            0xbcc38203cf5bb232,
            0xbcb35244a5a3b87f,
            0x3c4232a86ad50912,
            0x3c2cc2b9ddd41d78,
            0xbbbabc89e8bda64e,
            0xbba14a3e2ef0249e,
        ],
    },
    J1TaylorExtendedSeries {
        a0: (0x3c4123b2f0e7d16e, 0xbfb8b67a2481077d),
        a1: (0x36545dc273417a4f, 0x39c1f801b18baa58),
        a2: (0xbc4b4ca91bf04b2b, 0x3fa8b51f21068ea2),
        a3: (0x3bb25eeafc73f1c2, 0xbf2ed935c7aefa31),
        a4: (0xbbf1e975e030950d, 0xbf707522a5037f2d),
        a5: (0x3b9aa03389e30080, 0x3ef8a196061f8bbc),
        a6: (0x3bc503aa76bc4b7f, 0x3f21874a47e3c1e3),
        a7: (0x3b251fbc851c46c7, 0xbeac10cf34c04f17),
        a8: (0xbb59f6ffb0f70e81, 0xbec3fd6c2d4fa2a4),
        a9: (0x3ac67cc02537fb0e, 0x3e50906d55522785),
        a10: (0x3af1ee728d4aea8e, 0x3e5c5a1c124dfa00),
        a11: (0x3a60982e4eb6b3b4, 0xbde7f883b31a5ffe),
        a12: (0x3a88a2534ddc2cea, 0xbdeb668cf53018d9),
        a13: (0x3a1faa5bcdefcf38, 0x3d7775372d06b112),
        c: [
            0x3d7331d871c0f9b1,
            0xbd009010257f7a34,
            0xbcf461c3e17de103,
            0x3c81abfab1cdce3a,
            0x3c70f6ee6d44ced4,
            0xbbfd803a24de064e,
            0xbbe68b389cae2e8f,
            0x3b741889ffa2f3b7,
        ],
    },
    J1TaylorExtendedSeries {
        a0: (0xb5f042d8acd9105a, 0xb97e8eb4b1d9d307),
        a1: (0x3c5b1c9821973f29, 0x3fb86e51be0a9153),
        a2: (0xbbeaa49437fad6db, 0xbf465ed1b387e5da),
        a3: (0x3c09cfc1365ec073, 0xbf9046fc5a218a86),
        a4: (0x3bbf17fba4328a6e, 0x3f1dca617fefa913),
        a5: (0xbbe411dde5129af5, 0x3f4a0300221528a7),
        a6: (0xbb3582a678d942ff, 0xbed7c7618906f1e2),
        a7: (0x3b985c5b06c73308, 0xbef3c838897d0a1e),
        a8: (0xbb1fc5e635e13516, 0x3e820ede9f9dd7dd),
        a9: (0x3b3ce12d11cca1c4, 0x3e918a94165592bb),
        a10: (0x3a8d01b092ab12ec, 0xbe1ff76205118ef6),
        a11: (0xbab001f7bc4029fd, 0xbe24599dfeef019d),
        a12: (0x3a4dc19057439a90, 0x3db2803e59981db8),
        a13: (0xba5b26a4bac34660, 0x3db0a33202fea958),
        c: [
            0xbd3e2bff735241be,
            0xbd34329cf31ab378,
            0x3cc2424863fae66d,
            0x3cb2eba63ba85f56,
            0xbc410baa23d6c743,
            0xbc2c2d5dfc001e9c,
            0x3bb91058859d6546,
            0x3ba0f269aee99d7f,
        ],
    },
    J1TaylorExtendedSeries {
        a0: (0x3c48caabfef075c0, 0x3fb829d06fee9266),
        a1: (0x366edecc9a511d57, 0xb9c13c1b4d791243),
        a2: (0xbc46f7f24518288c, 0xbfa8289a526d7785),
        a3: (0xbbbd55f5846513b9, 0x3f2cd680355c9eb6),
        a4: (0x3c19dafae9c4b6a5, 0x3f7017d70f512861),
        a5: (0x3b94ef8a721a13eb, 0xbef707978e2a0db8),
        a6: (0xbbb84f09202d6ee6, 0xbf21247ce15e7385),
        a7: (0x3b3c3327ce08b2c0, 0x3eaa3f6125485ec1),
        a8: (0xbb6516288cde1493, 0x3ec38da848401be9),
        a9: (0x3ad1a5d469d74092, 0xbe4efe39d8db4cd9),
        a10: (0x3af1f141862d206a, 0xbe5bbd40f1db0e8d),
        a11: (0xba8dea568da6e73e, 0x3de66f7d4943701e),
        a12: (0xba7a56d406c3fc16, 0x3dead0e8148219f6),
        a13: (0xb9fd434e820bb74f, 0xbd75f78595142b60),
        c: [
            0xbd72ca9bb0ee0cc6,
            0x3cff09e49e54b5f9,
            0x3cf3f60ec71dec17,
            0xbc8091d795b72914,
            0xbc709f0c94180e06,
            0x3bfbae8413b37f33,
            0x3be618f7a8b101c1,
            0xbb72decbea01812e,
        ],
    },
    J1TaylorExtendedSeries {
        a0: (0xb5fef55730e32960, 0x3981d48bd2ba605f),
        a1: (0x3c21907f595a19be, 0xbfb7e656efb009ae),
        a2: (0x3bea538229f5f57e, 0x3f44f15066f3d876),
        a3: (0xbbef8f1c97b6beec, 0x3f8fd932c26aad94),
        a4: (0xbba4e6ba919146b6, 0xbf1be460dd86a0a4),
        a5: (0x3bdda15b624650e9, 0xbf49733b591879f8),
        a6: (0x3b7b0d78875ea531, 0x3ed64488c56022e0),
        a7: (0xbb74591ce60915aa, 0x3ef35ba58bf2f993),
        a8: (0xbb2bb204bd9709ff, 0xbe80ea47bceb9a8f),
        a9: (0xbb16e85b0b6b7954, 0xbe912b327055d0e9),
        a10: (0x3aab331f171cbb12, 0x3e1df42fc4e2481a),
        a11: (0x3abc8a24e460a229, 0x3e23ec3e76876db6),
        a12: (0xba4e369f4adcefec, 0xbdb1580393f46171),
        a13: (0x3a53083bdad618ac, 0xbdb04b0353c35235),
        c: [
            0x3d3c4ca30fcd08b2,
            0x3d33c947b92592b3,
            0xbcc122c71525057e,
            0xbcb28ac70d9916a0,
            0x3c4001f6f0738cf6,
            0x3c2ba00a1c83c8d3,
            0xbbb78dfab756468a,
            0xbba09f1dfbce0606,
        ],
    },
    J1TaylorExtendedSeries {
        a0: (0xbc2f952341a43833, 0xbfb7a62320798175),
        a1: (0x3675941f9242750c, 0x39f6e89cb4d5bbfd),
        a2: (0x3c4f5aadff795699, 0x3fa7a50ca4504bb8),
        a3: (0xbbc425ed6c52749b, 0xbf2b095ccb50a68c),
        a4: (0x3c0088e4445944f0, 0xbf6f80ef11daa37a),
        a5: (0xbb8b1aad717f7ca5, 0x3ef59822dc75b064),
        a6: (0x3bc22419101fc217, 0x3f20c7e6a7c66630),
        a7: (0x3b4c28f61e1125cc, 0xbea89e00b5c358d0),
        a8: (0x3b4efe98f8abd222, 0xbec324d5238b26d0),
        a9: (0xbae1107e66fe6981, 0x3e4d13a888d5dd35),
        a10: (0x3afdd0f6004f9d6e, 0x3e5b29f941a95afe),
        a11: (0x3a34979949562d8a, 0xbde50e6ebb3d234c),
        a12: (0xba8c58781239a0f6, 0xbdea4434c6b918ac),
        a13: (0xba1f27e279f12ae8, 0x3d74a03fe2cc7383),
        c: [
            0x3d7269628f178fe1,
            0xbcfd28cc17f55c66,
            0xbcf390703d8c346a,
            0x3c7f26d530308dba,
            0x3c704bf56a18e7ad,
            0xbbfa139e57b145a9,
            0xbbe5abf667fc5202,
            0x3b74746f6c21f668,
        ],
    },
];
