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
use crate::dyadic_float::{DyadicFloat128, DyadicSign};

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

degree = 23

print(f"pub(crate) static J1_COEFFS: [J1TaylorExtendedSeries; {len(intervals)}] = [")
for i in range(0, len(intervals)):
    interval = intervals[i]
    call_sollya_on_interval(interval[0], interval[1], interval[2], degree)
    coeffs = load_coefficients(f"coefficients.txt")
    print_remez_coeffs(coeffs)
print("];")
```
**/
pub(crate) static J1_COEFFS: [[(u64, u64); 24]; 47] = [
    [
        (0x3c61f1c324453b22, 0x3fe29ea3d19f035f),
        (0x35fa0d84854cf50d, 0xb9a948d98c79544d),
        (0xbc6e3631ae170265, 0xbfca41115c5df243),
        (0xbc18acc50ffb929f, 0x3f78d1448e6fed48),
        (0x3c0e85a7cffdb986, 0x3f8c441a2f9de22b),
        (0x3bd94ca1957aa906, 0xbf386671c18b088a),
        (0x3bd950ccbdd65736, 0xbf39e2504ddc7608),
        (0xbb61b3ae884e8323, 0x3ee34ccbca0c75d1),
        (0x3b7ab8f6e0e5328c, 0x3eda4973784d1087),
        (0x3b1f40e6cc67a3d5, 0xbe81045322aaab45),
        (0x3b06ea738f9ecf6a, 0xbe70fae0da6cdcef),
        (0xbabbaad3785f91ec, 0x3e13546cef5ed00e),
        (0x3a89c3da0f7d1467, 0x3dfe5ee82e6676f1),
        (0xba262fd170d5dbad, 0xbd9ec80cc8b644d5),
        (0xba23036bcf1b50d2, 0xbd83eb2e99627fad),
        (0xb9cc7616af2b9f56, 0x3d222bfcdb211968),
        (0xb9909fcb58b04d57, 0x3d03fb337676fc98),
        (0x3946bee5d8936909, 0xbca0901290ec992d),
        (0x391f5d570e3a4795, 0xbc7fa6a8aa2582e5),
        (0xb8b23101eae00b2a, 0x3c18078252545e8a),
        (0x386d13c8640abb7d, 0x3bf44b354a191eee),
        (0xb827a697c66ddb50, 0xbb8c75d1f25f35fc),
        (0x380d5a5801772dad, 0xbb6549fce4790a31),
        (0x378246c0eb5ba71c, 0x3affdbc62c23326e),
    ],
    [
        (0x35f9dfa8e2931f8b, 0xb95730bb59760ac5),
        (0x3c62de1143765a99, 0xbfd9c6cf582cbf7f),
        (0xbc46b8d73329a70b, 0x3faae8a39f51ad04),
        (0xbc47767d9698b3c4, 0x3fab589d1da13905),
        (0x3c0e65e41f715973, 0xbf7537544c331da7),
        (0x3bc117d1587f9833, 0xbf624b3409959064),
        (0x3bb8c81b4b053820, 0x3f26e4c2d5354224),
        (0x3b9859e946b98719, 0x3f083a06e30c4109),
        (0x3b51fd9c02e20b11, 0xbec9799d4c9f2549),
        (0xbb4ffc0b40a079e5, 0xbea33825cd2e2c16),
        (0x3aeb5032c6339f4e, 0x3e617069233e916c),
        (0x3a80c496c8db027f, 0x3e34569b22afc3c8),
        (0xba9a374b068a8d41, 0xbdf03b9e9651056f),
        (0x3a4c5798959dc28e, 0xbdbec62310af5fa2),
        (0xba11dc4cca029308, 0x3d75ec84e47b7d9e),
        (0x39d2cf2c8d7c612a, 0x3d417a40c942a964),
        (0x3971ff6b4d2284f1, 0xbcf67cb1d01a03c7),
        (0xb943da1ee0003f81, 0xbcbee7ff9372e125),
        (0xb8f9e098c8ea9826, 0x3c721fb8d23ce65d),
        (0x38d89765c4f16618, 0x3c35e053c1413792),
        (0xb8575b3ab1c49578, 0xbbe79265a43bea11),
        (0xb82c7bf56bef062f, 0xbba95d6cc94f810f),
        (0xb7f5bad532a4cc34, 0x3b59332b7e688450),
        (0x37ac610ae03b7f77, 0x3b1a2c1db24ae655),
    ],
    [
        (0xbc782d627160714a, 0xbfd626ee83500bf2),
        (0x35e215832cad36da, 0x39465d6b7a56edac),
        (0x3c6ae8952e6f1d78, 0x3fc55f6bec9ef962),
        (0x3c0d30f1a30544cb, 0xbf83d23336fd10e4),
        (0x3c2695c3589d6e1d, 0xbf88c77a983a0814),
        (0x3bea1200f42f4b09, 0x3f45cdc98db1cbe2),
        (0x3bd343248aee865b, 0x3f373576ff46ee3b),
        (0x3b8a8cd3c63aa60b, 0xbef2461447d7b423),
        (0x3b6a8a9edacbe6fe, 0xbed7b853456b6eaa),
        (0xbaff3625eaafdafc, 0x3e90abfc68274a98),
        (0xbaea64c013116245, 0x3e6ea7a1ee26124d),
        (0x3acef44f1141cf97, 0xbe235c0413e01419),
        (0xba7020655e19f397, 0xbdfb5c5d512fbafe),
        (0x3a2522345bdc9b62, 0x3daf4c5e26fd6f49),
        (0xba1a9c1d9fde6bb0, 0x3d81e4c43397bb8d),
        (0x39aead7c3865323c, 0xbd32addefc4e427f),
        (0x399ef826261673f0, 0xbd01e4fadc073884),
        (0xb9584f8fae96d4c1, 0x3cb12a0b89e91563),
        (0x3911cf3d91af7eef, 0x3c7c4246886f561c),
        (0xb8932b17cb520173, 0xbc290ee1bd0b60e5),
        (0x388f8ba7ff012875, 0xbbf21116fef61a90),
        (0x38093d7080102807, 0x3b9dcbd61b5fcfcc),
        (0xb7ddad88c94d63e2, 0x3b62f7833d3a6e0a),
        (0xb7a7e960c2811c48, 0xbb0f75a3322d6adb),
    ],
    [
        (0x35d923fac9298498, 0xb943ab00450c21a5),
        (0x3c7af22d033ee0a8, 0x3fd33518b3874e8a),
        (0xbc23b4d62984701e, 0xbf95e70dc60362bf),
        (0xbc476d8715df734e, 0xbfa80c83bdeee5b0),
        (0x3c03dac20aaf9d3e, 0x3f69a4b292e3de42),
        (0x3bd887eed7bf2f4a, 0x3f613fbc7d698217),
        (0xbbac52a914699396, 0xbf207358bbdbff91),
        (0xbb64764498326fe8, 0xbf0796a751f89051),
        (0xbb6df66e2a5b80ce, 0x3ec4255b015aded4),
        (0xbb3ae5bbd6b6b7b2, 0x3ea3026e0ce97ab9),
        (0xbaf53c3132b45416, 0xbe5d48dcdae92f2c),
        (0x3ac329ab58cdf5c8, 0xbe344639d7eeb0a6),
        (0xba7227361e36b492, 0x3dec62ccb4a32eb2),
        (0x3a436c156b981e3c, 0x3dbecae92e85448e),
        (0xba15900bb4077790, 0xbd73bb6898d7381b),
        (0x39ceeaf13516632f, 0xbd4183edbb8f3baa),
        (0xb946b19e901a3f2a, 0x3cf4ae3e7e418192),
        (0x395392c9d1b10508, 0x3cbefbb1e7e7f06a),
        (0xb91d78a8cbca8fc8, 0xbc70f28d91f2bb8c),
        (0xb8c1412abdb69981, 0xbc35ec327b001a36),
        (0x3866c77527da5dfe, 0x3be65536194079a9),
        (0x37d83c37616b601d, 0x3ba9644a2a48ce1f),
        (0xb7fc2f10d93b3bda, 0xbb58069a9c48e50d),
        (0x37b28c16423dadcd, 0xbb18f8666a0fbd52),
    ],
    [
        (0xbc7d48dbfa0ea1a4, 0x3fd17dbf09d40d25),
        (0xb5f12f0e1a0ea25f, 0xb95b2c846748c809),
        (0xbc61eb914d33c2b5, 0xbfc1404bf647c28f),
        (0x3c098a23a393f866, 0x3f74f4df2769f830),
        (0x3c24ae93dcb9bf99, 0x3f85c6285429b66d),
        (0x3bdddb259c934d19, 0xbf3d68ab722881bd),
        (0xbbd8cd731b94002b, 0xbf356acb6452d860),
        (0x3b60d2c0a7987aa7, 0x3eec10b47cf7ef69),
        (0xbb3cc8bfc9cf9d54, 0x3ed67eaae97bbc86),
        (0x3b021d7ef8f06ed1, 0xbe8bb6530c63f2df),
        (0xbaf7d003a6f8a38a, 0xbe6d87201e450edd),
        (0xba5df682e9056c25, 0x3e20f47c83ec550b),
        (0x3a784cc2fbcd5b3c, 0x3dfa98331f6ea797),
        (0x3a3af02acc5d445a, 0xbdac70414a236ade),
        (0xba1ac139214727ff, 0xbd817c057a5fc937),
        (0x39abeb4e52ce6eb7, 0x3d316fea14d29892),
        (0x39a6a5fe18553faf, 0x3d0189bd3eb9f3c6),
        (0x3946aa341fe1b254, 0xbcb05af44bb23c5f),
        (0x39010e85f5ece2bd, 0xbc7bbdb45d92e4fc),
        (0xb89b4dc8ed25771b, 0x3c284360a40ed89d),
        (0x38779d6e01f88769, 0x3bf1bfb967530933),
        (0xb82198845d1335e3, 0xbb9d37a860a66a8e),
        (0xb806974ac174343c, 0xbb62a0bbacab19c8),
        (0x37906170939ea2eb, 0x3b0e5619160fd338),
    ],
    [
        (0xb60ab4d43f876701, 0xb969264876fdf3dc),
        (0xbc5052a3a2541c57, 0xbfcff654544ebcd1),
        (0xbc01b402d42eae53, 0x3f89223ff2c0785b),
        (0xbc323a2755909c5f, 0x3fa4b0c5d5da6789),
        (0x3bb3ea54acd19ff8, 0xbf5f91a9ee0d2897),
        (0x3bfc41f7f5f8cde0, 0xbf5f51c2489b9e6f),
        (0x3bad7af9dcc67129, 0x3f16b4c9ca0f770d),
        (0xbbadd341169fa322, 0x3f063c5475439cb2),
        (0x3b59f6ed48ce32ee, 0xbebe3725daf69867),
        (0x3b21f00381a6cb34, 0xbea25c1238b32e59),
        (0xbaebcfa91cc88fb7, 0x3e57486f6b9aa951),
        (0x3adb600d25e5d0fc, 0x3e33e3bf248277ee),
        (0xba5d774620c03673, 0xbde78a38a73e7c07),
        (0x3a59e464c4822ecc, 0xbdbe844eb6b211d0),
        (0x3a19b0c1a7646942, 0x3d70e24abb406b2c),
        (0x39d0a87b33c9afed, 0x3d41797e5ead05d1),
        (0x397bece1bd73373b, 0xbcf21fc0f119f3a6),
        (0x39422a1731e7d087, 0xbcbf0c12e4c83411),
        (0x390bbdff1dae57a0, 0x3c6e44232fae84a9),
        (0x38d9844b6a7d039c, 0x3c360804e71ba342),
        (0x38824f2c3cf7a869, 0xbbe43e13e137bab9),
        (0xb8314bfff0c04984, 0xbba9903092d3c1c6),
        (0xb7e5d3b1cb26c861, 0x3b56091b8d2a60e0),
        (0x37988d9bfceddedd, 0x3b18e8b26041597f),
    ],
    [
        (0x3c50f8942d3f902b, 0xbfcddceb4ce1bf4a),
        (0xb62c35de62add74d, 0x3982dff6c1e15731),
        (0x3c433d5334b42b3d, 0x3fbda52116c0a640),
        (0xbbe72468a28b4443, 0xbf6a9da4603b67ea),
        (0xbc0c5e0831771de5, 0xbf8331e74ea59ab8),
        (0xbbdac4579b1b1f5c, 0x3f33e5cb6eba6eaa),
        (0x3bd0ad97399a9309, 0x3f33885fe9afa541),
        (0xbb74d732b8ca0df8, 0xbee494c0f4b0680b),
        (0x3b783e803ffeb597, 0xbed512b9d37762d7),
        (0x3b1382127dd60cd6, 0x3e85a861082bfb7f),
        (0xbb05de78d8bef340, 0x3e6c323ea0a042c3),
        (0xbaba616d9a2d0acd, 0xbe1bcc962f7b91eb),
        (0xba9907e5c13d4bf5, 0xbdf9bc94e2f29a4f),
        (0xba4dd160d698c980, 0x3da82bc6fcfa8ee1),
        (0xba0febc61d1c4bb8, 0x3d81141ce7b77ff8),
        (0xb9b71ce5b90d51b7, 0xbd2e79ccb186eec5),
        (0xb9729b55c285b0b8, 0xbd013e1fbe002f34),
        (0x394b7b4b328dbd37, 0x3cad36d12c86586f),
        (0x3913c4bb8b85fa48, 0x3c7b66070cb3ac2c),
        (0x38cf4dafb0ab1446, 0xbc260cbaacd645d2),
        (0xb892f6df683de3de, 0xbbf19605537cf4ab),
        (0xb806112791bea534, 0x3b9aed6b864dcd18),
        (0xb7dc24510173386c, 0x3b627eca9f946e65),
        (0xb79c60e9a1146a47, 0xbb0bfef86c0d63d6),
    ],
    [
        (0x3601c162477c9abb, 0x396159aca5db6572),
        (0x3c6c8c66d2e42063, 0x3fcbf3337873a7d8),
        (0x3c25e81c4baa875d, 0xbf80c83a2d7add33),
        (0x3c44192692d7c60a, 0xbfa251858011816b),
        (0x3bd0475c48fd4015, 0x3f559eb160bf72d8),
        (0x3b9e04d420af1ac6, 0x3f5c5bce33af2d77),
        (0x3b7eb4916c85fb39, 0xbf10413e306e0039),
        (0xbb9e3a42c4a3fc86, 0xbf04a6704d05ad0b),
        (0xbb5a4f57a7a75039, 0x3eb6c43eedfed6c9),
        (0x3b2d30c35715c7fa, 0x3ea16abd7815de74),
        (0xbae381d774667710, 0xbe5257f16f5d4346),
        (0xbad24f76ba3d2cac, 0xbe332db1b4b2ff8b),
        (0xba886b6bef46e0df, 0x3de33acccf7bfdcb),
        (0x3a40f79d1cb77e70, 0x3dbdc8f5682566ce),
        (0xb9f02059cbf2dc1b, 0xbd6c6513386da08e),
        (0xb9aee6c6f049704d, 0xbd413585a9b760a0),
        (0x396399b887813c07, 0x3cef322ea2ca7028),
        (0x3957137aaf1256dd, 0x3cbec749ba83ed87),
        (0xb90e8d0898b37a61, 0xbc6a8aba908ced81),
        (0xb8d53fcac77872af, 0xbc35f2abf33d3a06),
        (0xb881d148fd35588c, 0x3be206e2b69fc921),
        (0xb8490f14e83f785c, 0x3ba98eef781c5c36),
        (0xb7ea20a48a9f15a7, 0xbb53df43cc1537c9),
        (0x37b0c3f1e8cb94e8, 0xbb18d93abc409fbc),
    ],
    [
        (0xbc26397095487bbc, 0x3fca7f63fea81f26),
        (0x35c7769fcbedbba7, 0xb94670973020cb53),
        (0x3c4341d92ebaf230, 0xbfba60afb06640cf),
        (0x3c0aa0cf7ee02729, 0x3f62c1e930935d3c),
        (0x3c2253175b5c623d, 0x3f814506466d7f1f),
        (0x3bcf594f6bbbe6b1, 0xbf2cca8c0c0eaa3f),
        (0xbbc23306b101db0c, 0xbf31df821cc1377e),
        (0xbb6f426f813ba1cf, 0x3edee8814ed0ac45),
        (0x3b3530e090190f04, 0x3ed3a365a4199dd1),
        (0xbb253c7c31f63d20, 0xbe80ed2f9c3e458e),
        (0xbae6d6e3d35d7af7, 0xbe6ab3b37c5271b3),
        (0xbab5876ea3da2c96, 0x3e1684d6e62b5c66),
        (0xba61146a8712cc76, 0x3df8b105a5120ecf),
        (0x3a497e08482e4f36, 0xbda42dc5991b9d46),
        (0xba1c93b9ca69994f, 0xbd808d6405ffdf4d),
        (0xb9bb9c06669b183b, 0x3d2a152033665a46),
        (0x39a7ab21a3b60319, 0x3d00d7c118b17a82),
        (0x39456872094da304, 0xbca984ffd3ac0d7d),
        (0x3909c82687cc959e, 0xbc7aecb2b37b54c2),
        (0xb8cf42b62c996f29, 0x3c2397005ee8dd8b),
        (0xb886e35dc2991042, 0x3bf15c9c692b5bb4),
        (0xb832ecad9f92b6aa, 0xbb98431fbb78646d),
        (0x38062f0967c8bd95, 0xbb62529cf346132a),
        (0x37a9a4acee74819c, 0x3b0964837b2ece97),
    ],
    [
        (0x360120a78f21538d, 0xb97c8b54640b722a),
        (0x3c6e9557ccd17041, 0xbfc925c6fca08f55),
        (0x3c091bef692396e7, 0x3f786dd32e059b0e),
        (0x3c3dac1b118bbe70, 0x3fa09463bbd0367f),
        (0x3be231ff192e138d, 0xbf4fda0298c8768b),
        (0xbbf0185527f60d8b, 0xbf59f4be60758fb1),
        (0xbbac474076e2749b, 0x3f0877991af9d1bb),
        (0x3b85632e8431fcaa, 0x3f032cb00ee8c1f3),
        (0x3b46a07ad1a2663e, 0xbeb19d8ce8c35f58),
        (0x3b443944f6114d11, 0xbea06a042fbba455),
        (0xbaeb45186aece340, 0x3e4d3a689e677731),
        (0xbadf9a18bdd365cc, 0x3e325108c4ce2b64),
        (0xba6505bb4a19817a, 0xbddf7b8e9ab5314c),
        (0x3a5ad487a5bbd7ab, 0xbdbcc40d05652650),
        (0xb9dcf99ac25c15c5, 0x3d67cd76e2d7d943),
        (0xb9d024e4f22a308b, 0x3d40c58770232cfe),
        (0x3952f8f057a0b1c2, 0xbceaafec4cc364ff),
        (0x3933c2290c85768f, 0xbcbe36dd573e6a63),
        (0xb8dfd609484ae481, 0x3c67199e5ff69f74),
        (0x388f3830f4ff63d7, 0x3c35ab8bb9a71dac),
        (0x387ca89ca7c948bf, 0xbbdfd6f8b70e778b),
        (0xb844329e8362812f, 0xbba9598a6da0f2c3),
        (0x37fd2c2ab7de98e2, 0x3b51c436362641e5),
        (0x37b80a5354c74c14, 0x3b18a9fec4d5983d),
    ],
    [
        (0x3c2a5f1938003f60, 0xbfc810f50225b04b),
        (0x35f6ee82870d84bf, 0xb97d05564a508c39),
        (0x3c5462bc86c50e66, 0x3fb7fdf97ac36b1f),
        (0x3beaefb0d3cf3530, 0xbf5c3c256a8caa05),
        (0xbc191dbdbe49d9f5, 0xbf7f98feb7286b47),
        (0xbbcd0c6b4ea34036, 0x3f25f6559e5686e2),
        (0xbbd9de4a690cb448, 0x3f3080f57ac215af),
        (0xbb69ed0a7174b68e, 0xbed80c51397e5eba),
        (0x3ae06f2933360b4a, 0xbed256db543cd140),
        (0x3b1af76fcd5528a0, 0x3e7af7598a219824),
        (0x3b06515217062c3a, 0x3e69398226ca2305),
        (0xba795db0bb871535, 0xbe1260985d92587d),
        (0xba74a82a307b019b, 0xbdf792bb3eea6f67),
        (0x3a3817c8f7ec2b4c, 0x3da0d862695e1ab2),
        (0xba1b8749aa4dca22, 0x3d7fe52adc2ba19c),
        (0x39cda3a73a18530c, 0xbd263801377428c8),
        (0xb997a6607141db31, 0xbd005a1befd1620e),
        (0x394738fbaf9e72c3, 0x3ca620504263602e),
        (0x3919acd0b705176b, 0x3c7a4e14e3f85acf),
        (0x38bbce9a9cbe0c23, 0xbc213ef23a1a245a),
        (0x389387792bae4d33, 0xbbf10cae24ff432c),
        (0x383135b979266427, 0x3b95a41b5b0291e9),
        (0x380630df36241504, 0x3b62113e229c18af),
        (0x3779c024d092ed81, 0xbb06d663d6acddea),
    ],
    [
        (0xb5e79c1121733c06, 0xb97dd765a46686aa),
        (0x3c62da0057f84d3f, 0x3fc70c511227d5aa),
        (0xbbfb574e506cd46d, 0xbf72ccb0e97558da),
        (0x3c2e61277dedefa5, 0xbf9e7dc08e70e99a),
        (0x3ba77952d9976f3f, 0x3f48acdc5b058c0e),
        (0xbbf340f4df902276, 0x3f580503724ad30a),
        (0xbba115f68568f545, 0xbf032ee4ca1fcafb),
        (0x3ba6a66dcfa3c51f, 0xbf01e5d2836c8d99),
        (0xbb2b2738baeaafd8, 0x3eac129f077bb163),
        (0x3b33825b34f3c4e4, 0x3e9ef161591181a2),
        (0x3ae4e499d619d66a, 0xbe47b9bb07f19f82),
        (0x3ac1a055b5da7ef9, 0xbe316f3937595d96),
        (0xba7ad6f9f1cbf4a9, 0x3dda0bc8665b6876),
        (0xba4fc8c4f001bcc9, 0x3dbba135f99a9e23),
        (0xba00713e44285b86, 0xbd640d543d2cb545),
        (0xb9d8806fb4e057c8, 0xbd403d0592186aa2),
        (0x3976ebbceb30453a, 0x3ce6db5e22c00286),
        (0xb95cfcc8d365b2a3, 0x3cbd745a1f9b2cc0),
        (0xb8e84bda052a24a4, 0xbc64141c2da5de65),
        (0x38cd92252322f7fa, 0xbc353f4e8c94bba5),
        (0xb8715b904e388e03, 0x3bdc0988ad3d45d6),
        (0xb7f1437e75284c5c, 0x3ba8f98a7bde325e),
        (0xb7e5101387db5e8b, 0xbb4fa52f2ce39750),
        (0x37b59c36f27d26bd, 0xbb1859ec268d8e36),
    ],
    [
        (0xbc6b166d180d579d, 0x3fc633e7f7f05301),
        (0x35ffee1a88c4c7db, 0x39597a12d5b67dc0),
        (0xbc100659a075cf2f, 0xbfb6273784c1c06e),
        (0x3bfcb74bd087b3a0, 0x3f563ae94ade18d4),
        (0x3c02a45d712493c7, 0x3f7d4666536c88b9),
        (0x3b825e9ac0d7a01e, 0xbf216d528345ca11),
        (0x3bbb9405cb89e345, 0xbf2ec0dcdbb7c5fe),
        (0xbb7cb6decc20866e, 0x3ed34e966b0b09f8),
        (0xbb7a79f0d76f0cbb, 0x3ed135c64dc2d8d0),
        (0x3b01209dcbada35f, 0xbe75f7bc78b5fc2b),
        (0xbaf600f169603982, 0xbe67dc35b0764096),
        (0xba9e587e011af4be, 0x3e0e6d697361ea54),
        (0x3a9d43cb8e00cb1a, 0x3df679e3704987b1),
        (0xba3d2f784a54c5d4, 0xbd9c595f278a4dc9),
        (0x3a0e613d8bf8177d, 0xbd7ea36aef56e594),
        (0x39bad5b0fa90b31d, 0x3d22fd66b4e699f0),
        (0xb983e4c916a5d37c, 0x3cffa04c9f95e420),
        (0xb929ce3db1f0e449, 0xbca32f47f100cd0a),
        (0x38d1dc8a519c68bb, 0xbc7996582817ca81),
        (0x38bd0c28a3f651c5, 0x3c1e4c27a12867d7),
        (0x3893ddf51221f64b, 0x3bf0aafd1db63771),
        (0x383228e1d6fab232, 0xbb933b296a9fe2c2),
        (0xb8026ac12c9d03cc, 0xbb61bce573aebd8d),
        (0x37a6fcc6415460a7, 0x3b0476ec3eb1a5b8),
    ],
    [
        (0xb605d92a84ac636d, 0x396c0a10d4b8d049),
        (0x3c6a47ab4241a9f3, 0xbfc5664e13b70622),
        (0x3c04d78c254f378c, 0x3f6e16555e108dc6),
        (0x3c1fe75afd6ceb7e, 0x3f9c5e1ad9fb2f40),
        (0x3be099fe50ede362, 0xbf43d369f958e56a),
        (0xbbf3ed70de3ce58d, 0xbf566f4ec27a96e9),
        (0xbb9b30f0921d946f, 0x3eff0de0532652d5),
        (0xbb747f8e74699e57, 0x3f00cf264341409e),
        (0x3b4f302f945bfb6b, 0xbea6f46d51e5766f),
        (0xbb37090d5299f2bd, 0xbe9d407f7c248d45),
        (0xbaeee2a1b6297a03, 0x3e43a33cd9df6696),
        (0x3ad36edf4ce5ba78, 0x3e309901b0a816e5),
        (0xba69797c9a665948, 0xbdd5d856a58443f0),
        (0x3a581a0de548434c, 0xbdba7cbcd8fc075f),
        (0x39f4161aad5502dc, 0x3d610b62c2fd4020),
        (0x39b2429685721bb6, 0x3d3f56a09da1af17),
        (0xb986edb5e42e4df0, 0xbce3ae6849a19196),
        (0xb956860a8cce360d, 0xbcbc977524bd849e),
        (0x3904f25392c2b4b0, 0x3c617fc9c0a81a86),
        (0xb8b42258c6440740, 0x3c34bbe439c863fb),
        (0x385dba6e1d13451b, 0xbbd8b5326e15f6b7),
        (0xb84cacd442a85cbe, 0xbba87bd53addeb1f),
        (0xb7ed606f4edc8467, 0x3b4c2aa54d8d8127),
        (0xb788c3a332c0e10d, 0x3b17efec97ba8f61),
    ],
    [
        (0xbc4f6f339127993c, 0xbfc4b71d4ca2cc69),
        (0x35ca4f8d42ff3858, 0x3971bd20421445e5),
        (0x3c5422c1a1ae8e1c, 0x3fb4ae245697fba6),
        (0xbbf4ff572c18ea0c, 0xbf5215e4e1a5f1d6),
        (0xbbf259ec9aa6f76a, 0xbf7b633ed6d9cf61),
        (0x3bab33aa4933effd, 0x3f1c7f17b4b7dbbd),
        (0xbb79333bf7f0a288, 0x3f2ce01b8b6aa34c),
        (0x3b5f4e5959a7242e, 0xbecfced71b11e35b),
        (0x3b57f36535a55c97, 0xbed03c9d5823261d),
        (0xbb15ca3a5fd6e4c3, 0x3e724508091063b2),
        (0x3b0278797fb483cd, 0x3e66a2d20111e303),
        (0x3aa462d914cbe692, 0xbe0995a18f8e6888),
        (0xba92d2f68fe7bbf3, 0xbdf572d1a074f644),
        (0xba3d2febd6b25c97, 0x3d981df03c191241),
        (0xba1df51a16a51cab, 0x3d7d6895e48f3e4a),
        (0x3979fc3e189a5876, 0xbd205887f0b4463d),
        (0xb99ef66a10579222, 0xbcfe86703dca7086),
        (0xb93a92c2e3662797, 0x3ca0b3c70d522804),
        (0xb91294d0f373d163, 0x3c78d28b840bbf9e),
        (0xb8a3ac6e917d987d, 0xbc1aa928bdc7c69b),
        (0xb889057510b0376b, 0xbbf03e6bbda694dd),
        (0xb83da10a22fdbcd6, 0x3b9117421c5b8302),
        (0xb804a226000323af, 0x3b615b28f6d151e7),
        (0x378448f5a65efa22, 0xbb02540c2640944f),
    ],
    [
        (0x35f8435be9512d0c, 0x39596a809853afc4),
        (0x3c6316f8ffd294bc, 0x3fc40f90793605bb),
        (0xbbd411ad350e3915, 0xbf68c833077fbeae),
        (0x3c051eb6f09da299, 0xbf9aa0ce0421d1a8),
        (0xbbec0fe78ad65dee, 0x3f405fa598ef5d1d),
        (0x3bff2085596f93c4, 0x3f551d30d78ab526),
        (0x3b92b31b50a0ff60, 0xbef9c5807675c5f6),
        (0xbb7cb5b05cc545a3, 0xbeffc1bbf57e3ae2),
        (0x3b3cc24f7215eccc, 0x3ea32dfea2518ce6),
        (0x3b18c92945575b15, 0x3e9bc212085dcbc6),
        (0xbaecddf3334d6159, 0xbe408b946d64c5c2),
        (0x3ac0b1a7f29c826b, 0xbe2fa8f9d8da736a),
        (0xba756320d1b782d7, 0x3dd293fe14af1d0b),
        (0xba5e3bcddd76a2bf, 0x3db96544cb75a592),
        (0x39f0a916b1e33716, 0xbd5d4750748e1ec5),
        (0x3992aa2022ba01f9, 0xbd3e341812329adf),
        (0x396e5a994597409d, 0x3ce112aa494174e0),
        (0x395a9a5bb21cd5ce, 0x3cbbb1656dd4f875),
        (0x38fd92cb597fcae9, 0xbc5ea7b95ef8e08b),
        (0x38c2aa6384ba6bbe, 0xbc342cad8878d4a5),
        (0xb868a671070d345c, 0x3bd5d76ac6e68151),
        (0xb82311f5dacdef0c, 0x3ba7ec1ff36b8002),
        (0xb7db96a9ba77f2bd, 0xbb491cad22e089cd),
        (0xb7b0a5b54e7f6cac, 0xbb177445429e42ae),
    ],
    [
        (0x3c4f5ffd019535e1, 0x3fc37dfa8f5a550a),
        (0xb5aaaff7214acba6, 0xb93cb5c70d300ac5),
        (0xbc5c4cd2161ee66e, 0xbfb3775c1a04f09c),
        (0x3bd3d562913491bc, 0x3f4e2b4810a46c60),
        (0xbc1b976f331a69fb, 0x3f79d151a72b83a8),
        (0x3bbcef8b51459dfe, 0xbf17d8e5a090e4e6),
        (0xbbcb7ef1646f65dc, 0xbf2b49a6427386a0),
        (0x3b592486856820bb, 0x3ecac10957ddd2eb),
        (0x3b424317f186fd5b, 0x3ececa620745d3d3),
        (0xbaf34eef3c6c9553, 0xbe6eefc7e795dcdd),
        (0x3ae4b6734611daf2, 0xbe658c5d2a0da41d),
        (0x3aa97ce5a3d8c757, 0x3e05d4721f44a8f9),
        (0xba9b1942d77688ab, 0x3df481ce2314af57),
        (0x3a2ac654ac1c449a, 0xbd94c0d3279e9252),
        (0xba1f470c6e36c705, 0xbd7c3ea70752fc73),
        (0xb97829121cc317d9, 0x3d1c61d5b166a9e5),
        (0xb989ebd713c25b94, 0x3cfd72bf188712ae),
        (0xb9078152b026b94d, 0xbc9d427a9891c892),
        (0xb904e8f099d7cfab, 0xbc780c9634dc0029),
        (0x38b5c69bf0be9c48, 0x3c178e382c9b5ff2),
        (0x38536f9bd3face06, 0x3bef99f1b78427a0),
        (0xb828de3e89d5fd1e, 0xbb8e7352b4b18839),
        (0xb805ad010982fbed, 0xbb60f1a2b1d2a86e),
        (0x37abf1ea36b06bd9, 0x3b00708aca8b132c),
    ],
    [
        (0xb5c06f78c8d4a275, 0xb93d524b57168b9f),
        (0x3c689d1f48185c7e, 0xbfc2f2072e638cf4),
        (0x3c0f48257333a5e0, 0x3f64df208bbd44f1),
        (0xbc282c4cf012e4f5, 0x3f992bb5e1e159fc),
        (0x3bb5967313f39524, 0xbf3ba181c06897cd),
        (0xbbea6566cfb71c2e, 0xbf53fe9d5baa4a3d),
        (0xbb98b4ff32e89d05, 0x3ef5d17602b01cac),
        (0xbb84f1d77f1fffd7, 0x3efe26d3747fe829),
        (0xbb47896f2e82c323, 0xbea0509768ab6ecb),
        (0xbb30ede6356a7015, 0xbe9a70f232d9d06c),
        (0x3ade8fe4e5a1ed7e, 0x3e3c509252de33f9),
        (0x3ac6ae9c749edc67, 0x3e2e454fee07116e),
        (0x3a73b1f5bbe63372, 0xbdd0015b062ba122),
        (0xba44c8407e1ecc95, 0xbdb860e95adf8412),
        (0x39fddd4f18607329, 0x3d59691e90f7ccfc),
        (0x39d00120e60006cf, 0x3d3d1ce7997b42b4),
        (0x39682a77f2b9ac62, 0xbcdddcd1c52fea1b),
        (0xb95c20e8a1ea2f57, 0xbcbacd10a03d8921),
        (0xb8f1731302036655, 0x3c5b043a4e72efe2),
        (0xb8d6d22b4e2f1fae, 0x3c3399ba323db8d6),
        (0x385b209322c8d73d, 0xbbd364e661355e9a),
        (0x3838be7c1f6126ff, 0xbba75381cf058174),
        (0x37d1ce55ad6a8310, 0x3b46755328d7d302),
        (0xb7b75f604bab0893, 0x3b16ee50ed715cef),
    ],
    [
        (0xbc6b9fbd89653a0a, 0xbfc2768d29c69936),
        (0x35ca97e1f25db68b, 0xb93ffabfa86ef843),
        (0xbc592c5350d4d817, 0x3fb271811730b0ef),
        (0xbbcd42067c35395a, 0xbf49a8df96a1225e),
        (0xbbde8196594d5f87, 0xbf787c81cf1c6fc4),
        (0x3ba098686a079be4, 0x3f14549cdbb77978),
        (0xbbcc0c93813e02ce, 0x3f29ed2568116e19),
        (0xbb638a1fa72a2206, 0xbec6e4136f033ace),
        (0xbb692b62efc42c9a, 0xbecd53330316cde7),
        (0x3af835560f2b31cd, 0x3e6a983b5782dfca),
        (0x3b0d287b84d78d98, 0x3e64952ba7c5a1dc),
        (0xba72905636bd36a9, 0xbe02df3ad6f82e0d),
        (0x3a99c062ae844d27, 0xbdf3a70f9a89d2c0),
        (0x3a3c9ffd5cdbf957, 0x3d920e086c17f618),
        (0x3a1e25413f7cf79b, 0x3d7b29a554c10cda),
        (0x39b937b07d944fb1, 0xbd18dbe08f4b3b70),
        (0x399966ea8cf6d238, 0xbcfc6bd9fc31d25f),
        (0xb9310d7ab87b2019, 0x3c99ce5a745e4984),
        (0xb918ade94c5553fb, 0x3c774addf0c26059),
        (0x38be87d79b895c6e, 0xbc14eb89eaaa0820),
        (0x388e01622132e11d, 0xbbeeb605f9932592),
        (0x37f1d8f9c78e8f53, 0x3b8b3a602349f9f9),
        (0x37fda89d6c3f2718, 0x3b6084d788afe5ac),
        (0x379093629c8ee8ea, 0xbafd92b61b88c69e),
    ],
    [
        (0xb61e3a4d721e8324, 0xb970e12b7d27cf07),
        (0x3c51f9b16832f365, 0x3fc1ff5eec6a01cd),
        (0xbc0f89ce0d1cad55, 0xbf61e438b722c3b5),
        (0x3c39a4b7b3ed5aa3, 0xbf97ed5fffc1c774),
        (0xbbdc35d9a8eaece6, 0x3f37b7997babd9ca),
        (0xbbe39da40664c597, 0x3f53081def9612c5),
        (0x3b8551b9d43f8119, 0xbef2c5f5edafc4e9),
        (0xbb9dbae69d6d1983, 0xbefcc11a59e13739),
        (0x3b35152e3ec9bcce, 0x3e9c2c3a1b8014a3),
        (0xbb333a08049d2a47, 0x3e9946d1dab7bd01),
        (0x3ad9a43f93c7a3c8, 0xbe388db61946be64),
        (0x3ac51a4e1f57e228, 0xbe2d04d33be580e8),
        (0xba5f39a3631ac3ea, 0x3dcbe64386d2c5c9),
        (0xba50e72fbf713685, 0x3db77142e0e4497d),
        (0x39f05d87c3131a89, 0xbd56458476678b0a),
        (0x39dcadde2840dbd3, 0xbd3c15e96b25b19b),
        (0xb970bb70de53f974, 0x3cda545e6eacddb5),
        (0x39573a21f1dd7482, 0x3cb9f0a9b7519cb4),
        (0xb8f26c0a88ddf73c, 0xbc57f751c3434bc5),
        (0x38dba5906a1609ba, 0xbc33083bfa2b9853),
        (0xb876257782c7fb32, 0x3bd14f52f8c40156),
        (0xb8455e6398eb5536, 0x3ba6b866ba0cd5f9),
        (0xb7dd533054211520, 0xbb442a164ec53a73),
        (0xb7b5895e0b2acd11, 0xbb1663bd8ed78f9c),
    ],
    [
        (0x3c54fa3fb220c497, 0x3fc194eba75b32f9),
        (0x36084c365974fe27, 0xb97351f4ef13fe26),
        (0x3c59f5fdd12caab1, 0xbfb190f7dc27362b),
        (0xbbd96244746b5f38, 0x3f462bb47a5c5f7f),
        (0x3c106edbe0b8c444, 0x3f7756ef20f5d2e2),
        (0x3bb1231bc38d74fa, 0xbf1198b0ba97ecfb),
        (0xbbafcaf470009ca5, 0xbf28be8cf9358d55),
        (0x3b341ae8492fe909, 0x3ec3dd6f7c8cc3c0),
        (0xbb227b6134d4cb1b, 0x3ecc09c80ee7f9af),
        (0xbb0404bf1dcee0f4, 0xbe6728e46a451e32),
        (0xbb0524de73bc4a47, 0xbe63b91113508622),
        (0xba958725ca8bba0a, 0x3e0080fddad62bf8),
        (0xba979be80a70704c, 0x3df2e111e88dae1d),
        (0x3a16d9f11e8a571a, 0xbd8fbae88bdab281),
        (0x3a026df37d7f2f86, 0xbd7a2a4fbea86012),
        (0x39b683e401b5c962, 0x3d15f540e97728a3),
        (0xb95202b2fa6547b4, 0x3cfb74c0b254e833),
        (0xb938463548f90bb3, 0xbc96eb75a8c251ef),
        (0x391ab9d89901dd33, 0xbc76910c0073cf6f),
        (0x38b57304713e6ec7, 0x3c12af37349ec30d),
        (0x388e4a0ca7b6a403, 0x3bedd6cb83761ce6),
        (0xb82fbe9670a0583a, 0xbb88745a179736e3),
        (0x37f9ae92b9c8d8b8, 0xbb60180931dceaf3),
        (0xb78cba499e17a878, 0x3afab192c9241f99),
    ],
    [
        (0x35cbb600d4239a22, 0x39251b64cfa537e4),
        (0x3c6e71c482be67bd, 0xbfc12dd57bf18ada),
        (0xbbe9a8a827c4bbb7, 0x3f5f1e1e7f393e83),
        (0x3c3286f932bea35e, 0x3f96d9afe88301fa),
        (0xbbd360330de30bf2, 0xbf34a538a482979b),
        (0x3bef838ddd50c780, 0xbf52316250b4ae37),
        (0xbb94c13fc989fce7, 0x3ef05f11577b4627),
        (0x3b6771efd22cab1c, 0x3efb86bad42fc220),
        (0x3b1a6ca937a8bc13, 0xbe98a1b3a9e92749),
        (0x3b32d767937947d7, 0xbe983dcaf3f8fcc5),
        (0x3aba8d8d67ac9dda, 0x3e3589a7ca5fdcf1),
        (0x3ac9ad46d29f3be0, 0x3e2be3ee3298bb99),
        (0xba6e973d129d8f1f, 0xbdc8913f1d0ff123),
        (0x3a203a3c6881bd79, 0xbdb695c386660814),
        (0x39fb47463477c5c8, 0x3d53b25d3647586b),
        (0xb98bb89b4c085acf, 0x3d3b20c42e642e99),
        (0x397afedf65fe5bfe, 0xbcd76505329d2147),
        (0x3958669626edb388, 0xbcb91f514b70dc66),
        (0xb8ecf324e9269490, 0x3c55661a888ce787),
        (0xb8b1401ccd948cf2, 0x3c327b48d4454d52),
        (0xb8698f6cf1c5e01c, 0xbbcf10cd4224310e),
        (0x384b17d2f0b22cb7, 0xbba61f092d4fb132),
        (0xb7db3b4d370a111b, 0x3b422f16a2b51a21),
        (0x37ace9fba95b7f40, 0x3b15d89abfbbc58c),
    ],
    [
        (0xbc27736b1f56d6fe, 0xbfc0d0d36473e98c),
        (0xb5eeec74a121164b, 0xb94d417ee84a02e1),
        (0xbc517d3bbb8e77ff, 0x3fb0cda9974abe2b),
        (0xbbe2a43b589bb30c, 0xbf4367f38f201c25),
        (0xbc17d132300354c7, 0xbf7656b75e3c242e),
        (0x3b94b9f5aa61b16a, 0x3f0ed82abf7489f1),
        (0xbbcb0bedfa6ed1b8, 0x3f27b4e5b83eeb36),
        (0x3b0620c0f4833a16, 0xbec171fd0fb670e7),
        (0xbb638b5f7ce177c8, 0xbecae62b4ad017fb),
        (0x3b0b87c8b1ba0970, 0x3e64648495a7b49e),
        (0x3b05b3200be63707, 0x3e62f42a577135ad),
        (0xba91dadfc29e6554, 0xbdfd286e7fa32656),
        (0x3a9d8beeb2cb7298, 0xbdf22dbcbf76a1c5),
        (0x3a2a43a2ea6cef43, 0x3d8c222accc0d332),
        (0xba0327331125a6e3, 0x3d793fc7f6c8caef),
        (0xb9acf6791f2a444c, 0xbd138c7e1f6709fa),
        (0x399fac78a0046126, 0xbcfa8e4efc5ad38f),
        (0x3933d41876c7d322, 0x3c947e808f83bea6),
        (0xb91b783f58eeeb68, 0x3c75e0f0594a5008),
        (0xb86c4cbf1e35e546, 0xbc10c8678c595f1c),
        (0xb87ae60cd98137bf, 0xbbecffab1af4612b),
        (0x382a6f17570f961d, 0x3b86111d5963e0eb),
        (0x37f00fe3e9f49b7c, 0x3b5f5ac4d0ed3683),
        (0x377ece1c4a712e4e, 0xbaf8306e05270fd5),
    ],
    [
        (0xb60d733c7a7e52a7, 0x3970250757e24008),
        (0x3c61a13e2fee5687, 0x3fc076826cc2c191),
        (0xbbcb789ffb6667b2, 0xbf5b62885e0070c6),
        (0x3c35dbe9d7210bbe, 0xbf95e7f53001e4b1),
        (0xbbddb8eb1d2bf603, 0x3f322ebeb8dc2202),
        (0xbbf72618e8704270, 0x3f517444a7a04cd0),
        (0x3b7b01e86b4a3da0, 0xbeece06f1f1fcd7e),
        (0x3b8016817bdb5904, 0xbefa7006e6ad9cfe),
        (0x3b13ab2699f50e0a, 0x3e95c42f02cf15ca),
        (0x3b08b6ecfe623226, 0x3e9750ca5e1366b4),
        (0xbabee952793294d7, 0xbe3314982df7eaa2),
        (0x3abefa05a9bee7a7, 0xbe2aded75306b3b3),
        (0x3a67ece2ba1179c1, 0x3dc5d47847d8ebeb),
        (0xba5c6d6993547389, 0x3db5ccf44d287a22),
        (0x39e36e4501d8e34a, 0xbd518fce3e03f8ff),
        (0xb9a320cde8d9c336, 0xbd3a3d6bcad0c437),
        (0xb976652eddf4a16e, 0x3cd4efbd765ccd7d),
        (0x392c39469b09709e, 0x3cb85a4b163d5fb1),
        (0xb8fcb83bb2002a51, 0xbc5339ddcc1bcec2),
        (0xb8d45b1a219ef1df, 0xbc31f48b401dcae3),
        (0xb86e5c376fb9eb81, 0x3bcc067f698492f8),
        (0x384a2903db3f0c8f, 0x3ba58a055bf542b9),
        (0xb7e2ee806b9979fc, 0xbb4078cceb371d66),
        (0xb7bbbf36698a9adf, 0xbb154fa9d8bb2e05),
    ],
    [
        (0x3c61a6e02553980e, 0x3fc02455675ab6d2),
        (0xb5c5f31d846909a7, 0x3920fd28a7ab34ea),
        (0xbc5a58b3083e7da5, 0xbfb021c155a720df),
        (0xbbea19c1039cb49f, 0x3f412be56fc1449a),
        (0x3c01bad3e8f49c57, 0x3f75749d556ad61c),
        (0xbb75621379d06ca1, 0xbf0b51f1f9bea93e),
        (0x3bc81ccb322d1126, 0xbf26c96a07e236bd),
        (0xbb21804e8cca606d, 0x3ebef3a7abd5ac6b),
        (0xbb6f6eedeef435f9, 0x3ec9e207c257433a),
        (0xbaf0b5ea62a8de24, 0xbe6220b96eef8058),
        (0xbad2fa9552352869, 0xbe624317cb296737),
        (0xba7c0f987894ba0a, 0x3df9fc2f2cd3917f),
        (0xba7936059699e1c2, 0x3df18ae8347e8254),
        (0x3a1573486603ca10, 0xbd892540423f2e30),
        (0xba1e19de8fd09a84, 0xbd78687dcdc2db72),
        (0xb97684a35e3ea9dc, 0x3d1187909104433c),
        (0x39630ebfb454443d, 0x3cf9b8362860195d),
        (0x38c1923b5504cc9f, 0xbc92711ece62b404),
        (0xb903c12865011bca, 0xbc753b340e1ec1f7),
        (0xb8a22b4edd218aaf, 0x3c0e50ed1c32be11),
        (0x387ce89517a6eb64, 0x3bec3274fc0c19bf),
        (0x382cd7af1eec7b55, 0xbb8401a453695584),
        (0xb7f9f56cebafb7ea, 0xbb5e8c7520273be2),
        (0xb794336346140121, 0x3af602279247901e),
    ],
    [
        (0x35eae2ce139b755d, 0x396742f004f756a7),
        (0x3c5d7cc41717159f, 0xbfbfa8b41711c83a),
        (0x3bf6219a48a24bc6, 0x3f5857d3969997d1),
        (0xbc395ccf34fc8573, 0x3f9511c6dadaaa12),
        (0x3bd13cc55aed0677, 0xbf302c289dbdbd4f),
        (0xbbde8aacf938011b, 0xbf50cc2238d229f9),
        (0x3b831f1c6dd9cc15, 0x3ee9b64d5c63668f),
        (0xbb967f7ad298da49, 0x3ef976fb023f0f79),
        (0xbb26429ad31066a5, 0xbe93693ba0b5ba70),
        (0xbb30bc57284d2f2f, 0xbe967b952987350c),
        (0x3adfd4fe0a523682, 0x3e310cb79a2addab),
        (0x3abbd4987c3e4ca2, 0x3e29f2079f8e397f),
        (0xba6750593e1e01e8, 0xbdc38d957eaa53a8),
        (0x3a5878f2c117ca3c, 0xbdb51511e93ba74c),
        (0x39d6e8ba61764a78, 0x3d4f8bb4d9d2e233),
        (0x39d8bdf5fdad66d3, 0x3d396afe82155942),
        (0x39700aa17bd0529e, 0xbcd2dc3c5a2d6062),
        (0x39264d3a8a0c9e0a, 0xbcb7a1c8dec30b78),
        (0xb8fa02e0c4da574e, 0x3c51600c1fbbd4ea),
        (0xb8dd822a2730dcdd, 0x3c3174c6425ae3bb),
        (0x38597c34fcb10a17, 0xbbc969aad5156298),
        (0xb84606c1cf335ef0, 0xbba4fad8aa607bb2),
        (0xb78af3fa0ad5b88d, 0x3b3df9a69e01b8c9),
        (0x37a075924edc374d, 0x3b14caa46ecb07e4),
    ],
    [
        (0x3c50e4250a158a23, 0xbfbf161d0c28b48c),
        (0x3611a74b58e68914, 0x3973ee78b8c079ca),
        (0x3c421e360c4c6fb3, 0x3faf11d837aa6f64),
        (0x3bd18f48c3538dba, 0xbf3eab76da4d07a0),
        (0x3c1ebe3b989625cd, 0xbf74ab329f067aea),
        (0xbb5e80c38d8ec580, 0x3f086ada57bc1c51),
        (0x3bc03960e72796ae, 0x3f25f6e78f11ab9a),
        (0x3b5e771583ebe81e, 0xbebbb271f54c8965),
        (0x3b5cfd56b457272d, 0xbec8f85328c26cb7),
        (0x3b05a90e2d074d12, 0x3e603f82aebdeac1),
        (0x3b0054f935b9aa08, 0x3e61a3010279a195),
        (0x3a5b79af494102ab, 0xbdf75660809cdedf),
        (0x3a8b7a2894daf42e, 0xbdf0f6931774a05d),
        (0x3a13e36d17ff7a01, 0x3d86a2e612670759),
        (0x3a1478f3b027f12a, 0x3d77a2a8e0311432),
        (0x39a89f64fe6a7cac, 0xbd0fa4c8e3cec9ee),
        (0x399791d5ae2f24c9, 0xbcf8f1926f3bf5c9),
        (0xb90779dd9b4f4523, 0x3c90b165b307c656),
        (0x38c874530c0e61a2, 0x3c749fd30e4702fb),
        (0xb8a3d520095f9475, 0xbc0b864b8eb17f2d),
        (0xb887f6b35d01d945, 0xbbeb6fef448f9ec9),
        (0xb82c404eab8fd764, 0x3b8238c2939c0306),
        (0xb7cf45b246f68cc2, 0x3b5dc6a5afbef1d7),
        (0x37879b3e205ae80b, 0xbaf41b15e5934cd1),
    ],
    [
        (0x35e2402b29b452e6, 0xb95669e940f0e62d),
        (0x3c0020b4016594be, 0x3fbe8727daa3daed),
        (0xbbea4d873618607e, 0xbf55d353e2854a37),
        (0x3c3361836c5324f0, 0xbf94524d4813cc25),
        (0xbbb70735fac009ce, 0x3f2d037574e28370),
        (0x3bf7f1d7c0f3582e, 0x3f50356bb747a763),
        (0x3b73015c49ea72dd, 0xbee7156bfccef376),
        (0xbb99523f50a17202, 0xbef896d7dc819faf),
        (0x3b307ed5444cfe33, 0x3e9172c6dadf4149),
        (0xbb334d67f3854d6b, 0x3e95baae8efc2e31),
        (0x3abbd2db893ce2b0, 0xbe2eb347eb4d6941),
        (0x3acc3ffb82d08806, 0xbe291a60a72a20e0),
        (0x3a62e32cc76420f9, 0x3dc1a345a9a6a5f2),
        (0xba511d514fc93073, 0x3db46c56b01906be),
        (0x39eb5948642788cd, 0xbd4c84bb3767683a),
        (0xb9c6dbefdede6e96, 0xbd38a83e6e4c14a6),
        (0xb960551193c229eb, 0x3cd11796d0057ab1),
        (0xb952ad5a14be2d5e, 0x3cb6f5675ceddd1b),
        (0xb8d70078fb2fe673, 0xbc4f936ccfdc26d8),
        (0xb8bf16cd40d79f50, 0xbc30fc2f7536e2fb),
        (0xb82fad0a6a231690, 0x3bc727e83d15ee4b),
        (0x3807b3ed34c03c50, 0x3ba472426379c0dd),
        (0xb7d0fc2af29aaf50, 0xbb3b644d13049235),
        (0xb7a0b4b41cc19800, 0xbb144a95537eaf6c),
    ],
    [
        (0xbc5b4c98f0d3c4c3, 0x3fbe0357c158b119),
        (0xb60d0c75a77e2033, 0xb967816058366b9c),
        (0xbc410072ccb8850d, 0xbfadffc2fc1a91f5),
        (0x3bd769e0c0dcd2d3, 0x3f3b9b82ae07da44),
        (0xbc02d8dbfa1e178d, 0x3f73f64e05320ac6),
        (0x3bafe44f6d4ea208, 0xbf05fe4b66cf19d9),
        (0xbbc540e87c50de38, 0xbf2539518e1b00f5),
        (0xbb397c57bc8091bf, 0x3eb8f8d01c487905),
        (0xbb5e179d1aaec556, 0x3ec825045b97e2dc),
        (0xbad3613aacf5b7ad, 0xbe5d565f3bb61dea),
        (0x3af7612fa1765a21, 0xbe611186586f4f74),
        (0xba999c78eb7ba93a, 0x3df51a669158191c),
        (0x3a8fdeaddc76e697, 0x3df06ef52f6715a9),
        (0xb9f43ddc362f0509, 0xbd848215e95caabe),
        (0xb9e3ce8533da1d00, 0xbd76ec8422019334),
        (0xb977fdbacbb88faa, 0x3d0cbaf5fc1449cc),
        (0xb994adec13653643, 0x3cf8393ffa3d864a),
        (0xb920340fd1f6eaf4, 0xbc8e6232fcacc946),
        (0xb91a15a57d08156b, 0xbc740e69a1eec644),
        (0xb882efe890d1719f, 0x3c091cd31485cd1d),
        (0xb87204ec8daa0e68, 0x3beab83a7b21fd23),
        (0x37ee5c4a3769c7c2, 0xbb80ab43d4957dfa),
        (0x37d0433b7ef63b8f, 0xbb5d0a0f23f73309),
        (0xb76e179e9a438fab, 0x3af2707876868cba),
    ],
    [
        (0x3610ef51183cb002, 0x3972fe6576610a9d),
        (0xbc5cb1f28997ca3a, 0xbfbd8293aa55d18f),
        (0xbbd0e0b711c1383a, 0x3f53b6beb83f2596),
        (0x3c36c091c5e2bd45, 0x3f93a5ccbc12a67b),
        (0x3bc80bb5067d449d, 0xbf2a3765d26aa42b),
        (0x3bc464654b3effa9, 0xbf4f5ab33748c215),
        (0xbb86c3ad3a1d34af, 0x3ee4df6f1c257a5c),
        (0xbb9372510da194bb, 0x3ef7cbd49c315be0),
        (0x3b1dfd4ff2a51cd9, 0xbe8f96098cf07175),
        (0x3b0aec441511ad05, 0xbe950b37dd43531f),
        (0x3aa7073eee76aad1, 0x3e2bd2e6405c605d),
        (0xbaba14df19ca3ca5, 0x3e285530df0d4b70),
        (0xba5940f4b5dbe91a, 0xbdc0029e21930f20),
        (0x3a40c5ac889e279f, 0xbdb3d11aeba731a1),
        (0x39b3628671f6f9fe, 0x3d49ef077e064e5f),
        (0x39cb6ce5407798b1, 0x3d37f3d211d80a08),
        (0xb947488a9811c854, 0xbccf2617ceaee07b),
        (0x39261e70e580a2ad, 0xbcb654785f893e5c),
        (0xb8d177e15b92fb85, 0x3c4cd5d45e8ed84c),
        (0x38da56fe6bb5e9e1, 0x3c308aa980f020fb),
        (0x38655215f7c6e0cd, 0xbbc53213f4bfc10c),
        (0x3846a2e45c9c7024, 0xbba3f087bcfff455),
        (0xb7bf1d6ed9e53a8a, 0x3b3922460f936d70),
        (0x37bd5d3f90200af4, 0x3b13cfffa188113e),
    ],
    [
        (0xbc57ac02118ce034, 0xbfbd0b36e5737458),
        (0x35e5dbf97ada5b8e, 0x39455ca750ec6752),
        (0xbc494ec699987d83, 0x3fad082ce3c6b59b),
        (0x3bdbdb6b7d6ce172, 0xbf3905d00c5e6800),
        (0xbc14d8fddd6666b4, 0xbf7352b073fdac7b),
        (0x3ba8a390ec67540a, 0x3f03f1ccfec2fc88),
        (0xbb8f74934a18e861, 0x3f248d74583834bc),
        (0x3b536c48c738b669, 0xbeb6a9ef0d896bae),
        (0xbb686cf2b10eb75c, 0xbec764d9798d6a80),
        (0xbaf45b48c636943c, 0x3e5aa785d6736f5c),
        (0x3b032808a4ecc777, 0x3e608cae36118cdb),
        (0x3a954183d5d06233, 0xbdf332ddfb39cd01),
        (0xba85b57db0ba4e83, 0xbdefe502ff1a8f08),
        (0x3a0e117f9b9a2106, 0x3d82afc83348eef3),
        (0x3a1f9e55699c4976, 0x3d764468c0a30d64),
        (0xb9800b377972327f, 0xbd0a399e849ce64d),
        (0xb987e7e0201e1a08, 0xbcf78e09771972ad),
        (0xb92e0978fe137b3a, 0x3c8bc9dea78f9166),
        (0x3917220a53df58e4, 0x3c738663c7562da0),
        (0x3893066b32284ca1, 0xbc07042b3fa91e95),
        (0x3880dae944b873e0, 0xbbea0b14ac68d97e),
        (0x3816034e3ccdafa3, 0x3b7e9f946caccd4b),
        (0x37ee935573e4975e, 0x3b5c56e575651cdb),
        (0x378326c82431dace, 0xbaf0f996416a3127),
    ],
    [
        (0xb5cba20478c84ae1, 0xb941c202ce623e76),
        (0xbc49df1f0f8d2107, 0x3fbc96700bf039e2),
        (0x3bd298b7ed2d3ac6, 0xbf51ec0b5de4befe),
        (0xbc26397704521dc0, 0xbf93095734a24496),
        (0xbbbe43e7480a487e, 0x3f27d74e12285cb2),
        (0x3bbff97a4e7b98d7, 0x3f4e636fe259352c),
        (0x3b74153ba07b8f20, 0xbee2fe11972bc0c6),
        (0x3b833d8792db4991, 0xbef712e4d44c4a74),
        (0xbb0f1553b4c20b5d, 0x3e8cc3adabae0452),
        (0x3b2bb2d302a59ebb, 0x3e946ad2d9cbeb5c),
        (0x3abbcd480e0362ce, 0xbe295d81ae83f621),
        (0xbab3ff4398d1910c, 0xbe27a02aefea3d60),
        (0xba59673869d1d69f, 0x3dbd3a949a722395),
        (0xba4dc28fd0eb18e2, 0x3db341e0bb193b48),
        (0xb9c0a6ee415a5306, 0xbd47b550a4f76700),
        (0x3994e75cd4e25783, 0xbd374c654a26537e),
        (0x396bf6b9da90343a, 0x3ccc86173683ecb7),
        (0xb9343315e10774e3, 0x3cb5be2db45fdd99),
        (0xb8dd50a28dc1737f, 0xbc4a7432d3c802d5),
        (0xb89104cb13cbb909, 0xbc301fe8f4a0a59a),
        (0x3862fc2b2f292a7c, 0x3bc37bcc945caecc),
        (0x383f40a66926ee1f, 0x3ba375a2264aec44),
        (0x37876edaf3b8c254, 0xbb3726c4f48ef3c9),
        (0x37bfb858dba2f507, 0xbb135b197250609d),
    ],
    [
        (0x3c2c279ff462c3be, 0x3fbc29ae8400a320),
        (0x35cbc9db3681c2a9, 0x396a5c619f751287),
        (0xbc47ac9bcf3441f8, 0xbfac27138da31c2b),
        (0xbbbe5ab192ad8423, 0x3f36d141fcbea853),
        (0x3c12b724fd73605f, 0x3f72bdc71062acd6),
        (0x3baf5f2a5e20af51, 0xbf0231cf643ffc17),
        (0x3b6e53f73e07f7ec, 0xbf23f0bf3b3fe8be),
        (0xbb4cfe2786564eab, 0x3eb4b05e955de175),
        (0x3b6d19da721e7488, 0x3ec6b52b868fa5e2),
        (0x3afd657f133da2a0, 0xbe585a7aa3e84cc3),
        (0x3afd8305110d0415, 0xbe6012d3384c9164),
        (0x3a82a12fcbf74576, 0x3df18f8c4872544f),
        (0x3a83ddc6dcd20390, 0x3deeffc4029f2f01),
        (0x3a2a5f86153ccacb, 0xbd811d5a3daf7136),
        (0x3a1939ffd2f149a3, 0xbd75a8d84ce122b1),
        (0x39948d45aa4f471b, 0x3d080dfa27f2ad6a),
        (0x397ee5f9b8f0f308, 0x3cf6eec07f0d91c7),
        (0xb92b05e0bbac3b01, 0xbc8987d6f4b58a17),
        (0xb90c392890cae7e7, 0xbc7307194fcc3abc),
        (0xb894aeaf37740080, 0x3c052f0ff8dc5a44),
        (0x3880e377a138ed87, 0x3be96804c6d0c4c3),
        (0xb7fb165b04bc999e, 0xbb7c3d1f829825de),
        (0x37fc727aeddf2fd2, 0xbb5bad09acad456f),
        (0xb74fed0f6e3c7260, 0x3aef5d6f1dc95101),
    ],
    [
        (0xb5e9f9b32cc9c472, 0xb94f1c4dc0fd7ad0),
        (0x3c58fff4515190b5, 0xbfbbbf246914235f),
        (0xbbb59d638f72d376, 0x3f5062daee35411a),
        (0x3c1bef9e896a99ca, 0x3f927a96f174b6d1),
        (0x3bc79c688e87e02d, 0xbf25cdb5dea9c121),
        (0x3bd1f78082ed604a, 0xbf4d818348f98a0f),
        (0xbb8cf93ec0b60def, 0x3ee160aab829409d),
        (0xbb9dfe5ed21cbf06, 0x3ef6698d6ee99eb9),
        (0xbb0db19c59160d9d, 0xbe8a5633d8f0b3bf),
        (0x3b3148c171a390c6, 0xbe93d788d61154a7),
        (0xbabd69ba99ca8de9, 0x3e273ec2ae0084b9),
        (0x3a4e64bbbfd099bd, 0x3e26f958f6235deb),
        (0x3a5eece34fec9144, 0xbdbad0939c43a9f7),
        (0xba3e162f724391e6, 0xbdb2bd56309cf194),
        (0xb9d6bdc444caeccb, 0x3d45c709d717e64d),
        (0x39923c04d3fffb00, 0x3d36b0b8fe7370f7),
        (0x39668435dfe7dada, 0xbcca3cece6ce50aa),
        (0x3928e245f64dcddd, 0xbcb531b157c3eb03),
        (0x38d9c8cd5a159ad6, 0x3c485f346abae02a),
        (0xb8c42e265a4503fd, 0x3c2f771606e9b677),
        (0x3869cf65db30f118, 0xbbc1faf5bae17a4f),
        (0x3843ce2b3bffa99c, 0xbba3015dbdae4304),
        (0x37d01253f07a747a, 0x3b356721e926aabc),
        (0x37ab991b20a5ed5f, 0x3b12ebe23c0e019c),
    ],
    [
        (0xbc4b1bd5a08c4697, 0xbfbb5b8273b75055),
        (0xb5fc5d98089f691c, 0xb9617b443b19e784),
        (0xbc19066393ec8b4b, 0x3fab59418c36a684),
        (0x3bd4c6e4a4b5908a, 0xbf34eafeaa92aa79),
        (0xbc1e8c2463f5dff7, 0xbf7235801af9be44),
        (0x3b9bcf2d229387a5, 0x3f00af9747d0be92),
        (0xbbaf4068b21f5f13, 0x3f23611db0e1566f),
        (0xbb5fa296ae418fd9, 0xbeb2fbe414da1250),
        (0xbb6579a7a7d2b7a9, 0xbec613ccbb9cbe59),
        (0xbafb9e935c5d4727, 0x3e565cf274e84d31),
        (0xbae9f705bc9fc3a8, 0x3e5f452996e3dc2b),
        (0x3a88db40653fa6ad, 0xbdf023f5382da3ad),
        (0xba87d1ea787e96d5, 0xbdee2be24fbad63a),
        (0x39e1374c2f3e53a0, 0x3d7f7ed3740f8d64),
        (0x39f934955bbeefc4, 0x3d75187e998a123d),
        (0x39aeb3c518f70011, 0xbd06293c9e78a94a),
        (0xb980e7fd6e1aba60, 0xbcf65a49785b93de),
        (0xb925d4c892512e3e, 0x3c878dbf031e576d),
        (0xb912bf351c2e7bd2, 0x3c728fde13655c16),
        (0xb8aaed4885e1dc5d, 0xbc0392b9e7067906),
        (0x3814e81144aa2192, 0xbbe8ce75f24e6d2f),
        (0xb81d860830f44f3e, 0x3b7a224d38043629),
        (0xb7ef6d218b7f22e1, 0x3b5b0c2b72afb1a5),
        (0xb7213d15be1c6e9a, 0xbaed12d5bc6b57ab),
    ],
    [
        (0xb61ce6cbc1f04255, 0x397036fa5f6395cb),
        (0xbc5024304247ada4, 0x3fbaf9cb49c4f935),
        (0x3bd43675d81a335c, 0xbf4e1d930b513228),
        (0xbc26b1ae60058494, 0xbf91f7a8fec6eba8),
        (0xbbccaafb65fcb92e, 0x3f240a55310866fc),
        (0xbbe93eb318fd63c7, 0x3f4cb20c812fd3aa),
        (0x3b7bc11f9b0249b7, 0xbedff51953c6b6cc),
        (0x3b7ca6d618312011, 0xbef5cdc48f5d75eb),
        (0xbb1e5d2044572b21, 0x3e883b091952c721),
        (0x3b392f30805656b3, 0x3e934fb685e58ab7),
        (0x3abdef10942d29ba, 0xbe2566fc4369ab71),
        (0x3ac87311849ba940, 0xbe265f0f7de29720),
        (0x3a589755adec28c8, 0x3db8b61e5f9b79f9),
        (0xba46ab4cdf8eab9f, 0x3db24253069b78a8),
        (0xb9d3414dfa14e877, 0xbd441732d722a86f),
        (0xb9c20d9bc0b7d20e, 0xbd361fa985a1652f),
        (0xb962a44c7096efd1, 0x3cc83c1756587290),
        (0xb94098c5fbe9b36b, 0x3cb4ae329744cc53),
        (0xb8e5a6b655add69f, 0xbc468a7f4c14db50),
        (0xb8c263b696f21145, 0xbc2eba469216c470),
        (0x386c7b9219cb3b7a, 0x3bc0a74a6d599cb0),
        (0x383259264f03f94d, 0x3ba2936d1cb0c961),
        (0xb7de00e6c53ea913, 0xbb33da81fd8c724e),
        (0x37a15c8a04c7e7ef, 0xbb128235bc0c47ff),
    ],
    [
        (0x3c5ffacf3e2418f7, 0x3fba9e13a0db6429),
        (0x360c15fcefbc6cf7, 0xb96b692a38880539),
        (0xbc305f5aef32722d, 0xbfaa9c1ca2161b9b),
        (0xbbd5ab98ac975e2e, 0x3f3344a09efdc635),
        (0x3be8bf0e555b3029, 0x3f71b82c430a2381),
        (0xbb93a994567a3c7b, 0xbefebfb97bca01f2),
        (0xbb991369e36b1830, 0xbf22dcdb1bc1d038),
        (0x3b5ea31075526d31, 0x3eb180047f0b79ae),
        (0x3b46ada0a01da18a, 0x3ec57eeeee84d0d0),
        (0xbadb2a1ce7632789, 0xbe54a0c699c8318b),
        (0x3addd02fc0b7c908, 0xbe5e7594e8a2c760),
        (0xba80e5f270287d53, 0x3dedccbbb4c0ba7f),
        (0x3a7ca797a3565ecd, 0x3ded6766337b5c91),
        (0x3a1ab23ec0eb4599, 0xbd7d1a05cea18a69),
        (0x39e1f4b5b52d14cc, 0xbd74922fb50a22bf),
        (0x395024c9150c0f64, 0x3d047fa35bf6ab34),
        (0x3985f2bf77091074, 0x3cf5cfa08701ca2f),
        (0xb9266c2a76db875b, 0xbc85d0099455bd99),
        (0xb91f43baef384b49, 0xbc72200ae9780072),
        (0xb8a3a3083fb87713, 0x3c02265fdf769ade),
        (0x38836fad8505e902, 0x3be83dc86da8bc25),
        (0x37ea0aade21b140f, 0xbb78447b5e2dbf2f),
        (0xb7f62e5f79c1ca9c, 0xbb5a73de1878cbf9),
        (0x378acaea6a6d9658, 0x3aeb08df74d22d98),
    ],
    [
        (0xb5d780da95d79b82, 0x395404a4367c1acf),
        (0xbc55d35a88f1e0a3, 0xbfba4407e04298d1),
        (0x3bdf6dddc07aba4f, 0x3f4bcc9df0cf00b2),
        (0xbc3c3cb8ccc39d2a, 0x3f917f0266db2149),
        (0x3ba3d974ad41d7db, 0xbf2280a052234a05),
        (0x3bb902b66919d00f, 0xbf4bf2ada1f44071),
        (0xbb64ad3ffb27706b, 0x3edd83d58032b48d),
        (0x3b3d6efc3f5aa589, 0x3ef53dd972d8f232),
        (0xbb230bc40c1510e8, 0xbe8663c1fe202028),
        (0xbb26dba834dbae96, 0xbe92d1fbf2203ff6),
        (0x3ac7c87693a0820a, 0x3e23c9f0b759c5f9),
        (0x3ac4688bb194b145, 0x3e25cfe1b012696d),
        (0xba597597731528c4, 0xbdb6ddc0795781e2),
        (0x3a30fd07b7c74b96, 0xbdb1cfd495042669),
        (0xb9e56d22c0bed868, 0x3d429b7af52144d6),
        (0x39c7617697a68dba, 0x3d3598302ab3a84e),
        (0xb969ff73c4a2ed6f, 0xbcc677f6a5ebe777),
        (0x3900de58e64b147d, 0xbcb432ecc1192f3c),
        (0x38ee92d01dc42743, 0x3c44ec1126d1a78e),
        (0x38c2b1157c94ede8, 0x3c2e0887dbe88798),
        (0xb7d2d49d120bbb3a, 0xbbbef4037751e14b),
        (0x3837430b053ee391, 0xbba22b7609bcdbcb),
        (0x37d8cec91a9bdbe3, 0x3b327988001f32dc),
        (0x376c8d2e6d2182d4, 0x3b121ddd67c06d53),
    ],
    [
        (0x3c423404089aea02, 0xbfb9ef3bb2213b0b),
        (0xb5c87de2f8cebd53, 0xb962e8408ad4e1a5),
        (0x3c342c5cb51f294f, 0x3fa9ed82007a9a45),
        (0xbbad41b75b8ea526, 0xbf31d2fdeeb29f8a),
        (0x3c18994595321a31, 0xbf71446866ff1b83),
        (0x3b94bf0ff8b9b23d, 0x3efc73b684f93259),
        (0xbbc79edfbf5325fb, 0x3f22628de594b6c9),
        (0xbb2befa01f9526e8, 0xbeb03303c1427449),
        (0x3b691967f42c542c, 0xbec4f51007c51087),
        (0x3af5410acbffa9a8, 0x3e531adfa36f2213),
        (0xbafc56975a6a16e6, 0x3e5db4f306b19095),
        (0x3a79c0979a0204dd, 0xbdeb9e33598f899e),
        (0xba765a4f110b219b, 0xbdecb09ed6ecb892),
        (0xb9f410899a06c92c, 0x3d7afe1183148af8),
        (0x39f78f1dd4b75a7f, 0x3d7414e442aa864a),
        (0xb99a9b988c884db4, 0xbd0307c12941fdc5),
        (0xb990c58e388b6028, 0xbcf54ddabc97b3bb),
        (0x390b34a178291a11, 0x3c84455c63029495),
        (0x391979fc805b19bd, 0x3c71b7024cd44274),
        (0x389e56ba399caf6e, 0xbc00e2d078dd0efb),
        (0xb82f4378acefe136, 0xbbe7b55bb497c4cb),
        (0x381a7fc1f4e5a681, 0x3b769ad969d0625b),
        (0xb7931840ec89fbc9, 0x3b59e3a726e2743d),
        (0x378455d827026a82, 0xbae936884740c735),
    ],
    [
        (0x35fae2734bb315aa, 0x3968400e2d1b8167),
        (0x3c5728ab934a269f, 0x3fb99be744018c90),
        (0xbbdb8852614fe955, 0xbf49c3f52a2af724),
        (0xbc2f281d89ca125d, 0xbf910f5ca51f98b0),
        (0xbbc3c6681ad02312, 0x3f2126c8e8ca2766),
        (0x3babb471b2a0c061, 0x3f4b416f7d4fc313),
        (0xbb677b339729278d, 0xbedb5e2e5580e1ce),
        (0xbb988d465ffb2680, 0xbef4b862279de756),
        (0xbb2b440cbefcb86d, 0x3e84c5071b39dc13),
        (0xbb3de6305e87339f, 0x3e925d2fc3b19021),
        (0x3ac90d431ee8f3e6, 0xbe225df322279972),
        (0xbac3c7a2f5c16d4a, 0xbe254a971eb6fe3b),
        (0x3a515ee7ac48908a, 0x3db53cc6c9922f25),
        (0x3a5ec538d819665b, 0x3db164f95180d8bf),
        (0x39d2657db4413c0c, 0xbd414b9ef404e359),
        (0x39d5a070336f534a, 0xbd3519623b4e4c39),
        (0x396e222733381912, 0x3cc4e7271a5acddf),
        (0xb94d7cdd22f0028c, 0x3cb3bf2a1b7bb45e),
        (0xb8e550d76b23a690, 0xbc437bc0b629d672),
        (0x38b920c14a99cb5c, 0xbc2d61049bf3447b),
        (0x3855956d27abf79f, 0x3bbcdb091753dc1f),
        (0x384b894b04dc498b, 0x3ba1c91964954692),
        (0x37d5bcb82226b372, 0xbb313e10f59f02a2),
        (0x37bd81d4f7329d91, 0xbb11be93ca114b15),
    ],
    [
        (0xbc54096ec8637e04, 0x3fb94d3276914e51),
        (0x35b1f799a74fb76a, 0x395d690fadf22a83),
        (0x3c4baee1d6c0d48a, 0xbfa94bac1950e319),
        (0x3bdbfe3ad50b9727, 0x3f308d4ff8f2059e),
        (0xbc1cca1c2772d45e, 0x3f70d90d29bfeecd),
        (0xbb987fee6366d20f, 0xbefa6d56162f7fb4),
        (0xbbc7e7232a191e61, 0xbf21f107da23807d),
        (0x3b42cbe0cc45eb3d, 0x3eae1a626277437b),
        (0x3b665f948e1d7ee0, 0x3ec474eafd0cc642),
        (0x3acb746c7398995a, 0xbe51c27144f42d78),
        (0x3af77d2dd1ae1ba3, 0xbe5d01999b1a4fe9),
        (0x3a77ecdd37e45371, 0x3de9b014801da87a),
        (0x3a855ec04282c8cc, 0x3dec061740bc8a41),
        (0x3a1fc05e7d7c7991, 0xbd791f8c297254d7),
        (0x3a0ace8e68225aa3, 0xbd739fb558b7a363),
        (0x398b465ba48d7de4, 0x3d01b9f74d482556),
        (0x399084cac0e412c2, 0x3cf4d425556e6c52),
        (0xb91cfaa8a45d49ae, 0xbc82e617aee470c1),
        (0xb91bd25d6d983a31, 0xbc71543271ff906d),
        (0xb89586b26f4fbcba, 0x3bff8440e63c211d),
        (0x386772bb21a29d94, 0x3be734945f8177cb),
        (0xb80904837292f3d9, 0xbb751e13261a9dc1),
        (0x37f62a65371653ca, 0xbb595b06b5ab1c41),
        (0x378aee2000820cbf, 0x3ae79434d2d5aa8d),
    ],
    [
        (0x35fd3e35d9fc508f, 0x39692a03d61a6765),
        (0x3c5e213a1a4b3671, 0xbfb8ffc9bd24fe08),
        (0x3be5e8a5d70cddcc, 0x3f47f7d46ab33721),
        (0x3c2f532ddb23da23, 0x3f90a7a725d3fbc4),
        (0xbbb3f089aa77a72c, 0xbf1fea1728f216b4),
        (0xbbdf63ef331627c3, 0xbf4a9cac69f0ed64),
        (0x3b78f775a8392238, 0x3ed977f48ff1056b),
        (0x3b8a77447708dbe0, 0x3ef43c2d8e698c10),
        (0x3b293fa288b304e6, 0xbe8355d1a6765ea6),
        (0x3b17742d486c25ae, 0xbe91f0553501d121),
        (0xba8ae50dffeaa9a9, 0x3e211b47f6a44829),
        (0xbaa73bc70607579c, 0x3e24ce23303889a9),
        (0x3a5d5eb483ce8058, 0xbdb3ca98df62221a),
        (0x3a46a518fdfd8131, 0xbdb100fc746529d6),
        (0x39e722b8c5f7aad4, 0x3d4020f11e2dc07c),
        (0xb9d9956b0c21d835, 0x3d34a26ef221978f),
        (0x396850ab1f7d4d53, 0xbcc38203ff3d18df),
        (0x3945a4bf6ed3a626, 0xbcb35244c14991ae),
        (0xb8eeb070694cf5e7, 0x3c4232db9b60cb00),
        (0xb8a56af95593329b, 0x3c2cc2f2efd1e051),
        (0xb817a34eb41964e7, 0xbbbafa643c1fee7f),
        (0x38319c576c8452aa, 0xbba16bf81ee3b787),
        (0x37d076f627ec3422, 0x3b3022fbfe27424b),
        (0xb7a801683d487848, 0x3b116413b0cae3c2),
    ],
    [
        (0x3c4123b2f0e7c9dd, 0xbfb8b67a2481077d),
        (0xb5edd5946b3cb2a5, 0x3946a50da827f406),
        (0xbc4b4ca91be60c39, 0x3fa8b51f21068ea2),
        (0x3bb25eeb0c76260b, 0xbf2ed935c7aefa31),
        (0xbbf1e976299c9c65, 0xbf707522a5037f2d),
        (0x3b9aa02b2a63445a, 0x3ef8a196061f8bbc),
        (0x3bc503ddc4900b75, 0x3f21874a47e3c1e3),
        (0x3b255f883f2d82a2, 0xbeac10cf34c04f17),
        (0xbb5b1d29b71073a0, 0xbec3fd6c2d4fa2a4),
        (0xbaeb8c54608e9282, 0x3e50906d55522785),
        (0xbac0c29c0a96e3da, 0x3e5c5a1c124dfa08),
        (0xba7aa4d02da2e1f4, 0xbde7f883b31a5f59),
        (0xba7747971035692d, 0xbdeb668cf53028e0),
        (0x3a1b938795c32e6d, 0x3d7775372d05b2b0),
        (0xba0cfec68b9fa162, 0x3d7331d871d67118),
        (0xb9ae3aa9398e5a58, 0xbd0090102489f480),
        (0xb99e9cec5597a630, 0xbcf461c4066d0e4a),
        (0x3914deb7f36c5b05, 0x3c81abf9941a3423),
        (0xb8df6670bd950fcb, 0x3c70f715e31f4f53),
        (0xb88076dcbbb6bbc6, 0xbbfd7ed45d36754d),
        (0x3884befa4b49ff70, 0xbbe6badf59009f89),
        (0xb817484bd8d9cd12, 0x3b73c80adb5c7d1b),
        (0x37e6bcd3e344ffce, 0x3b58d97e797244fc),
        (0x377002e0099ec551, 0xbae61ba441542113),
    ],
    [
        (0xb5e15739920a6042, 0x3957411605b9cd8d),
        (0x3c5b1c9821974147, 0x3fb86e51be0a9153),
        (0xbbeaa494385da6b8, 0xbf465ed1b387e5da),
        (0x3c09cfc1363fac8e, 0xbf9046fc5a218a86),
        (0x3bbf17fc5592840d, 0x3f1dca617fefa913),
        (0xbbe411ddd3f3c7e1, 0x3f4a0300221528a7),
        (0xbb35922b096464c9, 0xbed7c7618906f1e2),
        (0x3b985c37a80f1f32, 0xbef3c838897d0a1e),
        (0xbb1cfd1e42294da1, 0x3e820ede9f9dd7dd),
        (0x3b3d337f0c931d76, 0x3e918a94165592bb),
        (0x3aba2064c2ef620a, 0xbe1ff76205118f09),
        (0x3ab62166ef13ed4e, 0xbe24599dfeef01a1),
        (0xba4104861aca6ff9, 0x3db2803e5998312f),
        (0x3a5eb79000c7671a, 0x3db0a33202feb041),
        (0x39ad6a6495675c7f, 0xbd3e2bff73866f0b),
        (0x39d02f8782470b05, 0xbd34329cf32bc6b1),
        (0xb94c1fab2f5d49d8, 0x3cc2424890e7103d),
        (0x395165f2cbeeb591, 0x3cb2eba6575a58fb),
        (0xb8ee488b45b19302, 0xbc410bda2b2fedfa),
        (0x38c169197f6744f6, 0xbc2c2d96bb6d554b),
        (0xb839df5ffc139231, 0x3bb94a6157935d73),
        (0x38438133f79f8aef, 0x3ba113b61477f5c4),
        (0xb7c400b9f75f31c1, 0xbb2e47f731863fb7),
        (0xb7bf6daddcc30386, 0xbb110e1251b8eda4),
    ],
    [
        (0x3c48caabfef07d2b, 0x3fb829d06fee9266),
        (0xb61c7c2ce51dc943, 0xb97019a7cd7017ec),
        (0xbc46f7f24522358e, 0xbfa8289a526d7785),
        (0xbbbd55f5936024f5, 0x3f2cd680355c9eb6),
        (0x3c19dafafbc6423c, 0x3f7017d70f512861),
        (0x3b94ef9247fbd18a, 0xbef707978e2a0db8),
        (0xbbb84f6dc795e259, 0xbf21247ce15e7385),
        (0x3b3c154e3f91cb1f, 0x3eaa3f6125485ec1),
        (0xbb6485e10f66ec36, 0x3ec38da848401be9),
        (0xbae822e833951880, 0xbe4efe39d8db4cd8),
        (0xbaf0a59f02f9da6b, 0xbe5bbd40f1db0e94),
        (0x3a88e5731b6d2f4f, 0x3de66f7d49436f83),
        (0xba3e79da115ba52a, 0x3dead0e8148229af),
        (0xb999a0303d870988, 0xbd75f78595133d59),
        (0x39d279bb13e85151, 0xbd72ca9bb1031b0b),
        (0x398dc73294489cd5, 0x3cff09e49c89467c),
        (0xb98c493ba49ccabf, 0x3cf3f60eeb57dab2),
        (0xb928107dd8db2eda, 0xbc8091d68a71f568),
        (0x391a43591422376a, 0xbc709f3347db6cf2),
        (0xb88bd7e4b3dfccd4, 0x3bfbad357dd86bf0),
        (0x386a0a018ebdb436, 0x3be647b364d0856e),
        (0x381f08d09391197b, 0xbb72939f68459b07),
        (0xb7ea807173f17da8, 0xbb585e93fb40bea8),
        (0xb78b806a2be014c6, 0x3ae4c7339928aa1f),
    ],
    [
        (0xb5dba30ce3053693, 0x3937ea2409abe46b),
        (0x3c21907f595a082a, 0xbfb7e656efb009ae),
        (0x3bea53822a52cdff, 0x3f44f15066f3d876),
        (0xbbef8f1c97361f96, 0x3f8fd932c26aad94),
        (0xbba4e6bbdef5d491, 0xbf1be460dd86a0a4),
        (0x3bdda15b3ef145a9, 0xbf49733b591879f8),
        (0x3b7b0e61e74c413c, 0x3ed64488c56022e0),
        (0xbb74588b74c5876f, 0x3ef35ba58bf2f993),
        (0xbb2d010d38a866aa, 0xbe80ea47bceb9a8f),
        (0xbb1839485029ea9f, 0xbe912b327055d0e9),
        (0xbab04b71fdf73b8a, 0x3e1df42fc4e2482c),
        (0xba732a2fe61d6156, 0x3e23ec3e76876dba),
        (0xba511e7db88fce7c, 0xbdb1580393f473be),
        (0xba5b5424c61c9156, 0xbdb04b0353c35933),
        (0xb9c9856b07421c8e, 0x3d3c4ca30ffe1933),
        (0x39c738eb74a7f61d, 0x3d33c947b936be26),
        (0xb9400509632244dc, 0xbcc122c73f63e006),
        (0x3916d7eae173a1b3, 0xbcb28ac72945c28e),
        (0xb8e73a0527a29464, 0x3c4002241bb36c01),
        (0x38cad2069d44143f, 0x3c2ba04276d6d12c),
        (0x385b64b543fd96d2, 0xbbb7c48f7753a6e6),
        (0xb848072abf17cd35, 0xbba0bffbbb75a502),
        (0xb7a49b7af751a348, 0x3b2c7ae054639d4d),
        (0x37b54394c52e23e8, 0x3b10bc481d251ae6),
    ],
    [
        (0xbc2f952341a4610c, 0xbfb7a62320798175),
        (0x3608d7aa0df33a13, 0x396cf7a4a97dd69b),
        (0x3c4f5aadff868840, 0x3fa7a50ca4504bb8),
        (0xbbc425ed1e39553c, 0xbf2b095ccb50a68c),
        (0x3c0088e41716f14b, 0xbf6f80ef11daa37a),
        (0xbb8b1b49b7c35e54, 0x3ef59822dc75b064),
        (0x3bc22455d3c1fa15, 0x3f20c7e6a7c66630),
        (0x3b4cb7a0033a7704, 0xbea89e00b5c358d0),
        (0x3b4c5ee49818e61f, 0xbec324d5238b26d0),
        (0x3ae1d17d1c1f9459, 0x3e4d13a888d5dd30),
        (0xbac771ff3649d30b, 0x3e5b29f941a95b07),
        (0xba65c3a8faacb7ca, 0xbde50e6ebb3d1df7),
        (0xba5e1fa2b117348a, 0xbdea4434c6b929d9),
        (0x3a1aa33d3910c73c, 0x3d74a03fe2c47ae3),
        (0x39fc6f3d0c6f0177, 0x3d7269628f2df267),
        (0x3998d07708e4c6bf, 0xbcfd28cc08f5cbae),
        (0xb95557f1f159116a, 0xbcf390706324d249),
        (0x3901d36b70d2e2f4, 0x3c7f26c3f4d4f4a6),
        (0x38f56e280758755e, 0x3c704c1cbbe93d3e),
        (0xb86ceb268d238b88, 0xbbfa089cde4d22c5),
        (0xb889788b18499f22, 0xbbe5da907334a8cc),
        (0x38130abe343bca1e, 0x3b7181ee2ddebf21),
        (0x37f04d8f909b3d76, 0x3b57e93ad9acbcc3),
        (0x378249417d24489f, 0xbae648d6e71bec37),
    ],
];

/**
J1 zeros and extremums.

Generated by Sage:
```python
mp.prec = 200

step = mpf("0.01")
epsilon = mpf("1e-35")
x = mpf("0.0")

def j1_prime(x):
    return diff(lambda t: besselj(1, t), x)

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

# Print results
for i, z in enumerate(j1_zeros):
    print(f"Zero {i+1}: x ≈ {z}")

print("Extrema (peaks/valleys) of J1(x):")
for e in j1_extrema:
    print(f"nExtrema: {e}")

j1_zeros.extend(j1_extrema)

j1_zeros = sorted(j1_zeros)

# Print results
for i, z in enumerate(j1_zeros):
    print(f"Peak or zero {i+1}: x ≈ {z}")

print("")

print("pub(crate) static J1_ZEROS: [DyadicFloat128; 48] = [")
print(f"DyadicFloat128 {
        sign: DyadicSign::Pos,
        exponent: 0,
        mantissa: 0x0u128,
    },")
for z in j1_zeros:
    print_dyadic(z)

print("];")
```
**/
pub(crate) static J1_ZEROS_RATIONAL: [DyadicFloat128; 48] = [
    DyadicFloat128 {
        sign: DyadicSign::Pos,
        exponent: 0,
        mantissa: 0x0u128,
    },
    DyadicFloat128 {
        sign: DyadicSign::Pos,
        exponent: -127,
        mantissa: 0xebabe8ff_6451d02c_2db0419f_b5cccccd_u128,
    },
    DyadicFloat128 {
        sign: DyadicSign::Pos,
        exponent: -126,
        mantissa: 0xf53aabad_7b78453f_d54ac5c9_b53851ec_u128,
    },
    DyadicFloat128 {
        sign: DyadicSign::Pos,
        exponent: -125,
        mantissa: 0xaa9b2de0_1923395c_646a75d7_53970a3d_u128,
    },
    DyadicFloat128 {
        sign: DyadicSign::Pos,
        exponent: -125,
        mantissa: 0xe07faf9d_a3927f26_ec9316de_3ec51eb8_u128,
    },
    DyadicFloat128 {
        sign: DyadicSign::Pos,
        exponent: -124,
        mantissa: 0x8894c078_5c45074e_0f409efa_e71ae148_u128,
    },
    DyadicFloat128 {
        sign: DyadicSign::Pos,
        exponent: -124,
        mantissa: 0xa2c68685_efe14a04_c214a2ac_56c8f5c3_u128,
    },
    DyadicFloat128 {
        sign: DyadicSign::Pos,
        exponent: -124,
        mantissa: 0xbb4bcbcb_f72d6335_09a74bfd_b831eb85_u128,
    },
    DyadicFloat128 {
        sign: DyadicSign::Pos,
        exponent: -124,
        mantissa: 0xd52dd798_872d112b_ce7fd18e_69347ae1_u128,
    },
    DyadicFloat128 {
        sign: DyadicSign::Pos,
        exponent: -124,
        mantissa: 0xedd14250_bd62c52d_9a8cf1ed_61feb852_u128,
    },
    DyadicFloat128 {
        sign: DyadicSign::Pos,
        exponent: -123,
        mantissa: 0x83c3d9b0_2846245a_b2e2eea9_9d166666_u128,
    },
    DyadicFloat128 {
        sign: DyadicSign::Pos,
        exponent: -123,
        mantissa: 0x901fcd12_732937bf_7425352a_c91f5c29_u128,
    },
    DyadicFloat128 {
        sign: DyadicSign::Pos,
        exponent: -123,
        mantissa: 0x9ced473a_0b651dbc_f9fcd02c_2f41eb85_u128,
    },
    DyadicFloat128 {
        sign: DyadicSign::Pos,
        exponent: -123,
        mantissa: 0xa950a125_0d4fd3d5_5fd5f1d8_35e851ec_u128,
    },
    DyadicFloat128 {
        sign: DyadicSign::Pos,
        exponent: -123,
        mantissa: 0xb614a71e_a6c55ee4_0cc50107_70eeb852_u128,
    },
    DyadicFloat128 {
        sign: DyadicSign::Pos,
        exponent: -123,
        mantissa: 0xc27d98ef_70b1a769_44625798_3a3ccccd_u128,
    },
    DyadicFloat128 {
        sign: DyadicSign::Pos,
        exponent: -123,
        mantissa: 0xcf3ab86e_7508311a_2686480d_8820a3d7_u128,
    },
    DyadicFloat128 {
        sign: DyadicSign::Pos,
        exponent: -123,
        mantissa: 0xdba80a21_3a8a6a17_bdc4e527_af21eb85_u128,
    },
    DyadicFloat128 {
        sign: DyadicSign::Pos,
        exponent: -123,
        mantissa: 0xe85fe7a3_8fe65f5e_98e309b4_482ae148_u128,
    },
    DyadicFloat128 {
        sign: DyadicSign::Pos,
        exponent: -123,
        mantissa: 0xf4d0bcfd_d29950ed_24e160a1_a41f5c29_u128,
    },
    DyadicFloat128 {
        sign: DyadicSign::Pos,
        exponent: -122,
        mantissa: 0x80c23b73_595f7e41_c5e32f99_499dc28f_u128,
    },
    DyadicFloat128 {
        sign: DyadicSign::Pos,
        exponent: -122,
        mantissa: 0x86fc175e_aa718d29_fe7db0c2_1ec6147b_u128,
    },
    DyadicFloat128 {
        sign: DyadicSign::Pos,
        exponent: -122,
        mantissa: 0x8d54486e_2f4bdc54_08801f75_f300f5c3_u128,
    },
    DyadicFloat128 {
        sign: DyadicSign::Pos,
        exponent: -122,
        mantissa: 0x938f58dc_0218717d_864bbf17_a30b3333_u128,
    },
    DyadicFloat128 {
        sign: DyadicSign::Pos,
        exponent: -122,
        mantissa: 0x99e6291e_ae5b48cf_57f65065_27e4cccd_u128,
    },
    DyadicFloat128 {
        sign: DyadicSign::Pos,
        exponent: -122,
        mantissa: 0xa0223f28_6d8c1d46_26a79fea_cbe0f5c3_u128,
    },
    DyadicFloat128 {
        sign: DyadicSign::Pos,
        exponent: -122,
        mantissa: 0xa677e78b_9a5b10a4_4deaab72_884e6666_u128,
    },
    DyadicFloat128 {
        sign: DyadicSign::Pos,
        exponent: -122,
        mantissa: 0xacb4de39_388419d0_fd96f29c_211bd70a_u128,
    },
    DyadicFloat128 {
        sign: DyadicSign::Pos,
        exponent: -122,
        mantissa: 0xb3098aeb_5899fa9e_e2de62f3_b20428f6_u128,
    },
    DyadicFloat128 {
        sign: DyadicSign::Pos,
        exponent: -122,
        mantissa: 0xb9474495_447e45d4_e9138245_e303851f_u128,
    },
    DyadicFloat128 {
        sign: DyadicSign::Pos,
        exponent: -122,
        mantissa: 0xbf9b1890_1456ac18_d5eaf7a4_90e1999a_u128,
    },
    DyadicFloat128 {
        sign: DyadicSign::Pos,
        exponent: -122,
        mantissa: 0xc5d97d10_1bef1b2e_cad77847_2dc75c29_u128,
    },
    DyadicFloat128 {
        sign: DyadicSign::Pos,
        exponent: -122,
        mantissa: 0xcc2c947c_b6a8f30b_af7b6617_5e28a3d7_u128,
    },
    DyadicFloat128 {
        sign: DyadicSign::Pos,
        exponent: -122,
        mantissa: 0xd26b8fe5_ab7c5d1c_0f53ea6c_1d04cccd_u128,
    },
    DyadicFloat128 {
        sign: DyadicSign::Pos,
        exponent: -122,
        mantissa: 0xd8be01c6_100c627a_83c08395_4d08a3d7_u128,
    },
    DyadicFloat128 {
        sign: DyadicSign::Pos,
        exponent: -122,
        mantissa: 0xdefd8375_bc854e19_2fb8e927_2df947ae_u128,
    },
    DyadicFloat128 {
        sign: DyadicSign::Pos,
        exponent: -122,
        mantissa: 0xe54f62d4_119259a1_39ce2cd0_8ac4cccd_u128,
    },
    DyadicFloat128 {
        sign: DyadicSign::Pos,
        exponent: -122,
        mantissa: 0xeb8f5cc3_417838c5_ae63da28_0e0dc28f_u128,
    },
    DyadicFloat128 {
        sign: DyadicSign::Pos,
        exponent: -122,
        mantissa: 0xf1e0b98e_b278efbb_46443469_ac7b3333_u128,
    },
    DyadicFloat128 {
        sign: DyadicSign::Pos,
        exponent: -122,
        mantissa: 0xf8211fcc_da5a984e_66ff290c_00880000_u128,
    },
    DyadicFloat128 {
        sign: DyadicSign::Pos,
        exponent: -122,
        mantissa: 0xfe72077d_88ab71fd_1ee82863_58a1999a_u128,
    },
    DyadicFloat128 {
        sign: DyadicSign::Pos,
        exponent: -121,
        mantissa: 0x825967e5_d8426804_01982d45_40bc51ec_u128,
    },
    DyadicFloat128 {
        sign: DyadicSign::Pos,
        exponent: -121,
        mantissa: 0x8581a6ef_3ada11d3_cacfc720_4188cccd_u128,
    },
    DyadicFloat128 {
        sign: DyadicSign::Pos,
        exponent: -121,
        mantissa: 0x88a237b0_78a2c052_e1df53db_64035c29_u128,
    },
    DyadicFloat128 {
        sign: DyadicSign::Pos,
        exponent: -121,
        mantissa: 0x8bca46db_1b3ade19_96cfb401_8b840000_u128,
    },
    DyadicFloat128 {
        sign: DyadicSign::Pos,
        exponent: -121,
        mantissa: 0x8eeb005b_a1d4db5a_4194eb14_e29ee148_u128,
    },
    DyadicFloat128 {
        sign: DyadicSign::Pos,
        exponent: -121,
        mantissa: 0x9212e3ee_566fb31d_695296c7_26c35c29_u128,
    },
    DyadicFloat128 {
        sign: DyadicSign::Pos,
        exponent: -121,
        mantissa: 0x9533c2cd_e320ecbc_d2eda12f_8f5d999a_u128,
    },
];

/**
Taylor coefficient at zero/extremum for J1.

Generated by SageMath:
```python
def print_taylor_coeffs_dyad(poly):
    print("[")
    for i in range(0, 24):
        coeff = poly[i]
        print_dyadic(coeff)
    print("],")

print(f"pub(crate) static J1_COEFFS_RATIONAL128: [DyadicFloat128; {len(j1_zeros)}] = [")

for i in range(0, len(j1_zeros)):
    k_range = j1_zeros[i]
    range_diff = k_range - prev_zero
    g_c = 1

    x0 = mp.mpf(k_range)
    from mpmath import mp, j1, taylor
    poly = taylor(lambda val: j1(val), x0, 25)
    print_taylor_coeffs_dyad(poly)
    prev_zero = j1_zeros[i]

print("];")
```
**/
pub(crate) static J1_COEFFS_RATIONAL128: [[DyadicFloat128; 24]; 47] = [
    [
        DyadicFloat128 {
            sign: DyadicSign::Pos,
            exponent: -128,
            mantissa: 0x94f51e8c_f81af88f_8e192229_d983c515_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Pos,
            exponent: -236,
            mantissa: 0x818dba08_700886d0_9cadc65b_432b1ef3_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Neg,
            exponent: -130,
            mantissa: 0xd2088ae2_ef921bc6_c635c2e2_4e08076c_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Pos,
            exponent: -135,
            mantissa: 0xc68a2473_7f6a3cea_675dfe29_e99fd66b_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Pos,
            exponent: -134,
            mantissa: 0xe220d17c_ef1158f4_2d3e85e0_4c781388_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Neg,
            exponent: -139,
            mantissa: 0xc3338e0c_58444cd6_6bc92ddb_6daa355f_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Neg,
            exponent: -139,
            mantissa: 0xcf12826e_e3b03cd5_e6762f6f_a5d1d805_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Pos,
            exponent: -144,
            mantissa: 0x9a665e50_63ae8772_5bf5e38d_fbb9c0a2_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Pos,
            exponent: -145,
            mantissa: 0xd24b9bc2_68843b57_40bc50fa_2df5de9f_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Neg,
            exponent: -150,
            mantissa: 0x88229915_555a2600_53514344_8d50d8c7_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Neg,
            exponent: -151,
            mantissa: 0x87d706d3_66e776c2_54517015_8c230d3d_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Pos,
            exponent: -157,
            mantissa: 0x9aa3677a_f680538e_3370d255_4e805daa_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Pos,
            exponent: -159,
            mantissa: 0xf2f74173_33b83e53_8a283a8f_faa9897a_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Neg,
            exponent: -165,
            mantissa: 0xf6406645_b1e1c963_b1d65111_e4c5a7d3_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Neg,
            exponent: -166,
            mantissa: 0x9f5974cb_14da2145_e6e27b14_90ac0242_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Pos,
            exponent: -172,
            mantissa: 0x915fe6d8_cab0e176_26e428b4_f5725260_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Pos,
            exponent: -174,
            mantissa: 0x9fd99bb5_2171a75e_f0fd978f_93a2a23b_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Neg,
            exponent: -180,
            mantissa: 0x8480943e_882d7097_12d63478_6406f825_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Neg,
            exponent: -183,
            mantissa: 0xfd354863_50619475_e812615c_e44da777_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Pos,
            exponent: -189,
            mantissa: 0xc03ba78b_e2435580_05b08a99_10eae691_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Pos,
            exponent: -191,
            mantissa: 0xa25bcabc_7e37daed_e6c4d840_525b5fb3_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Neg,
            exponent: -198,
            mantissa: 0xe35582ff_3b983b88_33390ef5_d5354d17_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Neg,
            exponent: -200,
            mantissa: 0xac02fde9_c710bc24_88e98162_11816735_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Pos,
            exponent: -207,
            mantissa: 0xdf6ade75_118d117f_411d416c_35b41dd9_u128,
        },
    ],
    [
        DyadicFloat128 {
            sign: DyadicSign::Neg,
            exponent: -241,
            mantissa: 0xcd48fcb9_f6181615_46a7547a_e4a9a4bc_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Neg,
            exponent: -129,
            mantissa: 0xce367ac1_65fbf6d2_1eebc89a_569f0b2b_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Pos,
            exponent: -132,
            mantissa: 0xd7451cfa_8d681d28_e5199aca_35ab411b_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Pos,
            exponent: -132,
            mantissa: 0xdac4e8ed_09c82511_304d2ce9_59432cf1_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Neg,
            exponent: -135,
            mantissa: 0xa9baa261_98ed3619_a1be074e_a9115a31_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Neg,
            exponent: -136,
            mantissa: 0x9259a04c_ac831fdd_d05d4df1_8fa9f556_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Pos,
            exponent: -140,
            mantissa: 0xb72616a9_aa11218c_81b03584_e4348fd2_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Pos,
            exponent: -142,
            mantissa: 0xc1d03718_62084985_9e8c0ae3_7ce6ab19_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Neg,
            exponent: -146,
            mantissa: 0xcbccea64_f92a46e0_19711b20_1d9dc2c2_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Neg,
            exponent: -148,
            mantissa: 0x99c12e69_7160b3ff_623fb4b7_de637daa_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Pos,
            exponent: -152,
            mantissa: 0x8b834919_f48b60c5_2e76253b_e35fa06b_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Pos,
            exponent: -155,
            mantissa: 0xa2b4d915_7e1e3f90_91258fda_ffa22bb9_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Neg,
            exponent: -159,
            mantissa: 0x81dcf4b2_882b4f4e_e4ad6f1f_099e7108_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Neg,
            exponent: -163,
            mantissa: 0xf6311885_7afa930a_4c0c34e6_33b45117_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Pos,
            exponent: -167,
            mantissa: 0xaf642723_db7a78c6_0409e80d_652b1dec_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Pos,
            exponent: -170,
            mantissa: 0x8bd2064a_1170b9f1_57bfc6ea_f253ba40_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Neg,
            exponent: -175,
            mantissa: 0xb3e58e80_1864dc48_1bcb2452_a3b3abeb_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Neg,
            exponent: -179,
            mantissa: 0xf73ffc8c_7cdb4beb_ad2e9462_17660328_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Pos,
            exponent: -183,
            mantissa: 0x90fdc5ee_26bfcf8b_fea2fef9_e563aef6_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Pos,
            exponent: -187,
            mantissa: 0xaf028bee_85ae4846_8d5834e7_a6499c65_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Neg,
            exponent: -192,
            mantissa: 0xbc92cad0_f7633dd5_575c6f58_00bb6bf5_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Neg,
            exponent: -196,
            mantissa: 0xcad34101_16e5c0cd_42a81f5a_5d74eba8_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Pos,
            exponent: -201,
            mantissa: 0xca1c6232_505ccac3_4d4f17c3_cc44bbd3_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Pos,
            exponent: -205,
            mantissa: 0xc3fae4a1_c12a1b49_353c12a2_62e6c6cb_u128,
        },
    ],
    [
        DyadicFloat128 {
            sign: DyadicSign::Neg,
            exponent: -129,
            mantissa: 0xb137741a_805f9305_ac4e2c0e_29393269_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Pos,
            exponent: -238,
            mantissa: 0xcc43e0ad_d96c321e_0816785d_5fe79885_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Pos,
            exponent: -130,
            mantissa: 0xaafb5f64_f7cb135d_12a5cde3_b31b622a_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Neg,
            exponent: -134,
            mantissa: 0x9e9199b7_e8871f16_7872e7cd_6a84672b_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Neg,
            exponent: -134,
            mantissa: 0xc63bd4c1_d0409d2d_4794ec64_6042b080_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Pos,
            exponent: -138,
            mantissa: 0xae6e4c6d_8e5f1342_401e710e_66320508_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Pos,
            exponent: -139,
            mantissa: 0xb9abb7fa_3771da68_64919d4a_5786c293_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Neg,
            exponent: -143,
            mantissa: 0x9230a23e_bda11657_32941247_7e49a2c3_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Neg,
            exponent: -145,
            mantissa: 0xbdc29a2b_5b754e57_56fd3781_94566688_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Pos,
            exponent: -149,
            mantissa: 0x855fe341_3a54bfc1_1b593185_00c0138a_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Pos,
            exponent: -152,
            mantissa: 0xf53d0f71_30926730_ee8e294c_8960d20d_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Neg,
            exponent: -156,
            mantissa: 0x9ae0209f_00a0c2ad_30fbbc4d_6c22811e_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Neg,
            exponent: -159,
            mantissa: 0xdae2ea89_7dd7fc29_2f317584_5c8c1d24_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Pos,
            exponent: -164,
            mantissa: 0xfa62f137_eb7481a4_481d08d2_198e3808_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Pos,
            exponent: -166,
            mantissa: 0x8f26219c_bdf285c3_9612f5ab_a93553fe_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Neg,
            exponent: -171,
            mantissa: 0x956ef7e2_6a8e25d6_793f6f5d_1e6c2184_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Neg,
            exponent: -174,
            mantissa: 0x8f27d6e0_721ae020_ba0bff98_dbd46d47_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Pos,
            exponent: -179,
            mantissa: 0x89505c42_8479b7b1_3eb99f7e_1f70fa25_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Pos,
            exponent: -183,
            mantissa: 0xe2123501_a171cb4b_584571e2_9799f0c9_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Neg,
            exponent: -188,
            mantissa: 0xc876f2c0_b79673cf_144f4b55_7e64c3d1_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Neg,
            exponent: -191,
            mantissa: 0x908983c8_3718a30d_a469a8d1_b6f770b2_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Pos,
            exponent: -197,
            mantissa: 0xee3dfb81_cb5940de_23344b87_3ede7a57_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Pos,
            exponent: -200,
            mantissa: 0x98b73f8a_c5132329_3f48dedb_aece19d6_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Neg,
            exponent: -206,
            mantissa: 0xeb11f6e1_d9256116_7d269458_30c7f691_u128,
        },
    ],
    [
        DyadicFloat128 {
            sign: DyadicSign::Neg,
            exponent: -238,
            mantissa: 0xfc5058ca_74954dc7_bfc4329d_12c83fad_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Pos,
            exponent: -129,
            mantissa: 0x99a8c59c_3a74535e_45a067dc_1490486d_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Neg,
            exponent: -133,
            mantissa: 0xaf386e30_1b15f93b_4d629847_2ae2b445_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Neg,
            exponent: -132,
            mantissa: 0xc0641def_772d82ed_b0e2bbed_0594a2d2_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Pos,
            exponent: -136,
            mantissa: 0xcd259497_1ef2127b_58415645_a9f8ffd7_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Pos,
            exponent: -136,
            mantissa: 0x89fde3eb_4c10b862_1fbb5b9a_20283175_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Neg,
            exponent: -140,
            mantissa: 0x839ac5de_dffc88e2_95492c88_73a85610_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Neg,
            exponent: -142,
            mantissa: 0xbcb53a8f_c4828828_ec7a04da_64c5861c_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Pos,
            exponent: -146,
            mantissa: 0xa12ad80a_d6f69c41_341c0456_e210e50d_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Pos,
            exponent: -148,
            mantissa: 0x98137067_4bd5c651_7e92d987_b02fd8cf_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Neg,
            exponent: -153,
            mantissa: 0xea46e6d7_497962af_67107731_f45c66ec_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Neg,
            exponent: -155,
            mantissa: 0xa231cebf_75852e5b_3af4a11f_d052017a_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Pos,
            exponent: -160,
            mantissa: 0xe31665a5_1975a404_e1968c26_a0819166_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Pos,
            exponent: -163,
            mantissa: 0xf6574974_2a22b5f9_a68f1d1c_138f6fbb_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Neg,
            exponent: -167,
            mantissa: 0x9ddb44c6_b9e6937d_63dafd5e_2bf6fc49_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Neg,
            exponent: -170,
            mantissa: 0x8c1f6ddc_77aaf1d5_420c79c5_0f5580a5_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Pos,
            exponent: -175,
            mantissa: 0xa571f3f2_6671eb4c_d620b745_ac979f34_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Pos,
            exponent: -179,
            mantissa: 0xf7dd8f38_03984fce_68da2b35_c3029032_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Neg,
            exponent: -183,
            mantissa: 0x87946d1f_7d1e3a8d_85e56ce1_a10738bd_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Neg,
            exponent: -187,
            mantissa: 0xaf618c6f_cda336b5_3ca47df9_2b5c240a_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Pos,
            exponent: -192,
            mantissa: 0xb2aad4a9_1a2fd37c_3227d051_8bc34b7c_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Pos,
            exponent: -196,
            mantissa: 0xcb19e3f8_99e493ff_331fd37e_6f7122a4_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Neg,
            exponent: -201,
            mantissa: 0xc18a2fc9_e3856cc0_346a5a5b_97cb6998_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Neg,
            exponent: -205,
            mantissa: 0xc40e5ccf_350ec490_bd55fbf4_e9ce6798_u128,
        },
    ],
    [
        DyadicFloat128 {
            sign: DyadicSign::Pos,
            exponent: -129,
            mantissa: 0x8bedf84e_a0692456_e480be2b_cb772a02_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Pos,
            exponent: -237,
            mantissa: 0x84c9625c_3d163144_27f709d1_f02da29d_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Neg,
            exponent: -130,
            mantissa: 0x8a025fb2_3e147a3d_7229a678_5f3b05d1_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Pos,
            exponent: -135,
            mantissa: 0xa7a6f93b_4fc18198_a23a3934_c8097716_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Pos,
            exponent: -134,
            mantissa: 0xae3142a1_4db36a95_d27b975c_85b206ea_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Neg,
            exponent: -139,
            mantissa: 0xeb455b91_440de444_9b4c52dc_ca4c4729_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Neg,
            exponent: -139,
            mantissa: 0xab565b22_96c30319_ae63ed7e_7feda252_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Pos,
            exponent: -144,
            mantissa: 0xe085a3e7_bf7b4886_95c7fadc_7a3b688f_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Pos,
            exponent: -145,
            mantissa: 0xb3f5574b_dde42fc6_702e4f9e_b632bd3a_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Neg,
            exponent: -150,
            mantissa: 0xddb29863_1f96f76e_78878cbf_e02dfa22_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Neg,
            exponent: -152,
            mantissa: 0xec3900f2_2876e983_f73085e8_89347958_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Pos,
            exponent: -156,
            mantissa: 0x87a3e41f_62a85708_0ede218c_59a2102a_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Pos,
            exponent: -159,
            mantissa: 0xd4c198fb_753ccb2b_20bb49f3_8ffb47f5_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Neg,
            exponent: -164,
            mantissa: 0xe3820a51_1b5335e9_7f8c19b1_7816b224_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Neg,
            exponent: -166,
            mantissa: 0x8be02bd2_fe69ea93_aab5dc96_21d6772d_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Pos,
            exponent: -171,
            mantissa: 0x8b7f50a6_8ff3e204_2c968df8_23e67d02_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Pos,
            exponent: -174,
            mantissa: 0x8c4de9f6_1b11a190_498a42b8_b9e17e5b_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Neg,
            exponent: -179,
            mantissa: 0x82d7a255_79559dda_05c46bf3_7cd3d738_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Neg,
            exponent: -183,
            mantissa: 0xddeda3d7_35ca670b_9408acc5_0ee9a842_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Pos,
            exponent: -188,
            mantissa: 0xc21af423_0f69b6b3_3a4e09c3_cdfa7763_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Pos,
            exponent: -191,
            mantissa: 0x8dfeb37d_d3d9b12c_f7a81930_22873bd0_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Neg,
            exponent: -197,
            mantissa: 0xe9a937d1_632f118c_ec539879_f2188377_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Neg,
            exponent: -200,
            mantissa: 0x960f159f_9046696b_be8d3df2_73084793_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Pos,
            exponent: -206,
            mantissa: 0xe8f897dc_d74be83d_91869579_e0ce2074_u128,
        },
    ],
    [
        DyadicFloat128 {
            sign: DyadicSign::Pos,
            exponent: -239,
            mantissa: 0xc66f76eb_844ecba5_9714a145_3cf4f092_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Neg,
            exponent: -130,
            mantissa: 0xffb2a2a2_75e68905_2a3a2541_c56b22da_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Pos,
            exponent: -134,
            mantissa: 0xc911ff96_03c2d772_5fe95e8a_d7bb699f_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Pos,
            exponent: -132,
            mantissa: 0xa5862eae_d33c46dc_5d8aa6f5_7a146e74_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Neg,
            exponent: -137,
            mantissa: 0xfc8d4f70_6944b7d8_2b56a72e_20fc2207_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Neg,
            exponent: -137,
            mantissa: 0xfa8e1244_dcf37477_c1013d2f_a7d881db_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Pos,
            exponent: -141,
            mantissa: 0xb5a64e50_7bb869d7_af9f28c7_680506cf_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Pos,
            exponent: -142,
            mantissa: 0xb1e2a3aa_1ce58c45_97d4bf4e_50973ce2_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Neg,
            exponent: -147,
            mantissa: 0xf1b92ed7_b4c334c1_27099933_231e6b9f_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Neg,
            exponent: -148,
            mantissa: 0x92e091c5_9972c770_6ac12352_aee7d434_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Pos,
            exponent: -153,
            mantissa: 0xba437b5c_d54a864c_a183a6c8_f9f6295a_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Pos,
            exponent: -155,
            mantissa: 0x9f1df924_13bf732b_87c32fbd_60a4b790_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Neg,
            exponent: -160,
            mantissa: 0xbc51c539_f3e0517b_9fdb9bee_257dc8f1_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Neg,
            exponent: -163,
            mantissa: 0xf42275b5_908d8180_ea277e7e_c500bc7f_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Pos,
            exponent: -167,
            mantissa: 0x871255da_038471c6_dae4342e_4f756866_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Pos,
            exponent: -170,
            mantissa: 0x8bcbf2f5_66efd540_d9d23c74_20b4d993_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Neg,
            exponent: -175,
            mantissa: 0x90fe0789_32d3dec7_82c86cb1_20f1c295_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Neg,
            exponent: -179,
            mantissa: 0xf8609722_2e5a0b71_55ab5d8e_855617d5_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Pos,
            exponent: -184,
            mantissa: 0xf2211aac_840593c8_5eb6aba7_04950881_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Pos,
            exponent: -187,
            mantissa: 0xb040231d_a0f8e9c1_d5a56b93_a77404d1_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Neg,
            exponent: -192,
            mantissa: 0xa1f1c593_1839f1db_c48deb40_e7d5a929_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Neg,
            exponent: -196,
            mantissa: 0xcc7d090e_40ae557a_b7ef33b2_9e64a687_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Pos,
            exponent: -201,
            mantissa: 0xb192e5d0_b7d31194_a11e9993_a13dc031_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Pos,
            exponent: -205,
            mantissa: 0xc59fea1c_6355987d_981d9deb_5b2d671a_u128,
        },
    ],
    [
        DyadicFloat128 {
            sign: DyadicSign::Neg,
            exponent: -130,
            mantissa: 0xeee75a67_0dfa4ef0_76bd2c06_fd55227d_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Pos,
            exponent: -238,
            mantissa: 0x84fb9283_30d5f7a2_77eadea0_911d00fd_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Pos,
            exponent: -131,
            mantissa: 0xed2908b6_05320133_d5334b42_c830a721_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Neg,
            exponent: -136,
            mantissa: 0xd4ed2301_db3f50b9_2345144a_17b447a1_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Neg,
            exponent: -134,
            mantissa: 0x998f3a75_2cd5c0e2_f0418be2_de042808_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Pos,
            exponent: -139,
            mantissa: 0x9f2e5b75_d3754ca7_750c8905_9c0fba59_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Pos,
            exponent: -139,
            mantissa: 0x9c42ff4d_7d2a0a15_b2e7bfda_f5126f47_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Neg,
            exponent: -144,
            mantissa: 0xa4a607a5_8340594d_72fe9613_abab0c74_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Neg,
            exponent: -145,
            mantissa: 0xa895ce9b_bb16b4f8_31e00806_e1016958_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Pos,
            exponent: -150,
            mantissa: 0xad430841_5fdbf937_aef41758_71fa2d1e_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Pos,
            exponent: -152,
            mantissa: 0xe191f505_0216154c_063ab178_1994038c_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Neg,
            exponent: -157,
            mantissa: 0xde64b17b_dc8f59eb_25b7edb5_f5c8682c_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Neg,
            exponent: -159,
            mantissa: 0xcde4a717_94d28f96_8d43d6e9_06e730da_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Pos,
            exponent: -164,
            mantissa: 0xc15e37e7_d4744a22_bef51c3b_46e51fa6_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Pos,
            exponent: -166,
            mantissa: 0x88a0e73d_bc230394_5d34f144_618144af_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Neg,
            exponent: -172,
            mantissa: 0xf3ce658c_306bcfa6_db9dede3_ea034323_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Neg,
            exponent: -174,
            mantissa: 0x89f0fdf0_52cb8786_8227f576_ab334875_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Pos,
            exponent: -180,
            mantissa: 0xe9b68958_6924e79e_f0c6c5c8_b9292406_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Pos,
            exponent: -183,
            mantissa: 0xdb30395d_efa3bcdd_3ac1db8d_ce90e5fd_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Neg,
            exponent: -188,
            mantissa: 0xb065c91d_1aa19408_95644e59_2ab1cb60_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Neg,
            exponent: -191,
            mantissa: 0x8cb11ba8_05926fdb_165be803_2a7a2efa_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Pos,
            exponent: -197,
            mantissa: 0xd75d0802_f0c182ab_29076509_0905c610_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Pos,
            exponent: -200,
            mantissa: 0x9503e8b2_63537b49_770b7196_57312fe4_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Neg,
            exponent: -206,
            mantissa: 0xd93d6097_82518839_8c0ff46c_0bc28c58_u128,
        },
    ],
    [
        DyadicFloat128 {
            sign: DyadicSign::Pos,
            exponent: -243,
            mantissa: 0xc570a566_01e6110e_66d0c71e_44c1677f_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Pos,
            exponent: -130,
            mantissa: 0xdf999bc3_9d3ec391_8cda5c84_0c47949c_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Neg,
            exponent: -134,
            mantissa: 0x8641d16b_d6e99542_fc768aaf_62eeb205_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Neg,
            exponent: -132,
            mantissa: 0x928c2c00_8c0b557c_db2da506_c4a06645_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Pos,
            exponent: -137,
            mantissa: 0xacf58b05_fb96c082_3ae248b5_3f8e8b99_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Pos,
            exponent: -137,
            mantissa: 0xe2de719d_796bb80f_026a0df8_6c3dae92_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Neg,
            exponent: -141,
            mantissa: 0x8209f183_7001c7c2_96de79af_b7107286_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Neg,
            exponent: -142,
            mantissa: 0xa5338268_2d6859e3_a426e412_d41c09b7_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Pos,
            exponent: -147,
            mantissa: 0xb621f76f_f6b644b6_199bfb2f_c9b9ae14_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Pos,
            exponent: -148,
            mantissa: 0x8b55ebc0_aef3a0e9_788869c6_91d938c1_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Neg,
            exponent: -153,
            mantissa: 0x92bf8b7a_ea1a3141_6eabc41c_a78d06bf_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Neg,
            exponent: -155,
            mantissa: 0x996d8da5_97fc5a20_86e8956c_25f60990_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Pos,
            exponent: -160,
            mantissa: 0x99d6667b_dfee6d1a_6b9c39c5_eb2a865b_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Pos,
            exponent: -163,
            mantissa: 0xee47ab41_2b35d010_b03336c2_218ddab3_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Neg,
            exponent: -168,
            mantissa: 0xe32899c3_6d5710ce_68ad6836_9c598760_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Neg,
            exponent: -170,
            mantissa: 0x89ac2d4d_ba39f67b_9952ad10_59ef2626_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Pos,
            exponent: -176,
            mantissa: 0xf9917517_106dafde_e6dea9bd_9ac29c73_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Pos,
            exponent: -179,
            mantissa: 0xf63a4dd1_8e876a5c_0ff46304_d651a7e0_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Neg,
            exponent: -184,
            mantissa: 0xd455d5a2_52d0ad6b_78a383ca_4128d7d1_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Neg,
            exponent: -187,
            mantissa: 0xaf955d13_e80aaf14_546c6c8b_01513535_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Pos,
            exponent: -192,
            mantissa: 0x903828b4_d04426ac_33193536_51bbff69_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Pos,
            exponent: -196,
            mantissa: 0xcc74e4fd_083252b1_0497684d_dc8f961d_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Neg,
            exponent: -201,
            mantissa: 0xa02abe44_eb2e496e_8cdd35ee_ba72f66e_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Neg,
            exponent: -205,
            mantissa: 0xc620e6f6_d869eff3_05ee263d_8635c9c5_u128,
        },
    ],
    [
        DyadicFloat128 {
            sign: DyadicSign::Pos,
            exponent: -130,
            mantissa: 0xd3fb1ff5_40f92fd3_8d1ed56f_088f447a_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Pos,
            exponent: -238,
            mantissa: 0x865f27f0_aefc0799_77082eda_fd5edd6e_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Neg,
            exponent: -131,
            mantissa: 0xd3057d83_320676cb_e26d1450_f1ac7c4b_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Pos,
            exponent: -136,
            mantissa: 0x960f4984_9ae9e354_19efdbf8_95e123e9_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Pos,
            exponent: -134,
            mantissa: 0x8a283233_6bf8fa4a_62eb6bb8_516adac9_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Neg,
            exponent: -140,
            mantissa: 0xe6546060_7551f414_d61269c3_f70d8302_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Neg,
            exponent: -139,
            mantissa: 0x8efc10e6_09bbf123_306ba391_f4394eef_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Pos,
            exponent: -145,
            mantissa: 0xf7440a76_8562260b_d8c14e72_b7a0abec_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Pos,
            exponent: -145,
            mantissa: 0x9d1b2d20_ccee882a_63c01816_3d227c34_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Neg,
            exponent: -150,
            mantissa: 0x87697ce1_f22c72a7_35d67a92_8d4e3fe1_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Neg,
            exponent: -152,
            mantissa: 0xd59d9be2_938d98be_e3e20c00_54c13ba4_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Pos,
            exponent: -157,
            mantissa: 0xb426b731_5ae32c39_efe97e28_6e3062cb_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Pos,
            exponent: -159,
            mantissa: 0xc5882d28_90768cfe_02e4a83e_793f28d1_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Neg,
            exponent: -164,
            mantissa: 0xa16e2cc8_dce8094a_fb054573_dded9a40_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Neg,
            exponent: -166,
            mantissa: 0x846b202f_ff1edd63_d34dd01e_f5d13ffa_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Pos,
            exponent: -172,
            mantissa: 0xd0a9019b_2d4efe15_c93fa626_599165ae_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Pos,
            exponent: -174,
            mantissa: 0x86be08c5_df59b85c_5d1dabb2_55047408_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Neg,
            exponent: -180,
            mantissa: 0xcc27fe94_2bda6516_8b00c64a_4c8bd7d5_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Neg,
            exponent: -183,
            mantissa: 0xd7659699_1349f884_86f8ab67_c76d852f_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Pos,
            exponent: -188,
            mantissa: 0x9cb7f969_23020ed6_074b2540_12812150_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Pos,
            exponent: -191,
            mantissa: 0x8ae5d722_0e511fb3_bcad3b67_14adf00d_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Neg,
            exponent: -197,
            mantissa: 0xc20df022_7d146b11_6d8c4715_7b32b7d0_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Neg,
            exponent: -200,
            mantissa: 0x93a33dbc_ce67e160_da5e67d9_192ceb8d_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Pos,
            exponent: -206,
            mantissa: 0xc6142bae_ba2d3495_2c871174_4d8fc4ec_u128,
        },
    ],
    [
        DyadicFloat128 {
            sign: DyadicSign::Neg,
            exponent: -238,
            mantissa: 0xcddbc882_c83e4210_818c8dce_6c7c2612_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Neg,
            exponent: -130,
            mantissa: 0xc92e37e5_047aa42d_550665d1_f82d50f7_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Pos,
            exponent: -135,
            mantissa: 0xc36e9970_2cd87191_bef69239_aefd8b88_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Pos,
            exponent: -132,
            mantissa: 0x84a31dde_81b3f9da_c1b118bb_9463ddb5_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Neg,
            exponent: -138,
            mantissa: 0xfed014c6_43b455b9_c01cdbb5_73c25cb5_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Neg,
            exponent: -137,
            mantissa: 0xcfa5f303_ac7d8a03_0aa4fd24_11c4974d_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Pos,
            exponent: -142,
            mantissa: 0xc3bcc8d7_ce8dd477_17f3973e_ff919fb6_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Pos,
            exponent: -142,
            mantissa: 0x99658077_460f98ab_197073a8_81b20c17_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Neg,
            exponent: -147,
            mantissa: 0x8cec6746_1afabe95_fc8e0a14_a226ce01_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Neg,
            exponent: -148,
            mantissa: 0x8350217d_dd22a578_ce24659c_9743aa63_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Pos,
            exponent: -154,
            mantissa: 0xe9d344f3_3bb984a8_9ed07056_aa99845f_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Pos,
            exponent: -155,
            mantissa: 0x92884626_715b1bf0_a3d08cc5_0e7943f3_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Neg,
            exponent: -161,
            mantissa: 0xfbdc74d5_a98a8dfe_72ad6ad3_71e24ff6_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Neg,
            exponent: -163,
            mantissa: 0xe620682b_29321004_e4e60502_216a147a_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Pos,
            exponent: -168,
            mantissa: 0xbe6bb716_bf165785_08b1c45f_0d0cbe24_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Pos,
            exponent: -170,
            mantissa: 0x862c3b81_18e04ebf_9f471ba9_13807e03_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Neg,
            exponent: -176,
            mantissa: 0xd57f6266_c90032df_b9aca45f_427c5c7b_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Neg,
            exponent: -179,
            mantissa: 0xf1b6eab8_446214a7_2b59fda1_9462b4ab_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Pos,
            exponent: -184,
            mantissa: 0xb8ccf405_dc7614a6_03aecf9b_7d448e97_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Pos,
            exponent: -187,
            mantissa: 0xad5c5c34_a8dd8d34_a58b4b29_db003b2d_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Neg,
            exponent: -193,
            mantissa: 0xfeb9bbdd_095f13a2_dfa5bda1_f337d54d_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Neg,
            exponent: -196,
            mantissa: 0xcacad8bd_7d3c69aa_7858e0cb_b983e197_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Pos,
            exponent: -201,
            mantissa: 0x8f367720_cd9a9541_97e383b4_81fc014c_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Pos,
            exponent: -205,
            mantissa: 0xc53a7c3e_813d1ebb_3bf34e5e_1db85a3e_u128,
        },
    ],
    [
        DyadicFloat128 {
            sign: DyadicSign::Neg,
            exponent: -130,
            mantissa: 0xc087a811_2d8257cb_41cd8fff_81467092_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Neg,
            exponent: -239,
            mantissa: 0xeab62441_c78a8e5b_1b932dab_73ab7aa1_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Pos,
            exponent: -131,
            mantissa: 0xbfefcbd6_1b58fa8c_5790d8a1_e14c4059_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Neg,
            exponent: -137,
            mantissa: 0xe1e12b54_65502651_04f2c2f9_3629be09_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Neg,
            exponent: -135,
            mantissa: 0xfcc7f5b9_435a3b23_b7b7c994_db2a6b6a_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Pos,
            exponent: -140,
            mantissa: 0xafb2acf2_b4370c5e_72961278_a506e519_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Pos,
            exponent: -139,
            mantissa: 0x8407abd6_10ad74c4_36b37460_d0f72a1c_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Neg,
            exponent: -145,
            mantissa: 0xc06289cb_f2f5d19e_d06d60ef_e89449d8_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Neg,
            exponent: -145,
            mantissa: 0x92b6daa1_e689fffe_fb14942e_0db56d2b_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Pos,
            exponent: -151,
            mantissa: 0xd7bacc51_0cc1235e_5b636109_7c953b6a_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Pos,
            exponent: -152,
            mantissa: 0xc9cc1136_51182ad2_753f471b_f40b3833_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Neg,
            exponent: -157,
            mantissa: 0x9304c2ec_92c3e750_54c4a177_99540bd4_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Neg,
            exponent: -159,
            mantissa: 0xbc95d9f7_537b4e2a_9ed41f4b_4facb103_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Pos,
            exponent: -164,
            mantissa: 0x86c3134a_f0d3d286_5b2bb2d9_7b4cecad_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Pos,
            exponent: -167,
            mantissa: 0xff2956e1_5d56782f_9a203b33_a9c94128_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Neg,
            exponent: -172,
            mantissa: 0xb1c009bb_9cc77e55_aba5e5a1_49d6315d_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Neg,
            exponent: -174,
            mantissa: 0x82d0df7e_df216000_d22c7870_f1341086_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Pos,
            exponent: -180,
            mantissa: 0xb102820b_9cd52ae7_eaac1f16_89bde6f6_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Pos,
            exponent: -183,
            mantissa: 0xd270a81d_b68720a9_e53947f4_880471c5_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Neg,
            exponent: -188,
            mantissa: 0x89f78a0f_b6b2076e_a361efa6_364aa0b6_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Neg,
            exponent: -191,
            mantissa: 0x886664b6_6a600fd7_b1688b7e_93ba8a99_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Pos,
            exponent: -197,
            mantissa: 0xad17f011_d7513122_5ebf160c_1433dbc7_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Pos,
            exponent: -200,
            mantissa: 0x9196be48_699de273_1fea9f50_846132e1_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Neg,
            exponent: -206,
            mantissa: 0xb2b21ab4_bf8afc61_de4ff2e2_58ddbfa0_u128,
        },
    ],
    [
        DyadicFloat128 {
            sign: DyadicSign::Pos,
            exponent: -239,
            mantissa: 0xded6e926_79c65a18_d3a3ab5b_320501b3_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Pos,
            exponent: -130,
            mantissa: 0xb8628891_3ead525b_400aff09_a7728431_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Neg,
            exponent: -135,
            mantissa: 0x9665874b_aac6d0da_ba728367_badb49ad_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Neg,
            exponent: -133,
            mantissa: 0xf3ee0473_874cce19_ed882120_8fad584a_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Pos,
            exponent: -138,
            mantissa: 0xc566e2d8_2c60702e_f2a5b485_8dd69ad9_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Pos,
            exponent: -137,
            mantissa: 0xc0281b92_56984d97_e1640cd8_c956e853_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Neg,
            exponent: -142,
            mantissa: 0x99772650_fe57da22_bed2e8b2_936413c0_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Neg,
            exponent: -142,
            mantissa: 0x8f2e941b_646cc52b_32437612_1aaf0cd4_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Pos,
            exponent: -148,
            mantissa: 0xe094f83b_dd8b1726_cdefc041_85357607_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Pos,
            exponent: -149,
            mantissa: 0xf78b0ac8_8c0d1270_3e77b1f6_96b8799f_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Neg,
            exponent: -154,
            mantissa: 0xbdcdd83f_8cfc0d73_24fceb74_979e01f8_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Neg,
            exponent: -155,
            mantissa: 0x8b79c9ba_caecaed2_6ad57b57_0df1fc52_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Pos,
            exponent: -161,
            mantissa: 0xd05e4332_db43d551_ca37f2d8_e79c0f1c_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Pos,
            exponent: -163,
            mantissa: 0xdd09afcc_d4f0cb10_52ba9242_ac1d5edd_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Neg,
            exponent: -168,
            mantissa: 0xa06aa1e9_65ef81f2_99934254_ca7a2a5e_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Neg,
            exponent: -170,
            mantissa: 0x81e82c90_c2f8dff6_8271e19f_bd436c98_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Pos,
            exponent: -176,
            mantissa: 0xb6daf116_9df3a7f5_1f5cf2eb_8066ed0a_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Pos,
            exponent: -179,
            mantissa: 0xeba2d0fb_bc2c7f32_19276816_1ee04e25_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Neg,
            exponent: -184,
            mantissa: 0xa0a0e25a_d2082983_14c949b4_d68bc1de_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Neg,
            exponent: -187,
            mantissa: 0xa9fa7365_29fb279f_6953c002_e3a76e4b_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Pos,
            exponent: -193,
            mantissa: 0xe04e0b99_749a23e1_76be66ec_ce1a3474_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Pos,
            exponent: -196,
            mantissa: 0xc7cb8fd0_f9ba52b5_1b83db32_bce4b6d2_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Neg,
            exponent: -202,
            mantissa: 0xff1ce11f_a8db60f1_8658f813_8bd95612_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Neg,
            exponent: -205,
            mantissa: 0xc317a0cf_15d29ccf_4f37317c_903b3fa7_u128,
        },
    ],
    [
        DyadicFloat128 {
            sign: DyadicSign::Pos,
            exponent: -130,
            mantissa: 0xb19f3fbf_8298049d_325cfe55_0c698a0a_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Pos,
            exponent: -239,
            mantissa: 0xc8c38c33_2bf6f875_22c98375_e37d1c59_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Neg,
            exponent: -131,
            mantissa: 0xb139bc26_0e037020_0cb340eb_b34e109b_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Pos,
            exponent: -137,
            mantissa: 0xb1d74a56_f0c6a396_e97a10e5_c46c6a9d_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Pos,
            exponent: -135,
            mantissa: 0xea33329b_6445c92a_45d712a2_f909558c_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Neg,
            exponent: -140,
            mantissa: 0x8b6a941a_2e5087db_42ca693f_ea2f1dda_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Neg,
            exponent: -140,
            mantissa: 0xf606e6dd_be2fee46_bfa473a7_2a94d250_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Pos,
            exponent: -145,
            mantissa: 0x9a74b358_584fbc69_23f6168d_03043c66_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Pos,
            exponent: -145,
            mantissa: 0x89ae326e_16c67cb0_c3ec85eb_05a48a3e_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Neg,
            exponent: -151,
            mantissa: 0xafbde3c5_afe156ed_7b39509a_125d229a_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Neg,
            exponent: -152,
            mantissa: 0xbee1ad83_b204b168_597192a5_3a3d8823_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Pos,
            exponent: -158,
            mantissa: 0xf36b4b9b_0f529c9e_e5c98845_8511d2c1_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Pos,
            exponent: -159,
            mantissa: 0xb3cf1b82_4c3da125_acf92327_a4dc06f9_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Neg,
            exponent: -165,
            mantissa: 0xe2caf93c_526b5ea5_e7f53b5d_764e1442_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Neg,
            exponent: -167,
            mantissa: 0xf51b577a_b77601b8_a094432f_bc20bdd1_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Pos,
            exponent: -172,
            mantissa: 0x97eb35a7_310c2b1a_f9206859_72ae5bbf_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Pos,
            exponent: -175,
            mantissa: 0xfd0264fd_56760462_034c089e_81b47dc2_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Neg,
            exponent: -180,
            mantissa: 0x997a3f81_c2beb502_9793a7a0_d69ffb1e_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Neg,
            exponent: -183,
            mantissa: 0xccb2c23c_f5d985a4_7e8c1175_bf084a94_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Pos,
            exponent: -189,
            mantissa: 0xf2613019_72ef264b_1a3eacc3_ec3ee632_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Pos,
            exponent: -191,
            mantissa: 0x8558da3f_ce3c2bf5_7d383b83_16611557_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Neg,
            exponent: -197,
            mantissa: 0x99d1e3c1_4abab37e_ddd8dda4_bb443190_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Neg,
            exponent: -200,
            mantissa: 0x8ef0cdca_b802e990_ac82d0d0_17968245_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Pos,
            exponent: -206,
            mantissa: 0xa0718860_2f6b3b6e_2c4eeae6_e38123e8_u128,
        },
    ],
    [
        DyadicFloat128 {
            sign: DyadicSign::Neg,
            exponent: -238,
            mantissa: 0x8784c6e1_c723fc6d_b65880e8_14a09d7c_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Neg,
            exponent: -130,
            mantissa: 0xab32709d_b8310cb7_0a97b7ca_c166298e_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Pos,
            exponent: -136,
            mantissa: 0xf0b2aaf0_846e329a_f184a9e8_8438f593_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Pos,
            exponent: -133,
            mantissa: 0xe2f0d6cf_d97a00ff_3ad7eb67_07cd99ff_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Neg,
            exponent: -138,
            mantissa: 0x9e9b4fca_c72b4dec_c035e379_39b76d4c_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Neg,
            exponent: -137,
            mantissa: 0xb37a7613_d4b74a7d_ae1bc6cc_d074dae2_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Pos,
            exponent: -143,
            mantissa: 0xf86f0299_3296a499_e1f1c631_e892c3a1_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Pos,
            exponent: -142,
            mantissa: 0x8679321a_0a04efae_01c4570c_fea94964_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Neg,
            exponent: -148,
            mantissa: 0xb7a36a8f_2bb3741a_0106d1d7_305a26c6_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Neg,
            exponent: -149,
            mantissa: 0xea03fbe1_246a2ae1_187f38f7_a6826b77_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Pos,
            exponent: -154,
            mantissa: 0x9d19e6ce_fb34ac31_e41538b4_9e218cfc_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Pos,
            exponent: -155,
            mantissa: 0x84c80d85_40b72a60_1ac38665_f5f841e9_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Neg,
            exponent: -161,
            mantissa: 0xaec2b52c_221fa660_20e2273e_a4d4ed75_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Neg,
            exponent: -163,
            mantissa: 0xd3e5e6c7_e03ac0f8_17af6046_da1c1df8_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Pos,
            exponent: -168,
            mantissa: 0x885b1617_ea3fae35_5a70837e_aef43b92_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Pos,
            exponent: -171,
            mantissa: 0xfab504ed_0cfb7f03_b48a5d59_dbb1321c_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Neg,
            exponent: -176,
            mantissa: 0x9d73424d_9b2131cb_0796a5c1_53c8c7f6_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Neg,
            exponent: -179,
            mantissa: 0xe4bba925_32193f3b_56789880_821bc8d1_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Pos,
            exponent: -184,
            mantissa: 0x8bfe4edb_a677eb5f_113e848f_64873e24_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Pos,
            exponent: -187,
            mantissa: 0xa5df2136_bbd143cd_ff57de00_5b4a5297_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Neg,
            exponent: -193,
            mantissa: 0xc5ab2cb1_45471498_effb0fdc_a767440c_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Neg,
            exponent: -196,
            mantissa: 0xc3de6117_108b89cd_df274f11_9b86d2ab_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Pos,
            exponent: -202,
            mantissa: 0xe3167fa1_31356b8a_b37ee74b_a8e3a0de_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Pos,
            exponent: -205,
            mantissa: 0xc005aefa_da957015_6c3911a7_8ec32066_u128,
        },
    ],
    [
        DyadicFloat128 {
            sign: DyadicSign::Neg,
            exponent: -130,
            mantissa: 0xa5b8ea65_166348fb_799c893c_c9e86f8c_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Neg,
            exponent: -239,
            mantissa: 0x9d00c7c7_bb43b47c_4b6aa33a_f5f1e14c_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Pos,
            exponent: -131,
            mantissa: 0xa57122b4_bfdd3284_583435d1_d8776d8c_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Neg,
            exponent: -137,
            mantissa: 0x90af270d_2f8eb29f_eae5830e_85422d66_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Neg,
            exponent: -135,
            mantissa: 0xdb19f6b6_ce7b0892_cf64d590_bd454fb3_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Pos,
            exponent: -141,
            mantissa: 0xe3f8bda5_bedde9b3_3aa46f3e_7290e5ed_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Pos,
            exponent: -140,
            mantissa: 0xe700dc5b_551a5fe6_ccc531d8_4972fe5e_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Neg,
            exponent: -146,
            mantissa: 0xfe76b8d8_8f1ad60b_1a17af84_37808bb9_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Neg,
            exponent: -145,
            mantissa: 0x81e4eac1_1930e740_66d94eb6_27ca3d10_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Pos,
            exponent: -151,
            mantissa: 0x92284048_831d8d46_4fad6e7d_489f7493_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Pos,
            exponent: -152,
            mantissa: 0xb5169008_8f181a57_462e73ef_b574a77e_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Neg,
            exponent: -158,
            mantissa: 0xccad0c7c_73443c2f_60ef7996_485b4bde_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Neg,
            exponent: -159,
            mantissa: 0xab968d03_a7b237a1_729a5cb2_276191b8_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Pos,
            exponent: -165,
            mantissa: 0xc0ef81e0_c88f84b9_926d5882_443edc96_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Pos,
            exponent: -167,
            mantissa: 0xeb44af24_7a3ae5f6_47486a6a_9a0da273_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Neg,
            exponent: -172,
            mantissa: 0x82c43f85_9efbbd92_22ac652c_02a2f5f0_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Neg,
            exponent: -175,
            mantissa: 0xf43381ee_f8ddd9de_d8dfba53_f0d742f6_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Pos,
            exponent: -180,
            mantissa: 0x859e3865_39ef934d_b03d36e8_bbc28626_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Pos,
            exponent: -183,
            mantissa: 0xc6945d19_460c44a8_aed61ce3_7a9e794e_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Neg,
            exponent: -189,
            mantissa: 0xd5493aeb_4414fe15_a4775820_0a8c6471_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Neg,
            exponent: -191,
            mantissa: 0x81f44bba_554bb9b4_724ad789_8dc3abf5_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Pos,
            exponent: -197,
            mantissa: 0x88b3c902_5251c2d5_9ece2f8e_68fe63da_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Pos,
            exponent: -200,
            mantissa: 0x8bde9b50_272a463d_361a5ac5_90e17580_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Neg,
            exponent: -206,
            mantissa: 0x8fe283a2_ec2b2dbd_da31091a_f8a33e87_u128,
        },
    ],
    [
        DyadicFloat128 {
            sign: DyadicSign::Pos,
            exponent: -240,
            mantissa: 0xc25aaeb3_8f1d59f9_e30ee0b6_fbaa61f1_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Pos,
            exponent: -130,
            mantissa: 0xa07c83c9_b02dda62_df1ffa52_977c0f2e_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Neg,
            exponent: -136,
            mantissa: 0xc641983b_fdf57050_46b4d439_b7ba3b92_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Neg,
            exponent: -133,
            mantissa: 0xd5067021_0e8d3fab_85243d89_397f5de2_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Pos,
            exponent: -138,
            mantissa: 0x82fd2cc7_7ae8e47e_030ea64b_72dccc32_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Pos,
            exponent: -137,
            mantissa: 0xa8e986bc_55a933e4_10ab2d5d_2afbac16_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Neg,
            exponent: -143,
            mantissa: 0xce2c03b3_ae2fada9_9c99910c_0701dc38_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Neg,
            exponent: -143,
            mantissa: 0xfe0ddfab_f1d710e5_ad804442_0fec1254_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Pos,
            exponent: -148,
            mantissa: 0x996ff512_8c6731cc_2b42d9f0_f1a4e032_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Pos,
            exponent: -149,
            mantissa: 0xde109042_ee5e30c6_42c8a604_353eec6b_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Neg,
            exponent: -154,
            mantissa: 0x845ca36b_262e13a8_94a863be_bc97ba0d_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Neg,
            exponent: -156,
            mantissa: 0xfd47cec6_d39b4dd6_8072858d_bb73c961_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Pos,
            exponent: -161,
            mantissa: 0x949ff0a5_78e87686_efbff859_9729dfb5_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Pos,
            exponent: -163,
            mantissa: 0xcb2a265b_ad2c686d_7bd2a28b_f7e3b639_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Neg,
            exponent: -169,
            mantissa: 0xea3a83a4_71674095_9083b21d_b20c39da_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Neg,
            exponent: -171,
            mantissa: 0xf1a0c091_94839306_d3d0fe03_ae01ada8_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Pos,
            exponent: -176,
            mantissa: 0x8895524a_8c3d238d_f3c8fe09_329427e9_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Pos,
            exponent: -179,
            mantissa: 0xdd8b2b6e_33b81de3_bb60ff0e_6b1b6063_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Neg,
            exponent: -185,
            mantissa: 0xf53dcc7a_3c64aee9_14e4c0fc_7dc0c0f8_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Neg,
            exponent: -187,
            mantissa: 0xa1656bf5_57b1cf31_8693766a_e63f0568_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Pos,
            exponent: -193,
            mantissa: 0xaebcc6c3_f56e9fed_9879b6ed_47f1c62f_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Pos,
            exponent: -196,
            mantissa: 0xbf610ccf_4e86ce8e_a363239b_4dcda5fe_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Neg,
            exponent: -202,
            mantissa: 0xca79ab9f_01a98a0a_19ee144e_1c342663_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Neg,
            exponent: -205,
            mantissa: 0xbc527577_c2960727_7d4763dd_05f499b9_u128,
        },
    ],
    [
        DyadicFloat128 {
            sign: DyadicSign::Pos,
            exponent: -130,
            mantissa: 0x9befd47a_d2a850fa_ffe80ca9_af10cab1_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Pos,
            exponent: -241,
            mantissa: 0xb8d09af7_b9af05fd_f9400ad2_651efa5b_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Neg,
            exponent: -131,
            mantissa: 0x9bbae0d0_2784e389_9a42c3dc_e256042b_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Pos,
            exponent: -138,
            mantissa: 0xf15a4085_2363013d_56291330_1ed54959_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Pos,
            exponent: -135,
            mantissa: 0xce8a8d39_5c1d3c8d_12199d0a_88d4528e_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Neg,
            exponent: -141,
            mantissa: 0xbec72d04_87272c62_0e95b814_be9a96ea_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Neg,
            exponent: -140,
            mantissa: 0xda4d3213_9c35036f_de2db396_0ae27738_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Pos,
            exponent: -146,
            mantissa: 0xd6084abe_ee975992_48209e49_4f3baa08_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Pos,
            exponent: -146,
            mantissa: 0xf653103a_2e9e9892_1cb704fb_a387663e_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Neg,
            exponent: -152,
            mantissa: 0xf77e3f3c_aee6e934_38d41f60_80196aff_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Neg,
            exponent: -152,
            mantissa: 0xac62e950_6d20e762_650a007e_2e73f76c_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Pos,
            exponent: -158,
            mantissa: 0xaea390fa_2547ca16_827bdba6_9f89bc3f_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Pos,
            exponent: -159,
            mantissa: 0xa40e7118_a57ac992_a275db60_361a2597_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Neg,
            exponent: -165,
            mantissa: 0xa606993c_f4906402_277925c2_0e2d5e2f_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Neg,
            exponent: -167,
            mantissa: 0xe1f5383a_982b14c3_f0bdac72_f8293e6e_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Pos,
            exponent: -173,
            mantissa: 0xe30ead8b_2fbed152_b925c059_e394d8d5_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Pos,
            exponent: -175,
            mantissa: 0xeb95f8c4_db3d9fe1_fa58aafe_e625bf65_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Neg,
            exponent: -181,
            mantissa: 0xea13d4bb_4f4f6884_39963e97_4565ab5d_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Neg,
            exponent: -183,
            mantissa: 0xc064b29b_82b4c3ef_0976df10_24853d54_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Pos,
            exponent: -189,
            mantissa: 0xbc71b7e0_15cefb1b_2c9f92cc_d857557b_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Pos,
            exponent: -192,
            mantissa: 0xfcd160b4_5b2212cb_f0e58089_34de6e6e_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Neg,
            exponent: -198,
            mantissa: 0xf38fc17b_78114595_4ec2cb80_bbe1023c_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Neg,
            exponent: -200,
            mantissa: 0x888d608f_140b0e44_70bdbbe2_a906cd98_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Pos,
            exponent: -206,
            mantissa: 0x812d495a_f915286d_7dd85ebb_ae3626a8_u128,
        },
    ],
    [
        DyadicFloat128 {
            sign: DyadicSign::Neg,
            exponent: -238,
            mantissa: 0x89115da2_5e789ead_3470760c_3aa96c49_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Neg,
            exponent: -130,
            mantissa: 0x97903973_1c679cec_5c16fcf4_70491c86_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Pos,
            exponent: -136,
            mantissa: 0xa6f9045d_ea278be9_04ae6675_9a94f076_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Pos,
            exponent: -133,
            mantissa: 0xc95daf0f_0acfde7d_3b30fed1_855d6c16_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Neg,
            exponent: -139,
            mantissa: 0xdd0c0e03_44be6753_4c67625b_8e5908dd_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Neg,
            exponent: -137,
            mantissa: 0x9ff4eadd_5251e9a6_566cfb06_bed69ba8_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Pos,
            exponent: -143,
            mantissa: 0xae8bb015_80e55ce9_601ceda6_f07e2847_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Pos,
            exponent: -143,
            mantissa: 0xf1369ba3_ff4146b0_e2862e39_730c82e8_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Neg,
            exponent: -148,
            mantissa: 0x8284bb45_5b765af1_3395851a_0dbdfd5b_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Neg,
            exponent: -149,
            mantissa: 0xd3879196_ce83621d_b831cb29_fb33fdb8_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Pos,
            exponent: -155,
            mantissa: 0xe2849296_f19fcbe9_2cfda968_0b7b4d14_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Pos,
            exponent: -156,
            mantissa: 0xf22a7f70_388b72c8_79725f5c_5d493c8c_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Neg,
            exponent: -161,
            mantissa: 0x800ad831_5d092b85_abea6df5_680a35e4_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Neg,
            exponent: -163,
            mantissa: 0xc3074ad6_fc207939_93ecbbb4_4aab7805_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Pos,
            exponent: -169,
            mantissa: 0xcb48f487_bece03a3_94f12d88_b68d6a8a_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Pos,
            exponent: -171,
            mantissa: 0xe8e73ccb_d9e068e6_5ba4ccd2_06bf6deb_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Neg,
            exponent: -177,
            mantissa: 0xeee68e2a_67729c4f_539fc3b8_5de19e28_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Neg,
            exponent: -179,
            mantissa: 0xd6688501_aa96e31e_1ab15fa1_0ec50d49_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Pos,
            exponent: -185,
            mantissa: 0xd821d3d0_44b0cdb2_da1fe30a_8527f3a0_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Pos,
            exponent: -187,
            mantissa: 0x9ccdd177_ca545048_588a0022_cb5d1590_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Neg,
            exponent: -193,
            mantissa: 0x9b287f60_387b3b64_4c19d482_de58dabf_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Neg,
            exponent: -196,
            mantissa: 0xba9c5870_22278b10_bd09c87d_bd92f00b_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Pos,
            exponent: -202,
            mantissa: 0xb516e73b_07c52af5_81cd3c60_c558a907_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Pos,
            exponent: -205,
            mantissa: 0xb83f51dd_11171a04_0a9f92ed_e0be5bd9_u128,
        },
    ],
    [
        DyadicFloat128 {
            sign: DyadicSign::Neg,
            exponent: -130,
            mantissa: 0x93b4694e_34c9b373_f7b12ca7_4151e266_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Neg,
            exponent: -240,
            mantissa: 0xe47f169b_0fa83b07_afddbe6d_ef31b0e0_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Pos,
            exponent: -131,
            mantissa: 0x938c08b9_858774da_7595e565_1155534a_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Neg,
            exponent: -138,
            mantissa: 0xcd46fcb5_0912f0ea_1033e194_0643e2be_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Neg,
            exponent: -135,
            mantissa: 0xc3e40e78_e37e203d_032cb2f0_ff90410c_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Pos,
            exponent: -141,
            mantissa: 0xa2a4e6dd_bbcbc109_86868518_e7524943_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Pos,
            exponent: -140,
            mantissa: 0xcf692b40_8b70c47e_6d90f8c6_1a25eebb_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Neg,
            exponent: -146,
            mantissa: 0xb7209b78_19d67271_43b5fe97_eb867383_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Neg,
            exponent: -146,
            mantissa: 0xea999818_b66f3b25_704385cd_e6ce450b_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Pos,
            exponent: -152,
            mantissa: 0xd4c1dabc_16fe5182_b5a8a9ea_4dffbc35_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Pos,
            exponent: -152,
            mantissa: 0xa4a95d3e_2d0ee3ad_02fd4a38_7350c9f5_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Neg,
            exponent: -158,
            mantissa: 0x96f9d6b7_c1706753_c43063c9_79d36f15_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Neg,
            exponent: -159,
            mantissa: 0x9d387cd4_4e96115b_e3593e2b_e87bca2c_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Pos,
            exponent: -165,
            mantissa: 0x90704360_bfaedda1_397de176_37ac4669_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Pos,
            exponent: -167,
            mantissa: 0xd94d2aa6_08acf708_a51d56c9_9bcaa22d_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Neg,
            exponent: -173,
            mantissa: 0xc6df047a_54fb4ad4_a31e29ed_409bce1c_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Neg,
            exponent: -175,
            mantissa: 0xe35ecfe2_2e1bff81_c516d61c_e44d118d_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Pos,
            exponent: -181,
            mantissa: 0xce72d39a_d9488868_85d76cd0_c88020ae_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Pos,
            exponent: -183,
            mantissa: 0xba56f075_dd8251e5_b428d8e3_c2593018_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Neg,
            exponent: -189,
            mantissa: 0xa75c4701_58047697_211e05ff_6a8d20e6_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Neg,
            exponent: -192,
            mantissa: 0xf5b1f932_038fd7ea_7b9458e5_e6795350_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Pos,
            exponent: -198,
            mantissa: 0xd9c98ced_34062905_020072f2_233864bf_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Pos,
            exponent: -200,
            mantissa: 0x852192d7_4675343c_94c87a09_d27cd582_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Neg,
            exponent: -207,
            mantissa: 0xe88810e6_6cc1935d_0125547d_73afcfb3_u128,
        },
    ],
    [
        DyadicFloat128 {
            sign: DyadicSign::Pos,
            exponent: -239,
            mantissa: 0xf83d39d8_3138bb40_b6f76d4c_dcdfa000_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Pos,
            exponent: -130,
            mantissa: 0x8ffaf763_500e691f_9b16832f_3625febf_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Neg,
            exponent: -136,
            mantissa: 0x8f21c5b9_161dabf1_39c1a397_07c1fa88_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Neg,
            exponent: -133,
            mantissa: 0xbf6afffe_0e3b9ccb_69098254_8c6b7f7d_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Pos,
            exponent: -139,
            mantissa: 0xbdbccbdd_5ece4c79_44cae46a_e0604b7c_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Pos,
            exponent: -137,
            mantissa: 0x9840ef7c_b09626c6_25bf9967_e1e390d5_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Neg,
            exponent: -143,
            mantissa: 0x962faf6d_7e2746aa_e465b684_a7e258b3_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Neg,
            exponent: -143,
            mantissa: 0xe608d2cf_09b9cbb7_5cd25d56_b864a03b_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Pos,
            exponent: -149,
            mantissa: 0xe161d0dc_00a51aa2_b0126100_7e9bddb7_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Pos,
            exponent: -149,
            mantissa: 0xca368ed5_bde80598_bbd68c48_4f33503c_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Neg,
            exponent: -155,
            mantissa: 0xc46db0ca_35f31ce0_7412facc_dd63b6af_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Neg,
            exponent: -156,
            mantissa: 0xe82699df_2c073d53_becf60d1_0fc68ad6_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Pos,
            exponent: -162,
            mantissa: 0xdf321c36_962e7c50_1c79e5d5_09237c05_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Pos,
            exponent: -163,
            mantissa: 0xbb8a1707_224bd674_816320c1_aa929496_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Neg,
            exponent: -169,
            mantissa: 0xb22c23b3_3cb4b437_a6b6047e_fc056b41_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Neg,
            exponent: -171,
            mantissa: 0xe0af4b59_2d6dce21_fde870d1_758b1197_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Pos,
            exponent: -177,
            mantissa: 0xd2a2f376_38e90628_0d7fa19a_6351930f_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Pos,
            exponent: -179,
            mantissa: 0xcf854dba_701a75d2_df9d3ba3_4c896fc0_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Neg,
            exponent: -185,
            mantissa: 0xbfba8f55_6aca1ad9_189963a5_2a291194_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Neg,
            exponent: -187,
            mantissa: 0x9841dfdd_48059730_63ecbcd9_6a71697f_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Pos,
            exponent: -193,
            mantissa: 0x8a7bc433_ccf2efca_5eba13a5_b7206ac6_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Pos,
            exponent: -196,
            mantissa: 0xb5c3ab51_a30f4d22_b2acc54f_89d2fb7b_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Neg,
            exponent: -202,
            mantissa: 0xa299e097_18a0807f_e3b6de14_d92b9198_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Neg,
            exponent: -205,
            mantissa: 0xb3fe00e9_363814d9_18a3a98b_2fd4a5b8_u128,
        },
    ],
    [
        DyadicFloat128 {
            sign: DyadicSign::Pos,
            exponent: -130,
            mantissa: 0x8ca75d3a_d997c94f_a3fb220c_4974174e_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Pos,
            exponent: -239,
            mantissa: 0xee4d1383_518482d1_df48df71_98d0d96e_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Neg,
            exponent: -131,
            mantissa: 0x8c87bee1_39b154c1_4045da6a_bddbf815_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Pos,
            exponent: -138,
            mantissa: 0xb15da3d2_e2fbf669_dbb8b935_e4e08cb8_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Pos,
            exponent: -135,
            mantissa: 0xbab77907_ae97120d_db7c176d_13bd6492_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Neg,
            exponent: -141,
            mantissa: 0x8cc585d4_bf67d5db_9c877612_15eef38c_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Neg,
            exponent: -140,
            mantissa: 0xc5f467c9_ac6aa8fe_57a49ad0_ba92198d_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Pos,
            exponent: -146,
            mantissa: 0x9eeb7be4_661e0050_6b697051_3c317311_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Pos,
            exponent: -146,
            mantissa: 0xe04e4077_3fcd77db_0d0f3af1_21cf2dcc_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Neg,
            exponent: -152,
            mantissa: 0xb9472352_28f19280_0a72ba93_c076b791_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Neg,
            exponent: -152,
            mantissa: 0x9dc8889a_843112ac_666330ae_d517e093_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Pos,
            exponent: -158,
            mantissa: 0x8407eed6_b15fbdcd_4b32d81e_4f67d341_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Pos,
            exponent: -159,
            mantissa: 0x97088f44_6d70f935_2e1bf38c_bc8398c5_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Neg,
            exponent: -166,
            mantissa: 0xfdd7445e_d590aa21_246dacfa_cbc8a11d_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Neg,
            exponent: -167,
            mantissa: 0xd1527df5_43453f11_cf8e9d1c_16461259_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Pos,
            exponent: -173,
            mantissa: 0xafaa074b_b4f45b1c_72c66902_4c878492_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Pos,
            exponent: -175,
            mantissa: 0xdba60593_436f9aed_4cdaf377_21f319bf_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Neg,
            exponent: -181,
            mantissa: 0xb75bad3e_e89d11c1_d434c751_62e14f50_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Neg,
            exponent: -183,
            mantissa: 0xb48860ee_432e1638_01eabc8d_a67f4621_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Pos,
            exponent: -189,
            mantissa: 0x9579b248_5ae386bc_0a86f27e_2fd4188d_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Pos,
            exponent: -192,
            mantissa: 0xeeb81b75_246aaadf_13a5c800_83ddda58_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Neg,
            exponent: -198,
            mantissa: 0xc39a78fe_e8834340_916abf22_66458a4b_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Neg,
            exponent: -200,
            mantissa: 0x81b579b4_61cfe07e_96e87e1e_1612578f_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Pos,
            exponent: -207,
            mantissa: 0xd1fec508_ef7a8cd8_497c5a10_6355bce1_u128,
        },
    ],
    [
        DyadicFloat128 {
            sign: DyadicSign::Neg,
            exponent: -241,
            mantissa: 0xbd4c8e7c_be71728b_69c35e37_0bbeec91_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Neg,
            exponent: -130,
            mantissa: 0x896eabdf_8c56cc31_c76fa833_085440a5_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Pos,
            exponent: -137,
            mantissa: 0xf8f0f3f9_c9f41665_757d83b5_d2c51a3d_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Pos,
            exponent: -133,
            mantissa: 0xb6cd7f44_180fd250_df2657d4_563ad2b2_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Neg,
            exponent: -139,
            mantissa: 0xa529c524_14bcda6c_0661bdff_5403b7f7_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Neg,
            exponent: -137,
            mantissa: 0x918b1285_a571b607_c7222abf_4ef3eb16_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Pos,
            exponent: -143,
            mantissa: 0x82f88abb_da313567_d80982a1_4a9750ab_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Pos,
            exponent: -143,
            mantissa: 0xdc35d6a1_7e11005d_c7be6484_6fe78a17_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Neg,
            exponent: -149,
            mantissa: 0xc50d9d4f_493a472c_a40d18a2_204155ae_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Neg,
            exponent: -149,
            mantissa: 0xc1ee579f_c7e625a5_10f6ec0a_4222bd5a_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Pos,
            exponent: -155,
            mantissa: 0xac4d3e52_fee788e7_77bd9bd2_0a543044_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Pos,
            exponent: -156,
            mantissa: 0xdf1f7194_c5dccb2f_fe892b2a_d471911d_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Neg,
            exponent: -162,
            mantissa: 0xc489f8e8_7f894d11_e9caca9f_46b3d1f5_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Neg,
            exponent: -163,
            mantissa: 0xb4ae1c33_304096cf_3e06b6da_d757ae51_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Pos,
            exponent: -169,
            mantissa: 0x9d92e9b2_3b17344e_15e704d1_6450f526_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Pos,
            exponent: -171,
            mantissa: 0xd9062173_21666533_65d3e5fc_39bd30ea_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Neg,
            exponent: -177,
            mantissa: 0xbb282995_a78efc87_3e4838e5_f795fbd2_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Neg,
            exponent: -179,
            mantissa: 0xc8fa8a5b_85a05040_a0f95b39_6d1aa4e2_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Pos,
            exponent: -185,
            mantissa: 0xab30d562_6e7ef395_c7c311d8_d2284b13_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Pos,
            exponent: -187,
            mantissa: 0x93da46ca_3c582f7a_933e8338_a0a117a5_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Neg,
            exponent: -194,
            mantissa: 0xf8888af8_f03aff5f_74b80e6e_1ee58d76_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Neg,
            exponent: -196,
            mantissa: 0xb0f8de81_0985665f_4a1bb50c_315e1ced_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Pos,
            exponent: -202,
            mantissa: 0x92a32378_7cbe3b69_c9e062ad_32ec9cdb_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Pos,
            exponent: -205,
            mantissa: 0xafb1df2b_1a66921d_4d548c92_59541de0_u128,
        },
    ],
    [
        DyadicFloat128 {
            sign: DyadicSign::Neg,
            exponent: -130,
            mantissa: 0x86869b23_9f4c602e_e6d63ead_ae0134e1_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Neg,
            exponent: -238,
            mantissa: 0x92418a2a_06587fc6_c4770ae8_4503a541_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Pos,
            exponent: -131,
            mantissa: 0x866d4cba_55f155d0_58888e31_137011b4_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Neg,
            exponent: -138,
            mantissa: 0x9b3f9c79_00e12a54_876b1365_3f16c747_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Neg,
            exponent: -135,
            mantissa: 0xb2b5baf1_e12172fa_264600bd_5cedb782_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Pos,
            exponent: -142,
            mantissa: 0xf6c155fb_a44f894b_9f5a7acb_7cc0c734_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Pos,
            exponent: -140,
            mantissa: 0xbda72dc1_f759ac9e_8241c6fd_adbe9b66_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Neg,
            exponent: -146,
            mantissa: 0x8b8fe87d_b38737f4_ef6dc4ec_a5f53278_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Neg,
            exponent: -146,
            mantissa: 0xd7315a56_80bfda71_6fac9650_3ed039ee_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Pos,
            exponent: -152,
            mantissa: 0xa32424ad_3da4f370_7ac34dfe_e439a1ea_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Pos,
            exponent: -152,
            mantissa: 0x97a152bb_89ad6abe_0425b6fc_b84a793c_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Neg,
            exponent: -159,
            mantissa: 0xe94373fd_1932b0b5_87eaa66c_ead384fe_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Neg,
            exponent: -159,
            mantissa: 0x916de5fb_b50e3808_6501588d_5d88ba4b_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Pos,
            exponent: -166,
            mantissa: 0xe1115666_069692fd_0646a471_744209b0_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Pos,
            exponent: -167,
            mantissa: 0xc9fe3fb6_469aa956_929761e3_d7760656_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Neg,
            exponent: -173,
            mantissa: 0x9c63f0fb_34759c7b_91dba889_0e031794_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Neg,
            exponent: -175,
            mantissa: 0xd47277e3_6f5b0a0e_086978fc_197e59f4_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Pos,
            exponent: -181,
            mantissa: 0xa3f40475_b9671ab0_2dfda08b_0be60dbf_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Pos,
            exponent: -183,
            mantissa: 0xaf0783af_b99c95f2_c636907c_485ae2d7_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Neg,
            exponent: -189,
            mantissa: 0x864335d2_1d75d9b6_81bd123a_b3439c32_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Neg,
            exponent: -192,
            mantissa: 0xe7ff0e07_01041934_4f73adf2_48659072_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Pos,
            exponent: -198,
            mantissa: 0xb0817cc4_70543a3b_e575af18_009c27a6_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Pos,
            exponent: -201,
            mantissa: 0xfcb52bf7_e84fd582_7419088f_9f9935b7_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Neg,
            exponent: -207,
            mantissa: 0xbe5d62fd_abd48442_1fcfec82_258b47e2_u128,
        },
    ],
    [
        DyadicFloat128 {
            sign: DyadicSign::Neg,
            exponent: -242,
            mantissa: 0x806de2ef_cd134f9b_9313fc0e_22bbf7ac_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Pos,
            exponent: -130,
            mantissa: 0x83b41366_160c8a34_27c5fdca_d0ed270e_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Neg,
            exponent: -137,
            mantissa: 0xdb1442f0_0386306d_e27fed99_fd60dd3d_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Neg,
            exponent: -133,
            mantissa: 0xaf3fa980_0f258544_82c51bde_7a3b670b_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Pos,
            exponent: -139,
            mantissa: 0x9175f5c6_e1100c48_e29c5bf9_a1a93426_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Pos,
            exponent: -137,
            mantissa: 0x8ba2253d_02667d1b_3ce2f1d5_60783cd7_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Neg,
            exponent: -144,
            mantissa: 0xe70378f8_fe6bee4f_e17e373c_8f2df889_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Neg,
            exponent: -143,
            mantissa: 0xd3803735_6ce7eefe_97e7b04d_5d100d9c_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Pos,
            exponent: -149,
            mantissa: 0xae217816_78ae509d_61b5cf7d_d9a6c8ae_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Pos,
            exponent: -149,
            mantissa: 0xba8652f0_9b35a062_da6da628_bcbb5cf0_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Neg,
            exponent: -155,
            mantissa: 0x98a4c16f_bf551108_a205d6fb_d27715b7_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Neg,
            exponent: -156,
            mantissa: 0xd6f6ba98_359d960d_367283a4_f4fac84f_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Pos,
            exponent: -162,
            mantissa: 0xaea3c23e_c75f87d5_5c4d3446_332fee7c_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Pos,
            exponent: -163,
            mantissa: 0xae67a269_43d1086d_1d7b943e_6bc8273b_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Neg,
            exponent: -169,
            mantissa: 0x8c7e71f0_20144f35_b7a9bdca_502ad000_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Neg,
            exponent: -171,
            mantissa: 0xd1eb5e56_861fe2a4_b982e0fb_a9b39611_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Pos,
            exponent: -177,
            mantissa: 0xa77debb3_93e02915_f1fc8118_8f78ec68_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Pos,
            exponent: -179,
            mantissa: 0xc2d258b1_fe602299_e307c98d_ddda3f32_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Neg,
            exponent: -185,
            mantissa: 0x99ceef65_3ce851b3_5201a650_447a3e60_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Neg,
            exponent: -187,
            mantissa: 0x8fa45a3d_e220cbe5_fb51cae9_e3120345_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Pos,
            exponent: -194,
            mantissa: 0xe035eb38_85a2d8eb_de2cf7fd_30cbea31_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Pos,
            exponent: -196,
            mantissa: 0xac50d6cf_32a244b8_f49424fb_c2979395_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Neg,
            exponent: -202,
            mantissa: 0x84d5f15c_5f0768a7_2912defa_0d9b6ce3_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Neg,
            exponent: -205,
            mantissa: 0xab72baec_b2f44342_68779f49_c2c22957_u128,
        },
    ],
    [
        DyadicFloat128 {
            sign: DyadicSign::Pos,
            exponent: -130,
            mantissa: 0x8122ab3a_d5b69234_dc04aa73_01d4dfff_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Pos,
            exponent: -240,
            mantissa: 0xdaf58313_5326a309_645674df_fa6b7260_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Neg,
            exponent: -131,
            mantissa: 0x810e0aad_3906fb4b_166107cf_c795a072_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Pos,
            exponent: -138,
            mantissa: 0x895f2b7e_0a24ccbc_c7df8c59_e85dedad_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Pos,
            exponent: -135,
            mantissa: 0xaba4eaab_56b0e11b_ad3e8f9a_bb97af16_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Neg,
            exponent: -142,
            mantissa: 0xda8f8fcd_f549f055_884dc051_f7f869e1_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Neg,
            exponent: -140,
            mantissa: 0xb64b503f_11b5e4fc_669ac92a_f5e6f414_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Pos,
            exponent: -147,
            mantissa: 0xf79d3d5e_ad6357b9_fe6c584d_4b50fc2a_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Pos,
            exponent: -146,
            mantissa: 0xcf103e12_ba19cc12_25ea2a26_8b7bad41_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Neg,
            exponent: -152,
            mantissa: 0x9105cb77_7c02c10a_ed16a933_26e6768f_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Neg,
            exponent: -152,
            mantissa: 0x9218be59_4b39b853_5f8f551b_cb0ea380_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Pos,
            exponent: -159,
            mantissa: 0xcfe17966_9c8bf5c1_180985ba_826e1367_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Pos,
            exponent: -159,
            mantissa: 0x8c5741a3_f412b280_7cd9f228_74895190_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Neg,
            exponent: -166,
            mantissa: 0xc92a0211_f96ecc1d_070557dc_7b4d6fb8_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Neg,
            exponent: -167,
            mantissa: 0xc343ee6e_171d46b5_bf7ffd70_a4d80eae_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Pos,
            exponent: -173,
            mantissa: 0x8c3c8488_1ea3a120_1bd34791_753786c3_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Pos,
            exponent: -175,
            mantissa: 0xcdc1b143_96192144_a97e12df_c962edf4_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Neg,
            exponent: -181,
            mantissa: 0x9388f66d_5764fc4b_5d714f18_fa3c3071_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Neg,
            exponent: -183,
            mantissa: 0xa9d9a151_240a6f28_b852aec6_6a42d019_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Pos,
            exponent: -190,
            mantissa: 0xf2875d17_723b00d6_b3b11f20_3344f526_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Pos,
            exponent: -192,
            mantissa: 0xe19552fa_cb530f29_70655986_8ca9f71a_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Neg,
            exponent: -198,
            mantissa: 0xa00678be_97e69aaf_678351f6_58d83494_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Neg,
            exponent: -201,
            mantissa: 0xf6377a2b_673fc5ab_6e596f5d_78b4291f_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Pos,
            exponent: -207,
            mantissa: 0xad417f85_56c5e40b_ecee2136_3237acd7_u128,
        },
    ],
    [
        DyadicFloat128 {
            sign: DyadicSign::Pos,
            exponent: -241,
            mantissa: 0xb51809c1_ffc46537_e970e3fa_bc1bd84d_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Neg,
            exponent: -131,
            mantissa: 0xfd45a0b8_8e41cc50_677d1d1d_4c0668f5_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Pos,
            exponent: -137,
            mantissa: 0xc2be9cb4_ccbe8ac4_3349144b_741bc83b_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Pos,
            exponent: -133,
            mantissa: 0xa88e36d6_d5508cd4_6619606f_48a264fc_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Neg,
            exponent: -139,
            mantissa: 0x816144ed_edea75d8_6754a3b7_92e72cfd_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Neg,
            exponent: -137,
            mantissa: 0x866111c6_914fc8f4_5567c9ab_a27a60f1_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Pos,
            exponent: -144,
            mantissa: 0xcdb26ae3_1b347a63_e3923a4e_ae7a7df1_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Pos,
            exponent: -143,
            mantissa: 0xcbb7d811_f87bc530_10a55ac0_8aa492a5_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Neg,
            exponent: -149,
            mantissa: 0x9b49dd05_add38164_31725df2_0effa425_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Neg,
            exponent: -149,
            mantissa: 0xb3dca94c_39a86217_8a404e29_061cfaba_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Pos,
            exponent: -155,
            mantissa: 0x8865bcd1_56ed5c0a_785d3e62_31e888c0_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Pos,
            exponent: -156,
            mantissa: 0xcf903cfc_71cbf9bc_10a636cc_a3b8df69_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Neg,
            exponent: -162,
            mantissa: 0x9c6cabf5_529d6be3_e6d0cf22_405018d9_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Neg,
            exponent: -163,
            mantissa: 0xa8a88f49_dd3a5cb7_3696541e_686fbc3f_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Pos,
            exponent: -170,
            mantissa: 0xfc5da6ce_979d1e89_34a03ae2_2a30dbfa_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Pos,
            exponent: -171,
            mantissa: 0xcb57f410_aad1e2be_d78620c6_2c92c720_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Neg,
            exponent: -177,
            mantissa: 0x96e1e2d2_097e4313_83a68ead_b6e956d7_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Neg,
            exponent: -179,
            mantissa: 0xbd0e46f6_3b7fb61a_5365a2bc_98bb2f7e_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Pos,
            exponent: -185,
            mantissa: 0x8b0061eb_bce6634c_5726f50a_bb24a780_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Pos,
            exponent: -187,
            mantissa: 0x8ba6325f_858b5022_2d147b10_7c49d72c_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Neg,
            exponent: -194,
            mantissa: 0xcb4f1ba9_277be060_34477226_3318b3f5_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Neg,
            exponent: -196,
            mantissa: 0xa7d781f1_e616a937_a7c5fb77_7c8eca75_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Pos,
            exponent: -203,
            mantissa: 0xf1bd31f1_3b49cb3f_236ec473_b921e21c_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Pos,
            exponent: -205,
            mantissa: 0xa74fc728_d4b98b03_fff3c104_c8974a9b_u128,
        },
    ],
    [
        DyadicFloat128 {
            sign: DyadicSign::Neg,
            exponent: -131,
            mantissa: 0xf8b0e861_45a45de3_7b5ebd4e_bbb05f9e_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Neg,
            exponent: -239,
            mantissa: 0xe6e2aa92_8c862674_d9d5165d_0e88c30c_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Pos,
            exponent: -132,
            mantissa: 0xf88ec1bd_537b2243_c6c1898e_1ba963b0_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Neg,
            exponent: -139,
            mantissa: 0xf55bb6d2_683cfdce_16e79570_7758a174_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Neg,
            exponent: -135,
            mantissa: 0xa55994f8_33d74c28_388ced8a_7652570b_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Pos,
            exponent: -142,
            mantissa: 0xc356d2bd_e0e287e1_7f3c4f29_09392051_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Pos,
            exponent: -140,
            mantissa: 0xafb73c78_8d5cd207_2c1dedcd_7804e24e_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Neg,
            exponent: -147,
            mantissa: 0xdd938faa_644b2431_1cfe7537_61489a9d_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Neg,
            exponent: -146,
            mantissa: 0xc7c29946_1365b630_2e2825a5_d12651d6_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Pos,
            exponent: -152,
            mantissa: 0x81fc1575_ef560ab4_bae20686_45fe7ba5_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Pos,
            exponent: -152,
            mantissa: 0x8d180813_cd0caa11_ea267610_1bdd9ec4_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Neg,
            exponent: -159,
            mantissa: 0xbab30404_e6f6f68b_949091f8_4900b782_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Neg,
            exponent: -159,
            mantissa: 0x87b498bb_a502f924_b566432a_314f2458_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Pos,
            exponent: -166,
            mantissa: 0xb5173093_383857b0_fd1253eb_43510c42_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Pos,
            exponent: -167,
            mantissa: 0xbd154701_88e1ccd5_e430be8b_6b8d06ba_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Neg,
            exponent: -174,
            mantissa: 0xfd26471e_700a0209_1efa3e56_4a42fdf0_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Neg,
            exponent: -175,
            mantissa: 0xc78c937a_719d6f61_acc86e90_73c9b620_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Pos,
            exponent: -181,
            mantissa: 0x858b2d93_0ade40fb_ddf01297_034082bd_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Pos,
            exponent: -183,
            mantissa: 0xa4fe994d_49e86e9c_f411dc61_bff51499_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Neg,
            exponent: -190,
            mantissa: 0xdc3251c9_acbde848_fcdec2be_cd9d8011_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Neg,
            exponent: -192,
            mantissa: 0xdb811b67_f35eb3d4_d6b171e4_ef696ae6_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Pos,
            exponent: -198,
            mantissa: 0x91c00de5_d4b56e00_648176d2_7c8cf3f8_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Pos,
            exponent: -201,
            mantissa: 0xeffe16f5_98848d73_e34e4300_64fa587c_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Neg,
            exponent: -207,
            mantissa: 0x9e4ff1dc_2532c47a_e0dbe107_efe8859f_u128,
        },
    ],
    [
        DyadicFloat128 {
            sign: DyadicSign::Pos,
            exponent: -241,
            mantissa: 0xba332764_df689c84_c10ddf0e_70688318_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Pos,
            exponent: -131,
            mantissa: 0xf4393ed5_1ed76810_20b40165_94ac1611_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Neg,
            exponent: -137,
            mantissa: 0xae9a9f14_2a51b9a4_d8736187_88ef2646_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Neg,
            exponent: -133,
            mantissa: 0xa2926a40_9e612593_cf92759b_5d857dec_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Pos,
            exponent: -140,
            mantissa: 0xe81baba7_141b7e8f_8ca05676_a9d4fdfc_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Pos,
            exponent: -137,
            mantissa: 0x81ab5dba_3d3b1afe_3af81e61_681e585e_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Neg,
            exponent: -144,
            mantissa: 0xb8ab5fe6_779baecf_ea3f80d2_6e5175f4_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Neg,
            exponent: -143,
            mantissa: 0xc4b6bee4_0cfd7b2a_47e9f366_c1287d28_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Pos,
            exponent: -149,
            mantissa: 0x8b9636d6_fa0a4a0f_e1c8a3d7_854ab6e0_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Pos,
            exponent: -149,
            mantissa: 0xadd57477_e1718596_52d980d1_c5233a51_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Neg,
            exponent: -156,
            mantissa: 0xf59a3f5a_6b4a065f_e28d3bac_efa15bc0_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Neg,
            exponent: -156,
            mantissa: 0xc8d30539_5106fc78_459e133e_35e42498_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Pos,
            exponent: -162,
            mantissa: 0x8d1a2d4d_352fb7f0_4b85775f_eab8178d_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Pos,
            exponent: -163,
            mantissa: 0xa362b580_c835f08b_9c496e4a_d32c92e3_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Neg,
            exponent: -170,
            mantissa: 0xe425d9bb_3bc1bee2_4c8db876_64b08092_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Neg,
            exponent: -171,
            mantissa: 0xc541f372_60b45c35_2d941369_cd69d179_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Pos,
            exponent: -177,
            mantissa: 0x88bcb680_bd295e76_7944f096_72b6bbe1_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Pos,
            exponent: -179,
            mantissa: 0xb7ab3ae7_9df569e3_5cdfe4ad_2df8b26f_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Neg,
            exponent: -186,
            mantissa: 0xfc9b6833_184c1218_708a184c_d6d8915d_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Neg,
            exponent: -187,
            mantissa: 0x87e17c02_164fcb46_b26d6285_5d6ae37f_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Pos,
            exponent: -194,
            mantissa: 0xb940e13b_7c86b325_ee7c4bc5_27ba167c_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Pos,
            exponent: -196,
            mantissa: 0xa392db9e_43c11e6f_f6afee11_8e7acdf3_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Neg,
            exponent: -203,
            mantissa: 0xdce9148c_6775d05b_0b570a62_11ee0e89_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Neg,
            exponent: -205,
            mantissa: 0xa3521b78_b6feff7f_0f2e5f07_303d562d_u128,
        },
    ],
    [
        DyadicFloat128 {
            sign: DyadicSign::Pos,
            exponent: -131,
            mantissa: 0xf01abe0a_c588c496_6ce1e587_67ab116d_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Neg,
            exponent: -240,
            mantissa: 0xfe3e9fbb_0644dd86_281eee32_9c04ef8b_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Neg,
            exponent: -132,
            mantissa: 0xeffe17e0_d48faa20_0e599710_c60aa110_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Pos,
            exponent: -139,
            mantissa: 0xdcdc1570_3ed222ed_3c181b7f_c38ee586_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Pos,
            exponent: -135,
            mantissa: 0x9fb27029_90562ed2_72405e6b_f69f9964_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Neg,
            exponent: -142,
            mantissa: 0xaff25b36_78cec403_76123610_b258852f_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Neg,
            exponent: -140,
            mantissa: 0xa9ca8c70_d807aaa8_1d108d1a_54ad4d59_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Pos,
            exponent: -147,
            mantissa: 0xc7c680e2_43c82734_1cf85a13_01997f9d_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Pos,
            exponent: -146,
            mantissa: 0xc12822dc_bf16de1e_89ad782a_a379926f_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Neg,
            exponent: -153,
            mantissa: 0xeab2f9dd_b0ef509a_4e93c89b_fde8e8c0_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Neg,
            exponent: -152,
            mantissa: 0x888c32c3_7a7b9e91_0e8dee51_ba4be617_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Pos,
            exponent: -159,
            mantissa: 0xa8d3348a_c0c8dbab_8ea60d1d_df68d71a_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Pos,
            exponent: -159,
            mantissa: 0x8377a97b_38ad5c6e_a177c50d_90e4846f_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Neg,
            exponent: -166,
            mantissa: 0xa410af4a_e553b720_e42c0a02_d4063964_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Neg,
            exponent: -167,
            mantissa: 0xb7642110_0cd86a55_d458aa8f_bb1c32a0_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Pos,
            exponent: -174,
            mantissa: 0xe5d7afe0_9c99c6a3_1bcd40b1_5bccca76_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Pos,
            exponent: -175,
            mantissa: 0xc1c9ffd2_7ad7e590_a3b54939_fc328c9b_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Neg,
            exponent: -182,
            mantissa: 0xf31197db_f049db5c_ff628043_381507d1_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Neg,
            exponent: -183,
            mantissa: 0xa0734de5_900640a7_a0c706ed_90a44794_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Pos,
            exponent: -190,
            mantissa: 0xc8e68ef0_175fab60_cb3e7c15_f23edcb3_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Pos,
            exponent: -192,
            mantissa: 0xd5c36b92_b06dc399_e45df0f7_3f4f0f5b_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Neg,
            exponent: -198,
            mantissa: 0x8554a4ee_b02a8bc1_635a38a0_7165620a_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Neg,
            exponent: -201,
            mantissa: 0xea0ed91c_a1c760cf_f806447a_5511f78c_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Pos,
            exponent: -207,
            mantissa: 0x91381947_6481ec80_3c5a5220_faeb37ef_u128,
        },
    ],
    [
        DyadicFloat128 {
            sign: DyadicSign::Neg,
            exponent: -238,
            mantissa: 0x82e0cb4e_a70bdb8b_c945c4ea_e09dad40_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Neg,
            exponent: -131,
            mantissa: 0xec149d52_ae8c7b96_3e5132f9_47100019_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Pos,
            exponent: -137,
            mantissa: 0x9db5f5c1_f92caf78_fa4771f8_7f2d8e91_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Pos,
            exponent: -133,
            mantissa: 0x9d2e65e0_9533dad8_1238bc57_a751d1c5_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Neg,
            exponent: -140,
            mantissa: 0xd1bb2e93_552154fe_895f329c_ab0a7549_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Neg,
            exponent: -138,
            mantissa: 0xfad599ba_4610a75c_dcd5a605_edbfa523_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Pos,
            exponent: -144,
            mantissa: 0xa6fb78e1_2bd2dd27_8a5c8758_24051819_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Pos,
            exponent: -143,
            mantissa: 0xbe5ea4e1_8adefd91_b5de5238_25231ee9_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Neg,
            exponent: -150,
            mantissa: 0xfcb04c67_838ba620_381d646b_50c23d5a_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Neg,
            exponent: -149,
            mantissa: 0xa859beea_1a98f794_4f2a6ca1_e3a9c376_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Pos,
            exponent: -156,
            mantissa: 0xde973202_e302e8d2_f7f46686_864f2769_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Pos,
            exponent: -156,
            mantissa: 0xc2a986f8_6a5b7e60_25031c6e_edc31335_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Neg,
            exponent: -162,
            mantissa: 0x8014f10c_98792427_f914d66e_e555e723_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Neg,
            exponent: -163,
            mantissa: 0x9e88d75d_398d0bed_7dcd540d_77d4756d_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Pos,
            exponent: -170,
            mantissa: 0xcf783bf0_32e8b439_85328c20_710bc9de_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Pos,
            exponent: -171,
            mantissa: 0xbf9e908e_c0652e89_a09daa34_c2f1fdea_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Neg,
            exponent: -178,
            mantissa: 0xf930be76_8277b9ff_eb5c7b60_6b6532c2_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Neg,
            exponent: -179,
            mantissa: 0xb2a3c2fc_823a6ba8_67bf0ab6_03334e29_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Pos,
            exponent: -186,
            mantissa: 0xe6aea485_d79b0828_fe053b7b_e91541c2_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Pos,
            exponent: -187,
            mantissa: 0x84554c68_c42fd2dd_b0730acc_15e1a190_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Neg,
            exponent: -194,
            mantissa: 0xa9921dc7_a5167012_00c557f4_0390db2f_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Neg,
            exponent: -196,
            mantissa: 0x9f850f02_707c244a_e5c12aab_0bbc8589_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Pos,
            exponent: -203,
            mantissa: 0xcab47d2a_a3cc01e2_a379bc1d_48040f74_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Pos,
            exponent: -205,
            mantissa: 0x9f7e9f08_547d7683_ff4d3bf6_f7dc4b03_u128,
        },
    ],
    [
        DyadicFloat128 {
            sign: DyadicSign::Neg,
            exponent: -131,
            mantissa: 0xe859b72b_9ba2c2f5_8042319c_06883bcd_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Pos,
            exponent: -241,
            mantissa: 0xafdbbe97_071358ed_697d0d20_02f10c7e_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Pos,
            exponent: -132,
            mantissa: 0xe841671e_35acd4d6_272cccf0_73164e3f_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Neg,
            exponent: -139,
            mantissa: 0xc82e8062_f33ffc84_9290524c_24fc847f_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Neg,
            exponent: -135,
            mantissa: 0x9a95839f_ed63da9b_1fbbad18_93ad26c0_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Pos,
            exponent: -142,
            mantissa: 0x9f8e67f6_17e44314_721d6f85_72303a30_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Pos,
            exponent: -140,
            mantissa: 0xa46ba2c1_c1a5dfc1_16da6920_5e7ddd43_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Neg,
            exponent: -147,
            mantissa: 0xb54f786c_4b5d6d92_76a39099_0cb762d0_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Neg,
            exponent: -146,
            mantissa: 0xbb26cbcc_6b54030d_a1c1a37e_a030f014_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Pos,
            exponent: -153,
            mantissa: 0xd53c2eb3_9b7add73_eb73bc3f_961b2668_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Pos,
            exponent: -152,
            mantissa: 0x846571b0_8c66da6b_fa7ed496_40220b19_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Neg,
            exponent: -159,
            mantissa: 0x9996efd9_ce68044f_50fcc29b_d14c64c4_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Neg,
            exponent: -160,
            mantissa: 0xff2817f8_d47866c8_05fa755d_408d1141_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Pos,
            exponent: -166,
            mantissa: 0x957e419a_47758fd4_edb3b769_1a7c894d_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Pos,
            exponent: -167,
            mantissa: 0xb2234605_18a88a57_d286e3c5_62b20429_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Neg,
            exponent: -174,
            mantissa: 0xd1ccf424_e1f94d23_6b5bb816_4aaf0287_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Neg,
            exponent: -175,
            mantissa: 0xbc704bb9_570f3d81_9f1580ef_cdb6d79c_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Pos,
            exponent: -182,
            mantissa: 0xde4ef533_d3a28f64_3915a3bb_6cd84638_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Pos,
            exponent: -183,
            mantissa: 0x9c331f0c_020bee0b_f3bfc2ff_d34270bc_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Neg,
            exponent: -190,
            mantissa: 0xb821511c_23032258_ab53bcf6_b768e760_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Neg,
            exponent: -192,
            mantissa: 0xd05a33f0_a8898f65_6b9f85fc_7f778221_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Pos,
            exponent: -199,
            mantissa: 0xf4f29f81_a7903183_d5af6007_8369ac18_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Pos,
            exponent: -201,
            mantissa: 0xe46b6e3b_155fc435_a61266bb_793dcf42_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Neg,
            exponent: -207,
            mantissa: 0x85b474ff_0a86a6e7_ba18a86f_7c7774b1_u128,
        },
    ],
    [
        DyadicFloat128 {
            sign: DyadicSign::Pos,
            exponent: -238,
            mantissa: 0x81b4fd38_dfb221b9_f4f66331_221c7f93_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Pos,
            exponent: -131,
            mantissa: 0xe4b3805f_81cf0e62_0e0f072d_ef844dc5_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Neg,
            exponent: -137,
            mantissa: 0x8f605aef_25f7ef6b_3a409697_49d947c8_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Neg,
            exponent: -133,
            mantissa: 0x984ab9a5_1224b163_97704521_ddbab4d7_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Pos,
            exponent: -140,
            mantissa: 0xbeba7091_42e58e1b_c18b8173_6484d106_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Pos,
            exponent: -138,
            mantissa: 0xf31b7f12_c9a9607f_e5e939fa_0238b016_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Neg,
            exponent: -144,
            mantissa: 0x97f08cb9_5e062ebe_ac4978b7_16c7413f_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Neg,
            exponent: -143,
            mantissa: 0xb89726a2_62539ecc_2786f7f5_4176ece1_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Pos,
            exponent: -150,
            mantissa: 0xe61d6d5d_70228f07_617d6fcd_b0b755bb_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Pos,
            exponent: -149,
            mantissa: 0xa35696ce_5f5ae1bb_2db97d78_03d4ea86_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Neg,
            exponent: -156,
            mantissa: 0xcaec0d74_1fb1065b_dcb6b9d9_eabc9529_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Neg,
            exponent: -156,
            mantissa: 0xbd01577f_51eb0142_5671ce47_34d58c87_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Pos,
            exponent: -163,
            mantissa: 0xe9d4a4d3_911ce4ae_683f14d3_7525e24c_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Pos,
            exponent: -163,
            mantissa: 0x9a0f05d8_c9da44ec_6378cc70_9b9cd3ce_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Neg,
            exponent: -170,
            mantissa: 0xbdaa8527_bba4b74f_871c88f6_c640fb25_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Neg,
            exponent: -171,
            mantissa: 0xba632a51_32b55f35_6dc671bb_05ae2b52_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Pos,
            exponent: -178,
            mantissa: 0xe430b9b5_165b336e_e6e0c601_24b06690_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Pos,
            exponent: -179,
            mantissa: 0xadf16da3_3e5938ff_e5a6b62a_1be1e048_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Neg,
            exponent: -186,
            mantissa: 0xd3a19810_dc0d4179_e659806e_4b2e192e_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Neg,
            exponent: -187,
            mantissa: 0x80ff480d_009675b6_b765eaf9_4c31c831_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Pos,
            exponent: -194,
            mantissa: 0x9bdfc574_4297c1ca_56d7686e_0591a74b_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Pos,
            exponent: -196,
            mantissa: 0x9bade85d_6c8d089f_0421ab21_e851b6b4_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Neg,
            exponent: -203,
            mantissa: 0xbab85634_82acff82_e7bd120b_2aacff16_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Neg,
            exponent: -205,
            mantissa: 0x9bd76cda_1c4f2e5c_a18e94c5_e0b2a67a_u128,
        },
    ],
    [
        DyadicFloat128 {
            sign: DyadicSign::Pos,
            exponent: -131,
            mantissa: 0xe14d7420_05190070_9e7fd18b_0f00dfb7_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Neg,
            exponent: -241,
            mantissa: 0xc918cc36_2ae1ce1f_1caa2e60_4dadca15_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Neg,
            exponent: -132,
            mantissa: 0xe1389c6d_18e15af5_9379e688_61955b54_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Pos,
            exponent: -139,
            mantissa: 0xb68a0fe5_f542970d_2a736a7f_652c4fb9_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Pos,
            exponent: -135,
            mantissa: 0x95ee3883_1566b256_e49faeb6_2a480c84_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Neg,
            exponent: -142,
            mantissa: 0x918e7b21_ffe0b414_1ab420e0_b2fa3d03_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Neg,
            exponent: -140,
            mantissa: 0x9f85f9d9_ff45eff0_d60558e2_0301cd4a_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Pos,
            exponent: -147,
            mantissa: 0xa582f4aa_ef0ba630_1d497c26_a9bf5e7d_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Pos,
            exponent: -146,
            mantissa: 0xb5a95c34_7d2f13a3_3ea6ff9f_13cfaeff_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Neg,
            exponent: -153,
            mantissa: 0xc2d3d51f_42661452_b2695af9_e35d137b_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Neg,
            exponent: -152,
            mantissa: 0x809699c2_648b1e2e_a2b89aec_adbed5c7_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Pos,
            exponent: -159,
            mantissa: 0x8c7c6243_92a27836_cbbe5c9c_865f67bd_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Pos,
            exponent: -160,
            mantissa: 0xf7fe2014_f9782dc5_f3efeb5c_2c62b965_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Neg,
            exponent: -166,
            mantissa: 0x88ead1ed_7b87cd6c_0690d28f_1406ba48_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Neg,
            exponent: -167,
            mantissa: 0xad46c267_095196f1_a037b342_820775b6_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Pos,
            exponent: -174,
            mantissa: 0xc06fd13f_909da358_1ee4f6a1_4bc17d87_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Pos,
            exponent: -175,
            mantissa: 0xb77603f8_f4ff6621_457401ba_43f999e8_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Neg,
            exponent: -182,
            mantissa: 0xcc3eb79d_b5cbb521_7e1b1ff8_a5ad82c8_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Neg,
            exponent: -183,
            mantissa: 0x9838cb4b_1e92c90c_69dbdeb6_cd6b6ad4_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Pos,
            exponent: -190,
            mantissa: 0xa978779d_2434b68f_33db5a00_396ce7e2_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Pos,
            exponent: -192,
            mantissa: 0xcb41abfe_d3423ff1_ce0346ce_11639b79_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Neg,
            exponent: -199,
            mantissa: 0xe1dfc83b_21222e9d_d8f2a389_6044728c_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Neg,
            exponent: -201,
            mantissa: 0xdf12e5b4_8755dae0_e9c7ca30_9f4a6cb0_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Pos,
            exponent: -208,
            mantissa: 0xf713ebf7_9ffdff0b_30f02d23_86c67af6_u128,
        },
    ],
    [
        DyadicFloat128 {
            sign: DyadicSign::Pos,
            exponent: -243,
            mantissa: 0xf7cc9c96_43fa98c7_720b42c8_e56338d1_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Neg,
            exponent: -131,
            mantissa: 0xddf92348_a11af4e0_0175d5cd_e96d7339_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Pos,
            exponent: -137,
            mantissa: 0x8316d771_aa08cfd4_c538e11b_02d42a6a_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Pos,
            exponent: -133,
            mantissa: 0x93d4b78b_a5b688df_7cf44b54_d24a489d_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Neg,
            exponent: -140,
            mantissa: 0xae6daef5_4e09050c_72ee30f4_18bef434_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Neg,
            exponent: -138,
            mantissa: 0xec0c1a47_cc5076e0_87f7d140_867d117c_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Pos,
            exponent: -144,
            mantissa: 0x8b0555c1_4a04e460_d82b2798_6bb9ce59_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Pos,
            exponent: -143,
            mantissa: 0xb34c6b77_4cf5c440_3425faf1_1282d2e9_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Neg,
            exponent: -150,
            mantissa: 0xd2b19ec7_859df8ed_98199b8d_0022cc3f_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Neg,
            exponent: -149,
            mantissa: 0x9ebc46b0_8aa535d6_e89983bc_e8386244_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Pos,
            exponent: -156,
            mantissa: 0xb9f61570_0425c640_4441c204_fe194251_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Pos,
            exponent: -156,
            mantissa: 0xb7cac7b1_1aef5806_ec2848c1_cd0e93c6_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Neg,
            exponent: -163,
            mantissa: 0xd6849ce2_1d4fef4a_d2660334_ec1194a2_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Neg,
            exponent: -163,
            mantissa: 0x95eab184_e78ca925_7e147ffe_a33866b0_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Pos,
            exponent: -170,
            mantissa: 0xae384eb8_bf971db6_ee2add39_11397ec1_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Pos,
            exponent: -171,
            mantissa: 0xb585c7f3_9ba4acd8_b21b12ed_a6cfe6b6_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Neg,
            exponent: -178,
            mantissa: 0xd1e76737_574e8d66_ddb38b2a_571ccdff_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Neg,
            exponent: -179,
            mantissa: 0xa98d8abe_6443ab55_8103825b_23344840_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Pos,
            exponent: -186,
            mantissa: 0xc2f9a4ad_2ab6c584_8cb6877e_2695ecfd_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Pos,
            exponent: -188,
            mantissa: 0xfbb8b111_4c330c15_0380e8a6_dabe119f_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Neg,
            exponent: -194,
            mantissa: 0x8fd8f4ab_570da3b4_4ee1a4bc_4eebe6ab_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Neg,
            exponent: -196,
            mantissa: 0x980bc932_59f1a692_17ba1d5c_df8d458a_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Pos,
            exponent: -203,
            mantissa: 0xac9ec599_8cd73452_ce1831f8_8654be3c_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Pos,
            exponent: -205,
            mantissa: 0x985ccdc3_e1dcb711_82224cd7_f1953918_u128,
        },
    ],
    [
        DyadicFloat128 {
            sign: DyadicSign::Neg,
            exponent: -131,
            mantissa: 0xdadc139d_ba82a9b1_bd5a08c4_697da595_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Neg,
            exponent: -239,
            mantissa: 0x983206c2_1d9f20c8_03ea78e4_f1725978_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Pos,
            exponent: -132,
            mantissa: 0xdaca0c61_b5341f9b_e671b04d_f4bb3f95_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Neg,
            exponent: -139,
            mantissa: 0xa757f554_9553c567_236b693a_be884c1f_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Neg,
            exponent: -135,
            mantissa: 0x91ac00d7_cdf223d1_848c7f04_8eef8f4e_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Pos,
            exponent: -142,
            mantissa: 0x857cba3e_85f491bc_f2d21045_06724b89_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Pos,
            exponent: -140,
            mantissa: 0x9b08ed87_0ab37705_fcbb61b5_a9b4773c_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Neg,
            exponent: -147,
            mantissa: 0x97df20a6_d09283f4_529c74d9_eb8873e1_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Neg,
            exponent: -146,
            mantissa: 0xb09e65dc_e5f2caaf_383bb021_fe2f2c3c_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Pos,
            exponent: -153,
            mantissa: 0xb2e793a7_4269848b_9c0c21b2_683b3c02_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Pos,
            exponent: -153,
            mantissa: 0xfa294cb7_1ee1566d_ec113c9b_a2a3c681_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Neg,
            exponent: -159,
            mantissa: 0x811fa9c1_6d1d6591_cd1f0b29_8ca34f88_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Neg,
            exponent: -160,
            mantissa: 0xf15f127d_d6b1f585_8be75085_667e66ce_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Pos,
            exponent: -167,
            mantissa: 0xfbf69ba0_7c67abc1_b3434849_e9bfc6b7_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Pos,
            exponent: -167,
            mantissa: 0xa8c3f4cc_50ccb4d5_76fb7ece_e83b6510_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Neg,
            exponent: -174,
            mantissa: 0xb149e4f3_c0dbb307_a513f561_90476c5f_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Neg,
            exponent: -175,
            mantissa: 0xb2d24bc3_6228a4db_ea37271e_2f8358fd_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Pos,
            exponent: -182,
            mantissa: 0xbc6df811_9a114db9_d5628f2c_28b14e72_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Pos,
            exponent: -183,
            mantissa: 0x947ef163_86edb7c9_195e113b_5ba6eb19_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Neg,
            exponent: -190,
            mantissa: 0x9c95c7b0_b38a9e57_95860e6c_4a5009d4_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Neg,
            exponent: -192,
            mantissa: 0xc6752cfa_251f1274_7f903f46_f71a5945_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Pos,
            exponent: -199,
            mantissa: 0xd109edd6_54931825_0ac145a1_92255f16_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Pos,
            exponent: -201,
            mantissa: 0xda02bb3b_2e046ffe_8d4aa7a1_2206666e_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Neg,
            exponent: -208,
            mantissa: 0xe50d9c93_251a5bdc_aa60ea19_9ce7f5a3_u128,
        },
    ],
    [
        DyadicFloat128 {
            sign: DyadicSign::Neg,
            exponent: -239,
            mantissa: 0xeb3d8c00_66154c1d_a116bcaf_184ab621_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Pos,
            exponent: -131,
            mantissa: 0xd7ce5a4e_27c9a5fb_79f7b70a_4b9c7a5f_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Neg,
            exponent: -138,
            mantissa: 0xf0ec985a_89913ebc_98a27e5c_6f3cf2e4_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Neg,
            exponent: -133,
            mantissa: 0x8fbd47f6_375d416b_1ae60058_4f21870d_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Pos,
            exponent: -140,
            mantissa: 0xa052a988_4337dc6a_a0934235_324b93e0_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Pos,
            exponent: -138,
            mantissa: 0xe5906409_7e9d4cd8_299ce072_dcbd6ec7_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Neg,
            exponent: -145,
            mantissa: 0xffa8ca9e_35b65c87_dc12a692_a25bf52f_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Neg,
            exponent: -143,
            mantissa: 0xae6e247a_ebaf571a_c94f90fc_4ec4e386_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Pos,
            exponent: -150,
            mantissa: 0xc1d848ca_9639061a_38670a8f_d3c9a3fa_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Pos,
            exponent: -149,
            mantissa: 0x9a7db42f_2c55bb25_e709c81e_c1d023e9_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Neg,
            exponent: -156,
            mantissa: 0xab37e21b_4d5b8636_5014bdfe_8de1d42b_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Neg,
            exponent: -156,
            mantissa: 0xb2f87bef_14b8fcf5_549c1315_c68851ee_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Pos,
            exponent: -163,
            mantissa: 0xc5b0f2fc_dbd0020b_2edf6e64_833cb7ec_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Pos,
            exponent: -163,
            mantissa: 0x92129834_dbc547eb_5767227d_bf7ffb0c_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Neg,
            exponent: -170,
            mantissa: 0xa0b996b9_15a10f92_9b5894e5_40405ab3_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Neg,
            exponent: -171,
            mantissa: 0xb0fd4c2d_0b4932cc_5615c45a_bde0607c_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Pos,
            exponent: -178,
            mantissa: 0xc1e0bab3_982cf7c8_04e46b92_9ab3fa37_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Pos,
            exponent: -179,
            mantissa: 0xa57194ba_6f8e2d20_67f1d6f4_a8ca66e6_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Neg,
            exponent: -186,
            mantissa: 0xb453fb9f_adf3e15a_8f989475_dc7b921d_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Neg,
            exponent: -188,
            mantissa: 0xf5d23572_2c5e623a_e401e1dc_8e05b2f4_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Pos,
            exponent: -194,
            mantissa: 0x853b831b_6ee404f2_b3951af3_ccbd6904_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Pos,
            exponent: -196,
            mantissa: 0x949c46cf_e5d142dc_2dbd9a0d_7352bf3a_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Neg,
            exponent: -203,
            mantissa: 0xa0206f6d_bdf5f352_24a64205_348716cf_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Neg,
            exponent: -205,
            mantissa: 0x950de419_ad98a608_87c3a7fd_99b63821_u128,
        },
    ],
    [
        DyadicFloat128 {
            sign: DyadicSign::Pos,
            exponent: -131,
            mantissa: 0xd4f09d06_db214bff_59e7c483_1eee3ccb_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Pos,
            exponent: -240,
            mantissa: 0x9a9844f4_c88b658d_de4fa4c9_11bb271d_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Neg,
            exponent: -132,
            mantissa: 0xd4e0e510_b0dcd905_f5aef327_44388e70_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Pos,
            exponent: -139,
            mantissa: 0x9a2504f7_ee31a54a_8cea6d00_a2b4a645_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Pos,
            exponent: -135,
            mantissa: 0x8dc16218_511c0862_fc3955b3_d2e0b091_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Neg,
            exponent: -143,
            mantissa: 0xf5fdcbde_500f9275_328aa108_e675e7b4_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Neg,
            exponent: -140,
            mantissa: 0x96e6d8de_0e81c064_4da87b5e_d7e50d23_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Pos,
            exponent: -147,
            mantissa: 0x8c0023f8_5bcd73d4_61d98c24_e2ab09c9_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Pos,
            exponent: -146,
            mantissa: 0xabf77774_268680b5_703a716c_f111d4a0_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Neg,
            exponent: -153,
            mantissa: 0xa50634ce_418c58d8_ca0dd870_5e067264_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Neg,
            exponent: -153,
            mantissa: 0xf3aca745_163aff1e_946545f0_110933f9_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Pos,
            exponent: -160,
            mantissa: 0xee65dda6_05d3f443_3bbbc5fe_26f83e93_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Pos,
            exponent: -160,
            mantissa: 0xeb3b319b_dae4ab9f_08fdcce6_86fd48fe_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Neg,
            exponent: -167,
            mantissa: 0xe8d02e75_0c501120_0bd41d3e_51f40245_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Neg,
            exponent: -167,
            mantissa: 0xa4917da8_514f8be6_c45916ed_ebd2f2c5_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Pos,
            exponent: -174,
            mantissa: 0xa3fd1adf_b13e6ebe_96a650e7_91aae11c_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Pos,
            exponent: -175,
            mantissa: 0xae7d0438_9114a137_29747e80_12480a85_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Neg,
            exponent: -182,
            mantissa: 0xae804c9b_dfad588a_21817ded_5082e637_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Neg,
            exponent: -183,
            mantissa: 0x9100580f_edfce685_74d2061a_53d45c3f_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Pos,
            exponent: -190,
            mantissa: 0x9132f802_60ea45eb_034d4bc0_971d6e66_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Pos,
            exponent: -192,
            mantissa: 0xc1efb8d7_455bd1e5_2805deaa_b210c3d8_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Neg,
            exponent: -199,
            mantissa: 0xc21c0000_0eefaff2_732b1072_afac0f23_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Neg,
            exponent: -201,
            mantissa: 0xd5378717_928a01b3_ebd0a111_56ec9ffc_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Pos,
            exponent: -208,
            mantissa: 0xd50243c6_28ad4980_a266a5c0_6bb15491_u128,
        },
    ],
    [
        DyadicFloat128 {
            sign: DyadicSign::Pos,
            exponent: -240,
            mantissa: 0xdec2262a_368ff1a1_2f346627_67b224ec_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Neg,
            exponent: -131,
            mantissa: 0xd2203f02_14c68aba_6b511e3c_145ee51c_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Pos,
            exponent: -138,
            mantissa: 0xde64ef86_780591f6_dddc07ad_d4be4e14_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Pos,
            exponent: -133,
            mantissa: 0x8bf81336_d90a4478_68e6678c_61d6e67a_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Neg,
            exponent: -140,
            mantissa: 0x94050291_1a502761_345a979e_f01185ef_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Neg,
            exponent: -138,
            mantissa: 0xdf956d0f_a203879b_f5265bbf_1b225221_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Pos,
            exponent: -145,
            mantissa: 0xec1eac01_95a466b5_2c05eaf6_e8ed103c_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Pos,
            exponent: -143,
            mantissa: 0xa9eecb96_c791900e_b77e8223_0935166c_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Neg,
            exponent: -150,
            mantissa: 0xb31e0ff1_01014261_8236c00d_122fda34_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Neg,
            exponent: -149,
            mantissa: 0x968fdf91_01ffb16d_bba4ce8b_d7d0e47d_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Pos,
            exponent: -156,
            mantissa: 0x9e4f85ba_ce2fcb0c_dcdc6916_f6e69b6e_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Pos,
            exponent: -156,
            mantissa: 0xae7f0d80_934b6a91_3fe67bf9_d1e6315b_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Neg,
            exponent: -163,
            mantissa: 0xb6ee03ca_bc0f4667_ecc85386_9a5d278c_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Neg,
            exponent: -163,
            mantissa: 0x8e7ea4a8_213351b1_3e076398_cc2c2abd_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Pos,
            exponent: -170,
            mantissa: 0x94dbd7a9_0a7de3b6_ecb48b9f_5a99761e_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Pos,
            exponent: -171,
            mantissa: 0xacc18155_9d6450ea_ca59f6f0_a2559bfc_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Neg,
            exponent: -178,
            mantissa: 0xb3bfb530_255a1286_de178d32_816fd258_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Neg,
            exponent: -179,
            mantissa: 0xa1976609_15de1b55_4fa4082b_a366d337_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Pos,
            exponent: -186,
            mantissa: 0xa7608a5f_db5fb85b_089f78e2_491291fc_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Pos,
            exponent: -188,
            mantissa: 0xf0443fc6_17abe988_22ee44d5_5d6130f4_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Neg,
            exponent: -195,
            mantissa: 0xf7a251bf_c2788d2c_638307bd_c0a279f8_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Neg,
            exponent: -196,
            mantissa: 0x915c8faf_a224d457_61d62c1a_d6853e46_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Pos,
            exponent: -203,
            mantissa: 0x9501fc01_80d50e91_ccecff8e_c4c5e08e_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Pos,
            exponent: -205,
            mantissa: 0x91e91fc7_bd697471_c1d8764c_596534b6_u128,
        },
    ],
    [
        DyadicFloat128 {
            sign: DyadicSign::Neg,
            exponent: -131,
            mantissa: 0xcf79dd91_09d856dc_bfbf7651_5fecb56d_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Pos,
            exponent: -240,
            mantissa: 0xe3d77c11_bdb690a5_999d70a8_2000423a_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Pos,
            exponent: -132,
            mantissa: 0xcf6c1003_d4d22942_c5cb51f2_b57616d0_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Neg,
            exponent: -139,
            mantissa: 0x8e97ef75_94fc5075_06dd6e2a_37cae0fc_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Neg,
            exponent: -135,
            mantissa: 0x8a234337_f8dc14ec_d74d5a02_5c34d435_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Pos,
            exponent: -143,
            mantissa: 0xe39db427_c992ca97_e1feec37_59c78b10_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Pos,
            exponent: -140,
            mantissa: 0x93146f2c_a5b6450c_2408fe80_745d0eb5_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Neg,
            exponent: -147,
            mantissa: 0x81981e0a_13a2486f_be4f199b_d2e6266f_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Neg,
            exponent: -146,
            mantissa: 0xa7a8803e_288434dc_d62652aa_75991cba_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Pos,
            exponent: -153,
            mantissa: 0x98d6fd1b_79109aa7_a3f4428f_dc6604a5_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Pos,
            exponent: -153,
            mantissa: 0xeda79835_8c84a481_ff41f406_6d0ecad9_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Neg,
            exponent: -160,
            mantissa: 0xdcf19acc_7c4cece1_1ad6c5d8_cece9434_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Neg,
            exponent: -160,
            mantissa: 0xe584f6b7_65c4b28a_f864e9d7_ab235393_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Pos,
            exponent: -167,
            mantissa: 0xd7f08c18_a454c54c_121b60ec_f23678ea_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Pos,
            exponent: -167,
            mantissa: 0xa0a72215_546aba37_f5f3f797_2be32eb5_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Neg,
            exponent: -174,
            mantissa: 0x983e094a_0c1cb438_eb154bfc_c4e5f8e4_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Neg,
            exponent: -175,
            mantissa: 0xaa6ed5e5_3db94bbc_caece01a_fc82e765_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Pos,
            exponent: -182,
            mantissa: 0xa22ae311_c0b3c60d_70b0a6e2_c0bf33f9_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Pos,
            exponent: -183,
            mantissa: 0x8db81326_d1c16581_c6ce2774_8d19ccfa_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Neg,
            exponent: -190,
            mantissa: 0x87167d4b_1ea62b02_5379e3ac_73f54dd9_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Neg,
            exponent: -192,
            mantissa: 0xbdac4b6f_1a2bdaf3_903cb72c_154285f6_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Pos,
            exponent: -199,
            mantissa: 0xb4cf7e52_e6c63ecf_de83e087_4c0799bd_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Pos,
            exponent: -201,
            mantissa: 0xd0ad70c4_676e310b_0e37cbd9_4399e7f7_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Neg,
            exponent: -208,
            mantissa: 0xc6ab9057_68353371_32e26d45_983cf0c3_u128,
        },
    ],
    [
        DyadicFloat128 {
            sign: DyadicSign::Neg,
            exponent: -239,
            mantissa: 0xc658525d_c958e344_2f9ab4bf_8b185392_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Pos,
            exponent: -131,
            mantissa: 0xccdf3a20_0c6482e5_15726944_d3f8a12a_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Neg,
            exponent: -138,
            mantissa: 0xce1fa951_57b921b8_852614fe_8587522f_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Neg,
            exponent: -133,
            mantissa: 0x887ae528_fcc581f2_81d89ca1_2e3a1dd2_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Pos,
            exponent: -140,
            mantissa: 0x89364746_513b2d87_32fca78c_496a300f_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Pos,
            exponent: -138,
            mantissa: 0xda0b7bea_7e189837_68e3656d_9e7e6f8b_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Neg,
            exponent: -145,
            mantissa: 0xdaf172ac_070e7177_b33eb16c_ea47359e_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Neg,
            exponent: -143,
            mantissa: 0xa5c3113c_ef3ab311_a8cc6eed_ee7df032_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Pos,
            exponent: -150,
            mantissa: 0xa62838d9_cee09497_877a101f_68532b74_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Pos,
            exponent: -149,
            mantissa: 0x92e97e1d_8c810443_3b3642ff_2518265b_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Neg,
            exponent: -156,
            mantissa: 0x92ef9911_3ccb8cf0_d8d844f8_1afc61df_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Neg,
            exponent: -156,
            mantissa: 0xaa54b8f5_b7f1da7d_8475add2_e1376223_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Pos,
            exponent: -163,
            mantissa: 0xa9e6364c_91795a09_02136e97_b1478395_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Pos,
            exponent: -163,
            mantissa: 0x8b27ca8c_06c606c9_ef13a6aa_b91a68a3_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Neg,
            exponent: -170,
            mantissa: 0x8a5cf7a0_276c45b9_338eb98d_5126b640_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Neg,
            exponent: -171,
            mantissa: 0xa8cb11da_72855d43_9ab2d11a_4918d713_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Pos,
            exponent: -178,
            mantissa: 0xa73938d3_8f91e358_584ac712_2c6575ee_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Pos,
            exponent: -179,
            mantissa: 0x9df950dc_2c8cbb27_c46e8514_d1471da5_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Neg,
            exponent: -186,
            mantissa: 0x9bde06c7_225c1ab5_3e948215_df383e52_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Neg,
            exponent: -188,
            mantissa: 0xeb0825ca_623ebf02_bdda7d9b_0470229a_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Pos,
            exponent: -195,
            mantissa: 0xe6da59aa_25e3fc84_d1399aae_a5765eef_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Pos,
            exponent: -196,
            mantissa: 0x8e49ab3a_3155cf84_fc3e7b98_a43ec863_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Neg,
            exponent: -203,
            mantissa: 0x8b11f794_9c03737b_a42a4fc8_7b61bb7f_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Neg,
            exponent: -205,
            mantissa: 0x8eec8bde_4299354f_f2f7dcae_a4ba58f4_u128,
        },
    ],
    [
        DyadicFloat128 {
            sign: DyadicSign::Pos,
            exponent: -131,
            mantissa: 0xca6993b4_8a72857e_d226f390_3f871be2_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Neg,
            exponent: -239,
            mantissa: 0xbe269d08_48e32334_a5cad7c8_a94dd482_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Neg,
            exponent: -132,
            mantissa: 0xca5d60ca_8718c48a_23c527e5_8ea7a528_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Pos,
            exponent: -139,
            mantissa: 0x846a7fc7_902cf37f_c75aa163_83160f27_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Pos,
            exponent: -135,
            mantissa: 0x86c8694d_ff766466_bc7b11e9_bb5d9b22_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Neg,
            exponent: -143,
            mantissa: 0xd36ab0b1_7bfda30f_fdcc44be_8233fc08_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Neg,
            exponent: -140,
            mantissa: 0x8f883ed1_1c03eafc_e466277f_0c0bb945_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Pos,
            exponent: -148,
            mantissa: 0xf0d31313_ba1bda59_7bbd61f5_9435fca0_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Pos,
            exponent: -146,
            mantissa: 0xa3a757e8_663212cb_f5a6e867_6f1de139_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Neg,
            exponent: -153,
            mantissa: 0x8e138a27_a16bbf91_b955775e_491a73cf_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Neg,
            exponent: -153,
            mantissa: 0xe80cccd8_d27f451c_ec540bd0_240ab349_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Pos,
            exponent: -160,
            mantissa: 0xcd80a400_ed43d015_f0e0348c_776a7d2f_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Pos,
            exponent: -160,
            mantissa: 0xe030ba05_e4522b2b_170819b5_0143955d_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Neg,
            exponent: -167,
            mantissa: 0xc8fc614b_92a3ed2f_151ee904_ee957211_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Neg,
            exponent: -167,
            mantissa: 0x9cfdaac5_bd52647f_aec8d1fe_8f89bfa5_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Pos,
            exponent: -174,
            mantissa: 0x8dcfba6a_3d9afb3a_d72b6b9c_a67d7b8e_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Pos,
            exponent: -175,
            mantissa: 0xa6a12aab_f0f8e259_e6421d02_7add91b2_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Neg,
            exponent: -182,
            mantissa: 0x9730bd71_3ca96c3c_9c8d95b2_bb3d2141_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Neg,
            exponent: -183,
            mantissa: 0x8aa1944c_614e7d8a_98659cdb_e721bc79_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Pos,
            exponent: -191,
            mantissa: 0xfc21fb1a_3cb1fa9f_d9ddf7e9_a73d3cce_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Pos,
            exponent: -192,
            mantissa: 0xb9a60987_e324735f_8632945b_60e3f704_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Neg,
            exponent: -199,
            mantissa: 0xa8e9caca_4d7de3d3_b527dd1b_07370901_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Neg,
            exponent: -201,
            mantissa: 0xcc607818_14780770_2404a8ba_81c35788_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Pos,
            exponent: -208,
            mantissa: 0xb9ce6234_252de047_c57503df_06742ec6_u128,
        },
    ],
    [
        DyadicFloat128 {
            sign: DyadicSign::Pos,
            exponent: -242,
            mantissa: 0xc2331c68_792b24e3_262c4446_73c3f141_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Neg,
            exponent: -131,
            mantissa: 0xc7fe4de9_27f03c3b_d8bcb699_31da9fde_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Pos,
            exponent: -138,
            mantissa: 0xbfbea355_99b90abd_14bae19e_aaac8265_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Pos,
            exponent: -133,
            mantissa: 0x853d392e_9fde21f5_32ddb23d_ab3c3d61_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Neg,
            exponent: -141,
            mantissa: 0xff50b947_90b5a27e_113551e5_5bd016f3_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Neg,
            exponent: -138,
            mantissa: 0xd4e5634f_876b21f6_3ef33193_27c6cd3e_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Pos,
            exponent: -145,
            mantissa: 0xcbbfa47f_882b5b1e_eeb9f0dd_38bc3b24_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Pos,
            exponent: -143,
            mantissa: 0xa1e16c73_4c6081a7_7447ea6f_e43604ee_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Neg,
            exponent: -150,
            mantissa: 0x9aae8d33_b2f52cd8_142dd54d_31c97b85_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Neg,
            exponent: -149,
            mantissa: 0x8f82a9a8_0e890744_5ff1862a_f41f3d55_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Pos,
            exponent: -156,
            mantissa: 0x88da3fb5_224147db_8b82fdfc_8c5d6d6d_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Pos,
            exponent: -156,
            mantissa: 0xa6711981_c44d474a_fe8002ea_50a8c447_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Neg,
            exponent: -163,
            mantissa: 0x9e54c6fb_1110f92b_03ad2eec_1743f215_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Neg,
            exponent: -163,
            mantissa: 0x8807e3a3_294eba15_fb98ac77_cf327d2c_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Pos,
            exponent: -170,
            mantissa: 0x810788f1_6e503c46_83f45de3_dcc29947_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Pos,
            exponent: -171,
            mantissa: 0xa5137791_0ce15cac_3316e87e_82baf0d9_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Neg,
            exponent: -178,
            mantissa: 0x9c101ffa_963a3686_eb3d6064_6eeb10f0_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Neg,
            exponent: -179,
            mantissa: 0x9a92260a_9d514944_1377b47c_fd951628_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Pos,
            exponent: -186,
            mantissa: 0x9196dddf_51afd41d_18f4d5de_97ab28b3_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Pos,
            exponent: -188,
            mantissa: 0xe617986b_f254317c_43340294_39008bb3_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Neg,
            exponent: -195,
            mantissa: 0xd7d51170_21c3d219_bd46ad4d_bcb827a6_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Neg,
            exponent: -196,
            mantissa: 0x8b60a105_0135dd54_6c56c8bd_73ae2652_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Pos,
            exponent: -203,
            mantissa: 0x82270bdb_83cbf1e3_6d7bd132_71591375_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Pos,
            exponent: -205,
            mantissa: 0x8c1601b7_b3722389_174cb7c2_c3fc5a68_u128,
        },
    ],
    [
        DyadicFloat128 {
            sign: DyadicSign::Neg,
            exponent: -131,
            mantissa: 0xc5b3d124_083be6ed_c4d0f183_62325eac_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Neg,
            exponent: -239,
            mantissa: 0x9f94cfad_f2be101b_c064c937_f38c3254_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Pos,
            exponent: -132,
            mantissa: 0xc5a8f908_34750c96_6adc833e_9835c32a_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Neg,
            exponent: -140,
            mantissa: 0xf6c9ae3d_77d186da_114f387f_2dc44c34_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Neg,
            exponent: -135,
            mantissa: 0x83a91528_1bf9688f_4bb14d27_de949662_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Pos,
            exponent: -143,
            mantissa: 0xc50cb030_fc5de354_056526dc_1abe1482_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Pos,
            exponent: -140,
            mantissa: 0x8c3a523f_1e0f1aa0_7bb97206_9bedfff7_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Neg,
            exponent: -148,
            mantissa: 0xe08679a6_0278b755_0367c2ed_0ae8d09c_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Neg,
            exponent: -146,
            mantissa: 0x9feb616a_7d1521b1_d5a1969a_1af33332_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Pos,
            exponent: -153,
            mantissa: 0x84836aaa_913c2646_cd3a6b82_9a392c6f_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Pos,
            exponent: -153,
            mantissa: 0xe2d0e092_6fd03fc9_4a642b8f_d23cdb1a_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Neg,
            exponent: -160,
            mantissa: 0xbfc41d98_d2fac858_7f2e6cfb_455834c0_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Neg,
            exponent: -160,
            mantissa: 0xdb3467a9_81472155_31c808ef_bb7d4129_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Pos,
            exponent: -167,
            mantissa: 0xbba9b968_2d92ea08_ec3d93de_67e405ef_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Pos,
            exponent: -167,
            mantissa: 0x998ec38e_b3befe0d_3057a5c9_c3ac0715_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Neg,
            exponent: -174,
            mantissa: 0x84808124_4c4e9aca_7344685c_75eecf2a_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Neg,
            exponent: -175,
            mantissa: 0xa30e2033_e3a16659_97d8d254_f49dc89f_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Pos,
            exponent: -182,
            mantissa: 0x8d5fcc9b_4b73bc50_20d5acdc_40bd6ead_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Pos,
            exponent: -183,
            mantissa: 0x87b8afd1_c1dc5261_226ab6fb_d316193e_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Neg,
            exponent: -191,
            mantissa: 0xebf69798_6d81bf52_ea1a35ff_e2a69712_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Neg,
            exponent: -192,
            mantissa: 0xb5d85a6d_3eb5bc6a_0019c93c_1cf93179_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Pos,
            exponent: -199,
            mantissa: 0x9e39f870_abf2df96_33aef844_900a2469_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Pos,
            exponent: -201,
            mantissa: 0xc84ca2e1_caca358e_6e6981d9_e0ff6618_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Neg,
            exponent: -208,
            mantissa: 0xae38e61f_0b6fecda_4da2a7bb_e0b54be1_u128,
        },
    ],
    [
        DyadicFloat128 {
            sign: DyadicSign::Neg,
            exponent: -240,
            mantissa: 0xf0ed1a5f_f5c6fedf_45847fbf_5af9f05b_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Pos,
            exponent: -131,
            mantissa: 0xc3728df0_548a9b63_930432e8_28fdb92b_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Neg,
            exponent: -138,
            mantissa: 0xb2f68d9c_3f2ed354_92870bb5_64ccaefd_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Neg,
            exponent: -133,
            mantissa: 0x8237e2d1_0c542f98_c0fb2701_57c779a0_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Pos,
            exponent: -141,
            mantissa: 0xee530bff_7d489be2_ff8ab511_72c3ff00_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Pos,
            exponent: -138,
            mantissa: 0xd0180110_a945357d_c44581bb_651fde74_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Neg,
            exponent: -145,
            mantissa: 0xbe3b0c48_378f102b_245aafa6_9e47b75a_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Neg,
            exponent: -143,
            mantissa: 0x9e41c44b_e850ecf4_790b8065_a6c2d959_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Pos,
            exponent: -150,
            mantissa: 0x9076f4fc_eebee630_361610f8_ab7e6d14_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Pos,
            exponent: -149,
            mantissa: 0x8c54a0b2_ac95dba6_715228c3_92562d91_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Neg,
            exponent: -157,
            mantissa: 0xffbb1028_8c7844dc_81a7ee18_310aa94b_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Neg,
            exponent: -156,
            mantissa: 0xa2cceff7_780d06a3_03df80cd_3297ff83_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Pos,
            exponent: -163,
            mantissa: 0x9401f2cc_c189a10b_7bcb17c1_41a89f3c_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Pos,
            exponent: -163,
            mantissa: 0x85199017_f58217c8_e6c0054d_f56a9bf0_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Neg,
            exponent: -171,
            mantissa: 0xf15ffb9c_3407bf22_6c231fd9_824f96ef_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Neg,
            exponent: -171,
            mantissa: 0xa194e799_5e5b735e_2ad75941_5dcdd8f9_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Pos,
            exponent: -178,
            mantissa: 0x92124487_db67d655_66aab46b_6a9b0526_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Pos,
            exponent: -179,
            mantissa: 0x975d32bb_24ecd68a_81049e0b_79ece402_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Neg,
            exponent: -186,
            mantissa: 0x885ed24d_f5f135ea_e955eb80_8e41148f_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Neg,
            exponent: -188,
            mantissa: 0xe16cb6ca_79a6b15e_1a905ee6_1d30997c_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Pos,
            exponent: -195,
            mantissa: 0xca54dc28_4eed232d_6aa14537_46a05dbd_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Pos,
            exponent: -196,
            mantissa: 0x889e902d_31df07b7_fa999ad2_ca3f6649_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Neg,
            exponent: -204,
            mantissa: 0xf43d14c8_7a4744a1_a238d3d6_bdd2f6fe_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Neg,
            exponent: -205,
            mantissa: 0x89634a34_260a5981_09d37f95_5d590c34_u128,
        },
    ],
    [
        DyadicFloat128 {
            sign: DyadicSign::Pos,
            exponent: -131,
            mantissa: 0xc14e837f_7493318c_aabfef07_d2c0833b_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Neg,
            exponent: -241,
            mantissa: 0xae4ac840_20ea03b2_bd44e8dd_5b81f5b2_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Neg,
            exponent: -132,
            mantissa: 0xc144d293_6bbc2ade_fe48a446_d0a00efd_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Pos,
            exponent: -140,
            mantissa: 0xe6b401aa_e4f5ae2a_a0a6c9de_efd3bcf8_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Pos,
            exponent: -135,
            mantissa: 0x80beb87a_89430b3b_5f5f790a_066e0987_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Neg,
            exponent: -143,
            mantissa: 0xb83cbc71_506dbd62_0db6dd54_d6407e55_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Neg,
            exponent: -140,
            mantissa: 0x8923e70a_f39c2984_f6dd5540_27ab8c2d_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Pos,
            exponent: -148,
            mantissa: 0xd1fb092a_42f609c1_549321f5_8962adf9_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Pos,
            exponent: -146,
            mantissa: 0x9c6d4242_00df456f_46d5ea1d_c3825c3c_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Neg,
            exponent: -154,
            mantissa: 0xf7f1cec6_da66c303_8fccb199_687f1999_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Neg,
            exponent: -153,
            mantissa: 0xddea078e_d874a220_ce4e2a5b_5a4510d9_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Pos,
            exponent: -160,
            mantissa: 0xb37bea4a_1b7c19e0_1f9c4e4c_0e275f51_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Pos,
            exponent: -160,
            mantissa: 0xd68740a4_114d972b_10d7c603_c3db5649_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Neg,
            exponent: -167,
            mantissa: 0xafbc2ca8_99e85877_a9aca636_72ce60f0_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Neg,
            exponent: -167,
            mantissa: 0x9654dd88_190d9600_25ed4eaf_a74d2e2e_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Pos,
            exponent: -175,
            mantissa: 0xf84f24e4_43f48544_f7cb4d60_4d0231d6_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Pos,
            exponent: -175,
            mantissa: 0x9fb0775b_37ba3589_2ca05766_0318f385_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Neg,
            exponent: -182,
            mantissa: 0x848eb44e_629e990a_7c26357c_c610e26a_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Neg,
            exponent: -183,
            mantissa: 0x84f99af4_3100cdaa_b16fcea9_9df33d40_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Pos,
            exponent: -191,
            mantissa: 0xdd69a154_703be3ae_94ae65ec_9cd980f4_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Pos,
            exponent: -192,
            mantissa: 0xb23ef438_fe574cc6_f18a86a3_4e834ae9_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Neg,
            exponent: -199,
            mantissa: 0x9497043c_c4090666_60ef818b_18e390c8_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Neg,
            exponent: -201,
            mantissa: 0xc46e185e_b21710ea_b5252925_059881e0_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Pos,
            exponent: -208,
            mantissa: 0xa3c0ff7f_661d581e_8c81645a_59c3eb0e_u128,
        },
    ],
    [
        DyadicFloat128 {
            sign: DyadicSign::Neg,
            exponent: -239,
            mantissa: 0x9e4e02fe_346eb5e3_8a88cd61_cc649ccb_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Neg,
            exponent: -131,
            mantissa: 0xbf32b77d_804d6fb9_be029a97_df59b9d2_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Pos,
            exponent: -138,
            mantissa: 0xa78a8337_9ec3b34a_70454a5b_2839a6a9_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Pos,
            exponent: -134,
            mantissa: 0xfec99613_556c9fc0_e1c6d193_d5e2308c_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Neg,
            exponent: -141,
            mantissa: 0xdf2306ec_3505214e_6bbdf1f5_0131deb3_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Neg,
            exponent: -138,
            mantissa: 0xcb99dac8_c3cfbe25_ea4c1123_05c9e0f0_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Pos,
            exponent: -145,
            mantissa: 0xb224462b_01170361_cc414086_29e91063_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Pos,
            exponent: -143,
            mantissa: 0x9add2c5f_97cc975d_3ba4e2e3_8a66ddd4_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Neg,
            exponent: -150,
            mantissa: 0x87523de7_5cd47ba0_2928ac92_544a66ca_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Neg,
            exponent: -149,
            mantissa: 0x89599382_ae8748c1_cbc3c4cd_38d7984f_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Pos,
            exponent: -157,
            mantissa: 0xefa17e27_12415e15_3493be42_ef808df7_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Pos,
            exponent: -156,
            mantissa: 0x9f61f3b4_3b6dcff2_210d2ee9_c7649beb_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Neg,
            exponent: -163,
            mantissa: 0x8ac01c9f_a39e19c4_c84b112b_e572a434_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Neg,
            exponent: -163,
            mantissa: 0x82581a9e_1ac9a7b5_67e24a71_3dfa7fdb_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Pos,
            exponent: -171,
            mantissa: 0xe265187f_f1508e69_7e516c1e_d3ac4b38_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Pos,
            exponent: -171,
            mantissa: 0x9e4a3dc9_b617e5c1_bf28e745_8bdf55c5_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Neg,
            exponent: -178,
            mantissa: 0x891639fb_b8530278_31fc0d44_5abbde3d_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Neg,
            exponent: -179,
            mantissa: 0x9456394a_8134fcdf_5d64d3d2_b661f6dc_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Pos,
            exponent: -186,
            mantissa: 0x801121c3_b5c84a9b_4d33aa7d_893e9463_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Pos,
            exponent: -188,
            mantissa: 0xdd0214a6_a92c210a_1100a850_e1a6cf9c_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Neg,
            exponent: -195,
            mantissa: 0xbe2631d2_b67aebc9_95e54187_b5715f7b_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Neg,
            exponent: -196,
            mantissa: 0x8600bc7a_56e2be2a_425840b0_4af2072b_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Pos,
            exponent: -204,
            mantissa: 0xe5b67812_d5ce81a4_59191582_3246858a_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Pos,
            exponent: -205,
            mantissa: 0x86d232dd_e6a22a54_372b42cf_980bf270_u128,
        },
    ],
    [
        DyadicFloat128 {
            sign: DyadicSign::Neg,
            exponent: -131,
            mantissa: 0xbd311903_cc0ba87e_548d0691_843dfb12_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Neg,
            exponent: -240,
            mantissa: 0x9aa210dc_a4219ce8_e719f21e_c99ad26b_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Pos,
            exponent: -132,
            mantissa: 0xbd286522_825dc3eb_55bff0d1_30cfe187_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Neg,
            exponent: -140,
            mantissa: 0xd84ae65a_85346284_bda3c5ff_70f81e7d_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Neg,
            exponent: -136,
            mantissa: 0xfc07788e_d51bcdee_e37d1dc8_49299b1d_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Pos,
            exponent: -143,
            mantissa: 0xacc116e3_ad831e4e_4b632001_3285067b_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Pos,
            exponent: -140,
            mantissa: 0x863f353e_33318244_8abb83ba_2927132d_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Neg,
            exponent: -148,
            mantissa: 0xc4f005ae_1ac67c69_08f07e5b_b9264ccc_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Neg,
            exponent: -146,
            mantissa: 0x9926a91c_59367f1d_0c56a392_e1ef0ee0_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Pos,
            exponent: -154,
            mantissa: 0xe89d4446_aee98232_bb38f883_14af9068_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Pos,
            exponent: -153,
            mantissa: 0xd94fca0d_4ad837af_f56e71d4_6b2add40_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Neg,
            exponent: -160,
            mantissa: 0xa87375d9_e8efad9d_99a05cfd_afff0cf7_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Neg,
            exponent: -160,
            mantissa: 0xd221a635_c94eeaf4_9a348e68_c309de1f_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Pos,
            exponent: -167,
            mantissa: 0xa501ff16_23c2128c_1a0bbedc_eadccef3_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Pos,
            exponent: -167,
            mantissa: 0x934b1479_6fcc6185_50d3cb8e_7ab20465_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Neg,
            exponent: -175,
            mantissa: 0xe9466047_7a13be35_0a224859_9c690315_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Neg,
            exponent: -175,
            mantissa: 0x9c838319_a55e18c3_6ce9cf94_0f2f001e_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Pos,
            exponent: -183,
            mantissa: 0xf9361f52_21f5ee33_0947bd5d_4e829739_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Pos,
            exponent: -183,
            mantissa: 0x8260e699_a224166e_ee24c3e5_9f81f87e_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Neg,
            exponent: -191,
            mantissa: 0xd044916e_75182b87_6260f632_cfd6a583_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Neg,
            exponent: -192,
            mantissa: 0xaed5dffb_615fd590_1e2171c4_80a6a9f0_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Pos,
            exponent: -199,
            mantissa: 0x8bde65b6_68169643_36c6067d_986e06ce_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Pos,
            exponent: -201,
            mantissa: 0xc0c130d1_10a70b68_16c33399_f4904abb_u128,
        },
        DyadicFloat128 {
            sign: DyadicSign::Neg,
            exponent: -208,
            mantissa: 0x9a42f7b3_0699d3f9_ae264dce_ad67c68f_u128,
        },
    ],
];
