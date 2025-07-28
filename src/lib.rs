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
// #![forbid(unsafe_code)]
#![deny(unreachable_pub)]
#![allow(
    clippy::excessive_precision,
    clippy::approx_constant,
    clippy::manual_range_contains
)]
#![deny(
    clippy::print_stdout,
    clippy::print_stderr,
    clippy::print_literal,
    clippy::print_in_format_impl
)]
mod acos;
mod acosf;
mod acospi;
mod acospif;
mod asin;
mod asin_eval_dyadic;
mod asinf;
mod asinpi;
mod asinpif;
mod bessel;
mod bits;
mod cbrt;
mod cbrtf;
mod ceil;
mod common;
mod compound;
mod compound_m1;
mod compound_m1f;
mod compoundf;
mod cosf;
mod csc;
mod cscf;
mod double_double;
mod dyadic_float;
mod dyadic_float256;
mod err;
mod exponents;
mod floor;
mod hyperbolic;
mod hypot;
mod hypotf;
mod logs;
mod polyeval;
mod pow;
mod pow_exec;
mod pow_tables;
mod powf;
mod powf_tables;
mod round;
mod round_ties_even;
mod sec;
mod secf;
mod shared_eval;
mod sin;
mod sin_helper;
mod sin_table;
mod sinc;
mod sincf;
mod sincos;
mod sincos_dyadic;
mod sincos_reduce;
mod sincos_reduce_tables;
mod sincosf;
mod sincospi;
mod sincospi_tables;
mod sincospif;
mod sinf;
mod sqrtf;
mod tangent;
mod triple_double;
mod trunc;
mod u256;

pub use acos::f_acos;
pub use acosf::f_acosf;
pub use acospi::f_acospi;
pub use acospif::f_acospif;
pub use asin::f_asin;
pub use asinf::f_asinf;
pub use asinpi::f_asinpi;
pub use asinpif::f_asinpif;
pub use bessel::{f_i0, f_i0f, f_i1, f_i1f, f_j0, f_j0f, f_j1, f_j1f, f_y0, f_y0f, f_y1, f_y1f};
pub use cbrt::f_cbrt;
pub use cbrtf::{cbrtf, f_cbrtf};
pub use ceil::{ceil, ceilf};
pub use common::{copysignfk, copysignk};
pub use compound::f_compound;
pub use compound_m1::f_compound_m1;
pub use compound_m1f::f_compound_m1f;
pub use compoundf::f_compoundf;
pub use cosf::f_cosf;
pub use csc::f_csc;
pub use cscf::f_cscf;
pub use err::{f_erf, f_erfc, f_erfcf, f_erff};
pub use exponents::{
    exp, expf, f_exp, f_exp2, f_exp2f, f_exp2m1, f_exp2m1f, f_exp10, f_exp10f, f_exp10m1,
    f_exp10m1f, f_expf, f_expm1, f_expm1f,
};
pub use floor::{floor, floorf};
pub use hyperbolic::{
    f_acosh, f_acoshf, f_asinh, f_asinhf, f_atanh, f_atanhf, f_cosh, f_coshf, f_sinh, f_sinhf,
    f_tanh, f_tanhf,
};
pub use hypot::f_hypot;
pub use hypotf::{f_hypot3f, f_hypotf};
pub use logs::{
    f_log, f_log1p, f_log1pf, f_log2, f_log2f, f_log2p1, f_log2p1f, f_log10, f_log10f, f_log10p1,
    f_log10p1f, f_logf, log, logf,
};
pub use pow::{f_pow, pow};
pub use powf::{dirty_powf, f_powf, powf};
pub use round::{round, roundf};
pub use round_ties_even::{round_ties_even, roundf_ties_even};
pub use sec::f_sec;
pub use secf::f_secf;
pub use sin::{f_cos, f_sin};
pub use sinc::f_sinc;
pub use sincf::f_sincf;
pub use sincos::f_sincos;
pub use sincosf::f_sincosf;
pub use sincospi::{f_cospi, f_sinpi};
pub use sincospif::{f_cospif, f_sinpif};
pub use sinf::f_sinf;
pub use sqrtf::sqrtf;
pub use tangent::{
    f_atan, f_atan2, f_atan2f, f_atan2pi, f_atan2pif, f_atanf, f_atanpi, f_atanpif, f_cot, f_cotf,
    f_tan, f_tanf, f_tanpi, f_tanpif,
};
pub use trunc::{trunc, truncf};
