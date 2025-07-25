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
mod acosh;
mod acoshf;
mod acospi;
mod acospif;
mod asin;
mod asin_eval_dyadic;
mod asinf;
mod asinh;
mod asinhf;
mod asinpi;
mod asinpif;
mod atan;
mod atan2;
mod atan2f;
mod atan2pi;
mod atan2pif;
mod atanf;
mod atanh;
mod atanhf;
mod atanpi;
mod atanpif;
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
mod cosh;
mod coshf;
mod cot;
mod cotf;
mod csc;
mod cscf;
mod double_double;
mod dyadic_float;
mod dyadic_float256;
mod erf;
mod erf_poly;
mod erfc;
mod erff;
mod erffc;
mod exp;
mod exp10;
mod exp10f;
mod exp10m1;
mod exp10m1f;
mod exp2;
mod exp2f;
mod exp2m1;
mod exp2m1f;
mod expf;
mod expm1;
mod expm1f;
mod floor;
mod hypot;
mod hypotf;
mod j0;
mod j0_coeffs;
mod j0f;
mod j0f_coeffs;
mod j1;
mod j1_coeffs;
mod j1f;
mod j1f_coeffs;
mod log;
mod log10;
mod log10_dyadic;
mod log10f;
mod log10p1;
mod log10p1_tables;
mod log10p1f;
mod log1p;
mod log1p_dd;
mod log1p_dyadic;
mod log1p_dyadic_tables;
mod log1pf;
mod log2;
mod log2_dyadic;
mod log2f;
mod log2p1;
mod log2p1_dyadic_tables;
mod log2p1_tables;
mod log2p1f;
mod log_dyadic;
mod log_range_reduction;
mod logf;
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
mod sinh;
mod sinhf;
mod sqrtf;
mod tan;
mod tanf;
mod tanh;
mod tanhf;
mod tanpi;
mod tanpif;
mod triple_double;
mod trunc;
mod u256;
mod y0f;
mod y0f_coeffs;

pub use acos::f_acos;
pub use acosf::f_acosf;
pub use acosh::f_acosh;
pub use acoshf::f_acoshf;
pub use acospi::f_acospi;
pub use acospif::f_acospif;
pub use asin::f_asin;
pub use asinf::f_asinf;
pub use asinh::f_asinh;
pub use asinhf::f_asinhf;
pub use asinpi::f_asinpi;
pub use asinpif::f_asinpif;
pub use atan::f_atan;
pub use atan2::f_atan2;
pub use atan2f::f_atan2f;
pub use atan2pi::f_atan2pi;
pub use atan2pif::f_atan2pif;
pub use atanf::f_atanf;
pub use atanh::f_atanh;
pub use atanhf::f_atanhf;
pub use atanpi::f_atanpi;
pub use atanpif::f_atanpif;
pub use cbrt::f_cbrt;
pub use cbrtf::{cbrtf, f_cbrtf};
pub use ceil::{ceil, ceilf};
pub use common::{copysignfk, copysignk};
pub use compound::f_compound;
pub use compound_m1::f_compound_m1;
pub use compound_m1f::f_compound_m1f;
pub use compoundf::f_compoundf;
pub use cosf::f_cosf;
pub use cosh::f_cosh;
pub use coshf::f_coshf;
pub use cot::f_cot;
pub use cotf::f_cotf;
pub use csc::f_csc;
pub use cscf::f_cscf;
pub use erf::f_erf;
pub use erfc::f_erfc;
pub use erff::f_erff;
pub use erffc::f_erfcf;
pub use exp::{exp, f_exp};
pub use exp2::f_exp2;
pub use exp2f::f_exp2f;
pub use exp2m1::f_exp2m1;
pub use exp2m1f::f_exp2m1f;
pub use exp10::f_exp10;
pub use exp10f::f_exp10f;
pub use exp10m1::f_exp10m1;
pub use exp10m1f::f_exp10m1f;
pub use expf::{expf, f_expf};
pub use expm1::f_expm1;
pub use expm1f::f_expm1f;
pub use floor::{floor, floorf};
pub use hypot::f_hypot;
pub use hypotf::{f_hypot3f, f_hypotf};
pub use j0::f_j0;
pub use j0f::f_j0f;
pub use j1::f_j1;
pub use j1f::f_j1f;
pub use log::{f_log, log};
pub use log1p::f_log1p;
pub use log1pf::f_log1pf;
pub use log2::f_log2;
pub use log2f::f_log2f;
pub use log2p1::f_log2p1;
pub use log2p1f::f_log2p1f;
pub use log10::f_log10;
pub use log10f::f_log10f;
pub use log10p1::f_log10p1;
pub use log10p1f::f_log10p1f;
pub use logf::f_logf;
pub use logf::logf;
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
pub use sinh::f_sinh;
pub use sinhf::f_sinhf;
pub use sqrtf::sqrtf;
pub use tan::f_tan;
pub use tanf::f_tanf;
pub use tanh::f_tanh;
pub use tanhf::f_tanhf;
pub use tanpi::f_tanpi;
pub use tanpif::f_tanpif;
pub use trunc::{trunc, truncf};
pub use y0f::f_y0f;
