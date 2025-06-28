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
#![forbid(unsafe_code)]
#![deny(unreachable_pub)]
#![allow(
    clippy::excessive_precision,
    clippy::approx_constant,
    clippy::manual_range_contains
)]
mod acos;
mod acosf;
mod acospif;
mod asin;
mod asin_eval_dyadic;
mod asinf;
mod asinpif;
mod atan;
mod atan2;
mod atan2f;
mod atanf;
mod bits;
mod cbrt;
mod cbrtf;
mod ceil;
mod common;
mod cosf;
mod coshf;
mod dekker;
mod dyadic_float;
mod exp;
mod exp10;
mod exp10f;
mod exp10m1f;
mod exp2;
mod exp2f;
mod exp2m1f;
mod expf;
mod expm1f;
mod floor;
mod hypot;
mod j1;
mod log;
mod log10;
mod log10_dyadic;
mod log10f;
mod log10p1f;
mod log1pf;
mod log2;
mod log2_dyadic;
mod log2f;
mod log2p1f;
mod log_dyadic;
mod log_range_reduction;
mod logf;
mod pow;
mod powf;
mod round;
mod round_ties_even;
mod sin;
mod sincos;
mod sincos_dyadic;
mod sincosf;
mod sincospi;
mod sincospi_tables;
mod sincospif;
mod sinf;
mod sinhf;
mod sqrtf;
mod tan;
mod tanf;
mod tanhf;
mod tanpi;
mod tanpif;
mod trunc;

pub use acos::f_acos;
pub use acosf::f_acosf;
pub use acospif::f_acospif;
pub use asin::f_asin;
pub use asinf::f_asinf;
pub use asinpif::f_asinpif;
pub use atan::f_atan;
pub use atan2::f_atan2;
pub use atan2f::f_atan2f;
pub use atanf::f_atanf;
pub use cbrt::f_cbrt;
pub use cbrtf::{cbrtf, f_cbrtf};
pub use ceil::{ceil, ceilf};
pub use common::{copysignfk, copysignk};
pub use cosf::f_cosf;
pub use coshf::f_coshf;
pub use exp::{exp, f_exp};
pub use exp2::f_exp2;
pub use exp2f::f_exp2f;
pub use exp2m1f::f_exp2m1f;
pub use exp10::f_exp10;
pub use exp10f::f_exp10f;
pub use exp10m1f::f_exp10m1f;
pub use expf::{expf, f_expf};
pub use expm1f::f_expm1f;
pub use floor::{floor, floorf};
pub use hypot::{f_hypot3f, f_hypotf};
pub use j1::f_j1;
pub use log::{f_log, log};
pub use log1pf::f_log1pf;
pub use log2::f_log2;
pub use log2f::f_log2f;
pub use log2p1f::f_log2p1f;
pub use log10::f_log10;
pub use log10f::f_log10f;
pub use log10p1f::f_log10p1f;
pub use logf::f_logf;
pub use logf::logf;
pub use pow::{f_pow, pow};
pub use powf::{dirty_powf, f_powf, powf};
pub use round::{round, roundf};
pub use round_ties_even::{round_ties_even, roundf_ties_even};
pub use sin::{f_cos, f_sin};
pub use sincos::f_sincos;
pub use sincosf::f_sincosf;
pub use sincospi::{f_cospi, f_sinpi};
pub use sincospif::{f_cospif, f_sinpif};
pub use sinf::f_sinf;
pub use sinhf::f_sinhf;
pub use sqrtf::sqrtf;
pub use tan::f_tan;
pub use tanf::f_tanf;
pub use tanhf::f_tanhf;
pub use tanpi::f_tanpi;
pub use tanpif::f_tanpif;
pub use trunc::{trunc, truncf};
