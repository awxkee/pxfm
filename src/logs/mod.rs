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
#![deny(unreachable_pub)]
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

pub(crate) use log::log_dyadic;
pub use log::{f_log, log};
pub use log1p::f_log1p;
pub(crate) use log1p::log1p_f64_dyadic;
pub(crate) use log1p_dd::log1p_f64_dd;
pub use log1pf::f_log1pf;
#[cfg(not(any(
    all(
        any(target_arch = "x86", target_arch = "x86_64"),
        target_feature = "fma"
    ),
    all(target_arch = "aarch64", target_feature = "neon")
)))]
pub(crate) use log2::LOG_CD;
pub use log2::f_log2;
pub(crate) use log2::{LOG_COEFFS, LOG_RANGE_REDUCTION};
pub use log2f::f_log2f;
pub(crate) use log2f::{LOG2_R, dirty_log2f};
pub use log2p1::f_log2p1;
pub use log2p1f::f_log2p1f;
pub(crate) use log10::LOG_R_DD;
pub use log10::f_log10;
pub use log10f::f_log10f;
pub use log10p1::f_log10p1;
pub use log10p1f::f_log10p1f;
#[allow(unused)]
pub(crate) use logf::LOG_REDUCTION_F32;
pub use logf::{f_logf, logf};
