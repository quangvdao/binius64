// Copyright 2026 The Binius Developers

//! Binius-native Keccak-f[1600] proving experiments.
//!
//! This crate is the implementation home for the Keccak-specific prover sketched
//! in `docs/keccak-prove`. The first implementation slices expose native trace
//! construction and the BitAnd-style chi residual columns that the proving path
//! will feed into an NTT-backed outer reduction.

#![warn(rustdoc::missing_crate_level_docs)]

pub mod bit_ntt;
pub mod constants;
pub mod layout;
pub mod operands;
pub mod round_message;
pub mod shift_claims;
pub mod shift_operands;
pub mod trace;
pub mod unrolled;
pub mod witness;
