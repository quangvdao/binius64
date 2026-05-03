// Copyright 2026 The Binius Developers

//! Binius-native Keccak-f[1600] proving experiments.
//!
//! This crate is the implementation home for the Keccak-specific prover sketched
//! in `docs/keccak-prove`. The first implementation slices expose native trace
//! construction and the BitAnd-style chi residual columns that the proving path
//! will feed into an NTT-backed outer reduction.

#![warn(rustdoc::missing_crate_level_docs)]

pub mod constants;
pub mod operands;
pub mod trace;
