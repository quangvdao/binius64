// Copyright 2026 The Binius Developers
//! Constraint Equivalence Checker (ceck)
//!
//! Proves or disproves equivalence between constraint systems using
//! random testing (randblast) and optionally Z3 SMT solving.

pub mod ast;
pub mod export;
pub mod parser;
pub mod randblast;
#[cfg(feature = "z3")]
pub mod smt_check;
pub mod translate;
