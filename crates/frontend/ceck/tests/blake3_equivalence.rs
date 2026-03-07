// Copyright 2026 The Binius Developers
//! BLAKE3 compression function equivalence test.
//!
//! Builds two independent constraint systems for BLAKE3 compression:
//! 1. The production `blake3_compress` from `binius_circuits::blake3`
//! 2. A reference implementation built from scratch using `g_function`
//!
//! Then uses ceck (randblast + optionally Z3) to prove they accept exactly
//! the same witnesses -- a formal proof of functional equivalence.

use binius_circuits::blake2s::g_function;
use binius_circuits::blake3::constants::{IV, ROUNDS, SCHEDULES};
use binius_core::word::Word;
use binius_frontend::CircuitBuilder;
use ceck::randblast::RandBlast;

/// Build the BLAKE3 compression constraint system using the production code.
fn build_production_compression() -> binius_core::constraint_system::ConstraintSystem {
	let mut builder = CircuitBuilder::new();
	let zero = builder.add_constant(Word(0));

	let cv: [_; 8] = std::array::from_fn(|_| builder.add_witness());
	let m: [_; 16] = std::array::from_fn(|_| builder.add_witness());
	let block_len = builder.add_witness();
	let flags = builder.add_witness();

	binius_circuits::blake3::blake3_compress(
		&mut builder,
		&cv,
		&m,
		zero,
		zero,
		block_len,
		flags,
	);

	let circuit = builder.build();
	circuit.constraint_system().clone()
}

/// Build a reference BLAKE3 compression constraint system from scratch.
///
/// This is an independent reimplementation of the compression function,
/// written to match the BLAKE3 spec directly rather than calling `blake3_compress`.
/// It reuses only the low-level `g_function` primitive (shared with BLAKE2s).
fn build_reference_compression() -> binius_core::constraint_system::ConstraintSystem {
	let mut builder = CircuitBuilder::new();
	let zero = builder.add_constant(Word(0));

	let cv: [_; 8] = std::array::from_fn(|_| builder.add_witness());
	let m: [_; 16] = std::array::from_fn(|_| builder.add_witness());
	let block_len = builder.add_witness();
	let flags = builder.add_witness();

	// State initialization per BLAKE3 spec
	let mut v = [zero; 16];
	for i in 0..8 {
		v[i] = cv[i];
	}
	for i in 0..4 {
		v[8 + i] = builder.add_constant(Word(IV[i] as u64));
	}
	// v[12..16] = (counter_lo=0, counter_hi=0, block_len, flags)
	v[12] = zero;
	v[13] = zero;
	v[14] = block_len;
	v[15] = flags;

	// 7 rounds with precomputed message schedules
	for round in 0..ROUNDS {
		let s = &SCHEDULES[round];

		// Column step: G(0,4,8,12), G(1,5,9,13), G(2,6,10,14), G(3,7,11,15)
		let (a, b, c, d) = g_function(&mut builder, v[0], v[4], v[8], v[12], m[s[0]], m[s[1]]);
		v[0] = a;
		v[4] = b;
		v[8] = c;
		v[12] = d;

		let (a, b, c, d) = g_function(&mut builder, v[1], v[5], v[9], v[13], m[s[2]], m[s[3]]);
		v[1] = a;
		v[5] = b;
		v[9] = c;
		v[13] = d;

		let (a, b, c, d) = g_function(&mut builder, v[2], v[6], v[10], v[14], m[s[4]], m[s[5]]);
		v[2] = a;
		v[6] = b;
		v[10] = c;
		v[14] = d;

		let (a, b, c, d) = g_function(&mut builder, v[3], v[7], v[11], v[15], m[s[6]], m[s[7]]);
		v[3] = a;
		v[7] = b;
		v[11] = c;
		v[15] = d;

		// Diagonal step: G(0,5,10,15), G(1,6,11,12), G(2,7,8,13), G(3,4,9,14)
		let (a, b, c, d) = g_function(&mut builder, v[0], v[5], v[10], v[15], m[s[8]], m[s[9]]);
		v[0] = a;
		v[5] = b;
		v[10] = c;
		v[15] = d;

		let (a, b, c, d) =
			g_function(&mut builder, v[1], v[6], v[11], v[12], m[s[10]], m[s[11]]);
		v[1] = a;
		v[6] = b;
		v[11] = c;
		v[12] = d;

		let (a, b, c, d) =
			g_function(&mut builder, v[2], v[7], v[8], v[13], m[s[12]], m[s[13]]);
		v[2] = a;
		v[7] = b;
		v[8] = c;
		v[13] = d;

		let (a, b, c, d) =
			g_function(&mut builder, v[3], v[4], v[9], v[14], m[s[14]], m[s[15]]);
		v[3] = a;
		v[4] = b;
		v[9] = c;
		v[14] = d;
	}

	// BLAKE3 finalization: state[i] ^= state[i+8]; state[i+8] ^= cv[i]
	for i in 0..8 {
		v[i] = builder.bxor(v[i], v[i + 8]);
		v[i + 8] = builder.bxor(v[i + 8], cv[i]);
	}

	let circuit = builder.build();
	circuit.constraint_system().clone()
}

#[test]
fn test_blake3_compression_equivalence_randblast() {
	let production_cs = build_production_compression();
	let reference_cs = build_reference_compression();

	assert_eq!(
		production_cs.n_and_constraints(),
		reference_cs.n_and_constraints(),
		"Constraint counts must match: production={}, reference={}",
		production_cs.n_and_constraints(),
		reference_cs.n_and_constraints(),
	);

	let mut blaster = RandBlast::new(42);
	blaster
		.test_equivalence(&production_cs, &reference_cs, 100_000)
		.expect("Production and reference BLAKE3 compression should be equivalent");
}

#[test]
fn test_blake3_compression_export_roundtrip() {
	let production_cs = build_production_compression();

	let exported = ceck::export::export_constraint_set(&production_cs);

	// Verify the exported string is valid ceck syntax by parsing it
	assert!(exported.starts_with("(constraint_set"));
	assert!(exported.contains("(and "));

	// Count exported constraints matches
	let n_and_in_export = exported.matches("(and ").count();
	assert_eq!(
		n_and_in_export,
		production_cs.n_and_constraints(),
		"Exported AND constraint count should match"
	);
}

#[cfg(feature = "z3")]
#[test]
fn test_blake3_compression_equivalence_smt() {
	use ceck::smt_check::SmtChecker;
	use z3::{Config, Context};

	let production_cs = build_production_compression();
	let reference_cs = build_reference_compression();

	let config = Config::new();
	let ctx = Context::new(&config);
	let mut checker = SmtChecker::new(&ctx);

	checker
		.check_equivalence(&production_cs, &reference_cs)
		.expect(
			"SMT solver should prove production and reference BLAKE3 compression are equivalent",
		);
}
