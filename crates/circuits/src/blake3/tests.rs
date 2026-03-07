// Copyright 2025 Irreducible Inc.
//! Comprehensive BLAKE3 tests: completeness, soundness, and isolation.

use binius_core::{verify::verify_constraints, word::Word};
use binius_frontend::CircuitBuilder;
use rand::{Rng, RngCore, SeedableRng, rngs::StdRng};
use rstest::rstest;

use super::{Blake3, constants::*, reference};

// ============================================================================
// Layer 1: Completeness tests (correct digest → constraints satisfied)
// ============================================================================

/// Build a circuit, populate with the correct digest, verify constraints pass.
fn test_circuit_with_message(message: &[u8]) -> [u8; 32] {
	let mut builder = CircuitBuilder::new();
	let blake3 = Blake3::new_witness(&mut builder, message.len());
	let circuit = builder.build();

	let expected_digest = reference::blake3_hash(message);

	let mut witness = circuit.new_witness_filler();
	blake3.populate_message(&mut witness, message);
	blake3.populate_digest(&mut witness, &expected_digest);

	circuit
		.populate_wire_witness(&mut witness)
		.expect("Circuit should accept valid witness");

	// Read back the digest from witness wires
	let circuit_digest: [u8; 32] = core::array::from_fn(|i| {
		let word_idx = i / 4;
		let byte_idx = i % 4;
		let word_val = witness[blake3.digest[word_idx]].0;
		((word_val >> (byte_idx * 8)) & 0xFF) as u8
	});

	let cs = circuit.constraint_system();
	verify_constraints(cs, &witness.into_value_vec()).expect("All constraints should be satisfied");

	circuit_digest
}

/// Validate circuit output matches the blake3 crate oracle.
fn validate_against_oracle(message: &[u8], test_name: &str) {
	let oracle_digest: [u8; 32] = blake3::hash(message).into();
	let reference_digest = reference::blake3_hash(message);
	assert_eq!(
		reference_digest, oracle_digest,
		"Reference mismatch for {test_name}: ref={reference_digest:02x?}, oracle={oracle_digest:02x?}"
	);

	let circuit_digest = test_circuit_with_message(message);
	assert_eq!(
		circuit_digest, oracle_digest,
		"Circuit mismatch for {test_name}: circuit={circuit_digest:02x?}, oracle={oracle_digest:02x?}"
	);
}

#[rstest]
#[case(0, "empty")]
#[case(1, "single_byte")]
#[case(63, "under_one_block")]
#[case(64, "exact_one_block")]
#[case(65, "over_one_block")]
#[case(127, "under_two_blocks")]
#[case(128, "exact_two_blocks")]
#[case(129, "over_two_blocks")]
#[case(1023, "under_one_chunk")]
#[case(1024, "exact_one_chunk")]
#[case(1025, "over_one_chunk")]
#[case(2048, "exact_two_chunks")]
#[case(2049, "over_two_chunks")]
#[case(3072, "exact_three_chunks")]
#[case(4096, "exact_four_chunks")]
#[case(4097, "over_four_chunks")]
#[case(8192, "exact_eight_chunks")]
fn test_blake3_boundary_lengths(#[case] len: usize, #[case] name: &str) {
	let mut rng = StdRng::seed_from_u64(len as u64);
	let mut message = vec![0u8; len];
	rng.fill_bytes(&mut message);
	validate_against_oracle(&message, name);
}

#[test]
fn test_blake3_known_vectors() {
	validate_against_oracle(b"", "empty_literal");
	validate_against_oracle(b"abc", "abc");
	validate_against_oracle(b"a", "single_a");
	validate_against_oracle(
		b"The quick brown fox jumps over the lazy dog",
		"pangram",
	);
}

#[test]
fn test_blake3_block_boundary_patterns() {
	validate_against_oracle(&[0x42u8; BLOCK_LEN - 1], "just_under_block");
	validate_against_oracle(&[0x42u8; BLOCK_LEN], "exact_block");
	validate_against_oracle(&[0x42u8; BLOCK_LEN + 1], "just_over_block");
	validate_against_oracle(&[0xAAu8; 2 * BLOCK_LEN], "two_blocks");
	validate_against_oracle(&vec![0xFFu8; CHUNK_LEN], "full_chunk");
}

#[test]
fn test_blake3_random_messages() {
	let mut rng = StdRng::seed_from_u64(42);

	for i in 0..100 {
		let size = rng.random_range(0..=1024);
		let mut message = vec![0u8; size];
		rng.fill(&mut message[..]);
		validate_against_oracle(&message, &format!("random_{i}_size_{size}"));
	}
}

#[test]
fn test_blake3_all_zeros() {
	for &len in &[0, 1, 64, 128, 1024, 2048, 4096] {
		validate_against_oracle(&vec![0u8; len], &format!("zeros_{len}"));
	}
}

#[test]
fn test_blake3_all_ones() {
	for &len in &[1, 64, 128, 1024, 2048, 4096] {
		validate_against_oracle(&vec![0xFFu8; len], &format!("ones_{len}"));
	}
}

#[test]
fn test_blake3_multi_chunk_non_power_of_two() {
	let mut rng = StdRng::seed_from_u64(77);
	for &len in &[3000, 5000, 6144, 7000] {
		let mut message = vec![0u8; len];
		rng.fill_bytes(&mut message);
		validate_against_oracle(&message, &format!("multi_chunk_{len}"));
	}
}

#[test]
fn test_blake3_random_multi_chunk_messages() {
	let mut rng = StdRng::seed_from_u64(99);

	for i in 0..20 {
		let size = rng.random_range(1025..=8192);
		let mut message = vec![0u8; size];
		rng.fill(&mut message[..]);
		validate_against_oracle(&message, &format!("random_multi_{i}_size_{size}"));
	}
}

// ============================================================================
// Layer 2: Soundness tests (wrong digest → constraints fail)
// ============================================================================

/// Helper: build circuit, populate with given (possibly wrong) digest, return result.
fn try_circuit_with_digest(message: &[u8], digest: &[u8; 32]) -> Result<(), String> {
	let mut builder = CircuitBuilder::new();
	let blake3 = Blake3::new_witness(&mut builder, message.len());
	let circuit = builder.build();

	let mut witness = circuit.new_witness_filler();
	blake3.populate_message(&mut witness, message);
	blake3.populate_digest(&mut witness, digest);

	circuit
		.populate_wire_witness(&mut witness)
		.map_err(|e| e.to_string())?;

	let cs = circuit.constraint_system();
	verify_constraints(cs, &witness.into_value_vec()).map_err(|e| e.to_string())
}

#[test]
fn test_wrong_digest_rejected() {
	let message = b"hello blake3";
	let wrong_digest = [0u8; 32];
	assert!(
		try_circuit_with_digest(message, &wrong_digest).is_err(),
		"All-zero digest should be rejected"
	);
}

#[test]
fn test_single_bit_flip_rejected() {
	let message = b"test bit flip";
	let mut digest: [u8; 32] = blake3::hash(message).into();
	digest[0] ^= 1; // Flip lowest bit of first byte
	assert!(
		try_circuit_with_digest(message, &digest).is_err(),
		"Single bit flip in digest should be rejected"
	);
}

#[test]
fn test_swapped_digest_words_rejected() {
	let message = b"test swapped words";
	let mut digest: [u8; 32] = blake3::hash(message).into();
	// Swap first and second 4-byte words
	let (a, b) = digest.split_at_mut(4);
	let mut tmp = [0u8; 4];
	tmp.copy_from_slice(&a[..4]);
	a[..4].copy_from_slice(&b[..4]);
	b[..4].copy_from_slice(&tmp);
	assert!(
		try_circuit_with_digest(message, &digest).is_err(),
		"Swapped digest words should be rejected"
	);
}

#[test]
fn test_wrong_message_content_rejected() {
	let message_a = b"message A";
	let message_b = b"message B";
	// Use digest of message_b but provide message_a
	let digest_b: [u8; 32] = blake3::hash(message_b.as_slice()).into();

	let mut builder = CircuitBuilder::new();
	let blake3_circuit = Blake3::new_witness(&mut builder, message_a.len());
	let circuit = builder.build();

	let mut witness = circuit.new_witness_filler();
	blake3_circuit.populate_message(&mut witness, message_a);
	blake3_circuit.populate_digest(&mut witness, &digest_b);

	let result = circuit.populate_wire_witness(&mut witness);
	if result.is_ok() {
		let cs = circuit.constraint_system();
		assert!(
			verify_constraints(cs, &witness.into_value_vec()).is_err(),
			"Wrong message content should be rejected"
		);
	}
}

// ============================================================================
// Layer 3: Isolation tests (compression function + constraint counts)
// ============================================================================

#[test]
fn test_wrong_digest_rejected_multi_chunk() {
	let message = vec![0xABu8; 2048];
	let wrong_digest = [0u8; 32];
	assert!(
		try_circuit_with_digest(&message, &wrong_digest).is_err(),
		"All-zero digest should be rejected for multi-chunk"
	);
}

#[test]
fn test_single_bit_flip_rejected_multi_chunk() {
	let message = vec![0xCDu8; 2048];
	let mut digest: [u8; 32] = blake3::hash(&message).into();
	digest[0] ^= 1;
	assert!(
		try_circuit_with_digest(&message, &digest).is_err(),
		"Single bit flip should be rejected for multi-chunk"
	);
}

// ============================================================================
// Layer 3: Isolation tests (compression function + constraint counts)
// ============================================================================

#[test]
fn test_reference_matches_oracle() {
	let mut rng = StdRng::seed_from_u64(123);
	for _ in 0..50 {
		let len = rng.random_range(0..=1024);
		let mut message = vec![0u8; len];
		rng.fill(&mut message[..]);

		let reference_digest = reference::blake3_hash(&message);
		let oracle_digest: [u8; 32] = blake3::hash(&message).into();
		assert_eq!(reference_digest, oracle_digest, "Reference/oracle mismatch at len={len}");
	}
}

#[test]
fn test_reference_matches_oracle_multi_chunk() {
	let mut rng = StdRng::seed_from_u64(456);
	for _ in 0..20 {
		let len = rng.random_range(1025..=8192);
		let mut message = vec![0u8; len];
		rng.fill(&mut message[..]);

		let reference_digest = reference::blake3_hash(&message);
		let oracle_digest: [u8; 32] = blake3::hash(&message).into();
		assert_eq!(
			reference_digest, oracle_digest,
			"Reference/oracle mismatch at len={len}"
		);
	}
}

#[test]
fn test_compression_traced() {
	let cv = IV;
	let block_words = [0u32; 16];
	let trace = reference::compress_traced(&cv, &block_words, 0, 0, CHUNK_START | CHUNK_END | ROOT);

	// Verify trace output matches non-traced version
	let output = reference::compress(&cv, &block_words, 0, 0, CHUNK_START | CHUNK_END | ROOT);
	assert_eq!(trace.output, output);

	// Verify 7 round states were captured
	assert_eq!(trace.round_states.len(), ROUNDS);
}

#[test]
fn test_constraint_count_single_compression() {
	let mut builder = CircuitBuilder::new();
	let zero = builder.add_constant(Word(0));
	let cv: [_; 8] = std::array::from_fn(|_| builder.add_witness());
	let m: [_; 16] = std::array::from_fn(|_| builder.add_witness());
	let block_len = builder.add_witness();
	let flags = builder.add_witness();

	super::blake3_compress(&mut builder, &cv, &m, zero, zero, block_len, flags);

	let circuit = builder.build();
	let cs = circuit.constraint_system();
	let n_and = cs.n_and_constraints();

	// 7 rounds x 8 G calls = 56 G calls total.
	// With optimized 32-bit shift gates: ~13.4 AND per G call.
	// Regression check: pin the exact count to detect accidental changes.
	assert_eq!(
		n_and, 752,
		"Single compression should have exactly 752 AND constraints, got {n_and}"
	);
}

#[test]
fn test_empty_message_flags() {
	// Empty message: single block with block_len=0, flags=CHUNK_START|CHUNK_END|ROOT
	let oracle_digest: [u8; 32] = blake3::hash(b"").into();
	let circuit_digest = test_circuit_with_message(b"");
	assert_eq!(circuit_digest, oracle_digest, "Empty message digest mismatch");
}

#[test]
fn test_constraint_count_parent_compression() {
	let mut builder = CircuitBuilder::new();
	let left: [_; 8] = std::array::from_fn(|_| builder.add_witness());
	let right: [_; 8] = std::array::from_fn(|_| builder.add_witness());

	super::blake3_parent(&mut builder, &left, &right, false);

	let circuit = builder.build();
	let cs = circuit.constraint_system();
	let n_and = cs.n_and_constraints();

	assert_eq!(
		n_and, 752,
		"Parent compression should have exactly 752 AND constraints, got {n_and}"
	);
}

#[test]
fn test_wire_readback_matches_reference() {
	let message = b"wire readback test";
	let mut builder = CircuitBuilder::new();
	let blake3 = Blake3::new_witness(&mut builder, message.len());
	let circuit = builder.build();

	let expected_digest = reference::blake3_hash(message);

	let mut witness = circuit.new_witness_filler();
	blake3.populate_message(&mut witness, message);
	blake3.populate_digest(&mut witness, &expected_digest);

	circuit
		.populate_wire_witness(&mut witness)
		.expect("Should accept valid witness");

	for i in 0..8 {
		let wire_val = witness[blake3.digest[i]].0 as u32;
		let expected_word =
			u32::from_le_bytes(expected_digest[i * 4..(i + 1) * 4].try_into().unwrap());
		assert_eq!(
			wire_val, expected_word,
			"Digest word {i} mismatch: wire={wire_val:#010x}, expected={expected_word:#010x}"
		);
	}
}

#[test]
fn test_wire_readback_multi_chunk() {
	let mut rng = StdRng::seed_from_u64(888);
	let mut message = vec![0u8; 2048];
	rng.fill_bytes(&mut message);

	let mut builder = CircuitBuilder::new();
	let blake3 = Blake3::new_witness(&mut builder, message.len());
	let circuit = builder.build();

	let expected_digest = reference::blake3_hash(&message);

	let mut witness = circuit.new_witness_filler();
	blake3.populate_message(&mut witness, &message);
	blake3.populate_digest(&mut witness, &expected_digest);

	circuit
		.populate_wire_witness(&mut witness)
		.expect("Should accept valid witness");

	for i in 0..8 {
		let wire_val = witness[blake3.digest[i]].0 as u32;
		let expected_word =
			u32::from_le_bytes(expected_digest[i * 4..(i + 1) * 4].try_into().unwrap());
		assert_eq!(
			wire_val, expected_word,
			"Digest word {i} mismatch: wire={wire_val:#010x}, expected={expected_word:#010x}"
		);
	}
}
