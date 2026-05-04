// Copyright 2025 Irreducible Inc.

use binius_frontend::{CircuitBuilder, Wire};

use crate::keccak::fixed_length;

use super::mldsa44;

/// Computes the fixed-shape ML-DSA-44 final challenge hash:
///
/// ```text
/// c_tilde_prime = SHAKE256(mu || w1Encode(w1_prime), 32)
/// ```
///
/// Inputs are packed as little-endian 64-bit wires and must cover exactly 832 bytes:
/// 64 bytes of public `mu` followed by 768 private `w1` bytes.
pub fn mldsa44_final_challenge_hash(
	builder: &CircuitBuilder,
	mu_and_w1_bytes: &[Wire],
) -> [Wire; mldsa44::C_TILDE_BYTES / 8] {
	assert_eq!(
		mu_and_w1_bytes.len(),
		mldsa44::FINAL_CHALLENGE_INPUT_BYTES.div_ceil(8),
		"ML-DSA-44 final challenge expects 832 bytes packed into 104 words",
	);

	let output = fixed_length::shake256(
		builder,
		mu_and_w1_bytes,
		mldsa44::FINAL_CHALLENGE_INPUT_BYTES,
		mldsa44::C_TILDE_BYTES,
	);
	output.try_into().unwrap()
}
