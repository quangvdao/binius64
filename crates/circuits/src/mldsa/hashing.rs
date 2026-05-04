// Copyright 2025 Irreducible Inc.

use binius_frontend::{CircuitBuilder, Wire};

use crate::keccak::fixed_length;

use super::MldsaParams;

/// Computes the fixed-shape ML-DSA final challenge hash:
///
/// ```text
/// c_tilde_prime = SHAKE256(mu || w1Encode(w1_prime), 32)
/// ```
///
/// Inputs are packed as little-endian 64-bit wires: public `mu` followed by private `w1` bytes.
pub fn final_challenge_hash_for<P: MldsaParams>(
	builder: &CircuitBuilder,
	mu_and_w1_bytes: &[Wire],
) -> Vec<Wire> {
	assert_eq!(
		mu_and_w1_bytes.len(),
		P::FINAL_CHALLENGE_INPUT_BYTES.div_ceil(8),
		"{} final challenge input has wrong packed length",
		P::label(),
	);

	fixed_length::shake256(
		builder,
		mu_and_w1_bytes,
		P::FINAL_CHALLENGE_INPUT_BYTES,
		P::C_TILDE_BYTES,
	)
}
