// Copyright 2025 Irreducible Inc.

//! Prototype ML-DSA bit-heavy circuit helpers.
//!
//! This module starts with the ML-DSA-44 hidden-signature target shape from the
//! top-level lattice-sig-aggregation executable spec. Public-only hashing and key
//! preprocessing are intentionally hoisted out of this circuit layer.

use binius_frontend::{CircuitBuilder, Wire};

use crate::keccak::fixed_length;

/// ML-DSA-44 constants for the first prototype target.
pub mod mldsa44 {
	/// Number of rows in the public matrix.
	pub const K: usize = 4;
	/// Number of columns in the public matrix.
	pub const L: usize = 4;
	pub const TAU: usize = 39;
	pub const BETA: u64 = 78;
	pub const GAMMA1: u64 = 1 << 17;
	pub const C_TILDE_BYTES: usize = 32;
	pub const MU_BYTES: usize = 64;
	pub const POLY_W1_PACKED_BYTES: usize = 192;
	pub const W1_ENCODE_BYTES: usize = K * POLY_W1_PACKED_BYTES;
	pub const FINAL_CHALLENGE_INPUT_BYTES: usize = MU_BYTES + W1_ENCODE_BYTES;
	pub const FINAL_CHALLENGE_KECCAK_F_CALLS: usize = 7;

	/// Packed `y = gamma1 - z` range equivalent to `|z| < gamma1 - beta`.
	pub const Z_NORM_PACKED_Y_MIN: u64 = BETA + 1;
	pub const Z_NORM_PACKED_Y_MAX: u64 = 2 * GAMMA1 - BETA - 1;
	pub const Z_COEFFICIENTS: usize = L * 256;
	pub const Z_BITS_PER_COEFF: usize = 18;
	pub const POLY_Z_PACKED_BYTES: usize = 576;
	pub const Z_PACKED_BYTES: usize = L * POLY_Z_PACKED_BYTES;
	pub const Z_PACKED_WORDS: usize = Z_PACKED_BYTES / 8;

	pub const SAMPLE_IN_BALL_SIGN_BYTES: usize = 8;
	pub const SAMPLE_IN_BALL_ONE_BLOCK_DRAW_CAP_BYTES: usize = 128;
	pub const SAMPLE_IN_BALL_ONE_BLOCK_STREAM_BYTES: usize =
		SAMPLE_IN_BALL_SIGN_BYTES + SAMPLE_IN_BALL_ONE_BLOCK_DRAW_CAP_BYTES;
	pub const SAMPLE_IN_BALL_ONE_BLOCK_KECCAK_F_CALLS: usize = 1;
}

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

/// Computes the fixed one-block ML-DSA-44 `SampleInBall` SHAKE stream.
///
/// This is the first fixed-cap sampler shape from the top-level plan:
///
/// ```text
/// stream = SHAKE256(c_tilde, 136)
/// sign bytes = stream[0..8]
/// position draws = stream[8..136]
/// ```
///
/// The Fisher-Yates rejection/update relation is layered on top of this stream. The cap is part of
/// the circuit identity, so repeated batches must not mix this helper with two-block sampler shapes.
pub fn mldsa44_sample_in_ball_one_block_stream(
	builder: &CircuitBuilder,
	c_tilde: &[Wire],
) -> [Wire; mldsa44::SAMPLE_IN_BALL_ONE_BLOCK_STREAM_BYTES / 8] {
	assert_eq!(
		c_tilde.len(),
		mldsa44::C_TILDE_BYTES / 8,
		"ML-DSA-44 c_tilde expects 32 bytes packed into 4 words",
	);

	let stream = fixed_length::shake256(
		builder,
		c_tilde,
		mldsa44::C_TILDE_BYTES,
		mldsa44::SAMPLE_IN_BALL_ONE_BLOCK_STREAM_BYTES,
	);
	stream.try_into().unwrap()
}

/// Enforces the ML-DSA-44 `z` norm check on decoded BitPack coefficients.
///
/// FIPS stores each `z` coefficient through `BitPack(z, -gamma1 + 1, gamma1)`.
/// Equivalently, the packed integer is:
///
/// ```text
/// y = gamma1 - z
/// ```
///
/// For ML-DSA-44, `||z||_infty < gamma1 - beta` is therefore:
///
/// ```text
/// 79 <= y <= 262065
/// ```
///
/// The inputs here are one 64-bit wire per decoded packed `y` coefficient.
pub fn assert_mldsa44_z_norm_from_packed_y(builder: &CircuitBuilder, packed_y_coeffs: &[Wire]) {
	assert_eq!(
		packed_y_coeffs.len(),
		mldsa44::Z_COEFFICIENTS,
		"ML-DSA-44 z has L * 256 = 1024 coefficients",
	);

	let min = builder.add_constant_64(mldsa44::Z_NORM_PACKED_Y_MIN);
	let max = builder.add_constant_64(mldsa44::Z_NORM_PACKED_Y_MAX);

	for (i, &y) in packed_y_coeffs.iter().enumerate() {
		let gte_min = builder.icmp_ule(min, y);
		let lte_max = builder.icmp_ule(y, max);
		builder.assert_true(format!("mldsa44_z_norm_gte_min[{i}]"), gte_min);
		builder.assert_true(format!("mldsa44_z_norm_lte_max[{i}]"), lte_max);
	}
}

/// Decodes ML-DSA-44 `z` from its fixed BitPack representation.
///
/// FIPS packs `z` coefficients as contiguous 18-bit little-endian integers:
///
/// ```text
/// y_i = gamma1 - z_i
/// ```
///
/// This helper keeps the signature hidden by taking the packed stream as witness wires. The output
/// is one 64-bit wire per decoded packed `y_i` coefficient.
pub fn mldsa44_decode_z_packed_y(builder: &CircuitBuilder, z_words: &[Wire]) -> Vec<Wire> {
	assert_eq!(
		z_words.len(),
		mldsa44::Z_PACKED_WORDS,
		"ML-DSA-44 z expects 2304 bytes packed into 288 words",
	);

	let mask = builder.add_constant_64((1 << mldsa44::Z_BITS_PER_COEFF) - 1);
	let mut coeffs = Vec::with_capacity(mldsa44::Z_COEFFICIENTS);

	for coeff_idx in 0..mldsa44::Z_COEFFICIENTS {
		let bit_offset = coeff_idx * mldsa44::Z_BITS_PER_COEFF;
		let word_idx = bit_offset / 64;
		let shift = bit_offset % 64;

		let coeff = if shift <= 64 - mldsa44::Z_BITS_PER_COEFF {
			builder.band(builder.shr(z_words[word_idx], shift as u32), mask)
		} else {
			let low_bits = 64 - shift;
			let high_bits = mldsa44::Z_BITS_PER_COEFF - low_bits;
			let high_mask = builder.add_constant_64((1 << high_bits) - 1);
			let low = builder.shr(z_words[word_idx], shift as u32);
			let high = builder.band(z_words[word_idx + 1], high_mask);
			let high_shifted = builder.shl(high, low_bits as u32);
			builder.band(builder.bor(low, high_shifted), mask)
		};

		coeffs.push(coeff);
	}

	coeffs
}

/// Decodes hidden ML-DSA-44 packed `z` bytes and enforces `||z||_infty < gamma1 - beta`.
pub fn assert_mldsa44_z_packed_bytes_norm(builder: &CircuitBuilder, z_words: &[Wire]) -> Vec<Wire> {
	let packed_y_coeffs = mldsa44_decode_z_packed_y(builder, z_words);
	assert_mldsa44_z_norm_from_packed_y(builder, &packed_y_coeffs);
	packed_y_coeffs
}

#[cfg(test)]
mod tests {
	use binius_core::{verify::verify_constraints, word::Word};
	use binius_frontend::CircuitBuilder;
	use rand::{RngCore, SeedableRng, rngs::StdRng};
	use sha3::{
		Shake256,
		digest::{ExtendableOutput, Update, XofReader},
	};

	use super::*;

	fn pack_mldsa44_z_y_coeffs(packed_y_coeff_values: &[u64]) -> Vec<u64> {
		assert_eq!(packed_y_coeff_values.len(), mldsa44::Z_COEFFICIENTS);

		let mut words = vec![0u64; mldsa44::Z_PACKED_WORDS];
		for (coeff_idx, &coeff) in packed_y_coeff_values.iter().enumerate() {
			assert!(coeff < (1 << mldsa44::Z_BITS_PER_COEFF));
			for bit in 0..mldsa44::Z_BITS_PER_COEFF {
				if (coeff >> bit) & 1 == 1 {
					let bit_idx = coeff_idx * mldsa44::Z_BITS_PER_COEFF + bit;
					words[bit_idx / 64] |= 1 << (bit_idx % 64);
				}
			}
		}

		words
	}

	fn verify_mldsa44_z_packed_bytes_norm_witness(packed_y_coeff_values: &[u64]) -> bool {
		let z_word_values = pack_mldsa44_z_y_coeffs(packed_y_coeff_values);
		let builder = CircuitBuilder::new();
		let z_words: Vec<_> = (0..z_word_values.len())
			.map(|_| builder.add_witness())
			.collect();
		assert_mldsa44_z_packed_bytes_norm(&builder, &z_words);

		let circuit = builder.build();
		let cs = circuit.constraint_system();
		let mut witness = circuit.new_witness_filler();
		for (&wire, &value) in z_words.iter().zip(z_word_values.iter()) {
			witness[wire] = Word(value);
		}
		if circuit.populate_wire_witness(&mut witness).is_err() {
			return false;
		}
		verify_constraints(cs, &witness.into_value_vec()).is_ok()
	}

	fn verify_mldsa44_z_decode_witness(packed_y_coeff_values: &[u64]) {
		let z_word_values = pack_mldsa44_z_y_coeffs(packed_y_coeff_values);
		let builder = CircuitBuilder::new();
		let z_words: Vec<_> = (0..z_word_values.len())
			.map(|_| builder.add_witness())
			.collect();
		let expected_coeffs: Vec<_> = (0..packed_y_coeff_values.len())
			.map(|_| builder.add_witness())
			.collect();

		let decoded = mldsa44_decode_z_packed_y(&builder, &z_words);
		for (i, (&decoded, &expected)) in decoded.iter().zip(expected_coeffs.iter()).enumerate() {
			builder.assert_eq(format!("mldsa44_z_decode[{i}]"), decoded, expected);
		}

		let circuit = builder.build();
		let cs = circuit.constraint_system();
		let mut witness = circuit.new_witness_filler();
		for (&wire, &value) in z_words.iter().zip(z_word_values.iter()) {
			witness[wire] = Word(value);
		}
		for (&wire, &value) in expected_coeffs.iter().zip(packed_y_coeff_values.iter()) {
			witness[wire] = Word(value);
		}

		circuit.populate_wire_witness(&mut witness).unwrap();
		verify_constraints(cs, &witness.into_value_vec())
			.expect("Circuit constraints should be satisfied");
	}

	fn verify_mldsa44_z_norm_witness(packed_y_coeff_values: &[u64]) -> bool {
		let builder = CircuitBuilder::new();
		let packed_y_coeffs: Vec<_> = (0..packed_y_coeff_values.len())
			.map(|_| builder.add_witness())
			.collect();
		assert_mldsa44_z_norm_from_packed_y(&builder, &packed_y_coeffs);

		let circuit = builder.build();
		let cs = circuit.constraint_system();
		let mut witness = circuit.new_witness_filler();
		for (&wire, &value) in packed_y_coeffs.iter().zip(packed_y_coeff_values.iter()) {
			witness[wire] = Word(value);
		}
		if circuit.populate_wire_witness(&mut witness).is_err() {
			return false;
		}
		verify_constraints(cs, &witness.into_value_vec()).is_ok()
	}

	#[test]
	fn mldsa44_final_challenge_hash_matches_shake256() {
		let mut rng = StdRng::seed_from_u64(0x4D4C4453413434);
		let mut input = vec![0u8; mldsa44::FINAL_CHALLENGE_INPUT_BYTES];
		rng.fill_bytes(&mut input);

		let mut hasher = Shake256::default();
		hasher.update(&input);
		let mut reader = hasher.finalize_xof();
		let mut expected = [0u8; mldsa44::C_TILDE_BYTES];
		reader.read(&mut expected);

		let builder = CircuitBuilder::new();
		let input_wires: Vec<_> = (0..input.len().div_ceil(8))
			.map(|_| builder.add_witness())
			.collect();
		let expected_wires: [Wire; mldsa44::C_TILDE_BYTES / 8] =
			std::array::from_fn(|_| builder.add_witness());

		let computed = mldsa44_final_challenge_hash(&builder, &input_wires);
		for i in 0..computed.len() {
			builder.assert_eq(format!("c_tilde_prime[{i}]"), computed[i], expected_wires[i]);
		}

		let circuit = builder.build();
		let cs = circuit.constraint_system();
		let mut witness = circuit.new_witness_filler();

		for (i, chunk) in input.chunks(8).enumerate() {
			let word = u64::from_le_bytes(chunk.try_into().unwrap());
			witness[input_wires[i]] = Word(word);
		}
		for (i, chunk) in expected.chunks(8).enumerate() {
			let word = u64::from_le_bytes(chunk.try_into().unwrap());
			witness[expected_wires[i]] = Word(word);
		}

		circuit.populate_wire_witness(&mut witness).unwrap();
		verify_constraints(cs, &witness.into_value_vec())
			.expect("Circuit constraints should be satisfied");
	}

	#[test]
	fn mldsa44_sample_in_ball_one_block_stream_matches_shake256() {
		let mut rng = StdRng::seed_from_u64(0x53414D504C453434);
		let mut c_tilde = [0u8; mldsa44::C_TILDE_BYTES];
		rng.fill_bytes(&mut c_tilde);

		let mut hasher = Shake256::default();
		hasher.update(&c_tilde);
		let mut reader = hasher.finalize_xof();
		let mut expected = [0u8; mldsa44::SAMPLE_IN_BALL_ONE_BLOCK_STREAM_BYTES];
		reader.read(&mut expected);

		let builder = CircuitBuilder::new();
		let c_tilde_wires: [Wire; mldsa44::C_TILDE_BYTES / 8] =
			std::array::from_fn(|_| builder.add_witness());
		let expected_wires: [Wire; mldsa44::SAMPLE_IN_BALL_ONE_BLOCK_STREAM_BYTES / 8] =
			std::array::from_fn(|_| builder.add_witness());

		let computed = mldsa44_sample_in_ball_one_block_stream(&builder, &c_tilde_wires);
		for i in 0..computed.len() {
			builder.assert_eq(
				format!("sample_in_ball_stream[{i}]"),
				computed[i],
				expected_wires[i],
			);
		}

		let circuit = builder.build();
		let cs = circuit.constraint_system();
		let mut witness = circuit.new_witness_filler();

		for (i, chunk) in c_tilde.chunks(8).enumerate() {
			let word = u64::from_le_bytes(chunk.try_into().unwrap());
			witness[c_tilde_wires[i]] = Word(word);
		}
		for (i, chunk) in expected.chunks(8).enumerate() {
			let word = u64::from_le_bytes(chunk.try_into().unwrap());
			witness[expected_wires[i]] = Word(word);
		}

		circuit.populate_wire_witness(&mut witness).unwrap();
		verify_constraints(cs, &witness.into_value_vec())
			.expect("Circuit constraints should be satisfied");
	}

	#[test]
	fn mldsa44_z_norm_accepts_boundary_values() {
		let mut packed_y = vec![mldsa44::Z_NORM_PACKED_Y_MIN; mldsa44::Z_COEFFICIENTS];
		packed_y[17] = mldsa44::Z_NORM_PACKED_Y_MAX;
		packed_y[999] = (mldsa44::Z_NORM_PACKED_Y_MIN + mldsa44::Z_NORM_PACKED_Y_MAX) / 2;

		assert!(verify_mldsa44_z_norm_witness(&packed_y));
	}

	#[test]
	fn mldsa44_decode_z_packed_y_matches_bitpack_layout() {
		let mut packed_y = vec![0u64; mldsa44::Z_COEFFICIENTS];
		for (i, coeff) in packed_y.iter_mut().enumerate() {
			*coeff = ((i as u64 * 65_537) + 0x12345) & ((1 << mldsa44::Z_BITS_PER_COEFF) - 1);
		}
		packed_y[0] = 0;
		packed_y[3] = (1 << mldsa44::Z_BITS_PER_COEFF) - 1;
		packed_y[7] = 0b10_1010_1111_0000_0101;
		packed_y[mldsa44::Z_COEFFICIENTS - 1] = 0x2_0001;

		verify_mldsa44_z_decode_witness(&packed_y);
	}

	#[test]
	fn mldsa44_z_packed_bytes_norm_accepts_boundary_values() {
		let mut packed_y = vec![mldsa44::Z_NORM_PACKED_Y_MIN; mldsa44::Z_COEFFICIENTS];
		packed_y[17] = mldsa44::Z_NORM_PACKED_Y_MAX;
		packed_y[999] = (mldsa44::Z_NORM_PACKED_Y_MIN + mldsa44::Z_NORM_PACKED_Y_MAX) / 2;

		assert!(verify_mldsa44_z_packed_bytes_norm_witness(&packed_y));
	}

	#[test]
	fn mldsa44_z_packed_bytes_norm_rejects_below_min() {
		let mut packed_y = vec![mldsa44::Z_NORM_PACKED_Y_MIN; mldsa44::Z_COEFFICIENTS];
		packed_y[41] = mldsa44::Z_NORM_PACKED_Y_MIN - 1;

		assert!(!verify_mldsa44_z_packed_bytes_norm_witness(&packed_y));
	}

	#[test]
	fn mldsa44_z_packed_bytes_norm_rejects_above_max() {
		let mut packed_y = vec![mldsa44::Z_NORM_PACKED_Y_MIN; mldsa44::Z_COEFFICIENTS];
		packed_y[271] = mldsa44::Z_NORM_PACKED_Y_MAX + 1;

		assert!(!verify_mldsa44_z_packed_bytes_norm_witness(&packed_y));
	}

	#[test]
	fn mldsa44_z_norm_rejects_below_min() {
		let mut packed_y = vec![mldsa44::Z_NORM_PACKED_Y_MIN; mldsa44::Z_COEFFICIENTS];
		packed_y[41] = mldsa44::Z_NORM_PACKED_Y_MIN - 1;

		assert!(!verify_mldsa44_z_norm_witness(&packed_y));
	}

	#[test]
	fn mldsa44_z_norm_rejects_above_max() {
		let mut packed_y = vec![mldsa44::Z_NORM_PACKED_Y_MIN; mldsa44::Z_COEFFICIENTS];
		packed_y[271] = mldsa44::Z_NORM_PACKED_Y_MAX + 1;

		assert!(!verify_mldsa44_z_norm_witness(&packed_y));
	}
}
