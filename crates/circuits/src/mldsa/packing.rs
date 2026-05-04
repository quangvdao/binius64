// Copyright 2025 Irreducible Inc.

use binius_frontend::{CircuitBuilder, Wire};

use super::MldsaParams;

/// Encodes ML-DSA `w1_prime` coefficients as `w1Encode(w1_prime)`.
///
/// ML-DSA-44 has `gamma2 = (q - 1) / 88`, so each `w1_prime` coefficient is in `[0, 43]` and
/// `SimpleBitPack` stores each coefficient in six little-endian bits. The output is the fixed
/// 768-byte `w1Encode` stream packed into 96 little-endian 64-bit words.
pub fn encode_w1_for<P: MldsaParams>(builder: &CircuitBuilder, w1_coeffs: &[Wire]) -> Vec<Wire> {
	encode_w1_impl_for::<P>(builder, w1_coeffs, true)
}

pub(crate) fn encode_w1_unchecked_for<P: MldsaParams>(
	builder: &CircuitBuilder,
	w1_coeffs: &[Wire],
) -> Vec<Wire> {
	encode_w1_impl_for::<P>(builder, w1_coeffs, false)
}

fn encode_w1_impl_for<P: MldsaParams>(
	builder: &CircuitBuilder,
	w1_coeffs: &[Wire],
	check_range: bool,
) -> Vec<Wire> {
	assert_eq!(w1_coeffs.len(), P::W1_COEFFICIENTS, "{} w1 coefficient count mismatch", P::label(),);

	let max_coeff = builder.add_constant_64(P::W1_COEFF_MAX);
	let mut words = vec![builder.add_constant_64(0); P::W1_ENCODE_WORDS];

	for (coeff_idx, &coeff) in w1_coeffs.iter().enumerate() {
		if check_range {
			builder.assert_true(
				format!("{}_w1_coeff_range[{coeff_idx}]", P::label()),
				builder.icmp_ule(coeff, max_coeff),
			);
		}

		let bit_offset = coeff_idx * P::W1_BITS_PER_COEFF;
		let word_idx = bit_offset / 64;
		let shift = bit_offset % 64;

		if shift <= 64 - P::W1_BITS_PER_COEFF {
			words[word_idx] = builder.bor(words[word_idx], builder.shl(coeff, shift as u32));
		} else {
			let low_bits = 64 - shift;
			words[word_idx] = builder.bor(words[word_idx], builder.shl(coeff, shift as u32));
			words[word_idx + 1] =
				builder.bor(words[word_idx + 1], builder.shr(coeff, low_bits as u32));
		}
	}

	words
}

/// Enforces the ML-DSA `z` norm check on decoded BitPack coefficients.
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
pub fn assert_z_norm_from_packed_y_for<P: MldsaParams>(
	builder: &CircuitBuilder,
	packed_y_coeffs: &[Wire],
) {
	assert_eq!(
		packed_y_coeffs.len(),
		P::Z_COEFFICIENTS,
		"{} z coefficient count mismatch",
		P::label(),
	);

	let min = builder.add_constant_64(P::Z_NORM_PACKED_Y_MIN);
	let max = builder.add_constant_64(P::Z_NORM_PACKED_Y_MAX);

	for (i, &y) in packed_y_coeffs.iter().enumerate() {
		let gte_min = builder.icmp_ule(min, y);
		let lte_max = builder.icmp_ule(y, max);
		builder.assert_true(format!("{}_z_norm_gte_min[{i}]", P::label()), gte_min);
		builder.assert_true(format!("{}_z_norm_lte_max[{i}]", P::label()), lte_max);
	}
}

/// Decodes ML-DSA `z` from its fixed BitPack representation.
///
/// FIPS packs `z` coefficients as contiguous 18-bit little-endian integers:
///
/// ```text
/// y_i = gamma1 - z_i
/// ```
///
/// This helper keeps the signature hidden by taking the packed stream as witness wires. The output
/// is one 64-bit wire per decoded packed `y_i` coefficient.
pub fn decode_z_packed_y_for<P: MldsaParams>(
	builder: &CircuitBuilder,
	z_words: &[Wire],
) -> Vec<Wire> {
	assert_eq!(z_words.len(), P::Z_PACKED_WORDS, "{} packed z word count mismatch", P::label(),);

	let mask = builder.add_constant_64((1 << P::Z_BITS_PER_COEFF) - 1);
	let mut coeffs = Vec::with_capacity(P::Z_COEFFICIENTS);

	for coeff_idx in 0..P::Z_COEFFICIENTS {
		let bit_offset = coeff_idx * P::Z_BITS_PER_COEFF;
		let word_idx = bit_offset / 64;
		let shift = bit_offset % 64;

		let coeff = if shift <= 64 - P::Z_BITS_PER_COEFF {
			builder.band(builder.shr(z_words[word_idx], shift as u32), mask)
		} else {
			let low_bits = 64 - shift;
			let high_bits = P::Z_BITS_PER_COEFF - low_bits;
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

/// Decodes hidden ML-DSA packed `z` bytes and enforces `||z||_infty < gamma1 - beta`.
pub fn assert_z_packed_bytes_norm_for<P: MldsaParams>(
	builder: &CircuitBuilder,
	z_words: &[Wire],
) -> Vec<Wire> {
	let packed_y_coeffs = decode_z_packed_y_for::<P>(builder, z_words);
	assert_z_norm_from_packed_y_for::<P>(builder, &packed_y_coeffs);
	packed_y_coeffs
}
