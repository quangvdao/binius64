// Copyright 2025 Irreducible Inc.

use binius_frontend::{CircuitBuilder, Wire};

use super::{
	mldsa44,
	util::{assert_true_cond, select_indexed_wire_unchecked, unpack_bytes_from_words},
};

/// Asserts each `h_i in {0, 1}` and that the total Hamming weight is at most `omega`.
///
/// Implemented as a running sum guarded by per-step monotonicity, so a malicious witness cannot
/// hide overflow inside the accumulator.
pub fn assert_mldsa44_h_bits_and_weight(builder: &CircuitBuilder, h_coeffs: &[Wire]) {
	assert_eq!(
		h_coeffs.len(),
		mldsa44::W1_COEFFICIENTS,
		"ML-DSA-44 h has K * 256 = 1024 coefficients",
	);

	let one = builder.add_constant_64(1);
	let mut h_weight = builder.add_constant_64(0);
	for (i, &h) in h_coeffs.iter().enumerate() {
		builder.assert_true(format!("mldsa44_h_bit[{i}]"), builder.icmp_ule(h, one));
		let (next_weight, _carry) = builder.iadd(h_weight, h);
		builder.assert_true(
			format!("mldsa44_h_weight_monotone[{i}]"),
			builder.icmp_ule(h_weight, next_weight),
		);
		h_weight = next_weight;
	}

	builder.assert_true(
		"mldsa44_h_weight_omega",
		builder.icmp_ule(h_weight, builder.add_constant_64(mldsa44::OMEGA)),
	);
}

fn assert_mldsa44_h_bits_and_exact_weight(
	builder: &CircuitBuilder,
	h_coeffs: &[Wire],
	weight: Wire,
) {
	assert_eq!(
		h_coeffs.len(),
		mldsa44::W1_COEFFICIENTS,
		"ML-DSA-44 h has K * 256 = 1024 coefficients",
	);

	let mut h_weight = builder.add_constant_64(0);
	for (i, &h) in h_coeffs.iter().enumerate() {
		let (next_weight, _carry) = builder.iadd(h_weight, h);
		builder.assert_true(
			format!("mldsa44_h_weight_exact_monotone[{i}]"),
			builder.icmp_ule(h_weight, next_weight),
		);
		h_weight = next_weight;
	}

	builder.assert_eq("mldsa44_h_weight_exact", h_weight, weight);
}

/// Canonically decodes the ML-DSA-44 compressed hint representation into expanded hint bits.
///
/// The 84-byte ML-DSA-44 hint encoding stores up to `omega = 80` position bytes followed by `k = 4`
/// endpoint bytes. This enforces:
///
/// - monotone endpoints in `[0, omega]`;
/// - strictly increasing positions inside each polynomial segment;
/// - zero unused position bytes after the final endpoint;
/// - decoded expanded bits match exactly the compressed positions.
pub fn mldsa44_decode_hint_canonical(builder: &CircuitBuilder, h_words: &[Wire]) -> Vec<Wire> {
	let h_bytes = unpack_bytes_from_words(builder, h_words, mldsa44::HINT_BYTES);
	let zero = builder.add_constant_64(0);
	let one = builder.add_constant_64(1);
	let omega = builder.add_constant_64(mldsa44::OMEGA);
	let poly_n = builder.add_constant_64(mldsa44::N as u64);

	let mut endpoints = Vec::with_capacity(mldsa44::K);
	let mut prev_endpoint = zero;
	for poly_idx in 0..mldsa44::K {
		let endpoint = h_bytes[mldsa44::OMEGA_USIZE + poly_idx];
		builder.assert_true(
			format!("mldsa44_hint_endpoint_monotone[{poly_idx}]"),
			builder.icmp_ule(prev_endpoint, endpoint),
		);
		builder.assert_true(
			format!("mldsa44_hint_endpoint_omega[{poly_idx}]"),
			builder.icmp_ule(endpoint, omega),
		);
		endpoints.push(endpoint);
		prev_endpoint = endpoint;
	}

	let mut segment_starts = Vec::with_capacity(mldsa44::K);
	segment_starts.push(zero);
	for poly_idx in 1..mldsa44::K {
		segment_starts.push(endpoints[poly_idx - 1]);
	}

	for pos_idx in 0..mldsa44::OMEGA_USIZE {
		let pos_idx_wire = builder.add_constant_64(pos_idx as u64);
		let used = builder.icmp_ult(pos_idx_wire, endpoints[mldsa44::K - 1]);
		assert_true_cond(
			builder,
			format!("mldsa44_hint_unused_zero[{pos_idx}]"),
			builder.icmp_eq(h_bytes[pos_idx], zero),
			builder.bnot(used),
		);
		assert_true_cond(
			builder,
			format!("mldsa44_hint_position_lt_n[{pos_idx}]"),
			builder.icmp_ult(h_bytes[pos_idx], poly_n),
			used,
		);
	}

	for poly_idx in 0..mldsa44::K {
		let start = segment_starts[poly_idx];
		let end = endpoints[poly_idx];
		for pos_idx in 0..mldsa44::OMEGA_USIZE {
			let pos_idx_wire = builder.add_constant_64(pos_idx as u64);
			let in_segment = builder
				.band(builder.icmp_ule(start, pos_idx_wire), builder.icmp_ult(pos_idx_wire, end));
			let has_previous = builder.icmp_ult(start, pos_idx_wire);
			let needs_order_check = builder.band(in_segment, has_previous);
			let prev_pos = h_bytes[pos_idx - usize::from(pos_idx > 0)];
			assert_true_cond(
				builder,
				format!("mldsa44_hint_strictly_increasing[{poly_idx}][{pos_idx}]"),
				builder.icmp_ult(prev_pos, h_bytes[pos_idx]),
				needs_order_check,
			);
		}
	}

	let mut h_coeffs = Vec::with_capacity(mldsa44::W1_COEFFICIENTS);
	for poly_idx in 0..mldsa44::K {
		let start = segment_starts[poly_idx];
		let end = endpoints[poly_idx];
		for coeff_idx in 0..mldsa44::N {
			let coeff_idx_wire = builder.add_constant_64(coeff_idx as u64);
			let mut coeff = zero;
			for (pos_idx, &pos_byte) in h_bytes[..mldsa44::OMEGA_USIZE].iter().enumerate() {
				let pos_idx_wire = builder.add_constant_64(pos_idx as u64);
				let in_segment = builder.band(
					builder.icmp_ule(start, pos_idx_wire),
					builder.icmp_ult(pos_idx_wire, end),
				);
				let position_matches = builder.icmp_eq(pos_byte, coeff_idx_wire);
				let is_set = builder.band(in_segment, position_matches);
				coeff = builder.select(is_set, one, coeff);
			}
			h_coeffs.push(coeff);
		}
	}

	h_coeffs
}

/// Enforces canonical ML-DSA-44 compressed hint decoding against expanded hidden hint bits.
///
/// This is semantically equivalent to [`mldsa44_decode_hint_canonical`] followed by equality to an
/// expanded witness, but is much cheaper for the full relation: the expanded `h` bits are supplied
/// redundantly as hidden witness and the compressed encoding proves their canonical origin.
pub fn assert_mldsa44_hint_canonical_matches_expanded(
	builder: &CircuitBuilder,
	h_words: &[Wire],
	h_coeffs: &[Wire],
) {
	assert_eq!(
		h_coeffs.len(),
		mldsa44::W1_COEFFICIENTS,
		"ML-DSA-44 h has K * 256 = 1024 coefficients",
	);

	let h_bytes = unpack_bytes_from_words(builder, h_words, mldsa44::HINT_BYTES);
	let zero = builder.add_constant_64(0);
	let one = builder.add_constant_64(1);
	let omega = builder.add_constant_64(mldsa44::OMEGA);
	let poly_n = builder.add_constant_64(mldsa44::N as u64);

	let mut endpoints = Vec::with_capacity(mldsa44::K);
	let mut prev_endpoint = zero;
	for poly_idx in 0..mldsa44::K {
		let endpoint = h_bytes[mldsa44::OMEGA_USIZE + poly_idx];
		builder.assert_true(
			format!("mldsa44_hint_match_endpoint_monotone[{poly_idx}]"),
			builder.icmp_ule(prev_endpoint, endpoint),
		);
		builder.assert_true(
			format!("mldsa44_hint_match_endpoint_omega[{poly_idx}]"),
			builder.icmp_ule(endpoint, omega),
		);
		endpoints.push(endpoint);
		prev_endpoint = endpoint;
	}

	let mut segment_starts = Vec::with_capacity(mldsa44::K);
	segment_starts.push(zero);
	for poly_idx in 1..mldsa44::K {
		segment_starts.push(endpoints[poly_idx - 1]);
	}

	for pos_idx in 0..mldsa44::OMEGA_USIZE {
		let pos_idx_wire = builder.add_constant_64(pos_idx as u64);
		let used = builder.icmp_ult(pos_idx_wire, endpoints[mldsa44::K - 1]);
		assert_true_cond(
			builder,
			format!("mldsa44_hint_match_unused_zero[{pos_idx}]"),
			builder.icmp_eq(h_bytes[pos_idx], zero),
			builder.bnot(used),
		);
		assert_true_cond(
			builder,
			format!("mldsa44_hint_match_position_lt_n[{pos_idx}]"),
			builder.icmp_ult(h_bytes[pos_idx], poly_n),
			used,
		);
	}

	for poly_idx in 0..mldsa44::K {
		let start = segment_starts[poly_idx];
		let end = endpoints[poly_idx];
		let poly_h = &h_coeffs[poly_idx * mldsa44::N..(poly_idx + 1) * mldsa44::N];
		for pos_idx in 0..mldsa44::OMEGA_USIZE {
			let pos_idx_wire = builder.add_constant_64(pos_idx as u64);
			let in_segment = builder
				.band(builder.icmp_ule(start, pos_idx_wire), builder.icmp_ult(pos_idx_wire, end));
			let has_previous = builder.icmp_ult(start, pos_idx_wire);
			let needs_order_check = builder.band(in_segment, has_previous);
			let prev_pos = h_bytes[pos_idx - usize::from(pos_idx > 0)];
			assert_true_cond(
				builder,
				format!("mldsa44_hint_match_strictly_increasing[{poly_idx}][{pos_idx}]"),
				builder.icmp_ult(prev_pos, h_bytes[pos_idx]),
				needs_order_check,
			);

			let selected_h = select_indexed_wire_unchecked(builder, poly_h, h_bytes[pos_idx]);
			assert_true_cond(
				builder,
				format!("mldsa44_hint_match_sets_h[{poly_idx}][{pos_idx}]"),
				builder.icmp_eq(selected_h, one),
				in_segment,
			);
		}
	}

	assert_mldsa44_h_bits_and_exact_weight(builder, h_coeffs, endpoints[mldsa44::K - 1]);
}
