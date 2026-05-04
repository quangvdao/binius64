// Copyright 2025 Irreducible Inc.

use binius_frontend::{CircuitBuilder, Wire};

use crate::keccak::fixed_length;

use super::{
	mldsa44,
	types::{Mldsa44SampleInBallOneBlock, Mldsa44SampleInBallOneBlockSparse},
	util::{assert_true_cond, iadd_wrapping, isub_one, select_indexed_wire_unchecked},
};

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

/// Proves the fixed one-block ML-DSA-44 `SampleInBall` rejection/update relation from a SHAKE
/// stream.
///
/// The SHAKE stream must be the 136-byte output of `SHAKE256(c_tilde, 136)`. The helper allocates
/// private draw-count witness wires and enforces:
///
/// - each iteration consumes at least one draw byte;
/// - skipped draws are `> i`;
/// - the final accepted draw is `<= i`;
/// - challenge coefficients follow the FIPS swap/update rule.
pub fn mldsa44_sample_in_ball_one_block_from_stream(
	builder: &CircuitBuilder,
	stream: &[Wire],
) -> Mldsa44SampleInBallOneBlock {
	let sparse = mldsa44_sample_in_ball_one_block_sparse_from_stream(builder, stream);
	let coeffs = mldsa44_expand_sparse_sample_in_ball(builder, &sparse.positions, &sparse.signs);
	Mldsa44SampleInBallOneBlock {
		coeffs,
		draw_counts: sparse.draw_counts,
	}
}

/// Proves the fixed one-block ML-DSA-44 `SampleInBall` rejection/update relation from a SHAKE
/// stream and keeps the challenge in its natural sparse representation.
pub fn mldsa44_sample_in_ball_one_block_sparse_from_stream(
	builder: &CircuitBuilder,
	stream: &[Wire],
) -> Mldsa44SampleInBallOneBlockSparse {
	assert_eq!(
		stream.len(),
		mldsa44::SAMPLE_IN_BALL_ONE_BLOCK_STREAM_BYTES / 8,
		"ML-DSA-44 one-block SampleInBall expects 136 stream bytes packed into 17 words",
	);

	let zero = builder.add_constant_64(0);
	let one = builder.add_constant_64(1);
	let mut sparse_positions = Vec::with_capacity(mldsa44::TAU);
	let mut sparse_signs = Vec::with_capacity(mldsa44::TAU);
	let draw_counts: [Wire; mldsa44::TAU] = std::array::from_fn(|_| builder.add_witness());

	let mut draw_bytes = Vec::with_capacity(mldsa44::SAMPLE_IN_BALL_ONE_BLOCK_DRAW_CAP_BYTES);
	for &word in &stream[1..] {
		for byte_idx in 0..8 {
			draw_bytes.push(builder.extract_byte(word, byte_idx));
		}
	}

	let mut cursor = zero;
	for (round, &draw_count) in draw_counts.iter().enumerate() {
		let i = mldsa44::N - mldsa44::TAU + round;
		let i_wire = builder.add_constant_64(i as u64);

		builder.assert_true(
			format!("sample_in_ball_draw_count_nonzero[{round}]"),
			builder.icmp_ule(one, draw_count),
		);

		let next_cursor = iadd_wrapping(builder, cursor, draw_count);
		let draw_cap =
			builder.add_constant_64(mldsa44::SAMPLE_IN_BALL_ONE_BLOCK_DRAW_CAP_BYTES as u64);
		builder.assert_true(
			format!("sample_in_ball_cursor_monotone[{round}]"),
			builder.icmp_ule(cursor, next_cursor),
		);
		builder.assert_true(
			format!("sample_in_ball_cursor_within_cap[{round}]"),
			builder.icmp_ule(next_cursor, draw_cap),
		);

		let accepted_pos = isub_one(builder, next_cursor);
		let accepted_draw = select_indexed_wire_unchecked(builder, &draw_bytes, accepted_pos);
		builder.assert_true(
			format!("sample_in_ball_accepted_draw_le_i[{round}]"),
			builder.icmp_ule(accepted_draw, i_wire),
		);

		for (draw_idx, &draw) in draw_bytes.iter().enumerate() {
			let draw_idx_wire = builder.add_constant_64(draw_idx as u64);
			let after_cursor = builder.icmp_ule(cursor, draw_idx_wire);
			let before_accepted = builder.icmp_ult(draw_idx_wire, accepted_pos);
			let is_skipped_pos = builder.band(after_cursor, before_accepted);

			assert_true_cond(
				builder,
				format!("sample_in_ball_skipped_draw_gt_i[{round}][{draw_idx}]"),
				builder.icmp_ugt(draw, i_wire),
				is_skipped_pos,
			);
		}

		let sign_bit =
			builder.band(builder.shr(stream[0], round as u32), builder.add_constant_64(1));
		let sign_is_negative = builder.shl(sign_bit, 63);
		let signed_coeff = builder.select(sign_is_negative, builder.add_constant_64(u64::MAX), one);

		for position in sparse_positions.iter_mut() {
			let is_j = builder.icmp_eq(*position, accepted_draw);
			*position = builder.select(is_j, i_wire, *position);
		}
		sparse_positions.push(accepted_draw);
		sparse_signs.push(signed_coeff);

		cursor = next_cursor;
	}

	Mldsa44SampleInBallOneBlockSparse {
		positions: sparse_positions.try_into().unwrap(),
		signs: sparse_signs.try_into().unwrap(),
		draw_counts,
	}
}

pub(crate) fn mldsa44_expand_sparse_sample_in_ball(
	builder: &CircuitBuilder,
	positions: &[Wire; mldsa44::TAU],
	signs: &[Wire; mldsa44::TAU],
) -> [Wire; mldsa44::N] {
	let zero = builder.add_constant_64(0);
	let mut state = vec![zero; mldsa44::N];
	for (coeff_idx, coeff) in state.iter_mut().enumerate() {
		let coeff_idx_wire = builder.add_constant_64(coeff_idx as u64);
		for (&position, &sign) in positions.iter().zip(signs.iter()) {
			let is_position = builder.icmp_eq(position, coeff_idx_wire);
			*coeff = builder.select(is_position, sign, *coeff);
		}
	}

	state.try_into().unwrap()
}

/// Computes and proves the fixed one-block ML-DSA-44 `SampleInBall` relation.
pub fn mldsa44_sample_in_ball_one_block(
	builder: &CircuitBuilder,
	c_tilde: &[Wire],
) -> Mldsa44SampleInBallOneBlock {
	let stream = mldsa44_sample_in_ball_one_block_stream(builder, c_tilde);
	mldsa44_sample_in_ball_one_block_from_stream(builder, &stream)
}

/// Computes and proves the fixed one-block ML-DSA-44 `SampleInBall` relation, returning sparse
/// challenge entries for the lattice bridge.
pub fn mldsa44_sample_in_ball_one_block_sparse(
	builder: &CircuitBuilder,
	c_tilde: &[Wire],
) -> Mldsa44SampleInBallOneBlockSparse {
	let stream = mldsa44_sample_in_ball_one_block_stream(builder, c_tilde);
	mldsa44_sample_in_ball_one_block_sparse_from_stream(builder, &stream)
}
