// Copyright 2025 Irreducible Inc.

use binius_frontend::{CircuitBuilder, Wire};

use crate::keccak::fixed_length;

use super::{
	MldsaParams,
	types::{MldsaSampleInBallFixedCap, MldsaSampleInBallFixedCapSparse},
	util::{assert_true_cond, iadd_wrapping, isub_one, select_indexed_wire_unchecked},
};

/// Computes the fixed-cap ML-DSA `SampleInBall` SHAKE stream.
///
/// This is the first fixed-cap sampler shape from the top-level plan:
///
/// ```text
/// stream = SHAKE256(c_tilde, sign_bytes + draw_cap_bytes)
/// sign bytes = stream[0..8]
/// position draws = stream[8..]
/// ```
///
/// The Fisher-Yates rejection/update relation is layered on top of this stream. The cap is part of
/// the circuit identity, so repeated batches must not mix different cap choices.
pub fn sample_in_ball_fixed_cap_stream_for<P: MldsaParams>(
	builder: &CircuitBuilder,
	c_tilde: &[Wire],
) -> Vec<Wire> {
	assert_eq!(
		c_tilde.len(),
		P::C_TILDE_BYTES / 8,
		"{} c_tilde packed word count mismatch",
		P::label(),
	);

	fixed_length::shake256(builder, c_tilde, P::C_TILDE_BYTES, P::SAMPLE_IN_BALL_STREAM_BYTES)
}

/// Proves the fixed-cap ML-DSA `SampleInBall` rejection/update relation from a SHAKE
/// stream.
///
/// The SHAKE stream must be the fixed-cap output of `SHAKE256(c_tilde, stream_bytes)`. The helper allocates
/// private draw-count witness wires and enforces:
///
/// - each iteration consumes at least one draw byte;
/// - skipped draws are `> i`;
/// - the final accepted draw is `<= i`;
/// - challenge coefficients follow the FIPS swap/update rule.
pub fn sample_in_ball_fixed_cap_from_stream_for<P: MldsaParams>(
	builder: &CircuitBuilder,
	stream: &[Wire],
) -> MldsaSampleInBallFixedCap {
	let sparse = sample_in_ball_fixed_cap_sparse_from_stream_for::<P>(builder, stream);
	let coeffs = expand_sparse_sample_in_ball_for::<P>(builder, &sparse.positions, &sparse.signs);
	MldsaSampleInBallFixedCap {
		coeffs,
		draw_counts: sparse.draw_counts,
	}
}

/// Proves the fixed-cap ML-DSA `SampleInBall` rejection/update relation from a SHAKE
/// stream and keeps the challenge in its natural sparse representation.
pub fn sample_in_ball_fixed_cap_sparse_from_stream_for<P: MldsaParams>(
	builder: &CircuitBuilder,
	stream: &[Wire],
) -> MldsaSampleInBallFixedCapSparse {
	assert_eq!(
		stream.len(),
		P::SAMPLE_IN_BALL_STREAM_BYTES.div_ceil(8),
		"{} fixed-cap SampleInBall packed stream length mismatch",
		P::label(),
	);

	let zero = builder.add_constant_64(0);
	let one = builder.add_constant_64(1);
	let mut sparse_positions = Vec::with_capacity(P::TAU);
	let mut sparse_signs = Vec::with_capacity(P::TAU);
	let draw_counts: Vec<Wire> = (0..P::TAU).map(|_| builder.add_witness()).collect();

	let mut draw_bytes = Vec::with_capacity(P::SAMPLE_IN_BALL_DRAW_CAP_BYTES);
	for &word in &stream[1..] {
		for byte_idx in 0..8 {
			if draw_bytes.len() < P::SAMPLE_IN_BALL_DRAW_CAP_BYTES {
				draw_bytes.push(builder.extract_byte(word, byte_idx));
			}
		}
	}
	assert_eq!(draw_bytes.len(), P::SAMPLE_IN_BALL_DRAW_CAP_BYTES);

	let mut cursor = zero;
	for (round, &draw_count) in draw_counts.iter().enumerate() {
		let i = P::N - P::TAU + round;
		let i_wire = builder.add_constant_64(i as u64);

		builder.assert_true(
			format!("{}_sample_in_ball_draw_count_nonzero[{round}]", P::label()),
			builder.icmp_ule(one, draw_count),
		);

		let next_cursor = iadd_wrapping(builder, cursor, draw_count);
		let draw_cap = builder.add_constant_64(P::SAMPLE_IN_BALL_DRAW_CAP_BYTES as u64);
		builder.assert_true(
			format!("{}_sample_in_ball_cursor_monotone[{round}]", P::label()),
			builder.icmp_ule(cursor, next_cursor),
		);
		builder.assert_true(
			format!("{}_sample_in_ball_cursor_within_cap[{round}]", P::label()),
			builder.icmp_ule(next_cursor, draw_cap),
		);

		let accepted_pos = isub_one(builder, next_cursor);
		let accepted_draw = select_indexed_wire_unchecked(builder, &draw_bytes, accepted_pos);
		builder.assert_true(
			format!("{}_sample_in_ball_accepted_draw_le_i[{round}]", P::label()),
			builder.icmp_ule(accepted_draw, i_wire),
		);

		let min_possible_draw_idx = round;
		let max_possible_draw_idx = P::SAMPLE_IN_BALL_DRAW_CAP_BYTES - P::TAU + round;
		for (draw_idx, &draw) in draw_bytes
			.iter()
			.enumerate()
			.skip(min_possible_draw_idx)
			.take(max_possible_draw_idx - min_possible_draw_idx + 1)
		{
			let draw_idx_wire = builder.add_constant_64(draw_idx as u64);
			let after_cursor = builder.icmp_ule(cursor, draw_idx_wire);
			let before_accepted = builder.icmp_ult(draw_idx_wire, accepted_pos);
			let is_skipped_pos = builder.band(after_cursor, before_accepted);

			assert_true_cond(
				builder,
				format!("{}_sample_in_ball_skipped_draw_gt_i[{round}][{draw_idx}]", P::label()),
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

	MldsaSampleInBallFixedCapSparse {
		positions: sparse_positions,
		signs: sparse_signs,
		draw_counts,
	}
}

pub(crate) fn expand_sparse_sample_in_ball_for<P: MldsaParams>(
	builder: &CircuitBuilder,
	positions: &[Wire],
	signs: &[Wire],
) -> Vec<Wire> {
	assert_eq!(positions.len(), P::TAU, "{} sparse position count mismatch", P::label());
	assert_eq!(signs.len(), P::TAU, "{} sparse sign count mismatch", P::label());
	let zero = builder.add_constant_64(0);
	let mut state = vec![zero; P::N];
	for (coeff_idx, coeff) in state.iter_mut().enumerate() {
		let coeff_idx_wire = builder.add_constant_64(coeff_idx as u64);
		for (&position, &sign) in positions.iter().zip(signs.iter()) {
			let is_position = builder.icmp_eq(position, coeff_idx_wire);
			*coeff = builder.select(is_position, sign, *coeff);
		}
	}

	state
}

/// Computes and proves the fixed-cap ML-DSA `SampleInBall` relation.
pub fn sample_in_ball_fixed_cap_for<P: MldsaParams>(
	builder: &CircuitBuilder,
	c_tilde: &[Wire],
) -> MldsaSampleInBallFixedCap {
	let stream = sample_in_ball_fixed_cap_stream_for::<P>(builder, c_tilde);
	sample_in_ball_fixed_cap_from_stream_for::<P>(builder, &stream)
}

/// Computes and proves the fixed-cap ML-DSA `SampleInBall` relation, returning sparse
/// challenge entries for the lattice bridge.
pub fn sample_in_ball_fixed_cap_sparse_for<P: MldsaParams>(
	builder: &CircuitBuilder,
	c_tilde: &[Wire],
) -> MldsaSampleInBallFixedCapSparse {
	let stream = sample_in_ball_fixed_cap_stream_for::<P>(builder, c_tilde);
	sample_in_ball_fixed_cap_sparse_from_stream_for::<P>(builder, &stream)
}
