// Copyright 2026 The Binius Developers

//! Upper-half sumcheck round-message accumulation for Keccak chi/iota constraints.

use std::iter;

use binius_field::{BinaryField, Field, PackedField};
use binius_utils::rayon::prelude::*;

use crate::{
	bit_ntt::{NttLookup, upper_half_residual_evals},
	constants::{LANE_BITS, N_LANES},
	trace::RoundTrace,
};

/// Accumulate Keccak chi/iota residual evaluations over round traces and lane weights.
///
/// Each weight corresponds to one `(round_trace, lane)` pair, in trace-major order. The returned
/// array is the 64-point upper-half extension-domain message.
pub fn upper_half_round_message<FChallenge, P>(
	lookup: &NttLookup<P>,
	round_traces: &[RoundTrace],
	eq_weights: &[FChallenge],
) -> [FChallenge; LANE_BITS]
where
	FChallenge: Field + From<P::Scalar>,
	P: PackedField,
	P::Scalar: BinaryField + Field,
{
	assert_eq!(eq_weights.len(), round_traces.len() * N_LANES);

	let mut acc = [FChallenge::ZERO; LANE_BITS];
	for (trace_idx, round_trace) in round_traces.iter().enumerate() {
		accumulate_round_trace(
			lookup,
			round_trace,
			trace_idx % crate::constants::N_ROUNDS,
			&eq_weights[trace_idx * N_LANES..(trace_idx + 1) * N_LANES],
			&mut acc,
		);
	}
	acc
}

/// Parallel variant of [`upper_half_round_message`].
pub fn par_upper_half_round_message<FChallenge, P>(
	lookup: &NttLookup<P>,
	round_traces: &[RoundTrace],
	eq_weights: &[FChallenge],
) -> [FChallenge; LANE_BITS]
where
	FChallenge: Field + From<P::Scalar> + Send + Sync,
	P: PackedField,
	P::Scalar: BinaryField + Field,
{
	assert_eq!(eq_weights.len(), round_traces.len() * N_LANES);

	round_traces
		.par_iter()
		.enumerate()
		.map(|(trace_idx, round_trace)| {
			let mut acc = [FChallenge::ZERO; LANE_BITS];
			accumulate_round_trace(
				lookup,
				round_trace,
				trace_idx % crate::constants::N_ROUNDS,
				&eq_weights[trace_idx * N_LANES..(trace_idx + 1) * N_LANES],
				&mut acc,
			);
			acc
		})
		.reduce(
			|| [FChallenge::ZERO; LANE_BITS],
			|mut lhs, rhs| {
				for (lhs_i, rhs_i) in iter::zip(&mut lhs, rhs) {
					*lhs_i += rhs_i;
				}
				lhs
			},
		)
}

/// Accumulate a round message with small-field lane weights, staying packed in the hot loop.
pub fn upper_half_round_message_small_weights<FChallenge, P>(
	lookup: &NttLookup<P>,
	round_traces: &[RoundTrace],
	eq_weights: &[P::Scalar],
) -> [FChallenge; LANE_BITS]
where
	FChallenge: Field + From<P::Scalar>,
	P: PackedField,
	P::Scalar: BinaryField + Field,
{
	assert_eq!(eq_weights.len(), round_traces.len() * N_LANES);

	let mut acc = [P::zero(); crate::bit_ntt::PACKED_EVALS];
	for (trace_idx, round_trace) in round_traces.iter().enumerate() {
		accumulate_round_trace_small_weights(
			lookup,
			round_trace,
			trace_idx % crate::constants::N_ROUNDS,
			&eq_weights[trace_idx * N_LANES..(trace_idx + 1) * N_LANES],
			&mut acc,
		);
	}

	std::array::from_fn(|i| FChallenge::from(P::iter_slice(&acc).nth(i).unwrap()))
}

/// Parallel variant of [`upper_half_round_message_small_weights`].
pub fn par_upper_half_round_message_small_weights<FChallenge, P>(
	lookup: &NttLookup<P>,
	round_traces: &[RoundTrace],
	eq_weights: &[P::Scalar],
) -> [FChallenge; LANE_BITS]
where
	FChallenge: Field + From<P::Scalar> + Send + Sync,
	P: PackedField,
	P::Scalar: BinaryField + Field,
{
	assert_eq!(eq_weights.len(), round_traces.len() * N_LANES);

	let packed_acc = round_traces
		.par_iter()
		.enumerate()
		.map(|(trace_idx, round_trace)| {
			let mut acc = [P::zero(); crate::bit_ntt::PACKED_EVALS];
			accumulate_round_trace_small_weights(
				lookup,
				round_trace,
				trace_idx % crate::constants::N_ROUNDS,
				&eq_weights[trace_idx * N_LANES..(trace_idx + 1) * N_LANES],
				&mut acc,
			);
			acc
		})
		.reduce(
			|| [P::zero(); crate::bit_ntt::PACKED_EVALS],
			|mut lhs, rhs| {
				for (lhs_i, rhs_i) in iter::zip(&mut lhs, rhs) {
					*lhs_i += rhs_i;
				}
				lhs
			},
		);

	std::array::from_fn(|i| FChallenge::from(P::iter_slice(&packed_acc).nth(i).unwrap()))
}

fn accumulate_round_trace<FChallenge, P>(
	lookup: &NttLookup<P>,
	round_trace: &RoundTrace,
	round: usize,
	eq_weights: &[FChallenge],
	acc: &mut [FChallenge; LANE_BITS],
) where
	FChallenge: Field + From<P::Scalar>,
	P: PackedField,
	P::Scalar: BinaryField + Field,
{
	debug_assert_eq!(eq_weights.len(), N_LANES);

	let residuals =
		upper_half_residual_evals::<P>(lookup, &round_trace.pre_chi, &round_trace.output, round);
	for (lane_residuals, &eq_weight) in iter::zip(&residuals, eq_weights) {
		for (acc_i, residual_i) in iter::zip(&mut *acc, P::iter_slice(lane_residuals)) {
			*acc_i += eq_weight * FChallenge::from(residual_i);
		}
	}
}

fn accumulate_round_trace_small_weights<P>(
	lookup: &NttLookup<P>,
	round_trace: &RoundTrace,
	round: usize,
	eq_weights: &[P::Scalar],
	acc: &mut [P; crate::bit_ntt::PACKED_EVALS],
) where
	P: PackedField,
	P::Scalar: BinaryField + Field,
{
	debug_assert_eq!(eq_weights.len(), N_LANES);

	let residuals =
		upper_half_residual_evals::<P>(lookup, &round_trace.pre_chi, &round_trace.output, round);
	for (lane_residuals, &eq_weight) in iter::zip(&residuals, eq_weights) {
		let eq_weight = P::broadcast(eq_weight);
		for (acc_i, residual_i) in iter::zip(&mut *acc, lane_residuals) {
			*acc_i += *residual_i * eq_weight;
		}
	}
}

#[cfg(test)]
mod tests {
	use binius_field::{AESTowerField8b, BinaryField128bGhash, PackedAESBinaryField16x8b, Random};
	use rand::{Rng, SeedableRng, rngs::StdRng};

	use crate::{
		bit_ntt::NttLookup,
		constants::N_ROUNDS,
		trace::{PermutationTrace, State},
	};

	use super::*;

	type B128 = BinaryField128bGhash;
	type P = PackedAESBinaryField16x8b;

	fn round_traces(trace: &PermutationTrace) -> Vec<RoundTrace> {
		trace.rounds.to_vec()
	}

	#[test]
	fn valid_keccak_trace_parallel_matches_sequential() {
		let mut rng = StdRng::seed_from_u64(12);
		let lookup = NttLookup::<P>::for_upper_half_domain();
		let trace = PermutationTrace::new(rng.random::<State>());
		let round_traces = round_traces(&trace);
		let eq_weights: Vec<_> = (0..round_traces.len() * N_LANES)
			.map(|_| B128::random(&mut rng))
			.collect();

		assert_eq!(
			par_upper_half_round_message(&lookup, &round_traces, &eq_weights),
			upper_half_round_message(&lookup, &round_traces, &eq_weights)
		);
	}

	#[test]
	fn parallel_round_message_matches_sequential() {
		let mut rng = StdRng::seed_from_u64(13);
		let lookup = NttLookup::<P>::for_upper_half_domain();
		let mut trace = PermutationTrace::new(rng.random::<State>());
		trace.rounds[7].output[3] ^= 0x0123_4567_89ab_cdef;

		let round_traces = round_traces(&trace);
		let eq_weights: Vec<_> = (0..N_ROUNDS * N_LANES)
			.map(|_| B128::random(&mut rng))
			.collect();

		assert_eq!(
			par_upper_half_round_message(&lookup, &round_traces, &eq_weights),
			upper_half_round_message(&lookup, &round_traces, &eq_weights)
		);
	}

	#[test]
	fn small_weight_round_message_matches_generic_round_message() {
		let mut rng = StdRng::seed_from_u64(14);
		let lookup = NttLookup::<P>::for_upper_half_domain();
		let mut trace = PermutationTrace::new(rng.random::<State>());
		trace.rounds[11].output[9] ^= 0xfedc_ba98_7654_3210;

		let round_traces = round_traces(&trace);
		let small_weights: Vec<_> = (0..N_ROUNDS * N_LANES)
			.map(|_| rng.random::<AESTowerField8b>())
			.collect();
		let big_weights: Vec<_> = small_weights.iter().copied().map(B128::from).collect();

		assert_eq!(
			upper_half_round_message_small_weights::<B128, P>(
				&lookup,
				&round_traces,
				&small_weights
			),
			upper_half_round_message(&lookup, &round_traces, &big_weights)
		);
		assert_eq!(
			par_upper_half_round_message_small_weights::<B128, P>(
				&lookup,
				&round_traces,
				&small_weights
			),
			upper_half_round_message(&lookup, &round_traces, &big_weights)
		);
	}
}
