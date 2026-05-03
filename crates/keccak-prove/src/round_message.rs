// Copyright 2026 The Binius Developers

//! Upper-half sumcheck round-message accumulation for Keccak chi/iota constraints.

use std::iter;

use binius_field::{BinaryField, Field, PackedField};
use binius_math::{
	BinarySubspace,
	univariate::{extrapolate_over_subspace, lagrange_evals_scalars},
};
use binius_utils::rayon::prelude::*;

use crate::{
	bit_ntt::{CHI_OPERAND_LANES, NttLookup, upper_half_residual_evals},
	constants::{LANE_BITS, N_LANES, ROUND_CONSTANTS},
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

/// Build the full first-round message evaluations expected by the verifier.
///
/// The first half is the original 64-point bit domain, where valid residuals are zero. The second
/// half is the shifted upper-half extension domain sent by the prover.
pub fn first_round_message_evals<FChallenge>(
	upper_half_message: &[FChallenge; LANE_BITS],
) -> Vec<FChallenge>
where
	FChallenge: Field,
{
	let mut message_evals = vec![FChallenge::ZERO; 2 * LANE_BITS];
	message_evals[LANE_BITS..2 * LANE_BITS].copy_from_slice(upper_half_message);
	message_evals
}

/// Extrapolate the first-round message at the verifier's univariate challenge.
pub fn next_sum_claim_from_upper_half<FChallenge>(
	upper_half_message: &[FChallenge; LANE_BITS],
	challenge: FChallenge,
	prover_message_domain: &BinarySubspace<FChallenge>,
) -> FChallenge
where
	FChallenge: BinaryField + Field,
{
	let first_round_message_evals = first_round_message_evals(upper_half_message);

	extrapolate_over_subspace(&prover_message_domain, &first_round_message_evals, challenge)
}

/// Production-shaped first-round flow for Keccak chi/iota constraints.
pub fn par_first_round_claim_small_weights<FChallenge, P>(
	lookup: &NttLookup<P>,
	round_traces: &[RoundTrace],
	eq_weights: &[P::Scalar],
	challenge: FChallenge,
) -> FChallenge
where
	FChallenge: BinaryField + Field + From<P::Scalar> + Send + Sync,
	P: PackedField,
	P::Scalar: BinaryField + Field,
{
	let upper_half_message =
		par_upper_half_round_message_small_weights(lookup, round_traces, eq_weights);
	let prover_message_domain =
		BinarySubspace::<P::Scalar>::with_dim(crate::bit_ntt::LOG_LANE_BITS + 1)
			.isomorphic::<FChallenge>();
	next_sum_claim_from_upper_half(&upper_half_message, challenge, &prover_message_domain)
}

/// Compute the post-first-challenge claim by directly folding each lane word.
///
/// This is a correctness-oriented bridge toward the remaining sumcheck rounds. It mirrors the
/// verifier's view after the first univariate challenge: the 64-bit lane axis has been folded to
/// the challenge, leaving the outer round/lane indexed relation.
pub fn folded_next_sum_claim_small_weights<FChallenge, P>(
	round_traces: &[RoundTrace],
	eq_weights: &[P::Scalar],
	challenge: FChallenge,
) -> FChallenge
where
	FChallenge: BinaryField + Field + From<P::Scalar>,
	P: PackedField,
	P::Scalar: BinaryField + Field,
{
	assert_eq!(eq_weights.len(), round_traces.len() * N_LANES);

	let input_domain = BinarySubspace::<P::Scalar>::with_dim(crate::bit_ntt::LOG_LANE_BITS + 1)
		.reduce_dim(crate::bit_ntt::LOG_LANE_BITS)
		.isomorphic::<FChallenge>();
	let lagrange_evals = lagrange_evals_scalars(&input_domain, challenge);

	let mut claim = FChallenge::ZERO;
	for (trace_idx, round_trace) in round_traces.iter().enumerate() {
		let round = trace_idx % crate::constants::N_ROUNDS;
		let weights = &eq_weights[trace_idx * N_LANES..(trace_idx + 1) * N_LANES];
		for (lane_idx, &weight) in weights.iter().enumerate() {
			let (p_lane, q_lane, r_lane) = CHI_OPERAND_LANES[lane_idx];
			let p = fold_word(!round_trace.pre_chi[p_lane], &lagrange_evals);
			let q = fold_word(round_trace.pre_chi[q_lane], &lagrange_evals);
			let r = fold_word(round_trace.pre_chi[r_lane], &lagrange_evals);
			let next = fold_word(round_trace.output[lane_idx], &lagrange_evals);
			let iota = fold_word(
				if lane_idx == 0 {
					ROUND_CONSTANTS[round]
				} else {
					0
				},
				&lagrange_evals,
			);
			claim += (p * q - r - next - iota) * FChallenge::from(weight);
		}
	}

	claim
}

fn fold_word<F>(word: u64, lagrange_evals: &[F]) -> F
where
	F: Field,
{
	lagrange_evals
		.iter()
		.enumerate()
		.filter_map(|(bit_idx, &eval)| (((word >> bit_idx) & 1) == 1).then_some(eval))
		.sum()
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

	#[test]
	fn first_round_claim_matches_manual_extrapolation() {
		let mut rng = StdRng::seed_from_u64(15);
		let lookup = NttLookup::<P>::for_upper_half_domain();
		let trace = PermutationTrace::new(rng.random::<State>());

		let round_traces = round_traces(&trace);
		let small_weights: Vec<_> = (0..N_ROUNDS * N_LANES)
			.map(|_| rng.random::<AESTowerField8b>())
			.collect();
		let challenge = B128::random(&mut rng);
		let prover_message_domain =
			BinarySubspace::<AESTowerField8b>::with_dim(crate::bit_ntt::LOG_LANE_BITS + 1)
				.isomorphic::<B128>();
		let upper_half = par_upper_half_round_message_small_weights::<B128, P>(
			&lookup,
			&round_traces,
			&small_weights,
		);
		let message_evals = first_round_message_evals(&upper_half);
		assert_eq!(message_evals.len(), 2 * LANE_BITS);
		assert_eq!(&message_evals[..LANE_BITS], [B128::ZERO; LANE_BITS]);
		assert_eq!(&message_evals[LANE_BITS..], upper_half);

		assert_eq!(
			par_first_round_claim_small_weights::<B128, P>(
				&lookup,
				&round_traces,
				&small_weights,
				challenge
			),
			next_sum_claim_from_upper_half(&upper_half, challenge, &prover_message_domain)
		);
		assert_eq!(
			folded_next_sum_claim_small_weights::<B128, P>(
				&round_traces,
				&small_weights,
				challenge
			),
			next_sum_claim_from_upper_half(&upper_half, challenge, &prover_message_domain)
		);
	}
}
