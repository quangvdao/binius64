// Copyright 2026 The Binius Developers

//! Upper-half sumcheck round-message accumulation for Keccak chi/iota constraints.

use std::{cell::UnsafeCell, iter, sync::Barrier, thread};

use binius_field::{
	BinaryField, Field, PackedField, UnderlierWithBitOps, WithUnderlier,
	linear_transformation::{
		BytewiseLookupTransformationFactory, LinearTransformationFactory,
		OutputWrappingTransformationFactory, Transformation,
	},
};
use binius_ip::{channel::IPVerifierChannel, mlecheck, sumcheck::RoundCoeffs};
use binius_ip_prover::channel::IPProverChannel;
use binius_ip_prover::sumcheck::{
	Error as SumcheckError, common::SumcheckProver, prove_single_mlecheck,
	quadratic_mle::QuadraticMleCheckProver,
};
use binius_math::{
	BinarySubspace, FieldBuffer,
	multilinear::eq::{
		eq_ind_partial_eval, eq_ind_partial_eval_scalars, eq_ind_truncate_low_inplace,
	},
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

/// Build the row equality weights for the padded `(round_trace, lane)` row space.
///
/// The active Keccak rows are a linear prefix of the padded Boolean hypercube used by the
/// remaining MLE-check. The prover's first-round message must use the same row point that the
/// verifier supplies to the MLE-check, so these weights are derived directly from those
/// challenges.
pub fn row_eq_weights<F>(row_count: usize, row_challenges: &[F]) -> Vec<F>
where
	F: Field,
{
	assert_eq!(row_challenges.len(), row_count.next_power_of_two().ilog2() as usize);
	let weights = eq_ind_partial_eval_scalars(row_challenges);
	weights[..row_count].to_vec()
}

/// Production-shaped first-round flow for verifier-field row challenges.
pub fn par_first_round_claim_from_row_challenges<FChallenge, P>(
	lookup: &NttLookup<P>,
	round_traces: &[RoundTrace],
	row_challenges: &[FChallenge],
	challenge: FChallenge,
) -> FChallenge
where
	FChallenge: BinaryField + Field + From<P::Scalar> + Send + Sync,
	P: PackedField,
	P::Scalar: BinaryField + Field,
{
	let row_count = round_traces.len() * N_LANES;
	let eq_weights = row_eq_weights(row_count, row_challenges);
	let upper_half_message = par_upper_half_round_message(lookup, round_traces, &eq_weights);
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

/// Folded chi/iota columns after the bit-axis univariate skip.
///
/// The rows are indexed by `(round_trace, lane)` and padded to the next power of two so the
/// remaining Spartan outer pass can use the standard Boolean-hypercube sumcheck machinery. This
/// mirrors production BitAnd's three-column shape by combining `R + next + iota` into the third
/// folded column.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct FoldedOuterColumns<F> {
	pub p: Vec<F>,
	pub q: Vec<F>,
	pub c: Vec<F>,
	pub log_rows: usize,
	pub row_count: usize,
}

/// Packed folded chi/iota columns after the bit-axis univariate skip.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct PackedFoldedOuterColumns<P: PackedField> {
	pub p: FieldBuffer<P>,
	pub q: FieldBuffer<P>,
	pub c: FieldBuffer<P>,
	pub log_rows: usize,
	pub row_count: usize,
}

impl<F> FoldedOuterColumns<F>
where
	F: Field,
{
	fn into_field_buffers(self) -> [FieldBuffer<F>; 3] {
		let log_rows = self.log_rows;
		[
			FieldBuffer::new(log_rows, self.p.into_boxed_slice()),
			FieldBuffer::new(log_rows, self.q.into_boxed_slice()),
			FieldBuffer::new(log_rows, self.c.into_boxed_slice()),
		]
	}

	fn into_packed_field_buffers<P>(self) -> [FieldBuffer<P>; 3]
	where
		P: PackedField<Scalar = F>,
	{
		[
			FieldBuffer::<P>::from_values(&self.p),
			FieldBuffer::<P>::from_values(&self.q),
			FieldBuffer::<P>::from_values(&self.c),
		]
	}
}

impl<P> PackedFoldedOuterColumns<P>
where
	P: PackedField,
{
	fn into_packed_field_buffers(self) -> [FieldBuffer<P>; 3] {
		[self.p, self.q, self.c]
	}
}

/// Pack folded outer columns once, outside the timed prover loop.
pub fn pack_folded_outer_columns<F, P>(
	columns: FoldedOuterColumns<F>,
) -> PackedFoldedOuterColumns<P>
where
	F: Field,
	P: PackedField<Scalar = F>,
{
	let log_rows = columns.log_rows;
	let row_count = columns.row_count;
	PackedFoldedOuterColumns {
		p: FieldBuffer::<P>::from_values(&columns.p),
		q: FieldBuffer::<P>::from_values(&columns.q),
		c: FieldBuffer::<P>::from_values(&columns.c),
		log_rows,
		row_count,
	}
}

/// Result of the Keccak chi/iota Spartan outer pass after the bit-axis univariate skip.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct SpartanOuterPass<F> {
	pub first_round_challenge: F,
	pub folded_claim: F,
	pub zerocheck_challenges: Vec<F>,
	pub round_messages: Vec<RoundCoeffs<F>>,
	pub sumcheck_challenges: Vec<F>,
	pub multilinear_evals: [F; 3],
	pub final_eval: F,
}

/// Transcripted output for the Keccak chi/iota outer segment.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct SpartanOuterTranscriptOutput<F> {
	pub p_eval: F,
	pub q_eval: F,
	pub c_eval: F,
	pub z_challenge: F,
	pub eval_point: Vec<F>,
}

/// Build the folded outer columns for the remaining chi/iota Spartan sumcheck.
pub fn folded_outer_columns<FChallenge, P>(
	round_traces: &[RoundTrace],
	challenge: FChallenge,
) -> FoldedOuterColumns<FChallenge>
where
	FChallenge: BinaryField + Field + From<P::Scalar> + WithUnderlier,
	FChallenge::Underlier: UnderlierWithBitOps,
	P: PackedField,
	P::Scalar: BinaryField + Field,
{
	assert!(!round_traces.is_empty());

	let row_count = round_traces.len() * N_LANES;
	let padded_row_count = row_count.next_power_of_two();
	let log_rows = padded_row_count.ilog2() as usize;
	let input_domain = BinarySubspace::<P::Scalar>::with_dim(crate::bit_ntt::LOG_LANE_BITS + 1)
		.reduce_dim(crate::bit_ntt::LOG_LANE_BITS)
		.isomorphic::<FChallenge>();
	let lagrange_evals = lagrange_evals_scalars(&input_domain, challenge);
	let transform = OutputWrappingTransformationFactory::new(BytewiseLookupTransformationFactory)
		.create(&lagrange_evals);

	let mut columns = FoldedOuterColumns {
		p: vec![FChallenge::ZERO; padded_row_count],
		q: vec![FChallenge::ZERO; padded_row_count],
		c: vec![FChallenge::ZERO; padded_row_count],
		log_rows,
		row_count,
	};

	columns.p[..row_count]
		.par_chunks_mut(N_LANES)
		.zip(columns.q[..row_count].par_chunks_mut(N_LANES))
		.zip(columns.c[..row_count].par_chunks_mut(N_LANES))
		.enumerate()
		.for_each(|(trace_idx, ((p_chunk, q_chunk), c_chunk))| {
			let round_trace = &round_traces[trace_idx];
			let round = trace_idx % crate::constants::N_ROUNDS;
			let folded_iota = transform.transform(&ROUND_CONSTANTS[round]);
			for lane_idx in 0..N_LANES {
				let (p_lane, q_lane, r_lane) = CHI_OPERAND_LANES[lane_idx];
				p_chunk[lane_idx] = transform.transform(&!round_trace.pre_chi[p_lane]);
				q_chunk[lane_idx] = transform.transform(&round_trace.pre_chi[q_lane]);
				let r = transform.transform(&round_trace.pre_chi[r_lane]);
				let next = transform.transform(&round_trace.output[lane_idx]);
				c_chunk[lane_idx] = r
					+ next + if lane_idx == 0 {
					folded_iota
				} else {
					FChallenge::ZERO
				};
			}
		});

	columns
}

#[cfg(test)]
fn folded_outer_columns_direct<FChallenge, P>(
	round_traces: &[RoundTrace],
	challenge: FChallenge,
) -> FoldedOuterColumns<FChallenge>
where
	FChallenge: BinaryField + Field + From<P::Scalar>,
	P: PackedField,
	P::Scalar: BinaryField + Field,
{
	assert!(!round_traces.is_empty());

	let row_count = round_traces.len() * N_LANES;
	let padded_row_count = row_count.next_power_of_two();
	let log_rows = padded_row_count.ilog2() as usize;
	let input_domain = BinarySubspace::<P::Scalar>::with_dim(crate::bit_ntt::LOG_LANE_BITS + 1)
		.reduce_dim(crate::bit_ntt::LOG_LANE_BITS)
		.isomorphic::<FChallenge>();
	let lagrange_evals = lagrange_evals_scalars(&input_domain, challenge);

	let mut columns = FoldedOuterColumns {
		p: Vec::with_capacity(padded_row_count),
		q: Vec::with_capacity(padded_row_count),
		c: Vec::with_capacity(padded_row_count),
		log_rows,
		row_count,
	};

	for (trace_idx, round_trace) in round_traces.iter().enumerate() {
		let round = trace_idx % crate::constants::N_ROUNDS;
		for lane_idx in 0..N_LANES {
			let (p_lane, q_lane, r_lane) = CHI_OPERAND_LANES[lane_idx];
			columns
				.p
				.push(fold_word(!round_trace.pre_chi[p_lane], &lagrange_evals));
			columns
				.q
				.push(fold_word(round_trace.pre_chi[q_lane], &lagrange_evals));
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
			columns.c.push(r + next + iota);
		}
	}

	columns.p.resize(padded_row_count, FChallenge::ZERO);
	columns.q.resize(padded_row_count, FChallenge::ZERO);
	columns.c.resize(padded_row_count, FChallenge::ZERO);
	columns
}

/// Evaluate the folded chi/iota relation against a Boolean-hypercube equality point.
pub fn folded_outer_claim<F>(columns: &FoldedOuterColumns<F>, zerocheck_challenges: &[F]) -> F
where
	F: Field,
{
	assert_eq!(zerocheck_challenges.len(), columns.log_rows);

	let eq_weights = eq_ind_partial_eval_scalars(zerocheck_challenges);
	(0..columns.p.len())
		.into_par_iter()
		.map(|i| (columns.p[i] * columns.q[i] - columns.c[i]) * eq_weights[i])
		.sum()
}

/// Run the remaining Spartan outer sumcheck after the bit-axis univariate skip.
///
/// This is the first full outer pass for the Keccak chi/iota relation: the caller supplies the
/// verifier's post-skip challenge, the outer zerocheck point over padded `(round_trace, lane)`
/// rows, and the per-round sumcheck challenges. The returned messages are the degree-2
/// sumcheck messages for all remaining rounds.
pub fn prove_spartan_outer_after_first_round<FChallenge, P>(
	round_traces: &[RoundTrace],
	first_round_challenge: FChallenge,
	zerocheck_challenges: Vec<FChallenge>,
	sumcheck_challenges: &[FChallenge],
) -> Result<SpartanOuterPass<FChallenge>, SumcheckError>
where
	FChallenge: BinaryField + Field + From<P::Scalar> + WithUnderlier,
	FChallenge::Underlier: UnderlierWithBitOps,
	P: PackedField,
	P::Scalar: BinaryField + Field,
{
	let columns = folded_outer_columns::<FChallenge, P>(round_traces, first_round_challenge);
	assert_eq!(zerocheck_challenges.len(), columns.log_rows);
	assert_eq!(sumcheck_challenges.len(), columns.log_rows);

	let folded_claim = folded_outer_claim(&columns, &zerocheck_challenges);
	prove_spartan_outer_from_folded_columns(
		columns,
		first_round_challenge,
		zerocheck_challenges,
		sumcheck_challenges,
		folded_claim,
	)
}

/// Prove the Keccak chi/iota outer segment with a transcript channel.
pub fn prove_spartan_outer_with_channel<FChallenge, P, Channel>(
	lookup: &NttLookup<P>,
	round_traces: &[RoundTrace],
	zerocheck_challenges: Vec<FChallenge>,
	channel: &mut Channel,
) -> Result<SpartanOuterTranscriptOutput<FChallenge>, SumcheckError>
where
	FChallenge: BinaryField + Field + From<P::Scalar> + Send + Sync + WithUnderlier,
	FChallenge::Underlier: UnderlierWithBitOps,
	P: PackedField,
	P::Scalar: BinaryField + Field,
	Channel: IPProverChannel<FChallenge>,
{
	let row_count = round_traces.len() * N_LANES;
	let eq_weights = row_eq_weights(row_count, &zerocheck_challenges);
	let upper_half_message =
		par_upper_half_round_message::<FChallenge, P>(lookup, round_traces, &eq_weights);
	channel.send_many(&upper_half_message);

	let z_challenge = channel.sample();
	let prover_message_domain =
		BinarySubspace::<P::Scalar>::with_dim(crate::bit_ntt::LOG_LANE_BITS + 1)
			.isomorphic::<FChallenge>();
	let folded_claim =
		next_sum_claim_from_upper_half(&upper_half_message, z_challenge, &prover_message_domain);
	let columns = folded_outer_columns::<FChallenge, P>(round_traces, z_challenge);
	assert_eq!(zerocheck_challenges.len(), columns.log_rows);

	let prover = QuadraticMleCheckProver::new(
		columns.into_field_buffers(),
		|[p, q, c]| p * q - c,
		|[p, q, _]| p * q,
		zerocheck_challenges,
		folded_claim,
	)?;
	let prove_output = prove_single_mlecheck(prover, channel)?;

	assert_eq!(prove_output.multilinear_evals.len(), 3);
	channel.send_many(&prove_output.multilinear_evals);

	let mut eval_point = prove_output.challenges;
	eval_point.reverse();

	Ok(SpartanOuterTranscriptOutput {
		p_eval: prove_output.multilinear_evals[0],
		q_eval: prove_output.multilinear_evals[1],
		c_eval: prove_output.multilinear_evals[2],
		z_challenge,
		eval_point,
	})
}

/// Verify the transcripted Keccak chi/iota outer segment.
pub fn verify_spartan_outer_with_channel<F, Channel>(
	zerocheck_challenges: &[Channel::Elem],
	channel: &mut Channel,
	round_message_univariate_domain: &BinarySubspace<Channel::Elem>,
) -> Result<SpartanOuterTranscriptOutput<Channel::Elem>, binius_ip::sumcheck::Error>
where
	F: BinaryField,
	Channel: IPVerifierChannel<F>,
	Channel::Elem: BinaryField,
{
	let upper_half_message = channel.recv_array::<LANE_BITS>()?;
	let z_challenge = channel.sample();
	let folded_claim = next_sum_claim_from_upper_half(
		&upper_half_message,
		z_challenge,
		round_message_univariate_domain,
	);

	let binius_ip::sumcheck::SumcheckOutput {
		eval,
		challenges: mut eval_point,
	} = mlecheck::verify(zerocheck_challenges, 2, folded_claim, channel)?;

	let [p_eval, q_eval, c_eval] = channel.recv_array()?;
	channel.assert_zero(p_eval * q_eval - c_eval - eval)?;
	eval_point.reverse();

	Ok(SpartanOuterTranscriptOutput {
		p_eval,
		q_eval,
		c_eval,
		z_challenge,
		eval_point,
	})
}

/// Run the post-skip Spartan outer pass when the folded claim is already known.
///
/// In the full BitAnd-style flow this claim is obtained by extrapolating the univariate skip
/// message at the verifier challenge, so the remaining MLE-check should not recompute it by
/// materializing a full equality tensor.
pub fn prove_spartan_outer_after_first_round_with_claim<FChallenge, P>(
	round_traces: &[RoundTrace],
	first_round_challenge: FChallenge,
	zerocheck_challenges: Vec<FChallenge>,
	sumcheck_challenges: &[FChallenge],
	folded_claim: FChallenge,
) -> Result<SpartanOuterPass<FChallenge>, SumcheckError>
where
	FChallenge: BinaryField + Field + From<P::Scalar> + WithUnderlier,
	FChallenge::Underlier: UnderlierWithBitOps,
	P: PackedField,
	P::Scalar: BinaryField + Field,
{
	let columns = folded_outer_columns::<FChallenge, P>(round_traces, first_round_challenge);
	assert_eq!(zerocheck_challenges.len(), columns.log_rows);
	assert_eq!(sumcheck_challenges.len(), columns.log_rows);

	prove_spartan_outer_from_folded_columns(
		columns,
		first_round_challenge,
		zerocheck_challenges,
		sumcheck_challenges,
		folded_claim,
	)
}

/// Experimental fused/adaptive post-skip outer pass local to Keccak v0.
///
/// This keeps production Binius MLE-check code untouched while testing the scheduler shape we want:
/// persistent workers per phase, fused bind-then-reduce rounds, adaptive worker shrinkage, and a
/// sequential tail when the live domain becomes too small.
pub fn prove_spartan_outer_after_first_round_with_claim_fused_adaptive<FChallenge, P>(
	round_traces: &[RoundTrace],
	first_round_challenge: FChallenge,
	zerocheck_challenges: Vec<FChallenge>,
	sumcheck_challenges: &[FChallenge],
	folded_claim: FChallenge,
) -> SpartanOuterPass<FChallenge>
where
	FChallenge: BinaryField + Field + From<P::Scalar> + WithUnderlier,
	FChallenge::Underlier: UnderlierWithBitOps,
	P: PackedField,
	P::Scalar: BinaryField + Field,
{
	let columns = folded_outer_columns::<FChallenge, P>(round_traces, first_round_challenge);
	assert_eq!(zerocheck_challenges.len(), columns.log_rows);
	assert_eq!(sumcheck_challenges.len(), columns.log_rows);

	prove_spartan_outer_from_folded_columns_fused_adaptive(
		columns,
		first_round_challenge,
		zerocheck_challenges,
		sumcheck_challenges,
		folded_claim,
	)
}

/// Run the generic post-skip outer pass from already folded columns.
///
/// This is primarily useful for benchmarking the remaining Spartan outer rounds separately from
/// folded-column construction.
pub fn prove_spartan_outer_from_folded_columns_with_claim<FChallenge>(
	columns: FoldedOuterColumns<FChallenge>,
	first_round_challenge: FChallenge,
	zerocheck_challenges: Vec<FChallenge>,
	sumcheck_challenges: &[FChallenge],
	folded_claim: FChallenge,
) -> Result<SpartanOuterPass<FChallenge>, SumcheckError>
where
	FChallenge: BinaryField + Field,
{
	prove_spartan_outer_from_folded_columns(
		columns,
		first_round_challenge,
		zerocheck_challenges,
		sumcheck_challenges,
		folded_claim,
	)
}

/// Run the experimental fused/adaptive post-skip outer pass from already folded columns.
pub fn prove_spartan_outer_from_folded_columns_with_claim_fused_adaptive<F>(
	columns: FoldedOuterColumns<F>,
	first_round_challenge: F,
	zerocheck_challenges: Vec<F>,
	sumcheck_challenges: &[F],
	folded_claim: F,
) -> SpartanOuterPass<F>
where
	F: BinaryField + Field,
{
	prove_spartan_outer_from_folded_columns_fused_adaptive(
		columns,
		first_round_challenge,
		zerocheck_challenges,
		sumcheck_challenges,
		folded_claim,
	)
}

/// Run the generic post-skip outer pass from folded columns using an explicitly packed field.
pub fn prove_spartan_outer_from_folded_columns_with_claim_packed<F, P>(
	columns: FoldedOuterColumns<F>,
	first_round_challenge: F,
	zerocheck_challenges: Vec<F>,
	sumcheck_challenges: &[F],
	folded_claim: F,
) -> Result<SpartanOuterPass<F>, SumcheckError>
where
	F: BinaryField + Field,
	P: PackedField<Scalar = F>,
{
	prove_spartan_outer_from_packed_buffers(
		columns.into_packed_field_buffers::<P>(),
		first_round_challenge,
		zerocheck_challenges,
		sumcheck_challenges,
		folded_claim,
	)
}

/// Run the fused bind/reduce post-skip outer pass from folded columns using packed buffers.
pub fn prove_spartan_outer_from_folded_columns_with_claim_packed_fused<F, P>(
	columns: FoldedOuterColumns<F>,
	first_round_challenge: F,
	zerocheck_challenges: Vec<F>,
	sumcheck_challenges: &[F],
	folded_claim: F,
) -> Result<SpartanOuterPass<F>, SumcheckError>
where
	F: BinaryField + Field,
	P: PackedField<Scalar = F>,
{
	prove_spartan_outer_from_packed_buffers_fused::<F, P>(
		columns.into_packed_field_buffers::<P>(),
		first_round_challenge,
		zerocheck_challenges,
		sumcheck_challenges,
		folded_claim,
	)
}

/// Run the generic post-skip outer pass from pre-packed folded columns.
pub fn prove_spartan_outer_from_packed_folded_columns_with_claim<F, P>(
	columns: PackedFoldedOuterColumns<P>,
	first_round_challenge: F,
	zerocheck_challenges: Vec<F>,
	sumcheck_challenges: &[F],
	folded_claim: F,
) -> Result<SpartanOuterPass<F>, SumcheckError>
where
	F: BinaryField + Field,
	P: PackedField<Scalar = F>,
{
	prove_spartan_outer_from_packed_buffers(
		columns.into_packed_field_buffers(),
		first_round_challenge,
		zerocheck_challenges,
		sumcheck_challenges,
		folded_claim,
	)
}

/// Run the fused bind/reduce post-skip outer pass from pre-packed folded columns.
pub fn prove_spartan_outer_from_packed_folded_columns_with_claim_fused<F, P>(
	columns: PackedFoldedOuterColumns<P>,
	first_round_challenge: F,
	zerocheck_challenges: Vec<F>,
	sumcheck_challenges: &[F],
	folded_claim: F,
) -> Result<SpartanOuterPass<F>, SumcheckError>
where
	F: BinaryField + Field,
	P: PackedField<Scalar = F>,
{
	prove_spartan_outer_from_packed_buffers_fused::<F, P>(
		columns.into_packed_field_buffers(),
		first_round_challenge,
		zerocheck_challenges,
		sumcheck_challenges,
		folded_claim,
	)
}

/// Run the persistent-worker fused post-skip outer pass from pre-packed folded columns.
pub fn prove_spartan_outer_from_packed_folded_columns_with_claim_persistent_fused<F, P>(
	columns: PackedFoldedOuterColumns<P>,
	first_round_challenge: F,
	zerocheck_challenges: Vec<F>,
	sumcheck_challenges: &[F],
	folded_claim: F,
) -> Result<SpartanOuterPass<F>, SumcheckError>
where
	F: BinaryField + Field,
	P: PackedField<Scalar = F> + Send + Sync,
{
	prove_spartan_outer_from_packed_buffers_persistent_fused::<F, P>(
		columns.into_packed_field_buffers(),
		first_round_challenge,
		zerocheck_challenges,
		sumcheck_challenges,
		folded_claim,
	)
}

fn prove_spartan_outer_from_folded_columns<FChallenge>(
	columns: FoldedOuterColumns<FChallenge>,
	first_round_challenge: FChallenge,
	zerocheck_challenges: Vec<FChallenge>,
	sumcheck_challenges: &[FChallenge],
	folded_claim: FChallenge,
) -> Result<SpartanOuterPass<FChallenge>, SumcheckError>
where
	FChallenge: BinaryField + Field,
{
	prove_spartan_outer_from_packed_buffers(
		columns.into_field_buffers(),
		first_round_challenge,
		zerocheck_challenges,
		sumcheck_challenges,
		folded_claim,
	)
}

fn prove_spartan_outer_from_packed_buffers<F, P>(
	packed_buffers: [FieldBuffer<P>; 3],
	first_round_challenge: F,
	zerocheck_challenges: Vec<F>,
	sumcheck_challenges: &[F],
	folded_claim: F,
) -> Result<SpartanOuterPass<F>, SumcheckError>
where
	F: BinaryField + Field,
	P: PackedField<Scalar = F>,
{
	let log_rows = packed_buffers[0].log_len();
	assert_eq!(packed_buffers[1].log_len(), log_rows);
	assert_eq!(packed_buffers[2].log_len(), log_rows);
	assert_eq!(zerocheck_challenges.len(), log_rows);
	assert_eq!(sumcheck_challenges.len(), log_rows);

	let mut prover = QuadraticMleCheckProver::new(
		packed_buffers,
		|[p, q, c]| p * q - c,
		|[p, q, _]| p * q,
		zerocheck_challenges.clone(),
		folded_claim,
	)?;

	let mut round_messages = Vec::with_capacity(log_rows);
	for &challenge in sumcheck_challenges {
		let mut coeffs = prover.execute()?;
		debug_assert_eq!(coeffs.len(), 1);
		round_messages.push(coeffs.remove(0));
		prover.fold(challenge)?;
	}

	let multilinear_evals: [F; 3] = prover.finish()?.try_into().expect("three folded columns");
	let final_eval = multilinear_evals[0] * multilinear_evals[1] - multilinear_evals[2];

	Ok(SpartanOuterPass {
		first_round_challenge,
		folded_claim,
		zerocheck_challenges,
		round_messages,
		sumcheck_challenges: sumcheck_challenges.to_vec(),
		multilinear_evals,
		final_eval,
	})
}

fn prove_spartan_outer_from_packed_buffers_fused<F, P>(
	mut packed_buffers: [FieldBuffer<P>; 3],
	first_round_challenge: F,
	zerocheck_challenges: Vec<F>,
	sumcheck_challenges: &[F],
	folded_claim: F,
) -> Result<SpartanOuterPass<F>, SumcheckError>
where
	F: BinaryField + Field,
	P: PackedField<Scalar = F>,
{
	let log_rows = packed_buffers[0].log_len();
	assert_eq!(packed_buffers[1].log_len(), log_rows);
	assert_eq!(packed_buffers[2].log_len(), log_rows);
	assert_eq!(zerocheck_challenges.len(), log_rows);
	assert_eq!(sumcheck_challenges.len(), log_rows);

	let full_zerocheck_challenges = zerocheck_challenges.clone();
	let mut eq = eq_ind_partial_eval::<P>(&zerocheck_challenges[..log_rows.saturating_sub(1)]);
	let mut round_messages = Vec::with_capacity(log_rows);
	let mut last_eval = folded_claim;
	let mut rounds_done = 0;
	let min_packed_reduce_words = keccak_packed_fused_min_reduce_words();

	while rounds_done < log_rows {
		let n_vars = log_rows - rounds_done;
		if n_vars <= P::LOG_WIDTH {
			break;
		}
		let reduce_words = 1 << (n_vars - 1 - P::LOG_WIDTH);
		if rounds_done > 0 && reduce_words < min_packed_reduce_words {
			break;
		}

		let alpha = zerocheck_challenges[n_vars - 1];
		let partial = if rounds_done == 0 {
			packed_reduce_round::<F, P>(&packed_buffers, &eq, n_vars)
		} else {
			packed_bind_then_reduce_round::<F, P>(
				&mut packed_buffers,
				&eq,
				n_vars,
				sumcheck_challenges[rounds_done - 1],
			)
		};
		if rounds_done > 0 {
			for buffer in &mut packed_buffers {
				buffer.truncate(n_vars);
			}
		}
		let coeffs = partial.interpolate_eq(last_eval, alpha);
		last_eval = coeffs.evaluate(sumcheck_challenges[rounds_done]);
		round_messages.push(coeffs);
		rounds_done += 1;

		if n_vars > 1 {
			eq_ind_truncate_low_inplace(&mut eq, n_vars - 2);
		}
	}

	let remaining_vars = log_rows - rounds_done;
	if rounds_done > 0 {
		packed_bind_round(
			&mut packed_buffers,
			remaining_vars + 1,
			sumcheck_challenges[rounds_done - 1],
		);
		for buffer in &mut packed_buffers {
			buffer.truncate(remaining_vars);
		}
	}

	if remaining_vars > 0 {
		let tail = prove_spartan_outer_from_packed_buffers(
			packed_buffers,
			first_round_challenge,
			zerocheck_challenges[..remaining_vars].to_vec(),
			&sumcheck_challenges[rounds_done..],
			last_eval,
		)?;
		round_messages.extend(tail.round_messages);
		return Ok(SpartanOuterPass {
			first_round_challenge,
			folded_claim,
			zerocheck_challenges: full_zerocheck_challenges,
			round_messages,
			sumcheck_challenges: sumcheck_challenges.to_vec(),
			multilinear_evals: tail.multilinear_evals,
			final_eval: tail.final_eval,
		});
	}

	let multilinear_evals = [
		packed_buffers[0].get(0),
		packed_buffers[1].get(0),
		packed_buffers[2].get(0),
	];
	let final_eval = multilinear_evals[0] * multilinear_evals[1] - multilinear_evals[2];

	Ok(SpartanOuterPass {
		first_round_challenge,
		folded_claim,
		zerocheck_challenges: full_zerocheck_challenges,
		round_messages,
		sumcheck_challenges: sumcheck_challenges.to_vec(),
		multilinear_evals,
		final_eval,
	})
}

fn prove_spartan_outer_from_packed_buffers_persistent_fused<F, P>(
	packed_buffers: [FieldBuffer<P>; 3],
	first_round_challenge: F,
	zerocheck_challenges: Vec<F>,
	sumcheck_challenges: &[F],
	folded_claim: F,
) -> Result<SpartanOuterPass<F>, SumcheckError>
where
	F: BinaryField + Field,
	P: PackedField<Scalar = F> + Send + Sync,
{
	let log_rows = packed_buffers[0].log_len();
	assert_eq!(packed_buffers[1].log_len(), log_rows);
	assert_eq!(packed_buffers[2].log_len(), log_rows);
	assert_eq!(zerocheck_challenges.len(), log_rows);
	assert_eq!(sumcheck_challenges.len(), log_rows);

	let full_zerocheck_challenges = zerocheck_challenges.clone();
	let mut state = PackedFusedOuterState::new(packed_buffers, &zerocheck_challenges);
	let mut round_messages = Vec::with_capacity(log_rows);
	let mut last_eval = folded_claim;
	let rounds_done = state.run_persistent_parallel_prefix(
		&zerocheck_challenges,
		sumcheck_challenges,
		&mut last_eval,
		&mut round_messages,
	);

	let remaining_vars = log_rows - rounds_done;
	let packed_buffers = state.into_truncated_packed_buffers(remaining_vars);
	if remaining_vars > 0 {
		let tail = prove_spartan_outer_from_packed_buffers(
			packed_buffers,
			first_round_challenge,
			zerocheck_challenges[..remaining_vars].to_vec(),
			&sumcheck_challenges[rounds_done..],
			last_eval,
		)?;
		round_messages.extend(tail.round_messages);
		return Ok(SpartanOuterPass {
			first_round_challenge,
			folded_claim,
			zerocheck_challenges: full_zerocheck_challenges,
			round_messages,
			sumcheck_challenges: sumcheck_challenges.to_vec(),
			multilinear_evals: tail.multilinear_evals,
			final_eval: tail.final_eval,
		});
	}

	let multilinear_evals = [
		packed_buffers[0].get(0),
		packed_buffers[1].get(0),
		packed_buffers[2].get(0),
	];
	let final_eval = multilinear_evals[0] * multilinear_evals[1] - multilinear_evals[2];

	Ok(SpartanOuterPass {
		first_round_challenge,
		folded_claim,
		zerocheck_challenges: full_zerocheck_challenges,
		round_messages,
		sumcheck_challenges: sumcheck_challenges.to_vec(),
		multilinear_evals,
		final_eval,
	})
}

fn packed_reduce_round<F, P>(
	buffers: &[FieldBuffer<P>; 3],
	eq: &FieldBuffer<P>,
	n_vars: usize,
) -> OuterPartial<F>
where
	F: Field,
	P: PackedField<Scalar = F>,
{
	debug_assert!(n_vars > P::LOG_WIDTH);
	debug_assert_eq!(buffers[0].log_len(), n_vars);
	debug_assert_eq!(eq.log_len(), n_vars - 1);

	let (p_0, p_1) = buffers[0].split_half_ref();
	let (q_0, q_1) = buffers[1].split_half_ref();
	let (_, c_1) = buffers[2].split_half_ref();
	let (y_1, y_inf) =
		(p_0.as_ref(), p_1.as_ref(), q_0.as_ref(), q_1.as_ref(), c_1.as_ref(), eq.as_ref())
			.into_par_iter()
			.map(|(p_0, p_1, q_0, q_1, c_1, eq_i)| {
				((*p_1 * *q_1 - *c_1) * *eq_i, ((*p_0 + *p_1) * (*q_0 + *q_1)) * *eq_i)
			})
			.reduce(
				|| (P::zero(), P::zero()),
				|(lhs_1, lhs_inf), (rhs_1, rhs_inf)| (lhs_1 + rhs_1, lhs_inf + rhs_inf),
			);

	OuterPartial {
		y_1: sum_packed(y_1),
		y_inf: sum_packed(y_inf),
	}
}

fn packed_bind_then_reduce_round<F, P>(
	buffers: &mut [FieldBuffer<P>; 3],
	eq: &FieldBuffer<P>,
	n_vars: usize,
	prev_challenge: F,
) -> OuterPartial<F>
where
	F: Field,
	P: PackedField<Scalar = F>,
{
	debug_assert!(n_vars > P::LOG_WIDTH);
	debug_assert_eq!(buffers[0].log_len(), n_vars + 1);
	debug_assert_eq!(eq.log_len(), n_vars - 1);

	let bind_offset = 1 << (n_vars - P::LOG_WIDTH);
	let reduce_offset = 1 << (n_vars - 1 - P::LOG_WIDTH);
	let challenge = P::broadcast(prev_challenge);
	let [p, q, c] = buffers;
	let (p_lo, p_hi, p_bind_lo, p_bind_hi) = split_bind_reduce_slices(p.as_mut(), bind_offset);
	let (q_lo, q_hi, q_bind_lo, q_bind_hi) = split_bind_reduce_slices(q.as_mut(), bind_offset);
	let (c_lo, c_hi, c_bind_lo, c_bind_hi) = split_bind_reduce_slices(c.as_mut(), bind_offset);
	debug_assert_eq!(p_lo.len(), reduce_offset);

	let columns = PackedBindReduceColumns {
		p_lo: p_lo.as_mut_ptr(),
		p_hi: p_hi.as_mut_ptr(),
		p_bind_lo: p_bind_lo.as_ptr(),
		p_bind_hi: p_bind_hi.as_ptr(),
		q_lo: q_lo.as_mut_ptr(),
		q_hi: q_hi.as_mut_ptr(),
		q_bind_lo: q_bind_lo.as_ptr(),
		q_bind_hi: q_bind_hi.as_ptr(),
		c_lo: c_lo.as_mut_ptr(),
		c_hi: c_hi.as_mut_ptr(),
		c_bind_lo: c_bind_lo.as_ptr(),
		c_bind_hi: c_bind_hi.as_ptr(),
		eq: eq.as_ref().as_ptr(),
	};
	let (y_1, y_inf) = (0..reduce_offset)
		.into_par_iter()
		.map(|i| columns.bind_then_reduce(i, challenge))
		.reduce(
			|| (P::zero(), P::zero()),
			|(lhs_1, lhs_inf), (rhs_1, rhs_inf)| (lhs_1 + rhs_1, lhs_inf + rhs_inf),
		);

	OuterPartial {
		y_1: sum_packed(y_1),
		y_inf: sum_packed(y_inf),
	}
}

struct PackedBindReduceColumns<P> {
	p_lo: *mut P,
	p_hi: *mut P,
	p_bind_lo: *const P,
	p_bind_hi: *const P,
	q_lo: *mut P,
	q_hi: *mut P,
	q_bind_lo: *const P,
	q_bind_hi: *const P,
	c_lo: *mut P,
	c_hi: *mut P,
	c_bind_lo: *const P,
	c_bind_hi: *const P,
	eq: *const P,
}

// SAFETY: callers partition the active low/high mutable prefixes into disjoint packed slots and
// only invoke `bind_then_reduce` with unique indices in a parallel iterator.
unsafe impl<P: Send + Sync> Sync for PackedBindReduceColumns<P> {}

impl<P> PackedBindReduceColumns<P>
where
	P: PackedField,
{
	fn bind_then_reduce(&self, i: usize, challenge: P) -> (P, P) {
		// SAFETY: `i` is provided by `0..reduce_offset`, so every pointer addition lands in bounds.
		// Low and high slots are disjoint slices, and each parallel task writes a unique index.
		unsafe {
			let p_0_slot = &mut *self.p_lo.add(i);
			let p_1_slot = &mut *self.p_hi.add(i);
			let q_0_slot = &mut *self.q_lo.add(i);
			let q_1_slot = &mut *self.q_hi.add(i);
			let c_0_slot = &mut *self.c_lo.add(i);
			let c_1_slot = &mut *self.c_hi.add(i);

			let p_0 = *p_0_slot + challenge * (*self.p_bind_lo.add(i) - *p_0_slot);
			let p_1 = *p_1_slot + challenge * (*self.p_bind_hi.add(i) - *p_1_slot);
			let q_0 = *q_0_slot + challenge * (*self.q_bind_lo.add(i) - *q_0_slot);
			let q_1 = *q_1_slot + challenge * (*self.q_bind_hi.add(i) - *q_1_slot);
			let c_0 = *c_0_slot + challenge * (*self.c_bind_lo.add(i) - *c_0_slot);
			let c_1 = *c_1_slot + challenge * (*self.c_bind_hi.add(i) - *c_1_slot);

			*p_0_slot = p_0;
			*p_1_slot = p_1;
			*q_0_slot = q_0;
			*q_1_slot = q_1;
			*c_0_slot = c_0;
			*c_1_slot = c_1;

			let eq_i = *self.eq.add(i);
			((p_1 * q_1 - c_1) * eq_i, ((p_0 + p_1) * (q_0 + q_1)) * eq_i)
		}
	}
}

#[derive(Clone, Copy)]
struct PackedOuterPartial<P> {
	y_1: P,
	y_inf: P,
}

impl<P> PackedOuterPartial<P>
where
	P: PackedField,
{
	fn zero() -> Self {
		Self {
			y_1: P::zero(),
			y_inf: P::zero(),
		}
	}

	fn add_assign(&mut self, rhs: Self) {
		self.y_1 += rhs.y_1;
		self.y_inf += rhs.y_inf;
	}

	fn to_outer_partial(self) -> OuterPartial<P::Scalar> {
		OuterPartial {
			y_1: sum_packed(self.y_1),
			y_inf: sum_packed(self.y_inf),
		}
	}
}

struct PackedPartialSlots<P>(Vec<UnsafeCell<PackedOuterPartial<P>>>);

// SAFETY: each worker writes only its own slot, and worker 0 reads slots only after the phase
// barrier proves every active writer has completed.
unsafe impl<P: Send> Sync for PackedPartialSlots<P> {}

struct SharedPackedOuterColumns<P: PackedField> {
	p: UnsafeCell<FieldBuffer<P>>,
	q: UnsafeCell<FieldBuffer<P>>,
	c: UnsafeCell<FieldBuffer<P>>,
	eq: UnsafeCell<FieldBuffer<P>>,
}

// SAFETY: all parallel methods partition packed-word ranges into disjoint windows. Each write lands
// in the low live prefix owned by that worker; high halves are read-only for the current round.
unsafe impl<P: PackedField + Send> Sync for SharedPackedOuterColumns<P> {}

struct PackedFusedOuterState<P: PackedField> {
	columns: SharedPackedOuterColumns<P>,
	log_rows: usize,
}

impl<P> PackedFusedOuterState<P>
where
	P: PackedField + Send + Sync,
	P::Scalar: BinaryField + Field,
{
	fn new(packed_buffers: [FieldBuffer<P>; 3], zerocheck_challenges: &[P::Scalar]) -> Self {
		let [p, q, c] = packed_buffers;
		let log_rows = p.log_len();
		let eq = eq_ind_partial_eval::<P>(&zerocheck_challenges[..log_rows.saturating_sub(1)]);
		Self {
			log_rows,
			columns: SharedPackedOuterColumns {
				p: UnsafeCell::new(p),
				q: UnsafeCell::new(q),
				c: UnsafeCell::new(c),
				eq: UnsafeCell::new(eq),
			},
		}
	}

	fn run_persistent_parallel_prefix(
		&mut self,
		zerocheck_challenges: &[P::Scalar],
		sumcheck_challenges: &[P::Scalar],
		last_eval: &mut P::Scalar,
		round_messages: &mut Vec<RoundCoeffs<P::Scalar>>,
	) -> usize {
		let plans = self.parallel_round_plans();
		let max_workers = plans.iter().map(|plan| plan.workers).max().unwrap_or(1);
		if max_workers < 2 {
			return 0;
		}

		let partials = PackedPartialSlots(
			(0..max_workers)
				.map(|_| UnsafeCell::new(PackedOuterPartial::zero()))
				.collect(),
		);
		let barrier = Barrier::new(max_workers);

		thread::scope(|scope| {
			for worker_idx in 1..max_workers {
				let partials = &partials;
				let barrier = &barrier;
				let state = &*self;
				let plans = &plans;
				scope.spawn(move || {
					for plan in plans {
						if worker_idx < plan.workers {
							let (lo, hi) =
								static_chunk_range(plan.reduce_words, plan.workers, worker_idx);
							let partial = if plan.abs_round == 0 {
								state.reduce_round(plan.reduce_words, lo, hi)
							} else {
								state.bind_then_reduce_round(
									plan.n_vars,
									plan.reduce_words,
									lo,
									hi,
									sumcheck_challenges[plan.abs_round - 1],
								)
							};
							// SAFETY: this worker is the unique writer for its slot.
							unsafe {
								*partials.0[worker_idx].get() = partial;
							}
						}
						barrier.wait();
						if worker_idx < plan.workers && plan.n_vars > 1 {
							let next_eq_words = plan.reduce_words >> 1;
							let (lo, hi) =
								static_chunk_range(next_eq_words, plan.workers, worker_idx);
							state.truncate_eq_round(next_eq_words, lo, hi);
						}
						barrier.wait();
					}

					if let Some(last_plan) = plans.last()
						&& worker_idx < last_plan.workers
					{
						let (lo, hi) = static_chunk_range(
							last_plan.reduce_words,
							last_plan.workers,
							worker_idx,
						);
						state.bind_round(
							last_plan.reduce_words,
							lo,
							hi,
							sumcheck_challenges[last_plan.abs_round],
						);
					}
				});
			}

			for plan in &plans {
				let (lo, hi) = static_chunk_range(plan.reduce_words, plan.workers, 0);
				let partial = if plan.abs_round == 0 {
					self.reduce_round(plan.reduce_words, lo, hi)
				} else {
					self.bind_then_reduce_round(
						plan.n_vars,
						plan.reduce_words,
						lo,
						hi,
						sumcheck_challenges[plan.abs_round - 1],
					)
				};
				// SAFETY: worker 0 is the unique writer for slot 0.
				unsafe {
					*partials.0[0].get() = partial;
				}
				barrier.wait();

				let mut sum = PackedOuterPartial::zero();
				for slot in &partials.0[..plan.workers] {
					// SAFETY: all workers reached the barrier after writing their slots.
					sum.add_assign(unsafe { *slot.get() });
				}
				let alpha = zerocheck_challenges[plan.n_vars - 1];
				let coeffs = sum.to_outer_partial().interpolate_eq(*last_eval, alpha);
				*last_eval = coeffs.evaluate(sumcheck_challenges[plan.abs_round]);
				round_messages.push(coeffs);

				if plan.n_vars > 1 {
					let next_eq_words = plan.reduce_words >> 1;
					let (lo, hi) = static_chunk_range(next_eq_words, plan.workers, 0);
					self.truncate_eq_round(next_eq_words, lo, hi);
				}
				barrier.wait();
			}

			if let Some(last_plan) = plans.last() {
				let (lo, hi) = static_chunk_range(last_plan.reduce_words, last_plan.workers, 0);
				self.bind_round(
					last_plan.reduce_words,
					lo,
					hi,
					sumcheck_challenges[last_plan.abs_round],
				);
			}
		});

		plans.len()
	}

	fn parallel_round_plans(&self) -> Vec<PackedParallelRoundPlan> {
		let mut plans = Vec::new();
		let mut workers = keccak_outer_max_workers().max(1);
		let min_reduce_words = keccak_packed_fused_min_reduce_words();
		let min_words_per_worker = keccak_packed_outer_min_words_per_worker();

		for abs_round in 0..self.log_rows {
			let n_vars = self.log_rows - abs_round;
			if n_vars <= P::LOG_WIDTH {
				break;
			}
			let reduce_words = 1 << (n_vars - 1 - P::LOG_WIDTH);
			if abs_round > 0 && reduce_words < min_reduce_words {
				break;
			}
			workers = workers.min(reduce_words);
			while workers >= 2 && reduce_words / workers < min_words_per_worker {
				workers /= 2;
			}
			if workers < 2 {
				break;
			}
			plans.push(PackedParallelRoundPlan {
				abs_round,
				n_vars,
				reduce_words,
				workers,
			});
		}

		plans
	}

	fn reduce_round(&self, reduce_words: usize, lo: usize, hi: usize) -> PackedOuterPartial<P> {
		// SAFETY: reduce is read-only over the current live prefix.
		let p = unsafe { (&*self.columns.p.get()).as_ref() };
		let q = unsafe { (&*self.columns.q.get()).as_ref() };
		let c = unsafe { (&*self.columns.c.get()).as_ref() };
		let eq = unsafe { (&*self.columns.eq.get()).as_ref() };
		let mut partial = PackedOuterPartial::zero();
		for i in lo..hi {
			let p_1 = p[i + reduce_words];
			let q_1 = q[i + reduce_words];
			let c_1 = c[i + reduce_words];
			let weight = eq[i];
			partial.y_1 += (p_1 * q_1 - c_1) * weight;
			partial.y_inf += ((p[i] + p_1) * (q[i] + q_1)) * weight;
		}
		partial
	}

	fn bind_then_reduce_round(
		&self,
		n_vars: usize,
		reduce_words: usize,
		lo: usize,
		hi: usize,
		prev_challenge: P::Scalar,
	) -> PackedOuterPartial<P> {
		debug_assert!(n_vars > P::LOG_WIDTH);
		let bind_offset = 1 << (n_vars - P::LOG_WIDTH);
		let challenge = P::broadcast(prev_challenge);
		let p = unsafe { (&mut *self.columns.p.get()).as_mut() };
		let q = unsafe { (&mut *self.columns.q.get()).as_mut() };
		let c = unsafe { (&mut *self.columns.c.get()).as_mut() };
		let eq = unsafe { (&*self.columns.eq.get()).as_ref() };
		let mut partial = PackedOuterPartial::zero();

		for i in lo..hi {
			let i_hi = i + bind_offset;
			let p_0 = p[i] + challenge * (p[i_hi] - p[i]);
			let q_0 = q[i] + challenge * (q[i_hi] - q[i]);
			let c_0 = c[i] + challenge * (c[i_hi] - c[i]);
			p[i] = p_0;
			q[i] = q_0;
			c[i] = c_0;

			let j = i + reduce_words;
			let j_hi = j + bind_offset;
			let p_1 = p[j] + challenge * (p[j_hi] - p[j]);
			let q_1 = q[j] + challenge * (q[j_hi] - q[j]);
			let c_1 = c[j] + challenge * (c[j_hi] - c[j]);
			p[j] = p_1;
			q[j] = q_1;
			c[j] = c_1;

			let weight = eq[i];
			partial.y_1 += (p_1 * q_1 - c_1) * weight;
			partial.y_inf += ((p_0 + p_1) * (q_0 + q_1)) * weight;
		}

		partial
	}

	fn bind_round(&self, reduce_words: usize, lo: usize, hi: usize, challenge: P::Scalar) {
		let challenge = P::broadcast(challenge);
		let p = unsafe { (&mut *self.columns.p.get()).as_mut() };
		let q = unsafe { (&mut *self.columns.q.get()).as_mut() };
		let c = unsafe { (&mut *self.columns.c.get()).as_mut() };
		for i in lo..hi {
			let p_lo = p[i];
			let q_lo = q[i];
			let c_lo = c[i];
			p[i] = p_lo + challenge * (p[i + reduce_words] - p_lo);
			q[i] = q_lo + challenge * (q[i + reduce_words] - q_lo);
			c[i] = c_lo + challenge * (c[i + reduce_words] - c_lo);
		}
	}

	fn truncate_eq_round(&self, next_eq_words: usize, lo: usize, hi: usize) {
		let eq = unsafe { (&mut *self.columns.eq.get()).as_mut() };
		for i in lo..hi {
			let high = eq[i + next_eq_words];
			eq[i] += high;
		}
	}

	fn into_truncated_packed_buffers(self, log_len: usize) -> [FieldBuffer<P>; 3] {
		let mut buffers = [
			self.columns.p.into_inner(),
			self.columns.q.into_inner(),
			self.columns.c.into_inner(),
		];
		for buffer in &mut buffers {
			buffer.truncate(log_len);
		}
		buffers
	}
}

#[derive(Clone, Copy)]
struct PackedParallelRoundPlan {
	abs_round: usize,
	n_vars: usize,
	reduce_words: usize,
	workers: usize,
}

fn keccak_packed_outer_min_words_per_worker() -> usize {
	std::env::var("KECCAK_PACKED_OUTER_MIN_WORDS_PER_WORKER")
		.ok()
		.and_then(|value| value.parse().ok())
		.unwrap_or(1 << 14)
}

fn split_bind_reduce_slices<P>(
	values: &mut [P],
	bind_offset: usize,
) -> (&mut [P], &mut [P], &[P], &[P])
where
	P: PackedField,
{
	let (current, bind_hi) = values.split_at_mut(bind_offset);
	let reduce_offset = bind_offset / 2;
	let (lo, hi) = current.split_at_mut(reduce_offset);
	let (bind_lo, bind_hi) = bind_hi.split_at(reduce_offset);
	(lo, hi, bind_lo, bind_hi)
}

fn packed_bind_round<F, P>(buffers: &mut [FieldBuffer<P>; 3], prev_log_len: usize, challenge: F)
where
	F: Field,
	P: PackedField<Scalar = F>,
{
	debug_assert!(prev_log_len > P::LOG_WIDTH);
	let offset = 1 << (prev_log_len - 1 - P::LOG_WIDTH);
	let challenge = P::broadcast(challenge);
	for buffer in buffers {
		let (lo, hi) = buffer.as_mut().split_at_mut(offset);
		(lo, &hi[..offset])
			.into_par_iter()
			.for_each(|(lo_i, hi_i)| *lo_i += challenge * (*hi_i - *lo_i));
	}
}

fn sum_packed<P>(packed: P) -> P::Scalar
where
	P: PackedField,
{
	packed
		.iter()
		.fold(P::Scalar::ZERO, |sum, value| sum + value)
}

fn keccak_packed_fused_min_reduce_words() -> usize {
	std::env::var("KECCAK_PACKED_FUSED_MIN_REDUCE_WORDS")
		.ok()
		.and_then(|value| value.parse().ok())
		.unwrap_or(1 << 16)
}

fn prove_spartan_outer_from_folded_columns_fused_adaptive<F>(
	columns: FoldedOuterColumns<F>,
	first_round_challenge: F,
	zerocheck_challenges: Vec<F>,
	sumcheck_challenges: &[F],
	folded_claim: F,
) -> SpartanOuterPass<F>
where
	F: BinaryField + Field,
{
	assert_eq!(zerocheck_challenges.len(), columns.log_rows);
	assert_eq!(sumcheck_challenges.len(), columns.log_rows);

	let log_rows = columns.log_rows;
	let mut state = FusedOuterState::new(columns, &zerocheck_challenges);
	let mut round_messages = Vec::with_capacity(log_rows);
	let mut last_eval = folded_claim;
	let (mut rounds_done, mut live_pairs) = state.run_persistent_parallel_prefix(
		1usize << (log_rows - 1),
		&zerocheck_challenges,
		sumcheck_challenges,
		&mut last_eval,
		&mut round_messages,
	);

	while rounds_done < log_rows {
		let coeffs = state.sequential_round(
			rounds_done,
			live_pairs,
			zerocheck_challenges[log_rows - rounds_done - 1],
			sumcheck_challenges[rounds_done],
			&mut last_eval,
		);
		round_messages.push(coeffs);
		rounds_done += 1;
		live_pairs >>= 1;
	}

	let multilinear_evals = state.final_evals();
	let final_eval = multilinear_evals[0] * multilinear_evals[1] - multilinear_evals[2];

	SpartanOuterPass {
		first_round_challenge,
		folded_claim,
		zerocheck_challenges,
		round_messages,
		sumcheck_challenges: sumcheck_challenges.to_vec(),
		multilinear_evals,
		final_eval,
	}
}

#[derive(Clone, Copy, Debug, Default)]
struct OuterPartial<F> {
	y_1: F,
	y_inf: F,
}

impl<F> OuterPartial<F>
where
	F: Field,
{
	fn add_assign(&mut self, rhs: Self) {
		self.y_1 += rhs.y_1;
		self.y_inf += rhs.y_inf;
	}

	fn interpolate_eq(self, sum: F, alpha: F) -> RoundCoeffs<F> {
		let y_0 = (sum - self.y_1 * alpha) * (F::ONE - alpha).invert_or_zero();
		let c_0 = y_0;
		let c_2 = self.y_inf;
		let c_1 = self.y_1 - c_0 - c_2;
		RoundCoeffs(vec![c_0, c_1, c_2])
	}
}

struct PartialSlots<F>(Vec<UnsafeCell<OuterPartial<F>>>);

// SAFETY: each worker writes only its own slot, and worker 0 reads slots only after the phase
// barrier proves every active writer has completed.
unsafe impl<F: Send> Sync for PartialSlots<F> {}

struct SharedOuterColumns<F> {
	p: UnsafeCell<Vec<F>>,
	q: UnsafeCell<Vec<F>>,
	c: UnsafeCell<Vec<F>>,
	eq: UnsafeCell<Vec<F>>,
}

// SAFETY: all parallel methods partition pair ranges into disjoint windows. Each write goes to the
// low half of the current live prefix, and each read is either from that worker's just-written
// range or from the read-only high half for the current round.
unsafe impl<F: Send> Sync for SharedOuterColumns<F> {}

struct FusedOuterState<F> {
	columns: SharedOuterColumns<F>,
	log_rows: usize,
}

impl<F> FusedOuterState<F>
where
	F: BinaryField + Field,
{
	fn new(columns: FoldedOuterColumns<F>, zerocheck_challenges: &[F]) -> Self {
		let eq_weights = eq_ind_partial_eval_scalars(
			&zerocheck_challenges[..columns.log_rows.saturating_sub(1)],
		);
		Self {
			log_rows: columns.log_rows,
			columns: SharedOuterColumns {
				p: UnsafeCell::new(columns.p),
				q: UnsafeCell::new(columns.q),
				c: UnsafeCell::new(columns.c),
				eq: UnsafeCell::new(eq_weights),
			},
		}
	}

	fn run_persistent_parallel_prefix(
		&mut self,
		live_pairs_start: usize,
		zerocheck_challenges: &[F],
		sumcheck_challenges: &[F],
		last_eval: &mut F,
		round_messages: &mut Vec<RoundCoeffs<F>>,
	) -> (usize, usize) {
		let mut plans = Vec::new();
		let mut rounds_done = 0usize;
		let mut live_pairs = live_pairs_start;
		let mut workers = keccak_outer_max_workers().min(live_pairs).max(1);
		let min_pairs_per_worker = keccak_outer_min_pairs_per_worker();

		while workers >= 2 && rounds_done < self.log_rows {
			while workers >= 2 && live_pairs / workers < min_pairs_per_worker {
				workers /= 2;
			}
			if workers < 2 {
				break;
			}

			plans.push(ParallelRoundPlan {
				abs_round: rounds_done,
				live_pairs,
				workers,
			});
			rounds_done += 1;
			live_pairs >>= 1;
		}

		let max_workers = plans.iter().map(|plan| plan.workers).max().unwrap_or(1);
		if max_workers < 2 {
			return (0, live_pairs_start);
		}

		let partials = PartialSlots(
			(0..max_workers)
				.map(|_| UnsafeCell::new(OuterPartial::default()))
				.collect(),
		);
		let barrier = Barrier::new(max_workers);

		thread::scope(|scope| {
			for worker_idx in 1..max_workers {
				let partials = &partials;
				let barrier = &barrier;
				let state = &*self;
				let plans = &plans;
				scope.spawn(move || {
					for plan in plans {
						if worker_idx < plan.workers {
							let (lo, hi) =
								static_chunk_range(plan.live_pairs, plan.workers, worker_idx);
							let partial = if plan.abs_round == 0 {
								state.reduce_round(plan.live_pairs, lo, hi)
							} else {
								state.bind_then_reduce_round(
									plan.live_pairs,
									lo,
									hi,
									sumcheck_challenges[plan.abs_round - 1],
								)
							};
							// SAFETY: this worker is the unique writer for its slot.
							unsafe {
								*partials.0[worker_idx].get() = partial;
							}
						}
						barrier.wait();
						if worker_idx < plan.workers && plan.abs_round + 1 < state.log_rows {
							let next_eq_len = plan.live_pairs >> 1;
							let (lo, hi) =
								static_chunk_range(next_eq_len, plan.workers, worker_idx);
							state.truncate_eq_round(next_eq_len, lo, hi);
						}
						barrier.wait();
					}

					if let Some(last_plan) = plans.last()
						&& worker_idx < last_plan.workers
					{
						let (lo, hi) =
							static_chunk_range(last_plan.live_pairs, last_plan.workers, worker_idx);
						state.bind_round(
							last_plan.live_pairs,
							lo,
							hi,
							sumcheck_challenges[last_plan.abs_round],
						);
					}
				});
			}

			for plan in &plans {
				let (lo, hi) = static_chunk_range(plan.live_pairs, plan.workers, 0);
				let partial = if plan.abs_round == 0 {
					self.reduce_round(plan.live_pairs, lo, hi)
				} else {
					self.bind_then_reduce_round(
						plan.live_pairs,
						lo,
						hi,
						sumcheck_challenges[plan.abs_round - 1],
					)
				};
				// SAFETY: worker 0 is the unique writer for slot 0.
				unsafe {
					*partials.0[0].get() = partial;
				}
				barrier.wait();

				let mut sum = OuterPartial::default();
				for slot in &partials.0[..plan.workers] {
					// SAFETY: all workers reached the barrier after writing their slots.
					sum.add_assign(unsafe { *slot.get() });
				}
				let alpha = zerocheck_challenges[self.log_rows - plan.abs_round - 1];
				let coeffs = sum.interpolate_eq(*last_eval, alpha);
				*last_eval = coeffs.evaluate(sumcheck_challenges[plan.abs_round]);
				round_messages.push(coeffs);

				if plan.abs_round + 1 < self.log_rows {
					let next_eq_len = plan.live_pairs >> 1;
					let (lo, hi) = static_chunk_range(next_eq_len, plan.workers, 0);
					self.truncate_eq_round(next_eq_len, lo, hi);
				}
				barrier.wait();
			}

			if let Some(last_plan) = plans.last() {
				let (lo, hi) = static_chunk_range(last_plan.live_pairs, last_plan.workers, 0);
				self.bind_round(
					last_plan.live_pairs,
					lo,
					hi,
					sumcheck_challenges[last_plan.abs_round],
				);
			}
		});

		(rounds_done, live_pairs)
	}

	fn sequential_round(
		&mut self,
		abs_round: usize,
		live_pairs: usize,
		alpha: F,
		challenge: F,
		last_eval: &mut F,
	) -> RoundCoeffs<F> {
		let partial = self.reduce_round(live_pairs, 0, live_pairs);
		let coeffs = partial.interpolate_eq(*last_eval, alpha);
		*last_eval = coeffs.evaluate(challenge);
		if abs_round + 1 < self.log_rows {
			self.truncate_eq_round(live_pairs >> 1, 0, live_pairs >> 1);
		}
		self.bind_round(live_pairs, 0, live_pairs, challenge);
		coeffs
	}

	fn reduce_round(&self, live_pairs: usize, lo: usize, hi: usize) -> OuterPartial<F> {
		// SAFETY: reduce is read-only over the current live prefix.
		let p = unsafe { &*self.columns.p.get() };
		let q = unsafe { &*self.columns.q.get() };
		let c = unsafe { &*self.columns.c.get() };
		let eq_weights = unsafe { &*self.columns.eq.get() };
		let mut partial = OuterPartial::default();
		for i in lo..hi {
			let p_1 = p[i + live_pairs];
			let q_1 = q[i + live_pairs];
			let c_1 = c[i + live_pairs];
			let p_inf = p[i] + p_1;
			let q_inf = q[i] + q_1;
			let weight = eq_weights[i];
			partial.y_1 += (p_1 * q_1 - c_1) * weight;
			partial.y_inf += (p_inf * q_inf) * weight;
		}
		partial
	}

	fn bind_then_reduce_round(
		&self,
		live_pairs: usize,
		lo: usize,
		hi: usize,
		prev_challenge: F,
	) -> OuterPartial<F> {
		let p = unsafe { &mut *self.columns.p.get() };
		let q = unsafe { &mut *self.columns.q.get() };
		let c = unsafe { &mut *self.columns.c.get() };
		let eq_weights = unsafe { &*self.columns.eq.get() };
		let mut partial = OuterPartial::default();
		for i in lo..hi {
			let i_hi = i + live_pairs * 2;
			let p_0 = p[i] + prev_challenge * (p[i_hi] - p[i]);
			let q_0 = q[i] + prev_challenge * (q[i_hi] - q[i]);
			let c_0 = c[i] + prev_challenge * (c[i_hi] - c[i]);
			p[i] = p_0;
			q[i] = q_0;
			c[i] = c_0;

			let j = i + live_pairs;
			let j_hi = j + live_pairs * 2;
			let p_1 = p[j] + prev_challenge * (p[j_hi] - p[j]);
			let q_1 = q[j] + prev_challenge * (q[j_hi] - q[j]);
			let c_1 = c[j] + prev_challenge * (c[j_hi] - c[j]);
			p[j] = p_1;
			q[j] = q_1;
			c[j] = c_1;

			let weight = eq_weights[i];
			partial.y_1 += (p_1 * q_1 - c_1) * weight;
			partial.y_inf += ((p_0 + p_1) * (q_0 + q_1)) * weight;
		}
		partial
	}

	fn bind_round(&self, live_pairs: usize, lo: usize, hi: usize, challenge: F) {
		let p = unsafe { &mut *self.columns.p.get() };
		let q = unsafe { &mut *self.columns.q.get() };
		let c = unsafe { &mut *self.columns.c.get() };
		for i in lo..hi {
			let p_lo = p[i];
			let q_lo = q[i];
			let c_lo = c[i];
			p[i] = p_lo + challenge * (p[i + live_pairs] - p_lo);
			q[i] = q_lo + challenge * (q[i + live_pairs] - q_lo);
			c[i] = c_lo + challenge * (c[i + live_pairs] - c_lo);
		}
	}

	fn truncate_eq_round(&self, next_eq_len: usize, lo: usize, hi: usize) {
		let eq_weights = unsafe { &mut *self.columns.eq.get() };
		for i in lo..hi {
			let high = eq_weights[i + next_eq_len];
			eq_weights[i] += high;
		}
	}

	fn final_evals(&self) -> [F; 3] {
		[
			unsafe { (&*self.columns.p.get())[0] },
			unsafe { (&*self.columns.q.get())[0] },
			unsafe { (&*self.columns.c.get())[0] },
		]
	}
}

#[derive(Clone, Copy)]
struct ParallelRoundPlan {
	abs_round: usize,
	live_pairs: usize,
	workers: usize,
}

fn static_chunk_range(live_pairs: usize, workers: usize, worker_idx: usize) -> (usize, usize) {
	let chunk = live_pairs.div_ceil(workers);
	let lo = (worker_idx * chunk).min(live_pairs);
	let hi = ((worker_idx + 1) * chunk).min(live_pairs);
	(lo, hi)
}

fn keccak_outer_max_workers() -> usize {
	std::env::var("KECCAK_OUTER_WORKERS")
		.or_else(|_| std::env::var("RAYON_NUM_THREADS"))
		.ok()
		.and_then(|value| value.parse().ok())
		.or_else(|| thread::available_parallelism().ok().map(usize::from))
		.unwrap_or(1)
}

fn keccak_outer_min_pairs_per_worker() -> usize {
	std::env::var("KECCAK_OUTER_MIN_PAIRS_PER_WORKER")
		.ok()
		.and_then(|value| value.parse().ok())
		.unwrap_or(1 << 13)
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
	use binius_math::multilinear::eq::eq_ind_partial_eval_scalars;
	use binius_prover::OptimalPackedB128;
	use binius_transcript::{ProverTranscript, fiat_shamir::HasherChallenger};
	use rand::{Rng, SeedableRng, rngs::StdRng};

	use crate::{
		bit_ntt::NttLookup,
		constants::N_ROUNDS,
		trace::{PermutationTrace, State},
	};

	use super::*;

	type B128 = BinaryField128bGhash;
	type P = PackedAESBinaryField16x8b;
	type StdChallenger = HasherChallenger<sha2::Sha256>;

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

	#[test]
	fn first_round_claim_from_row_challenges_matches_manual_weights() {
		let mut rng = StdRng::seed_from_u64(18);
		let lookup = NttLookup::<P>::for_upper_half_domain();
		let traces: Vec<_> = (0..2)
			.map(|_| PermutationTrace::new(rng.random::<State>()))
			.collect();
		let round_traces: Vec<_> = traces.iter().flat_map(|trace| trace.rounds).collect();
		let row_count = round_traces.len() * N_LANES;
		let log_rows = row_count.next_power_of_two().ilog2() as usize;
		let row_challenges: Vec<_> = (0..log_rows).map(|_| B128::random(&mut rng)).collect();
		let first_round_challenge = B128::random(&mut rng);
		let eq_weights = row_eq_weights(row_count, &row_challenges);
		let upper_half =
			par_upper_half_round_message::<B128, P>(&lookup, &round_traces, &eq_weights);
		let prover_message_domain =
			BinarySubspace::<AESTowerField8b>::with_dim(crate::bit_ntt::LOG_LANE_BITS + 1)
				.isomorphic::<B128>();

		assert_eq!(
			par_first_round_claim_from_row_challenges::<B128, P>(
				&lookup,
				&round_traces,
				&row_challenges,
				first_round_challenge,
			),
			next_sum_claim_from_upper_half(
				&upper_half,
				first_round_challenge,
				&prover_message_domain
			)
		);
	}

	#[test]
	fn spartan_outer_pass_runs_after_univariate_skip() {
		let mut rng = StdRng::seed_from_u64(16);
		let traces: Vec<_> = (0..2)
			.map(|_| PermutationTrace::new(rng.random::<State>()))
			.collect();
		let round_traces: Vec<_> = traces.iter().flat_map(|trace| trace.rounds).collect();
		let first_round_challenge = B128::random(&mut rng);

		let columns = folded_outer_columns::<B128, P>(&round_traces, first_round_challenge);
		assert_eq!(
			columns,
			folded_outer_columns_direct::<B128, P>(&round_traces, first_round_challenge)
		);
		assert_eq!(columns.row_count, round_traces.len() * N_LANES);
		assert!(columns.p.len().is_power_of_two());
		assert_eq!(columns.p.len(), 1 << columns.log_rows);

		let zerocheck_challenges: Vec<_> = (0..columns.log_rows)
			.map(|_| B128::random(&mut rng))
			.collect();
		let eq_weights = eq_ind_partial_eval_scalars(&zerocheck_challenges);
		let manual_claim: B128 = (0..columns.p.len())
			.map(|i| (columns.p[i] * columns.q[i] - columns.c[i]) * eq_weights[i])
			.sum();
		assert_eq!(folded_outer_claim(&columns, &zerocheck_challenges), manual_claim);

		let sumcheck_challenges: Vec<_> = (0..columns.log_rows)
			.map(|_| B128::random(&mut rng))
			.collect();
		let pass = prove_spartan_outer_after_first_round::<B128, P>(
			&round_traces,
			first_round_challenge,
			zerocheck_challenges,
			&sumcheck_challenges,
		)
		.unwrap();

		assert_eq!(pass.folded_claim, manual_claim);
		assert_eq!(pass.round_messages.len(), columns.log_rows);

		let mut claim = pass.folded_claim;
		for (round_idx, (round_message, &challenge)) in
			iter::zip(&pass.round_messages, &pass.sumcheck_challenges).enumerate()
		{
			let alpha = pass.zerocheck_challenges[columns.log_rows - 1 - round_idx];
			assert_eq!(
				claim,
				round_message.evaluate(B128::ZERO) * (B128::ONE - alpha)
					+ round_message.evaluate(B128::ONE) * alpha
			);
			claim = round_message.evaluate(challenge);
		}
		assert_eq!(claim, pass.final_eval);
	}

	#[test]
	fn fused_adaptive_spartan_outer_matches_generic_mlecheck() {
		let mut rng = StdRng::seed_from_u64(19);
		let traces: Vec<_> = (0..4)
			.map(|_| PermutationTrace::new(rng.random::<State>()))
			.collect();
		let round_traces: Vec<_> = traces.iter().flat_map(|trace| trace.rounds).collect();
		let first_round_challenge = B128::random(&mut rng);
		let columns = folded_outer_columns::<B128, P>(&round_traces, first_round_challenge);
		let zerocheck_challenges: Vec<_> = (0..columns.log_rows)
			.map(|_| B128::random(&mut rng))
			.collect();
		let folded_claim = folded_outer_claim(&columns, &zerocheck_challenges);
		let sumcheck_challenges: Vec<_> = (0..columns.log_rows)
			.map(|_| B128::random(&mut rng))
			.collect();

		let generic = prove_spartan_outer_after_first_round_with_claim::<B128, P>(
			&round_traces,
			first_round_challenge,
			zerocheck_challenges.clone(),
			&sumcheck_challenges,
			folded_claim,
		)
		.unwrap();
		let fused = prove_spartan_outer_after_first_round_with_claim_fused_adaptive::<B128, P>(
			&round_traces,
			first_round_challenge,
			zerocheck_challenges,
			&sumcheck_challenges,
			folded_claim,
		);

		assert_eq!(fused, generic);
	}

	#[test]
	fn packed_spartan_outer_matches_generic_mlecheck() {
		let mut rng = StdRng::seed_from_u64(20);
		let traces: Vec<_> = (0..8)
			.map(|_| PermutationTrace::new(rng.random::<State>()))
			.collect();
		let round_traces: Vec<_> = traces.iter().flat_map(|trace| trace.rounds).collect();
		let first_round_challenge = B128::random(&mut rng);
		let columns = folded_outer_columns::<B128, P>(&round_traces, first_round_challenge);
		let zerocheck_challenges: Vec<_> = (0..columns.log_rows)
			.map(|_| B128::random(&mut rng))
			.collect();
		let folded_claim = folded_outer_claim(&columns, &zerocheck_challenges);
		let sumcheck_challenges: Vec<_> = (0..columns.log_rows)
			.map(|_| B128::random(&mut rng))
			.collect();

		let generic = prove_spartan_outer_from_folded_columns_with_claim(
			columns.clone(),
			first_round_challenge,
			zerocheck_challenges.clone(),
			&sumcheck_challenges,
			folded_claim,
		)
		.unwrap();
		let packed_columns = pack_folded_outer_columns::<B128, OptimalPackedB128>(columns.clone());
		let packed_generic =
			prove_spartan_outer_from_folded_columns_with_claim_packed::<B128, OptimalPackedB128>(
				columns.clone(),
				first_round_challenge,
				zerocheck_challenges.clone(),
				&sumcheck_challenges,
				folded_claim,
			)
			.unwrap();
		let prepacked_generic =
			prove_spartan_outer_from_packed_folded_columns_with_claim::<B128, OptimalPackedB128>(
				packed_columns.clone(),
				first_round_challenge,
				zerocheck_challenges.clone(),
				&sumcheck_challenges,
				folded_claim,
			)
			.unwrap();
		let packed_fused = prove_spartan_outer_from_folded_columns_with_claim_packed_fused::<
			B128,
			OptimalPackedB128,
		>(
			columns,
			first_round_challenge,
			zerocheck_challenges.clone(),
			&sumcheck_challenges,
			folded_claim,
		)
		.unwrap();
		let prepacked_fused = prove_spartan_outer_from_packed_folded_columns_with_claim_fused::<
			B128,
			OptimalPackedB128,
		>(
			packed_columns.clone(),
			first_round_challenge,
			zerocheck_challenges.clone(),
			&sumcheck_challenges,
			folded_claim,
		)
		.unwrap();
		let prepacked_persistent_fused =
			prove_spartan_outer_from_packed_folded_columns_with_claim_persistent_fused::<
				B128,
				OptimalPackedB128,
			>(
				packed_columns,
				first_round_challenge,
				zerocheck_challenges,
				&sumcheck_challenges,
				folded_claim,
			)
			.unwrap();

		assert_eq!(packed_generic, generic);
		assert_eq!(prepacked_generic, generic);
		assert_eq!(packed_fused, generic);
		assert_eq!(prepacked_fused, generic);
		assert_eq!(prepacked_persistent_fused, generic);
	}

	#[test]
	fn transcripted_spartan_outer_verifier_replays() {
		let mut rng = StdRng::seed_from_u64(17);
		let lookup = NttLookup::<P>::for_upper_half_domain();
		let traces: Vec<_> = (0..2)
			.map(|_| PermutationTrace::new(rng.random::<State>()))
			.collect();
		let round_traces: Vec<_> = traces.iter().flat_map(|trace| trace.rounds).collect();
		let log_rows = (round_traces.len() * N_LANES).next_power_of_two().ilog2() as usize;
		let zerocheck_challenges: Vec<_> = (0..log_rows).map(|_| B128::random(&mut rng)).collect();

		let mut prover_transcript = ProverTranscript::<StdChallenger>::default();
		let prove_output = prove_spartan_outer_with_channel::<B128, P, _>(
			&lookup,
			&round_traces,
			zerocheck_challenges.clone(),
			&mut prover_transcript,
		)
		.unwrap();

		let mut verifier_transcript = prover_transcript.into_verifier();
		let round_message_domain =
			BinarySubspace::<AESTowerField8b>::with_dim(crate::bit_ntt::LOG_LANE_BITS + 1)
				.isomorphic::<B128>();
		let verify_output = verify_spartan_outer_with_channel(
			&zerocheck_challenges,
			&mut verifier_transcript,
			&round_message_domain,
		)
		.unwrap();

		assert_eq!(prove_output, verify_output);
	}
}
