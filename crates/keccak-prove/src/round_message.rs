// Copyright 2026 The Binius Developers

//! Upper-half sumcheck round-message accumulation for Keccak chi/iota constraints.

use std::iter;

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
	multilinear::eq::eq_ind_partial_eval_scalars,
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
	assert_eq!(zerocheck_challenges.len(), columns.log_rows);
	assert_eq!(sumcheck_challenges.len(), columns.log_rows);

	let log_rows = columns.log_rows;
	let mut prover = QuadraticMleCheckProver::new(
		columns.into_field_buffers(),
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

	let multilinear_evals: [FChallenge; 3] =
		prover.finish()?.try_into().expect("three folded columns");
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
