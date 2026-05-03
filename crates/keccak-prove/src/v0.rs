// Copyright 2026 The Binius Developers

//! End-to-end v0 production-path constraint system for committed Keccak traces.

use binius_core::{
	ShiftVariant,
	constraint_system::{
		AndConstraint, ConstraintSystem, ShiftedValueIndex, ValueIndex, ValueVec, ValueVecLayout,
	},
	word::Word,
};
use binius_field::{BinaryField, FieldOps, util::powers};
use binius_math::{
	BinarySubspace,
	inner_product::inner_product_scalars,
	multilinear::{eq::eq_ind_partial_eval_scalars, evaluate::evaluate_inplace_scalars},
	univariate::lagrange_evals_scalars,
};
use binius_verifier::{
	config::{LOG_WORD_SIZE_BITS, WORD_SIZE_BITS},
	protocols::shift::{BITAND_ARITY, OperatorData, SHIFT_VARIANT_COUNT, evaluate_h_op},
};

use crate::{
	bit_ntt::CHI_OPERAND_LANES,
	constants::{N_LANES, N_ROUNDS, RHO_OFFSETS, ROUND_CONSTANTS, lane_x},
	layout,
	shift_operands::{
		RHO_PI_PREIMAGE, d_correctness_operand_at, value_index_at, virtual_b_operand_at,
	},
	witness::CommittedKeccakWitness,
};

/// Public constant index for the all-one word used by `P = 1 + B[x+1,y]`.
pub const ALL_ONE_INDEX: usize = 0;

/// Public constant offset for the 24 Keccak round constants.
pub const ROUND_CONSTANTS_OFFSET: usize = ALL_ONE_INDEX + 1;

/// Number of explicit public constants in the v0 production constraint system.
pub const N_CONSTANTS: usize = ROUND_CONSTANTS_OFFSET + N_ROUNDS;

/// Padded public segment size required by production Binius64.
pub const PUBLIC_WORDS: usize = 32;

/// Number of BitAnd constraint slots per `(permutation, round)` block.
pub const CONSTRAINT_ROWS_PER_ROUND: usize = 32;

/// First slot of the five `D` correctness rows inside a round constraint block.
pub const D_CONSTRAINT_OFFSET: usize = N_LANES;

const _: () = assert!(N_CONSTANTS <= PUBLIC_WORDS);

/// Return the value-vector base offset of the committed Keccak `A`/`D` witness.
#[inline]
pub const fn committed_witness_base() -> usize {
	PUBLIC_WORDS
}

/// Return the public constants used by the v0 production constraint system.
pub fn constants() -> Vec<Word> {
	let mut constants = Vec::with_capacity(N_CONSTANTS);
	constants.push(Word::ALL_ONE);
	constants.extend(ROUND_CONSTANTS.into_iter().map(Word));
	constants
}

/// Value-vector layout used by the v0 production constraint system.
pub fn value_vec_layout(n_permutations: usize) -> ValueVecLayout {
	let witness_len = layout::witness_len(n_permutations);
	let committed_total_len = (committed_witness_base() + witness_len).next_power_of_two();
	let n_internal = committed_total_len - committed_witness_base() - witness_len;

	ValueVecLayout {
		n_const: N_CONSTANTS,
		n_inout: 0,
		n_witness: witness_len,
		n_internal,
		offset_inout: PUBLIC_WORDS,
		offset_witness: committed_witness_base(),
		committed_total_len,
		n_scratch: 0,
	}
}

/// Return a production `ValueVec` containing public constants and the committed Keccak witness.
pub fn value_vec(witness: &CommittedKeccakWitness) -> ValueVec {
	let layout = value_vec_layout(witness.n_permutations());

	let mut public = vec![Word::ZERO; committed_witness_base()];
	public[..N_CONSTANTS].copy_from_slice(&constants());

	let private_len = layout.committed_total_len - committed_witness_base();
	let mut private = vec![Word::ZERO; private_len];
	private[..witness.words().len()].copy_from_slice(witness.words());

	ValueVec::new_from_data(layout, public, private)
		.expect("constructed Keccak v0 value vector has matching layout")
}

/// Number of active chi constraints for `n_permutations` Keccak-f[1600] traces.
#[inline]
pub const fn chi_row_count(n_permutations: usize) -> usize {
	n_permutations * N_ROUNDS * N_LANES
}

/// Number of active chi plus `D` correctness rows for `n_permutations`.
#[inline]
pub const fn active_row_count(n_permutations: usize) -> usize {
	n_permutations * N_ROUNDS * (N_LANES + layout::D_WORDS_PER_BLOCK)
}

/// Number of v0 constraint slots for `n_permutations`, including two padding rows per round.
#[inline]
pub const fn constraint_row_count(n_permutations: usize) -> usize {
	n_permutations * N_ROUNDS * CONSTRAINT_ROWS_PER_ROUND
}

/// Build full chi/iota AND constraints over committed `A`/`D` witness terms.
///
/// Each row is the production BitAnd relation:
///
/// ```text
/// (1 + B[x+1,y]) & B[x+2,y] = B[x,y] + A_next[x,y] + iota
/// ```
///
/// where every virtual `B` lane is lowered to shifted committed `A` and `D` terms.
pub fn chi_constraints_at(base_index: usize, n_permutations: usize) -> Vec<AndConstraint> {
	let mut constraints = Vec::with_capacity(chi_row_count(n_permutations));

	for permutation in 0..n_permutations {
		for round in 0..N_ROUNDS {
			for lane_idx in 0..N_LANES {
				constraints.push(chi_constraint_at(base_index, permutation, round, lane_idx));
			}
		}
	}

	constraints
}

/// Build full chi/iota AND constraints at the canonical committed witness base.
#[inline]
pub fn chi_constraints(n_permutations: usize) -> Vec<AndConstraint> {
	chi_constraints_at(committed_witness_base(), n_permutations)
}

/// Build the v0 production constraint system for a batch of committed Keccak traces.
pub fn constraint_system(n_permutations: usize) -> ConstraintSystem {
	let mut and_constraints = Vec::with_capacity(constraint_row_count(n_permutations));
	for permutation in 0..n_permutations {
		for round in 0..N_ROUNDS {
			for lane_idx in 0..N_LANES {
				and_constraints.push(chi_constraint_at(
					committed_witness_base(),
					permutation,
					round,
					lane_idx,
				));
			}

			for x in 0..layout::D_WORDS_PER_BLOCK {
				and_constraints.push(AndConstraint {
					a: Vec::new(),
					b: Vec::new(),
					c: d_correctness_operand_at(committed_witness_base(), permutation, round, x),
				});
			}

			and_constraints.extend([AndConstraint::default(), AndConstraint::default()]);
		}
	}

	let mut constraint_system = ConstraintSystem::new(
		constants(),
		value_vec_layout(n_permutations),
		and_constraints,
		Vec::new(),
	);
	constraint_system
		.validate_and_prepare()
		.expect("constructed Keccak v0 production constraint system is valid");
	constraint_system
}

fn chi_constraint_at(
	base_index: usize,
	permutation: usize,
	round: usize,
	lane_idx: usize,
) -> AndConstraint {
	let (p_lane, q_lane, r_lane) = CHI_OPERAND_LANES[lane_idx];

	let mut a = vec![ShiftedValueIndex::plain(ValueIndex(ALL_ONE_INDEX as u32))];
	a.extend(virtual_b_operand_at(base_index, permutation, round, p_lane));

	let mut c = virtual_b_operand_at(base_index, permutation, round, r_lane);
	c.push(ShiftedValueIndex::plain(value_index_at(
		base_index,
		layout::a_index(permutation, round + 1, lane_idx),
	)));
	if lane_idx == 0 {
		c.push(ShiftedValueIndex::plain(ValueIndex((ROUND_CONSTANTS_OFFSET + round) as u32)));
	}

	AndConstraint {
		a,
		b: virtual_b_operand_at(base_index, permutation, round, q_lane),
		c,
	}
}

/// Evaluate the v0 BitAnd monster multilinear with Keccak's tensor-structured schema.
///
/// This is a verifier-side specialization of production Shift's generic matrix evaluation. It
/// avoids materializing operand vectors and streams the fixed Keccak chi and `D` correctness
/// pattern directly into the `(shift_variant, amount)` matrix evaluations.
pub fn structured_bitand_monster_eval<F, E>(
	n_permutations: usize,
	operator_data: &OperatorData<E, BITAND_ARITY>,
	subspace: &BinarySubspace<F>,
	lambda: E,
	r_j: &[E],
	r_s: &[E],
	r_y: &[E],
) -> E
where
	F: BinaryField,
	E: FieldOps<Scalar = F> + From<F>,
{
	assert_eq!(subspace.dim(), LOG_WORD_SIZE_BITS);

	let r_x_prime_tensor = eq_ind_partial_eval_scalars(&operator_data.r_x_prime);
	let r_y_tensor = eq_ind_partial_eval_scalars(r_y);
	let l_tilde = lagrange_evals_scalars(subspace, operator_data.r_zhat_prime.clone());
	let h_op_evals = evaluate_h_op(&l_tilde, r_j, r_s);
	let lambda_powers = powers(lambda)
		.skip(1)
		.take(BITAND_ARITY)
		.collect::<Vec<_>>();

	let matrix_evals = structured_bitand_matrix_evals(
		n_permutations,
		&lambda_powers,
		&r_x_prime_tensor,
		&r_y_tensor,
	);

	inner_product_scalars(
		matrix_evals
			.into_iter()
			.map(|mut evals_op| evaluate_inplace_scalars(&mut evals_op[..], r_s)),
		h_op_evals,
	)
}

fn structured_bitand_matrix_evals<F, E>(
	n_permutations: usize,
	operand_coeffs: &[E],
	r_x_prime_tensor: &[E],
	r_y_tensor: &[E],
) -> [[E; WORD_SIZE_BITS]; SHIFT_VARIANT_COUNT]
where
	F: BinaryField,
	E: FieldOps<Scalar = F> + From<F>,
{
	assert_eq!(operand_coeffs.len(), BITAND_ARITY);

	let mut evals =
		std::array::from_fn(|_| std::array::from_fn::<E, WORD_SIZE_BITS, _>(|_| E::zero()));

	for permutation in 0..n_permutations {
		for round in 0..N_ROUNDS {
			for lane_idx in 0..N_LANES {
				let row_weight =
					&r_x_prime_tensor[chi_constraint_row_index(permutation, round, lane_idx)];
				let (p_lane, q_lane, r_lane) = CHI_OPERAND_LANES[lane_idx];

				accumulate_plain_term::<F, E>(
					&mut evals,
					&operand_coeffs[0],
					row_weight,
					ALL_ONE_INDEX,
					r_y_tensor,
				);
				accumulate_virtual_b_terms::<F, E>(
					&mut evals,
					&operand_coeffs[0],
					row_weight,
					permutation,
					round,
					p_lane,
					r_y_tensor,
				);
				accumulate_virtual_b_terms::<F, E>(
					&mut evals,
					&operand_coeffs[1],
					row_weight,
					permutation,
					round,
					q_lane,
					r_y_tensor,
				);
				accumulate_virtual_b_terms::<F, E>(
					&mut evals,
					&operand_coeffs[2],
					row_weight,
					permutation,
					round,
					r_lane,
					r_y_tensor,
				);
				accumulate_plain_term::<F, E>(
					&mut evals,
					&operand_coeffs[2],
					row_weight,
					committed_witness_base() + layout::a_index(permutation, round + 1, lane_idx),
					r_y_tensor,
				);
				if lane_idx == 0 {
					accumulate_plain_term::<F, E>(
						&mut evals,
						&operand_coeffs[2],
						row_weight,
						ROUND_CONSTANTS_OFFSET + round,
						r_y_tensor,
					);
				}
			}

			for x in 0..layout::D_WORDS_PER_BLOCK {
				let row_weight =
					&r_x_prime_tensor[d_correctness_constraint_row_index(permutation, round, x)];
				accumulate_plain_term::<F, E>(
					&mut evals,
					&operand_coeffs[2],
					row_weight,
					committed_witness_base() + layout::d_index(permutation, round, x),
					r_y_tensor,
				);

				let left_x = (x + 4) % 5;
				for y in 0..5 {
					accumulate_plain_term::<F, E>(
						&mut evals,
						&operand_coeffs[2],
						row_weight,
						committed_witness_base()
							+ layout::a_index(permutation, round, left_x + 5 * y),
						r_y_tensor,
					);
				}

				let right_x = (x + 1) % 5;
				for y in 0..5 {
					accumulate_rotl_term::<F, E>(
						&mut evals,
						&operand_coeffs[2],
						row_weight,
						committed_witness_base()
							+ layout::a_index(permutation, round, right_x + 5 * y),
						1,
						r_y_tensor,
					);
				}
			}
		}
	}

	evals
}

#[inline]
fn chi_constraint_row_index(permutation: usize, round: usize, lane_idx: usize) -> usize {
	constraint_row_index(permutation, round, lane_idx)
}

#[inline]
fn d_correctness_constraint_row_index(permutation: usize, round: usize, x: usize) -> usize {
	constraint_row_index(permutation, round, D_CONSTRAINT_OFFSET + x)
}

#[inline]
fn constraint_row_index(permutation: usize, round: usize, slot: usize) -> usize {
	(permutation * N_ROUNDS + round) * CONSTRAINT_ROWS_PER_ROUND + slot
}

#[inline]
fn accumulate_virtual_b_terms<F, E>(
	evals: &mut [[E; WORD_SIZE_BITS]; SHIFT_VARIANT_COUNT],
	operand_coeff: &E,
	row_weight: &E,
	permutation: usize,
	round: usize,
	b_lane: usize,
	r_y_tensor: &[E],
) where
	F: BinaryField,
	E: FieldOps<Scalar = F> + From<F>,
{
	let source_lane = RHO_PI_PREIMAGE[b_lane];
	let source_x = lane_x(source_lane);
	let rho = RHO_OFFSETS[source_lane];
	accumulate_rotl_term::<F, E>(
		evals,
		operand_coeff,
		row_weight,
		committed_witness_base() + layout::a_index(permutation, round, source_lane),
		rho,
		r_y_tensor,
	);
	accumulate_rotl_term::<F, E>(
		evals,
		operand_coeff,
		row_weight,
		committed_witness_base() + layout::d_index(permutation, round, source_x),
		rho,
		r_y_tensor,
	);
}

#[inline]
fn accumulate_rotl_term<F, E>(
	evals: &mut [[E; WORD_SIZE_BITS]; SHIFT_VARIANT_COUNT],
	operand_coeff: &E,
	row_weight: &E,
	value_index: usize,
	amount: u32,
	r_y_tensor: &[E],
) where
	F: BinaryField,
	E: FieldOps<Scalar = F> + From<F>,
{
	let amount = (amount % 64) as usize;
	if amount == 0 {
		accumulate_plain_term::<F, E>(evals, operand_coeff, row_weight, value_index, r_y_tensor);
	} else {
		accumulate_shifted_term::<F, E>(
			evals,
			operand_coeff,
			row_weight,
			value_index,
			ShiftVariant::Rotr,
			64 - amount,
			r_y_tensor,
		);
	}
}

#[inline]
fn accumulate_plain_term<F, E>(
	evals: &mut [[E; WORD_SIZE_BITS]; SHIFT_VARIANT_COUNT],
	operand_coeff: &E,
	row_weight: &E,
	value_index: usize,
	r_y_tensor: &[E],
) where
	F: BinaryField,
	E: FieldOps<Scalar = F> + From<F>,
{
	accumulate_shifted_term::<F, E>(
		evals,
		operand_coeff,
		row_weight,
		value_index,
		ShiftVariant::Sll,
		0,
		r_y_tensor,
	);
}

#[inline]
fn accumulate_shifted_term<F, E>(
	evals: &mut [[E; WORD_SIZE_BITS]; SHIFT_VARIANT_COUNT],
	operand_coeff: &E,
	row_weight: &E,
	value_index: usize,
	shift_variant: ShiftVariant,
	amount: usize,
	r_y_tensor: &[E],
) where
	F: BinaryField,
	E: FieldOps<Scalar = F> + From<F>,
{
	let shift_id = match shift_variant {
		ShiftVariant::Sll => 0,
		ShiftVariant::Slr => 1,
		ShiftVariant::Sar => 2,
		ShiftVariant::Rotr => 3,
		ShiftVariant::Sll32 => 4,
		ShiftVariant::Srl32 => 5,
		ShiftVariant::Sra32 => 6,
		ShiftVariant::Rotr32 => 7,
	};
	evals[shift_id][amount] +=
		operand_coeff.clone() * row_weight.clone() * r_y_tensor[value_index].clone();
}

#[cfg(test)]
mod tests {
	use binius_core::{constraint_system::Operand, verify::verify_constraints};
	use binius_field::{BinaryField128bGhash, Random, arch::OptimalPackedB128};
	use binius_math::BinarySubspace;
	use binius_prover::{Prover, hash::parallel_compression::ParallelCompressionAdaptor};
	use binius_transcript::ProverTranscript;
	use binius_verifier::{
		Verifier,
		config::StdChallenger,
		hash::{StdCompression, StdDigest},
		protocols::shift::{
			OperatorData as VerifierOperatorData, evaluate_monster_multilinear_for_operation,
		},
	};
	use rand::{Rng, SeedableRng, rngs::StdRng};

	use crate::{
		trace::{PermutationTrace, State},
		witness::CommittedKeccakWitness,
	};

	use super::*;

	#[test]
	fn v0_constraint_system_satisfies_native_keccak_traces() {
		let mut rng = StdRng::seed_from_u64(33);
		let traces: Vec<_> = (0..2)
			.map(|_| PermutationTrace::new(rng.random::<State>()))
			.collect();
		let witness = CommittedKeccakWitness::from_traces(&traces);
		let cs = constraint_system(traces.len());
		let value_vec = value_vec(&witness);

		verify_constraints(&cs, &value_vec).unwrap();
		assert_eq!(value_vec.public()[ALL_ONE_INDEX], Word::ALL_ONE);
		for (round, &round_constant) in ROUND_CONSTANTS.iter().enumerate() {
			assert_eq!(value_vec.public()[ROUND_CONSTANTS_OFFSET + round], Word(round_constant));
		}
	}

	#[test]
	fn v0_production_prover_verifies_small_trace() {
		const LOG_INV_RATE: usize = 1;

		let mut rng = StdRng::seed_from_u64(34);
		let trace = PermutationTrace::new(rng.random::<State>());
		let witness = CommittedKeccakWitness::from_traces(std::slice::from_ref(&trace));
		let cs = constraint_system(1);
		let value_vec = value_vec(&witness);

		let verifier =
			Verifier::<StdDigest, _>::setup(cs, LOG_INV_RATE, StdCompression::default()).unwrap();
		let prover = Prover::<OptimalPackedB128, _, StdDigest>::setup(
			verifier.clone(),
			ParallelCompressionAdaptor::new(StdCompression::default()),
		)
		.unwrap();

		let mut prover_transcript = ProverTranscript::new(StdChallenger::default());
		prover
			.prove(value_vec.clone(), &mut prover_transcript)
			.unwrap();

		let mut verifier_transcript = prover_transcript.into_verifier();
		verifier
			.verify(value_vec.public(), &mut verifier_transcript)
			.unwrap();
		verifier_transcript.finalize().unwrap();
	}

	#[test]
	fn structured_monster_eval_matches_generic_shift_verifier() {
		type F = BinaryField128bGhash;

		let mut rng = StdRng::seed_from_u64(35);
		let n_permutations = 2;
		let cs = constraint_system(n_permutations);
		let subspace = BinarySubspace::<F>::with_dim(LOG_WORD_SIZE_BITS);
		let operator_data = VerifierOperatorData::new(
			F::random(&mut rng),
			(0..cs.and_constraints.len().ilog2())
				.map(|_| F::random(&mut rng))
				.collect(),
			std::array::from_fn(|_| F::random(&mut rng)),
		);
		let lambda = F::random(&mut rng);
		let r_j = (0..LOG_WORD_SIZE_BITS)
			.map(|_| F::random(&mut rng))
			.collect::<Vec<_>>();
		let r_s = (0..LOG_WORD_SIZE_BITS)
			.map(|_| F::random(&mut rng))
			.collect::<Vec<_>>();
		let r_y = (0..cs.value_vec_layout.committed_total_len.ilog2())
			.map(|_| F::random(&mut rng))
			.collect::<Vec<_>>();

		let mut a = Vec::<&Operand>::with_capacity(cs.and_constraints.len());
		let mut b = Vec::<&Operand>::with_capacity(cs.and_constraints.len());
		let mut c = Vec::<&Operand>::with_capacity(cs.and_constraints.len());
		for constraint in &cs.and_constraints {
			a.push(&constraint.a);
			b.push(&constraint.b);
			c.push(&constraint.c);
		}
		let generic = evaluate_monster_multilinear_for_operation(
			&[a, b, c],
			&operator_data,
			&subspace,
			lambda,
			&r_j,
			&r_s,
			&r_y,
		)
		.unwrap();
		let structured = structured_bitand_monster_eval(
			n_permutations,
			&operator_data,
			&subspace,
			lambda,
			&r_j,
			&r_s,
			&r_y,
		);

		assert_eq!(structured, generic);
	}
}
