// Copyright 2026 The Binius Developers

//! Constraint-system schemas that hand Keccak chi and theta-correction claims to Shift.

use binius_core::{
	constraint_system::{
		AndConstraint, ConstraintSystem, ShiftedValueIndex, ValueVec, ValueVecLayout,
	},
	consts::MIN_WORDS_PER_SEGMENT,
	word::Word,
};
use binius_field::{BinaryField, Field};
use binius_math::{
	BinarySubspace, multilinear::eq::eq_ind_partial_eval_scalars, univariate::lagrange_evals,
};

use crate::{
	bit_ntt::CHI_OPERAND_LANES,
	constants::{N_LANES, N_ROUNDS, ROUND_CONSTANTS},
	layout,
	shift_operands::{d_correctness_operand_at, value_index_at, virtual_b_operand_at},
	witness::CommittedKeccakWitness,
};

/// Number of public words kept before the committed Keccak witness segment.
///
/// Shift's production value-vector layout requires a power-of-two public segment of at least this
/// size. The current Keccak operand schemas treat both public words as zero constants and keep all
/// Keccak state data in the private witness segment beginning at [`committed_witness_base`].
pub const SHIFT_PUBLIC_WORDS: usize = MIN_WORDS_PER_SEGMENT;

/// Return the value-vector base offset of the committed Keccak `A`/`D` witness.
#[inline]
pub const fn committed_witness_base() -> usize {
	SHIFT_PUBLIC_WORDS
}

/// Number of active chi constraints for `n_permutations` Keccak-f[1600] traces.
#[inline]
pub const fn chi_row_count(n_permutations: usize) -> usize {
	n_permutations * N_ROUNDS * N_LANES
}

/// Number of active `D` correctness constraints for `n_permutations` Keccak-f[1600] traces.
#[inline]
pub const fn d_correctness_row_count(n_permutations: usize) -> usize {
	n_permutations * N_ROUNDS * layout::D_WORDS_PER_BLOCK
}

/// Value-vector layout used by the Keccak Shift schemas.
pub fn value_vec_layout(n_permutations: usize) -> ValueVecLayout {
	let witness_len = layout::witness_len(n_permutations);
	let committed_total_len = (committed_witness_base() + witness_len).next_power_of_two();
	let n_internal = committed_total_len - committed_witness_base() - witness_len;

	ValueVecLayout {
		n_const: SHIFT_PUBLIC_WORDS,
		n_inout: 0,
		n_witness: witness_len,
		n_internal,
		offset_inout: SHIFT_PUBLIC_WORDS,
		offset_witness: committed_witness_base(),
		committed_total_len,
		n_scratch: 0,
	}
}

/// Return a full power-of-two value vector containing public zeros and the committed Keccak
/// witness at [`committed_witness_base`].
pub fn shift_witness_words(witness: &CommittedKeccakWitness) -> Vec<Word> {
	let layout = value_vec_layout(witness.n_permutations());
	let mut words = vec![Word::ZERO; layout.committed_total_len];
	let start = committed_witness_base();
	words[start..start + witness.words().len()].copy_from_slice(witness.words());
	words
}

/// Return a production `ValueVec` containing the shifted Keccak witness.
pub fn shift_value_vec(witness: &CommittedKeccakWitness) -> ValueVec {
	let layout = value_vec_layout(witness.n_permutations());
	let words = shift_witness_words(witness);
	ValueVec::new_from_data(
		layout,
		words[..committed_witness_base()].to_vec(),
		words[committed_witness_base()..].to_vec(),
	)
	.expect("constructed Keccak Shift value vector has matching layout")
}

/// Build the BitAnd-shaped operand constraints for all chi/iota rows.
///
/// Each active row corresponds to `(permutation, round, lane)` in trace-major order. The operands
/// handed to Shift contain only committed-witness terms:
///
/// ```text
/// a = B[x+1,y]
/// b = B[x+2,y]
/// c = B[x,y] + A_next[x,y]
/// ```
///
/// The chi outer relation uses `P = 1 + B[x+1,y]` and
/// `C = B[x,y] + A_next[x,y] + iota`. Those all-one and iota contributions are transparent
/// corrections to the outer evaluation claims, not committed witness terms.
pub fn chi_constraints_at(base_index: usize, n_permutations: usize) -> Vec<AndConstraint> {
	let mut constraints = Vec::with_capacity(chi_row_count(n_permutations));

	for permutation in 0..n_permutations {
		for round in 0..N_ROUNDS {
			for lane_idx in 0..N_LANES {
				let (p_lane, q_lane, r_lane) = CHI_OPERAND_LANES[lane_idx];
				let mut c = virtual_b_operand_at(base_index, permutation, round, r_lane);
				c.push(ShiftedValueIndex::plain(value_index_at(
					base_index,
					layout::a_index(permutation, round + 1, lane_idx),
				)));

				constraints.push(AndConstraint {
					a: virtual_b_operand_at(base_index, permutation, round, p_lane),
					b: virtual_b_operand_at(base_index, permutation, round, q_lane),
					c,
				});
			}
		}
	}

	constraints
}

/// Build the BitAnd-shaped operand constraints for all chi/iota rows at the canonical witness base.
#[inline]
pub fn chi_constraints(n_permutations: usize) -> Vec<AndConstraint> {
	chi_constraints_at(committed_witness_base(), n_permutations)
}

/// Build degenerate `0 * 0 = D_correctness_operand` constraints.
///
/// For valid witnesses, the third operand evaluates to zero:
///
/// ```text
/// D_r[x] + sum_y A_r[x-1,y] + sum_y rotl_1(A_r[x+1,y]) = 0.
/// ```
pub fn d_correctness_constraints_at(
	base_index: usize,
	n_permutations: usize,
) -> Vec<AndConstraint> {
	let mut constraints = Vec::with_capacity(d_correctness_row_count(n_permutations));

	for permutation in 0..n_permutations {
		for round in 0..N_ROUNDS {
			for x in 0..layout::D_WORDS_PER_BLOCK {
				constraints.push(AndConstraint {
					a: Vec::new(),
					b: Vec::new(),
					c: d_correctness_operand_at(base_index, permutation, round, x),
				});
			}
		}
	}

	constraints
}

/// Build the `D` correctness constraints at the canonical witness base.
#[inline]
pub fn d_correctness_constraints(n_permutations: usize) -> Vec<AndConstraint> {
	d_correctness_constraints_at(committed_witness_base(), n_permutations)
}

/// Build and validate the Shift constraint system for chi operand pushback.
pub fn chi_constraint_system(n_permutations: usize) -> ConstraintSystem {
	let mut constraint_system = ConstraintSystem::new(
		vec![Word::ZERO; SHIFT_PUBLIC_WORDS],
		value_vec_layout(n_permutations),
		chi_constraints(n_permutations),
		Vec::new(),
	);
	constraint_system
		.validate_and_prepare()
		.expect("constructed Keccak chi Shift schema is valid");
	constraint_system
}

/// Build and validate the Shift constraint system for `D` correctness pushback.
pub fn d_correctness_constraint_system(n_permutations: usize) -> ConstraintSystem {
	let mut constraint_system = ConstraintSystem::new(
		vec![Word::ZERO; SHIFT_PUBLIC_WORDS],
		value_vec_layout(n_permutations),
		d_correctness_constraints(n_permutations),
		Vec::new(),
	);
	constraint_system
		.validate_and_prepare()
		.expect("constructed Keccak D-correctness Shift schema is valid");
	constraint_system
}

/// Evaluate the active-row selector for chi rows at the Shift row challenge.
///
/// This is the transparent correction that converts a full chi `P = 1 + B[x+1,y]` claim into the
/// witness-only `B[x+1,y]` operand claim consumed by Shift.
pub fn chi_active_selector_eval<F>(n_permutations: usize, r_x_prime: &[F]) -> F
where
	F: Field,
{
	active_row_selector_eval(chi_row_count(n_permutations), r_x_prime)
}

/// Evaluate the transparent iota column in the chi `C` operand at the Shift challenge.
///
/// This correction converts `C = B[x,y] + A_next[x,y] + iota` into the witness-only
/// `B[x,y] + A_next[x,y]` operand consumed by Shift.
pub fn chi_iota_eval<F>(
	n_permutations: usize,
	r_zhat_prime: F,
	r_x_prime: &[F],
	subspace: &BinarySubspace<F>,
) -> F
where
	F: BinaryField,
{
	assert_eq!(r_x_prime.len(), chi_constraint_system_log_rows(n_permutations));
	let row_weights = eq_ind_partial_eval_scalars(r_x_prime);
	let bit_evals = lagrange_evals(subspace, r_zhat_prime);
	let bit_evals = bit_evals.as_ref();

	let mut result = F::ZERO;
	for permutation in 0..n_permutations {
		for round in 0..N_ROUNDS {
			let row_index = chi_constraint_row_index(permutation, round, 0);
			result +=
				row_weights[row_index] * fold_word_with_lagrange(ROUND_CONSTANTS[round], bit_evals);
		}
	}
	result
}

/// Convert full chi-outer operand evaluations into the witness-only evaluations expected by Shift.
pub fn witness_only_chi_evals<F>(
	n_permutations: usize,
	r_zhat_prime: F,
	r_x_prime: &[F],
	subspace: &BinarySubspace<F>,
	p_eval: F,
	q_eval: F,
	c_eval: F,
) -> [F; 3]
where
	F: BinaryField,
{
	[
		p_eval + chi_active_selector_eval(n_permutations, r_x_prime),
		q_eval,
		c_eval + chi_iota_eval(n_permutations, r_zhat_prime, r_x_prime, subspace),
	]
}

/// Return the log row count of the prepared chi Shift constraint system.
#[inline]
pub fn chi_constraint_system_log_rows(n_permutations: usize) -> usize {
	chi_row_count(n_permutations)
		.max(1)
		.next_power_of_two()
		.ilog2() as usize
}

/// Return the log row count of the prepared `D` correctness Shift constraint system.
#[inline]
pub fn d_correctness_constraint_system_log_rows(n_permutations: usize) -> usize {
	d_correctness_row_count(n_permutations)
		.max(1)
		.next_power_of_two()
		.ilog2() as usize
}

#[inline]
fn chi_constraint_row_index(permutation: usize, round: usize, lane_idx: usize) -> usize {
	(permutation * N_ROUNDS + round) * N_LANES + lane_idx
}

fn active_row_selector_eval<F>(active_rows: usize, r_x_prime: &[F]) -> F
where
	F: Field,
{
	let row_weights = eq_ind_partial_eval_scalars(r_x_prime);
	assert!(active_rows <= row_weights.len());
	row_weights[..active_rows].iter().copied().sum()
}

fn fold_word_with_lagrange<F>(word: u64, lagrange_evals: &[F]) -> F
where
	F: Field,
{
	lagrange_evals
		.iter()
		.enumerate()
		.filter_map(|(bit_idx, &eval)| (((word >> bit_idx) & 1) == 1).then_some(eval))
		.sum()
}

#[cfg(test)]
mod tests {
	use binius_core::{
		constraint_system::{MulConstraint, Operand},
		verify::{eval_operand, eval_shifted_word, verify_constraints},
	};
	use binius_field::{AESTowerField8b, BinaryField128bGhash, PackedBinaryGhash2x128b, Random};
	use binius_math::{
		inner_product::{inner_product, inner_product_buffers},
		multilinear::eq::eq_ind_partial_eval,
		univariate::lagrange_evals,
	};
	use binius_prover::{
		fold_word::fold_words,
		protocols::shift::{
			OperatorData as ProverOperatorData, build_key_collection, prove as prove_shift,
		},
	};
	use binius_transcript::ProverTranscript;
	use binius_verifier::{
		config::{LOG_WORD_SIZE_BITS, StdChallenger},
		protocols::shift::{
			OperatorData as VerifierOperatorData, check_eval, verify as verify_shift,
		},
	};
	use rand::{Rng, SeedableRng, rngs::StdRng};

	use crate::{
		trace::{PermutationTrace, State},
		witness::CommittedKeccakWitness,
	};

	use super::*;

	type F = BinaryField128bGhash;
	type P = PackedBinaryGhash2x128b;

	fn eval_operand_words(words: &[Word], operand: &Operand) -> Word {
		operand.iter().fold(Word::ZERO, |acc, term| {
			acc ^ eval_shifted_word(
				words[term.value_index.0 as usize],
				term.shift_variant,
				term.amount,
			)
		})
	}

	fn compute_bitand_images(
		constraints: &[AndConstraint],
		value_vec: &ValueVec,
	) -> [Vec<Word>; 3] {
		let mut a_image = Vec::with_capacity(constraints.len());
		let mut b_image = Vec::with_capacity(constraints.len());
		let mut c_image = Vec::with_capacity(constraints.len());

		for constraint in constraints {
			a_image.push(eval_operand(value_vec, &constraint.a));
			b_image.push(eval_operand(value_vec, &constraint.b));
			c_image.push(eval_operand(value_vec, &constraint.c));
		}

		[a_image, b_image, c_image]
	}

	fn compute_intmul_images(
		constraints: &[MulConstraint],
		value_vec: &ValueVec,
	) -> [Vec<Word>; 4] {
		let mut a_image = Vec::with_capacity(constraints.len());
		let mut b_image = Vec::with_capacity(constraints.len());
		let mut lo_image = Vec::with_capacity(constraints.len());
		let mut hi_image = Vec::with_capacity(constraints.len());

		for constraint in constraints {
			a_image.push(eval_operand(value_vec, &constraint.a));
			b_image.push(eval_operand(value_vec, &constraint.b));
			lo_image.push(eval_operand(value_vec, &constraint.lo));
			hi_image.push(eval_operand(value_vec, &constraint.hi));
		}

		[a_image, b_image, lo_image, hi_image]
	}

	fn evaluate_image(
		subspace: &BinarySubspace<F>,
		image: &[Word],
		r_zhat_prime: F,
		r_x_prime_tensor: &[F],
	) -> F {
		let l_tilde = lagrange_evals(subspace, r_zhat_prime);
		let l_tilde = l_tilde.as_ref();
		let univariate = image
			.iter()
			.map(|&word| {
				(0..64)
					.filter(|&i| ((word.0 >> i) & 1) == 1)
					.map(|i| l_tilde[i as usize])
					.sum()
			})
			.collect::<Vec<_>>();
		inner_product(r_x_prime_tensor.iter().copied(), univariate.iter().copied())
	}

	fn evaluate_witness(words: &[Word], r_j: &[F], r_y: &[F]) -> F {
		let r_j_tensor = eq_ind_partial_eval::<F>(r_j);
		let r_y_tensor = eq_ind_partial_eval::<F>(r_y);

		let r_j_witness = fold_words::<_, F>(words, r_j_tensor.as_ref());

		inner_product_buffers(&r_j_witness, &r_y_tensor)
	}

	fn prove_and_verify_shift_system(cs: &ConstraintSystem, value_vec: &ValueVec) {
		let mut rng = StdRng::seed_from_u64(31);
		let subspace = BinarySubspace::<AESTowerField8b>::with_dim(LOG_WORD_SIZE_BITS).isomorphic();

		let r_x_prime_bitand = (0..cs.and_constraints.len().ilog2())
			.map(|_| F::random(&mut rng))
			.collect::<Vec<_>>();
		let r_x_prime_intmul = (0..cs.mul_constraints.len().ilog2())
			.map(|_| F::random(&mut rng))
			.collect::<Vec<_>>();
		let r_zhat_prime_bitand = F::random(&mut rng);
		let r_zhat_prime_intmul = F::random(&mut rng);

		let bitand_evals = compute_bitand_images(&cs.and_constraints, value_vec).map(|image| {
			evaluate_image(
				&subspace,
				&image,
				r_zhat_prime_bitand,
				eq_ind_partial_eval(&r_x_prime_bitand).as_ref(),
			)
		});
		let intmul_evals = compute_intmul_images(&cs.mul_constraints, value_vec).map(|image| {
			evaluate_image(
				&subspace,
				&image,
				r_zhat_prime_intmul,
				eq_ind_partial_eval(&r_x_prime_intmul).as_ref(),
			)
		});

		let key_collection = build_key_collection(cs);
		let mut prover_transcript = ProverTranscript::<StdChallenger>::default();
		let prover_output = prove_shift::<F, P, _>(
			&key_collection,
			value_vec.combined_witness(),
			ProverOperatorData {
				evals: bitand_evals.to_vec(),
				r_zhat_prime: r_zhat_prime_bitand,
				r_x_prime: r_x_prime_bitand.clone(),
			},
			ProverOperatorData {
				evals: intmul_evals.to_vec(),
				r_zhat_prime: r_zhat_prime_intmul,
				r_x_prime: r_x_prime_intmul.clone(),
			},
			&mut prover_transcript,
		)
		.unwrap();

		let mut verifier_transcript = prover_transcript.into_verifier();
		let verifier_bitand_data =
			VerifierOperatorData::new(r_zhat_prime_bitand, r_x_prime_bitand, bitand_evals);
		let verifier_intmul_data =
			VerifierOperatorData::new(r_zhat_prime_intmul, r_x_prime_intmul, intmul_evals);
		let verifier_output = verify_shift::<F, _>(
			cs,
			&verifier_bitand_data,
			&verifier_intmul_data,
			&mut verifier_transcript,
		)
		.unwrap();

		check_eval::<F, _>(
			cs,
			&verifier_bitand_data,
			&verifier_intmul_data,
			&subspace,
			&verifier_output,
			&mut verifier_transcript,
		)
		.unwrap();
		verifier_transcript.finalize().unwrap();

		let expected_eval = evaluate_witness(
			value_vec.combined_witness(),
			verifier_output.r_j(),
			verifier_output.r_y(),
		);
		assert_eq!(expected_eval, verifier_output.witness_eval);

		let eval_point = [verifier_output.r_j(), verifier_output.r_y()].concat();
		assert_eq!(prover_output.challenges, eval_point);
		assert_eq!(prover_output.eval, verifier_output.witness_eval);
	}

	#[test]
	fn shift_value_layout_places_committed_witness_after_public_segment() {
		let mut rng = StdRng::seed_from_u64(25);
		let trace = PermutationTrace::new(rng.random::<State>());
		let witness = CommittedKeccakWitness::from_traces(std::slice::from_ref(&trace));
		let shifted = shift_witness_words(&witness);
		let value_vec = shift_value_vec(&witness);

		assert_eq!(committed_witness_base(), 2);
		assert_eq!(shifted.len(), value_vec.size());
		assert_eq!(
			&shifted[committed_witness_base()..committed_witness_base() + witness.words().len()],
			witness.words()
		);
		assert!(
			shifted[..committed_witness_base()]
				.iter()
				.all(|&word| word == Word::ZERO)
		);
	}

	#[test]
	fn chi_constraints_lower_to_virtual_b_and_next_state() {
		let mut rng = StdRng::seed_from_u64(26);
		let trace = PermutationTrace::new(rng.random::<State>());
		let witness = CommittedKeccakWitness::from_traces(std::slice::from_ref(&trace));
		let shifted = shift_witness_words(&witness);
		let constraints = chi_constraints(1);

		assert_eq!(constraints.len(), N_ROUNDS * N_LANES);
		for round in 0..N_ROUNDS {
			for lane_idx in 0..N_LANES {
				let row = chi_constraint_row_index(0, round, lane_idx);
				let (p_lane, q_lane, r_lane) = CHI_OPERAND_LANES[lane_idx];
				let constraint = &constraints[row];

				assert_eq!(
					eval_operand_words(&shifted, &constraint.a),
					Word(trace.rounds[round].pre_chi[p_lane]),
					"round={round} lane={lane_idx} p",
				);
				assert_eq!(
					eval_operand_words(&shifted, &constraint.b),
					Word(trace.rounds[round].pre_chi[q_lane]),
					"round={round} lane={lane_idx} q",
				);
				assert_eq!(
					eval_operand_words(&shifted, &constraint.c),
					Word(
						trace.rounds[round].pre_chi[r_lane] ^ trace.rounds[round].output[lane_idx]
					),
					"round={round} lane={lane_idx} c",
				);
			}
		}
	}

	#[test]
	fn d_correctness_constraints_are_satisfied_as_degenerate_and_rows() {
		let mut rng = StdRng::seed_from_u64(27);
		let trace = PermutationTrace::new(rng.random::<State>());
		let witness = CommittedKeccakWitness::from_traces(std::slice::from_ref(&trace));
		let constraints = d_correctness_constraints(1);
		let value_vec = shift_value_vec(&witness);

		assert_eq!(constraints.len(), N_ROUNDS * layout::D_WORDS_PER_BLOCK);
		for constraint in &constraints {
			assert!(constraint.a.is_empty());
			assert!(constraint.b.is_empty());
			assert_eq!(eval_operand(&value_vec, &constraint.c), Word::ZERO);
		}
	}

	#[test]
	fn d_correctness_constraint_system_validates_and_verifies() {
		let mut rng = StdRng::seed_from_u64(28);
		let traces: Vec<_> = (0..2)
			.map(|_| PermutationTrace::new(rng.random::<State>()))
			.collect();
		let witness = CommittedKeccakWitness::from_traces(&traces);
		let cs = d_correctness_constraint_system(traces.len());
		let value_vec = shift_value_vec(&witness);

		verify_constraints(&cs, &value_vec).unwrap();
		assert_eq!(
			cs.and_constraints.len(),
			d_correctness_row_count(traces.len()).next_power_of_two()
		);
	}

	#[test]
	fn transparent_chi_corrections_match_full_operands() {
		let mut rng = StdRng::seed_from_u64(29);
		let trace = PermutationTrace::new(rng.random::<State>());
		let witness = CommittedKeccakWitness::from_traces(std::slice::from_ref(&trace));
		let value_vec = shift_value_vec(&witness);
		let cs = chi_constraint_system(1);
		let subspace = BinarySubspace::<AESTowerField8b>::with_dim(LOG_WORD_SIZE_BITS).isomorphic();
		let r_x_prime = (0..chi_constraint_system_log_rows(1))
			.map(|_| F::random(&mut rng))
			.collect::<Vec<_>>();
		let r_zhat_prime = F::random(&mut rng);
		let row_weights = eq_ind_partial_eval(&r_x_prime);
		let row_weights = row_weights.as_ref();

		let images = compute_bitand_images(&cs.and_constraints, &value_vec);
		let witness_only_evals = images
			.clone()
			.map(|image| evaluate_image(&subspace, &image, r_zhat_prime, row_weights));

		let mut full_p_image = images[0].clone();
		let full_q_image = images[1].clone();
		let mut full_c_image = images[2].clone();
		for round in 0..N_ROUNDS {
			for lane_idx in 0..N_LANES {
				let row = chi_constraint_row_index(0, round, lane_idx);
				full_p_image[row] = full_p_image[row] ^ Word::ALL_ONE;
				if lane_idx == 0 {
					full_c_image[row] = full_c_image[row] ^ Word(ROUND_CONSTANTS[round]);
				}
			}
		}

		let full_p_eval = evaluate_image(&subspace, &full_p_image, r_zhat_prime, row_weights);
		let full_q_eval = evaluate_image(&subspace, &full_q_image, r_zhat_prime, row_weights);
		let full_c_eval = evaluate_image(&subspace, &full_c_image, r_zhat_prime, row_weights);

		assert_eq!(
			witness_only_chi_evals(
				1,
				r_zhat_prime,
				&r_x_prime,
				&subspace,
				full_p_eval,
				full_q_eval,
				full_c_eval
			),
			witness_only_evals
		);
	}

	#[test]
	fn production_shift_reduces_keccak_chi_operand_claims() {
		let mut rng = StdRng::seed_from_u64(30);
		let traces: Vec<_> = (0..2)
			.map(|_| PermutationTrace::new(rng.random::<State>()))
			.collect();
		let witness = CommittedKeccakWitness::from_traces(&traces);
		let cs = chi_constraint_system(traces.len());
		let value_vec = shift_value_vec(&witness);

		prove_and_verify_shift_system(&cs, &value_vec);
	}

	#[test]
	fn production_shift_reduces_keccak_d_correctness_claims() {
		let mut rng = StdRng::seed_from_u64(32);
		let traces: Vec<_> = (0..2)
			.map(|_| PermutationTrace::new(rng.random::<State>()))
			.collect();
		let witness = CommittedKeccakWitness::from_traces(&traces);
		let cs = d_correctness_constraint_system(traces.len());
		let value_vec = shift_value_vec(&witness);

		prove_and_verify_shift_system(&cs, &value_vec);
	}
}
