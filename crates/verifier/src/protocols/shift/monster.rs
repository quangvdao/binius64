// Copyright 2025 Irreducible Inc.

use std::{array, iter};

use binius_core::{
	ShiftVariant,
	constraint_system::{Operand, ShiftedValueIndex},
};
use binius_field::{BinaryField, FieldOps, util::powers};
use binius_math::{
	BinarySubspace,
	inner_product::inner_product_scalars,
	multilinear::{
		eq::{eq_ind, eq_ind_partial_eval_scalars},
		evaluate::evaluate_inplace_scalars,
	},
	univariate::lagrange_evals_scalars,
};

use super::{
	SHIFT_VARIANT_COUNT,
	error::Error,
	shift_ind::{partial_eval_phi, partial_eval_sigmas, partial_eval_sigmas_transpose},
	verify::OperatorData,
};
use crate::config::{LOG_WORD_SIZE_BITS, WORD_SIZE_BITS};

/// Evaluates the three h multilinear polynomials (corresponding to SLL, SRL, SRA) at challenge
/// points.
///
/// This is the verifier's version of the h-parts evaluation - instead of building
/// full multilinear polynomials, it directly computes their evaluations.
pub fn evaluate_h_op<E: FieldOps>(l_tilde: &[E], r_j: &[E], r_s: &[E]) -> [E; SHIFT_VARIANT_COUNT] {
	assert_eq!(l_tilde.len(), WORD_SIZE_BITS);
	assert_eq!(r_j.len(), LOG_WORD_SIZE_BITS);
	assert_eq!(r_s.len(), LOG_WORD_SIZE_BITS);

	// Use helper functions to compute shift indicator helpers for 64-bit shifts
	let (sigma, sigma_prime) = partial_eval_sigmas(r_j, r_s);
	let sigma_transpose = partial_eval_sigmas_transpose(r_j, r_s);
	let phi = partial_eval_phi(r_s);
	let j_product: E = r_j.iter().cloned().product();

	// Use helper functions to compute shift indicator helpers for 32-bit shifts
	let (sigma32, sigma32_prime) = partial_eval_sigmas(&r_j[..5], &r_s[..5]);
	let sigma32_transpose = partial_eval_sigmas_transpose(&r_j[..5], &r_s[..5]);
	let phi32 = partial_eval_phi(&r_s[..5]);
	let j_product32: E = r_j[..5].iter().cloned().product();

	// Compute final results
	let sll = inner_product_scalars(l_tilde.iter().cloned(), sigma_transpose);
	let srl = inner_product_scalars(l_tilde.iter().cloned(), sigma.iter().cloned());
	// sra == ∑ᵢ L̃(i) ⋅ (srlᵢ + ∏ₖ rⱼ[k] ⋅ φᵢ)
	//     == ∑ᵢ L̃(i) ⋅ srlᵢ + ∏ₖ rⱼ[k] ⋅ [ ∑ᵢ L̃(i) ⋅ φᵢ ]
	//     == srl + ∏ₖ rⱼ[k] ⋅ [ ∑ᵢ L̃(i) ⋅ φᵢ ]
	let sra = srl.clone() + j_product * inner_product_scalars(l_tilde.iter().cloned(), phi);
	let rotr = inner_product_scalars(
		l_tilde.iter().cloned(),
		iter::zip(&sigma, &sigma_prime).map(|(s_i, s_prime_i)| s_i.clone() + s_prime_i),
	);

	let r_j_rest_tensor = eq_ind_partial_eval_scalars(&r_j[5..]);
	let chunk_size = 1 << 5; // 32

	let sll32 = inner_product_scalars(
		l_tilde.chunks(chunk_size).map(|chunk| {
			inner_product_scalars(chunk.iter().cloned(), sigma32_transpose.iter().cloned())
		}),
		r_j_rest_tensor.iter().cloned(),
	);
	let srl32 = inner_product_scalars(
		l_tilde
			.chunks(chunk_size)
			.map(|chunk| inner_product_scalars(chunk.iter().cloned(), sigma32.iter().cloned())),
		r_j_rest_tensor.iter().cloned(),
	);
	let sra32 = srl32.clone()
		+ inner_product_scalars(
			l_tilde.chunks(chunk_size).map(|chunk| {
				j_product32.clone()
					* inner_product_scalars(chunk.iter().cloned(), phi32.iter().cloned())
			}),
			r_j_rest_tensor.iter().cloned(),
		);
	let rotr32 = inner_product_scalars(
		l_tilde.chunks(chunk_size).map(|chunk| {
			inner_product_scalars(
				chunk.iter().cloned(),
				iter::zip(&sigma32, &sigma32_prime).map(|(s_i, s_prime_i)| s_i.clone() + s_prime_i),
			)
		}),
		r_j_rest_tensor,
	);

	[sll, srl, sra, rotr, sll32, srl32, sra32, rotr32]
}

/// Evaluates the monster multilinear polynomial for a constraint operation.
///
/// The monster multilinear encodes all constraints of a given type (AND or MUL) into
/// a single polynomial. For an operation with multiple operands, it computes:
///
/// $$
/// \sum_{\text{m_idx} \in \text{enumerate(operands)}}
///     \lambda^{\text{m_idx}+1}
///     \sum_{\text{op}} h_{\text{op}}(r_j, r_s) \cdot M_{\text{m}, \text{op}}(r_x', r_y, r_s)
/// $$
///
/// where:
/// - `m_idx` indexes the operand position (0 to ARITY-1)
/// - `op` ranges over the shift variants (logical/arithmetic shifts)
/// - `h_op` is the shift selector polynomial
/// - `M_m,op` is the multilinear extension of the operand values
///
/// # Arguments
///
/// * `operand_vecs` - Vector of operand vectors, one per operand position in the constraint
/// * `operator_data` - Contains the multilinear challenge `r_x'`, univariate challenge `r_zhat'`,
///   and evaluation claims for each operand
/// * `lambda` - Random coefficient for batching operand evaluations
/// * `r_j` - Challenge point for bit index variables (length `LOG_WORD_SIZE_BITS`)
/// * `r_s` - Challenge point for shift variables (length `LOG_WORD_SIZE_BITS`)
/// * `r_y` - Challenge point for word index variables
///
/// # Returns
///
/// The evaluation of the monster multilinear at the given challenge points, or an error
/// if the computation fails.
///
/// # Errors
///
/// Returns an error if the binary subspace construction fails.
pub fn evaluate_monster_multilinear_for_operation<F, E, const ARITY: usize>(
	operand_vecs: &[Vec<&Operand>],
	operator_data: &OperatorData<E, ARITY>,
	subspace: &BinarySubspace<F>,
	lambda: E,
	r_j: &[E],
	r_s: &[E],
	r_y: &[E],
) -> Result<E, Error>
where
	F: BinaryField,
	E: FieldOps<Scalar = F> + From<F>,
{
	assert_eq!(subspace.dim(), LOG_WORD_SIZE_BITS); // precondition

	let r_x_prime_tensor = eq_ind_partial_eval_scalars(&operator_data.r_x_prime);
	let r_y_tensor = eq_ind_partial_eval_scalars(r_y);

	let l_tilde = lagrange_evals_scalars(subspace, operator_data.r_zhat_prime.clone());
	let h_op_evals = evaluate_h_op(&l_tilde, r_j, r_s);

	let lambda_powers = powers(lambda).skip(1).take(ARITY).collect::<Vec<_>>();
	let evals = evaluate_matrices(operand_vecs, &lambda_powers, &r_x_prime_tensor, &r_y_tensor);

	let eval = inner_product_scalars(
		evals.map(|mut evals_op| evaluate_inplace_scalars(&mut evals_op[..], r_s)),
		h_op_evals,
	);

	Ok(eval)
}

/// Tensor layout for a batch of identical operation constraints.
///
/// This describes a flat materialization where both constraint and value indices are laid out as:
///
/// ```text
/// flat_index = instance * base_len + local_index
/// ```
///
/// Equivalently, the low multilinear variables address the local base circuit and the high
/// variables address the repeated instance axis.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct RepeatedMonsterLayout {
	/// Number of instance-axis multilinear variables.
	pub log_instances: usize,
	/// Number of constraints in the base operation after operation-specific padding.
	pub base_constraint_count: usize,
	/// Number of committed values in the base circuit after value-vector padding.
	pub base_value_count: usize,
}

impl RepeatedMonsterLayout {
	/// Creates a repeated monster layout.
	///
	/// # Panics
	///
	/// Panics if either base size is not a power of two.
	pub fn new(
		log_instances: usize,
		base_constraint_count: usize,
		base_value_count: usize,
	) -> Self {
		assert!(
			base_constraint_count.is_power_of_two(),
			"base constraint count must be a power of two"
		);
		assert!(base_value_count.is_power_of_two(), "base value count must be a power of two");
		Self {
			log_instances,
			base_constraint_count,
			base_value_count,
		}
	}

	fn base_log_constraints(self) -> usize {
		self.base_constraint_count.ilog2() as usize
	}

	fn base_log_values(self) -> usize {
		self.base_value_count.ilog2() as usize
	}
}

/// Evaluates the monster multilinear for a repeated identical operation without scanning every
/// materialized instance.
///
/// `base_operand_vecs` must contain local base-circuit value indices. The verifier challenges must
/// be the same challenges that would be used for the flat materialization with:
///
/// ```text
/// flat_constraint = instance * base_constraint_count + local_constraint
/// flat_value      = instance * base_value_count      + local_value
/// ```
///
/// The instance axis collapses via:
///
/// ```text
/// sum_u eq(r_x_instance, u) * eq(r_y_instance, u) = eq(r_x_instance, r_y_instance)
/// ```
pub fn evaluate_repeated_monster_multilinear_for_operation<F, E, const ARITY: usize>(
	base_operand_vecs: &[Vec<&Operand>],
	operator_data: &OperatorData<E, ARITY>,
	subspace: &BinarySubspace<F>,
	lambda: E,
	r_j: &[E],
	r_s: &[E],
	r_y: &[E],
	layout: RepeatedMonsterLayout,
) -> Result<E, Error>
where
	F: BinaryField,
	E: FieldOps<Scalar = F> + From<F>,
{
	assert_eq!(subspace.dim(), LOG_WORD_SIZE_BITS); // precondition
	assert_eq!(base_operand_vecs.len(), ARITY);
	for operands in base_operand_vecs {
		assert_eq!(operands.len(), layout.base_constraint_count);
	}

	let base_log_constraints = layout.base_log_constraints();
	let base_log_values = layout.base_log_values();

	assert_eq!(operator_data.r_x_prime.len(), base_log_constraints + layout.log_instances);
	assert_eq!(r_y.len(), base_log_values + layout.log_instances);

	let (r_x_local, r_x_instance) = operator_data.r_x_prime.split_at(base_log_constraints);
	let (r_y_local, r_y_instance) = r_y.split_at(base_log_values);

	let r_x_local_tensor = eq_ind_partial_eval_scalars(r_x_local);
	let r_y_local_tensor = eq_ind_partial_eval_scalars(r_y_local);
	let instance_factor = eq_ind(r_x_instance, r_y_instance);

	let l_tilde = lagrange_evals_scalars(subspace, operator_data.r_zhat_prime.clone());
	let h_op_evals = evaluate_h_op(&l_tilde, r_j, r_s);

	let lambda_powers = powers(lambda).skip(1).take(ARITY).collect::<Vec<_>>();
	let evals =
		evaluate_matrices(base_operand_vecs, &lambda_powers, &r_x_local_tensor, &r_y_local_tensor);

	let eval = inner_product_scalars(
		evals.map(|mut evals_op| evaluate_inplace_scalars(&mut evals_op[..], r_s)),
		h_op_evals,
	);

	Ok(instance_factor * eval)
}

/// Calculate a batched sum of the M_{\text{op}}(r'_x, r_y, s) matrices.
///
/// Computes one evaluation per shift op, shift amount pair. The matrix evaluations are scaled by
/// the operand coefficients (λ powers) and accumulated.
///
/// # Arguments
///
/// * `operands` - Three-dimensional array of operands where:
///   - dim 1: operation arity
///   - dim 2: constraints
///   - dim 3: shifted indices per term
/// * `operand_coeffs` - Coefficients (λ powers) for batching operand evaluations
/// * `r_x_prime_tensor` - Multilinear challenge tensor for constraint variables
/// * `r_y_tensor` - Challenge tensor for word index variables
fn evaluate_matrices<F: BinaryField, E: FieldOps<Scalar = F> + From<F>>(
	operands: &[Vec<&Operand>],
	operand_coeffs: &[E],
	r_x_prime_tensor: &[E],
	r_y_tensor: &[E],
) -> [[E; WORD_SIZE_BITS]; SHIFT_VARIANT_COUNT] {
	assert_eq!(operands.len(), operand_coeffs.len());

	let zero_evals = array::from_fn(|_| array::from_fn::<E, WORD_SIZE_BITS, _>(|_| E::zero()));

	iter::zip(operand_coeffs, operands)
		.map(|(coeff, constraint_operands)| {
			let mut evals: [[E; WORD_SIZE_BITS]; SHIFT_VARIANT_COUNT] =
				array::from_fn(|_| array::from_fn::<E, WORD_SIZE_BITS, _>(|_| E::zero()));
			for (operand_terms, constraint_eval) in iter::zip(constraint_operands, r_x_prime_tensor)
			{
				for ShiftedValueIndex {
					value_index,
					shift_variant,
					amount,
				} in *operand_terms
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
					evals[shift_id][*amount] +=
						constraint_eval.clone() * &r_y_tensor[value_index.0 as usize];
				}
			}

			// Scale all evaluations by the batching coefficient.
			for evals_op in &mut evals {
				for evals_op_s in &mut *evals_op {
					*evals_op_s *= coeff;
				}
			}

			evals
		})
		.fold(zero_evals, |mut a, b| {
			for (a_op, b_op) in iter::zip(&mut a, b) {
				for (a_op_s, b_op_s) in iter::zip(&mut *a_op, b_op) {
					*a_op_s += b_op_s;
				}
			}
			a
		})
}

#[cfg(test)]
mod tests {
	use std::{
		array,
		hint::black_box,
		time::{Duration, Instant},
	};

	use crate::protocols::shift::{BITAND_ARITY, INTMUL_ARITY};
	use binius_core::constraint_system::{AndConstraint, MulConstraint, Operand, ValueIndex};
	use binius_field::{BinaryField128bGhash, Field, Random};
	use binius_math::{
		BinarySubspace,
		test_utils::{index_to_hypercube_point, random_scalars},
		univariate::lagrange_evals_scalars,
	};
	use itertools::Itertools;
	use rand::{Rng, SeedableRng, rngs::StdRng};

	use super::*;

	fn average_elapsed(iterations: usize, mut f: impl FnMut()) -> Duration {
		let start = Instant::now();
		for _ in 0..iterations {
			f();
		}
		start.elapsed() / iterations as u32
	}

	fn test_term(value_index: usize, seed: usize) -> ShiftedValueIndex {
		let value_index = ValueIndex(value_index as u32);
		let amount = (seed * 7 + 3) % WORD_SIZE_BITS;
		match seed % SHIFT_VARIANT_COUNT {
			0 => ShiftedValueIndex::plain(value_index),
			1 => ShiftedValueIndex::sll(value_index, amount),
			2 => ShiftedValueIndex::srl(value_index, amount),
			3 => ShiftedValueIndex::sar(value_index, amount),
			4 => ShiftedValueIndex::rotr(value_index, amount),
			5 => ShiftedValueIndex::sll32(value_index, amount % 32),
			6 => ShiftedValueIndex::srl32(value_index, amount % 32),
			_ => ShiftedValueIndex::sra32(value_index, amount % 32),
		}
	}

	fn offset_operand(operand: &Operand, offset: usize) -> Operand {
		operand
			.iter()
			.map(|term| ShiftedValueIndex {
				value_index: ValueIndex(term.value_index.0 + offset as u32),
				shift_variant: term.shift_variant,
				amount: term.amount,
			})
			.collect()
	}

	fn make_base_and_constraints(
		base_constraint_count: usize,
		base_value_count: usize,
	) -> Vec<AndConstraint> {
		(0..base_constraint_count)
			.map(|row| AndConstraint {
				a: vec![
					test_term((row * 3 + 1) % base_value_count, row),
					test_term((row * 5 + 7) % base_value_count, row + 1),
				],
				b: vec![
					test_term((row * 11 + 13) % base_value_count, row + 2),
					test_term((row * 17 + 19) % base_value_count, row + 3),
				],
				c: vec![
					test_term((row * 23 + 29) % base_value_count, row + 4),
					test_term((row * 31 + 37) % base_value_count, row + 5),
				],
			})
			.collect()
	}

	fn repeat_and_constraints(
		base_constraints: &[AndConstraint],
		log_instances: usize,
		base_value_count: usize,
	) -> Vec<AndConstraint> {
		let instances = 1 << log_instances;
		let mut constraints = Vec::with_capacity(instances * base_constraints.len());
		for instance in 0..instances {
			let value_offset = instance * base_value_count;
			for constraint in base_constraints {
				constraints.push(AndConstraint {
					a: offset_operand(&constraint.a, value_offset),
					b: offset_operand(&constraint.b, value_offset),
					c: offset_operand(&constraint.c, value_offset),
				});
			}
		}
		constraints
	}

	fn make_base_mul_constraints(
		base_constraint_count: usize,
		base_value_count: usize,
	) -> Vec<MulConstraint> {
		(0..base_constraint_count)
			.map(|row| MulConstraint {
				a: vec![
					test_term((row * 3 + 2) % base_value_count, row),
					test_term((row * 5 + 11) % base_value_count, row + 1),
				],
				b: vec![
					test_term((row * 7 + 17) % base_value_count, row + 2),
					test_term((row * 13 + 23) % base_value_count, row + 3),
				],
				lo: vec![
					test_term((row * 19 + 31) % base_value_count, row + 4),
					test_term((row * 29 + 41) % base_value_count, row + 5),
				],
				hi: vec![
					test_term((row * 37 + 43) % base_value_count, row + 6),
					test_term((row * 47 + 53) % base_value_count, row + 7),
				],
			})
			.collect()
	}

	fn repeat_mul_constraints(
		base_constraints: &[MulConstraint],
		log_instances: usize,
		base_value_count: usize,
	) -> Vec<MulConstraint> {
		let instances = 1 << log_instances;
		let mut constraints = Vec::with_capacity(instances * base_constraints.len());
		for instance in 0..instances {
			let value_offset = instance * base_value_count;
			for constraint in base_constraints {
				constraints.push(MulConstraint {
					a: offset_operand(&constraint.a, value_offset),
					b: offset_operand(&constraint.b, value_offset),
					lo: offset_operand(&constraint.lo, value_offset),
					hi: offset_operand(&constraint.hi, value_offset),
				});
			}
		}
		constraints
	}

	#[test]
	fn repeated_monster_eval_matches_flat_large_and_materialization() {
		type F = BinaryField128bGhash;

		let mut rng = StdRng::seed_from_u64(1);
		let base_constraint_count = 1 << 8;
		let base_value_count = 1 << 9;
		let log_instances = 8;

		let base_constraints = make_base_and_constraints(base_constraint_count, base_value_count);
		let flat_constraints =
			repeat_and_constraints(&base_constraints, log_instances, base_value_count);
		assert_eq!(flat_constraints.len(), 1 << 16);

		let operator_data: OperatorData<F, BITAND_ARITY> = OperatorData::new(
			F::random(&mut rng),
			random_scalars(&mut rng, base_constraint_count.ilog2() as usize + log_instances),
			array::from_fn(|_| F::random(&mut rng)),
		);
		let subspace = BinarySubspace::<F>::with_dim(LOG_WORD_SIZE_BITS);
		let lambda = F::random(&mut rng);
		let r_j = random_scalars(&mut rng, LOG_WORD_SIZE_BITS);
		let r_s = random_scalars(&mut rng, LOG_WORD_SIZE_BITS);
		let r_y = random_scalars(&mut rng, base_value_count.ilog2() as usize + log_instances);

		let (base_a, base_b, base_c) = base_constraints
			.iter()
			.map(|AndConstraint { a, b, c }| (a, b, c))
			.multiunzip();
		let (flat_a, flat_b, flat_c) = flat_constraints
			.iter()
			.map(|AndConstraint { a, b, c }| (a, b, c))
			.multiunzip();

		let flat = evaluate_monster_multilinear_for_operation(
			&[flat_a, flat_b, flat_c],
			&operator_data,
			&subspace,
			lambda,
			&r_j,
			&r_s,
			&r_y,
		)
		.unwrap();
		let repeated = evaluate_repeated_monster_multilinear_for_operation(
			&[base_a, base_b, base_c],
			&operator_data,
			&subspace,
			lambda,
			&r_j,
			&r_s,
			&r_y,
			RepeatedMonsterLayout::new(log_instances, base_constraint_count, base_value_count),
		)
		.unwrap();

		assert_eq!(repeated, flat);
	}

	#[test]
	fn repeated_monster_eval_matches_flat_mul_materialization() {
		type F = BinaryField128bGhash;

		let mut rng = StdRng::seed_from_u64(2);
		let base_constraint_count = 1 << 7;
		let base_value_count = 1 << 8;
		let log_instances = 7;

		let base_constraints = make_base_mul_constraints(base_constraint_count, base_value_count);
		let flat_constraints =
			repeat_mul_constraints(&base_constraints, log_instances, base_value_count);
		assert_eq!(flat_constraints.len(), 1 << 14);

		let operator_data: OperatorData<F, INTMUL_ARITY> = OperatorData::new(
			F::random(&mut rng),
			random_scalars(&mut rng, base_constraint_count.ilog2() as usize + log_instances),
			array::from_fn(|_| F::random(&mut rng)),
		);
		let subspace = BinarySubspace::<F>::with_dim(LOG_WORD_SIZE_BITS);
		let lambda = F::random(&mut rng);
		let r_j = random_scalars(&mut rng, LOG_WORD_SIZE_BITS);
		let r_s = random_scalars(&mut rng, LOG_WORD_SIZE_BITS);
		let r_y = random_scalars(&mut rng, base_value_count.ilog2() as usize + log_instances);

		let (base_a, base_b, base_lo, base_hi) = base_constraints
			.iter()
			.map(|MulConstraint { a, b, lo, hi }| (a, b, lo, hi))
			.multiunzip();
		let (flat_a, flat_b, flat_lo, flat_hi) = flat_constraints
			.iter()
			.map(|MulConstraint { a, b, lo, hi }| (a, b, lo, hi))
			.multiunzip();

		let flat = evaluate_monster_multilinear_for_operation(
			&[flat_a, flat_b, flat_lo, flat_hi],
			&operator_data,
			&subspace,
			lambda,
			&r_j,
			&r_s,
			&r_y,
		)
		.unwrap();
		let repeated = evaluate_repeated_monster_multilinear_for_operation(
			&[base_a, base_b, base_lo, base_hi],
			&operator_data,
			&subspace,
			lambda,
			&r_j,
			&r_s,
			&r_y,
			RepeatedMonsterLayout::new(log_instances, base_constraint_count, base_value_count),
		)
		.unwrap();

		assert_eq!(repeated, flat);
	}

	#[test]
	#[ignore = "prints concrete verifier-side runtimes for flat vs repeated monster evaluation"]
	fn repeated_monster_eval_print_runtimes() {
		type F = BinaryField128bGhash;

		let base_constraint_count = 1 << 8;
		let base_value_count = 1 << 9;
		let base_constraints = make_base_and_constraints(base_constraint_count, base_value_count);
		let (base_a, base_b, base_c): (Vec<_>, Vec<_>, Vec<_>) = base_constraints
			.iter()
			.map(|AndConstraint { a, b, c }| (a, b, c))
			.multiunzip();
		let base_operands = [base_a, base_b, base_c];
		let subspace = BinarySubspace::<F>::with_dim(LOG_WORD_SIZE_BITS);

		println!("base_constraints={base_constraint_count}, base_values={base_value_count}");
		println!("log_instances,instances,flat_constraints,flat_ms,repeated_ms,speedup");

		for log_instances in [0usize, 4, 8, 10, 12] {
			let mut rng = StdRng::seed_from_u64(10_000 + log_instances as u64);
			let flat_constraints =
				repeat_and_constraints(&base_constraints, log_instances, base_value_count);
			let (flat_a, flat_b, flat_c): (Vec<_>, Vec<_>, Vec<_>) = flat_constraints
				.iter()
				.map(|AndConstraint { a, b, c }| (a, b, c))
				.multiunzip();
			let flat_operands = [flat_a, flat_b, flat_c];

			let operator_data: OperatorData<F, BITAND_ARITY> = OperatorData::new(
				F::random(&mut rng),
				random_scalars(&mut rng, base_constraint_count.ilog2() as usize + log_instances),
				array::from_fn(|_| F::random(&mut rng)),
			);
			let lambda = F::random(&mut rng);
			let r_j = random_scalars(&mut rng, LOG_WORD_SIZE_BITS);
			let r_s = random_scalars(&mut rng, LOG_WORD_SIZE_BITS);
			let r_y = random_scalars(&mut rng, base_value_count.ilog2() as usize + log_instances);

			let flat = evaluate_monster_multilinear_for_operation(
				&flat_operands,
				&operator_data,
				&subspace,
				lambda,
				&r_j,
				&r_s,
				&r_y,
			)
			.unwrap();

			let repeated_layout =
				RepeatedMonsterLayout::new(log_instances, base_constraint_count, base_value_count);
			let repeated = evaluate_repeated_monster_multilinear_for_operation(
				&base_operands,
				&operator_data,
				&subspace,
				lambda,
				&r_j,
				&r_s,
				&r_y,
				repeated_layout,
			)
			.unwrap();
			assert_eq!(repeated, flat);

			let iterations = match log_instances {
				0 | 4 => 256,
				8 => 64,
				10 => 32,
				_ => 8,
			};
			let flat_elapsed = average_elapsed(iterations, || {
				let eval = evaluate_monster_multilinear_for_operation(
					&flat_operands,
					&operator_data,
					&subspace,
					lambda,
					&r_j,
					&r_s,
					&r_y,
				)
				.unwrap();
				black_box(eval);
			});
			let repeated_elapsed = average_elapsed(iterations, || {
				let eval = evaluate_repeated_monster_multilinear_for_operation(
					&base_operands,
					&operator_data,
					&subspace,
					lambda,
					&r_j,
					&r_s,
					&r_y,
					repeated_layout,
				)
				.unwrap();
				black_box(eval);
			});

			let flat_again = evaluate_monster_multilinear_for_operation(
				&flat_operands,
				&operator_data,
				&subspace,
				lambda,
				&r_j,
				&r_s,
				&r_y,
			)
			.unwrap();
			let repeated_again = evaluate_repeated_monster_multilinear_for_operation(
				&base_operands,
				&operator_data,
				&subspace,
				lambda,
				&r_j,
				&r_s,
				&r_y,
				repeated_layout,
			)
			.unwrap();
			assert_eq!(flat_again, repeated_again);

			let flat_ms = flat_elapsed.as_secs_f64() * 1_000.0;
			let repeated_ms = repeated_elapsed.as_secs_f64() * 1_000.0;
			println!(
				"{log_instances},{},{},{flat_ms:.3},{repeated_ms:.3},{:.1}x",
				1usize << log_instances,
				flat_constraints.len(),
				flat_ms / repeated_ms,
			);
		}
	}

	#[test]
	fn test_evaluate_h_op_hypercube_vertices() {
		// Property-based test: for random i, j, s in {0..63}, with challenge being
		// the i-th element of the subspace, the outputs must match indicator relations
		// over integers:
		// - sll == 1 iff j + s == i
		// - srl == 1 iff i + s == j
		// - sra == 1 iff i + s == j || i + s >= 64 && j == 63
		// - rotr == 1 iff (i + s) % 64 == j
		let mut rng = StdRng::seed_from_u64(0);
		let subspace = BinarySubspace::<BinaryField128bGhash>::with_dim(LOG_WORD_SIZE_BITS);

		// Run a reasonable number of random trials
		for _trial in 0..1024 {
			let i = rng.random_range(0..64);
			let j = rng.random_range(0..64);
			let s = rng.random_range(0..64);

			let challenge = subspace.get(i);
			let l_tilde = lagrange_evals_scalars(&subspace, challenge);

			let r_j = index_to_hypercube_point::<BinaryField128bGhash>(LOG_WORD_SIZE_BITS, j);
			let r_s = index_to_hypercube_point::<BinaryField128bGhash>(LOG_WORD_SIZE_BITS, s);

			let [sll, srl, sra, rotr, sll32, srl32, sra32, rotr32] =
				evaluate_h_op(&l_tilde, &r_j, &r_s);

			let expected_sll = j + s == i;
			let expected_srl = i + s == j;
			let expected_sra = (i + s).min(63) == j;
			let expected_rotr = (i + s) % 64 == j;

			let i_hi = i / 32;
			let i_lo = i % 32;
			let j_hi = j / 32;
			let j_lo = j % 32;
			let s_lo = s % 32;

			let expected_sll32 = i_hi == j_hi && j_lo + s_lo == i_lo;
			let expected_srl32 = i_hi == j_hi && i_lo + s_lo == j_lo;
			let expected_sra32 = i_hi == j_hi && (i_lo + s_lo).min(31) == j_lo;
			let expected_rotr32 = i_hi == j_hi && (i_lo + s_lo) % 32 == j_lo;

			let to_field = |b: bool| {
				if b {
					BinaryField128bGhash::ONE
				} else {
					BinaryField128bGhash::ZERO
				}
			};

			assert_eq!(sll, to_field(expected_sll), "sll failed for i={i}, j={j}, s={s}");
			assert_eq!(srl, to_field(expected_srl), "srl failed for i={i}, j={j}, s={s}");
			assert_eq!(sra, to_field(expected_sra), "sra failed for i={i}, j={j}, s={s}");
			assert_eq!(rotr, to_field(expected_rotr), "rotr failed for i={i}, j={j}, s={s}");
			assert_eq!(sll32, to_field(expected_sll32), "sll32 failed for i={i}, j={j}, s={s}");
			assert_eq!(srl32, to_field(expected_srl32), "srl32 failed for i={i}, j={j}, s={s}");
			assert_eq!(sra32, to_field(expected_sra32), "sra32 failed for i={i}, j={j}, s={s}");
			assert_eq!(rotr32, to_field(expected_rotr32), "rotr32 failed for i={i}, j={j}, s={s}");
		}
	}

	#[test]
	fn test_evaluate_h_op_multilinearity() {
		// Test that the function is multilinear in each variable
		let mut rng = StdRng::seed_from_u64(0);

		// Generate random evaluation points
		let challenge = BinaryField128bGhash::random(&mut rng);
		let subspace = BinarySubspace::<BinaryField128bGhash>::with_dim(LOG_WORD_SIZE_BITS);
		let l_tilde = lagrange_evals_scalars(&subspace, challenge);
		let r_j = random_scalars::<BinaryField128bGhash>(&mut rng, LOG_WORD_SIZE_BITS);
		let r_s = random_scalars::<BinaryField128bGhash>(&mut rng, LOG_WORD_SIZE_BITS);

		// Check linearity in each variable
		for i in 0..LOG_WORD_SIZE_BITS {
			// Check r_j[i]
			let mut r_j_at_0 = r_j.clone();
			r_j_at_0[i] = BinaryField128bGhash::ZERO;
			let mut r_j_at_1 = r_j.clone();
			r_j_at_1[i] = BinaryField128bGhash::ONE;
			let [result_0, result_1, result_y] = [&r_j_at_0, &r_j_at_1, &r_j]
				.map(|r_j_variant| evaluate_h_op(&l_tilde, r_j_variant, &r_s));
			for variant in 0..SHIFT_VARIANT_COUNT {
				let expected = result_0[variant] * (BinaryField128bGhash::ONE - r_j[i])
					+ result_1[variant] * r_j[i];
				assert_eq!(result_y[variant], expected, "Not linear in r_j[{i}]");
			}

			// Check r_s[i]
			let mut r_s_at_0 = r_s.clone();
			r_s_at_0[i] = BinaryField128bGhash::ZERO;
			let mut r_s_at_1 = r_s.clone();
			r_s_at_1[i] = BinaryField128bGhash::ONE;
			let [result_0, result_1, result_y] = [&r_s_at_0, &r_s_at_1, &r_s]
				.map(|r_s_variant| evaluate_h_op(&l_tilde, &r_j, r_s_variant));
			for variant in 0..SHIFT_VARIANT_COUNT {
				let expected = result_0[variant] * (BinaryField128bGhash::ONE - r_s[i])
					+ result_1[variant] * r_s[i];
				assert_eq!(result_y[variant], expected, "Not linear in r_s[{i}]");
			}
		}
	}
}
