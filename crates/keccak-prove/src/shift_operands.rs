// Copyright 2026 The Binius Developers

//! Shift-compatible operands for the committed Keccak witness layout.

use binius_core::constraint_system::{Operand, ShiftedValueIndex, ValueIndex};

use crate::{
	constants::{N_LANES, RHO_OFFSETS, lane_x},
	layout,
};

/// Source lane of each virtual pre-chi `B` lane after rho+pi.
#[rustfmt::skip]
pub const RHO_PI_PREIMAGE: [usize; N_LANES] = [
	 0,  6, 12, 18, 24,
	 3,  9, 10, 16, 22,
	 1,  7, 13, 19, 20,
	 4,  5, 11, 17, 23,
	 2,  8, 14, 15, 21,
];

/// Convert a committed witness word index to a production `ValueIndex`.
#[inline]
pub fn value_index(index: usize) -> ValueIndex {
	assert!(u32::try_from(index).is_ok());
	ValueIndex(index as u32)
}

/// Return a shifted term for rotating a committed word left by `amount`.
#[inline]
pub fn rotl_term(index: usize, amount: u32) -> ShiftedValueIndex {
	let amount = (amount % 64) as usize;
	let index = value_index(index);
	if amount == 0 {
		ShiftedValueIndex::plain(index)
	} else {
		ShiftedValueIndex::rotr(index, 64 - amount)
	}
}

/// Return the virtual `B_r[b_lane]` operand lowered to committed `A` and `D` terms.
///
/// Keccak rho+pi gives:
///
/// ```text
/// B_r[pi(x,y)] = rotl_{rho[x,y]}(A_r[x,y] + D_r[x])
/// ```
///
/// Since `D` is committed and `B` is not, this lowers every virtual `B` reference to two shifted
/// committed terms.
pub fn virtual_b_operand(permutation: usize, round: usize, b_lane: usize) -> Operand {
	assert!(b_lane < N_LANES);
	let source_lane = RHO_PI_PREIMAGE[b_lane];
	let source_x = lane_x(source_lane);
	let rho = RHO_OFFSETS[source_lane];

	vec![
		rotl_term(layout::a_index(permutation, round, source_lane), rho),
		rotl_term(layout::d_index(permutation, round, source_x), rho),
	]
}

/// Return the linear operand proving `D_r[x]` is the theta correction for column `x`.
///
/// The relation is:
///
/// ```text
/// D_r[x] + sum_y A_r[x-1,y] + sum_y rotl_1(A_r[x+1,y]) = 0
/// ```
pub fn d_correctness_operand(permutation: usize, round: usize, x: usize) -> Operand {
	assert!(x < layout::D_WORDS_PER_BLOCK);
	let left_x = (x + 4) % 5;
	let right_x = (x + 1) % 5;

	let mut operand = Vec::with_capacity(11);
	operand.push(ShiftedValueIndex::plain(value_index(layout::d_index(permutation, round, x))));
	for y in 0..5 {
		operand.push(ShiftedValueIndex::plain(value_index(layout::a_index(
			permutation,
			round,
			left_x + 5 * y,
		))));
	}
	for y in 0..5 {
		operand.push(rotl_term(layout::a_index(permutation, round, right_x + 5 * y), 1));
	}
	operand
}

#[cfg(test)]
mod tests {
	use binius_core::{verify::eval_shifted_word, word::Word};
	use rand::{Rng, SeedableRng, rngs::StdRng};

	use crate::{
		constants::lane,
		trace::{PermutationTrace, State},
		witness::{CommittedKeccakWitness, theta_corrections},
	};

	use super::*;

	fn eval_operand(words: &[Word], operand: &Operand) -> Word {
		operand.iter().fold(Word::ZERO, |acc, term| {
			acc ^ eval_shifted_word(
				words[term.value_index.0 as usize],
				term.shift_variant,
				term.amount,
			)
		})
	}

	#[test]
	fn rho_pi_preimage_matches_keccak_formula() {
		for b_lane in 0..N_LANES {
			let source_lane = RHO_PI_PREIMAGE[b_lane];
			let x = source_lane % 5;
			let y = source_lane / 5;
			assert_eq!(b_lane, lane(y, (2 * x + 3 * y) % 5));
		}
	}

	#[test]
	fn virtual_b_operand_evaluates_pre_chi_lane() {
		let mut rng = StdRng::seed_from_u64(23);
		let trace = PermutationTrace::new(rng.random::<State>());
		let witness = CommittedKeccakWitness::from_traces(std::slice::from_ref(&trace));

		for round in 0..crate::constants::N_ROUNDS {
			for b_lane in 0..N_LANES {
				assert_eq!(
					eval_operand(witness.words(), &virtual_b_operand(0, round, b_lane)),
					Word(trace.rounds[round].pre_chi[b_lane]),
					"round={round} b_lane={b_lane}",
				);
			}
		}
	}

	#[test]
	fn d_correctness_operand_evaluates_to_zero() {
		let mut rng = StdRng::seed_from_u64(24);
		let trace = PermutationTrace::new(rng.random::<State>());
		let witness = CommittedKeccakWitness::from_traces(std::slice::from_ref(&trace));

		for round in 0..crate::constants::N_ROUNDS {
			assert_eq!(
				theta_corrections(&trace.rounds[round].input),
				std::array::from_fn(|x| witness.d(0, round, x).0)
			);
			for x in 0..layout::D_WORDS_PER_BLOCK {
				assert_eq!(
					eval_operand(witness.words(), &d_correctness_operand(0, round, x)),
					Word::ZERO,
					"round={round} x={x}",
				);
			}
		}
	}
}
