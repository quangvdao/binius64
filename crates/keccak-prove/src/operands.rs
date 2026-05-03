// Copyright 2026 The Binius Developers

//! BitAnd-style Keccak chi operands and residuals.

use crate::{
	constants::{ROUND_CONSTANTS, lane},
	trace::State,
	unrolled,
};

#[cfg(test)]
use crate::constants::N_LANES;

/// Virtual chi operands for one valid lane position.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct ChiOperands {
	/// `P = 1 + D[x+1,y]`, represented as bitwise NOT over a 64-bit lane.
	pub p: u64,
	/// `Q = D[x+2,y]`.
	pub q: u64,
	/// `R = D[x,y]`.
	pub r: u64,
}

/// Return the chi operand words for valid lane coordinates `(x, y)`.
#[inline]
pub fn chi_operands(pre_chi: &State, x: usize, y: usize) -> ChiOperands {
	debug_assert!(x < 5);
	debug_assert!(y < 5);

	ChiOperands {
		p: !pre_chi[lane((x + 1) % 5, y)],
		q: pre_chi[lane((x + 2) % 5, y)],
		r: pre_chi[lane(x, y)],
	}
}

/// Return the BitAnd-style residual word
///
/// ```text
/// H = P * Q - R - A_next - Iota
/// ```
///
/// over the base bit domain. In characteristic two this is
/// `(p & q) ^ r ^ next ^ iota`.
#[inline]
pub fn chi_iota_residual_word(
	pre_chi: &State,
	next: &State,
	round: usize,
	x: usize,
	y: usize,
) -> u64 {
	let ChiOperands { p, q, r } = chi_operands(pre_chi, x, y);
	let iota = if x == 0 && y == 0 {
		ROUND_CONSTANTS[round]
	} else {
		0
	};

	(p & q) ^ r ^ next[lane(x, y)] ^ iota
}

/// Compute all 25 residual words for one Keccak round.
pub fn chi_iota_residual_words(pre_chi: &State, next: &State, round: usize) -> State {
	unrolled::chi_iota_residual_words(pre_chi, next, round)
}

#[cfg(test)]
fn chi_iota_residual_words_reference(pre_chi: &State, next: &State, round: usize) -> State {
	let mut residual = [0u64; N_LANES];
	for y in 0..5 {
		for x in 0..5 {
			residual[lane(x, y)] = chi_iota_residual_word(pre_chi, next, round, x, y);
		}
	}
	residual
}

#[cfg(test)]
mod tests {
	use rand::{Rng, SeedableRng, rngs::StdRng};

	use crate::{
		constants::{N_ROUNDS, lane},
		trace::{PermutationTrace, State},
	};

	use super::*;

	#[test]
	fn keccak_bitand_residual_matches_native_round() {
		let mut rng = StdRng::seed_from_u64(1);

		for _ in 0..8 {
			let trace = PermutationTrace::new(rng.random::<State>());

			for round in 0..N_ROUNDS {
				let round_trace = trace.rounds[round];
				let residual =
					chi_iota_residual_words(&round_trace.pre_chi, &round_trace.output, round);

				assert_eq!(residual, [0u64; N_LANES], "round {round}");
				assert_eq!(
					residual,
					chi_iota_residual_words_reference(
						&round_trace.pre_chi,
						&round_trace.output,
						round
					)
				);
			}
		}
	}

	#[test]
	fn chi_operands_match_lane_neighbors() {
		let pre_chi: State = std::array::from_fn(|i| 1u64 << i);

		for y in 0..5 {
			for x in 0..5 {
				let operands = chi_operands(&pre_chi, x, y);
				assert_eq!(operands.p, !pre_chi[lane((x + 1) % 5, y)]);
				assert_eq!(operands.q, pre_chi[lane((x + 2) % 5, y)]);
				assert_eq!(operands.r, pre_chi[lane(x, y)]);
			}
		}
	}
}
