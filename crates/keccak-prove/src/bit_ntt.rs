// Copyright 2026 The Binius Developers

//! NTT lookup support for the Keccak chi/iota residual columns.

use binius_field::{BinaryField, BinaryField1b, Field, PackedBinaryField8x1b, PackedField};
use binius_math::{BinarySubspace, univariate::lagrange_evals_scalars};

use crate::{
	constants::{LANE_BITS, N_LANES, ROUND_CONSTANTS},
	trace::State,
};

/// Number of skipped row-index variables in one 64-bit Keccak lane.
pub const LOG_LANE_BITS: usize = 6;

const BYTE_BITS: usize = 8;
const BYTE_VALUES: usize = 1 << BYTE_BITS;
const LANE_BYTES: usize = LANE_BITS / BYTE_BITS;

/// Number of packed field elements needed to hold one 64-point lane evaluation.
pub const PACKED_EVALS: usize = LANE_BITS / 16;

#[rustfmt::skip]
const CHI_OPERAND_LANES: [(usize, usize, usize); N_LANES] = [
	( 1,  2,  0), ( 2,  3,  1), ( 3,  4,  2), ( 4,  0,  3), ( 0,  1,  4),
	( 6,  7,  5), ( 7,  8,  6), ( 8,  9,  7), ( 9,  5,  8), ( 5,  6,  9),
	(11, 12, 10), (12, 13, 11), (13, 14, 12), (14, 10, 13), (10, 11, 14),
	(16, 17, 15), (17, 18, 16), (18, 19, 17), (19, 15, 18), (15, 16, 19),
	(21, 22, 20), (22, 23, 21), (23, 24, 22), (24, 20, 23), (20, 21, 24),
];

/// A byte-chunked NTT lookup for 64 one-bit lane coefficients.
#[derive(Clone)]
pub struct NttLookup<P>(Box<[[[P; PACKED_EVALS]; BYTE_VALUES]; LANE_BYTES]>);

impl<P> NttLookup<P>
where
	P: PackedField,
	P::Scalar: BinaryField + Field,
{
	/// Precompute all byte contributions for a 64-point input domain and output domain.
	pub fn new(
		ntt_input_domain: &BinarySubspace<P::Scalar>,
		ntt_output_domain: &[P::Scalar],
	) -> Self {
		assert_eq!(P::WIDTH, 16);
		assert_eq!(ntt_input_domain.dim(), LOG_LANE_BITS);
		assert_eq!(ntt_output_domain.len(), LANE_BITS);

		let mut lookup = Box::new([[[P::zero(); PACKED_EVALS]; BYTE_VALUES]; LANE_BYTES]);
		let mut eval_point_lagrange_evals =
			vec![vec![P::Scalar::ZERO; LANE_BITS]; ntt_output_domain.len()];

		for (eval_point_idx, eval_point) in ntt_output_domain.iter().enumerate() {
			eval_point_lagrange_evals[eval_point_idx] =
				lagrange_evals_scalars(ntt_input_domain, *eval_point);
		}

		for byte_idx in 0..LANE_BYTES {
			for bit_idx in 0..BYTE_BITS {
				let one_hot_byte = 1 << bit_idx;
				let nonzero_lagrange_basis_coeffs: Vec<_> =
					PackedBinaryField8x1b::from_underlier(one_hot_byte)
						.iter()
						.collect();
				let mut lagrange_basis_coeffs = [BinaryField1b::ZERO; LANE_BITS];

				for (i, nonzero_lagrange_basis_coeff) in
					nonzero_lagrange_basis_coeffs.into_iter().enumerate()
				{
					lagrange_basis_coeffs[byte_idx * BYTE_BITS + i] = nonzero_lagrange_basis_coeff;
				}

				for eval_point_idx in 0..LANE_BITS {
					let mut result = P::Scalar::ZERO;
					for basis_point_idx in 0..LANE_BITS {
						result += eval_point_lagrange_evals[eval_point_idx][basis_point_idx]
							* lagrange_basis_coeffs[basis_point_idx];
					}

					let packed_idx = eval_point_idx / P::WIDTH;
					let scalar_idx = eval_point_idx % P::WIDTH;
					lookup[byte_idx][one_hot_byte as usize][packed_idx].set(scalar_idx, result);
				}
			}
		}

		for byte_idx in 0..LANE_BYTES {
			for byte_value in 0..BYTE_VALUES {
				let mut result = [P::zero(); PACKED_EVALS];
				for bit_idx in 0..BYTE_BITS {
					let one_hot_byte = byte_value & (1 << bit_idx);
					for packed_idx in 0..PACKED_EVALS {
						result[packed_idx] += lookup[byte_idx][one_hot_byte][packed_idx];
					}
				}
				lookup[byte_idx][byte_value] = result;
			}
		}

		Self(lookup)
	}

	/// Precompute the lookup for the upper half of the first prover-message domain.
	pub fn for_upper_half_domain() -> Self
	where
		P::Scalar: From<u8>,
	{
		let (input_domain, output_domain) = upper_half_domains::<P::Scalar>();
		Self::new(&input_domain, &output_domain)
	}

	/// Evaluate a 64-bit word, interpreted as 64 one-bit coefficients.
	#[inline]
	pub fn eval_word(&self, word: u64) -> [P; PACKED_EVALS] {
		let mut result = [P::zero(); PACKED_EVALS];

		for (byte_idx, byte_value) in word.to_le_bytes().into_iter().enumerate() {
			let row = &self.0[byte_idx][byte_value as usize];
			for packed_idx in 0..PACKED_EVALS {
				result[packed_idx] += row[packed_idx];
			}
		}

		result
	}
}

/// Return the 64-point input domain and shifted upper-half output domain.
pub fn upper_half_domains<F>() -> (BinarySubspace<F>, Vec<F>)
where
	F: BinaryField + From<u8>,
{
	let prover_message_domain = BinarySubspace::<F>::with_dim(LOG_LANE_BITS + 1);
	let shift = prover_message_domain.basis()[LOG_LANE_BITS];
	let input_domain = prover_message_domain.reduce_dim(LOG_LANE_BITS);
	let output_domain = input_domain.iter().map(|x| shift + x).collect();

	(input_domain, output_domain)
}

/// Evaluate all Keccak chi/iota residual columns on the shifted upper-half domain.
pub fn upper_half_residual_evals<P>(
	lookup: &NttLookup<P>,
	pre_chi: &State,
	next: &State,
	round: usize,
) -> [[P; PACKED_EVALS]; N_LANES]
where
	P: PackedField,
	P::Scalar: BinaryField + Field,
{
	std::array::from_fn(|lane_idx| {
		let (p_lane, q_lane, r_lane) = CHI_OPERAND_LANES[lane_idx];
		let p = lookup.eval_word(!pre_chi[p_lane]);
		let q = lookup.eval_word(pre_chi[q_lane]);
		let r = lookup.eval_word(pre_chi[r_lane]);
		let next = lookup.eval_word(next[lane_idx]);
		let iota = lookup.eval_word(if lane_idx == 0 {
			ROUND_CONSTANTS[round]
		} else {
			0
		});

		std::array::from_fn(|packed_idx| {
			p[packed_idx] * q[packed_idx] - r[packed_idx] - next[packed_idx] - iota[packed_idx]
		})
	})
}

#[cfg(test)]
mod tests {
	use binius_field::{AESTowerField8b, Field, PackedAESBinaryField16x8b, PackedField};
	use binius_math::{BinarySubspace, univariate::lagrange_evals_scalars};
	use rand::{Rng, SeedableRng, rngs::StdRng};

	use crate::{
		constants::{N_ROUNDS, ROUND_CONSTANTS},
		trace::{PermutationTrace, State},
	};

	use super::*;

	type P = PackedAESBinaryField16x8b;
	type F = AESTowerField8b;

	fn eval_word_direct(input_domain: &BinarySubspace<F>, eval_point: F, word: u64) -> F {
		let lagrange_evals = lagrange_evals_scalars(input_domain, eval_point);
		let mut result = F::ZERO;

		for (bit_idx, lagrange_eval) in lagrange_evals.into_iter().enumerate() {
			if (word >> bit_idx) & 1 == 1 {
				result += lagrange_eval;
			}
		}

		result
	}

	#[test]
	fn ntt_lookup_matches_direct_lagrange_eval() {
		let (input_domain, output_domain) = upper_half_domains::<F>();
		let lookup = NttLookup::<P>::new(&input_domain, &output_domain);
		let mut rng = StdRng::seed_from_u64(3);

		for _ in 0..16 {
			let word = rng.random::<u64>();
			let packed = lookup.eval_word(word);

			for (i, eval_point) in output_domain.iter().copied().enumerate() {
				assert_eq!(
					P::iter_slice(&packed).nth(i).unwrap(),
					eval_word_direct(&input_domain, eval_point, word)
				);
			}
		}
	}

	#[test]
	fn upper_half_residual_evals_match_direct_extension() {
		let (input_domain, output_domain) = upper_half_domains::<F>();
		let lookup = NttLookup::<P>::new(&input_domain, &output_domain);
		let mut rng = StdRng::seed_from_u64(4);

		for _ in 0..4 {
			let trace = PermutationTrace::new(rng.random::<State>());

			for round in 0..N_ROUNDS {
				let round_trace = trace.rounds[round];
				let residuals = upper_half_residual_evals::<P>(
					&lookup,
					&round_trace.pre_chi,
					&round_trace.output,
					round,
				);

				for lane_idx in 0..N_LANES {
					let (p_lane, q_lane, r_lane) = CHI_OPERAND_LANES[lane_idx];
					let iota = if lane_idx == 0 {
						ROUND_CONSTANTS[round]
					} else {
						0
					};

					for (i, eval_point) in output_domain.iter().copied().enumerate() {
						let p = eval_word_direct(
							&input_domain,
							eval_point,
							!round_trace.pre_chi[p_lane],
						);
						let q = eval_word_direct(
							&input_domain,
							eval_point,
							round_trace.pre_chi[q_lane],
						);
						let r = eval_word_direct(
							&input_domain,
							eval_point,
							round_trace.pre_chi[r_lane],
						);
						let next = eval_word_direct(
							&input_domain,
							eval_point,
							round_trace.output[lane_idx],
						);
						let iota = eval_word_direct(&input_domain, eval_point, iota);

						assert_eq!(
							P::iter_slice(&residuals[lane_idx]).nth(i).unwrap(),
							p * q - r - next - iota,
							"round {round}, lane {lane_idx}, eval {i}"
						);
					}
				}
			}
		}
	}
}
