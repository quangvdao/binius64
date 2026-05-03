// Copyright 2026 The Binius Developers

//! Committed `A` and `D` witness construction for Keccak-f[1600].

use binius_core::word::Word;

use crate::{
	constants::{N_LANES, N_ROUNDS, lane},
	layout,
	trace::{PermutationTrace, State},
};

/// Committed witness words for a batch of Keccak-f[1600] permutations.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct CommittedKeccakWitness {
	words: Vec<Word>,
	n_permutations: usize,
}

impl CommittedKeccakWitness {
	/// Build the committed `A` and `D` witness for a batch of permutation traces.
	pub fn from_traces(traces: &[PermutationTrace]) -> Self {
		let mut words = vec![Word::ZERO; layout::witness_len(traces.len())];

		for (permutation, trace) in traces.iter().enumerate() {
			for round in 0..N_ROUNDS {
				fill_state(&mut words, permutation, round, &trace.rounds[round].input);
				let d = theta_corrections(&trace.rounds[round].input);
				for (x, &d_x) in d.iter().enumerate() {
					words[layout::d_index(permutation, round, x)] = Word(d_x);
				}
			}

			fill_state(&mut words, permutation, N_ROUNDS, &trace.output());
		}

		Self {
			words,
			n_permutations: traces.len(),
		}
	}

	/// Return the committed witness words.
	pub fn words(&self) -> &[Word] {
		&self.words
	}

	/// Consume the witness and return the committed words.
	pub fn into_words(self) -> Vec<Word> {
		self.words
	}

	/// Return the number of Keccak-f permutations in the batch.
	pub fn n_permutations(&self) -> usize {
		self.n_permutations
	}

	/// Return one committed `A_r[lane]` word.
	pub fn a(&self, permutation: usize, round_block: usize, lane: usize) -> Word {
		self.words[layout::a_index(permutation, round_block, lane)]
	}

	/// Return one committed `D_r[x]` word.
	pub fn d(&self, permutation: usize, round: usize, x: usize) -> Word {
		self.words[layout::d_index(permutation, round, x)]
	}
}

/// Compute the five Keccak theta correction words for a round input state.
pub fn theta_corrections(state: &State) -> [u64; 5] {
	let c0 = state[lane(0, 0)]
		^ state[lane(0, 1)]
		^ state[lane(0, 2)]
		^ state[lane(0, 3)]
		^ state[lane(0, 4)];
	let c1 = state[lane(1, 0)]
		^ state[lane(1, 1)]
		^ state[lane(1, 2)]
		^ state[lane(1, 3)]
		^ state[lane(1, 4)];
	let c2 = state[lane(2, 0)]
		^ state[lane(2, 1)]
		^ state[lane(2, 2)]
		^ state[lane(2, 3)]
		^ state[lane(2, 4)];
	let c3 = state[lane(3, 0)]
		^ state[lane(3, 1)]
		^ state[lane(3, 2)]
		^ state[lane(3, 3)]
		^ state[lane(3, 4)];
	let c4 = state[lane(4, 0)]
		^ state[lane(4, 1)]
		^ state[lane(4, 2)]
		^ state[lane(4, 3)]
		^ state[lane(4, 4)];

	[
		c4 ^ c1.rotate_left(1),
		c0 ^ c2.rotate_left(1),
		c1 ^ c3.rotate_left(1),
		c2 ^ c4.rotate_left(1),
		c3 ^ c0.rotate_left(1),
	]
}

fn fill_state(words: &mut [Word], permutation: usize, round_block: usize, state: &State) {
	for (lane, &word) in state.iter().enumerate().take(N_LANES) {
		words[layout::a_index(permutation, round_block, lane)] = Word(word);
	}
}

#[cfg(test)]
mod tests {
	use rand::{Rng, SeedableRng, rngs::StdRng};

	use super::*;

	fn theta_corrections_reference(state: &State) -> [u64; 5] {
		let c = std::array::from_fn::<_, 5, _>(|x| {
			state[lane(x, 0)]
				^ state[lane(x, 1)]
				^ state[lane(x, 2)]
				^ state[lane(x, 3)]
				^ state[lane(x, 4)]
		});

		[
			c[4] ^ c[1].rotate_left(1),
			c[0] ^ c[2].rotate_left(1),
			c[1] ^ c[3].rotate_left(1),
			c[2] ^ c[4].rotate_left(1),
			c[3] ^ c[0].rotate_left(1),
		]
	}

	#[test]
	fn theta_corrections_match_reference() {
		let mut rng = StdRng::seed_from_u64(21);

		for _ in 0..32 {
			let state = rng.random::<State>();
			assert_eq!(theta_corrections(&state), theta_corrections_reference(&state));
		}
	}

	#[test]
	fn witness_materializes_a_and_d_with_padding_zero() {
		let mut rng = StdRng::seed_from_u64(22);
		let traces: Vec<_> = (0..2)
			.map(|_| PermutationTrace::new(rng.random::<State>()))
			.collect();
		let witness = CommittedKeccakWitness::from_traces(&traces);

		assert_eq!(witness.n_permutations(), traces.len());
		assert_eq!(witness.words().len(), layout::witness_len(traces.len()));

		for (permutation, trace) in traces.iter().enumerate() {
			for round in 0..N_ROUNDS {
				for lane in 0..N_LANES {
					assert_eq!(
						witness.a(permutation, round, lane),
						Word(trace.rounds[round].input[lane]),
					);
				}

				let d = theta_corrections(&trace.rounds[round].input);
				for (x, &d_x) in d.iter().enumerate() {
					assert_eq!(witness.d(permutation, round, x), Word(d_x));
				}
			}

			for lane in 0..N_LANES {
				assert_eq!(witness.a(permutation, N_ROUNDS, lane), Word(trace.output()[lane]),);
			}

			for round_block in 0..layout::BLOCKS_PER_PERMUTATION {
				for padding_slot in 0..layout::PADDING_WORDS_PER_BLOCK {
					assert_eq!(
						witness.words()
							[layout::padding_index(permutation, round_block, padding_slot)],
						Word::ZERO,
					);
				}
			}

			for x in 0..layout::D_WORDS_PER_BLOCK {
				assert_eq!(
					witness.words()
						[layout::block_slot_index(permutation, N_ROUNDS, layout::D_OFFSET + x,)],
					Word::ZERO,
				);
			}
		}
	}
}
