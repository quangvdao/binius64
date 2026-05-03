// Copyright 2026 The Binius Developers

//! Committed witness layout for the specialized Keccak prover.

use crate::constants::{N_LANES, N_ROUNDS};

/// Number of committed round-boundary blocks per Keccak-f[1600] permutation.
pub const BLOCKS_PER_PERMUTATION: usize = N_ROUNDS + 1;

/// Number of committed words per round-boundary block.
pub const WORDS_PER_BLOCK: usize = 32;

/// Offset of the first round-boundary state lane in a block.
pub const A_OFFSET: usize = 0;

/// Number of committed round-boundary state lanes in a block.
pub const A_WORDS_PER_BLOCK: usize = N_LANES;

/// Offset of the first theta correction word in a block.
pub const D_OFFSET: usize = A_OFFSET + A_WORDS_PER_BLOCK;

/// Number of theta correction words per non-final block.
pub const D_WORDS_PER_BLOCK: usize = 5;

/// Offset of the first padding word in a block.
pub const PADDING_OFFSET: usize = D_OFFSET + D_WORDS_PER_BLOCK;

/// Number of padding words per block.
pub const PADDING_WORDS_PER_BLOCK: usize = WORDS_PER_BLOCK - PADDING_OFFSET;

/// Number of committed word slots per Keccak-f[1600] permutation.
pub const WORDS_PER_PERMUTATION: usize = BLOCKS_PER_PERMUTATION * WORDS_PER_BLOCK;

/// Number of active committed state words per permutation.
pub const A_WORDS_PER_PERMUTATION: usize = BLOCKS_PER_PERMUTATION * A_WORDS_PER_BLOCK;

/// Number of active theta correction words per permutation.
pub const D_WORDS_PER_PERMUTATION: usize = N_ROUNDS * D_WORDS_PER_BLOCK;

/// Number of active committed words per permutation.
pub const ACTIVE_WORDS_PER_PERMUTATION: usize = A_WORDS_PER_PERMUTATION + D_WORDS_PER_PERMUTATION;

/// Number of inactive committed padding slots per permutation.
pub const PADDING_WORDS_PER_PERMUTATION: usize =
	WORDS_PER_PERMUTATION - ACTIVE_WORDS_PER_PERMUTATION;

/// Return the first committed word index for `(permutation, round_block)`.
#[inline]
pub fn block_index(permutation: usize, round_block: usize) -> usize {
	assert!(round_block < BLOCKS_PER_PERMUTATION);
	permutation * WORDS_PER_PERMUTATION + round_block * WORDS_PER_BLOCK
}

/// Return the committed word index for one block-local slot.
#[inline]
pub fn block_slot_index(permutation: usize, round_block: usize, slot: usize) -> usize {
	assert!(slot < WORDS_PER_BLOCK);
	block_index(permutation, round_block) + slot
}

/// Return the committed word index for `A_r[lane]`.
#[inline]
pub fn a_index(permutation: usize, round_block: usize, lane: usize) -> usize {
	assert!(lane < A_WORDS_PER_BLOCK);
	block_slot_index(permutation, round_block, A_OFFSET + lane)
}

/// Return the committed word index for `D_r[x]`.
#[inline]
pub fn d_index(permutation: usize, round: usize, x: usize) -> usize {
	assert!(round < N_ROUNDS);
	assert!(x < D_WORDS_PER_BLOCK);
	block_slot_index(permutation, round, D_OFFSET + x)
}

/// Return the committed word index for one padding slot.
#[inline]
pub fn padding_index(permutation: usize, round_block: usize, padding_slot: usize) -> usize {
	assert!(padding_slot < PADDING_WORDS_PER_BLOCK);
	block_slot_index(permutation, round_block, PADDING_OFFSET + padding_slot)
}

/// Return the committed witness length for a batch of permutations.
#[inline]
pub fn witness_len(n_permutations: usize) -> usize {
	n_permutations * WORDS_PER_PERMUTATION
}

#[cfg(test)]
mod tests {
	use super::*;

	#[test]
	fn constants_match_locked_layout() {
		assert_eq!(BLOCKS_PER_PERMUTATION, 25);
		assert_eq!(WORDS_PER_BLOCK, 32);
		assert_eq!(A_OFFSET, 0);
		assert_eq!(A_WORDS_PER_BLOCK, 25);
		assert_eq!(D_OFFSET, 25);
		assert_eq!(D_WORDS_PER_BLOCK, 5);
		assert_eq!(PADDING_OFFSET, 30);
		assert_eq!(PADDING_WORDS_PER_BLOCK, 2);
		assert_eq!(WORDS_PER_PERMUTATION, 800);
		assert_eq!(A_WORDS_PER_PERMUTATION, 625);
		assert_eq!(D_WORDS_PER_PERMUTATION, 120);
		assert_eq!(ACTIVE_WORDS_PER_PERMUTATION, 745);
		assert_eq!(PADDING_WORDS_PER_PERMUTATION, 55);
	}

	#[test]
	fn indices_follow_block_formula() {
		assert_eq!(block_index(0, 0), 0);
		assert_eq!(block_index(0, 1), 32);
		assert_eq!(block_index(1, 0), 800);
		assert_eq!(a_index(1, 2, 7), 800 + 2 * 32 + 7);
		assert_eq!(d_index(1, 2, 3), 800 + 2 * 32 + 25 + 3);
		assert_eq!(padding_index(1, 24, 1), 800 + 24 * 32 + 31);
		assert_eq!(witness_len(3), 2400);
	}
}
