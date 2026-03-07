// Copyright 2026 The Binius Developers
//! BLAKE3 constants
//!
//! Constants for the BLAKE3 hash function, including initialization vectors,
//! message permutation, precomputed round schedules, and domain separation flags.

/// BLAKE3 block size in bytes.
pub const BLOCK_LEN: usize = 64;

/// BLAKE3 chunk size in bytes (16 blocks per chunk).
pub const CHUNK_LEN: usize = 1024;

/// Number of rounds in the BLAKE3 compression function.
pub const ROUNDS: usize = 7;

/// BLAKE3 initialization vectors (same as SHA-256 / BLAKE2s).
pub const IV: [u32; 8] = [
	0x6A09E667, 0xBB67AE85, 0x3C6EF372, 0xA54FF53A, 0x510E527F, 0x9B05688C, 0x1F83D9AB,
	0x5BE0CD19,
];

/// Domain separation flags for BLAKE3 compression.
pub const CHUNK_START: u32 = 1 << 0;
pub const CHUNK_END: u32 = 1 << 1;
pub const PARENT: u32 = 1 << 2;
pub const ROOT: u32 = 1 << 3;

/// Fixed message word permutation applied between rounds.
///
/// Unlike BLAKE2s which uses a 10-row SIGMA table, BLAKE3 applies this single
/// permutation iteratively. Round `i` uses message words permuted by applying
/// this permutation `i` times to the identity ordering.
pub const MSG_PERMUTATION: [usize; 16] = [2, 6, 3, 10, 7, 0, 4, 13, 1, 11, 12, 5, 9, 14, 15, 8];

/// Precomputed message schedules for all 7 rounds.
///
/// `SCHEDULES[i]` gives the message word ordering for round `i`, computed by
/// applying `MSG_PERMUTATION` `i` times to the identity `[0, 1, ..., 15]`.
/// The 7th permutation result (after round 6) is discarded per the spec.
///
/// These are verified at compile time via `const` evaluation.
pub const SCHEDULES: [[usize; 16]; ROUNDS] = compute_schedules();

const fn compute_schedules() -> [[usize; 16]; ROUNDS] {
	let mut schedules = [[0usize; 16]; ROUNDS];

	// Round 0: identity
	let mut i = 0;
	while i < 16 {
		schedules[0][i] = i;
		i += 1;
	}

	// Rounds 1..6: iteratively apply MSG_PERMUTATION
	let mut r = 1;
	while r < ROUNDS {
		let mut j = 0;
		while j < 16 {
			schedules[r][j] = schedules[r - 1][MSG_PERMUTATION[j]];
			j += 1;
		}
		r += 1;
	}

	// Compile-time verification: each schedule must be a permutation of 0..16
	let mut r = 0;
	while r < ROUNDS {
		let mut seen = [false; 16];
		let mut j = 0;
		while j < 16 {
			assert!(schedules[r][j] < 16, "schedule index out of range");
			assert!(!seen[schedules[r][j]], "duplicate index in schedule");
			seen[schedules[r][j]] = true;
			j += 1;
		}
		r += 1;
	}

	schedules
}

#[cfg(test)]
mod tests {
	use super::*;

	#[test]
	fn verify_schedules_against_iterative_permute() {
		fn permute(m: &[usize; 16]) -> [usize; 16] {
			let mut out = [0; 16];
			for i in 0..16 {
				out[i] = m[MSG_PERMUTATION[i]];
			}
			out
		}

		let mut m: [usize; 16] = core::array::from_fn(|i| i);
		for r in 0..ROUNDS {
			assert_eq!(SCHEDULES[r], m, "schedule mismatch at round {r}");
			m = permute(&m);
		}
	}

	#[test]
	fn verify_iv_matches_blake2s() {
		let blake2s_iv = crate::blake2s::constants::IV;
		assert_eq!(IV, blake2s_iv);
	}
}
