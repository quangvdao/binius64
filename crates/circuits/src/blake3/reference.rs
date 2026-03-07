// Copyright 2025 Irreducible Inc.
//! BLAKE3 reference implementation for testing.
//!
//! Ported from the official BLAKE3 reference implementation at
//! <https://github.com/BLAKE3-team/BLAKE3/blob/master/reference_impl/reference_impl.rs>.
//!
//! Exposes intermediate compression state for isolation testing.

use super::constants::*;

fn g(state: &mut [u32; 16], a: usize, b: usize, c: usize, d: usize, mx: u32, my: u32) {
	state[a] = state[a].wrapping_add(state[b]).wrapping_add(mx);
	state[d] = (state[d] ^ state[a]).rotate_right(16);
	state[c] = state[c].wrapping_add(state[d]);
	state[b] = (state[b] ^ state[c]).rotate_right(12);
	state[a] = state[a].wrapping_add(state[b]).wrapping_add(my);
	state[d] = (state[d] ^ state[a]).rotate_right(8);
	state[c] = state[c].wrapping_add(state[d]);
	state[b] = (state[b] ^ state[c]).rotate_right(7);
}

fn round(state: &mut [u32; 16], m: &[u32; 16]) {
	g(state, 0, 4, 8, 12, m[0], m[1]);
	g(state, 1, 5, 9, 13, m[2], m[3]);
	g(state, 2, 6, 10, 14, m[4], m[5]);
	g(state, 3, 7, 11, 15, m[6], m[7]);
	g(state, 0, 5, 10, 15, m[8], m[9]);
	g(state, 1, 6, 11, 12, m[10], m[11]);
	g(state, 2, 7, 8, 13, m[12], m[13]);
	g(state, 3, 4, 9, 14, m[14], m[15]);
}

fn permute(m: &mut [u32; 16]) {
	let mut permuted = [0; 16];
	for i in 0..16 {
		permuted[i] = m[MSG_PERMUTATION[i]];
	}
	*m = permuted;
}

/// Result of a compression call, exposing intermediate round states.
pub struct CompressTrace {
	/// State after each of the 7 rounds (before finalization).
	pub round_states: [[u32; 16]; ROUNDS],
	/// Final 16-word output after finalization XORs.
	pub output: [u32; 16],
}

/// BLAKE3 compression function with intermediate state tracing.
pub fn compress_traced(
	chaining_value: &[u32; 8],
	block_words: &[u32; 16],
	counter: u64,
	block_len: u32,
	flags: u32,
) -> CompressTrace {
	let counter_low = counter as u32;
	let counter_high = (counter >> 32) as u32;

	#[rustfmt::skip]
	let mut state = [
		chaining_value[0], chaining_value[1], chaining_value[2], chaining_value[3],
		chaining_value[4], chaining_value[5], chaining_value[6], chaining_value[7],
		IV[0], IV[1], IV[2], IV[3],
		counter_low, counter_high, block_len, flags,
	];

	let mut block = *block_words;
	let mut round_states = [[0u32; 16]; ROUNDS];

	for r in 0..ROUNDS {
		round(&mut state, &block);
		round_states[r] = state;
		permute(&mut block);
	}

	for i in 0..8 {
		state[i] ^= state[i + 8];
		state[i + 8] ^= chaining_value[i];
	}

	CompressTrace {
		round_states,
		output: state,
	}
}

/// BLAKE3 compression function (simple version, returns output only).
pub fn compress(
	chaining_value: &[u32; 8],
	block_words: &[u32; 16],
	counter: u64,
	block_len: u32,
	flags: u32,
) -> [u32; 16] {
	compress_traced(chaining_value, block_words, counter, block_len, flags).output
}

fn first_8_words(output: &[u32; 16]) -> [u32; 8] {
	output[0..8].try_into().unwrap()
}

fn words_from_little_endian_bytes(bytes: &[u8], words: &mut [u32]) {
	for (four_bytes, word) in bytes.chunks_exact(4).zip(words) {
		*word = u32::from_le_bytes(four_bytes.try_into().unwrap());
	}
}

/// Compress a single chunk (up to 1024 bytes) and return its 8-word chaining value.
///
/// Processes all blocks in the chunk with the correct flags and counter.
/// If `is_root` is true, the last block gets `CHUNK_END | ROOT`.
pub fn chunk_cv(chunk_bytes: &[u8], counter: u64, is_root: bool) -> [u32; 8] {
	let chunk_len = chunk_bytes.len();
	let num_blocks = chunk_len.div_ceil(BLOCK_LEN).max(1);
	let mut cv = IV;

	for block_idx in 0..num_blocks {
		let is_first = block_idx == 0;
		let is_last = block_idx == num_blocks - 1;

		let block_start = block_idx * BLOCK_LEN;
		let block_end = (block_start + BLOCK_LEN).min(chunk_len);
		let actual_bytes = if chunk_bytes.is_empty() {
			0
		} else {
			block_end - block_start
		};

		let mut block_bytes = [0u8; BLOCK_LEN];
		if actual_bytes > 0 {
			block_bytes[..actual_bytes].copy_from_slice(&chunk_bytes[block_start..block_end]);
		}

		let mut block_words = [0u32; 16];
		words_from_little_endian_bytes(&block_bytes, &mut block_words);

		let mut flags = 0u32;
		if is_first {
			flags |= CHUNK_START;
		}
		if is_last {
			flags |= CHUNK_END;
			if is_root {
				flags |= ROOT;
			}
		}

		let output = compress(&cv, &block_words, counter, actual_bytes as u32, flags);
		cv = first_8_words(&output);
	}

	cv
}

/// Compute a parent chaining value from two child chaining values.
///
/// Uses `compress(IV, left || right, 0, BLOCK_LEN, PARENT [| ROOT])`.
pub fn parent_cv(left: &[u32; 8], right: &[u32; 8], is_root: bool) -> [u32; 8] {
	let mut block_words = [0u32; 16];
	block_words[..8].copy_from_slice(left);
	block_words[8..].copy_from_slice(right);

	let mut flags = PARENT;
	if is_root {
		flags |= ROOT;
	}

	let output = compress(&IV, &block_words, 0, BLOCK_LEN as u32, flags);
	first_8_words(&output)
}

/// Recursively merge chaining values in a left-leaning binary tree.
///
/// Left subtree gets the largest power-of-2 count that fits.
fn merge_cvs(cvs: &[[u32; 8]], is_root: bool) -> [u32; 8] {
	assert!(!cvs.is_empty());
	if cvs.len() == 1 {
		return cvs[0];
	}
	let split = cvs.len().next_power_of_two() >> 1;
	let left = merge_cvs(&cvs[..split], false);
	let right = merge_cvs(&cvs[split..], false);
	parent_cv(&left, &right, is_root)
}

/// Compute BLAKE3 hash of a message (arbitrary length).
///
/// Returns the 32-byte digest. Supports single-chunk and multi-chunk
/// tree hashing.
pub fn blake3_hash(message: &[u8]) -> [u8; 32] {
	let num_chunks = message.len().div_ceil(CHUNK_LEN).max(1);

	let chunk_cvs: Vec<[u32; 8]> = (0..num_chunks)
		.map(|chunk_idx| {
			let start = chunk_idx * CHUNK_LEN;
			let end = (start + CHUNK_LEN).min(message.len());
			let chunk_bytes = if message.is_empty() { &[] } else { &message[start..end] };
			let is_single_chunk = num_chunks == 1;
			chunk_cv(chunk_bytes, chunk_idx as u64, is_single_chunk)
		})
		.collect();

	let root_cv = if chunk_cvs.len() == 1 {
		chunk_cvs[0]
	} else {
		merge_cvs(&chunk_cvs, true)
	};

	let mut digest = [0u8; 32];
	for (i, word) in root_cv.iter().enumerate() {
		digest[i * 4..(i + 1) * 4].copy_from_slice(&word.to_le_bytes());
	}
	digest
}
