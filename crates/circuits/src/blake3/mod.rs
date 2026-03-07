// Copyright 2026 The Binius Developers
//! BLAKE3 hash function circuit implementation
//!
//! BLAKE3 is a cryptographic hash function that uses the same G mixing function
//! as BLAKE2s but with 7 rounds (vs 10), a fixed message permutation (vs SIGMA
//! table), and a Merkle tree structure for multi-chunk parallelism.
//!
//! This implementation supports arbitrary-length messages via BLAKE3's tree hashing
//! mode. Messages up to 1024 bytes use a single chunk; larger messages are split
//! into 1024-byte chunks, compressed independently, and merged in a left-leaning
//! binary Merkle tree using parent compressions.
//!
//! # Differences from BLAKE2s
//!
//! - **Rounds**: 7 instead of 10
//! - **Message schedule**: Fixed `MSG_PERMUTATION` applied between rounds
//! - **State init v\[12..16\]**: Direct assign `(counter_lo, counter_hi, block_len, flags)`
//!   -- no XOR with IV
//! - **Finalization**: `state[i] ^= state[i+8]; state[i+8] ^= cv[i]` producing
//!   16-word output (only first 8 used for chaining values)
//! - **Counter**: Chunk index, not cumulative byte count
//! - **Flags**: Bitfield with `CHUNK_START`, `CHUNK_END`, `PARENT`, `ROOT` domain
//!   separation
//!
//! # Tree Hashing
//!
//! For messages larger than one chunk (1024 bytes):
//! 1. Split message into 1024-byte chunks
//! 2. Compress each chunk independently (counter = chunk index)
//! 3. Merge chaining values in a left-leaning binary tree using parent compressions
//! 4. The root compression (final parent or sole chunk) gets the `ROOT` flag
//!
//! # Constraint Cost
//!
//! - Per G function: 12 AND constraints (same as BLAKE2s)
//! - Per round: 8 G calls = 96 AND constraints
//! - Per compression: 7 rounds = 672 AND constraints
//! - Per parent compression: same as one compression (single 64-byte block)
//! - For N chunks: (2N-1) compressions total (N chunks + N-1 parent merges)
//! - 30% fewer constraints per compression than BLAKE2s (672 vs 960)

pub mod constants;
pub mod reference;
#[cfg(test)]
mod tests;

use binius_core::word::Word;
use binius_frontend::{CircuitBuilder, Wire, WitnessFiller};

use crate::blake2s::g_function;
use constants::{
	BLOCK_LEN, CHUNK_END, CHUNK_LEN, CHUNK_START, IV, PARENT, ROOT, ROUNDS, SCHEDULES,
};

/// BLAKE3 compression function circuit.
///
/// Compresses a single 64-byte block with the given chaining value, counter,
/// block length, and domain separation flags. Returns the full 16-word state
/// (callers take `first_8_words` for chaining, or use all 16 for root output).
///
/// # Arguments
///
/// * `chaining_value` - 8 input chaining value words (IV for first block in chunk)
/// * `block_words` - 16 message words (one 64-byte block)
/// * `counter` - Chunk counter as two 32-bit words `(lo, hi)`
/// * `block_len` - Actual byte count in this block (as a wire)
/// * `flags` - Domain separation flags (as a wire)
///
/// # Constraint Cost
///
/// 7 rounds x 8 G calls x 12 AND = **672 AND constraints**
pub fn blake3_compress(
	builder: &mut CircuitBuilder,
	chaining_value: &[Wire; 8],
	block_words: &[Wire; 16],
	counter_lo: Wire,
	counter_hi: Wire,
	block_len: Wire,
	flags: Wire,
) -> [Wire; 16] {
	let mut v = [builder.add_constant(Word(0)); 16];

	// v[0..8] = chaining_value
	v[0..8].copy_from_slice(chaining_value);

	// v[8..12] = IV[0..4]  (NOT all 8 IVs -- BLAKE3 differs from BLAKE2s here)
	for i in 0..4 {
		v[8 + i] = builder.add_constant(Word(IV[i] as u64));
	}

	// v[12..16] = (counter_lo, counter_hi, block_len, flags)
	// Note: NO XOR with IV -- direct assignment (unlike BLAKE2s)
	v[12] = counter_lo;
	v[13] = counter_hi;
	v[14] = block_len;
	v[15] = flags;

	// 7 rounds of G function mixing with precomputed message schedules
	for round in 0..ROUNDS {
		let s = &SCHEDULES[round];

		// Column step
		let (v0, v4, v8, v12) = g_function(
			builder,
			v[0],
			v[4],
			v[8],
			v[12],
			block_words[s[0]],
			block_words[s[1]],
		);
		v[0] = v0;
		v[4] = v4;
		v[8] = v8;
		v[12] = v12;

		let (v1, v5, v9, v13) = g_function(
			builder,
			v[1],
			v[5],
			v[9],
			v[13],
			block_words[s[2]],
			block_words[s[3]],
		);
		v[1] = v1;
		v[5] = v5;
		v[9] = v9;
		v[13] = v13;

		let (v2, v6, v10, v14) = g_function(
			builder,
			v[2],
			v[6],
			v[10],
			v[14],
			block_words[s[4]],
			block_words[s[5]],
		);
		v[2] = v2;
		v[6] = v6;
		v[10] = v10;
		v[14] = v14;

		let (v3, v7, v11, v15) = g_function(
			builder,
			v[3],
			v[7],
			v[11],
			v[15],
			block_words[s[6]],
			block_words[s[7]],
		);
		v[3] = v3;
		v[7] = v7;
		v[11] = v11;
		v[15] = v15;

		// Diagonal step
		let (v0, v5, v10, v15) = g_function(
			builder,
			v[0],
			v[5],
			v[10],
			v[15],
			block_words[s[8]],
			block_words[s[9]],
		);
		v[0] = v0;
		v[5] = v5;
		v[10] = v10;
		v[15] = v15;

		let (v1, v6, v11, v12) = g_function(
			builder,
			v[1],
			v[6],
			v[11],
			v[12],
			block_words[s[10]],
			block_words[s[11]],
		);
		v[1] = v1;
		v[6] = v6;
		v[11] = v11;
		v[12] = v12;

		let (v2, v7, v8, v13) = g_function(
			builder,
			v[2],
			v[7],
			v[8],
			v[13],
			block_words[s[12]],
			block_words[s[13]],
		);
		v[2] = v2;
		v[7] = v7;
		v[8] = v8;
		v[13] = v13;

		let (v3, v4, v9, v14) = g_function(
			builder,
			v[3],
			v[4],
			v[9],
			v[14],
			block_words[s[14]],
			block_words[s[15]],
		);
		v[3] = v3;
		v[4] = v4;
		v[9] = v9;
		v[14] = v14;
	}

	// BLAKE3 finalization: state[i] ^= state[i+8]; state[i+8] ^= cv[i]
	for i in 0..8 {
		v[i] = builder.bxor(v[i], v[i + 8]);
		v[i + 8] = builder.bxor(v[i + 8], chaining_value[i]);
	}

	v
}

/// BLAKE3 chunk compression circuit.
///
/// Compresses a single chunk (up to 1024 bytes / 16 blocks) by iterating through
/// its blocks with the correct flags and chunk counter. Returns the 8-word
/// chaining value for this chunk.
///
/// # Arguments
///
/// * `message_qwords` - 64-bit packed message words for this chunk (each qword
///   holds two 32-bit BLAKE3 words)
/// * `chunk_byte_len` - Actual byte count in this chunk (at most `CHUNK_LEN`)
/// * `counter` - Chunk index (0-based)
/// * `is_root` - If true, the last block gets `CHUNK_END | ROOT`; otherwise just
///   `CHUNK_END`
pub fn blake3_chunk(
	builder: &mut CircuitBuilder,
	message_qwords: &[Wire],
	chunk_byte_len: usize,
	counter: u64,
	is_root: bool,
) -> [Wire; 8] {
	let num_blocks = chunk_byte_len.div_ceil(BLOCK_LEN).max(1);
	let zero = builder.add_constant(Word(0));

	let iv_wires: [Wire; 8] =
		std::array::from_fn(|i| builder.add_constant(Word(IV[i] as u64)));
	let mut cv = iv_wires;

	for block_idx in 0..num_blocks {
		let mut m = [zero; 16];

		for word_idx in 0..16 {
			let qword_idx = block_idx * 8 + word_idx / 2;
			let message_qword = message_qwords.get(qword_idx).copied().unwrap_or(zero);

			let message_dword = if word_idx % 2 == 0 {
				builder.band(message_qword, builder.add_constant_64(0xFFFF_FFFF))
			} else {
				builder.shr(message_qword, 32)
			};

			let first_byte_offset = block_idx * BLOCK_LEN + word_idx * 4;
			let padded_message_dword = if first_byte_offset + 4 > chunk_byte_len {
				if first_byte_offset < chunk_byte_len {
					let nonzero_bytes = (chunk_byte_len - first_byte_offset) as u32;
					builder.band(
						message_dword,
						builder.add_constant(Word::ALL_ONE >> (64 - nonzero_bytes * 8)),
					)
				} else {
					zero
				}
			} else {
				message_dword
			};

			m[word_idx] = padded_message_dword;
		}

		let is_first = block_idx == 0;
		let is_last = block_idx == num_blocks - 1;

		let mut flag_bits = 0u32;
		if is_first {
			flag_bits |= CHUNK_START;
		}
		if is_last {
			flag_bits |= CHUNK_END;
			if is_root {
				flag_bits |= ROOT;
			}
		}

		let bytes_in_block = if is_last {
			let rem = chunk_byte_len % BLOCK_LEN;
			if rem == 0 && chunk_byte_len > 0 {
				BLOCK_LEN
			} else if chunk_byte_len == 0 {
				0
			} else {
				rem
			}
		} else {
			BLOCK_LEN
		};

		let counter_lo = builder.add_constant(Word(counter & 0xFFFF_FFFF));
		let counter_hi = builder.add_constant(Word(counter >> 32));
		let block_len_wire = builder.add_constant(Word(bytes_in_block as u64));
		let flags_wire = builder.add_constant(Word(flag_bits as u64));

		let output = blake3_compress(
			builder,
			&cv,
			&m,
			counter_lo,
			counter_hi,
			block_len_wire,
			flags_wire,
		);

		cv = output[0..8].try_into().unwrap();
	}

	cv
}

/// BLAKE3 parent node compression circuit.
///
/// Merges two child chaining values by compressing `left_cv || right_cv` as a
/// single 64-byte block with IV as the chaining value and the `PARENT` flag.
///
/// # Arguments
///
/// * `left_cv` - 8-word chaining value from the left child
/// * `right_cv` - 8-word chaining value from the right child
/// * `is_root` - If true, adds `ROOT` to the flags (for the final merge)
pub fn blake3_parent(
	builder: &mut CircuitBuilder,
	left_cv: &[Wire; 8],
	right_cv: &[Wire; 8],
	is_root: bool,
) -> [Wire; 8] {
	let iv_wires: [Wire; 8] =
		std::array::from_fn(|i| builder.add_constant(Word(IV[i] as u64)));
	let zero = builder.add_constant(Word(0));

	let mut block_words = [zero; 16];
	block_words[..8].copy_from_slice(left_cv);
	block_words[8..].copy_from_slice(right_cv);

	let mut flags = PARENT;
	if is_root {
		flags |= ROOT;
	}

	let block_len_wire = builder.add_constant(Word(BLOCK_LEN as u64));
	let flags_wire = builder.add_constant(Word(flags as u64));

	let output = blake3_compress(
		builder,
		&iv_wires,
		&block_words,
		zero,
		zero,
		block_len_wire,
		flags_wire,
	);

	output[0..8].try_into().unwrap()
}

/// Recursively merge chaining values in a left-leaning binary tree.
///
/// The left subtree always gets the largest power-of-2 count that fits,
/// matching the official BLAKE3 tree structure.
fn merge_cvs(
	builder: &mut CircuitBuilder,
	cvs: &[[Wire; 8]],
	is_root: bool,
) -> [Wire; 8] {
	assert!(!cvs.is_empty());
	if cvs.len() == 1 {
		return cvs[0];
	}
	let split = cvs.len().next_power_of_two() >> 1;
	let left = merge_cvs(builder, &cvs[..split], false);
	let right = merge_cvs(builder, &cvs[split..], false);
	blake3_parent(builder, &left, &right, is_root)
}

/// BLAKE3 hash function circuit.
///
/// Verifies that a message of fixed `length` produces a specific 256-bit digest.
/// Supports arbitrary-length messages: single-chunk (up to 1024 bytes) uses
/// sequential block compression, multi-chunk uses tree hashing with parent merges.
pub struct Blake3 {
	/// Message size in bytes this circuit supports
	pub length: usize,
	/// Witness wires for the input message (little-endian packed into 64-bit words)
	pub message: Vec<Wire>,
	/// Witness wires for the expected 256-bit digest (8 x 32-bit words)
	pub digest: [Wire; 8],
}

impl Blake3 {
	/// Create a new BLAKE3 circuit with witness variables.
	///
	/// # Arguments
	///
	/// * `builder` - Circuit builder to add constraints to
	/// * `length` - Fixed message size in bytes (any length supported)
	pub fn new_witness(builder: &mut CircuitBuilder, length: usize) -> Self {
		let message: Vec<Wire> = (0..length.div_ceil(8))
			.map(|_| builder.add_witness())
			.collect();
		let digest = std::array::from_fn(|_| builder.add_witness());

		Self::build_circuit(builder, length, &message, digest);

		Self {
			length,
			message,
			digest,
		}
	}

	/// Build the BLAKE3 circuit constraints.
	///
	/// For single-chunk messages (up to 1024 bytes), compresses blocks
	/// sequentially within one chunk. For multi-chunk messages, compresses
	/// each chunk independently and merges chaining values via a left-leaning
	/// binary tree of parent compressions.
	fn build_circuit(
		builder: &mut CircuitBuilder,
		length: usize,
		message: &[Wire],
		expected_digest: [Wire; 8],
	) {
		let num_chunks = length.div_ceil(CHUNK_LEN).max(1);

		let mut chunk_cvs = Vec::with_capacity(num_chunks);

		for chunk_idx in 0..num_chunks {
			let chunk_byte_start = chunk_idx * CHUNK_LEN;
			let chunk_byte_len = if length == 0 {
				0
			} else {
				CHUNK_LEN.min(length - chunk_byte_start)
			};

			let qword_start = chunk_byte_start / 8;
			let num_qwords = chunk_byte_len.div_ceil(8);
			let chunk_qwords = if num_qwords > 0 {
				&message[qword_start..qword_start + num_qwords]
			} else {
				&[]
			};

			let is_single_chunk = num_chunks == 1;
			let cv = blake3_chunk(
				builder,
				chunk_qwords,
				chunk_byte_len,
				chunk_idx as u64,
				is_single_chunk,
			);
			chunk_cvs.push(cv);
		}

		let root_cv = if chunk_cvs.len() == 1 {
			chunk_cvs[0]
		} else {
			merge_cvs(builder, &chunk_cvs, true)
		};

		for i in 0..8 {
			builder.assert_eq("blake3_digest_match", root_cv[i], expected_digest[i]);
		}
	}

	/// Populate the message witness data.
	pub fn populate_message(&self, witness: &mut WitnessFiller, message: &[u8]) {
		assert!(
			message.len() == self.length,
			"Only messages of length {} supported, given {} bytes",
			self.length,
			message.len(),
		);

		for (i, bytes) in message.chunks(8).enumerate() {
			let mut le_bytes = [0; 8];
			le_bytes[..bytes.len()].copy_from_slice(bytes);
			witness[self.message[i]] = Word(u64::from_le_bytes(le_bytes));
		}
	}

	/// Populate the expected digest witness data.
	pub fn populate_digest(&self, witness: &mut WitnessFiller, digest: &[u8; 32]) {
		for i in 0..8 {
			let word_bytes = &digest[i * 4..(i + 1) * 4];
			let word = u32::from_le_bytes(word_bytes.try_into().unwrap());
			witness[self.digest[i]] = Word(word as u64);
		}
	}
}
