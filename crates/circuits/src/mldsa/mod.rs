// Copyright 2025 Irreducible Inc.

//! Prototype ML-DSA bit-heavy circuit helpers.
//!
//! Public-only hashing and key preprocessing are intentionally hoisted out of this circuit layer.

mod hashing;
mod hint;
mod packing;
mod params;
mod relations;
mod sample_in_ball;
#[cfg(test)]
mod tests;
mod types;
mod use_hint;
mod util;

pub use hashing::final_challenge_hash_for;
pub use hint::{
	assert_h_bits_and_weight_for, assert_hint_canonical_matches_expanded_for,
	decode_hint_canonical_for,
};
pub use packing::{
	assert_z_norm_from_packed_y_for, assert_z_packed_bytes_norm_for, decode_z_packed_y_for,
	encode_w1_for,
};
pub use params::{Mldsa44, Mldsa65, Mldsa87, MldsaParams, mldsa44, mldsa65, mldsa87};
pub use relations::{
	full_bit_heavy_one_block_canonical_hint_matched_relation_for,
	full_bit_heavy_one_block_canonical_hint_relation_for, full_bit_heavy_one_block_relation_for,
	one_block_canonical_hint_hash_relation_for, one_block_hidden_hash_relation_for,
	one_block_use_hint_hash_relation_for, one_block_w1encode_hash_relation_for,
};
pub use sample_in_ball::{
	sample_in_ball_one_block_for, sample_in_ball_one_block_from_stream_for,
	sample_in_ball_one_block_sparse_for, sample_in_ball_one_block_sparse_from_stream_for,
	sample_in_ball_one_block_stream_for,
};
pub use types::{MldsaBitHeavyCircuit, MldsaSampleInBallOneBlock, MldsaSampleInBallOneBlockSparse};
pub use use_hint::{high_bits_for, use_hint_coeff_for, use_hint_for};
