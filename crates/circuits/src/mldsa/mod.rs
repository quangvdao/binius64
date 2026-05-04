// Copyright 2025 Irreducible Inc.

//! Prototype ML-DSA bit-heavy circuit helpers.
//!
//! This module starts with the ML-DSA-44 hidden-signature target shape from the
//! top-level lattice-sig-aggregation executable spec. Public-only hashing and key
//! preprocessing are intentionally hoisted out of this circuit layer.

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

pub use hashing::mldsa44_final_challenge_hash;
pub use hint::{
	assert_mldsa44_h_bits_and_weight, assert_mldsa44_hint_canonical_matches_expanded,
	mldsa44_decode_hint_canonical,
};
pub use packing::{
	assert_mldsa44_z_norm_from_packed_y, assert_mldsa44_z_packed_bytes_norm,
	mldsa44_decode_z_packed_y, mldsa44_encode_w1,
};
pub use params::{Mldsa44, MldsaParams, mldsa44};
pub use relations::{
	mldsa44_full_bit_heavy_one_block_canonical_hint_matched_relation,
	mldsa44_full_bit_heavy_one_block_canonical_hint_relation,
	mldsa44_full_bit_heavy_one_block_relation, mldsa44_one_block_canonical_hint_hash_relation,
	mldsa44_one_block_hidden_hash_relation, mldsa44_one_block_use_hint_hash_relation,
	mldsa44_one_block_w1encode_hash_relation,
};
pub use sample_in_ball::{
	mldsa44_sample_in_ball_one_block, mldsa44_sample_in_ball_one_block_from_stream,
	mldsa44_sample_in_ball_one_block_sparse, mldsa44_sample_in_ball_one_block_sparse_from_stream,
	mldsa44_sample_in_ball_one_block_stream,
};
pub use types::{
	Mldsa44BitHeavyCircuit, Mldsa44SampleInBallOneBlock, Mldsa44SampleInBallOneBlockSparse,
};
pub use use_hint::{mldsa44_high_bits, mldsa44_use_hint, mldsa44_use_hint_coeff};
