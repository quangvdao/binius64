// Copyright 2025 Irreducible Inc.

use binius_frontend::{CircuitBuilder, Wire};

use super::{
	MldsaParams,
	hashing::final_challenge_hash_for,
	hint::{assert_hint_canonical_matches_expanded_for, decode_hint_canonical_for},
	packing::{assert_z_packed_bytes_norm_for, encode_w1_for, encode_w1_unchecked_for},
	sample_in_ball::{expand_sparse_sample_in_ball_for, sample_in_ball_one_block_sparse_for},
	types::{MldsaBitHeavyCircuit, MldsaSampleInBallOneBlock, MldsaSampleInBallOneBlockSparse},
	use_hint::{use_hint_checked_for, use_hint_for},
};

/// Binds hidden `c_tilde` to both `SampleInBall` and the ML-DSA-44 final challenge hash.
///
/// This is the first combined bit-heavy relation with the lattice layer still mocked:
///
/// ```text
/// c = SampleInBall(c_tilde)
/// c_tilde == SHAKE256(mu || w1Encode(w1_prime), 32)
/// ```
///
/// `mu_and_w1_bytes` is a fixed 832-byte packed wire slice. In the target aggregate relation, `mu`
/// is a fixed-width public/committed input and `w1Encode(w1_prime)` is derived from the hidden
/// lattice bridge plus `UseHint`; for this circuit phase both arrive as wires.
pub fn one_block_hidden_hash_relation_for<P: MldsaParams>(
	builder: &CircuitBuilder,
	c_tilde: &[Wire],
	mu_and_w1_bytes: &[Wire],
) -> MldsaSampleInBallOneBlock {
	let sparse = one_block_hidden_hash_sparse_relation_for::<P>(builder, c_tilde, mu_and_w1_bytes);
	let coeffs = expand_sparse_sample_in_ball_for::<P>(builder, &sparse.positions, &sparse.signs);
	MldsaSampleInBallOneBlock {
		coeffs,
		draw_counts: sparse.draw_counts,
	}
}

fn one_block_hidden_hash_sparse_relation_for<P: MldsaParams>(
	builder: &CircuitBuilder,
	c_tilde: &[Wire],
	mu_and_w1_bytes: &[Wire],
) -> MldsaSampleInBallOneBlockSparse {
	let c_tilde_prime = final_challenge_hash_for::<P>(builder, mu_and_w1_bytes);
	for (i, (&expected, &actual)) in c_tilde.iter().zip(c_tilde_prime.iter()).enumerate() {
		builder.assert_eq(format!("{}_c_tilde_final_hash[{i}]", P::label()), expected, actual);
	}

	sample_in_ball_one_block_sparse_for::<P>(builder, c_tilde)
}

/// Binds hidden `c_tilde` to `SampleInBall` and to a final hash whose `w1` bytes are encoded in
/// circuit from hidden `w1_prime` coefficient wires.
///
/// This keeps `UseHint` mocked by taking `w1_prime` directly as witness wires, but removes the
/// earlier shortcut where prepacked `w1` bytes were supplied by the witness.
pub fn one_block_w1encode_hash_relation_for<P: MldsaParams>(
	builder: &CircuitBuilder,
	c_tilde: &[Wire],
	mu_words: &[Wire],
	w1_coeffs: &[Wire],
) -> MldsaSampleInBallOneBlock {
	let sparse =
		one_block_w1encode_hash_sparse_relation_for::<P>(builder, c_tilde, mu_words, w1_coeffs);
	let coeffs = expand_sparse_sample_in_ball_for::<P>(builder, &sparse.positions, &sparse.signs);
	MldsaSampleInBallOneBlock {
		coeffs,
		draw_counts: sparse.draw_counts,
	}
}

fn one_block_w1encode_hash_sparse_relation_for<P: MldsaParams>(
	builder: &CircuitBuilder,
	c_tilde: &[Wire],
	mu_words: &[Wire],
	w1_coeffs: &[Wire],
) -> MldsaSampleInBallOneBlockSparse {
	one_block_w1encode_hash_sparse_relation_impl_for::<P>(
		builder, c_tilde, mu_words, w1_coeffs, true,
	)
}

fn one_block_w1encode_hash_sparse_relation_unchecked_w1_for<P: MldsaParams>(
	builder: &CircuitBuilder,
	c_tilde: &[Wire],
	mu_words: &[Wire],
	w1_coeffs: &[Wire],
) -> MldsaSampleInBallOneBlockSparse {
	one_block_w1encode_hash_sparse_relation_impl_for::<P>(
		builder, c_tilde, mu_words, w1_coeffs, false,
	)
}

fn one_block_w1encode_hash_sparse_relation_impl_for<P: MldsaParams>(
	builder: &CircuitBuilder,
	c_tilde: &[Wire],
	mu_words: &[Wire],
	w1_coeffs: &[Wire],
	check_w1_range: bool,
) -> MldsaSampleInBallOneBlockSparse {
	assert_eq!(mu_words.len(), P::MU_WORDS, "{} mu packed word count mismatch", P::label(),);

	let w1_words = if check_w1_range {
		encode_w1_for::<P>(builder, w1_coeffs)
	} else {
		encode_w1_unchecked_for::<P>(builder, w1_coeffs)
	};
	let mut hash_input = Vec::with_capacity(P::FINAL_CHALLENGE_INPUT_WORDS);
	hash_input.extend_from_slice(mu_words);
	hash_input.extend_from_slice(&w1_words);

	one_block_hidden_hash_sparse_relation_for::<P>(builder, c_tilde, &hash_input)
}

/// Binds hidden `c_tilde` to `SampleInBall` and to a final hash whose `w1` bytes are derived from
/// hidden `h` and `wApprox` through `UseHint`.
pub fn one_block_use_hint_hash_relation_for<P: MldsaParams>(
	builder: &CircuitBuilder,
	c_tilde: &[Wire],
	mu_words: &[Wire],
	h_coeffs: &[Wire],
	w_approx_coeffs: &[Wire],
) -> MldsaSampleInBallOneBlock {
	let w1_coeffs = use_hint_for::<P>(builder, h_coeffs, w_approx_coeffs);
	one_block_w1encode_hash_relation_for::<P>(builder, c_tilde, mu_words, &w1_coeffs)
}

/// Like [`one_block_use_hint_hash_relation_for`], but decodes canonical compressed hint bytes
/// instead of taking expanded hint bits as witness.
pub fn one_block_canonical_hint_hash_relation_for<P: MldsaParams>(
	builder: &CircuitBuilder,
	c_tilde: &[Wire],
	mu_words: &[Wire],
	h_words: &[Wire],
	w_approx_coeffs: &[Wire],
) -> MldsaSampleInBallOneBlock {
	let h_coeffs = decode_hint_canonical_for::<P>(builder, h_words);
	let w1_coeffs = use_hint_checked_for::<P>(builder, &h_coeffs, w_approx_coeffs);
	one_block_w1encode_hash_relation_for::<P>(builder, c_tilde, mu_words, &w1_coeffs)
}

/// Builds the full prototype ML-DSA-44 bit-heavy relation.
///
/// This circuit keeps the signature hidden and assumes public-only preprocessing has been hoisted:
///
/// - `mu` is supplied as fixed-width public/committed input wires by the caller;
/// - `Ahat`, `t1_ntt_shifted`, and the lattice arithmetic are outside this binary circuit;
/// - `wApprox` is the hidden bridge value produced by the lattice layer;
/// - `c_tilde`, packed `z`, and expanded hint bits `h` are hidden signature witness values.
///
/// The relation enforced here is:
///
/// ```text
/// z_packed_y = BitUnpack(z_bytes)
/// ||z||_infty < gamma1 - beta
/// h_i in {0,1}, sum(h_i) <= omega
/// w1_prime = UseHint(h, wApprox)
/// c_tilde = SHAKE256(mu || w1Encode(w1_prime), 32)
/// c = SampleInBall(c_tilde)
/// ```
pub fn full_bit_heavy_one_block_relation_for<P: MldsaParams>(
	builder: &CircuitBuilder,
	c_tilde: &[Wire],
	z_words: &[Wire],
	h_coeffs: &[Wire],
	mu_words: &[Wire],
	w_approx_coeffs: &[Wire],
) -> MldsaBitHeavyCircuit<P> {
	let z_packed_y = assert_z_packed_bytes_norm_for::<P>(builder, z_words);
	let w1_prime = use_hint_for::<P>(builder, h_coeffs, w_approx_coeffs);
	let sample_in_ball = one_block_w1encode_hash_sparse_relation_unchecked_w1_for::<P>(
		builder, c_tilde, mu_words, &w1_prime,
	);

	MldsaBitHeavyCircuit {
		z_packed_y,
		w1_prime,
		sample_in_ball,
		_params: std::marker::PhantomData,
	}
}

/// Full ML-DSA bit-heavy relation with canonical compressed hint decoding.
pub fn full_bit_heavy_one_block_canonical_hint_relation_for<P: MldsaParams>(
	builder: &CircuitBuilder,
	c_tilde: &[Wire],
	z_words: &[Wire],
	h_words: &[Wire],
	mu_words: &[Wire],
	w_approx_coeffs: &[Wire],
) -> MldsaBitHeavyCircuit<P> {
	let z_packed_y = assert_z_packed_bytes_norm_for::<P>(builder, z_words);
	let h_coeffs = decode_hint_canonical_for::<P>(builder, h_words);
	let w1_prime = use_hint_checked_for::<P>(builder, &h_coeffs, w_approx_coeffs);
	let sample_in_ball = one_block_w1encode_hash_sparse_relation_unchecked_w1_for::<P>(
		builder, c_tilde, mu_words, &w1_prime,
	);

	MldsaBitHeavyCircuit {
		z_packed_y,
		w1_prime,
		sample_in_ball,
		_params: std::marker::PhantomData,
	}
}

/// Full ML-DSA bit-heavy relation with canonical compressed hint bytes and redundant expanded
/// hidden hint bits.
pub fn full_bit_heavy_one_block_canonical_hint_matched_relation_for<P: MldsaParams>(
	builder: &CircuitBuilder,
	c_tilde: &[Wire],
	z_words: &[Wire],
	h_words: &[Wire],
	h_coeffs: &[Wire],
	mu_words: &[Wire],
	w_approx_coeffs: &[Wire],
) -> MldsaBitHeavyCircuit<P> {
	let z_packed_y = assert_z_packed_bytes_norm_for::<P>(builder, z_words);
	assert_hint_canonical_matches_expanded_for::<P>(builder, h_words, h_coeffs);
	let w1_prime = use_hint_checked_for::<P>(builder, h_coeffs, w_approx_coeffs);
	let sample_in_ball = one_block_w1encode_hash_sparse_relation_unchecked_w1_for::<P>(
		builder, c_tilde, mu_words, &w1_prime,
	);

	MldsaBitHeavyCircuit {
		z_packed_y,
		w1_prime,
		sample_in_ball,
		_params: std::marker::PhantomData,
	}
}
