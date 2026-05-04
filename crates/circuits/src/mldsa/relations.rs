// Copyright 2025 Irreducible Inc.

use binius_frontend::{CircuitBuilder, Wire};

use super::{
	hashing::mldsa44_final_challenge_hash,
	hint::{assert_mldsa44_hint_canonical_matches_expanded, mldsa44_decode_hint_canonical},
	mldsa44,
	packing::{assert_mldsa44_z_packed_bytes_norm, mldsa44_encode_w1, mldsa44_encode_w1_unchecked},
	sample_in_ball::{
		mldsa44_expand_sparse_sample_in_ball, mldsa44_sample_in_ball_one_block_sparse,
	},
	types::{
		Mldsa44BitHeavyCircuit, Mldsa44SampleInBallOneBlock, Mldsa44SampleInBallOneBlockSparse,
	},
	use_hint::{mldsa44_use_hint, mldsa44_use_hint_checked},
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
pub fn mldsa44_one_block_hidden_hash_relation(
	builder: &CircuitBuilder,
	c_tilde: &[Wire],
	mu_and_w1_bytes: &[Wire],
) -> Mldsa44SampleInBallOneBlock {
	let sparse = mldsa44_one_block_hidden_hash_sparse_relation(builder, c_tilde, mu_and_w1_bytes);
	let coeffs = mldsa44_expand_sparse_sample_in_ball(builder, &sparse.positions, &sparse.signs);
	Mldsa44SampleInBallOneBlock {
		coeffs,
		draw_counts: sparse.draw_counts,
	}
}

fn mldsa44_one_block_hidden_hash_sparse_relation(
	builder: &CircuitBuilder,
	c_tilde: &[Wire],
	mu_and_w1_bytes: &[Wire],
) -> Mldsa44SampleInBallOneBlockSparse {
	let c_tilde_prime = mldsa44_final_challenge_hash(builder, mu_and_w1_bytes);
	for (i, (&expected, &actual)) in c_tilde.iter().zip(c_tilde_prime.iter()).enumerate() {
		builder.assert_eq(format!("mldsa44_c_tilde_final_hash[{i}]"), expected, actual);
	}

	mldsa44_sample_in_ball_one_block_sparse(builder, c_tilde)
}

/// Binds hidden `c_tilde` to `SampleInBall` and to a final hash whose `w1` bytes are encoded in
/// circuit from hidden `w1_prime` coefficient wires.
///
/// This keeps `UseHint` mocked by taking `w1_prime` directly as witness wires, but removes the
/// earlier shortcut where prepacked `w1` bytes were supplied by the witness.
pub fn mldsa44_one_block_w1encode_hash_relation(
	builder: &CircuitBuilder,
	c_tilde: &[Wire],
	mu_words: &[Wire],
	w1_coeffs: &[Wire],
) -> Mldsa44SampleInBallOneBlock {
	let sparse =
		mldsa44_one_block_w1encode_hash_sparse_relation(builder, c_tilde, mu_words, w1_coeffs);
	let coeffs = mldsa44_expand_sparse_sample_in_ball(builder, &sparse.positions, &sparse.signs);
	Mldsa44SampleInBallOneBlock {
		coeffs,
		draw_counts: sparse.draw_counts,
	}
}

fn mldsa44_one_block_w1encode_hash_sparse_relation(
	builder: &CircuitBuilder,
	c_tilde: &[Wire],
	mu_words: &[Wire],
	w1_coeffs: &[Wire],
) -> Mldsa44SampleInBallOneBlockSparse {
	mldsa44_one_block_w1encode_hash_sparse_relation_impl(
		builder, c_tilde, mu_words, w1_coeffs, true,
	)
}

fn mldsa44_one_block_w1encode_hash_sparse_relation_unchecked_w1(
	builder: &CircuitBuilder,
	c_tilde: &[Wire],
	mu_words: &[Wire],
	w1_coeffs: &[Wire],
) -> Mldsa44SampleInBallOneBlockSparse {
	mldsa44_one_block_w1encode_hash_sparse_relation_impl(
		builder, c_tilde, mu_words, w1_coeffs, false,
	)
}

fn mldsa44_one_block_w1encode_hash_sparse_relation_impl(
	builder: &CircuitBuilder,
	c_tilde: &[Wire],
	mu_words: &[Wire],
	w1_coeffs: &[Wire],
	check_w1_range: bool,
) -> Mldsa44SampleInBallOneBlockSparse {
	assert_eq!(
		mu_words.len(),
		mldsa44::MU_WORDS,
		"ML-DSA-44 mu expects 64 bytes packed into 8 words",
	);

	let w1_words = if check_w1_range {
		mldsa44_encode_w1(builder, w1_coeffs)
	} else {
		mldsa44_encode_w1_unchecked(builder, w1_coeffs)
	};
	let mut hash_input = Vec::with_capacity(mldsa44::FINAL_CHALLENGE_INPUT_WORDS);
	hash_input.extend_from_slice(mu_words);
	hash_input.extend_from_slice(&w1_words);

	mldsa44_one_block_hidden_hash_sparse_relation(builder, c_tilde, &hash_input)
}

/// Binds hidden `c_tilde` to `SampleInBall` and to a final hash whose `w1` bytes are derived from
/// hidden `h` and `wApprox` through `UseHint`.
pub fn mldsa44_one_block_use_hint_hash_relation(
	builder: &CircuitBuilder,
	c_tilde: &[Wire],
	mu_words: &[Wire],
	h_coeffs: &[Wire],
	w_approx_coeffs: &[Wire],
) -> Mldsa44SampleInBallOneBlock {
	let w1_coeffs = mldsa44_use_hint(builder, h_coeffs, w_approx_coeffs);
	mldsa44_one_block_w1encode_hash_relation(builder, c_tilde, mu_words, &w1_coeffs)
}

/// Like [`mldsa44_one_block_use_hint_hash_relation`], but decodes canonical compressed hint bytes
/// instead of taking expanded hint bits as witness.
pub fn mldsa44_one_block_canonical_hint_hash_relation(
	builder: &CircuitBuilder,
	c_tilde: &[Wire],
	mu_words: &[Wire],
	h_words: &[Wire],
	w_approx_coeffs: &[Wire],
) -> Mldsa44SampleInBallOneBlock {
	let h_coeffs = mldsa44_decode_hint_canonical(builder, h_words);
	let w1_coeffs = mldsa44_use_hint_checked(builder, &h_coeffs, w_approx_coeffs);
	mldsa44_one_block_w1encode_hash_relation(builder, c_tilde, mu_words, &w1_coeffs)
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
pub fn mldsa44_full_bit_heavy_one_block_relation(
	builder: &CircuitBuilder,
	c_tilde: &[Wire],
	z_words: &[Wire],
	h_coeffs: &[Wire],
	mu_words: &[Wire],
	w_approx_coeffs: &[Wire],
) -> Mldsa44BitHeavyCircuit {
	let z_packed_y = assert_mldsa44_z_packed_bytes_norm(builder, z_words);
	let w1_prime = mldsa44_use_hint(builder, h_coeffs, w_approx_coeffs);
	let sample_in_ball = mldsa44_one_block_w1encode_hash_sparse_relation_unchecked_w1(
		builder, c_tilde, mu_words, &w1_prime,
	);

	Mldsa44BitHeavyCircuit {
		z_packed_y,
		w1_prime,
		sample_in_ball,
	}
}

/// Full ML-DSA-44 bit-heavy relation with canonical compressed hint decoding.
pub fn mldsa44_full_bit_heavy_one_block_canonical_hint_relation(
	builder: &CircuitBuilder,
	c_tilde: &[Wire],
	z_words: &[Wire],
	h_words: &[Wire],
	mu_words: &[Wire],
	w_approx_coeffs: &[Wire],
) -> Mldsa44BitHeavyCircuit {
	let z_packed_y = assert_mldsa44_z_packed_bytes_norm(builder, z_words);
	let h_coeffs = mldsa44_decode_hint_canonical(builder, h_words);
	let w1_prime = mldsa44_use_hint_checked(builder, &h_coeffs, w_approx_coeffs);
	let sample_in_ball = mldsa44_one_block_w1encode_hash_sparse_relation_unchecked_w1(
		builder, c_tilde, mu_words, &w1_prime,
	);

	Mldsa44BitHeavyCircuit {
		z_packed_y,
		w1_prime,
		sample_in_ball,
	}
}

/// Full ML-DSA-44 bit-heavy relation with canonical compressed hint bytes and redundant expanded
/// hidden hint bits.
pub fn mldsa44_full_bit_heavy_one_block_canonical_hint_matched_relation(
	builder: &CircuitBuilder,
	c_tilde: &[Wire],
	z_words: &[Wire],
	h_words: &[Wire],
	h_coeffs: &[Wire],
	mu_words: &[Wire],
	w_approx_coeffs: &[Wire],
) -> Mldsa44BitHeavyCircuit {
	let z_packed_y = assert_mldsa44_z_packed_bytes_norm(builder, z_words);
	assert_mldsa44_hint_canonical_matches_expanded(builder, h_words, h_coeffs);
	let w1_prime = mldsa44_use_hint_checked(builder, h_coeffs, w_approx_coeffs);
	let sample_in_ball = mldsa44_one_block_w1encode_hash_sparse_relation_unchecked_w1(
		builder, c_tilde, mu_words, &w1_prime,
	);

	Mldsa44BitHeavyCircuit {
		z_packed_y,
		w1_prime,
		sample_in_ball,
	}
}
