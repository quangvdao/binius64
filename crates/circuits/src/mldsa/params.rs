// Copyright 2025 Irreducible Inc.

/// Shared parameter interface for ML-DSA bit-heavy circuit shapes.
///
/// The current implementation exposes ML-DSA-44 concrete helpers, but keeping the
/// variant constants behind this trait makes the future ML-DSA-65/87 cutover less invasive.
pub trait MldsaParams {
	const N: usize;
	const Q: u64;
	const K: usize;
	const L: usize;
	const TAU: usize;
	const BETA: u64;
	const OMEGA: u64;
	const GAMMA1: u64;
	const GAMMA2: u64;
	const C_TILDE_BYTES: usize;
	const MU_BYTES: usize;
	const POLY_W1_PACKED_BYTES: usize;
	const W1_BITS_PER_COEFF: usize;
	const W1_COEFF_MAX: u64;
	const POLY_Z_PACKED_BYTES: usize;
	const Z_BITS_PER_COEFF: usize;
	const SAMPLE_IN_BALL_SIGN_BYTES: usize;
	const SAMPLE_IN_BALL_ONE_BLOCK_DRAW_CAP_BYTES: usize;
}

/// ML-DSA-44 constants for the first prototype target.
pub mod mldsa44 {
	pub const N: usize = 256;
	pub const Q: u64 = 8_380_417;
	/// Number of rows in the public matrix.
	pub const K: usize = 4;
	/// Number of columns in the public matrix.
	pub const L: usize = 4;
	pub const TAU: usize = 39;
	pub const BETA: u64 = 78;
	pub const OMEGA: u64 = 80;
	pub const OMEGA_USIZE: usize = OMEGA as usize;
	pub const GAMMA1: u64 = 1 << 17;
	pub const GAMMA2: u64 = (Q - 1) / 88;
	pub const TWO_GAMMA2: u64 = 2 * GAMMA2;
	pub const C_TILDE_BYTES: usize = 32;
	pub const MU_BYTES: usize = 64;
	pub const MU_WORDS: usize = MU_BYTES / 8;
	pub const POLY_W1_PACKED_BYTES: usize = 192;
	pub const W1_ENCODE_BYTES: usize = K * POLY_W1_PACKED_BYTES;
	pub const W1_ENCODE_WORDS: usize = W1_ENCODE_BYTES / 8;
	pub const W1_COEFFICIENTS: usize = K * N;
	pub const W1_BITS_PER_COEFF: usize = 6;
	pub const W1_M: u64 = (Q - 1) / TWO_GAMMA2;
	pub const W1_COEFF_MAX: u64 = 43;
	pub const HINT_BYTES: usize = OMEGA_USIZE + K;
	pub const HINT_WORDS: usize = HINT_BYTES.div_ceil(8);
	pub const FINAL_CHALLENGE_INPUT_BYTES: usize = MU_BYTES + W1_ENCODE_BYTES;
	pub const FINAL_CHALLENGE_INPUT_WORDS: usize = FINAL_CHALLENGE_INPUT_BYTES / 8;
	pub const FINAL_CHALLENGE_KECCAK_F_CALLS: usize = 7;

	/// Packed `y = gamma1 - z` range equivalent to `|z| < gamma1 - beta`.
	pub const Z_NORM_PACKED_Y_MIN: u64 = BETA + 1;
	pub const Z_NORM_PACKED_Y_MAX: u64 = 2 * GAMMA1 - BETA - 1;
	pub const Z_COEFFICIENTS: usize = L * 256;
	pub const Z_BITS_PER_COEFF: usize = 18;
	pub const POLY_Z_PACKED_BYTES: usize = 576;
	pub const Z_PACKED_BYTES: usize = L * POLY_Z_PACKED_BYTES;
	pub const Z_PACKED_WORDS: usize = Z_PACKED_BYTES / 8;

	pub const SAMPLE_IN_BALL_SIGN_BYTES: usize = 8;
	pub const SAMPLE_IN_BALL_ONE_BLOCK_DRAW_CAP_BYTES: usize = 128;
	pub const SAMPLE_IN_BALL_ONE_BLOCK_STREAM_BYTES: usize =
		SAMPLE_IN_BALL_SIGN_BYTES + SAMPLE_IN_BALL_ONE_BLOCK_DRAW_CAP_BYTES;
	pub const SAMPLE_IN_BALL_ONE_BLOCK_KECCAK_F_CALLS: usize = 1;
}

/// Marker type for the ML-DSA-44 parameter set.
pub struct Mldsa44;

impl MldsaParams for Mldsa44 {
	const N: usize = mldsa44::N;
	const Q: u64 = mldsa44::Q;
	const K: usize = mldsa44::K;
	const L: usize = mldsa44::L;
	const TAU: usize = mldsa44::TAU;
	const BETA: u64 = mldsa44::BETA;
	const OMEGA: u64 = mldsa44::OMEGA;
	const GAMMA1: u64 = mldsa44::GAMMA1;
	const GAMMA2: u64 = mldsa44::GAMMA2;
	const C_TILDE_BYTES: usize = mldsa44::C_TILDE_BYTES;
	const MU_BYTES: usize = mldsa44::MU_BYTES;
	const POLY_W1_PACKED_BYTES: usize = mldsa44::POLY_W1_PACKED_BYTES;
	const W1_BITS_PER_COEFF: usize = mldsa44::W1_BITS_PER_COEFF;
	const W1_COEFF_MAX: u64 = mldsa44::W1_COEFF_MAX;
	const POLY_Z_PACKED_BYTES: usize = mldsa44::POLY_Z_PACKED_BYTES;
	const Z_BITS_PER_COEFF: usize = mldsa44::Z_BITS_PER_COEFF;
	const SAMPLE_IN_BALL_SIGN_BYTES: usize = mldsa44::SAMPLE_IN_BALL_SIGN_BYTES;
	const SAMPLE_IN_BALL_ONE_BLOCK_DRAW_CAP_BYTES: usize =
		mldsa44::SAMPLE_IN_BALL_ONE_BLOCK_DRAW_CAP_BYTES;
}
