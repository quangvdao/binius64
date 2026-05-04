// Copyright 2025 Irreducible Inc.

use binius_frontend::Wire;

use super::MldsaParams;

/// ML-DSA one-block `SampleInBall` relation output plus private trace wires.
pub struct MldsaSampleInBallOneBlock {
	/// Challenge polynomial coefficients encoded as unsigned 64-bit words: `0`, `1`, or `u64::MAX`
	/// for `-1`.
	pub coeffs: Vec<Wire>,
	/// Private bounded rejection-loop counts. Entry `r` is the number of one-byte draws consumed
	/// for the `r`th nonzero challenge coefficient.
	pub draw_counts: Vec<Wire>,
}

/// Sparse ML-DSA one-block `SampleInBall` relation output plus private trace wires.
pub struct MldsaSampleInBallOneBlockSparse {
	/// Sparse challenge positions in update order.
	pub positions: Vec<Wire>,
	/// Sparse challenge signs encoded as unsigned 64-bit words: `1` or `u64::MAX` for `-1`.
	pub signs: Vec<Wire>,
	/// Private bounded rejection-loop counts. Entry `r` is the number of one-byte draws consumed
	/// for the `r`th nonzero challenge coefficient.
	pub draw_counts: Vec<Wire>,
}

/// Top-level ML-DSA bit-heavy circuit outputs.
pub struct MldsaBitHeavyCircuit<P: MldsaParams> {
	/// Decoded packed `y = gamma1 - z` coefficients from hidden signature `z`.
	pub z_packed_y: Vec<Wire>,
	/// Hidden `w1_prime = UseHint(h, wApprox)` coefficients.
	pub w1_prime: Vec<Wire>,
	/// Hidden sparse challenge polynomial plus private sampler trace.
	pub sample_in_ball: MldsaSampleInBallOneBlockSparse,
	pub(crate) _params: std::marker::PhantomData<P>,
}
