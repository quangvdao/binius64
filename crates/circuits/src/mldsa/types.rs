// Copyright 2025 Irreducible Inc.

use binius_frontend::Wire;

use super::mldsa44;

/// ML-DSA-44 one-block `SampleInBall` relation output plus private trace wires.
pub struct Mldsa44SampleInBallOneBlock {
	/// Challenge polynomial coefficients encoded as unsigned 64-bit words: `0`, `1`, or `u64::MAX`
	/// for `-1`.
	pub coeffs: [Wire; mldsa44::N],
	/// Private bounded rejection-loop counts. Entry `r` is the number of one-byte draws consumed
	/// for the `r`th nonzero challenge coefficient.
	pub draw_counts: [Wire; mldsa44::TAU],
}

/// Sparse ML-DSA-44 one-block `SampleInBall` relation output plus private trace wires.
pub struct Mldsa44SampleInBallOneBlockSparse {
	/// Sparse challenge positions in update order.
	pub positions: [Wire; mldsa44::TAU],
	/// Sparse challenge signs encoded as unsigned 64-bit words: `1` or `u64::MAX` for `-1`.
	pub signs: [Wire; mldsa44::TAU],
	/// Private bounded rejection-loop counts. Entry `r` is the number of one-byte draws consumed
	/// for the `r`th nonzero challenge coefficient.
	pub draw_counts: [Wire; mldsa44::TAU],
}

/// Top-level ML-DSA-44 bit-heavy circuit outputs.
pub struct Mldsa44BitHeavyCircuit {
	/// Decoded packed `y = gamma1 - z` coefficients from hidden signature `z`.
	pub z_packed_y: Vec<Wire>,
	/// Hidden `w1_prime = UseHint(h, wApprox)` coefficients.
	pub w1_prime: Vec<Wire>,
	/// Hidden sparse challenge polynomial plus private sampler trace.
	pub sample_in_ball: Mldsa44SampleInBallOneBlockSparse,
}
