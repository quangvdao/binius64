// Copyright 2026 The Binius Developers

//! Keccak-f[1600] constants and lane helpers.

/// Number of 64-bit lanes in a Keccak-f[1600] state.
pub const N_LANES: usize = 25;

/// Number of bits in each Keccak lane.
pub const LANE_BITS: usize = 64;

/// Number of Keccak-f[1600] rounds.
pub const N_ROUNDS: usize = 24;

/// Iota round constants.
pub const ROUND_CONSTANTS: [u64; N_ROUNDS] = [
	0x0000_0000_0000_0001,
	0x0000_0000_0000_8082,
	0x8000_0000_0000_808A,
	0x8000_0000_8000_8000,
	0x0000_0000_0000_808B,
	0x0000_0000_8000_0001,
	0x8000_0000_8000_8081,
	0x8000_0000_0000_8009,
	0x0000_0000_0000_008A,
	0x0000_0000_0000_0088,
	0x0000_0000_8000_8009,
	0x0000_0000_8000_000A,
	0x0000_0000_8000_808B,
	0x8000_0000_0000_008B,
	0x8000_0000_0000_8089,
	0x8000_0000_0000_8003,
	0x8000_0000_0000_8002,
	0x8000_0000_0000_0080,
	0x0000_0000_0000_800A,
	0x8000_0000_8000_000A,
	0x8000_0000_8000_8081,
	0x8000_0000_0000_8080,
	0x0000_0000_8000_0001,
	0x8000_0000_8000_8008,
];

/// Rho rotation offsets in lane order `x + 5*y`.
#[rustfmt::skip]
pub const RHO_OFFSETS: [u32; N_LANES] = [
	 0,  1, 62, 28, 27,
	36, 44,  6, 55, 20,
	 3, 10, 43, 25, 39,
	41, 45, 15, 21,  8,
	18,  2, 61, 56, 14,
];

/// Return the lane index for coordinates `(x, y)`.
#[inline(always)]
pub const fn lane(x: usize, y: usize) -> usize {
	x + 5 * y
}

/// Return the x-coordinate of a lane index.
#[inline(always)]
pub const fn lane_x(lane: usize) -> usize {
	lane % 5
}

/// Return the y-coordinate of a lane index.
#[inline(always)]
pub const fn lane_y(lane: usize) -> usize {
	lane / 5
}
