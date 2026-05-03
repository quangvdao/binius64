// Copyright 2026 The Binius Developers

//! Unrolled Keccak-f[1600] hot-path helpers.
//!
//! These functions intentionally avoid lane-coordinate loops, `% 5`, and `/ 5`
//! in the row-generation path used by the prover.

use crate::{constants::ROUND_CONSTANTS, trace::State};

/// Apply theta followed by rho+pi, returning the pre-chi lanes.
#[inline(always)]
pub fn theta_rho_pi(state: State) -> State {
	let a00 = state[0];
	let a10 = state[1];
	let a20 = state[2];
	let a30 = state[3];
	let a40 = state[4];
	let a01 = state[5];
	let a11 = state[6];
	let a21 = state[7];
	let a31 = state[8];
	let a41 = state[9];
	let a02 = state[10];
	let a12 = state[11];
	let a22 = state[12];
	let a32 = state[13];
	let a42 = state[14];
	let a03 = state[15];
	let a13 = state[16];
	let a23 = state[17];
	let a33 = state[18];
	let a43 = state[19];
	let a04 = state[20];
	let a14 = state[21];
	let a24 = state[22];
	let a34 = state[23];
	let a44 = state[24];

	let c0 = a00 ^ a01 ^ a02 ^ a03 ^ a04;
	let c1 = a10 ^ a11 ^ a12 ^ a13 ^ a14;
	let c2 = a20 ^ a21 ^ a22 ^ a23 ^ a24;
	let c3 = a30 ^ a31 ^ a32 ^ a33 ^ a34;
	let c4 = a40 ^ a41 ^ a42 ^ a43 ^ a44;

	let d0 = c4 ^ c1.rotate_left(1);
	let d1 = c0 ^ c2.rotate_left(1);
	let d2 = c1 ^ c3.rotate_left(1);
	let d3 = c2 ^ c4.rotate_left(1);
	let d4 = c3 ^ c0.rotate_left(1);

	let b00 = a00 ^ d0;
	let b10 = a10 ^ d1;
	let b20 = a20 ^ d2;
	let b30 = a30 ^ d3;
	let b40 = a40 ^ d4;
	let b01 = a01 ^ d0;
	let b11 = a11 ^ d1;
	let b21 = a21 ^ d2;
	let b31 = a31 ^ d3;
	let b41 = a41 ^ d4;
	let b02 = a02 ^ d0;
	let b12 = a12 ^ d1;
	let b22 = a22 ^ d2;
	let b32 = a32 ^ d3;
	let b42 = a42 ^ d4;
	let b03 = a03 ^ d0;
	let b13 = a13 ^ d1;
	let b23 = a23 ^ d2;
	let b33 = a33 ^ d3;
	let b43 = a43 ^ d4;
	let b04 = a04 ^ d0;
	let b14 = a14 ^ d1;
	let b24 = a24 ^ d2;
	let b34 = a34 ^ d3;
	let b44 = a44 ^ d4;

	[
		b00,
		b11.rotate_left(44),
		b22.rotate_left(43),
		b33.rotate_left(21),
		b44.rotate_left(14),
		b30.rotate_left(28),
		b41.rotate_left(20),
		b02.rotate_left(3),
		b13.rotate_left(45),
		b24.rotate_left(61),
		b10.rotate_left(1),
		b21.rotate_left(6),
		b32.rotate_left(25),
		b43.rotate_left(8),
		b04.rotate_left(18),
		b40.rotate_left(27),
		b01.rotate_left(36),
		b12.rotate_left(10),
		b23.rotate_left(15),
		b34.rotate_left(56),
		b20.rotate_left(62),
		b31.rotate_left(55),
		b42.rotate_left(39),
		b03.rotate_left(41),
		b14.rotate_left(2),
	]
}

/// Apply chi followed by iota to a pre-chi state.
#[inline(always)]
pub fn chi_iota(pre_chi: State, round: usize) -> State {
	let p00 = pre_chi[0];
	let p10 = pre_chi[1];
	let p20 = pre_chi[2];
	let p30 = pre_chi[3];
	let p40 = pre_chi[4];
	let p01 = pre_chi[5];
	let p11 = pre_chi[6];
	let p21 = pre_chi[7];
	let p31 = pre_chi[8];
	let p41 = pre_chi[9];
	let p02 = pre_chi[10];
	let p12 = pre_chi[11];
	let p22 = pre_chi[12];
	let p32 = pre_chi[13];
	let p42 = pre_chi[14];
	let p03 = pre_chi[15];
	let p13 = pre_chi[16];
	let p23 = pre_chi[17];
	let p33 = pre_chi[18];
	let p43 = pre_chi[19];
	let p04 = pre_chi[20];
	let p14 = pre_chi[21];
	let p24 = pre_chi[22];
	let p34 = pre_chi[23];
	let p44 = pre_chi[24];

	[
		p00 ^ ((!p10) & p20) ^ ROUND_CONSTANTS[round],
		p10 ^ ((!p20) & p30),
		p20 ^ ((!p30) & p40),
		p30 ^ ((!p40) & p00),
		p40 ^ ((!p00) & p10),
		p01 ^ ((!p11) & p21),
		p11 ^ ((!p21) & p31),
		p21 ^ ((!p31) & p41),
		p31 ^ ((!p41) & p01),
		p41 ^ ((!p01) & p11),
		p02 ^ ((!p12) & p22),
		p12 ^ ((!p22) & p32),
		p22 ^ ((!p32) & p42),
		p32 ^ ((!p42) & p02),
		p42 ^ ((!p02) & p12),
		p03 ^ ((!p13) & p23),
		p13 ^ ((!p23) & p33),
		p23 ^ ((!p33) & p43),
		p33 ^ ((!p43) & p03),
		p43 ^ ((!p03) & p13),
		p04 ^ ((!p14) & p24),
		p14 ^ ((!p24) & p34),
		p24 ^ ((!p34) & p44),
		p34 ^ ((!p44) & p04),
		p44 ^ ((!p04) & p14),
	]
}

/// Compute all 25 BitAnd-style chi+iota residual words.
#[inline(always)]
pub fn chi_iota_residual_words(pre_chi: &State, next: &State, round: usize) -> State {
	let p00 = pre_chi[0];
	let p10 = pre_chi[1];
	let p20 = pre_chi[2];
	let p30 = pre_chi[3];
	let p40 = pre_chi[4];
	let p01 = pre_chi[5];
	let p11 = pre_chi[6];
	let p21 = pre_chi[7];
	let p31 = pre_chi[8];
	let p41 = pre_chi[9];
	let p02 = pre_chi[10];
	let p12 = pre_chi[11];
	let p22 = pre_chi[12];
	let p32 = pre_chi[13];
	let p42 = pre_chi[14];
	let p03 = pre_chi[15];
	let p13 = pre_chi[16];
	let p23 = pre_chi[17];
	let p33 = pre_chi[18];
	let p43 = pre_chi[19];
	let p04 = pre_chi[20];
	let p14 = pre_chi[21];
	let p24 = pre_chi[22];
	let p34 = pre_chi[23];
	let p44 = pre_chi[24];

	[
		((!p10) & p20) ^ p00 ^ next[0] ^ ROUND_CONSTANTS[round],
		((!p20) & p30) ^ p10 ^ next[1],
		((!p30) & p40) ^ p20 ^ next[2],
		((!p40) & p00) ^ p30 ^ next[3],
		((!p00) & p10) ^ p40 ^ next[4],
		((!p11) & p21) ^ p01 ^ next[5],
		((!p21) & p31) ^ p11 ^ next[6],
		((!p31) & p41) ^ p21 ^ next[7],
		((!p41) & p01) ^ p31 ^ next[8],
		((!p01) & p11) ^ p41 ^ next[9],
		((!p12) & p22) ^ p02 ^ next[10],
		((!p22) & p32) ^ p12 ^ next[11],
		((!p32) & p42) ^ p22 ^ next[12],
		((!p42) & p02) ^ p32 ^ next[13],
		((!p02) & p12) ^ p42 ^ next[14],
		((!p13) & p23) ^ p03 ^ next[15],
		((!p23) & p33) ^ p13 ^ next[16],
		((!p33) & p43) ^ p23 ^ next[17],
		((!p43) & p03) ^ p33 ^ next[18],
		((!p03) & p13) ^ p43 ^ next[19],
		((!p14) & p24) ^ p04 ^ next[20],
		((!p24) & p34) ^ p14 ^ next[21],
		((!p34) & p44) ^ p24 ^ next[22],
		((!p44) & p04) ^ p34 ^ next[23],
		((!p04) & p14) ^ p44 ^ next[24],
	]
}
