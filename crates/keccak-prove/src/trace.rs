// Copyright 2026 The Binius Developers

//! Native Keccak-f[1600] trace construction.

use std::array;

use crate::{
	constants::{N_LANES, N_ROUNDS},
	unrolled,
};

#[cfg(test)]
use crate::constants::{RHO_OFFSETS, ROUND_CONSTANTS, lane};

/// One Keccak-f[1600] state represented as 25 little-endian 64-bit lanes.
pub type State = [u64; N_LANES];

/// Round-local native values needed by the Keccak-specific prover.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct RoundTrace {
	/// Input state at the start of the round.
	pub input: State,
	/// Pre-chi state after theta, rho, and pi.
	pub pre_chi: State,
	/// Output state after chi and iota.
	pub output: State,
}

/// Full 24-round Keccak-f[1600] trace.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct PermutationTrace {
	/// Per-round traces in forward order.
	pub rounds: [RoundTrace; N_ROUNDS],
}

impl PermutationTrace {
	/// Build a full trace from one input state.
	pub fn new(input: State) -> Self {
		let mut state = input;
		let rounds = array::from_fn(|round| {
			let input = state;
			let pre_chi = theta_rho_pi(input);
			let output = chi_iota(pre_chi, round);
			state = output;
			RoundTrace {
				input,
				pre_chi,
				output,
			}
		});

		Self { rounds }
	}

	/// Return the final permutation output.
	pub fn output(&self) -> State {
		self.rounds[N_ROUNDS - 1].output
	}
}

/// Apply one full Keccak-f[1600] permutation.
pub fn permutation(input: State) -> State {
	PermutationTrace::new(input).output()
}

/// Apply one full Keccak round.
pub fn round(input: State, round: usize) -> State {
	chi_iota(theta_rho_pi(input), round)
}

/// Apply theta followed by rho+pi, returning the pre-chi lanes.
pub fn theta_rho_pi(input: State) -> State {
	unrolled::theta_rho_pi(input)
}

/// Apply chi followed by iota to a pre-chi state.
pub fn chi_iota(pre_chi: State, round: usize) -> State {
	unrolled::chi_iota(pre_chi, round)
}

#[cfg(test)]
fn theta_rho_pi_reference(input: State) -> State {
	let mut state = input;
	theta_reference(&mut state);
	rho_pi_reference(&mut state);
	state
}

#[cfg(test)]
fn chi_iota_reference(mut pre_chi: State, round: usize) -> State {
	chi(&mut pre_chi);
	pre_chi[0] ^= ROUND_CONSTANTS[round];
	pre_chi
}

#[cfg(test)]
fn theta_reference(state: &mut State) {
	let c = array::from_fn::<_, 5, _>(|x| {
		state[lane(x, 0)]
			^ state[lane(x, 1)]
			^ state[lane(x, 2)]
			^ state[lane(x, 3)]
			^ state[lane(x, 4)]
	});

	let d = [
		c[4] ^ c[1].rotate_left(1),
		c[0] ^ c[2].rotate_left(1),
		c[1] ^ c[3].rotate_left(1),
		c[2] ^ c[4].rotate_left(1),
		c[3] ^ c[0].rotate_left(1),
	];

	for y in 0..5 {
		for x in 0..5 {
			state[lane(x, y)] ^= d[x];
		}
	}
}

#[cfg(test)]
fn rho_pi_reference(state: &mut State) {
	let mut output = [0u64; N_LANES];
	for y in 0..5 {
		for x in 0..5 {
			output[lane(y, (2 * x + 3 * y) % 5)] =
				state[lane(x, y)].rotate_left(RHO_OFFSETS[lane(x, y)]);
		}
	}
	*state = output;
}

#[cfg(test)]
fn chi(state: &mut State) {
	for y in 0..5 {
		let a0 = state[lane(0, y)];
		let a1 = state[lane(1, y)];
		let a2 = state[lane(2, y)];
		let a3 = state[lane(3, y)];
		let a4 = state[lane(4, y)];

		state[lane(0, y)] = a0 ^ ((!a1) & a2);
		state[lane(1, y)] = a1 ^ ((!a2) & a3);
		state[lane(2, y)] = a2 ^ ((!a3) & a4);
		state[lane(3, y)] = a3 ^ ((!a4) & a0);
		state[lane(4, y)] = a4 ^ ((!a0) & a1);
	}
}

#[cfg(test)]
mod tests {
	use rand::{Rng, SeedableRng, rngs::StdRng};

	use super::*;

	#[test]
	fn trace_rounds_chain() {
		let mut rng = StdRng::seed_from_u64(0);
		let input = rng.random::<State>();
		let trace = PermutationTrace::new(input);

		assert_eq!(trace.rounds[0].input, input);
		for round_idx in 0..N_ROUNDS {
			assert_eq!(
				trace.rounds[round_idx].output,
				round(trace.rounds[round_idx].input, round_idx)
			);
			if round_idx + 1 < N_ROUNDS {
				assert_eq!(trace.rounds[round_idx].output, trace.rounds[round_idx + 1].input);
			}
		}
	}

	#[test]
	fn unrolled_round_steps_match_reference() {
		let mut rng = StdRng::seed_from_u64(2);

		for _ in 0..32 {
			let input = rng.random::<State>();
			let pre_chi = theta_rho_pi(input);

			assert_eq!(pre_chi, theta_rho_pi_reference(input));
			for round_idx in 0..N_ROUNDS {
				assert_eq!(chi_iota(pre_chi, round_idx), chi_iota_reference(pre_chi, round_idx));
				assert_eq!(round(input, round_idx), chi_iota_reference(pre_chi, round_idx));
			}
		}
	}
}
