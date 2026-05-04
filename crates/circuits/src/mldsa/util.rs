// Copyright 2025 Irreducible Inc.

use binius_frontend::{CircuitBuilder, Wire};

use crate::multiplexer::single_wire_multiplex;

/// Conditional `assert_true`: asserts `pred` only when `cond` is the all-ones MSB-true wire.
///
/// Implemented by selecting between `pred` and a hard-coded "true" wire (`u64::MAX`) so that the
/// assertion gate sees an MSB-true value whenever `cond` is false. This trick keeps the assertion
/// graph fixed (the gate is always added) while letting the prover skip the actual check on
/// inactive iterations of a loop. Both branches must always be honestly populated; otherwise a
/// malicious prover could disable a real check by toggling `cond`. Callers must therefore drive
/// `cond` from a constraint already enforced elsewhere.
pub(crate) fn assert_true_cond(
	builder: &CircuitBuilder,
	name: impl Into<String>,
	pred: Wire,
	cond: Wire,
) {
	let true_wire = builder.add_constant_64(u64::MAX);
	let guarded_pred = builder.select(cond, pred, true_wire);
	builder.assert_true(name, guarded_pred);
}

pub(crate) fn iadd_wrapping(builder: &CircuitBuilder, a: Wire, b: Wire) -> Wire {
	let (sum, _carry) = builder.iadd(a, b);
	sum
}

pub(crate) fn isub_one(builder: &CircuitBuilder, x: Wire) -> Wire {
	let one = builder.add_constant_64(1);
	let zero = builder.add_constant_64(0);
	let (diff, _borrow) = builder.isub_bin_bout(x, one, zero);
	diff
}

/// Multiplies a wire by a constant using a shift-and-add ladder, instead of an `imul` gate.
///
/// **Correctness contract**: the caller must guarantee the product fits in 64 bits, i.e.
/// `x * c < 2^64`. Each `iadd` discards its high carry, so an overflowing multiplication is
/// silently truncated and would corrupt downstream constraints. Used by `high_bits_for` where
/// `x <= q-1 < 2^23` and `c = 11275 < 2^14`, leaving 27 bits of slack. Trades one `MUL` opcode for
/// `popcount(c)` shifts and `popcount(c) - 1` adds, which is cheaper for low-Hamming-weight
/// constants under the AND-reduction proof system.
pub(crate) fn mul_const_no_overflow(builder: &CircuitBuilder, x: Wire, c: u64) -> Wire {
	let mut acc = None;
	for bit in 0..u64::BITS {
		if (c >> bit) & 1 == 0 {
			continue;
		}

		let term = if bit == 0 { x } else { builder.shl(x, bit) };
		acc = Some(match acc {
			Some(acc) => builder.iadd(acc, term).0,
			None => term,
		});
	}

	acc.unwrap_or_else(|| builder.add_constant_64(0))
}

/// Balanced multiplexer for indices whose range has already been proved by surrounding constraints.
pub(crate) fn select_indexed_wire_unchecked(
	builder: &CircuitBuilder,
	values: &[Wire],
	index: Wire,
) -> Wire {
	single_wire_multiplex(builder, values, index)
}

pub(crate) fn unpack_bytes_from_words(
	builder: &CircuitBuilder,
	words: &[Wire],
	len_bytes: usize,
) -> Vec<Wire> {
	assert_eq!(words.len(), len_bytes.div_ceil(8), "word count must match byte length",);

	(0..len_bytes)
		.map(|byte_idx| builder.extract_byte(words[byte_idx / 8], (byte_idx % 8) as u32))
		.collect()
}
