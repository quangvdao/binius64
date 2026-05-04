// Copyright 2025 Irreducible Inc.

use binius_frontend::{CircuitBuilder, Wire};

use super::{
	MldsaParams,
	hint::assert_h_bits_and_weight_for,
	util::{isub_one, mul_const_no_overflow},
};

/// Computes ML-DSA `HighBits(r)` for `r in Z_q`.
///
/// This is the Barrett-style integer-rounding formula used verbatim by the Dilithium / ML-DSA
/// reference implementation. For `gamma2 = (q - 1) / 88` this is:
///
/// ```text
/// t1 = floor((r + 127) / 128)
/// t2 = t1 * 11275                  // ~ 2^31 / (2 * gamma2)
/// t3 = floor((t2 + 2^23) / 2^24)
/// r1 = (t3 > 43) ? 0 : t3
/// ```
///
/// For `gamma2 = (q - 1) / 32`, the generic implementation uses the corresponding reference
/// constants `1025`, `2^21`, and a low-four-bit mask.
///
/// The `r in [0, q-1]` precondition is enforced in-circuit. The Barrett correctness is sweep-tested
/// in `mldsa44_high_bits_matches_reference_*` for the ML-DSA-44 instantiation.
pub fn high_bits_for<P: MldsaParams>(builder: &CircuitBuilder, r: Wire) -> Wire {
	let q_minus_one = builder.add_constant_64(P::Q - 1);
	builder
		.assert_true(format!("{}_high_bits_r_in_zq", P::label()), builder.icmp_ule(r, q_minus_one));

	let (r_plus_127, _carry) = builder.iadd(r, builder.add_constant_64(127));
	let t = builder.shr(r_plus_127, 7);

	if P::GAMMA2 == (P::Q - 1) / 32 {
		let t = mul_const_no_overflow(builder, t, 1025);
		let (t, _carry) = builder.iadd(t, builder.add_constant_64(1 << 21));
		builder.band(builder.shr(t, 22), builder.add_constant_64(15))
	} else {
		let t = mul_const_no_overflow(builder, t, 11_275);
		let (t, _carry) = builder.iadd(t, builder.add_constant_64(1 << 23));
		let t = builder.shr(t, 24);
		let over_max = builder.icmp_ugt(t, builder.add_constant_64(P::W1_COEFF_MAX));
		builder.select(over_max, builder.add_constant_64(0), t)
	}
}

/// Returns the all-ones MSB-true wire iff the centered remainder `r0 := r - r1 * (2*gamma2) (mod q)`
/// lies strictly in `(0, (q-1)/2]`, i.e. is "positive" in FIPS Decompose's centered representative.
///
/// We deliberately do not reduce modulo `q`: the caller has already proved `r1 = HighBits(r)`, so
/// `r1 * (2*gamma2)` is in `[0, q-1]` and the integer subtraction `r - r1*(2*gamma2)` either gives
/// `r0` (when no borrow) or wraps when `r0` would be negative. The `(no_borrow && nonzero &&
/// in_upper_half)` triple captures exactly the "positive non-zero" branch, which is what `UseHint`
/// uses to decide whether to increment or decrement `r1`.
fn r0_is_positive_for<P: MldsaParams>(builder: &CircuitBuilder, r: Wire, r1: Wire) -> Wire {
	let r1_alpha = mul_const_no_overflow(builder, r1, P::TWO_GAMMA2);
	let zero = builder.add_constant_64(0);
	let (diff, _borrow_bits) = builder.isub_bin_bout(r, r1_alpha, zero);

	let no_borrow = builder.icmp_ule(r1_alpha, r);
	let diff_nonzero = builder.icmp_ne(diff, zero);
	let diff_lte_half = builder.icmp_ule(diff, builder.add_constant_64((P::Q - 1) / 2));
	builder.band(builder.band(no_borrow, diff_nonzero), diff_lte_half)
}

/// Computes ML-DSA `UseHint(h, r)` for one coefficient.
///
/// Per FIPS 204 algorithm `UseHint`, when `h = 0` the result is `r1 = HighBits(r)`. When `h = 1`
/// the result is `(r1 + 1) mod m` if `r0` is positive, else `(r1 - 1) mod m`, where
/// `m = (q - 1) / (2 * gamma2)`. The wrap-around arithmetic is done with explicit selects because
/// the value only ever increments or decrements by one. `h_i in {0,1}` is asserted here for safety;
/// the caller is *also* expected to enforce the global Hamming-weight bound.
pub fn use_hint_coeff_for<P: MldsaParams>(builder: &CircuitBuilder, h: Wire, r: Wire) -> Wire {
	let one = builder.add_constant_64(1);
	builder.assert_true(format!("{}_hint_bit", P::label()), builder.icmp_ule(h, one));
	use_hint_coeff_checked_for::<P>(builder, h, r)
}

pub(crate) fn use_hint_coeff_checked_for<P: MldsaParams>(
	builder: &CircuitBuilder,
	h: Wire,
	r: Wire,
) -> Wire {
	let zero = builder.add_constant_64(0);
	let one = builder.add_constant_64(1);

	let r1 = high_bits_for::<P>(builder, r);
	let r0_positive = r0_is_positive_for::<P>(builder, r, r1);

	let r1_is_max = builder.icmp_eq(r1, builder.add_constant_64(P::W1_COEFF_MAX));
	let (r1_plus_one, _carry) = builder.iadd(r1, one);
	let inc = builder.select(r1_is_max, zero, r1_plus_one);

	let r1_is_zero = builder.icmp_eq(r1, zero);
	let dec_wrapped = isub_one(builder, r1);
	let dec = builder.select(r1_is_zero, builder.add_constant_64(P::W1_COEFF_MAX), dec_wrapped);

	let hinted = builder.select(r0_positive, inc, dec);
	let h_is_one = builder.shl(h, 63);
	builder.select(h_is_one, hinted, r1)
}

/// Computes ML-DSA `UseHint(h, wApprox)` coefficient-wise.
///
/// In addition to the per-coefficient relation, this enforces `h_i in {0, 1}` and the canonical
/// weight bound `sum(h_i) <= omega`. Callers therefore do not need a separate weight gadget.
pub fn use_hint_for<P: MldsaParams>(
	builder: &CircuitBuilder,
	h_coeffs: &[Wire],
	w_approx_coeffs: &[Wire],
) -> Vec<Wire> {
	assert_eq!(h_coeffs.len(), P::W1_COEFFICIENTS, "{} h coefficient count mismatch", P::label(),);
	assert_eq!(
		w_approx_coeffs.len(),
		P::W1_COEFFICIENTS,
		"{} wApprox coefficient count mismatch",
		P::label(),
	);

	assert_h_bits_and_weight_for::<P>(builder, h_coeffs);

	use_hint_checked_for::<P>(builder, h_coeffs, w_approx_coeffs)
}

pub(crate) fn use_hint_checked_for<P: MldsaParams>(
	builder: &CircuitBuilder,
	h_coeffs: &[Wire],
	w_approx_coeffs: &[Wire],
) -> Vec<Wire> {
	assert_eq!(h_coeffs.len(), P::W1_COEFFICIENTS, "{} h coefficient count mismatch", P::label(),);
	assert_eq!(
		w_approx_coeffs.len(),
		P::W1_COEFFICIENTS,
		"{} wApprox coefficient count mismatch",
		P::label(),
	);

	h_coeffs
		.iter()
		.zip(w_approx_coeffs.iter())
		.map(|(&h, &r)| use_hint_coeff_checked_for::<P>(builder, h, r))
		.collect()
}
