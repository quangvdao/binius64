// Copyright 2025 Irreducible Inc.

use binius_core::{verify::verify_constraints, word::Word};
use binius_frontend::{CircuitBuilder, CircuitStat, Wire};
use rand::{RngCore, SeedableRng, rngs::StdRng};
use rstest::rstest;
use sha3::{
	Shake256,
	digest::{ExtendableOutput, Update, XofReader},
};

use super::*;

fn words_from_bytes(bytes: &[u8]) -> Vec<u64> {
	bytes
		.chunks(8)
		.map(|chunk| {
			let mut word_bytes = [0u8; 8];
			word_bytes[..chunk.len()].copy_from_slice(chunk);
			u64::from_le_bytes(word_bytes)
		})
		.collect()
}

fn pack_mldsa44_z_y_coeffs(packed_y_coeff_values: &[u64]) -> Vec<u64> {
	assert_eq!(packed_y_coeff_values.len(), mldsa44::Z_COEFFICIENTS);

	let mut words = vec![0u64; mldsa44::Z_PACKED_WORDS];
	for (coeff_idx, &coeff) in packed_y_coeff_values.iter().enumerate() {
		assert!(coeff < (1 << mldsa44::Z_BITS_PER_COEFF));
		for bit in 0..mldsa44::Z_BITS_PER_COEFF {
			if (coeff >> bit) & 1 == 1 {
				let bit_idx = coeff_idx * mldsa44::Z_BITS_PER_COEFF + bit;
				words[bit_idx / 64] |= 1 << (bit_idx % 64);
			}
		}
	}

	words
}

fn pack_mldsa44_w1_coeffs(w1_coeff_values: &[u64]) -> Vec<u64> {
	assert_eq!(w1_coeff_values.len(), mldsa44::W1_COEFFICIENTS);

	let mut words = vec![0u64; mldsa44::W1_ENCODE_WORDS];
	for (coeff_idx, &coeff) in w1_coeff_values.iter().enumerate() {
		assert!(coeff <= mldsa44::W1_COEFF_MAX);
		for bit in 0..mldsa44::W1_BITS_PER_COEFF {
			if (coeff >> bit) & 1 == 1 {
				let bit_idx = coeff_idx * mldsa44::W1_BITS_PER_COEFF + bit;
				words[bit_idx / 64] |= 1 << (bit_idx % 64);
			}
		}
	}

	words
}

fn pack_mldsa44_hint_bytes(poly_positions: &[Vec<u8>; mldsa44::K]) -> [u8; mldsa44::HINT_BYTES] {
	let mut h = [0u8; mldsa44::HINT_BYTES];
	let mut idx = 0usize;
	for (poly_idx, positions) in poly_positions.iter().enumerate() {
		for &pos in positions {
			h[idx] = pos;
			idx += 1;
		}
		h[mldsa44::OMEGA_USIZE + poly_idx] = idx as u8;
	}
	h
}

fn host_mldsa44_decode_hint(h_bytes: &[u8; mldsa44::HINT_BYTES]) -> Option<Vec<u64>> {
	let mut h = vec![0u64; mldsa44::W1_COEFFICIENTS];
	let mut index = 0usize;
	for poly_idx in 0..mldsa44::K {
		let endpoint = h_bytes[mldsa44::OMEGA_USIZE + poly_idx] as usize;
		if endpoint < index || endpoint > mldsa44::OMEGA_USIZE {
			return None;
		}
		let mut prev = None;
		for &pos in &h_bytes[index..endpoint] {
			if let Some(prev) = prev {
				if pos <= prev {
					return None;
				}
			}
			h[poly_idx * mldsa44::N + pos as usize] = 1;
			prev = Some(pos);
		}
		index = endpoint;
	}
	if h_bytes[index..mldsa44::OMEGA_USIZE].iter().any(|&x| x != 0) {
		return None;
	}
	Some(h)
}

fn host_mldsa44_high_bits(r: u64) -> u64 {
	assert!(r < mldsa44::Q);
	let mut r1 = (r + 127) >> 7;
	r1 = ((r1 * 11_275) + (1 << 23)) >> 24;
	if r1 > mldsa44::W1_COEFF_MAX { 0 } else { r1 }
}

fn host_mldsa44_r0_is_positive(r: u64, r1: u64) -> bool {
	let r1_alpha = r1 * mldsa44::TWO_GAMMA2;
	r >= r1_alpha && r != r1_alpha && r - r1_alpha <= (mldsa44::Q - 1) / 2
}

fn host_mldsa44_use_hint(h: u64, r: u64) -> u64 {
	assert!(h <= 1);
	let r1 = host_mldsa44_high_bits(r);
	if h == 0 {
		return r1;
	}

	if host_mldsa44_r0_is_positive(r, r1) {
		if r1 == mldsa44::W1_COEFF_MAX {
			0
		} else {
			r1 + 1
		}
	} else if r1 == 0 {
		mldsa44::W1_COEFF_MAX
	} else {
		r1 - 1
	}
}

fn verify_mldsa44_z_packed_bytes_norm_witness(packed_y_coeff_values: &[u64]) -> bool {
	let z_word_values = pack_mldsa44_z_y_coeffs(packed_y_coeff_values);
	let builder = CircuitBuilder::new();
	let z_words: Vec<_> = (0..z_word_values.len())
		.map(|_| builder.add_witness())
		.collect();
	assert_mldsa44_z_packed_bytes_norm(&builder, &z_words);

	let circuit = builder.build();
	let cs = circuit.constraint_system();
	let mut witness = circuit.new_witness_filler();
	for (&wire, &value) in z_words.iter().zip(z_word_values.iter()) {
		witness[wire] = Word(value);
	}
	if circuit.populate_wire_witness(&mut witness).is_err() {
		return false;
	}
	verify_constraints(cs, &witness.into_value_vec()).is_ok()
}

fn host_sample_in_ball_one_block(
	stream: &[u8],
) -> Option<([u64; mldsa44::N], [u64; mldsa44::TAU])> {
	assert_eq!(stream.len(), mldsa44::SAMPLE_IN_BALL_ONE_BLOCK_STREAM_BYTES);

	let mut signs = u64::from_le_bytes(stream[..8].try_into().unwrap());
	let mut draw_cursor = mldsa44::SAMPLE_IN_BALL_SIGN_BYTES;
	let mut coeffs = [0u64; mldsa44::N];
	let mut draw_counts = [0u64; mldsa44::TAU];

	for (round, i) in (mldsa44::N - mldsa44::TAU..mldsa44::N).enumerate() {
		let mut count = 0u64;
		let accepted_j = loop {
			if draw_cursor >= stream.len() {
				return None;
			}
			let draw = stream[draw_cursor] as usize;
			draw_cursor += 1;
			count += 1;
			if draw <= i {
				break draw;
			}
		};

		coeffs[i] = coeffs[accepted_j];
		coeffs[accepted_j] = if signs & 1 == 1 { u64::MAX } else { 1 };
		signs >>= 1;
		draw_counts[round] = count;
	}

	Some((coeffs, draw_counts))
}

fn deterministic_sample_in_ball_stream() -> [u8; mldsa44::SAMPLE_IN_BALL_ONE_BLOCK_STREAM_BYTES] {
	let mut stream = [0u8; mldsa44::SAMPLE_IN_BALL_ONE_BLOCK_STREAM_BYTES];
	stream[..8].copy_from_slice(&0x0000_0055_AA55_33CCu64.to_le_bytes());

	let mut pos = mldsa44::SAMPLE_IN_BALL_SIGN_BYTES;
	for (round, i) in (mldsa44::N - mldsa44::TAU..mldsa44::N).enumerate() {
		if round % 5 == 0 && i < 255 {
			stream[pos] = 255;
			pos += 1;
		}
		stream[pos] = ((i * 17 + round * 29) % (i + 1)) as u8;
		pos += 1;
	}

	stream
}

fn verify_mldsa44_sample_in_ball_stream_witness(
	stream: &[u8; mldsa44::SAMPLE_IN_BALL_ONE_BLOCK_STREAM_BYTES],
	draw_counts: &[u64; mldsa44::TAU],
	expected_coeffs: &[u64; mldsa44::N],
) -> bool {
	let stream_word_values = words_from_bytes(stream);
	let builder = CircuitBuilder::new();
	let stream_words: Vec<_> = (0..stream_word_values.len())
		.map(|_| builder.add_witness())
		.collect();
	let expected_coeff_wires: Vec<_> = (0..expected_coeffs.len())
		.map(|_| builder.add_witness())
		.collect();

	let sample = mldsa44_sample_in_ball_one_block_from_stream(&builder, &stream_words);
	for (i, (&computed, &expected)) in sample
		.coeffs
		.iter()
		.zip(expected_coeff_wires.iter())
		.enumerate()
	{
		builder.assert_eq(format!("sample_in_ball_coeff[{i}]"), computed, expected);
	}

	let circuit = builder.build();
	let cs = circuit.constraint_system();
	let mut witness = circuit.new_witness_filler();

	for (&wire, &value) in stream_words.iter().zip(stream_word_values.iter()) {
		witness[wire] = Word(value);
	}
	for (&wire, &value) in sample.draw_counts.iter().zip(draw_counts.iter()) {
		witness[wire] = Word(value);
	}
	for (&wire, &value) in expected_coeff_wires.iter().zip(expected_coeffs.iter()) {
		witness[wire] = Word(value);
	}

	if circuit.populate_wire_witness(&mut witness).is_err() {
		return false;
	}
	verify_constraints(cs, &witness.into_value_vec()).is_ok()
}

fn verify_mldsa44_one_block_hidden_hash_relation_witness(
	mu_and_w1_bytes: &[u8],
	c_tilde: &[u8; mldsa44::C_TILDE_BYTES],
	expected_coeffs: &[u64; mldsa44::N],
	draw_counts: &[u64; mldsa44::TAU],
) -> bool {
	assert_eq!(mu_and_w1_bytes.len(), mldsa44::FINAL_CHALLENGE_INPUT_BYTES);
	let mu_and_w1_word_values = words_from_bytes(mu_and_w1_bytes);
	let c_tilde_word_values = words_from_bytes(c_tilde);

	let builder = CircuitBuilder::new();
	let mu_and_w1_words: Vec<_> = (0..mu_and_w1_word_values.len())
		.map(|_| builder.add_witness())
		.collect();
	let c_tilde_words: Vec<_> = (0..c_tilde_word_values.len())
		.map(|_| builder.add_witness())
		.collect();
	let expected_coeff_wires: Vec<_> = (0..expected_coeffs.len())
		.map(|_| builder.add_witness())
		.collect();

	let sample = mldsa44_one_block_hidden_hash_relation(&builder, &c_tilde_words, &mu_and_w1_words);
	for (i, (&computed, &expected)) in sample
		.coeffs
		.iter()
		.zip(expected_coeff_wires.iter())
		.enumerate()
	{
		builder.assert_eq(format!("hidden_hash_sample_coeff[{i}]"), computed, expected);
	}

	let circuit = builder.build();
	let cs = circuit.constraint_system();
	let mut witness = circuit.new_witness_filler();

	for (&wire, &value) in mu_and_w1_words.iter().zip(mu_and_w1_word_values.iter()) {
		witness[wire] = Word(value);
	}
	for (&wire, &value) in c_tilde_words.iter().zip(c_tilde_word_values.iter()) {
		witness[wire] = Word(value);
	}
	for (&wire, &value) in sample.draw_counts.iter().zip(draw_counts.iter()) {
		witness[wire] = Word(value);
	}
	for (&wire, &value) in expected_coeff_wires.iter().zip(expected_coeffs.iter()) {
		witness[wire] = Word(value);
	}

	if circuit.populate_wire_witness(&mut witness).is_err() {
		return false;
	}
	verify_constraints(cs, &witness.into_value_vec()).is_ok()
}

fn verify_mldsa44_w1_encode_witness(w1_coeff_values: &[u64]) -> bool {
	let expected_word_values = pack_mldsa44_w1_coeffs(w1_coeff_values);
	let builder = CircuitBuilder::new();
	let w1_coeffs: Vec<_> = (0..w1_coeff_values.len())
		.map(|_| builder.add_witness())
		.collect();
	let expected_words: Vec<_> = (0..expected_word_values.len())
		.map(|_| builder.add_witness())
		.collect();

	let encoded = mldsa44_encode_w1(&builder, &w1_coeffs);
	for (i, (&computed, &expected)) in encoded.iter().zip(expected_words.iter()).enumerate() {
		builder.assert_eq(format!("mldsa44_w1_encode[{i}]"), computed, expected);
	}

	let circuit = builder.build();
	let cs = circuit.constraint_system();
	let mut witness = circuit.new_witness_filler();

	for (&wire, &value) in w1_coeffs.iter().zip(w1_coeff_values.iter()) {
		witness[wire] = Word(value);
	}
	for (&wire, &value) in expected_words.iter().zip(expected_word_values.iter()) {
		witness[wire] = Word(value);
	}

	if circuit.populate_wire_witness(&mut witness).is_err() {
		return false;
	}
	verify_constraints(cs, &witness.into_value_vec()).is_ok()
}

fn verify_mldsa44_w1encode_hash_relation_witness(
	mu: &[u8; mldsa44::MU_BYTES],
	w1_coeff_values: &[u64],
	c_tilde: &[u8; mldsa44::C_TILDE_BYTES],
	expected_coeffs: &[u64; mldsa44::N],
	draw_counts: &[u64; mldsa44::TAU],
) -> bool {
	let mu_word_values = words_from_bytes(mu);
	let c_tilde_word_values = words_from_bytes(c_tilde);

	let builder = CircuitBuilder::new();
	let mu_words: Vec<_> = (0..mu_word_values.len())
		.map(|_| builder.add_witness())
		.collect();
	let c_tilde_words: Vec<_> = (0..c_tilde_word_values.len())
		.map(|_| builder.add_witness())
		.collect();
	let w1_coeffs: Vec<_> = (0..w1_coeff_values.len())
		.map(|_| builder.add_witness())
		.collect();
	let expected_coeff_wires: Vec<_> = (0..expected_coeffs.len())
		.map(|_| builder.add_witness())
		.collect();

	let sample =
		mldsa44_one_block_w1encode_hash_relation(&builder, &c_tilde_words, &mu_words, &w1_coeffs);
	for (i, (&computed, &expected)) in sample
		.coeffs
		.iter()
		.zip(expected_coeff_wires.iter())
		.enumerate()
	{
		builder.assert_eq(format!("w1encode_hash_sample_coeff[{i}]"), computed, expected);
	}

	let circuit = builder.build();
	let cs = circuit.constraint_system();
	let mut witness = circuit.new_witness_filler();

	for (&wire, &value) in mu_words.iter().zip(mu_word_values.iter()) {
		witness[wire] = Word(value);
	}
	for (&wire, &value) in c_tilde_words.iter().zip(c_tilde_word_values.iter()) {
		witness[wire] = Word(value);
	}
	for (&wire, &value) in w1_coeffs.iter().zip(w1_coeff_values.iter()) {
		witness[wire] = Word(value);
	}
	for (&wire, &value) in sample.draw_counts.iter().zip(draw_counts.iter()) {
		witness[wire] = Word(value);
	}
	for (&wire, &value) in expected_coeff_wires.iter().zip(expected_coeffs.iter()) {
		witness[wire] = Word(value);
	}

	if circuit.populate_wire_witness(&mut witness).is_err() {
		return false;
	}
	verify_constraints(cs, &witness.into_value_vec()).is_ok()
}

fn verify_mldsa44_use_hint_witness(
	h_values: &[u64],
	r_values: &[u64],
	expected_w1: &[u64],
) -> bool {
	let builder = CircuitBuilder::new();
	let h_wires: Vec<_> = (0..h_values.len()).map(|_| builder.add_witness()).collect();
	let r_wires: Vec<_> = (0..r_values.len()).map(|_| builder.add_witness()).collect();
	let expected_wires: Vec<_> = (0..expected_w1.len())
		.map(|_| builder.add_witness())
		.collect();

	let w1 = mldsa44_use_hint(&builder, &h_wires, &r_wires);
	for (i, (&computed, &expected)) in w1.iter().zip(expected_wires.iter()).enumerate() {
		builder.assert_eq(format!("mldsa44_use_hint[{i}]"), computed, expected);
	}

	let circuit = builder.build();
	let cs = circuit.constraint_system();
	let mut witness = circuit.new_witness_filler();

	for (&wire, &value) in h_wires.iter().zip(h_values.iter()) {
		witness[wire] = Word(value);
	}
	for (&wire, &value) in r_wires.iter().zip(r_values.iter()) {
		witness[wire] = Word(value);
	}
	for (&wire, &value) in expected_wires.iter().zip(expected_w1.iter()) {
		witness[wire] = Word(value);
	}

	if circuit.populate_wire_witness(&mut witness).is_err() {
		return false;
	}
	verify_constraints(cs, &witness.into_value_vec()).is_ok()
}

fn verify_mldsa44_use_hint_hash_relation_witness(
	mu: &[u8; mldsa44::MU_BYTES],
	h_values: &[u64],
	r_values: &[u64],
	c_tilde: &[u8; mldsa44::C_TILDE_BYTES],
	expected_coeffs: &[u64; mldsa44::N],
	draw_counts: &[u64; mldsa44::TAU],
) -> bool {
	let mu_word_values = words_from_bytes(mu);
	let c_tilde_word_values = words_from_bytes(c_tilde);

	let builder = CircuitBuilder::new();
	let mu_words: Vec<_> = (0..mu_word_values.len())
		.map(|_| builder.add_witness())
		.collect();
	let c_tilde_words: Vec<_> = (0..c_tilde_word_values.len())
		.map(|_| builder.add_witness())
		.collect();
	let h_wires: Vec<_> = (0..h_values.len()).map(|_| builder.add_witness()).collect();
	let r_wires: Vec<_> = (0..r_values.len()).map(|_| builder.add_witness()).collect();
	let expected_coeff_wires: Vec<_> = (0..expected_coeffs.len())
		.map(|_| builder.add_witness())
		.collect();

	let sample = mldsa44_one_block_use_hint_hash_relation(
		&builder,
		&c_tilde_words,
		&mu_words,
		&h_wires,
		&r_wires,
	);
	for (i, (&computed, &expected)) in sample
		.coeffs
		.iter()
		.zip(expected_coeff_wires.iter())
		.enumerate()
	{
		builder.assert_eq(format!("use_hint_hash_sample_coeff[{i}]"), computed, expected);
	}

	let circuit = builder.build();
	let cs = circuit.constraint_system();
	let mut witness = circuit.new_witness_filler();

	for (&wire, &value) in mu_words.iter().zip(mu_word_values.iter()) {
		witness[wire] = Word(value);
	}
	for (&wire, &value) in c_tilde_words.iter().zip(c_tilde_word_values.iter()) {
		witness[wire] = Word(value);
	}
	for (&wire, &value) in h_wires.iter().zip(h_values.iter()) {
		witness[wire] = Word(value);
	}
	for (&wire, &value) in r_wires.iter().zip(r_values.iter()) {
		witness[wire] = Word(value);
	}
	for (&wire, &value) in sample.draw_counts.iter().zip(draw_counts.iter()) {
		witness[wire] = Word(value);
	}
	for (&wire, &value) in expected_coeff_wires.iter().zip(expected_coeffs.iter()) {
		witness[wire] = Word(value);
	}

	if circuit.populate_wire_witness(&mut witness).is_err() {
		return false;
	}
	verify_constraints(cs, &witness.into_value_vec()).is_ok()
}

fn verify_mldsa44_hint_decode_witness(
	h_bytes: &[u8; mldsa44::HINT_BYTES],
	expected_h: &[u64],
) -> bool {
	let h_word_values = words_from_bytes(h_bytes);
	let builder = CircuitBuilder::new();
	let h_words: Vec<_> = (0..h_word_values.len())
		.map(|_| builder.add_witness())
		.collect();
	let expected_wires: Vec<_> = (0..expected_h.len())
		.map(|_| builder.add_witness())
		.collect();

	let decoded = mldsa44_decode_hint_canonical(&builder, &h_words);
	for (i, (&computed, &expected)) in decoded.iter().zip(expected_wires.iter()).enumerate() {
		builder.assert_eq(format!("mldsa44_hint_decode[{i}]"), computed, expected);
	}

	let circuit = builder.build();
	let cs = circuit.constraint_system();
	let mut witness = circuit.new_witness_filler();

	for (&wire, &value) in h_words.iter().zip(h_word_values.iter()) {
		witness[wire] = Word(value);
	}
	for (&wire, &value) in expected_wires.iter().zip(expected_h.iter()) {
		witness[wire] = Word(value);
	}

	if circuit.populate_wire_witness(&mut witness).is_err() {
		return false;
	}
	verify_constraints(cs, &witness.into_value_vec()).is_ok()
}

fn verify_mldsa44_hint_matches_expanded_witness(
	h_bytes: &[u8; mldsa44::HINT_BYTES],
	h_values: &[u64],
) -> bool {
	let h_word_values = words_from_bytes(h_bytes);
	let builder = CircuitBuilder::new();
	let h_words: Vec<_> = (0..h_word_values.len())
		.map(|_| builder.add_witness())
		.collect();
	let h_coeffs: Vec<_> = (0..h_values.len()).map(|_| builder.add_witness()).collect();

	assert_mldsa44_hint_canonical_matches_expanded(&builder, &h_words, &h_coeffs);

	let circuit = builder.build();
	let cs = circuit.constraint_system();
	let mut witness = circuit.new_witness_filler();

	for (&wire, &value) in h_words.iter().zip(h_word_values.iter()) {
		witness[wire] = Word(value);
	}
	for (&wire, &value) in h_coeffs.iter().zip(h_values.iter()) {
		witness[wire] = Word(value);
	}

	if circuit.populate_wire_witness(&mut witness).is_err() {
		return false;
	}
	verify_constraints(cs, &witness.into_value_vec()).is_ok()
}

fn build_mldsa44_full_bit_heavy_one_block_circuit() -> binius_frontend::Circuit {
	let builder = CircuitBuilder::new();
	let c_tilde: Vec<_> = (0..mldsa44::C_TILDE_BYTES / 8)
		.map(|_| builder.add_witness())
		.collect();
	let z_words: Vec<_> = (0..mldsa44::Z_PACKED_WORDS)
		.map(|_| builder.add_witness())
		.collect();
	let h_coeffs: Vec<_> = (0..mldsa44::W1_COEFFICIENTS)
		.map(|_| builder.add_witness())
		.collect();
	let mu_words: Vec<_> = (0..mldsa44::MU_WORDS)
		.map(|_| builder.add_inout())
		.collect();
	let w_approx_coeffs: Vec<_> = (0..mldsa44::W1_COEFFICIENTS)
		.map(|_| builder.add_witness())
		.collect();

	mldsa44_full_bit_heavy_one_block_relation(
		&builder,
		&c_tilde,
		&z_words,
		&h_coeffs,
		&mu_words,
		&w_approx_coeffs,
	);

	builder.build()
}

fn build_mldsa44_full_bit_heavy_one_block_canonical_hint_circuit() -> binius_frontend::Circuit {
	let builder = CircuitBuilder::new();
	let c_tilde: Vec<_> = (0..mldsa44::C_TILDE_BYTES / 8)
		.map(|_| builder.add_witness())
		.collect();
	let z_words: Vec<_> = (0..mldsa44::Z_PACKED_WORDS)
		.map(|_| builder.add_witness())
		.collect();
	let h_words: Vec<_> = (0..mldsa44::HINT_WORDS)
		.map(|_| builder.add_witness())
		.collect();
	let mu_words: Vec<_> = (0..mldsa44::MU_WORDS)
		.map(|_| builder.add_inout())
		.collect();
	let w_approx_coeffs: Vec<_> = (0..mldsa44::W1_COEFFICIENTS)
		.map(|_| builder.add_witness())
		.collect();

	mldsa44_full_bit_heavy_one_block_canonical_hint_relation(
		&builder,
		&c_tilde,
		&z_words,
		&h_words,
		&mu_words,
		&w_approx_coeffs,
	);

	builder.build()
}

fn build_mldsa44_full_bit_heavy_one_block_canonical_hint_matched_circuit()
-> binius_frontend::Circuit {
	let builder = CircuitBuilder::new();
	let c_tilde: Vec<_> = (0..mldsa44::C_TILDE_BYTES / 8)
		.map(|_| builder.add_witness())
		.collect();
	let z_words: Vec<_> = (0..mldsa44::Z_PACKED_WORDS)
		.map(|_| builder.add_witness())
		.collect();
	let h_words: Vec<_> = (0..mldsa44::HINT_WORDS)
		.map(|_| builder.add_witness())
		.collect();
	let h_coeffs: Vec<_> = (0..mldsa44::W1_COEFFICIENTS)
		.map(|_| builder.add_witness())
		.collect();
	let mu_words: Vec<_> = (0..mldsa44::MU_WORDS)
		.map(|_| builder.add_inout())
		.collect();
	let w_approx_coeffs: Vec<_> = (0..mldsa44::W1_COEFFICIENTS)
		.map(|_| builder.add_witness())
		.collect();

	mldsa44_full_bit_heavy_one_block_canonical_hint_matched_relation(
		&builder,
		&c_tilde,
		&z_words,
		&h_words,
		&h_coeffs,
		&mu_words,
		&w_approx_coeffs,
	);

	builder.build()
}

fn verify_mldsa44_z_decode_witness(packed_y_coeff_values: &[u64]) {
	let z_word_values = pack_mldsa44_z_y_coeffs(packed_y_coeff_values);
	let builder = CircuitBuilder::new();
	let z_words: Vec<_> = (0..z_word_values.len())
		.map(|_| builder.add_witness())
		.collect();
	let expected_coeffs: Vec<_> = (0..packed_y_coeff_values.len())
		.map(|_| builder.add_witness())
		.collect();

	let decoded = mldsa44_decode_z_packed_y(&builder, &z_words);
	for (i, (&decoded, &expected)) in decoded.iter().zip(expected_coeffs.iter()).enumerate() {
		builder.assert_eq(format!("mldsa44_z_decode[{i}]"), decoded, expected);
	}

	let circuit = builder.build();
	let cs = circuit.constraint_system();
	let mut witness = circuit.new_witness_filler();
	for (&wire, &value) in z_words.iter().zip(z_word_values.iter()) {
		witness[wire] = Word(value);
	}
	for (&wire, &value) in expected_coeffs.iter().zip(packed_y_coeff_values.iter()) {
		witness[wire] = Word(value);
	}

	circuit.populate_wire_witness(&mut witness).unwrap();
	verify_constraints(cs, &witness.into_value_vec())
		.expect("Circuit constraints should be satisfied");
}

fn verify_mldsa44_z_norm_witness(packed_y_coeff_values: &[u64]) -> bool {
	let builder = CircuitBuilder::new();
	let packed_y_coeffs: Vec<_> = (0..packed_y_coeff_values.len())
		.map(|_| builder.add_witness())
		.collect();
	assert_mldsa44_z_norm_from_packed_y(&builder, &packed_y_coeffs);

	let circuit = builder.build();
	let cs = circuit.constraint_system();
	let mut witness = circuit.new_witness_filler();
	for (&wire, &value) in packed_y_coeffs.iter().zip(packed_y_coeff_values.iter()) {
		witness[wire] = Word(value);
	}
	if circuit.populate_wire_witness(&mut witness).is_err() {
		return false;
	}
	verify_constraints(cs, &witness.into_value_vec()).is_ok()
}

#[test]
fn mldsa44_final_challenge_hash_matches_shake256() {
	let mut rng = StdRng::seed_from_u64(0x4D4C4453413434);
	let mut input = vec![0u8; mldsa44::FINAL_CHALLENGE_INPUT_BYTES];
	rng.fill_bytes(&mut input);

	let mut hasher = Shake256::default();
	hasher.update(&input);
	let mut reader = hasher.finalize_xof();
	let mut expected = [0u8; mldsa44::C_TILDE_BYTES];
	reader.read(&mut expected);

	let builder = CircuitBuilder::new();
	let input_wires: Vec<_> = (0..input.len().div_ceil(8))
		.map(|_| builder.add_witness())
		.collect();
	let expected_wires: [Wire; mldsa44::C_TILDE_BYTES / 8] =
		std::array::from_fn(|_| builder.add_witness());

	let computed = mldsa44_final_challenge_hash(&builder, &input_wires);
	for i in 0..computed.len() {
		builder.assert_eq(format!("c_tilde_prime[{i}]"), computed[i], expected_wires[i]);
	}

	let circuit = builder.build();
	let cs = circuit.constraint_system();
	let mut witness = circuit.new_witness_filler();

	for (i, chunk) in input.chunks(8).enumerate() {
		let word = u64::from_le_bytes(chunk.try_into().unwrap());
		witness[input_wires[i]] = Word(word);
	}
	for (i, chunk) in expected.chunks(8).enumerate() {
		let word = u64::from_le_bytes(chunk.try_into().unwrap());
		witness[expected_wires[i]] = Word(word);
	}

	circuit.populate_wire_witness(&mut witness).unwrap();
	verify_constraints(cs, &witness.into_value_vec())
		.expect("Circuit constraints should be satisfied");
}

#[test]
fn mldsa44_sample_in_ball_one_block_stream_matches_shake256() {
	let mut rng = StdRng::seed_from_u64(0x53414D504C453434);
	let mut c_tilde = [0u8; mldsa44::C_TILDE_BYTES];
	rng.fill_bytes(&mut c_tilde);

	let mut hasher = Shake256::default();
	hasher.update(&c_tilde);
	let mut reader = hasher.finalize_xof();
	let mut expected = [0u8; mldsa44::SAMPLE_IN_BALL_ONE_BLOCK_STREAM_BYTES];
	reader.read(&mut expected);

	let builder = CircuitBuilder::new();
	let c_tilde_wires: [Wire; mldsa44::C_TILDE_BYTES / 8] =
		std::array::from_fn(|_| builder.add_witness());
	let expected_wires: [Wire; mldsa44::SAMPLE_IN_BALL_ONE_BLOCK_STREAM_BYTES / 8] =
		std::array::from_fn(|_| builder.add_witness());

	let computed = mldsa44_sample_in_ball_one_block_stream(&builder, &c_tilde_wires);
	for i in 0..computed.len() {
		builder.assert_eq(format!("sample_in_ball_stream[{i}]"), computed[i], expected_wires[i]);
	}

	let circuit = builder.build();
	let cs = circuit.constraint_system();
	let mut witness = circuit.new_witness_filler();

	for (i, chunk) in c_tilde.chunks(8).enumerate() {
		let word = u64::from_le_bytes(chunk.try_into().unwrap());
		witness[c_tilde_wires[i]] = Word(word);
	}
	for (i, chunk) in expected.chunks(8).enumerate() {
		let word = u64::from_le_bytes(chunk.try_into().unwrap());
		witness[expected_wires[i]] = Word(word);
	}

	circuit.populate_wire_witness(&mut witness).unwrap();
	verify_constraints(cs, &witness.into_value_vec())
		.expect("Circuit constraints should be satisfied");
}

#[test]
fn mldsa44_sample_in_ball_one_block_accepts_valid_trace() {
	let stream = deterministic_sample_in_ball_stream();
	let (expected_coeffs, draw_counts) =
		host_sample_in_ball_one_block(&stream).expect("stream should satisfy one-block cap");

	assert!(verify_mldsa44_sample_in_ball_stream_witness(&stream, &draw_counts, &expected_coeffs,));
}

#[test]
fn mldsa44_sample_in_ball_one_block_rejects_wrong_draw_count() {
	let stream = deterministic_sample_in_ball_stream();
	let (expected_coeffs, mut draw_counts) =
		host_sample_in_ball_one_block(&stream).expect("stream should satisfy one-block cap");
	let bad_round = draw_counts.iter().position(|&count| count > 1).unwrap();
	draw_counts[bad_round] -= 1;

	assert!(
		!verify_mldsa44_sample_in_ball_stream_witness(&stream, &draw_counts, &expected_coeffs,)
	);
}

#[test]
fn mldsa44_one_block_hidden_hash_relation_accepts_valid_binding() {
	let mut rng = StdRng::seed_from_u64(0x4841534852454C);
	let mut mu_and_w1_bytes = vec![0u8; mldsa44::FINAL_CHALLENGE_INPUT_BYTES];
	rng.fill_bytes(&mut mu_and_w1_bytes);

	let mut final_hasher = Shake256::default();
	final_hasher.update(&mu_and_w1_bytes);
	let mut final_reader = final_hasher.finalize_xof();
	let mut c_tilde = [0u8; mldsa44::C_TILDE_BYTES];
	final_reader.read(&mut c_tilde);

	let mut sample_hasher = Shake256::default();
	sample_hasher.update(&c_tilde);
	let mut sample_reader = sample_hasher.finalize_xof();
	let mut sample_stream = [0u8; mldsa44::SAMPLE_IN_BALL_ONE_BLOCK_STREAM_BYTES];
	sample_reader.read(&mut sample_stream);
	let (expected_coeffs, draw_counts) =
		host_sample_in_ball_one_block(&sample_stream).expect("stream should satisfy cap");

	assert!(verify_mldsa44_one_block_hidden_hash_relation_witness(
		&mu_and_w1_bytes,
		&c_tilde,
		&expected_coeffs,
		&draw_counts,
	));
}

#[test]
fn mldsa44_one_block_hidden_hash_relation_rejects_c_tilde_mutation() {
	let mut rng = StdRng::seed_from_u64(0x4841534852454D);
	let mut mu_and_w1_bytes = vec![0u8; mldsa44::FINAL_CHALLENGE_INPUT_BYTES];
	rng.fill_bytes(&mut mu_and_w1_bytes);

	let mut final_hasher = Shake256::default();
	final_hasher.update(&mu_and_w1_bytes);
	let mut final_reader = final_hasher.finalize_xof();
	let mut c_tilde = [0u8; mldsa44::C_TILDE_BYTES];
	final_reader.read(&mut c_tilde);

	let mut sample_hasher = Shake256::default();
	sample_hasher.update(&c_tilde);
	let mut sample_reader = sample_hasher.finalize_xof();
	let mut sample_stream = [0u8; mldsa44::SAMPLE_IN_BALL_ONE_BLOCK_STREAM_BYTES];
	sample_reader.read(&mut sample_stream);
	let (expected_coeffs, draw_counts) =
		host_sample_in_ball_one_block(&sample_stream).expect("stream should satisfy cap");

	c_tilde[0] ^= 1;

	assert!(!verify_mldsa44_one_block_hidden_hash_relation_witness(
		&mu_and_w1_bytes,
		&c_tilde,
		&expected_coeffs,
		&draw_counts,
	));
}

#[test]
fn mldsa44_w1_encode_matches_simple_bitpack_layout() {
	let mut w1 = vec![0u64; mldsa44::W1_COEFFICIENTS];
	for (i, coeff) in w1.iter_mut().enumerate() {
		*coeff = ((i as u64 * 17) + 5) % (mldsa44::W1_COEFF_MAX + 1);
	}
	w1[0] = 0;
	w1[1] = mldsa44::W1_COEFF_MAX;
	w1[10] = 1;
	w1[mldsa44::W1_COEFFICIENTS - 1] = 42;

	assert!(verify_mldsa44_w1_encode_witness(&w1));
}

#[test]
fn mldsa44_w1_encode_rejects_out_of_range_coeff() {
	let mut w1 = vec![0u64; mldsa44::W1_COEFFICIENTS];
	w1[123] = mldsa44::W1_COEFF_MAX + 1;

	let builder = CircuitBuilder::new();
	let w1_coeffs: Vec<_> = (0..w1.len()).map(|_| builder.add_witness()).collect();
	mldsa44_encode_w1(&builder, &w1_coeffs);

	let circuit = builder.build();
	let cs = circuit.constraint_system();
	let mut witness = circuit.new_witness_filler();
	for (&wire, &value) in w1_coeffs.iter().zip(w1.iter()) {
		witness[wire] = Word(value);
	}

	assert!(
		circuit.populate_wire_witness(&mut witness).is_err()
			|| verify_constraints(cs, &witness.into_value_vec()).is_err()
	);
}

#[test]
fn mldsa44_w1encode_hash_relation_accepts_valid_binding() {
	let mut rng = StdRng::seed_from_u64(0x5731454E434F4445);
	let mut mu = [0u8; mldsa44::MU_BYTES];
	rng.fill_bytes(&mut mu);
	let mut w1 = vec![0u64; mldsa44::W1_COEFFICIENTS];
	for coeff in &mut w1 {
		*coeff = (rng.next_u64() % (mldsa44::W1_COEFF_MAX + 1)) as u64;
	}

	let w1_words = pack_mldsa44_w1_coeffs(&w1);
	let mut hash_input = Vec::with_capacity(mldsa44::FINAL_CHALLENGE_INPUT_BYTES);
	hash_input.extend_from_slice(&mu);
	for word in w1_words {
		hash_input.extend_from_slice(&word.to_le_bytes());
	}

	let mut final_hasher = Shake256::default();
	final_hasher.update(&hash_input);
	let mut final_reader = final_hasher.finalize_xof();
	let mut c_tilde = [0u8; mldsa44::C_TILDE_BYTES];
	final_reader.read(&mut c_tilde);

	let mut sample_hasher = Shake256::default();
	sample_hasher.update(&c_tilde);
	let mut sample_reader = sample_hasher.finalize_xof();
	let mut sample_stream = [0u8; mldsa44::SAMPLE_IN_BALL_ONE_BLOCK_STREAM_BYTES];
	sample_reader.read(&mut sample_stream);
	let (expected_coeffs, draw_counts) =
		host_sample_in_ball_one_block(&sample_stream).expect("stream should satisfy cap");

	assert!(verify_mldsa44_w1encode_hash_relation_witness(
		&mu,
		&w1,
		&c_tilde,
		&expected_coeffs,
		&draw_counts,
	));
}

#[test]
fn mldsa44_use_hint_matches_reference_boundaries() {
	let mut h = vec![0u64; mldsa44::W1_COEFFICIENTS];
	let mut r = vec![0u64; mldsa44::W1_COEFFICIENTS];

	let cases = [
		(0, 0),
		(1, 0),
		(1, 1),
		(1, mldsa44::GAMMA2),
		(1, mldsa44::TWO_GAMMA2),
		(1, mldsa44::TWO_GAMMA2 - 1),
		(1, mldsa44::TWO_GAMMA2 + 1),
		(1, 43 * mldsa44::TWO_GAMMA2),
		(1, 43 * mldsa44::TWO_GAMMA2 + 1),
		(1, 43 * mldsa44::TWO_GAMMA2 + mldsa44::GAMMA2),
		(1, mldsa44::Q - 1),
		(0, mldsa44::Q - 1),
	];
	let initial_weight: u64 = cases.iter().map(|&(hint, _)| hint).sum();
	for (i, &(hint, value)) in cases.iter().enumerate() {
		h[i] = hint;
		r[i] = value;
	}
	// Sprinkle additional ones across the rest of the polynomial without exceeding `omega`.
	let mut remaining_ones = mldsa44::OMEGA - initial_weight;
	for i in cases.len()..mldsa44::W1_COEFFICIENTS {
		if remaining_ones > 0 && i % 13 == 0 {
			h[i] = 1;
			remaining_ones -= 1;
		}
		r[i] = ((i as u64 * 65_537) + 12_345) % mldsa44::Q;
	}

	let expected: Vec<_> = h
		.iter()
		.zip(r.iter())
		.map(|(&hint, &value)| host_mldsa44_use_hint(hint, value))
		.collect();

	assert!(verify_mldsa44_use_hint_witness(&h, &r, &expected));
}

#[test]
fn mldsa44_use_hint_rejects_out_of_range_r() {
	let h = vec![0u64; mldsa44::W1_COEFFICIENTS];
	let mut r = vec![0u64; mldsa44::W1_COEFFICIENTS];
	r[9] = mldsa44::Q;
	let expected = vec![0u64; mldsa44::W1_COEFFICIENTS];

	assert!(!verify_mldsa44_use_hint_witness(&h, &r, &expected));
}

#[test]
fn mldsa44_use_hint_hash_relation_accepts_valid_binding() {
	let mut rng = StdRng::seed_from_u64(0x55534548494E54);
	let mut mu = [0u8; mldsa44::MU_BYTES];
	rng.fill_bytes(&mut mu);
	let mut h = vec![0u64; mldsa44::W1_COEFFICIENTS];
	let mut r = vec![0u64; mldsa44::W1_COEFFICIENTS];
	// Pick `omega` distinct random positions to set h = 1; everything else stays 0. The
	// canonical relation rejects any witness whose Hamming weight exceeds omega.
	let mut hint_positions: Vec<usize> = (0..mldsa44::W1_COEFFICIENTS).collect();
	for swap_idx in 0..(mldsa44::OMEGA_USIZE) {
		let pick = swap_idx + (rng.next_u64() as usize) % (hint_positions.len() - swap_idx);
		hint_positions.swap(swap_idx, pick);
	}
	for &pos in &hint_positions[..mldsa44::OMEGA_USIZE] {
		h[pos] = 1;
	}
	for ri in r.iter_mut() {
		*ri = rng.next_u64() % mldsa44::Q;
	}
	let w1: Vec<_> = h
		.iter()
		.zip(r.iter())
		.map(|(&hint, &value)| host_mldsa44_use_hint(hint, value))
		.collect();

	let w1_words = pack_mldsa44_w1_coeffs(&w1);
	let mut hash_input = Vec::with_capacity(mldsa44::FINAL_CHALLENGE_INPUT_BYTES);
	hash_input.extend_from_slice(&mu);
	for word in w1_words {
		hash_input.extend_from_slice(&word.to_le_bytes());
	}

	let mut final_hasher = Shake256::default();
	final_hasher.update(&hash_input);
	let mut final_reader = final_hasher.finalize_xof();
	let mut c_tilde = [0u8; mldsa44::C_TILDE_BYTES];
	final_reader.read(&mut c_tilde);

	let mut sample_hasher = Shake256::default();
	sample_hasher.update(&c_tilde);
	let mut sample_reader = sample_hasher.finalize_xof();
	let mut sample_stream = [0u8; mldsa44::SAMPLE_IN_BALL_ONE_BLOCK_STREAM_BYTES];
	sample_reader.read(&mut sample_stream);
	let (expected_coeffs, draw_counts) =
		host_sample_in_ball_one_block(&sample_stream).expect("stream should satisfy cap");

	assert!(verify_mldsa44_use_hint_hash_relation_witness(
		&mu,
		&h,
		&r,
		&c_tilde,
		&expected_coeffs,
		&draw_counts,
	));
}

#[test]
fn mldsa44_hint_decode_accepts_canonical_encoding() {
	let h_bytes = pack_mldsa44_hint_bytes(&[vec![0, 9, 71], vec![3, 12], vec![], vec![4, 255]]);
	let expected = host_mldsa44_decode_hint(&h_bytes).unwrap();

	assert!(verify_mldsa44_hint_decode_witness(&h_bytes, &expected));
}

#[test]
fn mldsa44_hint_decode_rejects_non_increasing_segment() {
	let mut h_bytes = pack_mldsa44_hint_bytes(&[vec![0, 9, 71], vec![3, 12], vec![], vec![4, 255]]);
	h_bytes[1] = 0;
	let expected = vec![0u64; mldsa44::W1_COEFFICIENTS];

	assert!(!verify_mldsa44_hint_decode_witness(&h_bytes, &expected));
}

#[test]
fn mldsa44_hint_decode_rejects_nonzero_unused_byte() {
	let mut h_bytes = pack_mldsa44_hint_bytes(&[vec![0, 9], vec![], vec![], vec![]]);
	h_bytes[20] = 7;
	let expected = vec![0u64; mldsa44::W1_COEFFICIENTS];

	assert!(!verify_mldsa44_hint_decode_witness(&h_bytes, &expected));
}

#[test]
fn mldsa44_hint_match_accepts_canonical_expanded_pair() {
	let h_bytes = pack_mldsa44_hint_bytes(&[vec![0, 9, 71], vec![3, 12], vec![], vec![4, 255]]);
	let h = host_mldsa44_decode_hint(&h_bytes).unwrap();

	assert!(verify_mldsa44_hint_matches_expanded_witness(&h_bytes, &h));
}

#[test]
fn mldsa44_hint_match_rejects_missing_expanded_bit() {
	let h_bytes = pack_mldsa44_hint_bytes(&[vec![0, 9, 71], vec![3, 12], vec![], vec![4, 255]]);
	let mut h = host_mldsa44_decode_hint(&h_bytes).unwrap();
	h[9] = 0;

	assert!(!verify_mldsa44_hint_matches_expanded_witness(&h_bytes, &h));
}

#[test]
fn mldsa44_hint_match_rejects_extra_expanded_bit() {
	let h_bytes = pack_mldsa44_hint_bytes(&[vec![0, 9, 71], vec![3, 12], vec![], vec![4, 255]]);
	let mut h = host_mldsa44_decode_hint(&h_bytes).unwrap();
	h[1] = 1;

	assert!(!verify_mldsa44_hint_matches_expanded_witness(&h_bytes, &h));
}

#[test]
fn mldsa44_full_bit_heavy_one_block_circuit_builds_and_prints_stats() {
	let circuit = build_mldsa44_full_bit_heavy_one_block_circuit();
	let stat = CircuitStat::collect(&circuit);

	println!("{stat}");
	println!("Gate composition JSON:");
	println!("{}", circuit.simple_json_dump());

	assert_eq!(stat.n_inout, mldsa44::MU_WORDS);
	assert!(stat.n_gates > 0);
	assert!(stat.n_and_constraints > 0);
	assert_eq!(stat.n_mul_constraints, 0);
}

#[test]
fn mldsa44_full_bit_heavy_one_block_canonical_hint_circuit_builds_and_prints_stats() {
	let circuit = build_mldsa44_full_bit_heavy_one_block_canonical_hint_circuit();
	let stat = CircuitStat::collect(&circuit);

	println!("{stat}");
	println!("Canonical hint gate composition JSON:");
	println!("{}", circuit.simple_json_dump());

	assert_eq!(stat.n_inout, mldsa44::MU_WORDS);
	assert!(stat.n_gates > 0);
	assert!(stat.n_and_constraints > 0);
	assert_eq!(stat.n_mul_constraints, 0);
}

#[test]
fn mldsa44_full_bit_heavy_one_block_canonical_hint_matched_circuit_builds_and_prints_stats() {
	let circuit = build_mldsa44_full_bit_heavy_one_block_canonical_hint_matched_circuit();
	let stat = CircuitStat::collect(&circuit);

	println!("{stat}");
	println!("Canonical matched hint gate composition JSON:");
	println!("{}", circuit.simple_json_dump());

	assert_eq!(stat.n_inout, mldsa44::MU_WORDS);
	assert!(stat.n_gates > 0);
	assert!(stat.n_and_constraints > 0);
	assert_eq!(stat.n_mul_constraints, 0);
}

#[test]
fn mldsa44_z_norm_accepts_boundary_values() {
	let mut packed_y = vec![mldsa44::Z_NORM_PACKED_Y_MIN; mldsa44::Z_COEFFICIENTS];
	packed_y[17] = mldsa44::Z_NORM_PACKED_Y_MAX;
	packed_y[999] = (mldsa44::Z_NORM_PACKED_Y_MIN + mldsa44::Z_NORM_PACKED_Y_MAX) / 2;

	assert!(verify_mldsa44_z_norm_witness(&packed_y));
}

#[test]
fn mldsa44_decode_z_packed_y_matches_bitpack_layout() {
	let mut packed_y = vec![0u64; mldsa44::Z_COEFFICIENTS];
	for (i, coeff) in packed_y.iter_mut().enumerate() {
		*coeff = ((i as u64 * 65_537) + 0x12345) & ((1 << mldsa44::Z_BITS_PER_COEFF) - 1);
	}
	packed_y[0] = 0;
	packed_y[3] = (1 << mldsa44::Z_BITS_PER_COEFF) - 1;
	packed_y[7] = 0b10_1010_1111_0000_0101;
	packed_y[mldsa44::Z_COEFFICIENTS - 1] = 0x2_0001;

	verify_mldsa44_z_decode_witness(&packed_y);
}

#[test]
fn mldsa44_z_packed_bytes_norm_accepts_boundary_values() {
	let mut packed_y = vec![mldsa44::Z_NORM_PACKED_Y_MIN; mldsa44::Z_COEFFICIENTS];
	packed_y[17] = mldsa44::Z_NORM_PACKED_Y_MAX;
	packed_y[999] = (mldsa44::Z_NORM_PACKED_Y_MIN + mldsa44::Z_NORM_PACKED_Y_MAX) / 2;

	assert!(verify_mldsa44_z_packed_bytes_norm_witness(&packed_y));
}

#[test]
fn mldsa44_z_packed_bytes_norm_rejects_below_min() {
	let mut packed_y = vec![mldsa44::Z_NORM_PACKED_Y_MIN; mldsa44::Z_COEFFICIENTS];
	packed_y[41] = mldsa44::Z_NORM_PACKED_Y_MIN - 1;

	assert!(!verify_mldsa44_z_packed_bytes_norm_witness(&packed_y));
}

#[test]
fn mldsa44_z_packed_bytes_norm_rejects_above_max() {
	let mut packed_y = vec![mldsa44::Z_NORM_PACKED_Y_MIN; mldsa44::Z_COEFFICIENTS];
	packed_y[271] = mldsa44::Z_NORM_PACKED_Y_MAX + 1;

	assert!(!verify_mldsa44_z_packed_bytes_norm_witness(&packed_y));
}

#[test]
fn mldsa44_z_norm_rejects_below_min() {
	let mut packed_y = vec![mldsa44::Z_NORM_PACKED_Y_MIN; mldsa44::Z_COEFFICIENTS];
	packed_y[41] = mldsa44::Z_NORM_PACKED_Y_MIN - 1;

	assert!(!verify_mldsa44_z_norm_witness(&packed_y));
}

#[test]
fn mldsa44_z_norm_rejects_above_max() {
	let mut packed_y = vec![mldsa44::Z_NORM_PACKED_Y_MIN; mldsa44::Z_COEFFICIENTS];
	packed_y[271] = mldsa44::Z_NORM_PACKED_Y_MAX + 1;

	assert!(!verify_mldsa44_z_norm_witness(&packed_y));
}

// =====================================================================
// Hardening tests
//
// The block below intentionally targets:
//
// 1. The Barrett-style `mldsa44_high_bits` formula (`((r + 127) >> 7) * 11275 + (1 << 23)) >>
//    24`) plus the `t > 43 ? 0 : t` clamp that emulates FIPS Decompose's special case for
//    `r' - r0' = q - 1`. We sweep all interesting boundary `r` and a wide pseudo-random sample
//    against the host reference.
// 2. The implicit "SampleInBall trail can only end at the first valid byte" property by
//    mutating the witnessed `draw_counts` away from the FIPS-conformant trail.
// 3. The `||z||_infty < gamma1 - beta` range across the canonical 18-bit encoding extremes (0,
//    2^18 - 1) and across coefficient slots that exercise both the single-word and two-word
//    decode branches.
// 4. The FIPS `omega` weight bound on `h`, now enforced inside `mldsa44_use_hint`.
// 5. Rejection of structurally invalid hint encodings (endpoint > omega, non-monotone
//    endpoints, position >= N).
//
// These cases are not redundant with the existing tests; they each exercise a constraint that
// the prior tests do not isolate.
// =====================================================================

fn build_mldsa44_high_bits_circuit(
	r_count: usize,
) -> (binius_frontend::Circuit, Vec<Wire>, Vec<Wire>) {
	let builder = CircuitBuilder::new();
	let r_wires: Vec<_> = (0..r_count).map(|_| builder.add_witness()).collect();
	let expected_wires: Vec<_> = (0..r_count).map(|_| builder.add_witness()).collect();
	for (i, (&r, &expected)) in r_wires.iter().zip(expected_wires.iter()).enumerate() {
		let computed = mldsa44_high_bits(&builder, r);
		builder.assert_eq(format!("hb_eq[{i}]"), computed, expected);
	}
	(builder.build(), r_wires, expected_wires)
}

fn run_mldsa44_high_bits_check(r_values: &[u64]) -> bool {
	let expected: Vec<u64> = r_values
		.iter()
		.map(|&r| host_mldsa44_high_bits(r))
		.collect();
	let (circuit, r_wires, expected_wires) = build_mldsa44_high_bits_circuit(r_values.len());
	let cs = circuit.constraint_system();
	let mut witness = circuit.new_witness_filler();
	for (&w, &v) in r_wires.iter().zip(r_values.iter()) {
		witness[w] = Word(v);
	}
	for (&w, &v) in expected_wires.iter().zip(expected.iter()) {
		witness[w] = Word(v);
	}
	if circuit.populate_wire_witness(&mut witness).is_err() {
		return false;
	}
	verify_constraints(cs, &witness.into_value_vec()).is_ok()
}

#[test]
fn mldsa44_high_bits_matches_reference_at_boundaries() {
	// Boundaries that flex the Barrett formula and the special-case clamp.
	let alpha = mldsa44::TWO_GAMMA2;
	let half_alpha = mldsa44::GAMMA2;
	let q_minus_one = mldsa44::Q - 1;
	let mut r = vec![
		0,
		1,
		half_alpha - 1,
		half_alpha,
		half_alpha + 1,
		alpha - 1,
		alpha,
		alpha + 1,
		alpha + half_alpha,
		alpha + half_alpha + 1,
		(mldsa44::Q - 1) / 2,
		(mldsa44::Q - 1) / 2 + 1,
		// Boundary between normal r1=43 and the special case (r1 -> 0).
		q_minus_one - half_alpha - 1,
		q_minus_one - half_alpha,
		q_minus_one - half_alpha + 1,
		q_minus_one - 1,
		q_minus_one,
	];
	// All multiples of alpha across [0, q-1].
	for k in 0..=(q_minus_one / alpha) {
		r.push(k * alpha);
	}
	// All k*alpha + half_alpha (round-half-up boundary points).
	for k in 0..((q_minus_one - half_alpha) / alpha) {
		r.push(k * alpha + half_alpha);
		r.push(k * alpha + half_alpha + 1);
	}
	assert!(run_mldsa44_high_bits_check(&r));
}

#[test]
fn mldsa44_high_bits_matches_reference_random_sweep() {
	// 4096 pseudo-random values uniformly across [0, q-1]. Plenty of coverage given the
	// formula is monotone-ish and the clamp is the only nonlinear part.
	let mut rng = StdRng::seed_from_u64(0x4842_4954_5343_4845);
	let mut r = Vec::with_capacity(4096);
	for _ in 0..4096 {
		r.push(rng.next_u64() % mldsa44::Q);
	}
	assert!(run_mldsa44_high_bits_check(&r));
}

#[test]
fn mldsa44_high_bits_rejects_r_eq_q() {
	let mut r = vec![0u64; 16];
	r[7] = mldsa44::Q;
	// `host_mldsa44_high_bits(q)` would assertion-fail; so we don't compute it.
	// Just check that any consistent expected value fails because the in-circuit r-range
	// assertion fires.
	let expected = vec![0u64; 16];
	let (circuit, r_wires, expected_wires) = build_mldsa44_high_bits_circuit(r.len());
	let cs = circuit.constraint_system();
	let mut witness = circuit.new_witness_filler();
	for (&w, &v) in r_wires.iter().zip(r.iter()) {
		witness[w] = Word(v);
	}
	for (&w, &v) in expected_wires.iter().zip(expected.iter()) {
		witness[w] = Word(v);
	}
	assert!(
		circuit.populate_wire_witness(&mut witness).is_err()
			|| verify_constraints(cs, &witness.into_value_vec()).is_err()
	);
}

#[test]
fn mldsa44_use_hint_rejects_h_ge_two() {
	let mut h = vec![0u64; mldsa44::W1_COEFFICIENTS];
	let r = vec![0u64; mldsa44::W1_COEFFICIENTS];
	// h = 2 with everything else 0 — weight is 2 (tolerated), but per-coefficient bit check
	// must fire.
	h[100] = 2;
	let expected = vec![0u64; mldsa44::W1_COEFFICIENTS];
	assert!(!verify_mldsa44_use_hint_witness(&h, &r, &expected));
}

#[test]
fn mldsa44_use_hint_rejects_weight_above_omega() {
	let mut h = vec![0u64; mldsa44::W1_COEFFICIENTS];
	// Set omega + 1 ones, all valid 0/1 bits.
	for idx in 0..(mldsa44::OMEGA_USIZE + 1) {
		h[idx] = 1;
	}
	let r = vec![0u64; mldsa44::W1_COEFFICIENTS];
	// host_mldsa44_use_hint(1, 0) returns W1_COEFF_MAX = 43 (since r1=0, r0=0 -> h=1 dec to 43).
	let expected: Vec<_> = h
		.iter()
		.zip(r.iter())
		.map(|(&hh, &rr)| host_mldsa44_use_hint(hh, rr))
		.collect();
	assert!(!verify_mldsa44_use_hint_witness(&h, &r, &expected));
}

#[test]
fn mldsa44_use_hint_accepts_weight_eq_omega() {
	let mut h = vec![0u64; mldsa44::W1_COEFFICIENTS];
	for idx in 0..mldsa44::OMEGA_USIZE {
		h[idx * 11 % mldsa44::W1_COEFFICIENTS] = 1;
	}
	// Saturation may have collisions; recount and ensure exactly omega ones.
	let mut count: u64 = h.iter().sum();
	while count < mldsa44::OMEGA {
		let pos = (count as usize * 17 + 13) % mldsa44::W1_COEFFICIENTS;
		if h[pos] == 0 {
			h[pos] = 1;
			count += 1;
		} else {
			let mut p = (pos + 1) % mldsa44::W1_COEFFICIENTS;
			while h[p] == 1 {
				p = (p + 1) % mldsa44::W1_COEFFICIENTS;
			}
			h[p] = 1;
			count += 1;
		}
	}
	assert_eq!(h.iter().sum::<u64>(), mldsa44::OMEGA);

	let r: Vec<u64> = (0..mldsa44::W1_COEFFICIENTS)
		.map(|i| (i as u64 * 65_537) % mldsa44::Q)
		.collect();
	let expected: Vec<_> = h
		.iter()
		.zip(r.iter())
		.map(|(&hh, &rr)| host_mldsa44_use_hint(hh, rr))
		.collect();
	assert!(verify_mldsa44_use_hint_witness(&h, &r, &expected));
}

#[test]
fn mldsa44_sample_in_ball_rejects_zero_draw_count() {
	let stream = deterministic_sample_in_ball_stream();
	let (expected_coeffs, mut draw_counts) = host_sample_in_ball_one_block(&stream).unwrap();
	// Zero out the first round's draw count (must be >= 1).
	draw_counts[0] = 0;
	assert!(
		!verify_mldsa44_sample_in_ball_stream_witness(&stream, &draw_counts, &expected_coeffs,)
	);
}

#[test]
fn mldsa44_sample_in_ball_rejects_overcap_draw_counts() {
	let stream = deterministic_sample_in_ball_stream();
	let (expected_coeffs, mut draw_counts) = host_sample_in_ball_one_block(&stream).unwrap();
	// Inflate the last round's draw_count to push the cumulative cursor past 128.
	let total: u64 = draw_counts.iter().sum();
	let last = draw_counts.last_mut().unwrap();
	*last += mldsa44::SAMPLE_IN_BALL_ONE_BLOCK_DRAW_CAP_BYTES as u64 + 1 - total;
	assert!(
		!verify_mldsa44_sample_in_ball_stream_witness(&stream, &draw_counts, &expected_coeffs,)
	);
}

#[test]
fn mldsa44_sample_in_ball_rejects_wrong_sign_bit() {
	let mut stream = deterministic_sample_in_ball_stream();
	let (expected_coeffs, draw_counts) = host_sample_in_ball_one_block(&stream).unwrap();
	// Flip a sign bit so the host-derived `expected_coeffs` no longer match what the circuit
	// computes from the (mutated) stream.
	stream[0] ^= 0b1;
	assert!(
		!verify_mldsa44_sample_in_ball_stream_witness(&stream, &draw_counts, &expected_coeffs,)
	);
}

#[test]
fn mldsa44_sample_in_ball_rejects_swapped_skip_and_accept() {
	// Construct a stream where round 0 (i = 217) needs at least one rejection. The host then
	// accepts the second byte; we mutate the witness to swap the bytes so the trail looks
	// like "accept byte 0, skip byte 1". The skipped-draw constraint requires skipped > i,
	// so a draw <= i in a skipped slot fails.
	let mut stream = [0u8; mldsa44::SAMPLE_IN_BALL_ONE_BLOCK_STREAM_BYTES];
	stream[..8].copy_from_slice(&0u64.to_le_bytes());
	// Round 0: i = 217. Put 220 (>217), then 100 (<=217), then easy bytes for later rounds.
	stream[mldsa44::SAMPLE_IN_BALL_SIGN_BYTES] = 220;
	stream[mldsa44::SAMPLE_IN_BALL_SIGN_BYTES + 1] = 100;
	let mut byte_pos = mldsa44::SAMPLE_IN_BALL_SIGN_BYTES + 2;
	for i in (mldsa44::N - mldsa44::TAU + 1)..mldsa44::N {
		stream[byte_pos] = (i % (i + 1)) as u8;
		byte_pos += 1;
	}
	let (expected_coeffs, _draw_counts) =
		host_sample_in_ball_one_block(&stream).expect("constructed stream should satisfy cap");

	// Build a malicious draw_counts that pretends the first round consumed only 1 byte, so
	// `accepted_pos = 0` and `accepted_draw = 220 > 217`.
	let mut bad_draw_counts = [1u64; mldsa44::TAU];
	// The remaining cumulative draws need to match the host trail to keep [round 1..] sound,
	// otherwise we would also fail there. So just assert the swap fails.
	bad_draw_counts[0] = 1;
	// The rest of `bad_draw_counts` is whatever the host says.
	let host = host_sample_in_ball_one_block(&stream).unwrap().1;
	bad_draw_counts[1..].copy_from_slice(&host[1..]);
	assert!(!verify_mldsa44_sample_in_ball_stream_witness(
		&stream,
		&bad_draw_counts,
		&expected_coeffs,
	));
}

#[rstest]
#[case(0)]
#[case(3)] // straddles bit-offsets [54, 71], two-word decode branch
#[case(4)] // bit-offsets [72, 89], single-word in word_idx=1
#[case(255)] // last coefficient of the first polynomial
#[case(256)] // first coefficient of the second polynomial
#[case(341)] // bit-offset 6138 = word 95 + shift 58 (single-word, last allowed shift)
#[case(1023)] // last overall coefficient (bit-offset 18414, ends at bit 18431)
fn mldsa44_z_packed_bytes_norm_rejects_below_min_at_slot(#[case] slot: usize) {
	let mut packed_y = vec![mldsa44::Z_NORM_PACKED_Y_MIN; mldsa44::Z_COEFFICIENTS];
	packed_y[slot] = mldsa44::Z_NORM_PACKED_Y_MIN - 1;
	assert!(!verify_mldsa44_z_packed_bytes_norm_witness(&packed_y));
}

#[rstest]
#[case(0)]
#[case(3)]
#[case(4)]
#[case(255)]
#[case(256)]
#[case(341)]
#[case(1023)]
fn mldsa44_z_packed_bytes_norm_rejects_above_max_at_slot(#[case] slot: usize) {
	let mut packed_y = vec![mldsa44::Z_NORM_PACKED_Y_MIN; mldsa44::Z_COEFFICIENTS];
	packed_y[slot] = mldsa44::Z_NORM_PACKED_Y_MAX + 1;
	assert!(!verify_mldsa44_z_packed_bytes_norm_witness(&packed_y));
}

#[rstest]
#[case(0)] // canonical extreme y = 0 (z = gamma1)
#[case(78)] // y = beta, just below MIN
#[case(262_066)] // y = MAX + 1 = 2*gamma1 - beta
#[case(262_143)] // canonical extreme y = 2^18 - 1 (z = gamma1 + 1 - 2^18)
fn mldsa44_z_packed_bytes_norm_rejects_canonical_extreme(#[case] y: u64) {
	// y must fit in 18 bits, but can be outside the norm-allowed window.
	assert!(y < (1u64 << mldsa44::Z_BITS_PER_COEFF));
	let mut packed_y = vec![mldsa44::Z_NORM_PACKED_Y_MIN; mldsa44::Z_COEFFICIENTS];
	packed_y[42] = y;
	assert!(!verify_mldsa44_z_packed_bytes_norm_witness(&packed_y));
}

#[rstest]
#[case(0)] // shift = 0
#[case(10)] // shift = 60, two-word straddle branch
#[case(11)] // shift = 2, second word
#[case(255)] // end of poly 0
#[case(1023)] // last overall slot
fn mldsa44_w1_encode_rejects_out_of_range_at_slot(#[case] slot: usize) {
	let mut w1 = vec![0u64; mldsa44::W1_COEFFICIENTS];
	w1[slot] = mldsa44::W1_COEFF_MAX + 1;

	let builder = CircuitBuilder::new();
	let w1_coeffs: Vec<_> = (0..w1.len()).map(|_| builder.add_witness()).collect();
	mldsa44_encode_w1(&builder, &w1_coeffs);

	let circuit = builder.build();
	let cs = circuit.constraint_system();
	let mut witness = circuit.new_witness_filler();
	for (&wire, &value) in w1_coeffs.iter().zip(w1.iter()) {
		witness[wire] = Word(value);
	}
	assert!(
		circuit.populate_wire_witness(&mut witness).is_err()
			|| verify_constraints(cs, &witness.into_value_vec()).is_err()
	);
}

#[test]
fn mldsa44_w1_encode_round_trip_via_host_pack() {
	// Random valid `w1`. Host pack matches circuit pack iff both implement FIPS SimpleBitPack.
	let mut rng = StdRng::seed_from_u64(0x57_3145_4E43_4F44);
	let mut w1 = vec![0u64; mldsa44::W1_COEFFICIENTS];
	for coeff in w1.iter_mut() {
		*coeff = rng.next_u64() % (mldsa44::W1_COEFF_MAX + 1);
	}
	assert!(verify_mldsa44_w1_encode_witness(&w1));
}

fn build_canonical_hint_bytes_with(positions_per_poly: &[Vec<u8>]) -> [u8; mldsa44::HINT_BYTES] {
	assert_eq!(positions_per_poly.len(), mldsa44::K);
	let mut h = [0u8; mldsa44::HINT_BYTES];
	let mut idx = 0usize;
	for (poly_idx, positions) in positions_per_poly.iter().enumerate() {
		for &pos in positions {
			assert!(idx < mldsa44::OMEGA_USIZE);
			h[idx] = pos;
			idx += 1;
		}
		h[mldsa44::OMEGA_USIZE + poly_idx] = idx as u8;
	}
	h
}

fn expand_hint_bytes_to_coeffs(h_bytes: &[u8; mldsa44::HINT_BYTES]) -> Vec<u64> {
	// Permissive expander: do not fail on out-of-range positions; they just get ignored. Used
	// for negative test fixtures.
	let mut h = vec![0u64; mldsa44::W1_COEFFICIENTS];
	let mut idx = 0usize;
	for poly_idx in 0..mldsa44::K {
		let endpoint = h_bytes[mldsa44::OMEGA_USIZE + poly_idx] as usize;
		let end = endpoint.min(mldsa44::OMEGA_USIZE);
		if idx <= end {
			for &pos in &h_bytes[idx..end] {
				let pos = pos as usize;
				if pos < mldsa44::N {
					h[poly_idx * mldsa44::N + pos] = 1;
				}
			}
		}
		idx = end;
	}
	h
}

#[test]
fn mldsa44_hint_decode_rejects_endpoint_above_omega() {
	// Endpoint > omega should be rejected by the canonical-decoder structural check.
	let mut h = build_canonical_hint_bytes_with(&[
		(0..3).collect::<Vec<u8>>(),
		(10..15).collect::<Vec<u8>>(),
		(50..55).collect::<Vec<u8>>(),
		(100..105).collect::<Vec<u8>>(),
	]);
	// Forge the last endpoint to omega + 1.
	h[mldsa44::OMEGA_USIZE + mldsa44::K - 1] = (mldsa44::OMEGA + 1) as u8;
	let expected = expand_hint_bytes_to_coeffs(&h);
	assert!(!verify_mldsa44_hint_decode_witness(&h, &expected));
}

#[test]
fn mldsa44_hint_decode_rejects_endpoint_not_monotone() {
	let mut h = build_canonical_hint_bytes_with(&[
		(0..3).collect::<Vec<u8>>(),
		(10..15).collect::<Vec<u8>>(),
		(50..55).collect::<Vec<u8>>(),
		(100..105).collect::<Vec<u8>>(),
	]);
	// Swap two endpoints to break monotonicity.
	h.swap(mldsa44::OMEGA_USIZE + 1, mldsa44::OMEGA_USIZE + 2);
	let expected = expand_hint_bytes_to_coeffs(&h);
	assert!(!verify_mldsa44_hint_decode_witness(&h, &expected));
}

#[test]
fn mldsa44_hint_decode_rejects_position_at_n() {
	// position 255 is valid (< 256), 0 is valid; but a position that equals N must fail. Since
	// our hint bytes are u8, we use a value reachable from the byte witness.
	let mut h = build_canonical_hint_bytes_with(&[
		vec![0, 1, 2],
		vec![10, 11, 12],
		vec![100],
		vec![200, 201],
	]);
	// Replace one position with a synthetic byte; out-of-range bytes are always >= N for
	// ML-DSA-44 since N=256, but `u8::MAX = 255 < N`, so we instead flip the sign of that
	// invariant by overriding N to 256: the canonical decoder asserts `position < N = 256`.
	// Every u8 position trivially satisfies this for ML-DSA-44, so we can't construct an
	// out-of-range u8 position. We instead test the canonical-decoder's "unused position must
	// be zero" invariant by injecting a non-zero byte after the last endpoint.
	let last_endpoint = h[mldsa44::OMEGA_USIZE + mldsa44::K - 1] as usize;
	assert!(last_endpoint < mldsa44::OMEGA_USIZE);
	h[last_endpoint] = 7; // first unused position byte must be zero
	let expected = expand_hint_bytes_to_coeffs(&h);
	assert!(!verify_mldsa44_hint_decode_witness(&h, &expected));
}

/// Canonical decode + matched-expansion regression: the first poly being empty puts later
/// segments at start = 0, exercising the `start == 0, poly_idx > 0` corner of the
/// strictly-increasing check (where `prev_pos = h_bytes[pos_idx - 1]` is read for
/// `pos_idx > 0` while `start == 0` makes `has_previous` track every in-segment slot).
#[test]
fn mldsa44_hint_decode_accepts_empty_leading_segment() {
	let h_bytes = pack_mldsa44_hint_bytes(&[vec![], vec![1, 4, 17], vec![], vec![5, 200]]);
	let expected = host_mldsa44_decode_hint(&h_bytes).expect("canonical encoding");
	assert!(verify_mldsa44_hint_decode_witness(&h_bytes, &expected));
	assert!(verify_mldsa44_hint_matches_expanded_witness(&h_bytes, &expected));
}

/// Canonical decode + matched-expansion regression: only the trailing poly carries any
/// positions, exercising the path where every non-final segment is empty and the strictly-
/// increasing check must remain inert across multiple zero-length segments.
#[test]
fn mldsa44_hint_decode_accepts_only_trailing_segment() {
	let h_bytes = pack_mldsa44_hint_bytes(&[vec![], vec![], vec![], vec![0, 9, 71, 200]]);
	let expected = host_mldsa44_decode_hint(&h_bytes).expect("canonical encoding");
	assert!(verify_mldsa44_hint_decode_witness(&h_bytes, &expected));
	assert!(verify_mldsa44_hint_matches_expanded_witness(&h_bytes, &expected));
}

/// Differential coverage between the two `z`-norm gadgets:
/// `assert_mldsa44_z_norm_from_packed_y` consumes already-decoded packed-`y` coefficients,
/// while `assert_mldsa44_z_packed_bytes_norm` decodes from the BitPack byte stream first.
/// Their accept sets must agree for every coefficient slot we exercise here.
///
/// Each case picks one slot, pins every other slot to a value safely inside the norm window
/// and varies the chosen slot across the canonical 18-bit extremes plus the in-window
/// midpoint. Both gadgets must independently agree with the host-side decision
/// `Z_NORM_PACKED_Y_MIN <= y <= Z_NORM_PACKED_Y_MAX`.
#[rstest]
#[case(0, 0)] // canonical extreme y = 0 (z = gamma1)
#[case(0, mldsa44::Z_NORM_PACKED_Y_MIN - 1)] // y = 78
#[case(0, mldsa44::Z_NORM_PACKED_Y_MIN)] // boundary low
#[case(0, mldsa44::Z_NORM_PACKED_Y_MIN + 1)]
#[case(0, (mldsa44::Z_NORM_PACKED_Y_MIN + mldsa44::Z_NORM_PACKED_Y_MAX) / 2)]
#[case(0, mldsa44::Z_NORM_PACKED_Y_MAX - 1)]
#[case(0, mldsa44::Z_NORM_PACKED_Y_MAX)] // boundary high
#[case(0, mldsa44::Z_NORM_PACKED_Y_MAX + 1)] // y = 262066
#[case(0, (1u64 << mldsa44::Z_BITS_PER_COEFF) - 1)] // y = 262143
#[case(3, mldsa44::Z_NORM_PACKED_Y_MIN - 1)] // two-word straddle slot
#[case(3, mldsa44::Z_NORM_PACKED_Y_MIN)]
#[case(3, mldsa44::Z_NORM_PACKED_Y_MAX)]
#[case(3, mldsa44::Z_NORM_PACKED_Y_MAX + 1)]
#[case(255, mldsa44::Z_NORM_PACKED_Y_MIN)] // last coeff of poly 0
#[case(256, mldsa44::Z_NORM_PACKED_Y_MAX)] // first coeff of poly 1
#[case(1023, mldsa44::Z_NORM_PACKED_Y_MIN - 1)] // last overall coeff
#[case(1023, mldsa44::Z_NORM_PACKED_Y_MAX + 1)]
fn mldsa44_z_norm_packed_y_and_packed_bytes_gadgets_agree(#[case] slot: usize, #[case] y: u64) {
	assert!(y < (1u64 << mldsa44::Z_BITS_PER_COEFF), "y must fit in 18 bits");
	let mut packed_y = vec![mldsa44::Z_NORM_PACKED_Y_MIN; mldsa44::Z_COEFFICIENTS];
	packed_y[slot] = y;
	let host_accepts = y >= mldsa44::Z_NORM_PACKED_Y_MIN && y <= mldsa44::Z_NORM_PACKED_Y_MAX;
	let from_packed_y = verify_mldsa44_z_norm_witness(&packed_y);
	let from_packed_bytes = verify_mldsa44_z_packed_bytes_norm_witness(&packed_y);
	assert_eq!(from_packed_y, host_accepts, "_packed_y disagrees with host at slot {slot} y={y}",);
	assert_eq!(
		from_packed_bytes, host_accepts,
		"_packed_bytes_norm disagrees with host at slot {slot} y={y}",
	);
	assert_eq!(
		from_packed_y, from_packed_bytes,
		"the two z-norm gadgets disagree at slot {slot} y={y}",
	);
}
