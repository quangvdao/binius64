// Copyright 2025 Irreducible Inc.

use std::time::{Duration, Instant};

use binius_circuits::mldsa::{
	assert_mldsa44_hint_canonical_matches_expanded, assert_mldsa44_z_packed_bytes_norm, mldsa44,
	mldsa44_encode_w1, mldsa44_final_challenge_hash,
	mldsa44_full_bit_heavy_one_block_canonical_hint_matched_relation,
	mldsa44_full_bit_heavy_one_block_relation, mldsa44_sample_in_ball_one_block_from_stream,
	mldsa44_sample_in_ball_one_block_sparse_from_stream, mldsa44_sample_in_ball_one_block_stream,
	mldsa44_use_hint_coeff,
};
use binius_core::{constraint_system::ValueVec, word::Word};
use binius_field::arch::OptimalPackedB128;
use binius_frontend::{Circuit, CircuitBuilder, CircuitStat, Wire};
use binius_prover::{Prover, hash::parallel_compression::ParallelCompressionAdaptor};
use binius_transcript::{ProverTranscript, VerifierTranscript};
use binius_verifier::{
	Verifier,
	config::StdChallenger,
	hash::{StdCompression, StdDigest},
};
use sha3::{
	Shake256,
	digest::{ExtendableOutput, Update, XofReader},
};

const LOG_INV_RATE: usize = 1;

#[derive(Clone, Copy)]
enum MldsaCircuitKind {
	ExpandedHint,
	CanonicalMatchedHint,
}

impl MldsaCircuitKind {
	fn name(self) -> &'static str {
		match self {
			Self::ExpandedHint => "mldsa44_expanded_hint",
			Self::CanonicalMatchedHint => "mldsa44_canonical_matched_hint",
		}
	}
}

struct MldsaMeasurementCircuit {
	circuit: Circuit,
	c_tilde: Vec<Wire>,
	z_words: Vec<Wire>,
	h_words: Vec<Wire>,
	h_coeffs: Vec<Wire>,
	mu_words: Vec<Wire>,
	w_approx_coeffs: Vec<Wire>,
	draw_counts: [Wire; mldsa44::TAU],
}

struct MldsaMeasurementWitness {
	c_tilde_words: Vec<u64>,
	z_words: Vec<u64>,
	h_words: Vec<u64>,
	h_coeffs: Vec<u64>,
	mu_words: Vec<u64>,
	w_approx_coeffs: Vec<u64>,
	draw_counts: [u64; mldsa44::TAU],
}

struct Timed<T> {
	value: T,
	elapsed: Duration,
}

#[derive(Debug)]
struct MldsaProofMetrics {
	circuit_build: Duration,
	witness_generation: Duration,
	verifier_setup: Duration,
	prover_setup: Duration,
	prove: Duration,
	verify: Duration,
	proof_size_bytes: usize,
}

fn timed<T>(f: impl FnOnce() -> T) -> Timed<T> {
	let start = Instant::now();
	let value = f();
	Timed {
		value,
		elapsed: start.elapsed(),
	}
}

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

fn append_words_as_bytes(out: &mut Vec<u8>, words: &[u64]) {
	for word in words {
		out.extend_from_slice(&word.to_le_bytes());
	}
}

fn shake256(input: &[u8], out_bytes: usize) -> Vec<u8> {
	let mut hasher = Shake256::default();
	hasher.update(input);
	let mut reader = hasher.finalize_xof();
	let mut out = vec![0u8; out_bytes];
	reader.read(&mut out);
	out
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

fn host_mldsa44_decode_hint(h_bytes: &[u8; mldsa44::HINT_BYTES]) -> Vec<u64> {
	let mut h = vec![0u64; mldsa44::W1_COEFFICIENTS];
	let mut index = 0usize;
	for poly_idx in 0..mldsa44::K {
		let endpoint = h_bytes[mldsa44::OMEGA_USIZE + poly_idx] as usize;
		assert!(endpoint >= index);
		assert!(endpoint <= mldsa44::OMEGA_USIZE);
		for &pos in &h_bytes[index..endpoint] {
			h[poly_idx * mldsa44::N + pos as usize] = 1;
		}
		index = endpoint;
	}
	assert!(h_bytes[index..mldsa44::OMEGA_USIZE].iter().all(|&x| x == 0));
	h
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

fn build_mldsa44_circuit(kind: MldsaCircuitKind) -> MldsaMeasurementCircuit {
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

	let relation = match kind {
		MldsaCircuitKind::ExpandedHint => mldsa44_full_bit_heavy_one_block_relation(
			&builder,
			&c_tilde,
			&z_words,
			&h_coeffs,
			&mu_words,
			&w_approx_coeffs,
		),
		MldsaCircuitKind::CanonicalMatchedHint => {
			mldsa44_full_bit_heavy_one_block_canonical_hint_matched_relation(
				&builder,
				&c_tilde,
				&z_words,
				&h_words,
				&h_coeffs,
				&mu_words,
				&w_approx_coeffs,
			)
		}
	};
	let draw_counts = relation.sample_in_ball.draw_counts;
	let circuit = builder.build();

	MldsaMeasurementCircuit {
		circuit,
		c_tilde,
		z_words,
		h_words,
		h_coeffs,
		mu_words,
		w_approx_coeffs,
		draw_counts,
	}
}

fn make_valid_witness() -> MldsaMeasurementWitness {
	let h_bytes = pack_mldsa44_hint_bytes(&[
		vec![0, 9, 71],
		vec![3, 12, 90],
		vec![8, 44],
		vec![4, 128, 255],
	]);
	let h_coeffs = host_mldsa44_decode_hint(&h_bytes);
	let h_words = words_from_bytes(&h_bytes);

	let mut mu = [0u8; mldsa44::MU_BYTES];
	for (i, byte) in mu.iter_mut().enumerate() {
		*byte = (0xA5u8).wrapping_add((i as u8).wrapping_mul(17));
	}
	let mu_words = words_from_bytes(&mu);

	let z_packed_y = vec![mldsa44::Z_NORM_PACKED_Y_MIN + 123; mldsa44::Z_COEFFICIENTS];
	let z_words = pack_mldsa44_z_y_coeffs(&z_packed_y);

	for seed in 0..10_000u64 {
		let w_approx_coeffs: Vec<_> = (0..mldsa44::W1_COEFFICIENTS)
			.map(|i| ((i as u64 * 65_537) + seed * 1_299_721 + 12_345) % mldsa44::Q)
			.collect();
		let w1: Vec<_> = h_coeffs
			.iter()
			.zip(w_approx_coeffs.iter())
			.map(|(&h, &r)| host_mldsa44_use_hint(h, r))
			.collect();
		let w1_words = pack_mldsa44_w1_coeffs(&w1);

		let mut final_hash_input = Vec::with_capacity(mldsa44::FINAL_CHALLENGE_INPUT_BYTES);
		final_hash_input.extend_from_slice(&mu);
		append_words_as_bytes(&mut final_hash_input, &w1_words);
		assert_eq!(final_hash_input.len(), mldsa44::FINAL_CHALLENGE_INPUT_BYTES);

		let c_tilde = shake256(&final_hash_input, mldsa44::C_TILDE_BYTES);
		let sample_stream = shake256(&c_tilde, mldsa44::SAMPLE_IN_BALL_ONE_BLOCK_STREAM_BYTES);
		if let Some((_coeffs, draw_counts)) = host_sample_in_ball_one_block(&sample_stream) {
			return MldsaMeasurementWitness {
				c_tilde_words: words_from_bytes(&c_tilde),
				z_words,
				h_words,
				h_coeffs,
				mu_words,
				w_approx_coeffs,
				draw_counts,
			};
		}
	}

	panic!("failed to find a one-block SampleInBall witness");
}

fn populate_mldsa44_witness(
	built: &MldsaMeasurementCircuit,
	witness_values: &MldsaMeasurementWitness,
) -> ValueVec {
	let mut filler = built.circuit.new_witness_filler();

	for (&wire, &value) in built
		.c_tilde
		.iter()
		.zip(witness_values.c_tilde_words.iter())
	{
		filler[wire] = Word(value);
	}
	for (&wire, &value) in built.z_words.iter().zip(witness_values.z_words.iter()) {
		filler[wire] = Word(value);
	}
	for (&wire, &value) in built.h_words.iter().zip(witness_values.h_words.iter()) {
		filler[wire] = Word(value);
	}
	for (&wire, &value) in built.h_coeffs.iter().zip(witness_values.h_coeffs.iter()) {
		filler[wire] = Word(value);
	}
	for (&wire, &value) in built.mu_words.iter().zip(witness_values.mu_words.iter()) {
		filler[wire] = Word(value);
	}
	for (&wire, &value) in built
		.w_approx_coeffs
		.iter()
		.zip(witness_values.w_approx_coeffs.iter())
	{
		filler[wire] = Word(value);
	}
	for (&wire, &value) in built
		.draw_counts
		.iter()
		.zip(witness_values.draw_counts.iter())
	{
		filler[wire] = Word(value);
	}

	built.circuit.populate_wire_witness(&mut filler).unwrap();
	filler.into_value_vec()
}

fn prove_and_measure(kind: MldsaCircuitKind) -> MldsaProofMetrics {
	let built = timed(|| build_mldsa44_circuit(kind));
	let stat = CircuitStat::collect(&built.value.circuit);
	println!("{} circuit stats:\n{stat}", kind.name());

	let witness_values = make_valid_witness();
	let witness = timed(|| populate_mldsa44_witness(&built.value, &witness_values));
	let cs = built.value.circuit.constraint_system().clone();

	let verifier = timed(|| {
		Verifier::<StdDigest, _>::setup(cs, LOG_INV_RATE, StdCompression::default()).unwrap()
	});

	let prover = timed(|| {
		Prover::<OptimalPackedB128, _, StdDigest>::setup(
			verifier.value.clone(),
			ParallelCompressionAdaptor::new(StdCompression::default()),
		)
		.unwrap()
	});

	let proof = timed(|| {
		let mut prover_transcript = ProverTranscript::new(StdChallenger::default());
		prover
			.value
			.prove(witness.value.clone(), &mut prover_transcript)
			.unwrap();
		prover_transcript.finalize()
	});

	let proof_size_bytes = proof.value.len();
	let verify = timed(|| {
		let mut verifier_transcript =
			VerifierTranscript::new(StdChallenger::default(), proof.value.clone());
		verifier
			.value
			.verify(witness.value.public(), &mut verifier_transcript)
			.unwrap();
		verifier_transcript.finalize().unwrap();
	});

	MldsaProofMetrics {
		circuit_build: built.elapsed,
		witness_generation: witness.elapsed,
		verifier_setup: verifier.elapsed,
		prover_setup: prover.elapsed,
		prove: proof.elapsed,
		verify: verify.elapsed,
		proof_size_bytes,
	}
}

fn print_metrics(name: &str, metrics: &MldsaProofMetrics) {
	println!(
		"{name} proof metrics:\n\
		 circuit_build: {:?}\n\
		 witness_generation: {:?}\n\
		 verifier_setup: {:?}\n\
		 prover_setup: {:?}\n\
		 prove: {:?}\n\
		 verify: {:?}\n\
		 proof_size_bytes: {}",
		metrics.circuit_build,
		metrics.witness_generation,
		metrics.verifier_setup,
		metrics.prover_setup,
		metrics.prove,
		metrics.verify,
		metrics.proof_size_bytes
	);
}

fn print_circuit_stats(name: &str, build: impl FnOnce(&CircuitBuilder)) {
	let builder = CircuitBuilder::new();
	build(&builder);
	let circuit = builder.build();
	let stat = CircuitStat::collect(&circuit);
	println!("{name} stats:\n{stat}");
	println!("{name} gate composition JSON:");
	println!("{}", circuit.simple_json_dump());
}

#[test]
#[ignore = "stats-only circuit component breakdown"]
fn mldsa44_constraint_breakdown() {
	print_circuit_stats("z_packed_decode_and_norm", |builder| {
		let z_words: Vec<_> = (0..mldsa44::Z_PACKED_WORDS)
			.map(|_| builder.add_witness())
			.collect();
		assert_mldsa44_z_packed_bytes_norm(builder, &z_words);
	});

	print_circuit_stats("hint_canonical_matches_expanded", |builder| {
		let h_words: Vec<_> = (0..mldsa44::HINT_WORDS)
			.map(|_| builder.add_witness())
			.collect();
		let h_coeffs: Vec<_> = (0..mldsa44::W1_COEFFICIENTS)
			.map(|_| builder.add_witness())
			.collect();
		assert_mldsa44_hint_canonical_matches_expanded(builder, &h_words, &h_coeffs);
	});

	print_circuit_stats("use_hint_no_weight_check", |builder| {
		let h_coeffs: Vec<_> = (0..mldsa44::W1_COEFFICIENTS)
			.map(|_| builder.add_witness())
			.collect();
		let w_approx_coeffs: Vec<_> = (0..mldsa44::W1_COEFFICIENTS)
			.map(|_| builder.add_witness())
			.collect();
		for (&h, &r) in h_coeffs.iter().zip(w_approx_coeffs.iter()) {
			mldsa44_use_hint_coeff(builder, h, r);
		}
	});

	print_circuit_stats("w1_encode", |builder| {
		let w1_coeffs: Vec<_> = (0..mldsa44::W1_COEFFICIENTS)
			.map(|_| builder.add_witness())
			.collect();
		mldsa44_encode_w1(builder, &w1_coeffs);
	});

	print_circuit_stats("final_challenge_shake256_7f", |builder| {
		let input_words: Vec<_> = (0..mldsa44::FINAL_CHALLENGE_INPUT_WORDS)
			.map(|_| builder.add_witness())
			.collect();
		mldsa44_final_challenge_hash(builder, &input_words);
	});

	print_circuit_stats("sample_in_ball_shake256_1f", |builder| {
		let c_tilde: Vec<_> = (0..mldsa44::C_TILDE_BYTES / 8)
			.map(|_| builder.add_witness())
			.collect();
		mldsa44_sample_in_ball_one_block_stream(builder, &c_tilde);
	});

	print_circuit_stats("sample_in_ball_rejection_only", |builder| {
		let stream: Vec<_> = (0..mldsa44::SAMPLE_IN_BALL_ONE_BLOCK_STREAM_BYTES / 8)
			.map(|_| builder.add_witness())
			.collect();
		mldsa44_sample_in_ball_one_block_from_stream(builder, &stream);
	});

	print_circuit_stats("sample_in_ball_sparse_rejection_only", |builder| {
		let stream: Vec<_> = (0..mldsa44::SAMPLE_IN_BALL_ONE_BLOCK_STREAM_BYTES / 8)
			.map(|_| builder.add_witness())
			.collect();
		mldsa44_sample_in_ball_one_block_sparse_from_stream(builder, &stream);
	});
}

#[test]
#[ignore = "expensive end-to-end proof measurement"]
fn measure_mldsa44_canonical_matched_proof() {
	let metrics = prove_and_measure(MldsaCircuitKind::CanonicalMatchedHint);
	print_metrics(MldsaCircuitKind::CanonicalMatchedHint.name(), &metrics);
}

#[test]
#[ignore = "expensive end-to-end proof measurement"]
fn measure_mldsa44_expanded_hint_proof() {
	let metrics = prove_and_measure(MldsaCircuitKind::ExpandedHint);
	print_metrics(MldsaCircuitKind::ExpandedHint.name(), &metrics);
}
