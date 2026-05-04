// Copyright 2025 Irreducible Inc.

use std::time::{Duration, Instant};

use binius_circuits::mldsa::{
	Mldsa44, Mldsa65, Mldsa87, MldsaParams, assert_hint_canonical_matches_expanded_for,
	assert_z_packed_bytes_norm_for, encode_w1_for, final_challenge_hash_for,
	full_bit_heavy_fixed_cap_canonical_hint_matched_relation_for,
	full_bit_heavy_fixed_cap_relation_for, mldsa44, mldsa65, mldsa87,
	sample_in_ball_fixed_cap_from_stream_for, sample_in_ball_fixed_cap_sparse_from_stream_for,
	sample_in_ball_fixed_cap_stream_for, use_hint_coeff_for,
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
	fn suffix(self) -> &'static str {
		match self {
			Self::ExpandedHint => "expanded_hint",
			Self::CanonicalMatchedHint => "canonical_matched_hint",
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
	draw_counts: Vec<Wire>,
}

struct MldsaMeasurementWitness {
	c_tilde_words: Vec<u64>,
	z_words: Vec<u64>,
	h_words: Vec<u64>,
	h_coeffs: Vec<u64>,
	mu_words: Vec<u64>,
	w_approx_coeffs: Vec<u64>,
	draw_counts: Vec<u64>,
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

#[derive(Clone, Copy)]
struct ComponentStat {
	name: &'static str,
	n_gates: usize,
	n_evaluations: usize,
	n_and_constraints: usize,
	n_mul_constraints: usize,
	n_inout: usize,
	n_private: usize,
	total_committed: usize,
	scratch: usize,
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

fn pack_z_y_coeffs_for<P: MldsaParams>(packed_y_coeff_values: &[u64]) -> Vec<u64> {
	assert_eq!(packed_y_coeff_values.len(), P::Z_COEFFICIENTS);

	let mut words = vec![0u64; P::Z_PACKED_WORDS];
	for (coeff_idx, &coeff) in packed_y_coeff_values.iter().enumerate() {
		assert!(coeff < (1 << P::Z_BITS_PER_COEFF));
		for bit in 0..P::Z_BITS_PER_COEFF {
			if (coeff >> bit) & 1 == 1 {
				let bit_idx = coeff_idx * P::Z_BITS_PER_COEFF + bit;
				words[bit_idx / 64] |= 1 << (bit_idx % 64);
			}
		}
	}

	words
}

fn pack_w1_coeffs_for<P: MldsaParams>(w1_coeff_values: &[u64]) -> Vec<u64> {
	assert_eq!(w1_coeff_values.len(), P::W1_COEFFICIENTS);

	let mut words = vec![0u64; P::W1_ENCODE_WORDS];
	for (coeff_idx, &coeff) in w1_coeff_values.iter().enumerate() {
		assert!(coeff <= P::W1_COEFF_MAX);
		for bit in 0..P::W1_BITS_PER_COEFF {
			if (coeff >> bit) & 1 == 1 {
				let bit_idx = coeff_idx * P::W1_BITS_PER_COEFF + bit;
				words[bit_idx / 64] |= 1 << (bit_idx % 64);
			}
		}
	}

	words
}

fn pack_hint_bytes_for<P: MldsaParams>(poly_positions: &[Vec<u8>]) -> Vec<u8> {
	assert_eq!(poly_positions.len(), P::K);
	let mut h = vec![0u8; P::HINT_BYTES];
	let mut idx = 0usize;
	for (poly_idx, positions) in poly_positions.iter().enumerate() {
		for &pos in positions {
			assert!(idx < P::OMEGA_USIZE);
			h[idx] = pos;
			idx += 1;
		}
		h[P::OMEGA_USIZE + poly_idx] = idx as u8;
	}
	h
}

fn host_decode_hint_for<P: MldsaParams>(h_bytes: &[u8]) -> Vec<u64> {
	assert_eq!(h_bytes.len(), P::HINT_BYTES);
	let mut h = vec![0u64; P::W1_COEFFICIENTS];
	let mut index = 0usize;
	for poly_idx in 0..P::K {
		let endpoint = h_bytes[P::OMEGA_USIZE + poly_idx] as usize;
		assert!(endpoint >= index);
		assert!(endpoint <= P::OMEGA_USIZE);
		for &pos in &h_bytes[index..endpoint] {
			h[poly_idx * P::N + pos as usize] = 1;
		}
		index = endpoint;
	}
	assert!(h_bytes[index..P::OMEGA_USIZE].iter().all(|&x| x == 0));
	h
}

fn host_high_bits_for<P: MldsaParams>(r: u64) -> u64 {
	assert!(r < P::Q);
	let mut r1 = (r + 127) >> 7;
	if P::GAMMA2 == (P::Q - 1) / 32 {
		r1 = ((r1 * 1025) + (1 << 21)) >> 22;
		r1 & 15
	} else {
		r1 = ((r1 * 11_275) + (1 << 23)) >> 24;
		if r1 > P::W1_COEFF_MAX { 0 } else { r1 }
	}
}

fn host_r0_is_positive_for<P: MldsaParams>(r: u64, r1: u64) -> bool {
	let r1_alpha = r1 * P::TWO_GAMMA2;
	r >= r1_alpha && r != r1_alpha && r - r1_alpha <= (P::Q - 1) / 2
}

fn host_use_hint_for<P: MldsaParams>(h: u64, r: u64) -> u64 {
	assert!(h <= 1);
	let r1 = host_high_bits_for::<P>(r);
	if h == 0 {
		return r1;
	}

	if host_r0_is_positive_for::<P>(r, r1) {
		if r1 == P::W1_COEFF_MAX { 0 } else { r1 + 1 }
	} else if r1 == 0 {
		P::W1_COEFF_MAX
	} else {
		r1 - 1
	}
}

fn host_sample_in_ball_fixed_cap_for<P: MldsaParams>(
	stream: &[u8],
) -> Option<(Vec<u64>, Vec<u64>)> {
	assert_eq!(stream.len(), P::SAMPLE_IN_BALL_STREAM_BYTES);

	let mut signs = u64::from_le_bytes(stream[..8].try_into().unwrap());
	let mut draw_cursor = P::SAMPLE_IN_BALL_SIGN_BYTES;
	let mut coeffs = vec![0u64; P::N];
	let mut draw_counts = vec![0u64; P::TAU];

	for (round, i) in (P::N - P::TAU..P::N).enumerate() {
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

fn build_circuit_for<P: MldsaParams>(kind: MldsaCircuitKind) -> MldsaMeasurementCircuit {
	let builder = CircuitBuilder::new();
	let c_tilde: Vec<_> = (0..P::C_TILDE_BYTES / 8)
		.map(|_| builder.add_witness())
		.collect();
	let z_words: Vec<_> = (0..P::Z_PACKED_WORDS)
		.map(|_| builder.add_witness())
		.collect();
	let h_words: Vec<_> = (0..P::HINT_WORDS).map(|_| builder.add_witness()).collect();
	let h_coeffs: Vec<_> = (0..P::W1_COEFFICIENTS)
		.map(|_| builder.add_witness())
		.collect();
	let mu_words: Vec<_> = (0..P::MU_WORDS).map(|_| builder.add_inout()).collect();
	let w_approx_coeffs: Vec<_> = (0..P::W1_COEFFICIENTS)
		.map(|_| builder.add_witness())
		.collect();

	let relation = match kind {
		MldsaCircuitKind::ExpandedHint => full_bit_heavy_fixed_cap_relation_for::<P>(
			&builder,
			&c_tilde,
			&z_words,
			&h_coeffs,
			&mu_words,
			&w_approx_coeffs,
		),
		MldsaCircuitKind::CanonicalMatchedHint => {
			full_bit_heavy_fixed_cap_canonical_hint_matched_relation_for::<P>(
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

fn sparse_hint_positions_for<P: MldsaParams>() -> Vec<Vec<u8>> {
	let mut positions = Vec::with_capacity(P::K);
	for poly_idx in 0..P::K {
		let count = match poly_idx % 4 {
			0 => 3,
			1 => 2,
			2 => 1,
			_ => 3,
		};
		let mut poly = Vec::with_capacity(count);
		for j in 0..count {
			let pos = ((17 * poly_idx + 43 * j + 5) % P::N) as u8;
			poly.push(pos);
		}
		poly.sort_unstable();
		poly.dedup();
		positions.push(poly);
	}
	assert!(positions.iter().map(Vec::len).sum::<usize>() <= P::OMEGA_USIZE);
	positions
}

fn make_valid_witness_for<P: MldsaParams>() -> MldsaMeasurementWitness {
	let h_bytes = pack_hint_bytes_for::<P>(&sparse_hint_positions_for::<P>());
	let h_coeffs = host_decode_hint_for::<P>(&h_bytes);
	let h_words = words_from_bytes(&h_bytes);

	let mut mu = vec![0u8; P::MU_BYTES];
	for (i, byte) in mu.iter_mut().enumerate() {
		*byte = (0xA5u8).wrapping_add((i as u8).wrapping_mul(17));
	}
	let mu_words = words_from_bytes(&mu);

	let z_packed_y = vec![P::Z_NORM_PACKED_Y_MIN + 123; P::Z_COEFFICIENTS];
	let z_words = pack_z_y_coeffs_for::<P>(&z_packed_y);

	for seed in 0..10_000u64 {
		let w_approx_coeffs: Vec<_> = (0..P::W1_COEFFICIENTS)
			.map(|i| ((i as u64 * 65_537) + seed * 1_299_721 + 12_345) % P::Q)
			.collect();
		let w1: Vec<_> = h_coeffs
			.iter()
			.zip(w_approx_coeffs.iter())
			.map(|(&h, &r)| host_use_hint_for::<P>(h, r))
			.collect();
		let w1_words = pack_w1_coeffs_for::<P>(&w1);

		let mut final_hash_input = Vec::with_capacity(P::FINAL_CHALLENGE_INPUT_BYTES);
		final_hash_input.extend_from_slice(&mu);
		append_words_as_bytes(&mut final_hash_input, &w1_words);
		assert_eq!(final_hash_input.len(), P::FINAL_CHALLENGE_INPUT_BYTES);

		let c_tilde = shake256(&final_hash_input, P::C_TILDE_BYTES);
		let sample_stream = shake256(&c_tilde, P::SAMPLE_IN_BALL_STREAM_BYTES);
		if let Some((_coeffs, draw_counts)) = host_sample_in_ball_fixed_cap_for::<P>(&sample_stream)
		{
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

	panic!("failed to find a fixed-cap SampleInBall witness for {}", P::label());
}

fn populate_witness(
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

fn prove_and_measure_for<P: MldsaParams>(kind: MldsaCircuitKind) -> MldsaProofMetrics {
	let built = timed(|| build_circuit_for::<P>(kind));
	let stat = CircuitStat::collect(&built.value.circuit);
	println!("{}_{} circuit stats:\n{stat}", P::label(), kind.suffix());

	let witness_values = make_valid_witness_for::<P>();
	let witness = timed(|| populate_witness(&built.value, &witness_values));
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

fn stat_canonical_matched_for<P: MldsaParams>() -> CircuitStat {
	let circuit = build_circuit_for::<P>(MldsaCircuitKind::CanonicalMatchedHint).circuit;
	CircuitStat::collect(&circuit)
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

fn duration_ms(duration: Duration) -> f64 {
	duration.as_secs_f64() * 1000.0
}

fn percentile_duration(
	mut values: Vec<Duration>,
	numerator: usize,
	denominator: usize,
) -> Duration {
	assert!(!values.is_empty());
	values.sort_unstable();
	let idx = ((values.len() - 1) * numerator).div_ceil(denominator);
	values[idx]
}

fn print_repeated_metrics(name: &str, metrics: &[MldsaProofMetrics]) {
	assert!(!metrics.is_empty());
	let proof_sizes: Vec<_> = metrics
		.iter()
		.map(|metric| metric.proof_size_bytes)
		.collect();
	let all_same_proof_size = proof_sizes.iter().all(|&size| size == proof_sizes[0]);
	println!("{name} repeated proof metrics over {} runs:", metrics.len());
	for (label, values) in [
		(
			"circuit_build_ms",
			metrics
				.iter()
				.map(|metric| metric.circuit_build)
				.collect::<Vec<_>>(),
		),
		(
			"witness_generation_ms",
			metrics
				.iter()
				.map(|metric| metric.witness_generation)
				.collect::<Vec<_>>(),
		),
		(
			"verifier_setup_ms",
			metrics
				.iter()
				.map(|metric| metric.verifier_setup)
				.collect::<Vec<_>>(),
		),
		(
			"prover_setup_ms",
			metrics
				.iter()
				.map(|metric| metric.prover_setup)
				.collect::<Vec<_>>(),
		),
		(
			"prove_ms",
			metrics
				.iter()
				.map(|metric| metric.prove)
				.collect::<Vec<_>>(),
		),
		(
			"verify_ms",
			metrics
				.iter()
				.map(|metric| metric.verify)
				.collect::<Vec<_>>(),
		),
	] {
		let median = percentile_duration(values.clone(), 1, 2);
		let p95 = percentile_duration(values, 95, 100);
		println!("  {label}: median={:.3}, p95={:.3}", duration_ms(median), duration_ms(p95));
	}
	println!(
		"  proof_size_bytes: {}{}",
		proof_sizes[0],
		if all_same_proof_size {
			""
		} else {
			" (varied across runs)"
		}
	);
}

fn component_stat(name: &'static str, build: impl FnOnce(&CircuitBuilder)) -> ComponentStat {
	let builder = CircuitBuilder::new();
	build(&builder);
	let circuit = builder.build();
	let stat = CircuitStat::collect(&circuit);
	ComponentStat {
		name,
		n_gates: stat.n_gates,
		n_evaluations: stat.n_eval_insn,
		n_and_constraints: stat.n_and_constraints,
		n_mul_constraints: stat.n_mul_constraints,
		n_inout: stat.n_inout,
		n_private: stat.n_witness,
		total_committed: stat.value_vec_len,
		scratch: stat.n_scratch,
	}
}

fn component_stats_for<P: MldsaParams>() -> Vec<ComponentStat> {
	vec![
		component_stat("z_packed_decode_and_norm", |builder| {
			let z_words: Vec<_> = (0..P::Z_PACKED_WORDS)
				.map(|_| builder.add_witness())
				.collect();
			assert_z_packed_bytes_norm_for::<P>(builder, &z_words);
		}),
		component_stat("hint_canonical_matches_expanded", |builder| {
			let h_words: Vec<_> = (0..P::HINT_WORDS).map(|_| builder.add_witness()).collect();
			let h_coeffs: Vec<_> = (0..P::W1_COEFFICIENTS)
				.map(|_| builder.add_witness())
				.collect();
			assert_hint_canonical_matches_expanded_for::<P>(builder, &h_words, &h_coeffs);
		}),
		component_stat("use_hint_no_weight_check", |builder| {
			let h_coeffs: Vec<_> = (0..P::W1_COEFFICIENTS)
				.map(|_| builder.add_witness())
				.collect();
			let w_approx_coeffs: Vec<_> = (0..P::W1_COEFFICIENTS)
				.map(|_| builder.add_witness())
				.collect();
			for (&h, &r) in h_coeffs.iter().zip(w_approx_coeffs.iter()) {
				use_hint_coeff_for::<P>(builder, h, r);
			}
		}),
		component_stat("w1_encode", |builder| {
			let w1_coeffs: Vec<_> = (0..P::W1_COEFFICIENTS)
				.map(|_| builder.add_witness())
				.collect();
			encode_w1_for::<P>(builder, &w1_coeffs);
		}),
		component_stat("final_challenge_shake256", |builder| {
			let input_words: Vec<_> = (0..P::FINAL_CHALLENGE_INPUT_WORDS)
				.map(|_| builder.add_witness())
				.collect();
			final_challenge_hash_for::<P>(builder, &input_words);
		}),
		component_stat("sample_in_ball_shake256", |builder| {
			let c_tilde: Vec<_> = (0..P::C_TILDE_BYTES / 8)
				.map(|_| builder.add_witness())
				.collect();
			sample_in_ball_fixed_cap_stream_for::<P>(builder, &c_tilde);
		}),
		component_stat("sample_in_ball_rejection_dense", |builder| {
			let stream: Vec<_> = (0..P::SAMPLE_IN_BALL_STREAM_BYTES.div_ceil(8))
				.map(|_| builder.add_witness())
				.collect();
			sample_in_ball_fixed_cap_from_stream_for::<P>(builder, &stream);
		}),
		component_stat("sample_in_ball_rejection_sparse", |builder| {
			let stream: Vec<_> = (0..P::SAMPLE_IN_BALL_STREAM_BYTES.div_ceil(8))
				.map(|_| builder.add_witness())
				.collect();
			sample_in_ball_fixed_cap_sparse_from_stream_for::<P>(builder, &stream);
		}),
		component_stat("full_canonical_matched_fixed_cap", |builder| {
			let c_tilde: Vec<_> = (0..P::C_TILDE_BYTES / 8)
				.map(|_| builder.add_witness())
				.collect();
			let z_words: Vec<_> = (0..P::Z_PACKED_WORDS)
				.map(|_| builder.add_witness())
				.collect();
			let h_words: Vec<_> = (0..P::HINT_WORDS).map(|_| builder.add_witness()).collect();
			let h_coeffs: Vec<_> = (0..P::W1_COEFFICIENTS)
				.map(|_| builder.add_witness())
				.collect();
			let mu_words: Vec<_> = (0..P::MU_WORDS).map(|_| builder.add_inout()).collect();
			let w_approx_coeffs: Vec<_> = (0..P::W1_COEFFICIENTS)
				.map(|_| builder.add_witness())
				.collect();
			full_bit_heavy_fixed_cap_canonical_hint_matched_relation_for::<P>(
				builder,
				&c_tilde,
				&z_words,
				&h_words,
				&h_coeffs,
				&mu_words,
				&w_approx_coeffs,
			);
		}),
	]
}

fn print_component_breakdown_for<P: MldsaParams>() {
	let stats = component_stats_for::<P>();
	println!("{} component breakdown:", P::label());
	println!(
		"| component | gates | evaluations | AND | MUL | public | private | committed | scratch |"
	);
	println!("|---|---:|---:|---:|---:|---:|---:|---:|---:|");
	for stat in stats {
		println!(
			"| {} | {} | {} | {} | {} | {} | {} | {} | {} |",
			stat.name,
			stat.n_gates,
			stat.n_evaluations,
			stat.n_and_constraints,
			stat.n_mul_constraints,
			stat.n_inout,
			stat.n_private,
			stat.total_committed,
			stat.scratch,
		);
	}
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
		assert_z_packed_bytes_norm_for::<Mldsa44>(builder, &z_words);
	});

	print_circuit_stats("hint_canonical_matches_expanded", |builder| {
		let h_words: Vec<_> = (0..mldsa44::HINT_WORDS)
			.map(|_| builder.add_witness())
			.collect();
		let h_coeffs: Vec<_> = (0..mldsa44::W1_COEFFICIENTS)
			.map(|_| builder.add_witness())
			.collect();
		assert_hint_canonical_matches_expanded_for::<Mldsa44>(builder, &h_words, &h_coeffs);
	});

	print_circuit_stats("use_hint_no_weight_check", |builder| {
		let h_coeffs: Vec<_> = (0..mldsa44::W1_COEFFICIENTS)
			.map(|_| builder.add_witness())
			.collect();
		let w_approx_coeffs: Vec<_> = (0..mldsa44::W1_COEFFICIENTS)
			.map(|_| builder.add_witness())
			.collect();
		for (&h, &r) in h_coeffs.iter().zip(w_approx_coeffs.iter()) {
			use_hint_coeff_for::<Mldsa44>(builder, h, r);
		}
	});

	print_circuit_stats("w1_encode", |builder| {
		let w1_coeffs: Vec<_> = (0..mldsa44::W1_COEFFICIENTS)
			.map(|_| builder.add_witness())
			.collect();
		encode_w1_for::<Mldsa44>(builder, &w1_coeffs);
	});

	print_circuit_stats("final_challenge_shake256_7f", |builder| {
		let input_words: Vec<_> = (0..mldsa44::FINAL_CHALLENGE_INPUT_WORDS)
			.map(|_| builder.add_witness())
			.collect();
		final_challenge_hash_for::<Mldsa44>(builder, &input_words);
	});

	print_circuit_stats("sample_in_ball_shake256_1f", |builder| {
		let c_tilde: Vec<_> = (0..mldsa44::C_TILDE_BYTES / 8)
			.map(|_| builder.add_witness())
			.collect();
		sample_in_ball_fixed_cap_stream_for::<Mldsa44>(builder, &c_tilde);
	});

	print_circuit_stats("sample_in_ball_rejection_only", |builder| {
		let stream: Vec<_> = (0..mldsa44::SAMPLE_IN_BALL_STREAM_BYTES.div_ceil(8))
			.map(|_| builder.add_witness())
			.collect();
		sample_in_ball_fixed_cap_from_stream_for::<Mldsa44>(builder, &stream);
	});

	print_circuit_stats("sample_in_ball_sparse_rejection_only", |builder| {
		let stream: Vec<_> = (0..mldsa44::SAMPLE_IN_BALL_STREAM_BYTES.div_ceil(8))
			.map(|_| builder.add_witness())
			.collect();
		sample_in_ball_fixed_cap_sparse_from_stream_for::<Mldsa44>(builder, &stream);
	});
}

#[test]
#[ignore = "stats-only circuit component breakdown for all parameter sets"]
fn mldsa_variant_constraint_breakdown() {
	print_component_breakdown_for::<Mldsa44>();
	print_component_breakdown_for::<Mldsa65>();
	print_component_breakdown_for::<Mldsa87>();
}

#[test]
#[ignore = "expensive end-to-end proof measurement"]
fn measure_mldsa44_canonical_matched_proof() {
	let metrics = prove_and_measure_for::<Mldsa44>(MldsaCircuitKind::CanonicalMatchedHint);
	print_metrics("mldsa44_canonical_matched_hint", &metrics);
}

#[test]
#[ignore = "expensive end-to-end proof measurement"]
fn measure_mldsa44_expanded_hint_proof() {
	let metrics = prove_and_measure_for::<Mldsa44>(MldsaCircuitKind::ExpandedHint);
	print_metrics("mldsa44_expanded_hint", &metrics);
}

#[test]
#[ignore = "stats-only full circuit variants"]
fn mldsa_variant_canonical_matched_stats() {
	let stat44 = stat_canonical_matched_for::<Mldsa44>();
	let stat65 = stat_canonical_matched_for::<Mldsa65>();
	let stat87 = stat_canonical_matched_for::<Mldsa87>();

	println!("mldsa44_canonical_matched_hint stats:\n{stat44}");
	println!("mldsa65_canonical_matched_hint stats:\n{stat65}");
	println!("mldsa87_canonical_matched_hint stats:\n{stat87}");

	assert_eq!(stat44.n_inout, mldsa44::MU_WORDS);
	assert_eq!(stat65.n_inout, mldsa65::MU_WORDS);
	assert_eq!(stat87.n_inout, mldsa87::MU_WORDS);
	assert_eq!(stat44.n_mul_constraints, 0);
	assert_eq!(stat65.n_mul_constraints, 0);
	assert_eq!(stat87.n_mul_constraints, 0);
}

#[test]
#[ignore = "expensive end-to-end proof measurement"]
fn measure_mldsa65_canonical_matched_proof() {
	let metrics = prove_and_measure_for::<Mldsa65>(MldsaCircuitKind::CanonicalMatchedHint);
	print_metrics("mldsa65_canonical_matched_hint", &metrics);
}

#[test]
#[ignore = "expensive end-to-end proof measurement"]
fn measure_mldsa87_canonical_matched_proof() {
	let metrics = prove_and_measure_for::<Mldsa87>(MldsaCircuitKind::CanonicalMatchedHint);
	print_metrics("mldsa87_canonical_matched_hint", &metrics);
}

#[test]
#[ignore = "expensive repeated end-to-end proof measurement"]
fn measure_mldsa_variants_canonical_matched_proof_repeated() {
	const RUNS: usize = 7;

	let metrics44: Vec<_> = (0..RUNS)
		.map(|_| prove_and_measure_for::<Mldsa44>(MldsaCircuitKind::CanonicalMatchedHint))
		.collect();
	print_repeated_metrics("mldsa44_canonical_matched_hint", &metrics44);

	let metrics65: Vec<_> = (0..RUNS)
		.map(|_| prove_and_measure_for::<Mldsa65>(MldsaCircuitKind::CanonicalMatchedHint))
		.collect();
	print_repeated_metrics("mldsa65_canonical_matched_hint", &metrics65);

	let metrics87: Vec<_> = (0..RUNS)
		.map(|_| prove_and_measure_for::<Mldsa87>(MldsaCircuitKind::CanonicalMatchedHint))
		.collect();
	print_repeated_metrics("mldsa87_canonical_matched_hint", &metrics87);
}
