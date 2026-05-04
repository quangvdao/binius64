// Copyright 2026 The Binius Developers

use std::time::{Duration, Instant};

use binius_circuits::keccak::{N_WORDS_PER_DIGEST, fixed_length::keccak256};
use binius_core::{
	constraint_system::{ConstraintSystem, ValueVec},
	verify::verify_constraints,
	word::Word,
};
use binius_field::arch::OptimalPackedB128;
use binius_frontend::{Circuit, CircuitBuilder, Wire};
use binius_prover::{Prover, hash::parallel_compression::ParallelCompressionAdaptor};
use binius_transcript::ProverTranscript;
use binius_verifier::{
	RepeatedConstraintSystem, Verifier,
	config::StdChallenger,
	hash::{StdCompression, StdDigest},
};
use sha3::{Digest, Keccak256};

const LOG_INV_RATE: usize = 1;
const SMOKE_LOG_INSTANCES: usize = 1;
const MESSAGE_LEN_BYTES: usize = 8;

struct FixedKeccakCircuit {
	circuit: Circuit,
	message: Vec<Wire>,
	digest: [Wire; N_WORDS_PER_DIGEST],
}

impl FixedKeccakCircuit {
	fn new(len_bytes: usize) -> Self {
		let builder = CircuitBuilder::new();
		let message = (0..len_bytes.div_ceil(8))
			.map(|_| builder.add_witness())
			.collect::<Vec<_>>();
		let digest = std::array::from_fn(|_| builder.add_witness());
		let computed_digest = keccak256(&builder, &message, len_bytes);

		for (index, (actual, expected)) in computed_digest.into_iter().zip(digest).enumerate() {
			builder.assert_eq(format!("keccak_digest[{index}]"), actual, expected);
		}

		let circuit = builder.build();
		Self {
			circuit,
			message,
			digest,
		}
	}

	fn constraint_system(&self) -> ConstraintSystem {
		let mut constraint_system = self.circuit.constraint_system().clone();
		constraint_system
			.validate_and_prepare()
			.expect("fixed Keccak constraint system prepares for proving");
		constraint_system
	}

	fn value_vec(&self, message: &[u8]) -> ValueVec {
		let mut witness = self.circuit.new_witness_filler();
		for (wire, chunk) in self.message.iter().zip(message.chunks(8)) {
			let mut word_bytes = [0u8; 8];
			word_bytes[..chunk.len()].copy_from_slice(chunk);
			witness[*wire] = Word::from_u64(u64::from_le_bytes(word_bytes));
		}

		let digest: [u8; 32] = Keccak256::digest(message).into();
		for (wire, chunk) in self.digest.iter().zip(digest.chunks(8)) {
			witness[*wire] = Word::from_u64(u64::from_le_bytes(chunk.try_into().unwrap()));
		}

		self.circuit
			.populate_wire_witness(&mut witness)
			.expect("fixed Keccak witness satisfies the frontend circuit");
		witness.into_value_vec()
	}
}

fn message_for_instance(instance: usize) -> [u8; MESSAGE_LEN_BYTES] {
	let mut message = [0u8; MESSAGE_LEN_BYTES];
	for (index, byte) in message.iter_mut().enumerate() {
		*byte = (0x42 + 17 * instance as u8 + index as u8) ^ ((index as u8) << 3);
	}
	message
}

fn repeated_keccak_fixture(
	log_instances: usize,
) -> (RepeatedConstraintSystem, ConstraintSystem, ValueVec) {
	let keccak = FixedKeccakCircuit::new(MESSAGE_LEN_BYTES);
	let base_constraint_system = keccak.constraint_system();
	let repeated = RepeatedConstraintSystem::new(base_constraint_system, log_instances);
	let flat_constraint_system = repeated.to_flat_constraint_system();
	let instances = (0..1usize << log_instances)
		.map(|instance| keccak.value_vec(&message_for_instance(instance)))
		.collect::<Vec<_>>();
	let flat_value_vec = repeated
		.to_flat_value_vec(&instances)
		.expect("Keccak instance value vectors flatten into repeated layout");

	(repeated, flat_constraint_system, flat_value_vec)
}

fn elapsed_for<T>(f: impl FnOnce() -> T) -> (T, Duration) {
	let start = Instant::now();
	let value = f();
	(value, start.elapsed())
}

fn average_elapsed(iterations: usize, mut f: impl FnMut()) -> Duration {
	let start = Instant::now();
	for _ in 0..iterations {
		f();
	}
	start.elapsed() / iterations as u32
}

#[test]
fn repeated_keccak_descriptor_flattens_builder_witnesses() {
	let (repeated, flat_constraint_system, flat_value_vec) =
		repeated_keccak_fixture(SMOKE_LOG_INSTANCES);

	assert_eq!(repeated.base().value_vec_layout.n_inout, 0);
	assert!(repeated.matches_flat_constraint_system(&flat_constraint_system));
	verify_constraints(&flat_constraint_system, &flat_value_vec)
		.expect("flattened repeated Keccak witness satisfies the expanded circuit");
}

#[test]
#[ignore = "runs an end-to-end repeated Keccak proof"]
fn repeated_keccak_prove_verify() {
	let (repeated, flat_constraint_system, flat_value_vec) =
		repeated_keccak_fixture(SMOKE_LOG_INSTANCES);
	verify_constraints(&flat_constraint_system, &flat_value_vec)
		.expect("flattened repeated Keccak witness satisfies the expanded circuit");

	let verifier = Verifier::<StdDigest, _>::setup_repeated(
		&repeated,
		LOG_INV_RATE,
		StdCompression::default(),
	)
	.expect("repeated verifier setup succeeds");
	let prover = Prover::<OptimalPackedB128, _, StdDigest>::setup_repeated(
		verifier.clone(),
		ParallelCompressionAdaptor::new(StdCompression::default()),
		repeated.clone(),
	)
	.expect("repeated Keccak prover setup succeeds");

	let mut prover_transcript = ProverTranscript::new(StdChallenger::default());
	prover
		.prove_repeated(&repeated, flat_value_vec.clone(), &mut prover_transcript)
		.expect("repeated Keccak prover succeeds");

	let mut verifier_transcript = prover_transcript.into_verifier();
	verifier
		.verify_repeated(flat_value_vec.public(), &repeated, &mut verifier_transcript)
		.expect("repeated Keccak verifier accepts");
	verifier_transcript
		.finalize()
		.expect("repeated Keccak transcript is exhausted");
}

#[test]
#[ignore = "prints flat vs repeated verifier runtimes for repeated Keccak"]
fn repeated_keccak_verifier_print_runtimes() {
	println!(
		"End-to-end Keccak verifier timing with hidden message/digest witness values and structured repeated Shift verification."
	);
	println!(
		"log_instances,instances,base_and_constraints,flat_and_constraints,prove_repeated_ms,flat_verify_ms,repeated_verify_ms,speedup"
	);

	for log_instances in [0usize, 1, 2, 4] {
		let (repeated, flat_constraint_system, flat_value_vec) =
			repeated_keccak_fixture(log_instances);
		verify_constraints(&flat_constraint_system, &flat_value_vec)
			.expect("flattened repeated Keccak witness satisfies the expanded circuit");

		let verifier = Verifier::<StdDigest, _>::setup_repeated(
			&repeated,
			LOG_INV_RATE,
			StdCompression::default(),
		)
		.expect("repeated verifier setup succeeds");
		let flat_prover = Prover::<OptimalPackedB128, _, StdDigest>::setup(
			verifier.clone(),
			ParallelCompressionAdaptor::new(StdCompression::default()),
		)
		.expect("flat Keccak prover setup succeeds");
		let repeated_prover = Prover::<OptimalPackedB128, _, StdDigest>::setup_repeated(
			verifier.clone(),
			ParallelCompressionAdaptor::new(StdCompression::default()),
			repeated.clone(),
		)
		.expect("repeated Keccak prover setup succeeds");

		let flat_prover_transcript = {
			let mut prover_transcript = ProverTranscript::new(StdChallenger::default());
			flat_prover
				.prove(flat_value_vec.clone(), &mut prover_transcript)
				.expect("flat Keccak prover succeeds");
			prover_transcript
		};
		let (repeated_prover_transcript, prove_repeated_elapsed) = elapsed_for(|| {
			let mut prover_transcript = ProverTranscript::new(StdChallenger::default());
			repeated_prover
				.prove_repeated(&repeated, flat_value_vec.clone(), &mut prover_transcript)
				.expect("repeated Keccak prover succeeds");
			prover_transcript
		});

		let verify_iterations = 64;
		let flat_verify_elapsed = average_elapsed(verify_iterations, || {
			let mut verifier_transcript = flat_prover_transcript.clone().into_verifier();
			verifier
				.verify(flat_value_vec.public(), &mut verifier_transcript)
				.expect("flat Keccak verifier accepts");
			verifier_transcript
				.finalize()
				.expect("flat Keccak transcript is exhausted");
		});
		let repeated_verify_elapsed = average_elapsed(verify_iterations, || {
			let mut verifier_transcript = repeated_prover_transcript.clone().into_verifier();
			verifier
				.verify_repeated(flat_value_vec.public(), &repeated, &mut verifier_transcript)
				.expect("repeated Keccak verifier accepts");
			verifier_transcript
				.finalize()
				.expect("repeated Keccak transcript is exhausted");
		});

		let prove_repeated_ms = prove_repeated_elapsed.as_secs_f64() * 1_000.0;
		let flat_verify_ms = flat_verify_elapsed.as_secs_f64() * 1_000.0;
		let repeated_verify_ms = repeated_verify_elapsed.as_secs_f64() * 1_000.0;
		println!(
			"{log_instances},{},{},{},{prove_repeated_ms:.3},{flat_verify_ms:.3},{repeated_verify_ms:.3},{:.2}x",
			1usize << log_instances,
			repeated.base().and_constraints.len(),
			flat_constraint_system.and_constraints.len(),
			flat_verify_ms / repeated_verify_ms,
		);
	}
}

#[test]
#[ignore = "prints flat vs repeated prover setup/key/prove runtimes for repeated Keccak"]
fn repeated_keccak_prover_key_materialization_print_runtimes() {
	println!(
		"Keccak prover timing with flat Shift keys for every instance vs compact repeated Shift keys for the base circuit."
	);
	println!(
		"log_instances,instances,flat_key_words,repeated_key_words,flat_keys,repeated_keys,flat_setup_ms,repeated_setup_ms,flat_repeated_prove_ms,compact_repeated_prove_ms"
	);

	for log_instances in [0usize, 1, 2, 4] {
		let (repeated, flat_constraint_system, flat_value_vec) =
			repeated_keccak_fixture(log_instances);
		verify_constraints(&flat_constraint_system, &flat_value_vec)
			.expect("flattened repeated Keccak witness satisfies the expanded circuit");

		let verifier = Verifier::<StdDigest, _>::setup_repeated(
			&repeated,
			LOG_INV_RATE,
			StdCompression::default(),
		)
		.expect("repeated verifier setup succeeds");

		let (flat_prover, flat_setup_elapsed) = elapsed_for(|| {
			Prover::<OptimalPackedB128, _, StdDigest>::setup(
				verifier.clone(),
				ParallelCompressionAdaptor::new(StdCompression::default()),
			)
			.expect("flat Keccak prover setup succeeds")
		});
		let (repeated_prover, repeated_setup_elapsed) = elapsed_for(|| {
			Prover::<OptimalPackedB128, _, StdDigest>::setup_repeated(
				verifier.clone(),
				ParallelCompressionAdaptor::new(StdCompression::default()),
				repeated.clone(),
			)
			.expect("repeated Keccak prover setup succeeds")
		});

		let (_, flat_repeated_prove_elapsed) = elapsed_for(|| {
			let mut prover_transcript = ProverTranscript::new(StdChallenger::default());
			flat_prover
				.prove_repeated(&repeated, flat_value_vec.clone(), &mut prover_transcript)
				.expect("flat-key repeated Keccak prover succeeds");
			prover_transcript
		});
		let (compact_transcript, compact_repeated_prove_elapsed) = elapsed_for(|| {
			let mut prover_transcript = ProverTranscript::new(StdChallenger::default());
			repeated_prover
				.prove_repeated(&repeated, flat_value_vec.clone(), &mut prover_transcript)
				.expect("compact repeated Keccak prover succeeds");
			prover_transcript
		});

		let mut verifier_transcript = compact_transcript.into_verifier();
		verifier
			.verify_repeated(flat_value_vec.public(), &repeated, &mut verifier_transcript)
			.expect("repeated Keccak verifier accepts compact-key proof");
		verifier_transcript
			.finalize()
			.expect("compact repeated Keccak transcript is exhausted");

		println!(
			"{log_instances},{},{},{},{},{},{:.3},{:.3},{:.3},{:.3}",
			1usize << log_instances,
			flat_prover.key_collection().key_ranges.len(),
			repeated_prover.key_collection().key_ranges.len(),
			flat_prover.key_collection().keys.len(),
			repeated_prover.key_collection().keys.len(),
			flat_setup_elapsed.as_secs_f64() * 1_000.0,
			repeated_setup_elapsed.as_secs_f64() * 1_000.0,
			flat_repeated_prove_elapsed.as_secs_f64() * 1_000.0,
			compact_repeated_prove_elapsed.as_secs_f64() * 1_000.0,
		);
	}
}
