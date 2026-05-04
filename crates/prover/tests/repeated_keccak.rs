// Copyright 2026 The Binius Developers

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
const LOG_INSTANCES: usize = 1;
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

fn repeated_keccak_fixture() -> (RepeatedConstraintSystem, ConstraintSystem, ValueVec) {
	let keccak = FixedKeccakCircuit::new(MESSAGE_LEN_BYTES);
	let base_constraint_system = keccak.constraint_system();
	let repeated = RepeatedConstraintSystem::new(base_constraint_system, LOG_INSTANCES);
	let flat_constraint_system = repeated.to_flat_constraint_system();
	let instances = (0..1usize << LOG_INSTANCES)
		.map(|instance| keccak.value_vec(&message_for_instance(instance)))
		.collect::<Vec<_>>();
	let flat_value_vec = repeated
		.to_flat_value_vec(&instances)
		.expect("Keccak instance value vectors flatten into repeated layout");

	(repeated, flat_constraint_system, flat_value_vec)
}

#[test]
fn repeated_keccak_descriptor_flattens_builder_witnesses() {
	let (repeated, flat_constraint_system, flat_value_vec) = repeated_keccak_fixture();

	assert_eq!(repeated.base().value_vec_layout.n_inout, 0);
	assert!(repeated.matches_flat_constraint_system(&flat_constraint_system));
	verify_constraints(&flat_constraint_system, &flat_value_vec)
		.expect("flattened repeated Keccak witness satisfies the expanded circuit");
}

#[test]
#[ignore = "runs an end-to-end repeated Keccak proof"]
fn repeated_keccak_prove_verify() {
	let (repeated, flat_constraint_system, flat_value_vec) = repeated_keccak_fixture();
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
