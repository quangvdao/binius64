// Copyright 2025 Irreducible Inc.

use binius_core::{
	constraint_system::{AndConstraint, ConstraintSystem, MulConstraint, ValueVec},
	verify::eval_operand,
	word::Word,
};
use binius_field::{
	AESTowerField8b as B8, BinaryField, ExtensionField, PackedAESBinaryField16x8b, PackedExtension,
	PackedField, UnderlierWithBitOps, WithUnderlier,
};
use binius_iop_prover::{
	basefold_channel::BaseFoldProverChannel, basefold_compiler::BaseFoldProverCompiler,
	channel::IOPProverChannel,
};
use binius_math::{
	BinarySubspace, FieldBuffer, FieldSlice,
	inner_product::inner_product,
	multilinear::{eq::eq_ind_partial_eval, evaluate::evaluate},
	ntt::{NeighborsLastMultiThread, domain_context::GenericPreExpanded},
	univariate::lagrange_evals,
};
use binius_transcript::{ProverTranscript, fiat_shamir::Challenger};
use binius_utils::{SerializeBytes, checked_arithmetics::checked_log_2, rayon::prelude::*};
use binius_verifier::{
	IOPVerifier, RepeatedConstraintSystem, Verifier,
	config::{
		B1, B128, LOG_WORD_SIZE_BITS, LOG_WORDS_PER_ELEM, PROVER_SMALL_FIELD_ZEROCHECK_CHALLENGES,
	},
	protocols::{bitand::AndCheckOutput, intmul::IntMulOutput, sumcheck::SumcheckOutput},
};
use digest::{Digest, FixedOutputReset, Output, block_api::BlockSizeUser};

use super::error::Error;
use crate::{
	and_reduction::prover::OblongZerocheckProver,
	hash::{ParallelDigest, parallel_compression::ParallelPseudoCompression},
	merkle_tree::prover::BinaryMerkleTreeProver,
	protocols::{
		intmul::{prove::IntMulProver, witness::Witness as IntMulWitness},
		shift::{
			KeyCollection, OperatorData, ShiftKeySource, build_key_collection,
			prove_with_key_source,
		},
	},
	ring_switch,
};

/// Type alias for the prover NTT parameterized by field.
type ProverNTT<F> = NeighborsLastMultiThread<GenericPreExpanded<F>>;

/// Type alias for the prover Merkle tree prover parameterized by field.
type ProverMerkleProver<F, ParallelMerkleHasher, ParallelMerkleCompress> =
	BinaryMerkleTreeProver<F, ParallelMerkleHasher, ParallelMerkleCompress>;

/// IOP prover for a particular constraint system.
///
/// This struct encapsulates the constraint system and pre-computed keys,
/// providing the core proving logic independent of the specific IOP compilation strategy.
/// Most users should use [`Prover`] instead, which wraps this with a BaseFold compiler.
#[derive(Debug)]
pub struct IOPProver {
	constraint_system: ConstraintSystem,
	log_public_words: usize,
	log_witness_elems: usize,
	shift_keys: ShiftKeyMaterialization,
}

#[derive(Debug)]
enum ShiftKeyMaterialization {
	Flat(KeyCollection),
	Repeated {
		descriptor: RepeatedConstraintSystem,
		base_key_collection: KeyCollection,
	},
}

impl ShiftKeyMaterialization {
	fn key_collection(&self) -> &KeyCollection {
		match self {
			Self::Flat(key_collection) => key_collection,
			Self::Repeated {
				base_key_collection,
				..
			} => base_key_collection,
		}
	}
}

impl IOPProver {
	/// Constructs an IOP prover from an IOP verifier and pre-computed keys.
	pub fn new(iop_verifier: IOPVerifier, key_collection: KeyCollection) -> Self {
		let log_public_words = iop_verifier.log_public_words();
		let log_witness_elems = iop_verifier.log_witness_elems();
		let constraint_system = iop_verifier.into_constraint_system();
		Self {
			constraint_system,
			log_public_words,
			log_witness_elems,
			shift_keys: ShiftKeyMaterialization::Flat(key_collection),
		}
	}

	/// Constructs an IOP prover for a flat repeated circuit using only the base key collection.
	pub fn new_repeated(
		iop_verifier: IOPVerifier,
		repeated: RepeatedConstraintSystem,
		base_key_collection: KeyCollection,
	) -> Result<Self, Error> {
		let log_public_words = iop_verifier.log_public_words();
		let log_witness_elems = iop_verifier.log_witness_elems();
		let constraint_system = iop_verifier.into_constraint_system();
		if !repeated.matches_flat_constraint_system(&constraint_system) {
			return Err(Error::ArgumentError {
				arg: "repeated".to_string(),
				msg: "repeated descriptor does not match the flat constraint system".to_string(),
			});
		}
		if base_key_collection.key_ranges.len()
			!= repeated.base().value_vec_layout.committed_total_len
		{
			return Err(Error::ArgumentError {
				arg: "base_key_collection".to_string(),
				msg: "base key collection does not match repeated base value count".to_string(),
			});
		}
		Ok(Self {
			constraint_system,
			log_public_words,
			log_witness_elems,
			shift_keys: ShiftKeyMaterialization::Repeated {
				descriptor: repeated,
				base_key_collection,
			},
		})
	}

	/// Returns the constraint system.
	pub fn constraint_system(&self) -> &ConstraintSystem {
		&self.constraint_system
	}

	/// Returns a reference to the KeyCollection.
	///
	/// This can be used to serialize the KeyCollection for later use.
	pub fn key_collection(&self) -> &KeyCollection {
		self.shift_keys.key_collection()
	}

	/// Proves using an IOP channel interface.
	///
	/// This is the core proving logic, independent of the specific IOP compilation strategy.
	/// For most users, [`Prover::prove`] is the simpler interface.
	pub fn prove<P, Channel>(&self, witness: ValueVec, channel: Channel) -> Result<(), Error>
	where
		P: PackedField<Scalar = B128>
			+ PackedExtension<B128>
			+ PackedExtension<B1>
			+ WithUnderlier<Underlier: UnderlierWithBitOps>,
		Channel: IOPProverChannel<P>,
	{
		let shift_key_source = match &self.shift_keys {
			ShiftKeyMaterialization::Flat(key_collection) => ShiftKeySource::flat(key_collection),
			ShiftKeyMaterialization::Repeated { .. } => {
				return Err(Error::ArgumentError {
					arg: "prover".to_string(),
					msg: "prover was set up with repeated Shift keys; use prove_repeated"
						.to_string(),
				});
			}
		};
		self.prove_with_shift_keys::<P, _>(witness, channel, shift_key_source)
	}

	/// Proves a repeated circuit using the compact repeated Shift key source.
	pub fn prove_repeated<P, Channel>(
		&self,
		repeated: &RepeatedConstraintSystem,
		witness: ValueVec,
		mut channel: Channel,
	) -> Result<(), Error>
	where
		P: PackedField<Scalar = B128>
			+ PackedExtension<B128>
			+ PackedExtension<B1>
			+ WithUnderlier<Underlier: UnderlierWithBitOps>,
		Channel: IOPProverChannel<P>,
	{
		if !repeated.matches_flat_shape(&self.constraint_system) {
			return Err(Error::ArgumentError {
				arg: "repeated".to_string(),
				msg: "repeated descriptor does not match the flat constraint system shape"
					.to_string(),
			});
		}

		let shift_key_source = match &self.shift_keys {
			ShiftKeyMaterialization::Flat(key_collection) => ShiftKeySource::flat(key_collection),
			ShiftKeyMaterialization::Repeated {
				descriptor,
				base_key_collection,
			} => {
				if descriptor.binding_scalars() != repeated.binding_scalars() {
					return Err(Error::ArgumentError {
						arg: "repeated".to_string(),
						msg: "repeated descriptor differs from prover setup".to_string(),
					});
				}
				ShiftKeySource::repeated(base_key_collection, repeated)
			}
		};

		channel.observe_many(&repeated.binding_scalars());
		self.prove_with_shift_keys::<P, _>(witness, channel, shift_key_source)
	}

	fn prove_with_shift_keys<P, Channel>(
		&self,
		witness: ValueVec,
		mut channel: Channel,
		shift_key_source: ShiftKeySource<'_>,
	) -> Result<(), Error>
	where
		P: PackedField<Scalar = B128>
			+ PackedExtension<B128>
			+ PackedExtension<B1>
			+ WithUnderlier<Underlier: UnderlierWithBitOps>,
		Channel: IOPProverChannel<P>,
	{
		let cs = &self.constraint_system;

		let _prove_guard = tracing::info_span!(
			"Prove",
			operation = "prove",
			perfetto_category = "operation",
			n_witness_words = cs.value_vec_layout.committed_total_len,
			n_bitand = cs.and_constraints.len(),
			n_intmul = cs.mul_constraints.len(),
		)
		.entered();

		// [phase] Setup - initialization and constraint system setup
		let setup_guard =
			tracing::info_span!("[phase] Setup", phase = "setup", perfetto_category = "phase")
				.entered();
		let witness_packed = pack_witness::<P>(self.log_witness_elems, &witness)?;
		drop(setup_guard);

		// Observe the public input as B128 elements (includes it in Fiat-Shamir).
		let n_public_elems = 1 << (self.log_public_words - LOG_WORDS_PER_ELEM);
		let public_elems = witness_packed
			.iter_scalars()
			.take(n_public_elems)
			.collect::<Vec<_>>();
		channel.observe_many(&public_elems);

		// [phase] Witness Commit - witness generation and commitment
		let witness_commit_guard = tracing::info_span!(
			"[phase] Witness Commit",
			phase = "witness_commit",
			perfetto_category = "phase"
		)
		.entered();

		// Commit witness via channel
		let trace_oracle = channel.send_oracle(witness_packed.to_ref());

		drop(witness_commit_guard);

		// [phase] IntMul Reduction - multiplication constraint reduction
		let intmul_guard = tracing::info_span!(
			"[phase] IntMul Reduction",
			phase = "intmul_reduction",
			perfetto_category = "phase",
			n_constraints = cs.mul_constraints.len()
		)
		.entered();
		let mul_witness = build_intmul_witness(&cs.mul_constraints, &witness);
		let intmul_output = prove_intmul_reduction::<_, P, _>(mul_witness, &mut channel)?;
		drop(intmul_guard);

		// [phase] BitAnd Reduction - AND constraint reduction
		let bitand_guard = tracing::info_span!(
			"[phase] BitAnd Reduction",
			phase = "bitand_reduction",
			perfetto_category = "phase",
			n_constraints = cs.and_constraints.len()
		)
		.entered();
		let bitand_claim = {
			let bitand_witness = build_bitand_witness(&cs.and_constraints, &witness);
			let AndCheckOutput {
				a_eval,
				b_eval,
				c_eval,
				z_challenge,
				eval_point,
			} = prove_bitand_reduction::<B128, _>(bitand_witness, &mut channel)?;
			OperatorData {
				evals: vec![a_eval, b_eval, c_eval],
				r_zhat_prime: z_challenge,
				r_x_prime: eval_point,
			}
		};
		drop(bitand_guard);

		// Build `OperatorData` for IntMul using the same `r_zhat_prime`
		// challenge as in BitAnd. Sharing this univariate challenge
		// improves ShiftReduction perf.
		let intmul_claim = {
			let IntMulOutput {
				eval_point,
				a_evals,
				b_evals,
				c_lo_evals,
				c_hi_evals,
			} = intmul_output;

			let r_zhat_prime = bitand_claim.r_zhat_prime;
			let subspace = BinarySubspace::<B8>::with_dim(LOG_WORD_SIZE_BITS).isomorphic();
			let l_tilde = lagrange_evals(&subspace, r_zhat_prime);
			let make_final_claim = |evals| inner_product(evals, l_tilde.iter_scalars());
			OperatorData {
				evals: vec![
					make_final_claim(a_evals),
					make_final_claim(b_evals),
					make_final_claim(c_lo_evals),
					make_final_claim(c_hi_evals),
				],
				r_zhat_prime,
				r_x_prime: eval_point,
			}
		};

		// [phase] Shift Reduction - shift operations
		let shift_guard = tracing::info_span!(
			"[phase] Shift Reduction",
			phase = "shift_reduction",
			perfetto_category = "phase"
		)
		.entered();
		let SumcheckOutput {
			challenges: eval_point,
			eval: _,
		} = prove_with_key_source::<_, P, _>(
			&shift_key_source,
			witness.combined_witness(),
			bitand_claim,
			intmul_claim,
			&mut channel,
		)?;
		drop(shift_guard);

		// [phase] Ring-Switching + PCS Opening
		let pcs_guard = tracing::info_span!(
			"[phase] PCS Opening",
			phase = "pcs_opening",
			perfetto_category = "phase"
		)
		.entered();

		// Ring-switching reduction
		let ring_switch::RingSwitchOutput {
			rs_eq_ind,
			sumcheck_claim,
		} = ring_switch::prove(&witness_packed, &eval_point, &mut channel);

		// Public input check batched with ring-switch
		let log_packing = <B128 as ExtensionField<B1>>::LOG_DEGREE;

		let log_public_elems = self.log_public_words - LOG_WORDS_PER_ELEM;
		let pubcheck_point = &eval_point[log_packing..][..log_public_elems];
		let pubcheck_claim = {
			let public_elems_buf = FieldSlice::from_slice(log_public_elems, &public_elems);
			evaluate(&public_elems_buf, pubcheck_point)
		};

		let batch_coeff: B128 = channel.sample();
		let batched_claim = sumcheck_claim + batch_coeff * pubcheck_claim;

		// Batch the pubcheck transparent with the ring-switch transparent
		let batched_transparent =
			compute_batched_transparent(rs_eq_ind, pubcheck_point, batch_coeff);

		// Prove oracle relations via channel (runs BaseFold internally)
		channel.prove_oracle_relations([(
			trace_oracle,
			witness_packed,
			batched_transparent,
			batched_claim,
		)]);

		drop(pcs_guard);

		Ok(())
	}
}

/// Struct for proving instances of a particular constraint system.
///
/// The [`Self::setup`] constructor pre-processes reusable structures for proving instances of the
/// given constraint system. Then [`Self::prove`] is called one or more times with individual
/// instances.
pub struct Prover<P, ParallelMerkleCompress, ParallelMerkleHasher>
where
	P: PackedField<Scalar = B128>,
	ParallelMerkleHasher: ParallelDigest,
	ParallelMerkleHasher::Digest: Digest + BlockSizeUser + FixedOutputReset,
	ParallelMerkleCompress: ParallelPseudoCompression<Output<ParallelMerkleHasher::Digest>, 2>,
{
	iop_prover: IOPProver,
	#[allow(clippy::type_complexity)]
	basefold_compiler: BaseFoldProverCompiler<
		P,
		ProverNTT<B128>,
		ProverMerkleProver<B128, ParallelMerkleHasher, ParallelMerkleCompress>,
	>,
}

impl<P, MerkleHash, ParallelMerkleCompress, ParallelMerkleHasher>
	Prover<P, ParallelMerkleCompress, ParallelMerkleHasher>
where
	P: PackedField<Scalar = B128>
		+ PackedExtension<B128>
		+ PackedExtension<B1>
		+ WithUnderlier<Underlier: UnderlierWithBitOps>,
	MerkleHash: Digest + BlockSizeUser + FixedOutputReset,
	ParallelMerkleHasher: ParallelDigest<Digest = MerkleHash>,
	ParallelMerkleCompress: ParallelPseudoCompression<Output<MerkleHash>, 2>,
	Output<MerkleHash>: SerializeBytes,
{
	/// Constructs a prover corresponding to a constraint system verifier.
	///
	/// See [`Prover`] struct documentation for details.
	pub fn setup(
		verifier: Verifier<MerkleHash, ParallelMerkleCompress::Compression>,
		compression: ParallelMerkleCompress,
	) -> Result<Self, Error> {
		let key_collection = build_key_collection(verifier.constraint_system());
		Self::setup_with_key_collection(verifier, compression, key_collection)
	}

	/// Constructs a prover for a flat repeated circuit while materializing Shift keys only for
	/// the base circuit.
	pub fn setup_repeated(
		verifier: Verifier<MerkleHash, ParallelMerkleCompress::Compression>,
		compression: ParallelMerkleCompress,
		repeated: RepeatedConstraintSystem,
	) -> Result<Self, Error> {
		let base_key_collection = build_key_collection(repeated.base());
		Self::setup_repeated_with_key_collection(
			verifier,
			compression,
			repeated,
			base_key_collection,
		)
	}

	/// Constructs a prover with a pre-built KeyCollection.
	///
	/// This allows loading a previously serialized KeyCollection to avoid
	/// the expensive key building phase during setup.
	pub fn setup_with_key_collection(
		verifier: Verifier<MerkleHash, ParallelMerkleCompress::Compression>,
		compression: ParallelMerkleCompress,
		key_collection: KeyCollection,
	) -> Result<Self, Error> {
		// Get max subspace from verifier's IOP compiler (reuses FRI params)
		let subspace = verifier.iop_compiler().max_subspace();
		let domain_context = GenericPreExpanded::generate_from_subspace(subspace);
		// FIXME TODO For mobile phones, the number of shares should potentially be more than the
		// number of threads, because the threads/cores have different performance (but in the NTT
		// each share has the same amount of work)
		let log_num_shares = binius_utils::rayon::current_num_threads().ilog2() as usize;
		let ntt = NeighborsLastMultiThread::new(domain_context, log_num_shares);

		let merkle_prover = BinaryMerkleTreeProver::<_, ParallelMerkleHasher, _>::new(compression);

		// Create prover compiler from verifier compiler (reuses FRI params and oracle specs)
		let basefold_compiler = BaseFoldProverCompiler::from_verifier_compiler(
			verifier.iop_compiler(),
			ntt,
			merkle_prover,
		);

		let iop_prover = IOPProver::new(verifier.into_iop_verifier(), key_collection);

		Ok(Prover {
			iop_prover,
			basefold_compiler,
		})
	}

	/// Constructs a repeated prover with a pre-built base KeyCollection.
	pub fn setup_repeated_with_key_collection(
		verifier: Verifier<MerkleHash, ParallelMerkleCompress::Compression>,
		compression: ParallelMerkleCompress,
		repeated: RepeatedConstraintSystem,
		base_key_collection: KeyCollection,
	) -> Result<Self, Error> {
		// Get max subspace from verifier's IOP compiler (reuses FRI params)
		let subspace = verifier.iop_compiler().max_subspace();
		let domain_context = GenericPreExpanded::generate_from_subspace(subspace);
		let log_num_shares = binius_utils::rayon::current_num_threads().ilog2() as usize;
		let ntt = NeighborsLastMultiThread::new(domain_context, log_num_shares);

		let merkle_prover = BinaryMerkleTreeProver::<_, ParallelMerkleHasher, _>::new(compression);

		let basefold_compiler = BaseFoldProverCompiler::from_verifier_compiler(
			verifier.iop_compiler(),
			ntt,
			merkle_prover,
		);

		let iop_prover =
			IOPProver::new_repeated(verifier.into_iop_verifier(), repeated, base_key_collection)?;

		Ok(Prover {
			iop_prover,
			basefold_compiler,
		})
	}

	/// Returns a reference to the IOP prover.
	pub fn iop_prover(&self) -> &IOPProver {
		&self.iop_prover
	}

	/// Returns a reference to the KeyCollection.
	///
	/// This can be used to serialize the KeyCollection for later use.
	pub fn key_collection(&self) -> &KeyCollection {
		self.iop_prover.key_collection()
	}

	pub fn prove<Challenger_: Challenger>(
		&self,
		witness: ValueVec,
		transcript: &mut ProverTranscript<Challenger_>,
	) -> Result<(), Error> {
		// Create channel and delegate to IOPProver::prove
		let channel = BaseFoldProverChannel::from_compiler(&self.basefold_compiler, transcript);
		self.iop_prover.prove::<P, _>(witness, channel)
	}

	/// Proves a flat repeated-circuit witness while binding the repeated descriptor into the
	/// transcript before public inputs.
	pub fn prove_repeated<Challenger_: Challenger>(
		&self,
		repeated: &RepeatedConstraintSystem,
		witness: ValueVec,
		transcript: &mut ProverTranscript<Challenger_>,
	) -> Result<(), Error> {
		let channel = BaseFoldProverChannel::from_compiler(&self.basefold_compiler, transcript);
		self.iop_prover
			.prove_repeated::<P, _>(repeated, witness, channel)
	}
}

/// Batches the pubcheck transparent polynomial with the ring-switch equality indicator.
///
/// Computes `rs_eq_ind + batch_coeff * eq(pubcheck_point || 0, ·)`, adding the scaled
/// pubcheck equality indicator to the first `2^log_public_elems` entries of `rs_eq_ind`.
fn compute_batched_transparent<P: PackedField<Scalar = B128>>(
	mut rs_eq_ind: FieldBuffer<P>,
	pubcheck_point: &[B128],
	batch_coeff: B128,
) -> FieldBuffer<P> {
	let log_public_elems = pubcheck_point.len();
	let pubcheck_eq = eq_ind_partial_eval::<P>(pubcheck_point);
	let mut chunk = rs_eq_ind.chunk_mut(log_public_elems, 0);
	let mut chunk_data = chunk.get();
	let batch = P::broadcast(batch_coeff);
	for (dst, src) in std::iter::zip(chunk_data.as_mut(), pubcheck_eq.as_ref()) {
		*dst += *src * batch;
	}
	drop(chunk);
	rs_eq_ind
}

fn pack_witness<P: PackedField<Scalar = B128>>(
	log_witness_elems: usize,
	witness: &ValueVec,
) -> Result<FieldBuffer<P>, Error> {
	// The number of field elements that constitute the packed witness.
	let n_witness_elems = witness.size().div_ceil(1 << LOG_WORDS_PER_ELEM);
	if n_witness_elems > 1 << log_witness_elems {
		return Err(Error::ArgumentError {
			arg: "witness".to_string(),
			msg: "witness element count is incompatible with the constraint system".to_string(),
		});
	}

	let len = 1 << log_witness_elems.saturating_sub(P::LOG_WIDTH);
	let mut padded_witness_elems = Vec::<P>::with_capacity(len);

	let combined_witness = witness.combined_witness();
	padded_witness_elems
		.spare_capacity_mut()
		.into_par_iter()
		.enumerate()
		.for_each(|(i, dst)| {
			// Pack B128 elements into packed elements
			let offset = i << (P::LOG_WIDTH + 1);
			let value = P::from_fn(|j| {
				let word_0 = combined_witness[offset + 2 * j];
				let word_1 = combined_witness[offset + 2 * j + 1];
				B128::new(((word_1.0 as u128) << 64) | (word_0.0 as u128))
			});

			dst.write(value);
		});

	// SAFETY: We just initialized all elements
	unsafe {
		padded_witness_elems.set_len(len);
	};

	let padded_witness_elems =
		FieldBuffer::new(log_witness_elems, padded_witness_elems.into_boxed_slice());
	Ok(padded_witness_elems)
}

fn prove_bitand_reduction<F, Channel>(
	witness: AndCheckWitness,
	channel: &mut Channel,
) -> Result<AndCheckOutput<F>, Error>
where
	F: BinaryField + From<B8>,
	Channel: binius_ip_prover::channel::IPProverChannel<F>,
{
	let prover_message_domain = BinarySubspace::<B8>::with_dim(LOG_WORD_SIZE_BITS + 1);
	let AndCheckWitness { a, b, c } = witness;

	let log_constraint_count = checked_log_2(a.len());

	let mut small_field_zerocheck_challenges = PROVER_SMALL_FIELD_ZEROCHECK_CHALLENGES.to_vec();
	small_field_zerocheck_challenges.truncate(log_constraint_count);

	let big_field_zerocheck_challenges =
		channel.sample_many(log_constraint_count - small_field_zerocheck_challenges.len());

	let prover = OblongZerocheckProver::<_, PackedAESBinaryField16x8b>::new(
		a,
		b,
		c,
		big_field_zerocheck_challenges,
		small_field_zerocheck_challenges,
		prover_message_domain.isomorphic(),
	);

	Ok(prover.prove_with_channel(channel)?)
}

fn prove_intmul_reduction<F, P, Channel>(
	witness: MulCheckWitness,
	channel: &mut Channel,
) -> Result<IntMulOutput<F>, Error>
where
	F: BinaryField,
	P: PackedField<Scalar = F>,
	Channel: binius_ip_prover::channel::IPProverChannel<F>,
{
	let MulCheckWitness { a, b, lo, hi } = witness;

	let mut mulcheck_prover = IntMulProver::new(0, channel);

	// Words must be converted to u64 because
	// `Bitwise` requires `From<u8>` and `Shr<usize>`
	// We could implement these for `Word` in the future.
	let convert_to_u64 = |w: Vec<Word>| w.into_iter().map(|w| w.0).collect::<Vec<u64>>();
	let [a_u64, b_u64, lo_u64, hi_u64] = [a, b, lo, hi].map(convert_to_u64);
	let intmul_witness =
		IntMulWitness::<P, _, _>::new(LOG_WORD_SIZE_BITS, &a_u64, &b_u64, &lo_u64, &hi_u64)?;

	Ok(mulcheck_prover.prove(intmul_witness)?)
}

struct AndCheckWitness {
	a: Vec<Word>,
	b: Vec<Word>,
	c: Vec<Word>,
}

struct MulCheckWitness {
	a: Vec<Word>,
	b: Vec<Word>,
	lo: Vec<Word>,
	hi: Vec<Word>,
}

#[tracing::instrument(skip_all, "Build BitAnd witness", level = "debug")]
fn build_bitand_witness(and_constraints: &[AndConstraint], witness: &ValueVec) -> AndCheckWitness {
	let n_constraints = and_constraints.len();

	let mut a = Vec::with_capacity(n_constraints);
	let mut b = Vec::with_capacity(n_constraints);
	let mut c = Vec::with_capacity(n_constraints);

	(and_constraints, a.spare_capacity_mut(), b.spare_capacity_mut(), c.spare_capacity_mut())
		.into_par_iter()
		.for_each(|(constraint, a_i, b_i, c_i)| {
			a_i.write(eval_operand(witness, &constraint.a));
			b_i.write(eval_operand(witness, &constraint.b));
			c_i.write(eval_operand(witness, &constraint.c));
		});

	// Safety: all entries in a, b, c are initialized in the parallel loop above.
	unsafe {
		a.set_len(n_constraints);
		b.set_len(n_constraints);
		c.set_len(n_constraints);
	}

	AndCheckWitness { a, b, c }
}

#[tracing::instrument(skip_all, "Build IntMul witness", level = "debug")]
fn build_intmul_witness(mul_constraints: &[MulConstraint], witness: &ValueVec) -> MulCheckWitness {
	let n_constraints = mul_constraints.len();

	let mut a = Vec::with_capacity(n_constraints);
	let mut b = Vec::with_capacity(n_constraints);
	let mut lo = Vec::with_capacity(n_constraints);
	let mut hi = Vec::with_capacity(n_constraints);

	(
		mul_constraints,
		a.spare_capacity_mut(),
		b.spare_capacity_mut(),
		lo.spare_capacity_mut(),
		hi.spare_capacity_mut(),
	)
		.into_par_iter()
		.for_each(|(constraint, a_i, b_i, lo_i, hi_i)| {
			a_i.write(eval_operand(witness, &constraint.a));
			b_i.write(eval_operand(witness, &constraint.b));
			lo_i.write(eval_operand(witness, &constraint.lo));
			hi_i.write(eval_operand(witness, &constraint.hi));
		});

	// Safety: all entries in a, b, lo, hi are initialized in the parallel loop above.
	unsafe {
		a.set_len(n_constraints);
		b.set_len(n_constraints);
		lo.set_len(n_constraints);
		hi.set_len(n_constraints);
	}

	MulCheckWitness { a, b, lo, hi }
}

#[cfg(test)]
mod tests {
	use std::time::{Duration, Instant};

	use binius_core::{
		ShiftVariant,
		constraint_system::{
			AndConstraint, ConstraintSystem, ShiftedValueIndex, ValueIndex, ValueVec,
			ValueVecLayout,
		},
		verify::verify_constraints,
		word::Word,
	};
	use binius_field::arch::OptimalPackedB128;
	use binius_transcript::ProverTranscript;
	use binius_verifier::{
		RepeatedConstraintSystem, Verifier,
		config::StdChallenger,
		hash::{StdCompression, StdDigest},
	};

	use crate::{Prover, hash::parallel_compression::ParallelCompressionAdaptor};

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

	fn test_term(value_index: usize, seed: usize) -> ShiftedValueIndex {
		let value_index = ValueIndex(value_index as u32);
		let amount = (seed * 7 + 3) % 64;
		if amount == 0 {
			return ShiftedValueIndex::plain(value_index);
		}
		let shift_variant = match seed % 4 {
			0 => ShiftVariant::Sll,
			1 => ShiftVariant::Slr,
			2 => ShiftVariant::Sar,
			_ => ShiftVariant::Rotr,
		};
		ShiftedValueIndex {
			value_index,
			shift_variant,
			amount,
		}
	}

	fn value_vec_layout(value_count: usize) -> ValueVecLayout {
		assert!(value_count.is_power_of_two());
		assert!(value_count >= 2);
		ValueVecLayout {
			n_const: 0,
			n_inout: 2,
			n_witness: value_count - 2,
			n_internal: 0,
			offset_inout: 0,
			offset_witness: 2,
			committed_total_len: value_count,
			n_scratch: 0,
		}
	}

	fn make_base_constraint_system(
		base_constraint_count: usize,
		base_value_count: usize,
	) -> ConstraintSystem {
		let and_constraints = (0..base_constraint_count)
			.map(|row| AndConstraint {
				a: vec![
					test_term((row * 3 + 1) % base_value_count, row),
					test_term((row * 5 + 7) % base_value_count, row + 1),
				],
				b: vec![
					test_term((row * 11 + 13) % base_value_count, row + 2),
					test_term((row * 17 + 19) % base_value_count, row + 3),
				],
				c: vec![
					test_term((row * 23 + 29) % base_value_count, row + 4),
					test_term((row * 31 + 37) % base_value_count, row + 5),
				],
			})
			.collect();

		let mut constraint_system = ConstraintSystem::new(
			Vec::new(),
			value_vec_layout(base_value_count),
			and_constraints,
			Vec::new(),
		);
		constraint_system
			.validate_and_prepare()
			.expect("constructed base constraint system is valid");
		constraint_system
	}

	fn make_base_constraint_system_with_shared_constant(
		base_constraint_count: usize,
		base_value_count: usize,
	) -> ConstraintSystem {
		assert!(base_value_count >= 4);
		let layout = ValueVecLayout {
			n_const: 1,
			n_inout: 1,
			n_witness: base_value_count - 2,
			n_internal: 0,
			offset_inout: 1,
			offset_witness: 2,
			committed_total_len: base_value_count,
			n_scratch: 0,
		};
		let and_constraints = (0..base_constraint_count)
			.map(|row| {
				let witness = test_term(2 + (row * 3 % (base_value_count - 2)), row);
				AndConstraint {
					a: vec![ShiftedValueIndex::plain(ValueIndex(0))],
					b: vec![witness],
					c: vec![witness],
				}
			})
			.collect();

		let mut constraint_system =
			ConstraintSystem::new(vec![Word::ALL_ONE], layout, and_constraints, Vec::new());
		constraint_system
			.validate_and_prepare()
			.expect("constructed base constraint system with constants is valid");
		constraint_system
	}

	fn zero_value_vec(constraint_system: &ConstraintSystem) -> ValueVec {
		ValueVec::new_from_data(
			constraint_system.value_vec_layout.clone(),
			vec![Word::ZERO; constraint_system.value_vec_layout.offset_witness],
			vec![
				Word::ZERO;
				constraint_system.value_vec_layout.committed_total_len
					- constraint_system.value_vec_layout.offset_witness
			],
		)
		.expect("zero value vector has matching layout")
	}

	fn setup_repeated_fixture(
		base_constraint_count: usize,
		base_value_count: usize,
		log_instances: usize,
	) -> (
		RepeatedConstraintSystem,
		ConstraintSystem,
		ValueVec,
		Verifier<StdDigest, StdCompression>,
		Prover<OptimalPackedB128, ParallelCompressionAdaptor<StdCompression>, StdDigest>,
	) {
		const LOG_INV_RATE: usize = 1;

		let base_constraint_system =
			make_base_constraint_system(base_constraint_count, base_value_count);
		let repeated = RepeatedConstraintSystem::new(base_constraint_system.clone(), log_instances);
		let flat_constraint_system = repeated.to_flat_constraint_system();
		let value_vec = zero_value_vec(&flat_constraint_system);
		verify_constraints(&flat_constraint_system, &value_vec)
			.expect("zero witness satisfies the repeated toy circuit");

		let verifier = Verifier::<StdDigest, _>::setup_repeated(
			&repeated,
			LOG_INV_RATE,
			StdCompression::default(),
		)
		.expect("repeated verifier setup succeeds");
		let prover = Prover::<OptimalPackedB128, _, StdDigest>::setup(
			verifier.clone(),
			ParallelCompressionAdaptor::new(StdCompression::default()),
		)
		.expect("flat prover setup succeeds");

		(repeated, flat_constraint_system, value_vec, verifier, prover)
	}

	#[test]
	fn repeated_verifier_accepts_bound_flat_proof() {
		let (repeated, _, value_vec, verifier, prover) = setup_repeated_fixture(1 << 4, 1 << 5, 2);

		let mut prover_transcript = ProverTranscript::new(StdChallenger::default());
		prover
			.prove_repeated(&repeated, value_vec.clone(), &mut prover_transcript)
			.expect("repeated-bound flat prover succeeds");

		let mut verifier_transcript = prover_transcript.into_verifier();
		verifier
			.verify_repeated(value_vec.public(), &repeated, &mut verifier_transcript)
			.expect("repeated verifier accepts");
		verifier_transcript
			.finalize()
			.expect("repeated transcript is exhausted");
	}

	#[test]
	fn repeated_prover_uses_base_shift_keys() {
		let (repeated, flat_constraint_system, value_vec, verifier, flat_prover) =
			setup_repeated_fixture(1 << 4, 1 << 5, 2);
		let repeated_prover = Prover::<OptimalPackedB128, _, StdDigest>::setup_repeated(
			verifier.clone(),
			ParallelCompressionAdaptor::new(StdCompression::default()),
			repeated.clone(),
		)
		.expect("repeated prover setup succeeds");

		assert_eq!(
			flat_prover.key_collection().key_ranges.len(),
			flat_constraint_system.value_vec_layout.committed_total_len
		);
		assert_eq!(
			repeated_prover.key_collection().key_ranges.len(),
			repeated.base().value_vec_layout.committed_total_len
		);
		assert!(
			repeated_prover.key_collection().key_ranges.len()
				< flat_prover.key_collection().key_ranges.len()
		);

		let mut prover_transcript = ProverTranscript::new(StdChallenger::default());
		repeated_prover
			.prove_repeated(&repeated, value_vec.clone(), &mut prover_transcript)
			.expect("compact repeated prover succeeds");

		let mut verifier_transcript = prover_transcript.into_verifier();
		verifier
			.verify_repeated(value_vec.public(), &repeated, &mut verifier_transcript)
			.expect("repeated verifier accepts compact-key proof");
		verifier_transcript
			.finalize()
			.expect("repeated transcript is exhausted");
	}

	#[test]
	fn repeated_prover_accepts_shared_base_constants() {
		const LOG_INV_RATE: usize = 1;
		let base_constraint_system =
			make_base_constraint_system_with_shared_constant(1 << 4, 1 << 5);
		let repeated = RepeatedConstraintSystem::new(base_constraint_system, 2);
		let flat_constraint_system = repeated.to_flat_constraint_system();
		let base_layout = repeated.base().value_vec_layout.clone();
		let instances = (0..1usize << repeated.log_instances())
			.map(|instance| {
				ValueVec::new_from_data(
					base_layout.clone(),
					vec![Word::ALL_ONE, Word::from_u64(instance as u64)],
					vec![
						Word::from_u64(instance as u64);
						base_layout.committed_total_len - base_layout.offset_witness
					],
				)
				.expect("base instance value vec has matching layout")
			})
			.collect::<Vec<_>>();
		let value_vec = repeated
			.to_flat_value_vec(&instances)
			.expect("instance value vecs flatten");
		verify_constraints(&flat_constraint_system, &value_vec)
			.expect("witness satisfies repeated circuit with shared constants");

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
		.expect("repeated prover setup with shared constants succeeds");

		let mut prover_transcript = ProverTranscript::new(StdChallenger::default());
		prover
			.prove_repeated(&repeated, value_vec.clone(), &mut prover_transcript)
			.expect("compact repeated prover with shared constants succeeds");

		let mut verifier_transcript = prover_transcript.into_verifier();
		verifier
			.verify_repeated(value_vec.public(), &repeated, &mut verifier_transcript)
			.expect("repeated verifier accepts shared constants");
		verifier_transcript
			.finalize()
			.expect("repeated transcript is exhausted");
	}

	#[test]
	fn repeated_prover_rejects_same_shape_different_flat_circuit() {
		const LOG_INV_RATE: usize = 1;
		let (repeated, mut flat_constraint_system, _, _, _) =
			setup_repeated_fixture(1 << 4, 1 << 5, 2);
		let first_second_instance_row = repeated.base().and_constraints.len();
		flat_constraint_system.and_constraints[first_second_instance_row]
			.a
			.swap(0, 1);

		assert!(repeated.matches_flat_shape(&flat_constraint_system));
		assert!(!repeated.matches_flat_constraint_system(&flat_constraint_system));

		let verifier = Verifier::<StdDigest, _>::setup(
			flat_constraint_system,
			LOG_INV_RATE,
			StdCompression::default(),
		)
		.expect("same-shape verifier setup succeeds");

		assert!(
			Prover::<OptimalPackedB128, _, StdDigest>::setup_repeated(
				verifier,
				ParallelCompressionAdaptor::new(StdCompression::default()),
				repeated,
			)
			.is_err()
		);
	}

	#[test]
	fn repeated_verifier_rejects_wrong_log_instances() {
		let (repeated, _, value_vec, verifier, prover) = setup_repeated_fixture(1 << 4, 1 << 5, 2);
		let wrong_repeated =
			RepeatedConstraintSystem::new(repeated.base().clone(), repeated.log_instances() + 1);

		let mut prover_transcript = ProverTranscript::new(StdChallenger::default());
		prover
			.prove_repeated(&repeated, value_vec.clone(), &mut prover_transcript)
			.expect("repeated-bound flat prover succeeds");

		let mut verifier_transcript = prover_transcript.into_verifier();
		assert!(
			verifier
				.verify_repeated(value_vec.public(), &wrong_repeated, &mut verifier_transcript)
				.is_err()
		);
	}

	#[test]
	fn repeated_verifier_rejects_wrong_base_shape_binding() {
		let (repeated, _, value_vec, verifier, prover) = setup_repeated_fixture(1 << 4, 1 << 5, 2);
		let wrong_repeated = {
			let mut base = repeated.base().clone();
			base.and_constraints[0].a.swap(0, 1);
			RepeatedConstraintSystem::new(base, repeated.log_instances())
		};

		let mut prover_transcript = ProverTranscript::new(StdChallenger::default());
		prover
			.prove_repeated(&repeated, value_vec.clone(), &mut prover_transcript)
			.expect("repeated-bound flat prover succeeds");

		let mut verifier_transcript = prover_transcript.into_verifier();
		assert!(
			verifier
				.verify_repeated(value_vec.public(), &wrong_repeated, &mut verifier_transcript)
				.is_err()
		);
	}

	#[test]
	#[ignore = "prints end-to-end flat verifier vs repeated verifier runtimes"]
	fn repeated_e2e_verifier_print_runtimes() {
		const LOG_INV_RATE: usize = 1;
		let base_constraint_count = 1 << 8;
		let base_value_count = 1 << 9;
		let base_constraint_system =
			make_base_constraint_system(base_constraint_count, base_value_count);

		println!(
			"End-to-end verifier timing with flat prover, flat PCS, flat public IO, and structured repeated Shift monster check."
		);
		println!(
			"Base shape: {base_constraint_count} AND rows, {} MUL row, {base_value_count} values",
			base_constraint_system.mul_constraints.len()
		);
		println!(
			"log_instances,instances,flat_and_constraints,prove_ms,flat_verify_ms,repeated_verify_ms,speedup"
		);

		for log_instances in [0usize, 4, 8, 10] {
			let repeated =
				RepeatedConstraintSystem::new(base_constraint_system.clone(), log_instances);
			let flat_constraint_system = repeated.to_flat_constraint_system();
			let value_vec = zero_value_vec(&flat_constraint_system);
			verify_constraints(&flat_constraint_system, &value_vec)
				.expect("zero witness satisfies the repeated toy circuit");

			let verifier = Verifier::<StdDigest, _>::setup_repeated(
				&repeated,
				LOG_INV_RATE,
				StdCompression::default(),
			)
			.expect("repeated verifier setup succeeds");
			let prover = Prover::<OptimalPackedB128, _, StdDigest>::setup(
				verifier.clone(),
				ParallelCompressionAdaptor::new(StdCompression::default()),
			)
			.expect("flat prover setup succeeds");

			let (flat_prover_transcript, prove_elapsed) = elapsed_for(|| {
				let mut prover_transcript = ProverTranscript::new(StdChallenger::default());
				prover
					.prove(value_vec.clone(), &mut prover_transcript)
					.expect("flat prover succeeds");
				prover_transcript
			});
			let repeated_prover_transcript = {
				let mut prover_transcript = ProverTranscript::new(StdChallenger::default());
				prover
					.prove_repeated(&repeated, value_vec.clone(), &mut prover_transcript)
					.expect("repeated-bound flat prover succeeds");
				prover_transcript
			};

			let verify_iterations = if log_instances >= 8 { 8 } else { 32 };
			let flat_verify_elapsed = average_elapsed(verify_iterations, || {
				let mut verifier_transcript = flat_prover_transcript.clone().into_verifier();
				verifier
					.verify(value_vec.public(), &mut verifier_transcript)
					.expect("flat verifier accepts");
				verifier_transcript
					.finalize()
					.expect("flat transcript is exhausted");
			});

			let repeated_verify_elapsed = average_elapsed(verify_iterations, || {
				let mut verifier_transcript = repeated_prover_transcript.clone().into_verifier();
				verifier
					.verify_repeated(value_vec.public(), &repeated, &mut verifier_transcript)
					.expect("repeated verifier accepts");
				verifier_transcript
					.finalize()
					.expect("repeated transcript is exhausted");
			});

			let prove_ms = prove_elapsed.as_secs_f64() * 1_000.0;
			let flat_verify_ms = flat_verify_elapsed.as_secs_f64() * 1_000.0;
			let repeated_verify_ms = repeated_verify_elapsed.as_secs_f64() * 1_000.0;
			println!(
				"{log_instances},{},{},{prove_ms:.3},{flat_verify_ms:.3},{repeated_verify_ms:.3},{:.2}x",
				1usize << log_instances,
				flat_constraint_system.and_constraints.len(),
				flat_verify_ms / repeated_verify_ms,
			);
		}
	}

	#[test]
	#[ignore = "prints flat vs repeated prover setup/key/prove runtimes"]
	fn repeated_prover_key_materialization_print_runtimes() {
		const LOG_INV_RATE: usize = 1;
		let base_constraint_count = 1 << 8;
		let base_value_count = 1 << 9;
		let base_constraint_system =
			make_base_constraint_system(base_constraint_count, base_value_count);

		println!(
			"Flat prover materializes Shift keys for every instance; repeated prover materializes the base keys once and applies instance offsets during Shift proving."
		);
		println!(
			"Base shape: {base_constraint_count} AND rows, {} MUL row, {base_value_count} values",
			base_constraint_system.mul_constraints.len()
		);
		println!(
			"log_instances,instances,flat_key_words,repeated_key_words,flat_keys,repeated_keys,flat_setup_ms,repeated_setup_ms,flat_repeated_prove_ms,compact_repeated_prove_ms"
		);

		for log_instances in [0usize, 4, 8, 10] {
			let repeated =
				RepeatedConstraintSystem::new(base_constraint_system.clone(), log_instances);
			let flat_constraint_system = repeated.to_flat_constraint_system();
			let value_vec = zero_value_vec(&flat_constraint_system);
			verify_constraints(&flat_constraint_system, &value_vec)
				.expect("zero witness satisfies the repeated toy circuit");

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
				.expect("flat prover setup succeeds")
			});
			let (repeated_prover, repeated_setup_elapsed) = elapsed_for(|| {
				Prover::<OptimalPackedB128, _, StdDigest>::setup_repeated(
					verifier.clone(),
					ParallelCompressionAdaptor::new(StdCompression::default()),
					repeated.clone(),
				)
				.expect("compact repeated prover setup succeeds")
			});

			let (_, flat_repeated_prove_elapsed) = elapsed_for(|| {
				let mut prover_transcript = ProverTranscript::new(StdChallenger::default());
				flat_prover
					.prove_repeated(&repeated, value_vec.clone(), &mut prover_transcript)
					.expect("flat-key repeated prover succeeds");
				prover_transcript
			});
			let (compact_transcript, compact_repeated_prove_elapsed) = elapsed_for(|| {
				let mut prover_transcript = ProverTranscript::new(StdChallenger::default());
				repeated_prover
					.prove_repeated(&repeated, value_vec.clone(), &mut prover_transcript)
					.expect("compact repeated prover succeeds");
				prover_transcript
			});

			let mut verifier_transcript = compact_transcript.into_verifier();
			verifier
				.verify_repeated(value_vec.public(), &repeated, &mut verifier_transcript)
				.expect("repeated verifier accepts compact-key proof");
			verifier_transcript
				.finalize()
				.expect("compact repeated transcript is exhausted");

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
}
