// Copyright 2026 The Binius Developers

use std::{env, hint::black_box, thread, time::Duration};

use binius_core::{constraint_system::Operand, word::Word};
use binius_field::{AESTowerField8b, Field, PackedAESBinaryField16x8b, Random};
use binius_ip::sumcheck::RoundCoeffs;
use binius_keccak_prove::{
	bit_ntt::{NttLookup, upper_half_domains, upper_half_residual_evals},
	round_message::{
		PackedFoldedOuterColumns, folded_outer_claim, folded_outer_columns,
		pack_folded_outer_columns, par_first_round_claim_small_weights,
		par_upper_half_round_message, par_upper_half_round_message_small_weights,
		prove_spartan_outer_after_first_round_with_claim,
		prove_spartan_outer_after_first_round_with_claim_fused_adaptive,
		prove_spartan_outer_from_folded_columns_with_claim,
		prove_spartan_outer_from_folded_columns_with_claim_fused_adaptive,
		prove_spartan_outer_from_folded_columns_with_claim_packed,
		prove_spartan_outer_from_folded_columns_with_claim_packed_fused,
		prove_spartan_outer_from_packed_folded_columns_with_claim,
		prove_spartan_outer_from_packed_folded_columns_with_claim_fused,
		prove_spartan_outer_from_packed_folded_columns_with_claim_persistent_fused,
		upper_half_round_message, upper_half_round_message_small_weights,
	},
	trace::{PermutationTrace, RoundTrace, State},
	v0,
};
use binius_math::{
	BinarySubspace, FieldBuffer,
	multilinear::eq::{eq_ind_partial_eval, eq_ind_partial_eval_scalars},
	univariate::lagrange_evals_scalars,
};
use binius_prover::{
	OptimalPackedB128, Prover,
	and_reduction::{
		prover_setup::ntt_lookup_from_prover_message_domain,
		sumcheck_round_messages::univariate_round_message_extension_domain,
	},
	hash::parallel_compression::ParallelCompressionAdaptor,
};
use binius_transcript::ProverTranscript;
use binius_verifier::{
	Verifier,
	config::{B128, LOG_WORD_SIZE_BITS, StdChallenger},
	hash::{StdCompression, StdDigest},
	protocols::{
		bitand::SKIPPED_VARS,
		shift::{OperatorData as VerifierOperatorData, evaluate_monster_multilinear_for_operation},
	},
};
use criterion::{BenchmarkId, Criterion, Throughput, criterion_group, criterion_main};
use rand::{Rng, SeedableRng, rngs::StdRng};

const KECCAK_BENCH_PERMS: usize = 128;
const KECCAK_ROUNDS_PER_PERM: usize = 24;
const KECCAK_LANES_PER_ROUND: usize = binius_keccak_prove::constants::N_LANES;
const KECCAK_SCALE_CHUNK_PERMS: usize = 65_536;

#[derive(Clone)]
struct SpartanOuterScaleInstance {
	packed_columns: PackedFoldedOuterColumns<OptimalPackedB128>,
	first_round_challenge: B128,
	zerocheck_challenges: Vec<B128>,
	sumcheck_challenges: Vec<B128>,
	folded_claim: B128,
}

#[derive(Clone)]
struct SpartanOuterOneFsChunkInstance {
	packed_columns: PackedFoldedOuterColumns<OptimalPackedB128>,
	folded_claim: B128,
}

fn eval_word_direct(
	input_domain: &BinarySubspace<AESTowerField8b>,
	eval_points: &[AESTowerField8b],
	word: u64,
) -> Vec<AESTowerField8b> {
	eval_points
		.iter()
		.map(|&eval_point| {
			lagrange_evals_scalars(input_domain, eval_point)
				.into_iter()
				.enumerate()
				.filter_map(|(bit_idx, lagrange_eval)| {
					(((word >> bit_idx) & 1) == 1).then_some(lagrange_eval)
				})
				.sum()
		})
		.collect()
}

fn bench_ntt_lookup(c: &mut Criterion) {
	let (input_domain, output_domain) = upper_half_domains::<AESTowerField8b>();
	let lookup = NttLookup::<PackedAESBinaryField16x8b>::new(&input_domain, &output_domain);

	let mut rng = StdRng::seed_from_u64(10);
	let words: Vec<_> = (0..1024).map(|_| rng.random::<u64>()).collect();

	let mut group = c.benchmark_group("keccak_ntt_lookup");

	group.throughput(Throughput::Elements(1));
	group.bench_function("direct_lagrange_word", |bench| {
		let mut i = 0;
		bench.iter(|| {
			let word = words[i & (words.len() - 1)];
			i += 1;
			black_box(eval_word_direct(&input_domain, &output_domain, word))
		});
	});

	group.throughput(Throughput::Elements(1));
	group.bench_function("byte_lookup_word", |bench| {
		let mut i = 0;
		bench.iter(|| {
			let word = words[i & (words.len() - 1)];
			i += 1;
			black_box(lookup.eval_word(word))
		});
	});

	group.throughput(Throughput::Elements(1));
	group.bench_function("keccak_lookup_precompute", |bench| {
		bench.iter(|| {
			black_box(NttLookup::<PackedAESBinaryField16x8b>::new(&input_domain, &output_domain))
		});
	});

	group.throughput(Throughput::Elements(1));
	group.bench_function("production_bitand_lookup_precompute", |bench| {
		bench.iter(|| {
			let prover_message_domain = BinarySubspace::<AESTowerField8b>::with_dim(
				binius_keccak_prove::bit_ntt::LOG_LANE_BITS + 1,
			);
			black_box(ntt_lookup_from_prover_message_domain::<PackedAESBinaryField16x8b>(
				prover_message_domain,
			))
		});
	});
}

fn bench_keccak_residuals(c: &mut Criterion) {
	let (input_domain, output_domain) = upper_half_domains::<AESTowerField8b>();
	let lookup = NttLookup::<PackedAESBinaryField16x8b>::new(&input_domain, &output_domain);

	let mut rng = StdRng::seed_from_u64(11);
	let traces: Vec<_> = (0..KECCAK_BENCH_PERMS)
		.map(|_| PermutationTrace::new(rng.random::<State>()))
		.collect();
	let round_traces: Vec<_> = traces.iter().flat_map(|trace| trace.rounds).collect();
	let eq_weights: Vec<_> = (0..round_traces.len() * KECCAK_LANES_PER_ROUND)
		.map(|_| B128::from(rng.random::<AESTowerField8b>()))
		.collect();
	let small_eq_weights: Vec<_> = (0..round_traces.len() * KECCAK_LANES_PER_ROUND)
		.map(|_| rng.random::<AESTowerField8b>())
		.collect();
	let first_round_challenge = B128::random(&mut rng);

	let mut group = c.benchmark_group("keccak_residuals");

	for rounds_per_iter in [1, 24] {
		group.throughput(Throughput::Elements((rounds_per_iter * KECCAK_LANES_PER_ROUND) as u64));
		group.bench_function(
			BenchmarkId::new("upper_half_residual_evals", rounds_per_iter),
			|bench| {
				let mut trace_idx = 0;
				let mut round = 0;
				bench.iter(|| {
					let mut out = None;
					for _ in 0..rounds_per_iter {
						let round_trace = traces[trace_idx].rounds[round];
						out = Some(upper_half_residual_evals::<PackedAESBinaryField16x8b>(
							&lookup,
							&round_trace.pre_chi,
							&round_trace.output,
							round,
						));

						round += 1;
						if round == 24 {
							round = 0;
							trace_idx = (trace_idx + 1) & (traces.len() - 1);
						}
					}
					black_box(out)
				});
			},
		);
	}

	let accumulator_constraints =
		KECCAK_BENCH_PERMS * KECCAK_ROUNDS_PER_PERM * KECCAK_LANES_PER_ROUND;
	assert_eq!(accumulator_constraints, round_traces.len() * KECCAK_LANES_PER_ROUND);
	group.throughput(Throughput::Elements(accumulator_constraints as u64));
	group.bench_function("upper_half_round_message_seq/128_perms", |bench| {
		bench.iter(|| {
			black_box(upper_half_round_message::<B128, PackedAESBinaryField16x8b>(
				&lookup,
				&round_traces,
				&eq_weights,
			))
		});
	});

	group.throughput(Throughput::Elements(accumulator_constraints as u64));
	group.bench_function("upper_half_round_message_par/128_perms", |bench| {
		bench.iter(|| {
			black_box(par_upper_half_round_message::<B128, PackedAESBinaryField16x8b>(
				&lookup,
				&round_traces,
				&eq_weights,
			))
		});
	});

	group.throughput(Throughput::Elements(accumulator_constraints as u64));
	group.bench_function("upper_half_round_message_small_seq/128_perms", |bench| {
		bench.iter(|| {
			black_box(upper_half_round_message_small_weights::<B128, PackedAESBinaryField16x8b>(
				&lookup,
				&round_traces,
				&small_eq_weights,
			))
		});
	});

	group.throughput(Throughput::Elements(accumulator_constraints as u64));
	group.bench_function("upper_half_round_message_small_par/128_perms", |bench| {
		bench.iter(|| {
			black_box(
				par_upper_half_round_message_small_weights::<B128, PackedAESBinaryField16x8b>(
					&lookup,
					&round_traces,
					&small_eq_weights,
				),
			)
		});
	});

	group.throughput(Throughput::Elements(accumulator_constraints as u64));
	group.bench_function("first_round_claim_small_par/128_perms", |bench| {
		bench.iter(|| {
			black_box(par_first_round_claim_small_weights::<B128, PackedAESBinaryField16x8b>(
				&lookup,
				&round_traces,
				&small_eq_weights,
				first_round_challenge,
			))
		});
	});
}

fn bench_production_bitand_round_message(c: &mut Criterion) {
	for log_num_rows in [12, 22] {
		bench_production_bitand_round_message_size(c, log_num_rows);
	}
}

fn bench_keccak_v0_production_path(c: &mut Criterion) {
	let mut group = c.benchmark_group("keccak_v0_production_path");
	group.sample_size(10);
	group.measurement_time(Duration::from_secs(6));

	for n_permutations in v0_production_path_perm_counts() {
		let mut rng = StdRng::seed_from_u64(20 + n_permutations as u64);
		let traces: Vec<_> = (0..n_permutations)
			.map(|_| PermutationTrace::new(rng.random::<State>()))
			.collect();
		let witness = binius_keccak_prove::witness::CommittedKeccakWitness::from_traces(&traces);
		let constraint_system = v0::constraint_system(n_permutations);
		let value_vec = v0::value_vec(&witness);
		let verifier =
			Verifier::<StdDigest, _>::setup(constraint_system, 1, StdCompression::default())
				.unwrap();
		let prover = Prover::<OptimalPackedB128, _, StdDigest>::setup(
			verifier.clone(),
			ParallelCompressionAdaptor::new(StdCompression::default()),
		)
		.unwrap();

		let mut proof_transcript = ProverTranscript::new(StdChallenger::default());
		prover
			.prove(value_vec.clone(), &mut proof_transcript)
			.unwrap();
		let proof = proof_transcript.finalize();
		let constraint_rows =
			n_permutations * KECCAK_ROUNDS_PER_PERM * (KECCAK_LANES_PER_ROUND + 5);

		group.throughput(Throughput::Elements(constraint_rows as u64));
		group.bench_function(
			BenchmarkId::new("prove_full_production_path", n_permutations),
			|bench| {
				bench.iter(|| {
					let mut prover_transcript = ProverTranscript::new(StdChallenger::default());
					prover
						.prove(value_vec.clone(), &mut prover_transcript)
						.unwrap();
					black_box(prover_transcript)
				});
			},
		);

		group.throughput(Throughput::Elements(constraint_rows as u64));
		group.bench_function(
			BenchmarkId::new("verify_full_production_path", n_permutations),
			|bench| {
				bench.iter(|| {
					let mut verifier_transcript = binius_transcript::VerifierTranscript::new(
						StdChallenger::default(),
						proof.clone(),
					);
					verifier
						.verify(value_vec.public(), &mut verifier_transcript)
						.unwrap();
					black_box(verifier_transcript.finalize().unwrap())
				});
			},
		);
	}
}

fn v0_production_path_perm_counts() -> Vec<usize> {
	env::var("KECCAK_V0_PRODUCTION_PERMS")
		.ok()
		.map(|value| {
			value
				.split(',')
				.map(str::trim)
				.filter(|part| !part.is_empty())
				.map(|part| {
					part.parse::<usize>()
						.expect("KECCAK_V0_PRODUCTION_PERMS must be comma-separated usize values")
				})
				.collect()
		})
		.unwrap_or_else(|| vec![1, 16, 128])
}

fn bench_keccak_v0_structured_verifier(c: &mut Criterion) {
	let mut group = c.benchmark_group("keccak_v0_structured_verifier");
	group.sample_size(10);

	for n_permutations in [128, 1024] {
		let mut rng = StdRng::seed_from_u64(21 + n_permutations as u64);
		let constraint_system = v0::constraint_system(n_permutations);
		let subspace = BinarySubspace::<B128>::with_dim(LOG_WORD_SIZE_BITS);
		let operator_data = VerifierOperatorData::new(
			B128::random(&mut rng),
			(0..constraint_system.and_constraints.len().ilog2())
				.map(|_| B128::random(&mut rng))
				.collect(),
			std::array::from_fn(|_| B128::random(&mut rng)),
		);
		let lambda = B128::random(&mut rng);
		let r_j = (0..LOG_WORD_SIZE_BITS)
			.map(|_| B128::random(&mut rng))
			.collect::<Vec<_>>();
		let r_s = (0..LOG_WORD_SIZE_BITS)
			.map(|_| B128::random(&mut rng))
			.collect::<Vec<_>>();
		let r_y = (0..constraint_system
			.value_vec_layout
			.committed_total_len
			.ilog2())
			.map(|_| B128::random(&mut rng))
			.collect::<Vec<_>>();

		let mut a = Vec::<&Operand>::with_capacity(constraint_system.and_constraints.len());
		let mut b = Vec::<&Operand>::with_capacity(constraint_system.and_constraints.len());
		let mut c = Vec::<&Operand>::with_capacity(constraint_system.and_constraints.len());
		for constraint in &constraint_system.and_constraints {
			a.push(&constraint.a);
			b.push(&constraint.b);
			c.push(&constraint.c);
		}
		let operand_vecs = [a, b, c];
		let constraint_rows =
			n_permutations * KECCAK_ROUNDS_PER_PERM * (KECCAK_LANES_PER_ROUND + 5);
		group.throughput(Throughput::Elements(constraint_rows as u64));

		group.bench_function(
			BenchmarkId::new("generic_shift_monster_eval", n_permutations),
			|bench| {
				bench.iter(|| {
					black_box(
						evaluate_monster_multilinear_for_operation(
							&operand_vecs,
							&operator_data,
							&subspace,
							lambda,
							&r_j,
							&r_s,
							&r_y,
						)
						.unwrap(),
					)
				});
			},
		);

		group.bench_function(
			BenchmarkId::new("structured_shift_monster_eval", n_permutations),
			|bench| {
				bench.iter(|| {
					black_box(v0::structured_bitand_monster_eval(
						n_permutations,
						&operator_data,
						&subspace,
						lambda,
						&r_j,
						&r_s,
						&r_y,
					))
				});
			},
		);
	}
}

fn bench_keccak_first_round_claim_scale(c: &mut Criterion) {
	let (input_domain, output_domain) = upper_half_domains::<AESTowerField8b>();
	let lookup = NttLookup::<PackedAESBinaryField16x8b>::new(&input_domain, &output_domain);

	let mut group = c.benchmark_group("keccak_first_round_claim_scale");
	group.sample_size(10);
	group.measurement_time(Duration::from_secs(6));

	for total_perms in [2048, 4096, 8192, 16384, 32768, 65536, 131072, 196608] {
		let chunks = scale_chunks(total_perms);
		let mut rng = StdRng::seed_from_u64(17 + total_perms as u64);
		let first_round_challenge = B128::random(&mut rng);
		let total_constraints = total_perms * KECCAK_ROUNDS_PER_PERM * KECCAK_LANES_PER_ROUND;
		group.throughput(Throughput::Elements(total_constraints as u64));
		group.bench_function(
			BenchmarkId::new("first_round_claim_small_par_distinct", total_perms),
			|bench| {
				bench.iter(|| {
					let mut claim = B128::ZERO;
					for (round_traces, small_eq_weights) in &chunks {
						claim += par_first_round_claim_small_weights::<
							B128,
							PackedAESBinaryField16x8b,
						>(
							&lookup, round_traces, small_eq_weights, first_round_challenge
						);
					}
					black_box(claim)
				});
			},
		);
	}
}

fn bench_keccak_spartan_outer(c: &mut Criterion) {
	let mut rng = StdRng::seed_from_u64(18);
	let traces: Vec<_> = (0..KECCAK_BENCH_PERMS)
		.map(|_| PermutationTrace::new(rng.random::<State>()))
		.collect();
	let round_traces: Vec<_> = traces.iter().flat_map(|trace| trace.rounds).collect();
	let first_round_challenge = B128::random(&mut rng);
	let columns = folded_outer_columns::<B128, PackedAESBinaryField16x8b>(
		&round_traces,
		first_round_challenge,
	);
	let zerocheck_challenges: Vec<_> = (0..columns.log_rows)
		.map(|_| B128::random(&mut rng))
		.collect();
	let folded_claim = folded_outer_claim(&columns, &zerocheck_challenges);
	let sumcheck_challenges: Vec<_> = (0..columns.log_rows)
		.map(|_| B128::random(&mut rng))
		.collect();
	let total_constraints = KECCAK_BENCH_PERMS * KECCAK_ROUNDS_PER_PERM * KECCAK_LANES_PER_ROUND;

	let mut group = c.benchmark_group("keccak_spartan_outer");
	group.sample_size(10);

	group.throughput(Throughput::Elements(total_constraints as u64));
	group.bench_function("folded_outer_columns/128_perms", |bench| {
		bench.iter(|| {
			black_box(folded_outer_columns::<B128, PackedAESBinaryField16x8b>(
				&round_traces,
				first_round_challenge,
			))
		});
	});

	group.throughput(Throughput::Elements(total_constraints as u64));
	group.bench_function("folded_outer_claim/128_perms", |bench| {
		bench.iter(|| black_box(folded_outer_claim(&columns, &zerocheck_challenges)));
	});

	group.throughput(Throughput::Elements(total_constraints as u64));
	group.bench_function("prove_after_univariate_skip/128_perms", |bench| {
		bench.iter(|| {
			black_box(
				prove_spartan_outer_after_first_round_with_claim::<
					B128,
					PackedAESBinaryField16x8b,
				>(
					&round_traces,
					first_round_challenge,
					zerocheck_challenges.clone(),
					&sumcheck_challenges,
					folded_claim,
				)
				.unwrap(),
			)
		});
	});

	group.throughput(Throughput::Elements(total_constraints as u64));
	group.bench_function("prove_after_univariate_skip_fused_adaptive/128_perms", |bench| {
		bench.iter(|| {
			black_box(prove_spartan_outer_after_first_round_with_claim_fused_adaptive::<
				B128,
				PackedAESBinaryField16x8b,
			>(
				&round_traces,
				first_round_challenge,
				zerocheck_challenges.clone(),
				&sumcheck_challenges,
				folded_claim,
			))
		});
	});
}

fn bench_keccak_spartan_outer_scale(c: &mut Criterion) {
	let mut group = c.benchmark_group("keccak_spartan_outer_scale");
	group.sample_size(10);
	group.measurement_time(Duration::from_secs(6));

	for total_perms in spartan_outer_scale_perm_counts() {
		let mut rng = StdRng::seed_from_u64(19 + total_perms as u64);
		let mut round_traces = Vec::with_capacity(total_perms * KECCAK_ROUNDS_PER_PERM);
		for _ in 0..total_perms {
			round_traces.extend(PermutationTrace::new(rng.random::<State>()).rounds);
		}
		let first_round_challenge = B128::random(&mut rng);
		let log_rows = (round_traces.len() * KECCAK_LANES_PER_ROUND)
			.next_power_of_two()
			.ilog2() as usize;
		let zerocheck_challenges: Vec<_> = (0..log_rows).map(|_| B128::random(&mut rng)).collect();
		let columns = folded_outer_columns::<B128, PackedAESBinaryField16x8b>(
			&round_traces,
			first_round_challenge,
		);
		let packed_columns = pack_folded_outer_columns::<B128, OptimalPackedB128>(columns.clone());
		let folded_claim = folded_outer_claim(&columns, &zerocheck_challenges);
		let sumcheck_challenges: Vec<_> = (0..log_rows).map(|_| B128::random(&mut rng)).collect();
		let total_constraints = total_perms * KECCAK_ROUNDS_PER_PERM * KECCAK_LANES_PER_ROUND;

		group.throughput(Throughput::Elements(total_constraints as u64));
		group.bench_function(
			BenchmarkId::new("prove_after_univariate_skip_distinct", total_perms),
			|bench| {
				bench.iter(|| {
					black_box(
						prove_spartan_outer_after_first_round_with_claim::<
							B128,
							PackedAESBinaryField16x8b,
						>(
							&round_traces,
							first_round_challenge,
							zerocheck_challenges.clone(),
							&sumcheck_challenges,
							folded_claim,
						)
						.unwrap(),
					)
				});
			},
		);

		group.throughput(Throughput::Elements(total_constraints as u64));
		group.bench_function(
			BenchmarkId::new("prove_after_univariate_skip_fused_adaptive_distinct", total_perms),
			|bench| {
				bench.iter(|| {
					black_box(prove_spartan_outer_after_first_round_with_claim_fused_adaptive::<
						B128,
						PackedAESBinaryField16x8b,
					>(
						&round_traces,
						first_round_challenge,
						zerocheck_challenges.clone(),
						&sumcheck_challenges,
						folded_claim,
					))
				});
			},
		);

		group.throughput(Throughput::Elements(total_constraints as u64));
		group.bench_function(
			BenchmarkId::new("prove_from_folded_columns_distinct", total_perms),
			|bench| {
				bench.iter(|| {
					black_box(
						prove_spartan_outer_from_folded_columns_with_claim(
							columns.clone(),
							first_round_challenge,
							zerocheck_challenges.clone(),
							&sumcheck_challenges,
							folded_claim,
						)
						.unwrap(),
					)
				});
			},
		);

		group.throughput(Throughput::Elements(total_constraints as u64));
		group.bench_function(
			BenchmarkId::new("prove_from_folded_columns_fused_adaptive_distinct", total_perms),
			|bench| {
				bench.iter(|| {
					black_box(prove_spartan_outer_from_folded_columns_with_claim_fused_adaptive(
						columns.clone(),
						first_round_challenge,
						zerocheck_challenges.clone(),
						&sumcheck_challenges,
						folded_claim,
					))
				});
			},
		);

		group.throughput(Throughput::Elements(total_constraints as u64));
		group.bench_function(
			BenchmarkId::new("prove_from_folded_columns_packed_distinct", total_perms),
			|bench| {
				bench.iter(|| {
					black_box(
						prove_spartan_outer_from_folded_columns_with_claim_packed::<
							B128,
							OptimalPackedB128,
						>(
							columns.clone(),
							first_round_challenge,
							zerocheck_challenges.clone(),
							&sumcheck_challenges,
							folded_claim,
						)
						.unwrap(),
					)
				});
			},
		);

		group.throughput(Throughput::Elements(total_constraints as u64));
		group.bench_function(
			BenchmarkId::new("prove_from_folded_columns_packed_fused_distinct", total_perms),
			|bench| {
				bench.iter(|| {
					black_box(
						prove_spartan_outer_from_folded_columns_with_claim_packed_fused::<
							B128,
							OptimalPackedB128,
						>(
							columns.clone(),
							first_round_challenge,
							zerocheck_challenges.clone(),
							&sumcheck_challenges,
							folded_claim,
						)
						.unwrap(),
					)
				});
			},
		);

		group.throughput(Throughput::Elements(total_constraints as u64));
		group.bench_function(
			BenchmarkId::new("prove_from_packed_folded_columns_distinct", total_perms),
			|bench| {
				bench.iter(|| {
					black_box(
						prove_spartan_outer_from_packed_folded_columns_with_claim::<
							B128,
							OptimalPackedB128,
						>(
							packed_columns.clone(),
							first_round_challenge,
							zerocheck_challenges.clone(),
							&sumcheck_challenges,
							folded_claim,
						)
						.unwrap(),
					)
				});
			},
		);

		group.throughput(Throughput::Elements(total_constraints as u64));
		group.bench_function(
			BenchmarkId::new("prove_from_packed_folded_columns_fused_distinct", total_perms),
			|bench| {
				bench.iter(|| {
					black_box(
						prove_spartan_outer_from_packed_folded_columns_with_claim_fused::<
							B128,
							OptimalPackedB128,
						>(
							packed_columns.clone(),
							first_round_challenge,
							zerocheck_challenges.clone(),
							&sumcheck_challenges,
							folded_claim,
						)
						.unwrap(),
					)
				});
			},
		);

		group.throughput(Throughput::Elements(total_constraints as u64));
		group.bench_function(
			BenchmarkId::new(
				"prove_from_packed_folded_columns_persistent_fused_distinct",
				total_perms,
			),
			|bench| {
				bench.iter(|| {
					black_box(
						prove_spartan_outer_from_packed_folded_columns_with_claim_persistent_fused::<
							B128,
							OptimalPackedB128,
						>(
							packed_columns.clone(),
							first_round_challenge,
							zerocheck_challenges.clone(),
							&sumcheck_challenges,
							folded_claim,
						)
						.unwrap(),
					)
				});
			},
		);

		if let Some(chunk_perms) = spartan_outer_coarse_chunk_perms()
			&& total_perms > chunk_perms
		{
			let instances = spartan_outer_coarse_instances(total_perms, chunk_perms);
			let jobs = spartan_outer_coarse_jobs().min(instances.len()).max(1);
			group.throughput(Throughput::Elements(total_constraints as u64));
			group.bench_function(
				BenchmarkId::new(
					format!(
						"prove_packed_persistent_fused_coarse_{}_jobs_{}_per_chunk",
						jobs, chunk_perms
					),
					total_perms,
				),
				|bench| {
					bench.iter(|| black_box(prove_coarse_persistent_batch(&instances, jobs)));
				},
			);
		}

		if let Some(chunk_perms) = spartan_outer_one_fs_chunk_perms()
			&& total_perms > chunk_perms
		{
			let one_fs =
				spartan_outer_one_fs_chunks(total_perms, chunk_perms, first_round_challenge);
			if spartan_outer_one_fs_verify() {
				verify_one_fs_chunked_matches_permuted_generic(&one_fs);
			}
			let jobs = spartan_outer_one_fs_jobs()
				.min(one_fs.instances.len())
				.max(1);
			group.throughput(Throughput::Elements(total_constraints as u64));
			group.bench_function(
				BenchmarkId::new(
					format!(
						"prove_packed_persistent_fused_one_fs_{}_jobs_{}_per_chunk",
						jobs, chunk_perms
					),
					total_perms,
				),
				|bench| {
					bench.iter(|| black_box(prove_one_fs_chunked_global_batch(&one_fs, jobs)));
				},
			);
		}
	}
}

struct SpartanOuterOneFsChunks {
	instances: Vec<SpartanOuterOneFsChunkInstance>,
	first_round_challenge: B128,
	chunk_zerocheck_challenges: Vec<B128>,
	local_zerocheck_challenges: Vec<B128>,
	local_sumcheck_challenges: Vec<B128>,
	chunk_sumcheck_challenges: Vec<B128>,
	folded_claim: B128,
}

fn spartan_outer_one_fs_chunk_perms() -> Option<usize> {
	env::var("KECCAK_SPARTAN_OUTER_ONE_FS_CHUNK_PERMS")
		.ok()
		.map(|value| {
			value
				.parse::<usize>()
				.expect("KECCAK_SPARTAN_OUTER_ONE_FS_CHUNK_PERMS must be a usize")
		})
}

fn spartan_outer_one_fs_verify() -> bool {
	env::var("KECCAK_SPARTAN_OUTER_ONE_FS_VERIFY")
		.ok()
		.is_some_and(|value| value == "1" || value.eq_ignore_ascii_case("true"))
}

fn spartan_outer_one_fs_jobs() -> usize {
	env::var("KECCAK_SPARTAN_OUTER_ONE_FS_JOBS")
		.ok()
		.and_then(|value| value.parse().ok())
		.or_else(|| thread::available_parallelism().ok().map(usize::from))
		.unwrap_or(1)
}

fn spartan_outer_one_fs_chunks(
	total_perms: usize,
	chunk_perms: usize,
	first_round_challenge: B128,
) -> SpartanOuterOneFsChunks {
	assert_eq!(total_perms % chunk_perms, 0);
	let chunk_count = total_perms / chunk_perms;
	assert!(chunk_count.is_power_of_two());
	let log_chunks = chunk_count.ilog2() as usize;

	let mut rng = StdRng::seed_from_u64(43 + total_perms as u64 * 17 + chunk_perms as u64 * 31);
	let mut instances = Vec::with_capacity(chunk_count);
	let mut local_log_rows = None;
	for chunk_idx in 0..chunk_count {
		let mut round_traces = Vec::with_capacity(chunk_perms * KECCAK_ROUNDS_PER_PERM);
		for _ in 0..chunk_perms {
			round_traces.extend(PermutationTrace::new(rng.random::<State>()).rounds);
		}
		let columns = folded_outer_columns::<B128, PackedAESBinaryField16x8b>(
			&round_traces,
			first_round_challenge,
		);
		local_log_rows = Some(*local_log_rows.get_or_insert(columns.log_rows));
		assert_eq!(columns.log_rows, local_log_rows.expect("local log rows set"));
		let packed_columns = pack_folded_outer_columns::<B128, OptimalPackedB128>(columns.clone());
		instances.push((chunk_idx, columns, packed_columns));
	}

	let local_log_rows = local_log_rows.expect("at least one chunk");
	let chunk_zerocheck_challenges: Vec<_> =
		(0..log_chunks).map(|_| B128::random(&mut rng)).collect();
	let local_zerocheck_challenges: Vec<_> = (0..local_log_rows)
		.map(|_| B128::random(&mut rng))
		.collect();
	let local_sumcheck_challenges: Vec<_> = (0..local_log_rows)
		.map(|_| B128::random(&mut rng))
		.collect();
	let chunk_sumcheck_challenges: Vec<_> =
		(0..log_chunks).map(|_| B128::random(&mut rng)).collect();
	let chunk_eq_weights = eq_ind_partial_eval_scalars(&chunk_zerocheck_challenges);

	let mut folded_claim = B128::ZERO;
	let instances = instances
		.into_iter()
		.map(|(chunk_idx, columns, packed_columns)| {
			let chunk_claim = folded_outer_claim(&columns, &local_zerocheck_challenges);
			folded_claim += chunk_eq_weights[chunk_idx] * chunk_claim;
			SpartanOuterOneFsChunkInstance {
				packed_columns,
				folded_claim: chunk_claim,
			}
		})
		.collect();

	SpartanOuterOneFsChunks {
		instances,
		first_round_challenge,
		chunk_zerocheck_challenges,
		local_zerocheck_challenges,
		local_sumcheck_challenges,
		chunk_sumcheck_challenges,
		folded_claim,
	}
}

fn prove_one_fs_chunked_global_batch(one_fs: &SpartanOuterOneFsChunks, jobs: usize) -> B128 {
	let chunk_eq_weights = eq_ind_partial_eval_scalars(&one_fs.chunk_zerocheck_challenges);
	let local_log_rows = one_fs.local_zerocheck_challenges.len();

	let mut worker_outputs = thread::scope(|scope| {
		let mut handles = Vec::new();
		for job_idx in 1..jobs {
			let chunk_eq_weights = &chunk_eq_weights;
			handles.push(scope.spawn(move || {
				prove_one_fs_chunked_global_worker(one_fs, chunk_eq_weights, jobs, job_idx)
			}));
		}

		let mut outputs = vec![prove_one_fs_chunked_global_worker(
			one_fs,
			&chunk_eq_weights,
			jobs,
			0,
		)];
		for handle in handles {
			outputs.push(handle.join().expect("one-FS chunk worker should not panic"));
		}
		outputs
	});

	let mut local_round_messages = vec![RoundCoeffs(vec![B128::ZERO; 3]); local_log_rows];
	let mut chunk_p = vec![B128::ZERO; one_fs.instances.len()];
	let mut chunk_q = vec![B128::ZERO; one_fs.instances.len()];
	let mut chunk_c = vec![B128::ZERO; one_fs.instances.len()];
	let mut tail_claim = B128::ZERO;
	for output in worker_outputs.drain(..) {
		for (acc, coeffs) in local_round_messages
			.iter_mut()
			.zip(output.local_round_messages)
		{
			*acc += &coeffs;
		}
		for (chunk_idx, evals, final_eval) in output.chunk_evals {
			chunk_p[chunk_idx] = evals[0];
			chunk_q[chunk_idx] = evals[1];
			chunk_c[chunk_idx] = evals[2];
			tail_claim += chunk_eq_weights[chunk_idx] * final_eval;
		}
	}

	let chunk_columns = PackedFoldedOuterColumns {
		p: FieldBuffer::<OptimalPackedB128>::from_values(&chunk_p),
		q: FieldBuffer::<OptimalPackedB128>::from_values(&chunk_q),
		c: FieldBuffer::<OptimalPackedB128>::from_values(&chunk_c),
		log_rows: one_fs.chunk_zerocheck_challenges.len(),
		row_count: one_fs.instances.len(),
	};
	let tail =
		prove_spartan_outer_from_packed_folded_columns_with_claim::<B128, OptimalPackedB128>(
			chunk_columns,
			one_fs.first_round_challenge,
			one_fs.chunk_zerocheck_challenges.clone(),
			&one_fs.chunk_sumcheck_challenges,
			tail_claim,
		)
		.unwrap();
	black_box(one_fs.folded_claim);
	tail.final_eval
}

fn verify_one_fs_chunked_matches_permuted_generic(one_fs: &SpartanOuterOneFsChunks) {
	let log_chunks = one_fs.chunk_zerocheck_challenges.len();
	let log_local = one_fs.local_zerocheck_challenges.len();
	let chunk_count = one_fs.instances.len();
	let local_rows = 1usize << log_local;
	let total_rows = local_rows * chunk_count;
	let mut p = vec![B128::ZERO; total_rows];
	let mut q = vec![B128::ZERO; total_rows];
	let mut c = vec![B128::ZERO; total_rows];

	for (chunk_idx, instance) in one_fs.instances.iter().enumerate() {
		for local_idx in 0..local_rows {
			let global_idx = chunk_idx + (local_idx << log_chunks);
			p[global_idx] = instance.packed_columns.p.get(local_idx);
			q[global_idx] = instance.packed_columns.q.get(local_idx);
			c[global_idx] = instance.packed_columns.c.get(local_idx);
		}
	}

	let packed_columns = PackedFoldedOuterColumns {
		p: FieldBuffer::<OptimalPackedB128>::from_values(&p),
		q: FieldBuffer::<OptimalPackedB128>::from_values(&q),
		c: FieldBuffer::<OptimalPackedB128>::from_values(&c),
		log_rows: log_local + log_chunks,
		row_count: total_rows,
	};
	let zerocheck_challenges = one_fs
		.chunk_zerocheck_challenges
		.iter()
		.chain(&one_fs.local_zerocheck_challenges)
		.copied()
		.collect::<Vec<_>>();
	let sumcheck_challenges = one_fs
		.local_sumcheck_challenges
		.iter()
		.chain(&one_fs.chunk_sumcheck_challenges)
		.copied()
		.collect::<Vec<_>>();

	let generic =
		prove_spartan_outer_from_packed_folded_columns_with_claim::<B128, OptimalPackedB128>(
			packed_columns,
			one_fs.first_round_challenge,
			zerocheck_challenges,
			&sumcheck_challenges,
			one_fs.folded_claim,
		)
		.unwrap();
	let one_fs_eval = prove_one_fs_chunked_global_batch(one_fs, spartan_outer_one_fs_jobs());
	assert_eq!(one_fs_eval, generic.final_eval);
}

struct OneFsWorkerOutput {
	local_round_messages: Vec<RoundCoeffs<B128>>,
	chunk_evals: Vec<(usize, [B128; 3], B128)>,
}

fn prove_one_fs_chunked_global_worker(
	one_fs: &SpartanOuterOneFsChunks,
	chunk_eq_weights: &[B128],
	jobs: usize,
	job_idx: usize,
) -> OneFsWorkerOutput {
	let mut local_round_messages =
		vec![RoundCoeffs(vec![B128::ZERO; 3]); one_fs.local_zerocheck_challenges.len()];
	let mut chunk_evals = Vec::new();
	for chunk_idx in (job_idx..one_fs.instances.len()).step_by(jobs) {
		let instance = &one_fs.instances[chunk_idx];
		let proof = prove_spartan_outer_from_packed_folded_columns_with_claim_persistent_fused::<
			B128,
			OptimalPackedB128,
		>(
			instance.packed_columns.clone(),
			one_fs.first_round_challenge,
			one_fs.local_zerocheck_challenges.clone(),
			&one_fs.local_sumcheck_challenges,
			instance.folded_claim,
		)
		.unwrap();
		let weight = chunk_eq_weights[chunk_idx];
		for (acc, coeffs) in local_round_messages.iter_mut().zip(proof.round_messages) {
			*acc += &(coeffs * weight);
		}
		chunk_evals.push((chunk_idx, proof.multilinear_evals, proof.final_eval));
	}

	OneFsWorkerOutput {
		local_round_messages,
		chunk_evals,
	}
}

fn spartan_outer_coarse_chunk_perms() -> Option<usize> {
	env::var("KECCAK_SPARTAN_OUTER_COARSE_CHUNK_PERMS")
		.ok()
		.map(|value| {
			value
				.parse::<usize>()
				.expect("KECCAK_SPARTAN_OUTER_COARSE_CHUNK_PERMS must be a usize")
		})
}

fn spartan_outer_coarse_jobs() -> usize {
	env::var("KECCAK_SPARTAN_OUTER_COARSE_JOBS")
		.ok()
		.and_then(|value| value.parse().ok())
		.or_else(|| thread::available_parallelism().ok().map(usize::from))
		.unwrap_or(1)
}

fn spartan_outer_coarse_instances(
	total_perms: usize,
	chunk_perms: usize,
) -> Vec<SpartanOuterScaleInstance> {
	let mut remaining_perms = total_perms;
	let mut chunk_idx = 0;
	let mut instances = Vec::new();

	while remaining_perms > 0 {
		let current_chunk_perms = remaining_perms.min(chunk_perms);
		let mut rng = StdRng::seed_from_u64(
			29 + total_perms as u64 * 17 + chunk_perms as u64 * 31 + chunk_idx as u64,
		);
		let mut round_traces = Vec::with_capacity(current_chunk_perms * KECCAK_ROUNDS_PER_PERM);
		for _ in 0..current_chunk_perms {
			round_traces.extend(PermutationTrace::new(rng.random::<State>()).rounds);
		}
		let first_round_challenge = B128::random(&mut rng);
		let columns = folded_outer_columns::<B128, PackedAESBinaryField16x8b>(
			&round_traces,
			first_round_challenge,
		);
		let zerocheck_challenges: Vec<_> = (0..columns.log_rows)
			.map(|_| B128::random(&mut rng))
			.collect();
		let folded_claim = folded_outer_claim(&columns, &zerocheck_challenges);
		let sumcheck_challenges: Vec<_> = (0..columns.log_rows)
			.map(|_| B128::random(&mut rng))
			.collect();
		let packed_columns = pack_folded_outer_columns::<B128, OptimalPackedB128>(columns);
		instances.push(SpartanOuterScaleInstance {
			packed_columns,
			first_round_challenge,
			zerocheck_challenges,
			sumcheck_challenges,
			folded_claim,
		});

		remaining_perms -= current_chunk_perms;
		chunk_idx += 1;
	}

	instances
}

fn prove_coarse_persistent_batch(instances: &[SpartanOuterScaleInstance], jobs: usize) -> B128 {
	thread::scope(|scope| {
		let mut handles = Vec::new();
		for job_idx in 1..jobs {
			handles.push(
				scope.spawn(move || prove_coarse_persistent_worker(instances, jobs, job_idx)),
			);
		}

		let mut acc = prove_coarse_persistent_worker(instances, jobs, 0);
		for handle in handles {
			acc += handle.join().expect("coarse worker should not panic");
		}
		acc
	})
}

fn prove_coarse_persistent_worker(
	instances: &[SpartanOuterScaleInstance],
	jobs: usize,
	job_idx: usize,
) -> B128 {
	let mut acc = B128::ZERO;
	for instance_idx in (job_idx..instances.len()).step_by(jobs) {
		let instance = &instances[instance_idx];
		let proof = prove_spartan_outer_from_packed_folded_columns_with_claim_persistent_fused::<
			B128,
			OptimalPackedB128,
		>(
			instance.packed_columns.clone(),
			instance.first_round_challenge,
			instance.zerocheck_challenges.clone(),
			&instance.sumcheck_challenges,
			instance.folded_claim,
		)
		.unwrap();
		acc += proof.final_eval;
	}
	acc
}

fn spartan_outer_scale_perm_counts() -> Vec<usize> {
	env::var("KECCAK_SPARTAN_OUTER_SCALE_PERMS")
		.ok()
		.map(|value| {
			value
				.split(',')
				.map(str::trim)
				.filter(|part| !part.is_empty())
				.map(|part| {
					part.parse::<usize>().expect(
						"KECCAK_SPARTAN_OUTER_SCALE_PERMS must be comma-separated usize values",
					)
				})
				.collect()
		})
		.unwrap_or_else(|| {
			vec![
				128, 256, 512, 1024, 2048, 4096, 8192, 12288, 16384, 24576, 28672, 32768, 55924,
				65536,
			]
		})
}

fn scale_chunks(total_perms: usize) -> Vec<(Vec<RoundTrace>, Vec<AESTowerField8b>)> {
	let mut remaining_perms = total_perms;
	let mut chunk_idx = 0;
	let mut chunks = Vec::new();

	while remaining_perms > 0 {
		let chunk_perms = remaining_perms.min(KECCAK_SCALE_CHUNK_PERMS);
		chunks.push(scale_chunk(total_perms, chunk_idx, chunk_perms));
		remaining_perms -= chunk_perms;
		chunk_idx += 1;
	}

	chunks
}

fn scale_chunk(
	total_perms: usize,
	chunk_idx: usize,
	chunk_perms: usize,
) -> (Vec<RoundTrace>, Vec<AESTowerField8b>) {
	let mut rng = StdRng::seed_from_u64(13 + total_perms as u64 * 17 + chunk_idx as u64);
	let mut round_traces = Vec::with_capacity(chunk_perms * KECCAK_ROUNDS_PER_PERM);
	for _ in 0..chunk_perms {
		round_traces.extend(PermutationTrace::new(rng.random::<State>()).rounds);
	}
	let small_eq_weights = (0..round_traces.len() * KECCAK_LANES_PER_ROUND)
		.map(|_| rng.random::<AESTowerField8b>())
		.collect();

	(round_traces, small_eq_weights)
}

fn bench_production_bitand_round_message_size(c: &mut Criterion, log_num_rows: usize) {
	let log_num_words = log_num_rows - SKIPPED_VARS;
	let mut rng = StdRng::seed_from_u64(12);
	let small_field_zerocheck_challenges = [
		AESTowerField8b::new(2),
		AESTowerField8b::new(4),
		AESTowerField8b::new(16),
	];
	let big_field_zerocheck_challenges =
		vec![
			B128::random(&mut rng);
			log_num_rows - SKIPPED_VARS - small_field_zerocheck_challenges.len()
		];
	let eq_ind_big_field_challenges = eq_ind_partial_eval(&big_field_zerocheck_challenges);

	let first_col: Vec<_> = (0..1 << log_num_words)
		.map(|_| Word(rng.random()))
		.collect();
	let second_col: Vec<_> = (0..1 << log_num_words)
		.map(|_| Word(rng.random()))
		.collect();
	let third_col: Vec<_> = first_col
		.iter()
		.zip(&second_col)
		.map(|(&a, &b)| a & b)
		.collect();

	let prover_message_domain = BinarySubspace::with_dim(SKIPPED_VARS + 1);
	let ntt_lookup =
		ntt_lookup_from_prover_message_domain::<PackedAESBinaryField16x8b>(prover_message_domain);

	let mut group =
		c.benchmark_group(format!("production_bitand_reference/log_rows={log_num_rows}"));
	if log_num_rows >= 22 {
		group.sample_size(10);
	}
	group.throughput(Throughput::Elements(1 << log_num_words));
	group.bench_function(
		BenchmarkId::new("univariate_round_message_extension_domain", log_num_rows),
		|bench| {
			bench.iter(|| {
				black_box(univariate_round_message_extension_domain::<
					B128,
					PackedAESBinaryField16x8b,
				>(
					&first_col,
					&second_col,
					&third_col,
					&eq_ind_big_field_challenges,
					&ntt_lookup,
					&small_field_zerocheck_challenges,
				))
			});
		},
	);
}

criterion_group!(
	keccak_ntt,
	bench_ntt_lookup,
	bench_keccak_residuals,
	bench_keccak_first_round_claim_scale,
	bench_keccak_spartan_outer,
	bench_keccak_spartan_outer_scale,
	bench_keccak_v0_production_path,
	bench_keccak_v0_structured_verifier,
	bench_production_bitand_round_message
);
criterion_main!(keccak_ntt);
