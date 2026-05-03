// Copyright 2026 The Binius Developers

use std::hint::black_box;

use binius_core::word::Word;
use binius_field::{AESTowerField8b, PackedAESBinaryField16x8b, Random};
use binius_keccak_prove::{
	bit_ntt::{NttLookup, upper_half_domains, upper_half_residual_evals},
	round_message::{par_upper_half_round_message, upper_half_round_message},
	trace::{PermutationTrace, State},
};
use binius_math::{
	BinarySubspace, multilinear::eq::eq_ind_partial_eval, univariate::lagrange_evals_scalars,
};
use binius_prover::and_reduction::{
	prover_setup::ntt_lookup_from_prover_message_domain,
	sumcheck_round_messages::univariate_round_message_extension_domain,
};
use binius_verifier::{config::B128, protocols::bitand::SKIPPED_VARS};
use criterion::{BenchmarkId, Criterion, Throughput, criterion_group, criterion_main};
use rand::{Rng, SeedableRng, rngs::StdRng};

const KECCAK_BENCH_PERMS: usize = 128;
const KECCAK_ROUNDS_PER_PERM: usize = 24;
const KECCAK_LANES_PER_ROUND: usize = binius_keccak_prove::constants::N_LANES;

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
}

fn bench_production_bitand_round_message(c: &mut Criterion) {
	for log_num_rows in [12, 22] {
		bench_production_bitand_round_message_size(c, log_num_rows);
	}
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
	bench_production_bitand_round_message
);
criterion_main!(keccak_ntt);
