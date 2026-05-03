// Copyright 2026 The Binius Developers

use std::hint::black_box;

use binius_field::{AESTowerField8b, PackedAESBinaryField16x8b};
use binius_keccak_prove::{
	bit_ntt::{NttLookup, upper_half_domains, upper_half_residual_evals},
	trace::{PermutationTrace, State},
};
use binius_math::{BinarySubspace, univariate::lagrange_evals_scalars};
use binius_prover::and_reduction::prover_setup::ntt_lookup_from_prover_message_domain;
use criterion::{BenchmarkId, Criterion, Throughput, criterion_group, criterion_main};
use rand::{Rng, SeedableRng, rngs::StdRng};

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
	group.throughput(Throughput::Elements(words.len() as u64));

	group.bench_function("direct_lagrange_word", |bench| {
		let mut i = 0;
		bench.iter(|| {
			let word = words[i & (words.len() - 1)];
			i += 1;
			black_box(eval_word_direct(&input_domain, &output_domain, word))
		});
	});

	group.bench_function("byte_lookup_word", |bench| {
		let mut i = 0;
		bench.iter(|| {
			let word = words[i & (words.len() - 1)];
			i += 1;
			black_box(lookup.eval_word(word))
		});
	});

	group.bench_function("keccak_lookup_precompute", |bench| {
		bench.iter(|| {
			black_box(NttLookup::<PackedAESBinaryField16x8b>::new(&input_domain, &output_domain))
		});
	});

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
	let traces: Vec<_> = (0..128)
		.map(|_| PermutationTrace::new(rng.random::<State>()))
		.collect();

	let mut group = c.benchmark_group("keccak_residuals");
	group.throughput(Throughput::Elements((traces.len() * 24) as u64));

	for rounds_per_iter in [1, 24] {
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
}

criterion_group!(keccak_ntt, bench_ntt_lookup, bench_keccak_residuals);
criterion_main!(keccak_ntt);
