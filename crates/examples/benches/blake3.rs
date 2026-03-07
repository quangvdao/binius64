// Copyright 2026 The Binius Developers

mod utils;

use std::alloc::System;

use binius_examples::circuits::blake3::{Blake3Example, Instance, Params};
use criterion::{Criterion, Throughput, criterion_group, criterion_main};
use peakmem_alloc::PeakAlloc;
use utils::{ExampleBenchmark, HashBenchConfig, print_benchmark_header, run_cs_benchmark};

#[global_allocator]
static BLAKE3_PEAK_ALLOC: PeakAlloc<System> = PeakAlloc::new(System);

struct Blake3Benchmark {
	max_bytes: usize,
	count: usize,
	log_inv_rate: usize,
}

impl Blake3Benchmark {
	fn new(max_bytes: usize, count: usize) -> Self {
		let config = HashBenchConfig::from_env();
		Self {
			max_bytes,
			count,
			log_inv_rate: config.log_inv_rate,
		}
	}
}

impl ExampleBenchmark for Blake3Benchmark {
	type Params = Params;
	type Instance = Instance;
	type Example = Blake3Example;

	fn create_params(&self) -> Self::Params {
		Params {
			max_bytes: Some(self.max_bytes),
			count: self.count,
		}
	}

	fn create_instance(&self) -> Self::Instance {
		Instance {
			message_len: Some(self.max_bytes),
			message_string: None,
		}
	}

	fn bench_name(&self) -> String {
		let chunks = self.max_bytes.div_ceil(1024).max(1);
		format!("{}x_{chunks}ch_{}B", self.count, self.max_bytes)
	}

	fn throughput(&self) -> Throughput {
		Throughput::Bytes((self.count * self.max_bytes) as u64)
	}

	fn proof_description(&self) -> String {
		let chunks = self.max_bytes.div_ceil(1024).max(1);
		format!(
			"{} invocations x {} bytes ({} chunks)",
			self.count, self.max_bytes, chunks
		)
	}

	fn log_inv_rate(&self) -> usize {
		self.log_inv_rate
	}

	fn print_params(&self) {
		let chunks = self.max_bytes.div_ceil(1024).max(1);
		let blocks = self.max_bytes.div_ceil(64);
		let tree_depth = if chunks > 1 {
			(chunks as f64).log2().ceil() as usize
		} else {
			0
		};
		let compressions = if chunks > 1 {
			chunks + (chunks - 1)
		} else {
			blocks
		};
		let mut params_list = vec![
			("Invocations".to_string(), format!("{}", self.count)),
			(
				"Message per invocation".to_string(),
				format!(
					"{} bytes ({} chunks, {} blocks)",
					self.max_bytes, chunks, blocks
				),
			),
		];
		if chunks > 1 {
			params_list.push((
				"Tree hashing".to_string(),
				format!(
					"depth {}, {} compressions ({} chunk + {} parent)",
					tree_depth, compressions, chunks, chunks - 1
				),
			));
		}
		params_list.push((
			"Log inverse rate".to_string(),
			self.log_inv_rate.to_string(),
		));
		print_benchmark_header("BLAKE3", &params_list);
	}
}

fn bench_blake3_tree_hash(c: &mut Criterion) {
	// Tree hashing: sweep message sizes (each > 1 chunk triggers tree mode)
	let message_sizes = [2048, 4096, 8192, 16384, 32768, 65536];

	for &max_bytes in &message_sizes {
		let chunks = max_bytes / 1024;
		let benchmark = Blake3Benchmark::new(max_bytes, 1);
		run_cs_benchmark(
			c,
			benchmark,
			&format!("blake3_tree_{chunks}ch"),
			&BLAKE3_PEAK_ALLOC,
		);
	}

	// Batch tree hashing: 4 KiB messages (4 chunks) with increasing batch counts
	let counts = [1, 16, 64, 256];
	for &count in &counts {
		let benchmark = Blake3Benchmark::new(4096, count);
		run_cs_benchmark(
			c,
			benchmark,
			&format!("blake3_tree_4ch_{count}x"),
			&BLAKE3_PEAK_ALLOC,
		);
	}
}

criterion_group!(benches, bench_blake3_tree_hash);
criterion_main!(benches);
