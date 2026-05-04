# Keccak v0 Performance Notes

These notes record benchmark observations for the Binius-native Keccak v0 path.
They are intentionally operational: what was measured, what configuration mattered, and which
numbers are currently meaningful.

## Snapshot: 2026-05-03

The current full-path benchmark compares:

- v0 specialized path:

  ```text
  cargo bench -p binius-keccak-prove --bench keccak_ntt -- keccak_v0_production_path
  ```

- production generic Keccak path:

  ```text
  HASH_MAX_BYTES=17408 LOG_INV_RATE=1 cargo bench -p binius-examples --bench keccak -- keccak_proof
  ```

The 17,408-byte production benchmark corresponds to 128 Keccak-f permutations at 136 bytes per
permutation.

## Benchmark Inventory

The local Keccak benchmark suite is:

```text
cargo bench -p binius-keccak-prove --bench keccak_ntt
```

It tracks these paths while the prover is still being assembled:

- `direct_lagrange_word`: slow, obviously correct evaluation of one 64-bit lane on the shifted
  upper-half domain.
- `byte_lookup_word` and `keccak_lookup_precompute`: the Keccak-local byte NTT lookup path.
- `upper_half_round_message_seq` and `upper_half_round_message_par`: Keccak chi/iota
  extension-domain accumulation with big-field lane weights.
- `upper_half_round_message_small_seq` and `upper_half_round_message_small_par`: the
  production-shaped variant that keeps lane weights packed in the NTT field and widens only after
  accumulation.
- `first_round_claim_small_par`: the production-shaped first-round flow, including upper-half
  accumulation, full 128-point message construction with zero base-domain values, and
  extrapolation at a verifier challenge.
- `keccak_first_round_claim_scale`: a distinct-data batch sweep up to 196,608 effective
  permutations.
- `keccak_spartan_outer`: folded-column construction, folded outer claim, and the remaining
  post-skip Spartan outer rounds for 128 permutations.
- `keccak_spartan_outer_scale`: a distinct-data batch sweep up to 65,536 permutations for the
  post-skip outer pass.
- `production_bitand_lookup_precompute`: the existing Binius64 BitAnd lookup setup, using the same
  domain shape.
- `production_bitand_reference`: the existing Binius64 BitAnd univariate round-message hot path.
- `keccak_v0_production_path`: the end-to-end v0 production proof and verifier path.
- `keccak_v0_structured_verifier`: generic Shift monster evaluation versus the first
  tensor-structured Keccak verifier prototype.

Criterion throughput should be interpreted as constraints processed per iteration:

- one 64-bit lane for word lookup benchmarks;
- 25 lane constraints per Keccak round residual benchmark;
- `128 * 24 * 25` lane constraints for the Keccak round-message accumulator;
- `2^(log_num_rows - 6)` word constraints for the production BitAnd reference.

For broader production comparisons, also run:

```text
cargo bench -p binius-prover --bench and_reduction
cargo bench -p binius-examples --bench keccak
```

The first is the optimized BitAnd outer-reduction path we are adapting. The second is the current
generic Keccak circuit proving baseline.

## Implementation Performance Lessons

- The NTT lookup setup matches production BitAnd closely. Lookup precompute has repeatedly measured
  in the same band as `binius-prover` BitAnd lookup precompute.
- The first generic `B128`-weighted accumulator was not the right hot-loop shape. Keeping weights
  packed in the small NTT field made the Keccak accumulator production-shaped and removed the
  apparent round-message bottleneck.
- The first folded-column builder was also the wrong hot-loop shape: it directly folded every
  64-bit word against all 64 Lagrange values. Switching to the same bytewise lookup transform used
  by production BitAnd reduced folded-column construction by roughly an order of magnitude on the
  128-permutation benchmark.
- Avoiding an extra scalar-column-to-`FieldBuffer` copy before `QuadraticMleCheckProver` matters.
  The post-skip outer pass now consumes the folded vectors directly into `FieldBuffer`s, matching
  the BitAnd reduction shape more closely.
- Packed experiments must also avoid scalar-to-packed conversion inside the timed prover loop. The
  current fair packed benchmark materializes `PackedFoldedOuterColumns<OptimalPackedB128>` once and
  times cloned packed `FieldBuffer`s. Timings that call
  `prove_from_folded_columns_packed*_distinct` still include packing overhead and should only be
  used to quantify that overhead.
- The post-skip outer pass should use the same three-column shape as production BitAnd. Combining
  `R + next + iota` into a single folded `C` column lets the quadratic MLE-check prove `P * Q - C`,
  rather than carrying five columns through every remaining round.
- The folded-column builder should be parallel like production word folding. Preallocating the
  padded columns and filling `(round_trace, lane)` chunks with Rayon removed the serial fold
  bottleneck.
- Microbenchmarks must report constraints processed per iteration. Earlier raw timings compared
  very different workloads and overstated the gap to production BitAnd.
- The transcripted segment now links the first-round message to the post-skip MLE-check using the
  verifier's full row point. The packed small-field path is still useful as a fast
  benchmark/reference, but it is no longer a correctness restriction for transcript replay.
- The full v0 production path is now at least parity with the generic Keccak circuit benchmark at
  the 128-permutation checkpoint on the initial implementation machine. This comparison is not the
  final apples-to-apples story for SHA3/SHAKE wrappers, but it is the right first full-path sanity
  check because both sides run production proof generation and verification.
- Verifier specialization needs a real hook into Shift verification. The first tensor evaluator
  proves the math and is modestly faster at 1,024 permutations, but production verification still
  calls the generic evaluator. The next useful step is to expose a structured monster-evaluation
  hook or a v0 verifier path that replaces only `shift::check_eval`.

## Hardware and Raw-Power Expectation

The two useful machines in this checkpoint were:

| Machine | CPU | Threads | Relevant ISA | Expectation |
|---|---:|---:|---|---|
| local laptop | Apple M4 Max | 16 logical, 12 performance + 4 efficiency | NEON, AES, SHA3, AMX available to the platform | Strong single-machine baseline |
| `leopard` | AMD Ryzen 9 9950X | 16 cores, 32 threads | AVX2, AVX512, GFNI, VAES, VPCLMULQDQ | More raw throughput if compiled and scheduled well |

Raw arithmetic microbenchmarks support that expectation. With `RUSTFLAGS="-C target-cpu=native"`,
`leopard` exposes the x86 GFNI/VPCLMUL/AVX512 paths and is much faster on CLMUL-style kernels than
the loaded local laptop run. For example, the `ghash_google` kernel measured about 952M elements/s
on `leopard` native using `__m256i`, versus about 266M elements/s on the local run's aarch64 CLMUL
path.

The local laptop was not idle during the 2026-05-03 investigation: load average was about 118, with
an unrelated debug `setup` process using about 645% CPU. Today's local timings should therefore not
be used as a clean local baseline. The earlier clean local checkpoint remains the meaningful local
comparison for now.

## Native x86 Build Is Mandatory

On `leopard`, a default x86 build was misleading. The production benchmark's feature report showed:

```text
CPU: generic
Available CPU Instructions: fxsr, sse, sse2
```

That build does not enable the important Binius x86 paths. Native compilation showed:

```text
CPU: native
Available CPU Instructions: aes, avx, avx2, avx512*, gfni, pclmulqdq, vaes, vpclmulqdq, ...
```

Use this for meaningful `leopard` benchmarks:

```text
CARGO_TARGET_DIR=target-native RUSTFLAGS="-C target-cpu=native" ...
```

The separate target directory avoids mixing native and portable artifacts.

## Rayon Thread Count Is the Main Surprise

The first `leopard` result looked slow because the benchmark used Rayon defaults, which chose the
32-thread global pool. For the 128-permutation Keccak proof size, 32 threads is a severe slowdown on
the Ryzen 9950X. The likely cause is overhead and cross-CCD/cache traffic, with SMT adding more
contention than useful parallelism.

`leopard` is a 16-core / 32-logical-thread Ryzen 9950X with SMT enabled and two L3 caches:

```text
CPUs 0-7 and 16-23 share L3 0.
CPUs 8-15 and 24-31 share L3 1.
CPUs 16-31 are SMT siblings of CPUs 0-15.
```

With native compilation, v0 128-permutation proving measured:

| `RAYON_NUM_THREADS` | v0 prove time |
|---:|---:|
| 1 | 20.5 ms |
| 8 | 10.6 ms |
| 16 | 21.3 ms |
| 32 | 104.8 ms |

Production generic Keccak showed the same shape:

| `RAYON_NUM_THREADS` | production prove time |
|---:|---:|
| 1 | 21.7 ms |
| 8 | 11.3 ms |
| 16 | 21.6 ms |
| 32 | 104.4 ms |

So the apparent server regression was not v0-specific. It is a production-path scheduling issue on
this workload size.

An experimental diagnostic decoupled the BaseFold/NTT share count from the Rayon pool size by
forcing the NTT to use 8, 16, or 32 logical shares while keeping the Rayon pool fixed. That did not
materially change the full-proof timings:

| Rayon threads | NTT logical shares | v0 128-perm prove |
|---:|---:|---:|
| 8 | implicit / 8 / 16 / 32 | 10.5-10.8 ms |
| 16 | implicit / 8 / 16 / 32 | 20.5-21.1 ms |
| 32 | implicit / 8 / 16 / 32 | 103.8-105.0 ms |

This falsifies the initial hypothesis that the cliff was primarily from the
`current_num_threads().ilog2()` NTT-share heuristic. The issue is more general Rayon pool behavior
on this prover path.

The phase-level sweep was more nuanced:

| Benchmark | 8 threads | 16 threads | 32 threads | Interpretation |
|---|---:|---:|---:|---|
| Keccak first-round claim, 128 perms | 185 us | 123 us | 130 us | small round-message work can use 16 threads |
| folded outer claim, 128 perms | 388 us | 516 us | 594 us | outer folding overhead grows with thread count |
| post-skip Spartan outer, 128 perms | 1.14 ms | 2.59 ms | 3.75 ms | the outer pass is the first clear regression |
| full v0 proof, 128 perms | 10.7 ms | 20.6 ms | 102.3 ms | full proof amplifies the scheduling cliff |

At larger full-proof sizes, 16 physical cores can help, but SMT still does not:

| Keccak-f permutations | 1 thread | 2 threads | 4 threads | 8 threads | 12 threads | 16 threads | 20 threads | 24 threads | 28 threads | 32 threads |
|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| 2,048 | 310.6 ms | 170.8 ms | 105.6 ms | 82.4 ms | 84.1 ms | 90.6 ms | 97.8 ms | 110.7 ms | 133.6 ms | 245.0 ms |
| 8,192 | 1.276 s | 704.4 ms | 450.8 ms | 349.6 ms | 327.8 ms | 327.7 ms | 338.1 ms | 355.6 ms | 384.7 ms | 553.9 ms |

The 8k-32k target range changes the practical recommendation. A focused 2026-05-03 native
`leopard` sweep used:

```text
CARGO_TARGET_DIR=target-native RUSTFLAGS="-C target-cpu=native"
KECCAK_V0_PRODUCTION_PERMS=8192,16384,32768
```

and measured full v0 proving:

| Keccak-f permutations | 8 threads | 12 threads | 16 threads | 24 threads | 32 threads |
|---:|---:|---:|---:|---:|---:|
| 8,192 | 349.7 ms | 327.0 ms | 326.7 ms | 356.3 ms | 554.2 ms |
| 16,384 | 720.4 ms | 664.1 ms | 655.0 ms | 691.2 ms | 932.8 ms |
| 32,768 | 1.465 s | 1.354 s | 1.316 s | 1.352 s | 1.616 s |

So for the actual 8k-32k target range, the best default is 12-16 physical cores, not 8. The
remaining wall is specifically SMT-heavy scheduling. 24 logical workers already regresses at
8k/16k, and 32 logical workers is bad across the range.

Affinity checks at 8,192 permutations showed:

| Config | v0 prove time |
|---|---:|
| 8 threads pinned to CPUs 0-7, one physical CCD | 349.1 ms |
| 8 threads pinned to CPUs 8-15, one physical CCD | 355.6 ms |
| 16 threads pinned to CPUs 0-15, all physical cores | 334.7 ms |
| 16 threads pinned to CPUs 0-7 and 16-23, one CCD plus SMT siblings | 406.3 ms |
| 32 threads pinned to CPUs 0-31, all logical threads | 543.5 ms |

At 32,768 permutations, affinity checks reinforced the same diagnosis:

| Config | v0 prove time |
|---|---:|
| 16 workers pinned to CPUs 0-15, all physical cores | 1.334 s |
| 16 workers pinned to CPUs 0-7 and 16-23, one CCD plus SMT siblings | 1.522 s |
| 16 workers pinned to CPUs 0-11, 12 physical cores | 1.376 s |

This is the strongest current evidence for the root cause. Sumcheck and the surrounding proof are
parallelizable, but this implementation streams large packed buffers and synchronizes each round.
Physical cores help once the instance is large enough; SMT siblings mostly compete for the same
execution ports, private caches, and memory bandwidth while adding scheduler traffic. The result is
a saturation point around 12-16 physical cores on `leopard`, not a fundamental sumcheck limit.

The practical rule for small checkpoints remains:

```text
CARGO_TARGET_DIR=target-native RUSTFLAGS="-C target-cpu=native" RAYON_NUM_THREADS=8 ...
```

for 128- and 2,048-permutation checkpoints. For the 8k-32k target range, use:

```text
CARGO_TARGET_DIR=target-native RUSTFLAGS="-C target-cpu=native" RAYON_NUM_THREADS=12 ...
```

or:

```text
CARGO_TARGET_DIR=target-native RUSTFLAGS="-C target-cpu=native" RAYON_NUM_THREADS=16 ...
```

Use 12 threads around 8k if minimizing tail latency; use 16 threads by 16k-32k. Avoid 24 and 32
threads for this path unless a future implementation change specifically proves otherwise.

The machine-level follow-up is to rerun a small subset with the CPU governor set to `performance`
if we get root access. During this sweep all CPUs reported the `powersave` governor. That may affect
absolute timings, but it does not explain why 32 logical threads perform worse than 12-16 on the
same machine and build.

## Packed Fused Outer Experiment

The scalar fused/adaptive experiment was the wrong implementation shape. It fused bind and reduce,
but it gave up the packed `FieldBuffer` representation that production Binius paths rely on. A fair
packed experiment now uses:

```text
pack_folded_outer_columns::<B128, OptimalPackedB128>(...)
prove_spartan_outer_from_packed_folded_columns_with_claim_fused(...)
```

This materializes the packed columns once, outside the timed prover loop. The older
`prove_from_folded_columns_packed*_distinct` benchmark names still repack scalar columns each
iteration and should not be used as the headline result.

The packed fused pass uses packed `FieldBuffer<P>` columns, packed `eq` weights, a fused bind-then-
reduce prefix, and then switches back to the generic packed quadratic prover for the small tail.
The default tail switchover is:

```text
KECCAK_PACKED_FUSED_MIN_REDUCE_WORDS=65536
```

This can be overridden for threshold sweeps.

On `leopard` with native x86 codegen and the 8-thread setting that works best for small-to-medium
Keccak proof sizes:

```text
CARGO_TARGET_DIR=target-native RUSTFLAGS="-C target-cpu=native" RAYON_NUM_THREADS=8 KECCAK_OUTER_WORKERS=8
```

the folded-column outer-only comparison measured:

| Keccak-f permutations | Generic folded columns | Prepacked fused packed columns | Result |
|---:|---:|---:|---|
| 128 | 1.152 ms | 0.766 ms | fused packed is ~1.5x faster |
| 8,192 | 131.1 ms | 126.4 ms | fused packed is ~3.6% faster |

The larger target-range sweep used the new bench filter:

```text
KECCAK_SPARTAN_OUTER_SCALE_PERMS=8192,16384,32768
```

At these sizes, prepacked fused is a small but consistent improvement over prepacked generic. The
thread-count effect is larger than the fused-vs-generic effect:

| Threads | 8,192 generic | 8,192 fused | 16,384 generic | 16,384 fused | 32,768 generic | 32,768 fused |
|---:|---:|---:|---:|---:|---:|---:|
| 8 | 125.6 ms | 124.9 ms | 260.9 ms | 255.3 ms | 501.4 ms | 483.7 ms |
| 12 | 133.9 ms | 130.6 ms | 253.3 ms | 242.1 ms | 493.0 ms | 472.2 ms |
| 16 | 127.6 ms | 126.4 ms | 252.1 ms | 243.7 ms | 496.2 ms | 475.7 ms |
| 24 | 133.3 ms | 128.9 ms | 254.1 ms | 243.7 ms | 500.4 ms | 479.6 ms |
| 32 | 131.9 ms | 127.8 ms | 267.0 ms | 254.1 ms | 503.5 ms | 483.3 ms |

For the outer pass alone, 8 threads is still best at 8,192, while 12-16 threads catch up around
16,384 and 32,768. The full proof benefits more clearly from 12-16 threads because other phases have
enough work to use the extra physical cores.

On the local laptop during a high-load run, timings were much noisier. A representative 8,192-perm
sample showed prepacked fused packed columns at about 38.7 ms median, but the generic folded-column
measurement varied enough that `leopard` should be treated as the cleaner checkpoint for this
specific experiment.

The production BitAnd reference hot path has a similar, though less severe, size-dependent shape:

| BitAnd reference size | 8 threads | 16 threads | 32 threads |
|---|---:|---:|---:|
| `log_rows=12` | 8.39 us | 19.6 us | 45.3 us |
| `log_rows=22` | 317 us | 243 us | 276 us |

This reinforces the same conclusion: larger production-shaped kernels can benefit from 16 physical
cores, but small kernels and SMT-heavy runs lose to scheduling/parallel overhead.

### Deeper Sumcheck Parallelism Diagnosis

There is no apparent mathematical limitation preventing sumcheck from using all physical cores.
The current limitation is the shape of the production Rayon implementation.

The current `QuadraticMleCheckProver` loop does the following every round:

1. `execute`: split multilinears, compute the round message with a fresh Rayon
   `into_par_iter().reduce()`.
2. Fiat-Shamir: send coefficients and sample the challenge.
3. `fold`: for each multilinear, call `fold_highest_var_inplace`, which itself uses Rayon again.

So each round pays multiple global Rayon joins, and the data is streamed in separate reduce and
bind passes. The live vector halves every round, but the Rayon pool size remains fixed. This is
exactly the regime where per-round dispatch, barriers, cache-line movement, and memory bandwidth
can dominate even though the first few rounds contain abundant parallel work.

The standalone `binius-prover` MLE-check bench shows the same shape without the rest of the Keccak
proof:

```text
cargo bench -p binius-prover --bench sumcheck -- "A.B-C/n_vars=..."
```

On `leopard` with native compilation:

| MLE-check size | 1 thread | 2 threads | 4 threads | 8 threads | 12 threads | 16 threads | 24 threads | 32 threads |
|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| `n_vars=12` | 91.9 us | 120 us | 150 us | 196 us | 406 us | 513 us | 591 us | 779 us |
| `n_vars=16` | 554 us | 441 us | 507 us | 563 us | 928 us | 1.27 ms | 1.47 ms | 1.93 ms |
| `n_vars=20` | 8.82 ms | 6.50 ms | 3.57 ms | 3.31 ms | 4.55 ms | 5.00 ms | 5.73 ms | pathological in Criterion |

The 32-thread `n_vars=20` Criterion run did not produce a result after several minutes and was
terminated. The useful takeaway is the monotone regression after 8 threads, not the missing exact
number.

This matches the findings in `/Users/quang.dao/Documents/Research/monomial-sumcheck-benchmarks`:

- `rayon::scope` and `par_iter` have a large fixed dispatch floor relative to small sumcheck
  rounds.
- A persistent pinned worker pool with one broadcast per sumcheck phase, static contiguous
  partitions, and adaptive active-worker shrinkage can move the crossover much earlier.
- The right shape is a parallel prefix plus a sequential tail: use many workers while
  `pairs_per_worker` is large, halve active workers as the live domain shrinks, then finish the
  tiny tail sequentially.
- A fused `bind_then_reduce_chunk` path is important because the next round's reduce reads exactly
  what the previous bind wrote. A split Rayon implementation writes the folded buffers, returns to
  the driver, then re-enters Rayon and rereads them.

The likely production direction is therefore not to tune Rayon harder, but to add a specialized
sumcheck execution backend:

1. Keep Rayon for broad embarrassingly parallel construction tasks where its work-stealing is a
   good fit.
2. Add an opt-in pinned/persistent sumcheck backend for hot MLE-check paths.
3. Start with `QuadraticMleCheckProver` or `BatchQuadraticMleCheckProver`, because Keccak v0's
   post-skip outer pass and production BitAnd both go through this family.
4. Use a trait shape similar to `sumcheck-parallel`'s `SumcheckRound`: `reduce_chunk`,
   `bind_chunk`, optional `bind_then_reduce_chunk`, `combine`, and a transcript/observe hook.
5. Pick active workers from problem size and working set, not from global logical CPU count:
   8 workers for small/medium, 12-16 physical-core workers for larger domains on `leopard`, and
   no SMT unless a kernel is proved compute-bound enough to benefit.

This should let us test the real theoretical parallelism of sumcheck. If a persistent pinned
backend still fails to use 16 physical cores on large domains, then we should suspect memory
bandwidth or field-operation throughput as the fundamental limiter. The current Rayon results do
not yet justify that conclusion.

### Keccak-Local Fused/Adaptive Sumcheck Experiment

The first fused-bind-then-reduce experiment is intentionally confined to
`crates/keccak-prove`. It leaves the production Binius64 MLE-check and BitAnd code untouched.

The experiment adds a Keccak-local post-skip Spartan outer prover with three scheduling changes:

- fuse the previous round's bind with the next round's reduction;
- carry the equality table forward by truncating its low/high halves instead of recomputing
  `eq_ind_partial_eval_scalars` each round;
- shrink the active worker count as the live domain halves, then finish the small tail
  sequentially.

The first local Mac result looked very promising, but it was not the fair conclusion. The
non-scale benchmark was run with `KECCAK_OUTER_WORKERS=8` but without setting
`RAYON_NUM_THREADS=8`, so the generic path paid the bad default Rayon configuration during
folded-column construction. That produced an apparent 128-permutation win:

| Path | Median time |
|---|---:|
| generic `QuadraticMleCheckProver` post-skip outer | 8.33 ms |
| Keccak-local fused/adaptive post-skip outer | 1.41 ms |

The corrected `leopard` native runs set both `RAYON_NUM_THREADS=8` and
`KECCAK_OUTER_WORKERS=8`, and added isolated `prove_from_folded_columns` benchmark entries so the
remaining outer sumcheck can be measured without folded-column construction.

With the first scoped-phase implementation, the fair `leopard` outer-only numbers were:

| Keccak-f permutations | Generic outer-only | Fused/adaptive outer-only |
|---:|---:|---:|
| 128 | 1.09 ms | 1.28 ms |
| 8,192 | 129.7 ms | 140.2 ms |

This falsifies the early interpretation: the current Keccak-local fused loop does not beat the
production-shaped generic MLE-check when measured fairly. The likely reason is that the generic
path keeps more of the field work in Binius's packed `FieldBuffer` machinery, while the experiment
uses scalar `B128` loops plus manual barriers.

Replacing per-phase scoped thread creation with a single persistent scoped worker team improved
the large case but made the small case worse:

| Keccak-f permutations | Generic outer-only | Persistent fused/adaptive outer-only |
|---:|---:|---:|
| 128 | 1.13 ms | 1.40 ms |
| 8,192 | 128.2 ms | 135.1 ms |

Worker and affinity checks on the persistent fused 8,192-permutation outer-only path showed:

| Config | Median time |
|---|---:|
| 8 workers | 141.2 ms in a later repeat |
| 12 workers | 140.2 ms |
| 16 workers | 140.5 ms |
| 8 workers, process pinned to CPUs 0-7 | 134.2 ms |

Pinning helps the persistent experiment somewhat, but not enough to beat the generic packed path.
The next useful direction is therefore not more process-level affinity tuning. It is to move the
fused bind/reduce idea into a packed-field backend, or into the production `QuadraticMleCheckProver`
shape where the bind and reduce can be fused without abandoning `FieldBuffer` packing.

### Packed Persistent-Worker Outer Experiment

The next experiment moved the persistent-worker idea back onto packed `FieldBuffer` columns:

```text
prove_spartan_outer_from_packed_folded_columns_with_claim_persistent_fused(...)
```

This path keeps the same algebraic shape as the fair packed fused prover. The only Fiat-Shamir
barrier is between round messages: once the previous challenge is available, each worker streams a
static packed-word range, binds the previous challenge, and accumulates the next round's message
from the just-bound values in the same pass. Round 0 is reduce-only, because there is no previous
challenge yet. After the persistent prefix, the prover binds the last produced challenge, truncates
the packed columns, and hands the small tail back to the generic packed quadratic prover.

Two knobs control the current experiment:

```text
KECCAK_PACKED_FUSED_MIN_REDUCE_WORDS=65536
KECCAK_PACKED_OUTER_MIN_WORDS_PER_WORKER=16384
```

The first is the tail switchover threshold shared with the Rayon fused packed path. The second is
the adaptive-worker shrinkage threshold for the persistent path. A `leopard` 32-worker threshold
sweep at 32,768 permutations measured:

| `KECCAK_PACKED_OUTER_MIN_WORDS_PER_WORKER` | Persistent packed fused, 32 workers |
|---:|---:|
| 8,192 | 495.3 ms |
| 16,384 | 482.3 ms |
| 32,768 | 498.6 ms |
| 65,536 | 493.3 ms |

So the default is currently 16,384 packed words per active worker. This confirms that adaptive
worker shrinkage matters, but it does not create a 32-thread breakthrough by itself.

On `leopard` native with `KECCAK_PACKED_OUTER_MIN_WORDS_PER_WORKER=4096`, the fair comparison was:

| Threads | 8,192 Rayon fused | 8,192 persistent | 16,384 Rayon fused | 16,384 persistent | 32,768 Rayon fused | 32,768 persistent |
|---:|---:|---:|---:|---:|---:|---:|
| 12 | 128.1 ms | 129.2 ms | 256.8 ms | 253.0 ms | 488.6 ms | 492.5 ms |
| 16 | 131.4 ms | 131.7 ms | 247.3 ms | 248.1 ms | 478.2 ms | 478.7 ms |
| 32 | 131.6 ms | 134.4 ms | 253.4 ms | 252.7 ms | 495.9 ms | 494.9 ms |

With the tuned 16,384-word threshold and 16 persistent workers, the target sweep measured:

| Keccak-f permutations | Persistent packed fused |
|---:|---:|
| 8,192 | 131.6 ms |
| 16,384 | 249.8 ms |
| 32,768 | 485.7 ms |

The result is useful but sobering: persistent workers plus fused bind/reduce reaches parity with
Rayon on the outer pass, and sometimes wins a few percent, but it does not solve the full 32-logical
thread utilization problem. The remaining likely bottleneck is not the Fiat-Shamir barrier itself;
it is the per-round streaming kernel's bandwidth/cache behavior. Every round reads and writes large
packed buffers, then the live domain halves. More logical workers increase L3/SMT contention and
barrier traffic before they add useful arithmetic throughput.

The strongest next implementation direction is therefore a larger grain of parallelism, not just a
different per-round scheduler. The safe target is still one Fiat-Shamir challenge per logical
sumcheck round, but with row variables ordered as `(chunk, local_row)` so high-to-low sumcheck binds
the local-row variables first. Chunk-local workers can then compute local round messages, aggregate
them into one global round message, sample one challenge, and bind every chunk with that same
challenge. Within one monolithic table layout, the next lower-level experiment would be a
hand-pinned worker pool with explicit CCD partitioning and hardware-counter measurements under a
lower `perf_event_paranoid` setting, to distinguish memory bandwidth from arithmetic-port
saturation.

### Coarse-Grain Multi-Instance Experiment

A follow-up benchmark tests that coarse-grain hypothesis without changing the production proof yet.
It is enabled only when:

```text
KECCAK_SPARTAN_OUTER_COARSE_CHUNK_PERMS=<chunk size>
KECCAK_SPARTAN_OUTER_COARSE_JOBS=<parallel chunk jobs>
```

The benchmark name is:

```text
prove_packed_persistent_fused_coarse_<jobs>_jobs_<chunk>_per_chunk
```

It splits the total permutation count into distinct chunk instances, proves each chunk with the
packed persistent fused outer prover, and runs several chunks concurrently. This is not yet an
integrated aggregation protocol; it is a performance probe for whether independent sumcheck
instances use the machine better than one large synchronized instance.

The answer is yes. On `leopard` native at 32,768 total permutations:

| Shape | Env sketch | Median time |
|---|---|---:|
| One synchronized instance | `RAYON_NUM_THREADS=16 KECCAK_OUTER_WORKERS=16`, 16,384 min words/worker | 490.1 ms |
| 4 chunks of 8,192 | 4 coarse jobs, 4 inner workers/chunk, 4,096 min words/worker | 308.5 ms |
| 8 chunks of 4,096 | 8 coarse jobs, 2 inner workers/chunk, 4,096 min words/worker | 297.9 ms |
| 8 chunks of 4,096 | 8 coarse jobs, 4 inner workers/chunk, 4,096 min words/worker | 304.9 ms |
| 16 chunks of 2,048 | 8 coarse jobs, 2 inner workers/chunk, 2,048 min words/worker | 276.7 ms |
| 16 chunks of 2,048 | 16 coarse jobs, 2 inner workers/chunk, 2,048 min words/worker | 303.7 ms |
| 32 chunks of 1,024 | 16 coarse jobs, 2 inner workers/chunk, 1,024 min words/worker | 286.4 ms |

The best measured shape in this sweep was therefore:

```text
KECCAK_SPARTAN_OUTER_SCALE_PERMS=32768
KECCAK_SPARTAN_OUTER_COARSE_CHUNK_PERMS=2048
KECCAK_SPARTAN_OUTER_COARSE_JOBS=8
KECCAK_PACKED_FUSED_MIN_REDUCE_WORDS=65536
KECCAK_PACKED_OUTER_MIN_WORDS_PER_WORKER=2048
RAYON_NUM_THREADS=8
KECCAK_OUTER_WORKERS=2
```

That is about 1.77x faster than the single 32,768-permutation outer instance in the same checkpoint.
It also explains why simply raising the worker count inside one instance was disappointing: the
machine wants more independent work, not more threads contending on the same shrinking buffers and
round barriers.

The implementation implication is important. For large Keccak batches, the better v0/v1 direction
is likely:

1. Partition the trace into proof chunks sized around the 1k-4k permutation range.
2. Prove chunk-local post-skip outer claims concurrently.
3. Aggregate the chunk claims with a small outer batching layer, rather than forcing one monolithic
   sumcheck instance over the entire padded row domain.

This would increase transcript and aggregation design work, but it is now the most concrete route
to using 16-32 hardware threads effectively.

### One-Fiat-Shamir Chunked Global Experiment

The independent multi-instance experiment above is a useful machine probe, but it is not the
protocol shape we want. The safe benchmark is:

```text
KECCAK_SPARTAN_OUTER_ONE_FS_CHUNK_PERMS=<chunk size>
KECCAK_SPARTAN_OUTER_ONE_FS_JOBS=<parallel chunk jobs>
```

with benchmark name:

```text
prove_packed_persistent_fused_one_fs_<jobs>_jobs_<chunk>_per_chunk
```

This benchmark keeps a single logical Fiat-Shamir timeline. It uses the variable order:

```text
global row index bits = [chunk bits | local-row bits]
sumcheck order        = local-row bits first, then chunk bits
```

For each local-row round:

```text
chunk 0 local round message
chunk 1 local round message
...
chunk m local round message
        |
        v
weighted sum by eq(chunk; r_chunk)
        |
        v
one global round message -> one challenge -> broadcast to all chunks
```

After all local-row variables are bound, the benchmark constructs a small chunk-axis
`P * Q - C` instance from the per-chunk multilinear evaluations and proves the remaining chunk
variables with the same single transcript structure.

There is an opt-in sanity check:

```text
KECCAK_SPARTAN_OUTER_ONE_FS_VERIFY=1
```

For a small 512-permutation run split into 128-permutation chunks, this check compared the one-FS
chunked final claim against a generic packed prover over the same permuted global table and passed.

On `leopard` native at 32,768 total permutations, the safe one-FS chunked benchmark measured:

| Shape | Env sketch | Median time |
|---|---|---:|
| One synchronized instance | `RAYON_NUM_THREADS=16 KECCAK_OUTER_WORKERS=16`, 16,384 min words/worker | 486.9 ms |
| One-FS chunked global | 16 chunks of 2,048; 8 jobs; 2 inner workers/chunk; 2,048 min words/worker | 276.0 ms |

The exact command shape was:

```text
KECCAK_SPARTAN_OUTER_SCALE_PERMS=32768
KECCAK_SPARTAN_OUTER_ONE_FS_CHUNK_PERMS=2048
KECCAK_SPARTAN_OUTER_ONE_FS_JOBS=8
KECCAK_PACKED_FUSED_MIN_REDUCE_WORDS=65536
KECCAK_PACKED_OUTER_MIN_WORDS_PER_WORKER=2048
RAYON_NUM_THREADS=8
KECCAK_OUTER_WORKERS=2
```

This is the result we wanted: almost the same speed as the independent multi-instance probe, but
without giving chunks independent Fiat-Shamir challenges. The implementation is still benchmark
infrastructure, not a finished verifier-integrated protocol path, but it strongly supports the
chunk-local table layout with one global transcript.

## Full-Path Checkpoint

For 128 Keccak-f permutations:

| Machine/config | v0 prove | production prove | v0 verify | production verify |
|---|---:|---:|---:|---:|
| local clean prior run | 56.1 ms | 61.1 ms | 2.69 ms | 2.83 ms |
| local loaded 2026-05-03 run | 125.2 ms | 219.7 ms | 5.10 ms | 5.79 ms |
| `leopard`, x86 generic, default Rayon | 118.2 ms | 121.2 ms | 26.1 ms | 39.2 ms |
| `leopard`, native, default Rayon | 104.1 ms | 105.3 ms | 5.15 ms | 9.18 ms |
| `leopard`, native, `RAYON_NUM_THREADS=8` | 10.6 ms | 11.3 ms | not measured in this sweep | not measured in this sweep |

The fair `leopard` proving comparison is the last row. Under that configuration, `leopard` is
substantially faster than the earlier clean local baseline, and v0 remains at production parity.

## Initial Local Scale Checkpoints

These numbers were captured on the initial implementation machine before the `leopard` scheduling
investigation. They remain useful as historical local checkpoints, but cross-machine runs should
follow the configuration rules above.

### First-Round Claim Scale

Command:

```text
cargo bench -p binius-keccak-prove --bench keccak_ntt -- keccak_first_round_claim_scale
```

With Criterion `sample_size(10)`, `first_round_claim_small_par_distinct` measured:

| Effective permutations | Median time |
|---:|---:|
| 2,048 | 8.78 ms |
| 4,096 | 21.1 ms |
| 8,192 | 35.9 ms |
| 16,384 | 76.0 ms |
| 32,768 | 155 ms |
| 65,536 | 411 ms |
| 131,072 | 983 ms |
| 196,608 | 1.50 s |

The first measured size crossing roughly 500 ms by median for this first-round-claim path was
131,072 Keccak-f permutations.

The scale benchmark allocates distinct traces and weights for each measured permutation count. For
larger totals it uses distinct 65,536-permutation chunks and accumulates the chunk claims in one
measured iteration, so the benchmark no longer depends on replaying the same 2,048-permutation
batch. It is still a first-round benchmark rather than an end-to-end proof benchmark.

### Post-Skip Spartan Outer

Command:

```text
cargo bench -p binius-keccak-prove --bench keccak_ntt -- keccak_spartan_outer
```

On 128 Keccak-f permutations, the checkpoint measured:

| Step | Median time |
|---|---:|
| Folded outer columns | 0.412 ms |
| Folded outer claim | 0.586 ms |
| Prove after univariate skip | 8.16 ms |

This benchmark includes the remaining degree-2 outer rounds after the bit-axis univariate skip, but
does not include transcript serialization, verifier replay, or boundary-opening reductions.

The scale command is:

```text
cargo bench -p binius-keccak-prove --bench keccak_ntt -- keccak_spartan_outer_scale
```

With Criterion `sample_size(10)`, `prove_after_univariate_skip_distinct` measured:

| Keccak-f permutations | Median time |
|---:|---:|
| 128 | 11.5 ms |
| 256 | 18.6 ms |
| 512 | 19.4 ms |
| 1,024 | 28.5 ms |
| 2,048 | 44.9 ms |
| 4,096 | 85.0 ms |
| 8,192 | 159 ms |
| 12,288 | 179 ms |
| 16,384 | 262 ms |
| 24,576 | 366 ms |
| 32,768 | 252 ms |
| 55,924 | 264 ms |
| 65,536 | 524 ms |

The first measured size crossing roughly 500 ms by median for this post-skip outer pass was 65,536
Keccak-f permutations.

### Production BitAnd Comparison

Command:

```text
cargo bench -p binius-prover --bench and_reduction -- "full zerocheck"
```

Production BitAnd full zerocheck measured about 24.8 ms at `2^27` rows, or about 84.5M word
constraints/s. The Keccak post-skip outer pass at 55,924 permutations processes almost exactly one
full padded row domain, `55,924 * 24 * 25 = 33.55M` folded lane constraints, in about 264 ms, or
about 126.9M constraints/s. Just after the next padding cliff, 65,536 permutations processes
`39.32M` folded lane constraints in about 524 ms, or about 75.1M constraints/s.

This says the hot path is now close to production BitAnd; the remaining visible cliff is largely the
padded Boolean row domain. It still excludes transcript serialization, verifier replay,
linear-layer pushback, and boundary openings.

### Structured Verifier Prototype

Command:

```text
cargo bench -p binius-keccak-prove --bench keccak_ntt -- keccak_v0_structured_verifier
```

This compares generic Shift monster matrix evaluation with the v0 tensor evaluator:

| Permutations | Generic | Structured |
|---:|---:|---:|
| 128 | 1.71 ms | 2.13 ms |
| 1,024 | 18.2 ms | 17.1 ms |

The structured evaluator is therefore not a universal drop-in win yet, but it validates the tensor
approach and starts to win at larger batches.

## Revised 500 ms Scale Point

The original `leopard` 500 ms crossing was measured with the bad default 32-thread pool and landed
near 5,462 permutations. That is no longer the right operational number.

With native compilation and `RAYON_NUM_THREADS=8`, the v0 full production path measured:

| Keccak-f permutations | v0 prove time |
|---:|---:|
| 2,048 | 81.9 ms |
| 5,462 | 325.5 ms |
| 8,192 | 351.1 ms |
| 10,000 | 366.0 ms |
| 10,485 | 372.5 ms |
| 10,486 | 533.1 ms |
| 10,922 | 535.2 ms |
| 10,923 | 668.3 ms |
| 16,384 | 717.5 ms |

The first clean point past 500 ms under the good `leopard` configuration is therefore 10,486
permutations.

The jump at 10,486 is a padding cliff. The committed witness has:

```text
32 public words + 800 committed words per permutation
```

The witness committed length crosses the `2^23` power-of-two boundary between 10,485 and 10,486
permutations:

```text
32 + 800 * 10,485 = 8,388,032 <= 2^23
32 + 800 * 10,486 = 8,388,832 >  2^23
```

The next constraint-row padding boundary is between 10,922 and 10,923 permutations because v0 has
768 constraint rows per permutation:

```text
768 * 10,922 = 8,388,096 <= 2^23
768 * 10,923 = 8,388,864 >  2^23
```

This explains the second jump around 10,923 permutations.

## Practical Benchmark Rules

- Always report compiler target and Rayon thread count for cross-machine numbers.
- On x86 servers, use `RUSTFLAGS="-C target-cpu=native"` unless the goal is portable-binary
  performance.
- On `leopard`, start with `RAYON_NUM_THREADS=8` for 128-2,048 permutation checkpoints.
- On `leopard`, use `RAYON_NUM_THREADS=12` or `16` for the 8k-32k target range. Avoid 24 and 32
  unless a new implementation demonstrates a win.
- Treat local laptop numbers as invalid if load average is high or another CPU-heavy job is active.
- When reporting scale sweeps, annotate power-of-two cliffs for both committed witness length and
  constraint-row length.
- Use `KECCAK_SPARTAN_OUTER_SCALE_PERMS` and `KECCAK_V0_PRODUCTION_PERMS` to restrict Criterion to
  the sizes being studied so large sweeps do not build unrelated instances.
- Continue to benchmark v0 against the production Keccak path, not only microbenchmarks, because
  the thread-count cliff affected both paths.
