# binius-keccak-prove

Implementation home for the Binius-native Keccak-f[1600] proving path.

The v0 target is a performance-oriented prover that adapts production Binius64 machinery rather than starting from a serial toy protocol:

- BitAnd-style $K=3$ outer relation for chi;
- bit-axis additive-NTT lookup from the start;
- Keccak-specialized replacement for generic shift reduction;
- parallel or segmented-parallel proving across batches and rounds;
- early benchmarks against the existing generic Keccak circuit path.

Protocol notes live in `../../docs/keccak-prove/`.

The first implementation milestone is a residual-zero NTT test for one Keccak round segment:

```text
cargo test -p binius-keccak-prove keccak_bitand_residual_matches_native_round
```

Run the full crate sanity suite with:

```text
cargo test -p binius-keccak-prove
```

The local microbenchmarks are:

```text
cargo bench -p binius-keccak-prove --bench keccak_ntt
```

Current implementation chain:

- committed `A`/`D` witness in the locked 32-word block layout;
- chi/iota Spartan outer in the production BitAnd shape;
- witness-only chi `P/Q/C` claims lowered through virtual `B` operands into committed `A`/`D`;
- `D` correctness encoded as degenerate AND rows;
- production Shift prover/verifier tests for both Keccak schemas.
- `v0` production path builds a normal Binius64 `ConstraintSystem` with public all-one/iota
  constants, proves it with the production `Prover`, and verifies through BitAnd, Shift,
  ring-switching, and PCS.
- v0 constraint rows are laid out as 32 slots per `(permutation, round)`:
  25 chi rows, 5 `D` correctness rows, and 2 padding rows. This keeps the row space tensor-shaped
  for verifier specialization without changing the current power-of-two proving cliffs.

They track three baselines while the prover path is still being assembled:

- `direct_lagrange_word`: slow, obviously correct evaluation of one 64-bit lane on the shifted upper-half domain.
- `byte_lookup_word` and `keccak_lookup_precompute`: the Keccak-local byte NTT lookup path.
- `upper_half_round_message_seq` and `upper_half_round_message_par`: Keccak chi/iota extension-domain accumulation with big-field lane weights.
- `upper_half_round_message_small_seq` and `upper_half_round_message_small_par`: the production-shaped variant that keeps lane weights packed in the NTT field and widens only after accumulation.
- `first_round_claim_small_par`: the production-shaped first-round flow, including upper-half accumulation, full 128-point message construction with zero base-domain values, and extrapolation at a verifier challenge.
- `keccak_first_round_claim_scale`: a distinct-data batch sweep up to 196,608 effective permutations with Criterion `sample_size(10)` to track where first-round work crosses roughly 500 ms.
- `keccak_spartan_outer`: the folded-column construction, folded outer claim, and remaining post-skip Spartan outer rounds for 128 permutations.
- `keccak_spartan_outer_scale`: a distinct-data batch sweep up to 65,536 permutations to track where the post-skip outer pass crosses roughly 500 ms.
- `production_bitand_lookup_precompute`: the existing Binius64 BitAnd lookup setup, using the same domain shape.
- `production_bitand_reference`: the existing Binius64 BitAnd univariate round-message hot path.
- `keccak_v0_production_path`: the end-to-end v0 production proof and verifier path.
- `keccak_v0_structured_verifier`: generic Shift monster evaluation versus the first
  tensor-structured Keccak verifier prototype.

Criterion throughput is reported as constraints processed per iteration:

- one 64-bit lane for word lookup benchmarks;
- 25 lane constraints per Keccak round residual benchmark;
- `128 * 24 * 25` lane constraints for the Keccak round-message accumulator;
- `2^(log_num_rows - 6)` word constraints for the production BitAnd reference.

The scale benchmark can be run directly with:

```text
cargo bench -p binius-keccak-prove --bench keccak_ntt -- keccak_first_round_claim_scale
```

On the initial implementation machine, the current `first_round_claim_small_par_distinct` path crossed roughly 500 ms by median at 131,072 Keccak-f permutations. The scale benchmark allocates distinct traces and weights for each measured permutation count, using distinct 65,536-permutation chunks for larger totals.

The initial post-skip outer benchmark is:

```text
cargo bench -p binius-keccak-prove --bench keccak_ntt -- keccak_spartan_outer
```

On 128 permutations, `prove_after_univariate_skip` measured about 8.16 ms median after switching folded-column construction to the same bytewise lookup transform used by production BitAnd, combining `R + next + iota` into the BitAnd-shaped third column, filling folded columns in parallel, and avoiding an extra scalar-column copy into `FieldBuffer`s. This covers the remaining degree-2 Spartan outer rounds after the bit-axis univariate skip, but not transcript serialization, verifier replay, or boundary-opening reductions.

The post-skip outer scale benchmark is:

```text
cargo bench -p binius-keccak-prove --bench keccak_ntt -- keccak_spartan_outer_scale
```

On the initial implementation machine, `prove_after_univariate_skip_distinct` crossed roughly 500 ms by median at 65,536 Keccak-f permutations. A well-packed point just below the previous padded-row cliff, 55,924 permutations, measured about 264 ms median.

For production BitAnd comparison:

```text
cargo bench -p binius-prover --bench and_reduction -- "full zerocheck"
```

On the same machine, production BitAnd full zerocheck measured about 24.8 ms at `2^27` rows, or about 84.5M word constraints/s. The Keccak post-skip outer pass measured about 126.9M folded lane constraints/s at 55,924 permutations and about 75.1M folded lane constraints/s at 65,536 permutations after the next padding cliff.

For broader production comparisons, also run:

```text
cargo bench -p binius-prover --bench and_reduction
cargo bench -p binius-examples --bench keccak
```

The first is the optimized BitAnd outer-reduction path we are adapting. The second is the current generic Keccak circuit proving baseline.

Current full-path checkpoint on this machine:

```text
cargo bench -p binius-keccak-prove --bench keccak_ntt -- keccak_v0_production_path
HASH_MAX_BYTES=17408 LOG_INV_RATE=1 cargo bench -p binius-examples --bench keccak -- keccak_proof
```

For 128 Keccak-f permutations, v0 measured about 56.1 ms median proving time and 2.69 ms
verification time. The generic Keccak circuit benchmark at the matching 17,408-byte workload
measured about 61.1 ms proving time and 2.83 ms verification time.

The first structured verifier prototype is correctness-checked against the generic Shift verifier
matrix evaluation. It is not yet wired into the production verifier, but as a microbench it measured
about 17.1 ms versus 18.2 ms for generic Shift monster evaluation at 1,024 permutations. At 128
permutations it is still slower than the generic prebuilt-operand microbench, so the next verifier
step is a production hook plus more tensor-specific accumulation rather than dropping it in blindly.
