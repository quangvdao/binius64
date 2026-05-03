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

They track three baselines while the prover path is still being assembled:

- `direct_lagrange_word`: slow, obviously correct evaluation of one 64-bit lane on the shifted upper-half domain.
- `byte_lookup_word` and `keccak_lookup_precompute`: the Keccak-local byte NTT lookup path.
- `upper_half_round_message_seq` and `upper_half_round_message_par`: Keccak chi/iota extension-domain accumulation with big-field lane weights.
- `upper_half_round_message_small_seq` and `upper_half_round_message_small_par`: the production-shaped variant that keeps lane weights packed in the NTT field and widens only after accumulation.
- `first_round_claim_small_par`: the production-shaped first-round flow, including upper-half accumulation, full 128-point message construction with zero base-domain values, and extrapolation at a verifier challenge.
- `keccak_first_round_claim_scale`: a distinct-data batch sweep up to 196,608 effective permutations with Criterion `sample_size(10)` to track where first-round work crosses roughly 500 ms.
- `keccak_spartan_outer`: the folded-column construction, folded outer claim, and remaining post-skip Spartan outer rounds for 128 permutations.
- `production_bitand_lookup_precompute`: the existing Binius64 BitAnd lookup setup, using the same domain shape.
- `production_bitand_reference`: the existing Binius64 BitAnd univariate round-message hot path.

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

On 128 permutations, `prove_after_univariate_skip` measured about 68 ms median. This covers the remaining degree-2 Spartan outer rounds after the bit-axis univariate skip, but not transcript serialization, verifier replay, or boundary-opening reductions.

For broader production comparisons, also run:

```text
cargo bench -p binius-prover --bench and_reduction
cargo bench -p binius-examples --bench keccak
```

The first is the optimized BitAnd outer-reduction path we are adapting. The second is the current generic Keccak circuit proving baseline.
