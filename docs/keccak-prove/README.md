# Keccak Proving Notes

These notes track the Keccak-specific proving work for Binius64.

- `full-protocol.md` is the current working protocol note.
- `design.md` preserves the earlier fused-round design.

The implementation home is `crates/keccak-prove`.

The intended first implementation target is Binius-native and performance-oriented:

- prove many Keccak-f[1600] permutations in parallel;
- adapt the production BitAnd and shift-reduction machinery rather than rebuilding a toy protocol;
- use the bit-axis NTT lookup path from the start;
- specialize the generic shifted-word machinery to Keccak's fixed theta, rho, pi, chi, and iota structure;
- benchmark against the existing generic Keccak circuit path early.

Serial claim-reduction tests are still useful as algebraic guardrails, but they are no longer the v0 protocol target.

## Implementation checkpoint

Current implementation lives in `crates/keccak-prove` on the `quang/keccak-prove` branch.

What has been implemented so far:

- native Keccak-f[1600] trace construction;
- fully unrolled theta/rho/pi and chi/iota hot paths, avoiding `%5` and `/5` in round execution;
- BitAnd-style chi/iota residual words for all 25 lanes;
- Keccak-local byte NTT lookup over the 64-bit lane axis, matching the production BitAnd lookup-domain shape;
- upper-half extension-domain residual evaluation;
- sequential and Rayon-parallel first-round accumulation;
- a production-shaped packed small-field accumulator that keeps lane weights in the NTT field and widens to `B128` only after accumulation;
- first-round message construction with the original 64-point domain set to zero and the shifted upper-half domain filled by the prover;
- first-round claim extrapolation at the verifier challenge;
- a folded post-first-challenge claim check that directly folds each lane word and matches the verifier's extrapolated next claim;
- explicit folded `P`, `Q`, `R`, next-state, and iota columns over padded `(round_trace, lane)` rows;
- a first full Spartan outer pass after the bit-axis univariate skip, using the production `QuadraticMleCheckProver` for the remaining degree-2 MLE-check rounds;
- transcripted prover and verifier replay for the chi/iota segment, following the production BitAnd channel flow;
- verifier-field row challenge support for the transcripted first-round message, so the row weights are derived from the same MLE-check point used by the verifier.

Important lessons from the initial implementation:

- The NTT lookup setup matches production BitAnd closely. Lookup precompute has repeatedly measured in the same band as `binius-prover` BitAnd lookup precompute.
- The first generic `B128`-weighted accumulator was not the right hot-loop shape. Keeping weights packed in the small NTT field made the Keccak accumulator production-shaped and removed the apparent round-message bottleneck.
- The first folded-column builder was also the wrong hot-loop shape: it directly folded every 64-bit word against all 64 Lagrange values. Switching to the same bytewise lookup transform used by production BitAnd reduced folded-column construction by roughly an order of magnitude on the 128-permutation benchmark.
- Avoiding an extra scalar-column-to-`FieldBuffer` copy before `QuadraticMleCheckProver` matters. The post-skip outer pass now consumes the folded vectors directly into `FieldBuffer`s, matching the BitAnd reduction shape more closely.
- The post-skip outer pass should use the same three-column shape as production BitAnd. Combining `R + next + iota` into a single folded `C` column lets the quadratic MLE-check prove `P * Q - C`, rather than carrying five columns through every remaining round.
- The folded-column builder should be parallel like production word folding. Preallocating the padded columns and filling `(round_trace, lane)` chunks with Rayon removed the serial fold bottleneck.
- Microbenchmarks must report constraints processed per iteration. Earlier raw timings compared very different workloads and overstated the gap to production BitAnd.
- The current first-round claim path is still only the first sumcheck-message layer, not an end-to-end Keccak proof. It does not yet include folding through subsequent sumcheck rounds, transcript integration, committed input/output boundary openings, or full proof serialization.
- The post-skip outer pass uses the MLE-check identity from production Binius, so each remaining round message is tied to the current zerocheck coordinate by `(1 - alpha) r(0) + alpha r(1)`, not by the vanilla `r(0) + r(1)` sumcheck identity.
- The transcripted segment now links the first-round message to the post-skip MLE-check using the verifier's full row point. The packed small-field path is still useful as a fast benchmark/reference, but it is no longer a correctness restriction for transcript replay.
- Production Shift reduces BitAnd and IntMul claims through a two-phase protocol: a `g * h` sumcheck over bit/shift variables, followed by a bivariate product against the folded committed witness and the monster multilinear. The Keccak linear-layer pushback should be expressed in that shape. A custom direct claim rewrite over the current `(round_trace, lane)` rows would be a different, less production-faithful design.

Observed benchmark checkpoint on this machine:

```text
cargo bench -p binius-keccak-prove --bench keccak_ntt -- keccak_first_round_claim_scale
```

With Criterion `sample_size(10)`, the current `first_round_claim_small_par_distinct`
path measured approximately:

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

Thus the first measured size crossing roughly 500 ms by median for the current first-round-claim path is **131,072 Keccak-f permutations**.

The scale benchmark now allocates distinct traces and weights for each measured permutation count. For larger totals it uses distinct 65,536-permutation chunks and accumulates the chunk claims in one measured iteration, so the benchmark no longer depends on replaying the same 2048-permutation batch.
It is still a first-round benchmark rather than an end-to-end proof benchmark.

The first benchmark for the post-skip Spartan outer pass is:

```text
cargo bench -p binius-keccak-prove --bench keccak_ntt -- keccak_spartan_outer
```

On 128 Keccak-f permutations, the initial checkpoint measured approximately:

| Step | Median time |
|---|---:|
| Folded outer columns | 0.412 ms |
| Folded outer claim | 0.586 ms |
| Prove after univariate skip | 8.16 ms |

This benchmark includes the remaining degree-2 outer rounds after the bit-axis univariate skip, but still does not include transcript serialization, verifier replay, or boundary-opening reductions.

The post-skip outer pass scale benchmark is:

```text
cargo bench -p binius-keccak-prove --bench keccak_ntt -- keccak_spartan_outer_scale
```

With Criterion `sample_size(10)`, the current `prove_after_univariate_skip_distinct`
path measured approximately:

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

Thus the first measured size crossing roughly 500 ms by median for the current post-skip outer pass is **65,536 Keccak-f permutations**.

For a closer production BitAnd comparison:

```text
cargo bench -p binius-prover --bench and_reduction -- "full zerocheck"
```

The production BitAnd full-zerocheck benchmark measured approximately **24.8 ms** at `2^27` rows, reported as about **84.5M word constraints/s**. The Keccak post-skip outer pass at 55,924 permutations processes almost exactly one full padded row domain, `55,924 * 24 * 25 = 33.55M` folded lane constraints, in about **264 ms**, or about **126.9M constraints/s**. Just after the next padding cliff, 65,536 permutations processes `39.32M` folded lane constraints in about **524 ms**, or about **75.1M constraints/s**. This says the hot path is now close to production BitAnd; the remaining visible cliff is largely the padded Boolean row domain. It still excludes transcript serialization, verifier replay, linear-layer pushback, and boundary openings.

## Next steps

The next implementation milestone is to connect the transcripted chi/iota segment to the committed-witness path:

1. Decide the committed-witness layout for Keccak state words so input, round-boundary, and output lanes have stable word indices.
2. Represent the Keccak linear layer (`theta`, `rho`, `pi`) as production Shift-compatible shifted operands or an equivalent `KeyCollection`-style relation, then reuse the Shift two-phase reduction for the folded `P`, `Q`, and `C` claims.
3. Integrate boundary openings for committed input/output states through the same ring-switching and PCS opening path used after production Shift.
4. Add an end-to-end benchmark against the generic `binius-examples` Keccak circuit path, while keeping the current microbenches as regression tripwires.
