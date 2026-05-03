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
- a first full Spartan outer pass after the bit-axis univariate skip, using the production `QuadraticMleCheckProver` for the remaining degree-2 MLE-check rounds.

Important lessons from the initial implementation:

- The NTT lookup setup matches production BitAnd closely. Lookup precompute has repeatedly measured in the same band as `binius-prover` BitAnd lookup precompute.
- The first generic `B128`-weighted accumulator was not the right hot-loop shape. Keeping weights packed in the small NTT field made the Keccak accumulator production-shaped and removed the apparent round-message bottleneck.
- The first folded-column builder was also the wrong hot-loop shape: it directly folded every 64-bit word against all 64 Lagrange values. Switching to the same bytewise lookup transform used by production BitAnd reduced folded-column construction by roughly an order of magnitude on the 128-permutation benchmark.
- Avoiding an extra scalar-column-to-`FieldBuffer` copy before `QuadraticMleCheckProver` matters. The post-skip outer pass now consumes the folded vectors directly into `FieldBuffer`s, matching the BitAnd reduction shape more closely.
- Microbenchmarks must report constraints processed per iteration. Earlier raw timings compared very different workloads and overstated the gap to production BitAnd.
- The current first-round claim path is still only the first sumcheck-message layer, not an end-to-end Keccak proof. It does not yet include folding through subsequent sumcheck rounds, transcript integration, committed input/output boundary openings, or full proof serialization.
- The post-skip outer pass uses the MLE-check identity from production Binius, so each remaining round message is tied to the current zerocheck coordinate by `(1 - alpha) r(0) + alpha r(1)`, not by the vanilla `r(0) + r(1)` sumcheck identity.

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
| Folded outer columns | 1.54 ms |
| Folded outer claim | 1.23 ms |
| Prove after univariate skip | 12.5 ms |

This benchmark includes the remaining degree-2 outer rounds after the bit-axis univariate skip, but still does not include transcript serialization, verifier replay, or boundary-opening reductions.

## Next steps

The next implementation milestone is to turn the current prover-driven pass into a transcripted, verifier-replayed segment:

1. Add transcript plumbing around the univariate message and the post-skip quadratic MLE-check rounds.
2. Add verifier-side replay for this one Keccak chi/iota segment, including the MLE-check round identity.
3. Integrate boundary openings for committed input/output states.
4. Push folded `P`, `Q`, and `R` claims backward through the Keccak linear layer (`pi`, `rho`, `theta`) instead of treating pre-chi lanes as terminal columns.
5. Add an end-to-end benchmark against the generic `binius-examples` Keccak circuit path, while keeping the current microbenches as regression tripwires.
