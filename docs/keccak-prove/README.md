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
- a folded post-first-challenge claim check that directly folds each lane word and matches the verifier's extrapolated next claim.

Important lessons from the initial implementation:

- The NTT lookup setup matches production BitAnd closely. Lookup precompute has repeatedly measured in the same band as `binius-prover` BitAnd lookup precompute.
- The first generic `B128`-weighted accumulator was not the right hot-loop shape. Keeping weights packed in the small NTT field made the Keccak accumulator production-shaped and removed the apparent round-message bottleneck.
- Microbenchmarks must report constraints processed per iteration. Earlier raw timings compared very different workloads and overstated the gap to production BitAnd.
- The current first-round claim path is still only the first sumcheck-message layer, not an end-to-end Keccak proof. It does not yet include folding through subsequent sumcheck rounds, transcript integration, committed input/output boundary openings, or full proof serialization.

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

## Next steps

The next implementation milestone is to move from "first-round message and next claim" to "full first sumcheck segment":

1. Turn the folded post-first-challenge claim into explicit folded operand columns for `P`, `Q`, `R`, next-state, and iota, rather than only a direct consistency check.
2. Run the remaining quadratic sumcheck rounds for the Keccak chi/iota relation using the next claim produced by `par_first_round_claim_small_weights`.
3. Add transcript plumbing and verifier-side replay for this one Keccak segment.
4. Integrate boundary openings for committed input/output states.
5. Add an end-to-end benchmark against the generic `binius-examples` Keccak circuit path, while keeping the current microbenches as regression tripwires.
