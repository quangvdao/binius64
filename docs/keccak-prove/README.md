# Keccak Proving Notes

These notes track the Keccak-specific proving work for Binius64.

- `full-protocol.md` is the current working protocol note.
- `design.md` preserves the earlier fused-round design.
- `performance.md` records benchmark configuration, machine comparisons, and scaling cliffs.

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
- verifier-field row challenge support for the transcripted first-round message, so the row weights are derived from the same MLE-check point used by the verifier;
- committed `A`/`D` witness layout helpers for the locked 32-word block layout;
- committed witness construction from native Keccak traces, including materialized theta correction words and zero padding;
- Shift-compatible operand helpers for virtual `B` references and `D` correctness relations.
- production-Shift constraint-system schemas for chi operand pushback and `D` correctness:
  - chi rows lower `P/Q/C` witness terms to committed `A` and `D` words through virtual `B`;
  - `D` rows are encoded as degenerate AND constraints `0 * 0 = D_correctness_operand`;
  - tests run the production Shift prover/verifier on both Keccak schemas.
- a full v0 production-path constraint system:
  - public constants hold the all-one word and 24 iota round constants;
  - committed private witness holds the locked `A`/`D` layout;
  - full chi rows prove `(1 + B[x+1,y]) & B[x+2,y] = B[x,y] + A_next[x,y] + iota`;
  - production `Prover`/`Verifier` tests cover BitAnd, Shift, ring-switching, and PCS end to end.
- a tensor-shaped v0 row layout with 32 rows per `(permutation, round)`:
  - 25 chi rows;
  - 5 `D` correctness rows;
  - 2 padding rows.
- a first structured verifier prototype for v0 Shift monster evaluation, checked against the generic Shift verifier evaluator.

Performance lessons, benchmark commands, and measured checkpoints are collected in
`performance.md`.

## Current Proof Chain

The current v0 chain is committed-round-boundary, not the older serial transparent GKR chain.

1. Commit the locked `A`/`D` witness:
   - every round-boundary `A_r[x,y]`;
   - every theta correction `D_r[x]`;
   - no committed pre-chi `B`.
2. Run the chi/iota Spartan outer segment in the production BitAnd shape:

   ```text
   P * Q - C = 0
   ```

   where

   ```text
   P = 1 + B[x+1,y]
   Q =     B[x+2,y]
   C =     B[x,y] + A_next[x,y] + iota
   ```

3. In the custom transcripted segment, convert final outer evaluations into witness-only Shift claims:
   - `P_witness = B[x+1,y]`;
   - `Q_witness = B[x+2,y]`;
   - `C_witness = B[x,y] + A_next[x,y]`;
   - the all-one term in `P` is a transparent active-row selector correction;
   - iota is a transparent correction to the `C` claim.
4. Use production Shift to reduce the chi operand claims to committed `A`/`D` witness evaluations. Each virtual pre-chi lane is lowered as:

   ```text
   B_r[pi(x,y)] = rotl_{rho[x,y]}(A_r[x,y] + D_r[x])
   ```

5. Use production Shift again for `D` correctness, encoded as degenerate AND rows:

   ```text
   0 * 0 = D_r[x] + sum_y A_r[x-1,y] + sum_y rotl_1(A_r[x+1,y])
   ```

6. The remaining integration work is to batch these Shift output claims with the boundary claims on `A_0` and `A_24`, then discharge them through the same ring-switching and PCS opening path used after production Shift.

The older serial GKR-style transparent kernel pushback remains useful as a mathematical guardrail, but it is not the hot implementation path now that `A` and `D` are committed and all round rows can be handled in parallel.

The implemented end-to-end v0 path currently takes the most production-faithful version of this
chain: build a normal Binius64 `ConstraintSystem` containing full chi rows plus `D` rows, then let
the production `Prover` run BitAnd, Shift, ring-switching, and PCS. That path is in
`crates/keccak-prove/src/v0.rs`.

## Locked-in committed witness layout

The committed witness layout for the specialized Keccak prover is:

```text
per permutation, per round-boundary block:
  slots 0..24   A_r[x,y]   round-boundary state lanes
  slots 25..29  D_r[x]     theta correction words
  slots 30..31  padding
```

Each block is exactly 32 committed words. A full Keccak-f[1600] permutation uses 25 blocks:

```text
blocks r = 0..24
A_r exists in every block
D_r exists only for r = 0..23
block 24 stores the final A_24 state and padding only
```

The stable word indices are:

```text
block(p, r) = p * 25 * 32 + r * 32
A(p, r, lane) = block(p, r) + lane
D(p, r, x)    = block(p, r) + 25 + x
```

Active data per permutation:

```text
A: 25 states * 25 lanes = 625 words
D: 24 rounds * 5 words  = 120 words
active total            = 745 words
padding                 = 55 words
committed total         = 800 words
```

This intentionally commits the theta correction words `D` but not the pre-chi state `B`.
The virtual pre-chi lanes are:

```text
B_r[y, 2x + 3y] = rotl_{rho[x,y]}(A_r[x,y] + D_r[x])
```

Equivalently, every `B` reference is lowered into two shifted committed terms:

```text
B_r[pi(x,y)] = shifted(A_r[x,y]) + shifted(D_r[x])
```

The `D` correctness relations are linear and Shift-compatible:

```text
D_r[x] + sum_y A_r[x-1,y] + sum_y rotl_1(A_r[x+1,y]) = 0
```

These can be encoded as degenerate AND constraints with `0 * 0 = linear_operand`, allowing the
production Shift reduction to reduce them to committed witness openings. Keeping `B` virtual avoids
an additional 25 committed words per round while preserving small enough shifted operands for the
chi/iota claims.

See `performance.md` for benchmark commands, measured checkpoints, and scale results.

## Next steps

The next implementation milestone is verifier/prover specialization on top of the full v0 path:

1. Expose a production Shift verification hook so v0 can replace generic monster evaluation with the structured tensor evaluator.
2. Continue optimizing the structured evaluator so it wins at 128 permutations, not only at larger batches.
3. Add a v0 setup/key-building benchmark, because generic `KeyCollection` construction is likely avoidable from the 32-slot tensor schema.
4. Scale the full-path comparison beyond 128 permutations and keep matching against `binius-examples` Keccak at equivalent rate-permutation counts.
