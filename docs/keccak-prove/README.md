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
