# binius-keccak-prove

Implementation home for the Binius-native Keccak-f[1600] proving path.

The v0 target is a performance-oriented prover that adapts production Binius64 machinery rather than starting from a serial toy protocol:

- BitAnd-style $K=3$ outer relation for chi;
- bit-axis additive-NTT lookup from the start;
- Keccak-specialized replacement for generic shift reduction;
- parallel or segmented-parallel proving across batches and rounds;
- early benchmarks against the existing generic Keccak circuit path.

Protocol notes live in `../../docs/keccak-prove/`. Benchmark setup and cross-machine performance
observations are tracked in `../../docs/keccak-prove/performance.md`.

The first implementation milestone is a residual-zero NTT test for one Keccak round segment:

```text
cargo test -p binius-keccak-prove keccak_bitand_residual_matches_native_round
```

Run the full crate sanity suite with:

```text
cargo test -p binius-keccak-prove
```

Benchmark commands, throughput conventions, and current results are tracked in
`../../docs/keccak-prove/performance.md`.

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

See `../../docs/keccak-prove/performance.md` for the benchmark matrix, production comparison
commands, scale sweeps, and cross-machine notes.
