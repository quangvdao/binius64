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

After that, wire one segment into the Binius transcript and sumcheck machinery, then benchmark against `binius_circuits` Keccak.
