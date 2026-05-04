# Repeated Constraint-System Batching

This branch prototypes proving many copies of the same Binius64 constraint system with different
witness values. The verifier binds a compact repeated-circuit descriptor and uses the repeated Shift
monster structure instead of treating the Shift verifier as a fully flat circuit.

## Current Descriptor Flow

The current repeated descriptor is `RepeatedConstraintSystem`.

1. Build and prepare one base `ConstraintSystem`.
2. Wrap it with `RepeatedConstraintSystem::new(base, log_instances)`.
3. Build one base-layout `ValueVec` per instance.
4. Flatten those instance witnesses with `RepeatedConstraintSystem::to_flat_value_vec`.
5. Set up with `Verifier::setup_repeated`.
6. Prove with `Prover::setup_repeated` and `prove_repeated`.
7. Verify with `verify_repeated`.

The v0 value layout is instance-major for non-constants and shared for base constants. Public input
handling is intentionally minimal: only the first flat public prefix remains public. The Keccak
example below uses no public `inout` values at all, so its message words, digest words, and internal
Keccak state are hidden witness data.

## Keccak Example

The integration test in `crates/prover/tests/repeated_keccak.rs` uses the frontend fixed-length
Keccak-256 gadget as the base circuit. It then repeats that same prepared circuit over distinct
hidden witness instances.

Fast correctness smoke test:

```sh
cargo test -p binius-prover --test repeated_keccak -- --nocapture
```

End-to-end repeated Keccak proof:

```sh
cargo test --release -p binius-prover --test repeated_keccak repeated_keccak_prove_verify -- --ignored --nocapture
```

Verifier timing harness:

```sh
cargo test --release -p binius-prover --test repeated_keccak repeated_keccak_verifier_print_runtimes -- --ignored --nocapture
```

Prover key-materialization timing harness:

```sh
cargo test --release -p binius-prover --test repeated_keccak repeated_keccak_prover_key_materialization_print_runtimes -- --ignored --nocapture
```

## What The Keccak Numbers Mean

The verifier harness compares ordinary flat verification against repeated verification over the same
hidden repeated Keccak witnesses. On small instance counts, repeated verification can lose to fixed
overhead and timing noise. In the local run that introduced the harness, the repeated verifier became
faster by 16 Keccak instances.

The prover-key harness compares a flat prover that materializes Shift keys for every flat instance
against a repeated prover that materializes Shift keys for the base circuit only. In the local run
that introduced the harness, `flat_key_words` grew from `1024` to `16384` between one and sixteen
Keccak instances, while `repeated_key_words` stayed at `1024`.

Total proving time does not yet scale like the key counts because the witness oracle, PCS, and
non-Shift work still operate on the flat repeated witness. The current implementation demonstrates
the repeated Shift path and descriptor binding, not full tensorization of the entire prover.

## Remaining Caveats

- The hot verifier path trusts the descriptor bound into the transcript; structural expansion checks
  are setup-side guardrails.
- The repeated implementation currently has one public binding mode, `FlatPublicInputs`.
- Keccak here is fixed-length Keccak-256. It is a concrete repeated frontend-circuit example, not a
  general SHAKE/sponge circuit yet.
- Broader ML-DSA proving still needs the sponge/absorbing relation, bit-heavy ML-DSA relations, and
  the binary/lattice bridge over hidden shared values.
