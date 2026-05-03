# `keccak-prove`: design (v0)

**Status:** initial scaffold.
**Target:** standalone, batched Keccak-f[1600] proof sub-protocol on top of `binius64` IOP primitives.
**Out of scope:** sponge framing, padding, suffix injection, signature-scheme glue.

## What this is, and is not

This crate is a clean re-implementation of the 1-sumcheck-per-round Keccak protocol that the autoresearch worktree `binius64-keccak-prove-chi-iota-stable` converged on.
It is a rebuild, not a new algorithm.
The protocol shape is the one already documented at `docs/keccak-gkr-design.md` of the relevant branch.

The motivation is implementation hygiene, not algorithmic improvement.
The autoresearch worktree converged on a good design but accumulated scaffolding for variants that did not pay off.
We rebuild without that scaffolding, supervised, one PR at a time.

What is dropped:

- two-sumcheck-per-round prototype,
- `oblong_round` and packed-suffix variants,
- explicit pre-chi point-opening layer,
- generic Jolt-style virtual-claim DAG.

What is kept and rebuilt cleanly:

- single fused chi+iota+(theta+rho+pi) MLE-check per Keccak round,
- rotation predicate and its eq-vector cyclic-shift implementation,
- round-chained transparent-kernel hand-off across the 24 rounds.

## Reference: one round of Keccak-f[1600], fully unrolled

Before any sumcheck machinery, this is the round we are proving.
Every lane is an explicit `u64`, no loops, no array indexing into the state during the round body.
The five FIPS-202 steps theta, rho, pi, chi, iota appear in order.
The caller supplies the round constant; the full permutation is 24 calls of this function with `RC[0..24]`.

Conventions (FIPS 202, also in `~/Documents/Research/lattice-sig-aggregation/keccak-description.md`):

- lane index = `x + 5*y`,
- bit `z = 0` is the LSB of the `u64`,
- rho uses left rotations by the FIPS offsets,
- pi sends source `(x, y)` to destination `(y, (2x + 3y) mod 5)`,
- chi uses row neighbors `(x+1, y)` and `(x+2, y)`,
- iota only modifies lane `(0, 0)`.

```rust
/// One round of Keccak-f[1600], fully unrolled on a 25-lane u64 state.
#[inline(always)]
pub fn keccak_round(state: [u64; 25], rc: u64) -> [u64; 25] {
    // Unpack 25 input lanes. a_xy is lane (x, y).
    let a00 = state[ 0]; let a10 = state[ 1]; let a20 = state[ 2]; let a30 = state[ 3]; let a40 = state[ 4];
    let a01 = state[ 5]; let a11 = state[ 6]; let a21 = state[ 7]; let a31 = state[ 8]; let a41 = state[ 9];
    let a02 = state[10]; let a12 = state[11]; let a22 = state[12]; let a32 = state[13]; let a42 = state[14];
    let a03 = state[15]; let a13 = state[16]; let a23 = state[17]; let a33 = state[18]; let a43 = state[19];
    let a04 = state[20]; let a14 = state[21]; let a24 = state[22]; let a34 = state[23]; let a44 = state[24];

    // ---- theta, part 1: column parities. C[x] = XOR_y A[x, y].
    let c0 = a00 ^ a01 ^ a02 ^ a03 ^ a04;
    let c1 = a10 ^ a11 ^ a12 ^ a13 ^ a14;
    let c2 = a20 ^ a21 ^ a22 ^ a23 ^ a24;
    let c3 = a30 ^ a31 ^ a32 ^ a33 ^ a34;
    let c4 = a40 ^ a41 ^ a42 ^ a43 ^ a44;

    // ---- theta, part 2: column corrections. D[x] = C[x-1] XOR rotl(C[x+1], 1).
    let d0 = c4 ^ c1.rotate_left(1);
    let d1 = c0 ^ c2.rotate_left(1);
    let d2 = c1 ^ c3.rotate_left(1);
    let d3 = c2 ^ c4.rotate_left(1);
    let d4 = c3 ^ c0.rotate_left(1);

    // ---- theta, part 3: apply correction. B[x, y] = A[x, y] XOR D[x].
    let b00 = a00 ^ d0; let b10 = a10 ^ d1; let b20 = a20 ^ d2; let b30 = a30 ^ d3; let b40 = a40 ^ d4;
    let b01 = a01 ^ d0; let b11 = a11 ^ d1; let b21 = a21 ^ d2; let b31 = a31 ^ d3; let b41 = a41 ^ d4;
    let b02 = a02 ^ d0; let b12 = a12 ^ d1; let b22 = a22 ^ d2; let b32 = a32 ^ d3; let b42 = a42 ^ d4;
    let b03 = a03 ^ d0; let b13 = a13 ^ d1; let b23 = a23 ^ d2; let b33 = a33 ^ d3; let b43 = a43 ^ d4;
    let b04 = a04 ^ d0; let b14 = a14 ^ d1; let b24 = a24 ^ d2; let b34 = a34 ^ d3; let b44 = a44 ^ d4;

    // ---- rho + pi, fused. P[y, (2x + 3y) mod 5] = rotl(B[x, y], r[x, y]).
    // The destination lane on the LHS receives the rotated source lane.
    // From source (x, y):                    destination     rho offset
    let p00 = b00;                         // (0, 0) -> (0, 0)        0
    let p13 = b01.rotate_left(36);         // (0, 1) -> (1, 3)       36
    let p21 = b02.rotate_left( 3);         // (0, 2) -> (2, 1)        3
    let p34 = b03.rotate_left(41);         // (0, 3) -> (3, 4)       41
    let p42 = b04.rotate_left(18);         // (0, 4) -> (4, 2)       18

    let p02 = b10.rotate_left( 1);         // (1, 0) -> (0, 2)        1
    let p10 = b11.rotate_left(44);         // (1, 1) -> (1, 0)       44
    let p23 = b12.rotate_left(10);         // (1, 2) -> (2, 3)       10
    let p31 = b13.rotate_left(45);         // (1, 3) -> (3, 1)       45
    let p44 = b14.rotate_left( 2);         // (1, 4) -> (4, 4)        2

    let p04 = b20.rotate_left(62);         // (2, 0) -> (0, 4)       62
    let p12 = b21.rotate_left( 6);         // (2, 1) -> (1, 2)        6
    let p20 = b22.rotate_left(43);         // (2, 2) -> (2, 0)       43
    let p33 = b23.rotate_left(15);         // (2, 3) -> (3, 3)       15
    let p41 = b24.rotate_left(61);         // (2, 4) -> (4, 1)       61

    let p01 = b30.rotate_left(28);         // (3, 0) -> (0, 1)       28
    let p14 = b31.rotate_left(55);         // (3, 1) -> (1, 4)       55
    let p22 = b32.rotate_left(25);         // (3, 2) -> (2, 2)       25
    let p30 = b33.rotate_left(21);         // (3, 3) -> (3, 0)       21
    let p43 = b34.rotate_left(56);         // (3, 4) -> (4, 3)       56

    let p03 = b40.rotate_left(27);         // (4, 0) -> (0, 3)       27
    let p11 = b41.rotate_left(20);         // (4, 1) -> (1, 1)       20
    let p24 = b42.rotate_left(39);         // (4, 2) -> (2, 4)       39
    let p32 = b43.rotate_left( 8);         // (4, 3) -> (3, 2)        8
    let p40 = b44.rotate_left(14);         // (4, 4) -> (4, 0)       14

    // ---- chi, lane by lane. E[x, y] = P[x, y] XOR ((NOT P[x+1, y]) AND P[x+2, y]).
    // The only step with an AND. Quadratic in P over GF(2), one product per lane.
    let e00 = p00 ^ (!p10 & p20);
    let e10 = p10 ^ (!p20 & p30);
    let e20 = p20 ^ (!p30 & p40);
    let e30 = p30 ^ (!p40 & p00);
    let e40 = p40 ^ (!p00 & p10);

    let e01 = p01 ^ (!p11 & p21);
    let e11 = p11 ^ (!p21 & p31);
    let e21 = p21 ^ (!p31 & p41);
    let e31 = p31 ^ (!p41 & p01);
    let e41 = p41 ^ (!p01 & p11);

    let e02 = p02 ^ (!p12 & p22);
    let e12 = p12 ^ (!p22 & p32);
    let e22 = p22 ^ (!p32 & p42);
    let e32 = p32 ^ (!p42 & p02);
    let e42 = p42 ^ (!p02 & p12);

    let e03 = p03 ^ (!p13 & p23);
    let e13 = p13 ^ (!p23 & p33);
    let e23 = p23 ^ (!p33 & p43);
    let e33 = p33 ^ (!p43 & p03);
    let e43 = p43 ^ (!p03 & p13);

    let e04 = p04 ^ (!p14 & p24);
    let e14 = p14 ^ (!p24 & p34);
    let e24 = p24 ^ (!p34 & p44);
    let e34 = p34 ^ (!p44 & p04);
    let e44 = p44 ^ (!p04 & p14);

    // ---- iota. Only lane (0, 0) changes; every other output lane equals E.
    let o00 = e00 ^ rc;

    [
        o00, e10, e20, e30, e40,
        e01, e11, e21, e31, e41,
        e02, e12, e22, e32, e42,
        e03, e13, e23, e33, e43,
        e04, e14, e24, e34, e44,
    ]
}
```

Op count for one round on `u64` lanes:

| Step                          | XOR | rotl | NOT | AND |
|-------------------------------|----:|-----:|----:|----:|
| theta column parities `C[x]`  |  20 |    0 |   0 |   0 |
| theta corrections `D[x]`      |   5 |    5 |   0 |   0 |
| theta apply `B[x,y]`          |  25 |    0 |   0 |   0 |
| rho + pi fused `P[x,y]`       |   0 |   24 |   0 |   0 |
| chi `E[x,y]`                  |  25 |    0 |  25 |  25 |
| iota (lane (0,0) only)        |   1 |    0 |   0 |   0 |
| **per round**                 |  **76** | **29** | **25** | **25** |

Lane `(0,0)` has rho offset `0`, so its rotate is a no-op and counted as 0 above.

What this tells the protocol:

- Theta + rho + pi is a fixed sparse linear map from the 25 input lanes to the 25 pre-chi lanes `P[x,y]`.
  No new witness; in the sumcheck, it becomes 25 transparent kernels carrying rotation predicates.
- Chi contributes exactly **one quadratic term** per output lane, `P[x+1,y] · P[x+2,y]`.
  Total per round: 25 quadratic monomials over the pre-chi state, plus 50 linear terms.
  This is the only source of round-polynomial degree in the sumcheck integrand.
- Iota contributes one transparent additive offset on lane `(0,0)`.
  No witness dependence.

The integrand of the per-round MLE-check is exactly the chi+iota formula in the unrolled code, with each `P[x,y]` substituted by its linear expansion in the input lanes:

$$O^{(t)}_{i,j}(k) = P^{(t)}_{i,j}(k) \;+\; P^{(t)}_{i+2,j}(k) \;+\; P^{(t)}_{i+1,j}(k)\cdot P^{(t)}_{i+2,j}(k) \;+\; \delta_{i0}\delta_{j0}\,\mathrm{RC}_t(k).$$

The carried claim is $\sum_{i,j} W^{(t)}_{i,j}(k)\, O^{(t)}_{i,j}(k)$.
After substitution and folding the eq factor into the round polynomial, the sumcheck sees a degree-2 polynomial in the input-lane values, on a hypercube of size $2^{6 + \log h}$ (= $64h$ bits).

## Dependencies

`binius64` is depended on as a git source pinned to `binius-zk/binius64@115108d5` (public main).
No symlinks to local sibling clones.
No nested workspaces.

The minimal dependency surface is `binius_field`, `binius_math`, `binius_ip`, `binius_ip_prover`, `binius_transcript`.
We do not depend on `binius_iop_prover`, `binius_prover`, `binius_verifier`, `binius_circuits`, or `binius_keccak_check` (which the upstream main does not have).
`binius_iop` is added when the boundary-opening shim is wired up.

## Witness layout

For a batch of $h$ keccak-f instances, with $\ell = \log_2 h$ and $n = 6 + \ell$:

- 25 lane MLEs $A^{(t)}_{x,y}$ in $n$ variables, for $t = 0, \ldots, 24$.
- The 6 low variables index the 64 bit positions of a lane; the $\ell$ high variables index the instance.
- $A^{(0)}$ (input) and $A^{(24)} = O^{(23)}$ (output) are committed to the parent IOP layer.
- $A^{(1)}, \ldots, A^{(23)}$ are virtual.

## Boundary claim

At each end the parent IOP supplies a per-lane transparent linear relation $\langle A_{x,y},\, K_{x,y}\rangle$ where $K_{x,y}$ is verifier-computable.
The usual case is $K_{x,y}(b) = \beta_{x,y}\,\operatorname{eq}(\alpha, b)$ for verifier-chosen $\alpha \in \mathbb{F}^n$ and $\beta \in \mathbb{F}^{25}$.
The protocol carries this 25-tuple of kernels through the round chain.

## One sumcheck per round

For each Keccak round $t = 23, 22, \ldots, 0$:

1. Carry: 25 lanewise transparent kernels $W^{(t)}_{i,j}$ on the **output** of round $t$.
2. One fused MLE-check sumcheck of degree 2, with $\ell$ rounds.
   The integrand is the chi+iota formula

   $$O^{(t)}_{i,j}(k) = P^{(t)}_{i,j}(k) + P^{(t)}_{i+2,j}(k) + P^{(t)}_{i+1,j}(k)\cdot P^{(t)}_{i+2,j}(k) + \delta_{i0}\delta_{j0}\,\mathrm{RC}_t(k),$$

   with each pre-chi lane $P^{(t)}_{i,j}$ inlined as the explicit linear combination of rotated input lanes (theta+rho+pi composed and pushed through theta's rotation by 1).
   The eq factor against the boundary point is folded in by `MleToSumCheckDecorator`.
3. The bit-index axis is collapsed at the boundary point, not as separate sumcheck rounds.
   The first MLE-check round on the instance axis runs over u64 lane data so that the bit-index Lagrange pack is amortized into the heaviest fold; the remaining $\ell - 1$ rounds run on $25 \times 64$ extension-field blocks.
4. Hand-off: 25 carried kernels $K^{(t)}_{u,v}$ on the **input** lanes $A^{(t)} = O^{(t-1)}$, computed from the round's challenge point and the fixed rotation offsets.
   These become $W^{(t-1)}$ for round $t-1$.

After 24 rounds the kernels on $A^{(0)}$ close against the committed input via the parent IOP's `OracleLinearRelation`.
The kernels on $A^{(24)}$ close the same way against the committed output.

## Rotation predicate

Identical to `keccak-gkr-design.md:46-127`:

$$\operatorname{rot}_k(\alpha, b) = \widetilde{\mathbf{1}[\langle \alpha\rangle \equiv \langle b\rangle + k \pmod{64}]}, \qquad \alpha, b \in \mathbb{F}^6.$$

At sumcheck time the 64-entry vector $\{\operatorname{rot}_k(\alpha, b)\}_{b\in\{0,1\}^6}$ is the cyclic shift of the standard eq vector by $k$ positions, costing $O(64)$ per lane.

Because the bit-index axis is collapsed into the boundary point before sumcheck, the rotation predicate is consumed at boundary-claim construction time.
It does not appear inside any sumcheck round polynomial.

## Field choice

Sumcheck challenges live in $\mathrm{GF}(2^{128})$ via `binius_field::OptimalB128`.
The witness lives in $\mathrm{GF}(2)$.
No intermediate $\mathrm{GF}(2^k)$ packing for $1 < k < 128$ is used in v0; the Hashcaster-style Frobenius/dual-basis construction is explicitly out of scope.

## Univariate-skip parameter, and why it is fixed for v0

The bit-index Lagrange pack is, structurally, univariate skip with $m = 64$ on a 6-dimensional binary subspace of $\mathrm{GF}(2^{128})$.
That is fixed.
The independent knob is whether to also apply univariate skip to the first $c$ instance rounds.
For v0 the choice is $c_{\text{inst}} = 1$, matching the prior worktree's "first MLE-check round on native u64".

Going to $c_{\text{inst}} > 1$ trades transcript and verifier-side interpolation for prover-side bit-index amortization.
The transcript blows up from $2$ field elements per round to $\sim 3 \cdot 2^{c_{\text{inst}}} - 2$ field elements for the fused round.
This is a benchmark question, not a design question, and is left for the experiments section below.

### The bit-axis binary subspace, made precise

In the language of `~/Documents/Research/speeding-up-sumcheck/sections/7_univariate_skip.tex`, the bit-axis pack chooses:

- $m = 64$ components to pack (the 64 bit positions inside one 64-bit lane);
- an interpolation domain $D_0 \subset \mathrm{GF}(2^{128})$ of size 64.

We take $D_0$ to be the $\mathbb{F}_2$-linear subspace

$$D_0 \;=\; \mathrm{span}_{\mathbb{F}_2}\{1, \beta, \beta^2, \ldots, \beta^5\} \;\subset\; \mathrm{GF}(2^{128}),$$

where $\beta$ is a generator of an embedded `AESTowerField8b` inside `B128`.
Concretely, this is what `binius_math::BinarySubspace::<B8>::with_dim(6).isomorphic::<B128>()` returns.
Bit position $j \in \{0, \ldots, 63\}$ identifies with the element $u_j \in D_0$ in standard order (`nat2(j_5, \ldots, j_0)`).

Why this specific subspace, instead of any 6-dim $\mathbb{F}_2$-subspace of $\mathrm{GF}(2^{128})$:

1. Any $\mathbb{F}_2$-linear $D_0$ would work for the bit-axis Lagrange pack itself, because the pack is just one inner product $\sum_j L_j(\alpha_{\text{bit}}) \cdot \text{bit}_j$ for the single boundary point $\alpha_{\text{bit}} \in \mathrm{GF}(2^{128})$.
2. Picking the AES-tower $\mathbb{F}_2$-basis is what lets us **reuse binius64's NTT-lookup machinery** unchanged.
   The byte tables in `crates/prover/src/and_reduction/{ntt_lookup.rs, fold_lookup.rs}` are tied to that specific basis.
3. The same basis is also what gets used downstream by ring switching and BaseFold in the underlying binius64 IOP layer, so keeping one canonical basis avoids extra change-of-basis work at the layer boundary.

The bit-axis multilinear restricted to one lane row is then a function $D_0 \to \mathbb{F}_2$, and its multilinear extension on the 6 bit-axis variables coincides with the unique univariate $f \in \mathrm{GF}(2^{128})[Z]$ of degree at most 63 with $f(u_j) = \text{bit } j$.
Lagrange evaluation of $f$ at $\alpha_{\text{bit}}$ is the bit-axis Lagrange pack.

### How the NTT lookup implements the pack

For one fixed boundary point $\alpha_{\text{bit}}$, the bit-axis pack is the linear map

$$\mathbb{F}_2^{64} \longrightarrow \mathrm{GF}(2^{128}), \qquad
(c_0, \ldots, c_{63}) \;\longmapsto\; \sum_{j=0}^{63} L_j(\alpha_{\text{bit}})\,c_j.$$

Two structural facts make this cheap:

1. The map is $\mathbb{F}_2$-linear in the input. So the 64-bit input vector splits into 8 independent **bytes** and the partial outputs add:

   $$\mathrm{pack}(c) \;=\; \sum_{b=0}^{7} \mathrm{pack}\big(\text{byte } b \text{ of } c\big).$$

2. Each per-byte map has only $256$ inputs, so it can be precomputed once per boundary point as a 256-entry table of `B128` scalars.

This is exactly the **Method of Four Russians** instantiation that binius64 uses for its `FoldLookup` and (in a slightly larger output-domain form) for `ntt_lookup`.
For our setting:

- **Setup, once per boundary point:** compute $L_j(\alpha_{\text{bit}})$ for $j \in \{0, \ldots, 63\}$ by `lagrange_evals_scalars` on $D_0$ (one O(64) pass), then build 8 tables of 256 `B128` entries each; total `8 * 256 * 16 B = 32 KiB`.
- **Per word:** 8 byte-indexed table loads + 7 `B128` adds + one outer `B128` multiply by the carried-kernel scalar (the boundary kernel $K_{x,y}(b) = \beta_{x,y}\,\mathrm{eq}(\alpha, b)$ contributes that scalar after the bit-axis is collapsed).

This is the same cost shape as binius64's `FoldLookup`, just with $D_0$ instead of binius64's larger $D$ (size 128 in their case, because Phase 1 needs a degree-2 round message; here we only need a single point evaluation).
In code, the implementation can lift a slimmed copy of `binius64::and_reduction::fold_lookup::FoldLookup` (which already does the byte-table construction over a `BinarySubspace`).
The relevant existing helpers are:

- `binius_math::BinarySubspace::with_dim(6)` to construct $D_0$,
- `binius_math::univariate::lagrange_evals_scalars(&D_0, alpha_bit)` to get the 64 weights,
- the byte-table construction in `crates/prover/src/and_reduction/fold_lookup.rs:38-89`.

### Relationship to binius64's BitAnd Phase 1

The same construction sits at the heart of binius64's BitAnd zerocheck (`hashcaster-vs-binius64.md` §3.2, §3.9), with one upgrade: BitAnd Phase 1 evaluates a degree-2 round message $A(Z) B(Z) - C(Z)$ on the bit axis, so it needs an output domain of size $2 \cdot 64 = 128$ rather than 64.
binius64 picks the doubled domain $D = \mathrm{span}_{\mathbb{F}_2}\{1, \beta, \ldots, \beta^6\}$ with $D_0$ as its lower half, and the prover sends only the upper-half evaluations because $R_0$ vanishes on $D_0$ where the AND constraints live.

`keccak-prove` does not run a sumcheck round on the bit axis: the bit axis is collapsed at boundary-claim construction time via Lagrange evaluation at one point.
So we stop at $D_0$ and the simpler 64-entry output suffices.
This is the pure "Lagrange pack" view of univariate skip, restricted to evaluating the packed univariate at the verifier's challenge directly, rather than to running a univariate sumcheck round followed by a fresh challenge.

## Open knobs to settle by benchmark, not by argument

- `c_inst` on the instance axis. Initial value 1.
- 1-sumcheck-fused vs 2-sumcheck-per-round. Initial choice 1-sumcheck-fused. The 2-sumcheck variant has lighter per-cube-point cost but worse transcript and inversion counts. Worth re-running once both exist.
- Whether to commit at every $R$ rounds. Initial $R = 24$ (Option A, only input and output committed). Option B (flattened trace) is not implemented in v0.
- Switch from $\mathrm{GF}(2^{128})$ challenges to $\mathrm{GF}(2^{64})$. Saves one factor of 2 per multiplication if soundness budget allows. Out of v0.

## Implementation plan, one PR at a time

Each step lands as a self-contained commit with a passing test against a reference computation.
No step writes code that depends on a later step.

1. **Scaffold.** `Cargo.toml`, `rust-toolchain.toml`, `lib.rs`, this design doc. *(this commit)*
2. **Trace.** `[u64; 25]` input → 25 lane MLEs at all 25 round boundaries, no extension-field expansion. Adapt the relevant logic from the autoresearch worktree's `crates/keccak-check/src/trace.rs`, do not bring over the whole module.
3. **Boundary claim.** A `MixedClaim` type, evaluation at a random $\alpha \in \mathbb{F}^n$ with the bit-index axis collapsed via Lagrange pack. Round-trip test against direct lane evaluation.
4. **Rotation predicate.** Cyclic-shift eq vector, helpers for the 25 per-lane carried kernels. Round-trip test against the explicit indicator.
5. **One fused round.** chi+iota with the linear pre-chi map inlined. Sumcheck over $\ell$ instance variables, first round specialized to u64. Test against a reference table-evaluation of the same round map.
6. **Round chain.** 24 rounds end-to-end, virtual-claim hand-off between rounds. Test against `keccak_f1600` from any standard implementation.
7. **Boundary opening shim.** Translate the carried lane kernels at $t=0, t=24$ to `OracleLinearRelation` claims against the parent IOP.
8. **Bench.** Prover wall-clock against the upstream `binius64::circuits::keccak` path and against the autoresearch worktree, at $h \in \{2^{10}, 2^{14}, 2^{17}\}$.
9. **Tune.** Once all of the above passes, run the experiments listed in the previous section.

Each step is reviewed before the next is started.
No AI-generated code lands without a human-supervised diff.

## References

- `~/Documents/SNARKs/binius64-keccak-prove-chi-iota-stable/docs/keccak-gkr-design.md` — prior protocol-level design note.
- `~/Documents/SNARKs/binius64-keccak-prove-chi-iota-stable/crates/keccak-check/` — prior implementation. Kept as legacy reference. Not depended on.
- `~/Documents/Research/speeding-up-sumcheck/sections/7_univariate_skip.tex` — univariate skip definition.
- `~/Documents/Research/lattice-sig-aggregation/keccak-description.md` — FIPS-202 reference for this repo.
- KeccakCheck (eprint 2025/1764) — original prime-field protocol.
- Hashcaster (`https://hackmd.io/@levs57/SJ4fuZMD0`, `https://github.com/morgana-proofs/hashcaster`) — Frobenius approach explicitly not used here.
