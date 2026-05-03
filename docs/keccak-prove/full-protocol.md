# `keccak-prove`: full protocol working note

**Status:** Binius-native v0 planning note.
**Goal:** implement a parallel Keccak-f[1600] proving path inside `crates/keccak-prove`.
**Current design target:** adapt production BitAnd, shift-reduction, and bit-axis NTT lookup machinery to Keccak's fixed round structure.

This file is intentionally separate from `design.md`.
`design.md` records the earlier one-fused-sumcheck-per-round design.

Earlier sections of this note develop a conservative serial claim-propagation invariant.
That invariant remains valuable for tests, but it is not the implementation target anymore.
Now that this work lives inside the Binius64 tree, v0 should be more ambitious: parallel, NTT-backed, and closely derived from the production BitAnd and shift code.

**Current implementation update.** After locking in the committed `A`/`D` layout, the near-term v0
path uses production Shift directly. The older sections that describe a Keccak-aware transparent
kernel replacement for Shift should be read as a possible later specialization, not as the current
implementation plan. The current chain is:

1. commit all round-boundary `A_r[x,y]` words and theta correction `D_r[x]` words;
2. prove chi/iota with the BitAnd-shaped Spartan outer relation `P * Q - C = 0`;
3. convert the final outer `P/Q/C` evaluations into witness-only Shift claims by applying the
   transparent all-one and iota corrections;
4. run production Shift on the chi operand schema, where every virtual `B` reference is lowered to
   shifted committed `A` and `D` terms;
5. run production Shift on degenerate AND rows for `D` correctness;
6. batch the resulting witness-evaluation claims with the boundary claims and discharge them
   through the production opening path.

## 1. Keccak state notation

We use the same lane coordinates as `canvas/keccak-shake-permutation.canvas.tsx`.

- A Keccak-f[1600] state has lanes indexed by $(x,y) \in \{0,1,2,3,4\}^2$.
- The lane index is `x + 5*y`.
- Each lane is a 64-bit word.
- Bit $z = 0$ is the least significant bit of the word.
- A full bit coordinate is $(x,y,z)$ with $z \in \{0,\ldots,63\}$.

For a batch of $h = 2^\ell$ Keccak permutations, let $u \in \{0,1\}^\ell$ index the batch item.
The bit-level state at round boundary $t$ is

$$A^{(t)}_{x,y}(z,u) \in \mathbb{F}_2,$$

where $t = 0,\ldots,24$.
The initial state is $A^{(0)}$ and the final state is $A^{(24)}$.

## 2. Round-step notation

Within one round $t$, we use the following symbolic layers:

$$A^{(t)} \xrightarrow{\theta} B^{(t)} \xrightarrow{\rho+\pi} D^{(t)} \xrightarrow{\chi} E^{(t)} \xrightarrow{\iota} A^{(t+1)}.$$

This matches the canvas terminology where the nonlinear step is written

$$E[x,y] = D[x,y] \oplus ((\neg D[x+1,y]) \wedge D[x+2,y]).$$

All $x$ coordinates are taken modulo 5.
The $y$ coordinate is unchanged by chi.

Over $\mathbb{F}_2$, $\neg q = 1 + q$.
So the chi equation becomes

$$E^{(t)}_{x,y} = D^{(t)}_{x,y} + (1 + D^{(t)}_{x+1,y})D^{(t)}_{x+2,y}.$$

Iota is affine:

$$A^{(t+1)}_{0,0}(z,u) = E^{(t)}_{0,0}(z,u) + \mathrm{RC}_t(z),$$

and for $(x,y) \ne (0,0)$,

$$A^{(t+1)}_{x,y}(z,u) = E^{(t)}_{x,y}(z,u).$$

## 3. Candidate Spartan-outer chi relation

The intended $K=3$ outer relation should isolate the quadratic part of chi using three virtual operand columns.
For each round $t$ and output lane $(x,y)$, define the three Spartan operands

$$P^{(t)}_{x,y} = 1 + D^{(t)}_{x+1,y},$$

$$Q^{(t)}_{x,y} = D^{(t)}_{x+2,y},$$

$$R^{(t)}_{x,y} = D^{(t)}_{x,y}.$$

Then chi+iota is equivalent to

$$P^{(t)}_{x,y} \cdot Q^{(t)}_{x,y} - R^{(t)}_{x,y} = A^{(t+1)}_{x,y} + \delta_{x,0}\delta_{y,0}\mathrm{RC}_t.$$

Since the field has characteristic two, the minus sign is the same operation as plus.
We keep the minus sign because this is the Spartan/R1CS shape:

$$P \cdot Q - R = \mathrm{output\ claim}.$$

This is better than putting $A^{(t+1)}$ and iota inside $R$.
The output lane value $A^{(t+1)}_{x,y}$ should be treated as an **input claim** to this round proof, either supplied by the committed output boundary or carried from the next step of the chain.
The iota term is a transparent correction to that input claim.

The constant $1$ in $P = 1 + D_{x+1,y}$ can be folded into the virtual operand $P$.
If the constraint domain is padded, this constant should be a transparent validity selector rather than the all-one polynomial on padded rows.

We use $(P,Q,R)$ for the Spartan operands so that $A^{(t)}$ remains reserved for Keccak round-boundary state.

## 4. Constraint-domain indexing

A single chi constraint is indexed by

$$\kappa = (t,x,y,z,u).$$

The dimensions are:

- $t \in \{0,\ldots,23\}$ for the 24 rounds,
- $(x,y) \in \{0,\ldots,4\}^2$ for the 25 lanes,
- $z \in \{0,\ldots,63\}$ for the bit position,
- $u \in \{0,1\}^\ell$ for the batch item.

So the raw chi constraint count is

$$24 \cdot 25 \cdot 64 \cdot h.$$

We need decide how this index is represented as a hypercube.
Natural options:

1. Pad the $(t,x,y)$ axis to a power of two and run one global Spartan-outer over all rounds and lanes.
2. Run one Spartan-outer per round and include only $(x,y,z,u)$ in the constraint domain.
3. Run one Spartan-outer per row $y$ or per round-row pair, to keep the lane-neighbor structure smaller.

The design target is likely option 1 or 2, but this should be settled after the inner reduction notation is clear.

## 5. Outer batch weight

The Spartan-outer proves a batched claim of the form

$$S_{\mathrm{outer}} = S_{\mathrm{next}}.$$

The left-hand side is

$$S_{\mathrm{outer}} = \sum_{\kappa \in \mathcal{C}} \beta(\kappa)\left(P(\kappa)Q(\kappa)-R(\kappa)\right).$$

The right-hand side is

$$S_{\mathrm{next}} = \sum_{\kappa \in \mathcal{C}} \beta(\kappa)\left(A_{\mathrm{next}}(\kappa)+\mathrm{Iota}(\kappa)\right).$$

Here $\kappa$ is a constraint coordinate.
At minimum it contains $(t,x,y,z,u)$: round, lane, bit, and batch item.
The right-hand side is not a witness polynomial introduced by the outer proof.
It is an input claim on the next round-boundary state, plus a transparent iota correction.

The batch weight should factor by axis:

$$\beta(t,x,y,z,u) = \beta_{\mathrm{round}}(t)\beta_x(x)\beta_y(y)\beta_{\mathrm{bit}}(z)\beta_{\mathrm{batch}}(u).$$

For the bit axis, the univariate-skip/Lagrange-pack weight is

$$\beta_{\mathrm{bit}}(z) = L_z(\zeta),$$

where $\zeta$ is the verifier's bit-axis challenge and $z$ ranges over the size-64 binary subspace $D_0$.
This is the same bit-axis pack described in `design.md`, but now it weights the outer chi constraints.

For the batch axis,

$$\beta_{\mathrm{batch}}(u) = \mathrm{eq}(\rho,u).$$

For the round and lane axes we have two implementation choices:

1. **Padded Boolean hypercube.**
   Encode $t$, $x$, and $y$ in binary, pad to the next power of two, and multiply by transparent validity selectors.
   This matches existing Boolean-hypercube sumcheck infrastructure and is still the default for the serial test harness.
   The Binius-native v0 target in §13 may benchmark exact mixed domains after the padded path works.
2. **Exact mixed domains.**
   Use a size-24 round domain and size-5 lane-coordinate domains, or factor $24$ and $25$ as small domains such as $3 \times 8$ and $5 \times 5$.
   This may reduce padding and transcript waste, but it requires more custom domain code.

The first implementation should use padded Boolean domains unless benchmarks show the padding cost dominates.
The notation below should not assume exact mixed domains.

### Padded Boolean domains for v0

For v0, use padded Boolean encodings for the non-bit structural axes.
In the global form:

- the round index $t \in \{0,\ldots,23\}$ is encoded by $5$ Boolean variables, padded to $32$ values;
- the lane coordinates $x,y \in \{0,\ldots,4\}$ are each encoded by $3$ Boolean variables, padded to $8$ values each;
- the batch index $u \in \{0,1\}^{\ell}$ is already Boolean;
- the bit index $z$ is handled by the size-64 binary-subspace Lagrange pack, not by six ordinary sumcheck rounds.

Let

$$V_{\mathrm{round}}(t) = \mathbf{1}[t < 24], \qquad V_5(a) = \mathbf{1}[a < 5].$$

In the padded Boolean implementation, the actual batch weight is

$$\beta_{\mathrm{pad}}(t,x,y,z,u) = V_{\mathrm{round}}(t)V_5(x)V_5(y)\beta(t,x,y,z,u).$$

This means invalid padded rows contribute zero to both sides of the outer equation.
With this convention, the constant term in $P = 1 + D_{x+1,y}$ may be the all-one polynomial, because invalid rows are killed by $\beta_{\mathrm{pad}}$.
Equivalently, an implementation may fold the validity selector into that constant term; the algebraic claim is the same.

The serial per-round form drops the round variables and uses only $x,y,z,u$.
Then the padded weight is

$$\beta_{\mathrm{pad}}(x,y,z,u) = V_5(x)V_5(y)\beta_x(x)\beta_y(y)\beta_{\mathrm{bit}}(z)\beta_{\mathrm{batch}}(u).$$

This is the likely v0 path.
It keeps one round's proof small and avoids committing intermediate states.

## 6. Keccak-aware Spartan-inner target

The Spartan-outer produces evaluation claims about $P,Q,R$.
The inner reduction must reduce these claims to claims about the underlying round-boundary lane polynomials.

The three virtual operand columns are:

$$P^{(t)}_{x,y} = 1 + D^{(t)}_{x+1,y},$$

$$Q^{(t)}_{x,y} = D^{(t)}_{x+2,y},$$

$$R^{(t)}_{x,y} = D^{(t)}_{x,y}.$$

The next-state value $A^{(t+1)}_{x,y}$ is not reduced by this inner step as a witness operand.
It is the input claim on the right-hand side of the outer equation.
In a serial chain, this claim is the carried claim from the next round.
At a segment boundary, it is a committed boundary claim.

The inner reduction therefore only has to explain the three $D$-derived virtual columns in terms of $A^{(t)}$.
This is where we take binius64's valuable idea, but simplify it for Keccak.

### What binius64 does

binius64's Shift Reduction is a generic Spartan-inner for arbitrary shifted-word constraints.
It handles eight shift variants:

$$\mathrm{SLL}, \mathrm{SRL}, \mathrm{SRA}, \mathrm{ROTR}, \mathrm{SLL32}, \mathrm{SRL32}, \mathrm{SRA32}, \mathrm{ROTR32}.$$

The prover builds a `KeyCollection` keyed by witness word, operand index, shift variant, and shift amount.
Then it runs two bivariate-product sumchecks:

1. Phase 1 sums over bit position and shift amount $(j,s)$.
   It proves a claim of the form $\sum_{j,s,\mathrm{op}} g_{\mathrm{op}}(j,s)h_{\mathrm{op}}(j,s)$.
   The $h_{\mathrm{op}}$ factors are transparent shift indicators.
   The $g_{\mathrm{op}}$ factors are witness/constraint dependent and are built from the `KeyCollection`.
2. Phase 2 sums over witness-word index $y$.
   It proves $\sum_y W(r_j,y)M(r_s,y)$, where $W(r_j,y)$ is the witness folded along its bit axis and $M$ is the monster multilinear built from the key collection.

This is powerful because binius64 supports arbitrary shifted operands in arbitrary constraints.
It is also more general than Keccak needs.

### What Keccak needs

Keccak has no general shifts in this inner map.
It only has rotations by fixed FIPS offsets inside $\rho$, plus fixed lane permutations and XOR-linear theta.
Every linear operation is public, fixed, and reused for all batch items:

- theta: fixed column-parity map on the $5 \times 5$ lane grid;
- rho: fixed rotate-left offset per lane;
- pi: fixed lane permutation;
- chi-neighbor selection: fixed maps $(x,y) \mapsto (x+1,y)$, $(x+2,y)$, and $(x,y)$.

So the inner map is much more structured than binius64's generic shifted-word map.
There is no need for an eight-variant shift table, no `KeyCollection`, and no monster multilinear indexed by arbitrary witness-word adjacency.

The natural Keccak-aware inner is a transparent-kernel pushback:

1. Start from the batched claims on $P,Q,R$ at the outer challenge.
2. Convert them into claims on three $D$ lanes:
   - $P$ contributes to $D^{(t)}_{x+1,y}$ and a transparent constant term;
   - $Q$ contributes to $D^{(t)}_{x+2,y}$;
   - $R$ contributes to $D^{(t)}_{x,y}$.
3. Push those $D$-lane kernels backward through the transposes of $\pi$, $\rho$, and $\theta$.
   For $\pi$ and $\rho$, the transpose is the inverse lane permutation or inverse bit rotation.
   For $\theta$, it is the transpose of the fixed column-parity linear map, not a generic shift reduction.
4. Output claims on the underlying $A^{(t)}_{x,y}$ lane polynomials.

The prover-side computation should be small and structured:

- lane-neighbor selection is just permutation of 25 lane kernels;
- pi transpose is another permutation of 25 lane kernels;
- rho transpose is a cyclic shift of the 64-entry bit-axis kernel by the fixed FIPS offset;
- theta transpose is a fixed sparse XOR map among the 25 lane kernels, using column parities and one-bit rotations.

This is closer to the old transparent-kernel idea than to binius64's generic Shift Reduction.
The difference is that the outer relation is now Spartan-like with $K=3$, while the inner is still Keccak-aware and transparent.

### Claim and kernel objects

The object carried between layers is a linear claim with an explicit transparent kernel descriptor:

$$\mathsf{Claim}(A^{(t)}, K^{(t)}, s^{(t)}) \quad\text{means}\quad s^{(t)}=\langle K^{(t)}, A^{(t)}\rangle.$$

Expanded over lanes, bits, and batch items,

$$s^{(t)}=\sum_{x,y,z,u} K^{(t)}_{x,y}(z,u)A^{(t)}_{x,y}(z,u).$$

For the v0 serial protocol, $K^{(t)}$ is not an oracle and is not committed.
It is verifier-computable from transcript challenges and the fixed Keccak linear maps.
The prover only needs it to compute the next claimed scalar and to prepare the final boundary opening.

The Spartan-outer ends with claimed evaluations of the three virtual columns.
After batching them with fresh challenges $\lambda_P,\lambda_Q,\lambda_R$, the inner target is one claim

$$
\lambda_P P_\zeta(\alpha_x,\alpha_y,\alpha_u)
+ \lambda_Q Q_\zeta(\alpha_x,\alpha_y,\alpha_u)
+ \lambda_R R_\zeta(\alpha_x,\alpha_y,\alpha_u),
$$

where $\zeta$ is the bit-axis challenge and $(\alpha_x,\alpha_y,\alpha_u)$ is the remaining outer sumcheck point.
This expression induces a transparent kernel $K_D$ on $D^{(t)}$ plus a transparent constant from the $1$ in $P$.

Use a canonical padded definition of the virtual columns on Boolean rows:

$$
P_{x,y}=V_5(x)V_5(y)(1+D_{x+1,y}),\qquad
Q_{x,y}=V_5(x)V_5(y)D_{x+2,y},\qquad
R_{x,y}=V_5(x)V_5(y)D_{x,y}.
$$

The neighbor indices are only interpreted on valid rows, and invalid padded rows are zero.
With this convention, the $P$ evaluation contributes the transparent constant

$$\lambda_P\cdot (V_5(\alpha_x)V_5(\alpha_y)),$$

and the nonconstant part contributes to the shifted $D$ lane.
The induced $D$ kernel is:

$$
K_D(a,b,z,u)=L_z(\zeta)\operatorname{eq}(\alpha_u,u)
\cdot
\sum_{\star\in\{P,Q,R\}}\lambda_\star K^\star_D(a,b;\alpha_x,\alpha_y),
$$

where

$$
K^P_D(a,b;\alpha_x,\alpha_y)=
\sum_{x,y:\,a=x+1,\,b=y}
\operatorname{eq}(\alpha_x,x)\operatorname{eq}(\alpha_y,y)V_5(x)V_5(y),
$$

$$
K^Q_D(a,b;\alpha_x,\alpha_y)=
\sum_{x,y:\,a=x+2,\,b=y}
\operatorname{eq}(\alpha_x,x)\operatorname{eq}(\alpha_y,y)V_5(x)V_5(y),
$$

and

$$
K^R_D(a,b;\alpha_x,\alpha_y)=
\sum_{x,y:\,a=x,\,b=y}
\operatorname{eq}(\alpha_x,x)\operatorname{eq}(\alpha_y,y)V_5(x)V_5(y).
$$

All $x$-coordinates in the neighbor equations are modulo $5$ after restricting to valid rows.
In characteristic two the sign of $R$ is the same as plus, so $\lambda_R$ is added without a minus.

This is the Keccak replacement for binius64's inner/shift sumcheck.
binius64 must prove that arbitrary shifted operands evaluate to one committed witness MLE evaluation.
Here the operand-to-state map is fixed and public, so the verifier can update the kernel directly.
The only prover work in this inner step is bookkeeping: compute the new scalar claim and carry the new transparent kernel descriptor.

### Bit-kernel rank caveat

There is one important implementation caveat.
The 128-evaluation bit-axis univariate skip is directly compatible with the Spartan outer only when the carried bit kernel is a point-evaluation functional, or another rank-one form that lets the prover reduce

$$\sum_z K(z)(P(z)Q(z)-R(z))$$

to a single evaluation

$$P(\zeta)Q(\zeta)-R(\zeta).$$

The first boundary claim can be chosen in this form.
After one Keccak-aware inner pushback, the bit kernel has been permuted by rho rotations and mixed by theta transpose.
The result is still transparent and cheap to evaluate, but it is generally a 64-entry lane-bit vector, not one point-evaluation vector.

For this reason, the serial test harness should support a general carried lane-bit kernel.
The exact path is an ordinary Spartan sumcheck over the bit variables, lane variables, and batch variables, with the carried kernel as a transparent multiplicative factor.
The rank-one univariate-skip path is useful as a special case, but the Binius-native v0 target in §13 should instead use the residual-zero NTT shape from the production BitAnd path.

### Linear map convention

We define the pre-chi layer in the normal Keccak execution order:

$$D^{(t)} = (\rho+\pi)(\theta(A^{(t)})).$$

This is the clean protocol definition.
Although rotations distribute over XOR, we should not define the protocol by commuting rotations through theta.
The offsets are tied to source lanes and then pi moves lanes, so the normal order is less error-prone and matches the FIPS implementation.

The inner reduction works in the reverse direction on verifier kernels.
If a claim has the form

$$\sum_{x,y,z,u} K_D(x,y,z,u)D^{(t)}_{x,y}(z,u),$$

then the inner step computes a kernel $K_A$ such that

$$\langle K_D, D^{(t)} \rangle = \langle K_A, A^{(t)} \rangle.$$

This is applying the transpose of the fixed linear map $(\rho+\pi)\circ\theta$:

$$K_A = \theta^T(\rho+\pi)^T K_D.$$

For $\pi$ and $\rho$, this is just inverse lane movement plus inverse bit rotation on the kernel.
For $\theta$, it is a small fixed sparse formula on the $5 \times 5$ lane kernels.

### Rho-pi transpose

Let rho+pi send source lane $(a,b)$ to destination lane

$$\pi(a,b) = (b, 2a+3b \bmod 5),$$

with rho offset $r[a,b]$.
Equivalently, for destination $(x,y)$ the source is

$$a = x + 3y \bmod 5, \qquad b = x.$$

If

$$D_{x,y}(z,u)=B_{a,b}(z-r[a,b],u),$$

then the transpose sends a destination kernel $K_D$ to a source kernel $K_B$ by

$$K_B(a,b)(z,u) \;{+}{=}\; K_D(b,2a+3b)(z+r[a,b],u).$$

All bit indices are modulo $64$.
Thus the prover/verifier implementation is just one lane permutation plus one cyclic shift of the 64-entry bit kernel per lane.

### Theta transpose

Theta is:

$$C_x(z,u)=\sum_{y=0}^{4} A_{x,y}(z,u),$$

$$T_x(z,u)=C_{x-1}(z,u)+C_{x+1}(z-1,u),$$

$$B_{x,y}(z,u)=A_{x,y}(z,u)+T_x(z,u).$$

Here $z-1$ is modulo $64$, corresponding to the one-bit rotate-left in theta.

Given a kernel $K_B$ on the theta output lanes, define the column-kernel sums

$$U_x(z,u)=\sum_{y=0}^{4} K_B(x,y,z,u).$$

Then theta transpose is:

$$K_A(a,v,z,u)=K_B(a,v,z,u)+U_{a+1}(z,u)+U_{a-1}(z+1,u).$$

The $x$ indices are modulo $5$ and the $z$ index is modulo $64$.
The last term is the transpose of the one-bit rotate-left in theta, i.e. a one-bit rotate-right on the kernel.

This formula is the main reason the Keccak inner is simpler than binius64 Shift Reduction.
There are only 25 lane kernels and two fixed cyclic shifts, not a general table of shifted witness operands.

For the serial harness, this inner step does not need its own sumcheck.
It is a deterministic verifier-side kernel update plus matching prover bookkeeping.
If a later globalized protocol commits all intermediate trace layers and wants to batch many linear-layer checks at once, we can reintroduce a matrix or lincheck-style sumcheck.
That is a different mode from the serial claim-propagation design.

## 7. Serial one-round claim flow

This is the likely v0 protocol for one round $t$, run backward from $t=23$ to $t=0$.

Input to round $t$:

- a carried claim $S_{\mathrm{next}}$ on $A^{(t+1)}$ with transparent kernel $K^{(t+1)}$;
- the transparent iota correction for round $t$;
- access to the current round-boundary state polynomials through the eventual commitment/opening layer.

Step 1: compute the transparent correction

$$I_t=\sum_{z,u} K^{(t+1)}_{0,0}(z,u)\mathrm{RC}_t(z).$$

Step 2: run the Spartan-outer sumcheck for

$$S_{\mathrm{outer}}=\sum_{x,y,z,u}K^{(t+1)}_{x,y}(z,u)(P Q-R).$$

The verifier checks that

$$S_{\mathrm{outer}} = S_{\mathrm{next}} + I_t.$$

Step 3: after the outer reduces to claims on $P,Q,R$ at its challenge point, batch those three claims with fresh lambdas.
This gives one linear claim on the three $D$-derived columns.

Step 4: convert the $P,Q,R$ claim into a kernel on $D^{(t)}$:

- the $P$ claim contributes to lane $(x+1,y)$ and also contributes a transparent constant term;
- the $Q$ claim contributes to lane $(x+2,y)$;
- the $R$ claim contributes to lane $(x,y)$ with the appropriate sign, which is the same as plus in characteristic two.

Step 5: push that $D$ kernel backward through $(\rho+\pi)^T$ and $\theta^T$ using the formulas above.
The result is the carried claim on $A^{(t)}$.

After 24 rounds, the final carried claim is on the committed input state $A^{(0)}$.
The original input claim at $t=23$ is on the committed output state $A^{(24)}$.

## 8. Prover message computation target

The outer prover should materialize or stream three virtual columns $P,Q,R$ for the current round.
The most direct v0 implementation keeps the right-hand-side claim on $A^{(t+1)}$ outside the outer composition, as in §3 and §7.

### Important difference from binius64 BitAnd

binius64 BitAnd proves a residual polynomial that is zero on all 64 base bit positions:

$$A(z,x)B(z,x)-C(z,x)=0.$$

Because of that, its first univariate-skip message only sends the 64 upper-half evaluations on the size-128 domain.

Here, by design, the outer composition is only

$$P(z,x,y,u)Q(z,x,y,u)-R(z,x,y,u).$$

This is **not** zero on the 64 base bit positions.
It equals the next-state value plus iota:

$$P Q - R = A_{\mathrm{next}}+\mathrm{Iota}.$$

Therefore, if we keep $A_{\mathrm{next}}+\mathrm{Iota}$ only as a **batched input claim** on the right-hand side, the bit-axis univariate message for the left-hand side should be treated as a general degree-126 univariate.
The prover should send enough evaluations to bind it, naturally the full 128 evaluations on the doubled binary subspace.

There is an important nuance.
On the 64 base-domain points, the lower evaluations of the left-hand side are semantically equal to $A_{\mathrm{next}}+\mathrm{Iota}$.
But the verifier does not generally know those 64 values individually.
In serial mode it knows a carried claim, i.e. one weighted aggregate of them at the bit challenge, not the whole vector of 64 lower evaluations.
So those lower values are not free interpolation data unless we change the protocol object being carried.

There is an alternative:
include $A_{\mathrm{next}}+\mathrm{Iota}$ inside the residual polynomial

$$H = P Q - R - A_{\mathrm{next}}-\mathrm{Iota}.$$

Then $H$ vanishes on the 64 base bit positions and we recover binius64's "send only upper 64" optimization.
But this also pulls $A_{\mathrm{next}}$ into the outer composition and final evaluation interface.
For v0, prefer the cleaner claim flow:

$$\text{outer proves } S_{\mathrm{outer}}, \qquad \text{verifier checks } S_{\mathrm{outer}} = S_{\mathrm{next}} + I_t.$$

The cost is at most 64 extra field elements in the first bit-skip message per round.
This is small compared with the implementation simplification and keeps the next-state value as a carried input claim.

A third possibility is to carry a stronger object: the whole 64-vector of lower bit-axis claims for $A_{\mathrm{next}}+\mathrm{Iota}$.
Then the verifier could use those values as the lower-half interpolation data and ask the prover only for the upper half of $P Q - R$.
This saves 64 field elements in the first message but makes the carried state much larger, so it is not the default design.

### Tradeoff for upper-half-only univariate skip

There are two ways to make the verifier know the 64 lower evaluations in every layer.
Both change the layer invariant.

The first option is the residual-zero invariant:

$$H = P Q - R - A_{\mathrm{next}}-\mathrm{Iota}.$$

Then $H$ is pointwise zero on the 64 base bit positions.
The verifier knows the lower half for free and the prover only sends the upper 64 evaluations, exactly as in binius64 BitAnd.
The cost is that $A_{\mathrm{next}}$ becomes part of the outer composition and therefore part of the final evaluation interface of the outer sumcheck.
This is clean when every layer boundary is committed or otherwise externally bound.
It is less clean in the serial virtual-chain design, where $A_{\mathrm{next}}$ is intended to appear only through the carried output-side claim.

The second option is a 64-vector carried invariant.
Instead of carrying one value

$$S_{\mathrm{next}}=\sum_{x,y,z,u} K(x,y,z,u)A_{\mathrm{next}}(x,y,z,u),$$

the verifier carries all 64 base-bit claims

$$S_{\mathrm{next},z}=\sum_{x,y,u} K_z(x,y,u)A_{\mathrm{next}}(x,y,z,u).$$

Then the lower half of the next outer univariate is known pointwise.
But for this invariant to compose across GKR layers, each layer must also output a new 64-vector for $A^{(t)}$.
That means either 64 parallel pushbacks or a polynomial-valued carried kernel.
This likely gives back most of the 64 field elements saved in the first message, and it makes the transcript and verifier state less uniform.

For v0, the recommended invariant is therefore a scalar claim plus its transparent kernel descriptor.
The verifier knows the scalar value at the sampled bit challenge, not the 64 base-domain values.
This keeps every layer interface uniform and keeps $A_{\mathrm{next}}$ outside the current outer witness columns.
The price is that the first bit-skip message sends the full 128 evaluations.

### Why an aggregate constant is not enough for upper-half-only

There is a useful algebraic normalization.
Let

$$C_{\mathrm{next}}=\sum_{\kappa}\beta(\kappa)(A_{\mathrm{next}}(\kappa)+\mathrm{Iota}(\kappa)).$$

If $\sum_{\kappa}\beta(\kappa)$ is nonzero, then

$$\sum_{\kappa}\beta(\kappa)(P Q-R-C_{\mathrm{next}}/\sum_{\kappa}\beta(\kappa)) = 0.$$

This is valid as a scalar sumcheck claim.
It lets us turn the outer into a zero-sum claim without putting $A_{\mathrm{next}}$ into the integrand.

However, it still does **not** make the 64 lower bit-axis evaluations known pointwise.
It only enforces one aggregate linear relation among them.
For upper-half-only interpolation, the verifier needs 64 point values on the base half (as in binius64, where they are all zero), not just one weighted sum of those values.
So the aggregate-constant trick can simplify the initial claim, but it does not recover the binius64 64-evaluation first message.

This is probably the best v0 invariant:

- a layer proof consumes one carried scalar claim on $A^{(t+1)}$, together with the transparent kernel descriptor defining that scalar;
- it uses that scalar claim to set the starting claim of the Spartan-outer;
- the outer and Keccak-aware inner output one new carried scalar claim on $A^{(t)}$, again with its transparent kernel descriptor.

In this invariant, every layer has the same interface.
The carried object is a scalar linear claim plus transparent kernel data, not a whole lower-half vector and not a new witness column.

### Historical rank-one bit-kernel path

When the carried bit kernel is a point-evaluation kernel, one round can use the 128-evaluation bit-axis skip.
This is the right special case for an initial boundary claim or a segmented mode that resets the kernel at committed boundaries.
It is not, by itself, the complete serial protocol after arbitrary Keccak inner pushback.
It is also not the Binius-native v0 target; §13 prefers the residual-zero BitAnd-style NTT path.

For one round $t$ in this fast path:

1. Build the 25 pre-chi lanes $D_{x,y}$ for all batch items by executing theta, rho, and pi over `u64` words.
2. Form virtual operand streams:
   - $P_{x,y} = 1 + D_{x+1,y}$;
   - $Q_{x,y} = D_{x+2,y}$;
   - $R_{x,y} = D_{x,y}$.
3. Embed the valid 25 lane positions into the padded $8 \times 8$ lane grid.
   Invalid lane rows may contain arbitrary zeros because the validity selector in $\beta_{\mathrm{pad}}$ kills them.
4. Run the bit-axis univariate-skip first round for

   $$G_0(Z)=\sum_{x,y,u}\widetilde K^{(t+1)}_{x,y}(Z,u)\left(P(Z,x,y,u)Q(Z,x,y,u)-R(Z,x,y,u)\right).$$

   Here $\widetilde K^{(t+1)}$ is the carried kernel evaluated on the bit-axis extension point $Z$.
   In the initial factorized case it specializes to the earlier product of lane, bit, batch, and validity weights.
   This is the left-hand side only.
   Since it is not known to vanish on the base half, v0 sends the full 128 evaluations on the doubled binary subspace.
5. Verifier samples the bit-axis challenge $\zeta$ and evaluates $G_0(\zeta)$ from the sent univariate.
   It checks this value against the carried right-hand-side claim:

   $$G_0(\zeta)=S_{\mathrm{next}}+I_t.$$

6. Prover folds each word of $P,Q,R$ along the bit axis at $\zeta$ using the same byte-table/Four-Russians idea as binius64 `FoldLookup`.
   This produces three scalar tables over the padded lane grid and batch axis:

   $$P_\zeta(x,y,u), \qquad Q_\zeta(x,y,u), \qquad R_\zeta(x,y,u).$$

7. Run the remaining degree-2 MLE-check over the $3+3+\ell$ Boolean variables for padded $(x,y,u)$:

   $$G_0(\zeta)=\sum_{x,y,u}\widetilde K^{(t+1)}_{x,y}(\zeta,u)\left(P_\zeta Q_\zeta-R_\zeta\right).$$

   This is a standard three-polynomial composition, with the carried kernel treated as a transparent factor.
   At the end the prover supplies evaluations of $P_\zeta,Q_\zeta,R_\zeta$ at the sumcheck point.
8. Batch those final three evaluations with fresh lambdas and push the resulting $D$-kernel back to $A^{(t)}$ via the Keccak-aware inner formulas in §6.

This is intentionally simple.
It computes theta before rotation, exactly like Keccak execution.
Later optimizations may avoid materializing all $D$ lanes by pushing rotations through the XOR-linear theta computation, but the protocol should not depend on that optimization.

The inner prover does not need binius64's generic `KeyCollection`.
It only needs to maintain 25 transparent kernels and apply fixed permutations and cyclic shifts.

### Data layout for v0

For one round, store pre-chi lanes as:

```text
D[batch_index][lane_index] : u64
lane_index = x + 5*y, 0 <= lane_index < 25
```

For the outer sumcheck, expose padded lane slots:

```text
padded_lane = x + 8*y, 0 <= x < 8, 0 <= y < 8
```

The mapping for valid slots is:

```text
valid if x < 5 and y < 5
lane_index = x + 5*y
```

The virtual operands can be accessed without materializing all three arrays:

```text
P_word(x,y,u) = !D_word((x+1) mod 5, y, u)   # bitwise 1 + D
Q_word(x,y,u) =  D_word((x+2) mod 5, y, u)
R_word(x,y,u) =  D_word(x, y, u)
```

For invalid padded slots, return zero.
The validity selector also zeroes their contribution.

### First-round univariate computation

The bit-axis first round is where we should reuse the binius64 idea most directly.
For each `(x,y,u)` row, the prover has three 64-bit words:

```text
p = P_word(x,y,u)
q = Q_word(x,y,u)
r = R_word(x,y,u)
```

The byte-lookup NTT maps each word from 64 base-domain bits to 128 output-domain evaluations.
For each output-domain point $Z$ and row $(x,y,u)$, accumulate

$$\widetilde K^{(t+1)}_{x,y}(Z,u)\left(P(Z)Q(Z)-R(Z)\right).$$

The lower 64 evaluations are especially cheap.
On a base-domain bit $z$, $P(z),Q(z),R(z)\in\{0,1\}$, so the row contribution is the bit of

```text
h_word = (p & q) ^ r
```

at position $z$, multiplied by the carried kernel value $\widetilde K^{(t+1)}_{x,y}(z,u)$.
Thus a direct implementation can compute the lower half by walking the valid rows and adding the appropriate kernel value into 64 accumulators selected by the set bits of `h_word`.
When the kernel factors into a row scalar times a 64-entry bit kernel, this is just row-scalar multiplication followed by bytewise table accumulation.
More generally, the carried kernel descriptor should expose the 64 base-bit values for each active lane/batch row, or expose a byte-table evaluator for them.
This is the same 256-entry Four-Russians pattern as `FoldLookup`, but applied to weighted bit accumulation rather than a single folded scalar.
This path avoids interpolating the low half through field arithmetic.

The upper 64 evaluations need the actual degree-2 univariate extension.
For them, use the binius64 NTT-lookup shape:

1. expand `p`, `q`, and `r` from the 64-bit base domain to the upper half of the doubled binary subspace;
2. evaluate the carried kernel descriptor $\widetilde K^{(t+1)}_{x,y}(Z,u)$ on the same upper points;
3. compute `P(Z) * Q(Z) - R(Z)` at each upper point;
4. multiply by the kernel value and accumulate into the upper-half message.

In v0 the prover sends `lower || upper`.
In the residual-zero variant, the lower array is known to the verifier as all zeros, so only the upper array is sent.

This is exactly binius64's Phase-1 arithmetic shape, except:

- there are only Keccak lane rows, not arbitrary constraint rows;
- the row weight is the carried transparent kernel, with the product of lane, bit, batch, and validity factors as the simplest special case;
- v0 sends all 128 evaluations, not only the upper 64, because the left-hand side alone does not vanish on the base domain;
- the lower 64 evaluations can be computed directly from `h_word = (p & q) ^ r`, while the upper 64 use the NTT lookup.

Once this is implemented, the obvious optimization experiment is the residual variant:
include $A_{\mathrm{next}}+\mathrm{Iota}$ in the same NTT loop and send only the upper 64 evaluations.
That may be worth it if transcript size or first-round interpolation becomes visible in benchmarks.

## 9. Serial claim invariant for tests

This section records a useful algebraic test harness, not the v0 implementation target.
It prioritizes a uniform serial claim invariant over the bit-axis univariate-skip optimization, so it is good for catching sign, rotation, and bit-order bugs.
The production implementation should instead follow the Binius-native target in §13.

### Statement

The parent protocol commits to, or otherwise binds, the boundary states $A^{(0)}$ and $A^{(24)}$.
The Keccak sub-protocol consumes one output-boundary claim

$$\mathsf{Claim}(A^{(24)},K^{(24)},s^{(24)})$$

and reduces it to one input-boundary claim

$$\mathsf{Claim}(A^{(0)},K^{(0)},s^{(0)}).$$

The parent opening layer checks both boundary claims against the committed boundary oracles.
Equivalently, the Keccak sub-protocol proves that the output claim is consistent with applying 24 Keccak rounds to the input state.

### Kernel representation

For the serial harness, represent the carried kernel as

```text
Kernel {
    batch_point: [F; ell],
    lane_bit: [[F; 64]; 25],
}
```

It evaluates as

$$K_{x,y}(z,u)=K^{\mathrm{lane}}_{x,y}[z]\cdot \operatorname{eq}(\rho,u),$$

where $\rho$ is `batch_point`.
This representation is closed under the Keccak-aware inner pushback.
Rho and pi permute the 25 lane-bit arrays and cyclically shift their 64 entries.
Theta transpose forms column sums of lane-bit arrays and adds two shifted column sums.
The batch point changes only when the outer sumcheck samples a new batch-axis point.

The initial output claim should be converted into this form before the first round.
For example, a parent claim with lane coefficients $\gamma_{x,y}$, bit point $\zeta$, and batch point $\rho$ has

$$K^{\mathrm{lane}}_{x,y}[z]=\gamma_{x,y}L_z(\zeta).$$

### Transcript prelude

The verifier and prover initialize the transcript with:

1. a domain separator for `keccak-prove.serial-v0`;
2. the batch size $h=2^\ell$;
3. the boundary commitment identifiers for $A^{(0)}$ and $A^{(24)}$;
4. the initial output claim scalar $s^{(24)}$ and a canonical encoding of $K^{(24)}$.

The transcript should also bind the choice of padded Boolean lane domains and the use of serial mode $s_{\mathrm{seg}}=1$.

### One round

Rounds are processed backward, for $t=23,\ldots,0$.
At the start of round $t$, the verifier holds

$$\mathsf{Claim}(A^{(t+1)},K^{(t+1)},s^{(t+1)}).$$

The prover has witness access to the batch of states at boundary $t$ and can compute the pre-chi words $D^{(t)}$.
The transcript absorbs a round separator containing $t$, the carried scalar $s^{(t+1)}$, and a canonical encoding or digest of the carried kernel descriptor $K^{(t+1)}$.

First, both sides compute the transparent iota correction

$$I_t=\sum_{z,u}K^{(t+1)}_{0,0}(z,u)\mathrm{RC}_t(z).$$

The starting outer claim is

$$c_0=s^{(t+1)}+I_t.$$

The outer sumcheck proves

$$
c_0=
\sum_{z,x,y,u}
K^{(t+1)}_{x,y}(z,u)
\left(P^{(t)}_{x,y}(z,u)Q^{(t)}_{x,y}(z,u)-R^{(t)}_{x,y}(z,u)\right).
$$

The sum is over the padded Boolean domain for $z,x,y,u$:
six bit variables for $z$, three variables for $x$, three variables for $y$, and $\ell$ variables for $u$.
The virtual columns $P,Q,R$ are zero on invalid padded lane rows.
The round order should be fixed as

```text
z[0..6), x[0..3), y[0..3), u[0..ell)
```

where `z[0]` is the least significant bit of the bit index.

For each sumcheck round $i$:

1. the prover sends the univariate round polynomial $g_i(X)$ and the transcript absorbs its coefficients;
2. the verifier checks $g_i(0)+g_i(1)=c_i$;
3. the verifier squeezes challenge $r_i$;
4. both set $c_{i+1}=g_i(r_i)$.

The round polynomial has degree at most three because it contains the transparent kernel factor and the quadratic Spartan composition.
At the end, the verifier has the point

$$r=(r_z,r_x,r_y,r_u).$$

The prover sends

$$p=P^{(t)}(r),\qquad q=Q^{(t)}(r),\qquad r_{\mathrm{sp}}=R^{(t)}(r).$$

The transcript absorbs $(p,q,r_{\mathrm{sp}})$ before any batching challenge is sampled.

The verifier evaluates the transparent factors at $r$ and checks

$$c_{\mathrm{final}}
=
K^{(t+1)}(r)(p q-r_{\mathrm{sp}}).$$

### Inner pushback

After the outer check, the verifier squeezes fresh batching challenges

$$\lambda_P,\lambda_Q,\lambda_R.$$

The batched virtual-column claim is

$$s_{PQR}=\lambda_P p+\lambda_Q q+\lambda_R r_{\mathrm{sp}}.$$

The constant part of $P=1+D_{x+1,y}$ contributes

$$c_P=V_5(r_x)V_5(r_y).$$

Therefore the induced claim on $D^{(t)}$ has scalar

$$s_D=s_{PQR}+\lambda_P c_P.$$

The induced kernel $K_D$ is the point-evaluation kernel at $(r_z,r_u)$, multiplied by the lane coefficients determined by $(r_x,r_y)$ and $(\lambda_P,\lambda_Q,\lambda_R)$:

$$
K_D(a,b,z,u)=\operatorname{eq}(r_z,z)\operatorname{eq}(r_u,u)
\cdot\left(
\lambda_P K_D^P(a,b;r_x,r_y)
+\lambda_Q K_D^Q(a,b;r_x,r_y)
+\lambda_R K_D^R(a,b;r_x,r_y)
\right).
$$

Here $K_D^P,K_D^Q,K_D^R$ are the valid-row neighbor kernels from §6.
The verifier computes

$$K^{(t)}=\theta^T(\rho+\pi)^T K_D$$

using the fixed formulas in §6.
The carried scalar is

$$s^{(t)}=s_D.$$

This completes round $t$.

### Protocol output

After the loop reaches $t=0$, the Keccak sub-protocol outputs

$$\mathsf{Claim}(A^{(0)},K^{(0)},s^{(0)}).$$

The parent protocol batches this input-boundary claim with the original output-boundary claim and any other oracle-opening claims.
In a standalone test harness, the verifier can check both claims directly by evaluating the explicit input and output traces.

### Serial test checklist

If we build the serial harness, it should land in this order:

1. `trace`: compute all 25 pre-chi words $D^{(t)}$ from one round input.
2. `kernel`: implement `Kernel`, rho-pi transpose, theta transpose, and direct claim evaluation against a trace.
3. `outer`: implement the ordinary degree-3 sumcheck over `(z,x,y,u)` with virtual accessors for $P,Q,R$.
4. `round`: implement one backward round, including transcript order, final $P,Q,R$ eval batching, and kernel pushback.
5. `chain`: run 24 rounds and compare the final input claim against a reference Keccak permutation trace.
6. compare its carried claims against the parallel v0 implementation as a debugging oracle.

The main tests should be:

- transpose identities for rho-pi and theta;
- kernel representation closure across one inner pushback;
- one-round sumcheck against a direct table sum on small batches;
- 24-round claim reduction against direct Keccak-f[1600];
- transcript determinism between prover and verifier.

## 10. Round segmentation parameter

We should support a segment length parameter $s_{\mathrm{seg}}$ from the start.
This controls how many Keccak rounds are virtualized between committed or externally supplied boundary claims.

- **Serial mode:** $s_{\mathrm{seg}} = 1$.
  Each round consumes an input claim on $A^{(t+1)}$, proves the $K=3$ outer relation for round $t$, and pushes claims back to $A^{(t)}$.
  No intermediate commitments are needed.
- **Fully parallel mode:** $s_{\mathrm{seg}} = 24$.
  One global outer covers all rounds.
  This likely requires committing to, or otherwise binding, all intermediate $A^{(t)}$ states, because all right-hand-side $A^{(t+1)}$ values appear as outer input claims at once.
- **Segmented mode:** $1 < s_{\mathrm{seg}} < 24$.
  Commit or externally bind only segment boundaries, and virtualize the rounds inside each segment by repeated claim pushback.

The protocol should expose this knob even if v0 only implements one value.
The old conservative v0 choice was serial mode, because it avoided committing intermediate layers and kept the inner step as transparent claim propagation.
That is no longer the right default now that the implementation lives inside Binius64.
The Binius-native v0 should start with a parallel or segmented-parallel shape and only use serial mode as a correctness harness.
Benchmarks should tune $s_{\mathrm{seg}}$ after the NTT-backed implementation is working.

## 11. Relation to binius64

binius64 has the same two-level architecture:

1. **BitAnd / Spartan-outer:** prove $A \cdot B - C = 0$ with only three operand polynomials.
2. **Shift Reduction / Spartan-inner:** reduce claims about operand polynomials $A,B,C$ to claims about the underlying witness words.

This note aims for the same shape. The current implementation uses production Shift directly:

- binius64 Shift handles arbitrary shifted-word constraints;
- Keccak chi operand rows instantiate a small fixed subset of that machinery;
- committed `D` words turn theta into separate linear correctness rows rather than an uncommitted pushback through theta.

Future specialization can exploit Keccak's fixed maps after the production-faithful chain is benchmarked.
The outer should stay Spartan-like with a small number of operand polynomials.

## 12. Remaining decisions before protocol text

- How large the first round segment should be: all 24 rounds in one global outer, or smaller segments such as 4, 6, 8, or 12 rounds.
- Whether the residual-zero variant should be the default because it recovers BitAnd's upper-half-only first message.
- How much generic Shift overhead remains visible after the production-faithful path is complete and benchmarked.
- Whether exact mixed domains for $(t,x,y)$ are worth implementing after padded Boolean domains are measured.

## 13. Binius-native v0 implementation target

The v0 implementation should live in `crates/keccak-prove`.
It should be ambitious because it is now inside the Binius64 workspace and can reuse production prover internals directly.

The target is not a slow serial teaching implementation.
The target is a Keccak-specific prover that adapts the production BitAnd and shift-reduction architecture:

1. **BitAnd-style outer.**
   Keep the $K=3$ operand shape

   $$
   P = 1 + D_{x+1,y}, \qquad Q = D_{x+2,y}, \qquad R = D_{x,y},
   $$

   and prove the chi relation with the same broad machinery used for word-level AND constraints.
2. **NTT lookup from the start.**
   The bit axis should use the byte-table additive-NTT/Four-Russians path immediately.
   This is not a later fast path.
   The implementation should adapt `crates/prover/src/and_reduction/ntt_lookup.rs` and the surrounding BitAnd first-round message code.
3. **Production Shift inner first.**
   Reuse the production shift-reduction implementation directly for v0.
   Keccak has fixed rho offsets, a fixed pi permutation, fixed chi neighbors, and committed theta corrections, so every virtual `B` reference can be represented as a small Shift operand over committed `A` and `D` words.
   A later Keccak-specialized inner may remove some generic `KeyCollection` overhead, but that should be an optimization after the production-faithful chain is complete and benchmarked.
4. **Parallel or segmented-parallel proving.**
   The first prover should batch across many Keccak permutations and should aim to batch across rounds as much as the claim interface permits.
   A segmented mode is acceptable if one global 24-round outer forces too much boundary machinery too early.
5. **Early benchmarkability.**
   The crate should benchmark against the existing generic Keccak circuit path in `binius_circuits` early, even before every optimization is final.

### Crate layout

Use a new workspace crate:

```text
crates/keccak-prove/
  Cargo.toml
  src/
    lib.rs
    constants.rs
    trace.rs
    operands.rs
    bit_ntt.rs
    round_message.rs
    layout.rs
    witness.rs
    shift_operands.rs
    shift_claims.rs
```

`constants.rs` owns the FIPS constants, rho offsets, lane-index helpers, and padded-domain constants.

`trace.rs` owns native Keccak execution and compact word-level traces.
It should be explicit and testable, but it is not the main proving abstraction.

`operands.rs` owns virtual accessors for $P,Q,R$ from pre-chi words.
It should support padded lane rows and validity selectors without materializing unnecessary arrays.

`bit_ntt.rs` owns the Keccak-adapted bit-axis first-round logic.
It should start by copying the shape of the production BitAnd NTT lookup code, then generalize only where Keccak requires row weights, round/lane selectors, or iota residual terms.

`round_message.rs` owns the Spartan/BitAnd-style outer relation over the chosen segment domain.
It should be written so the segment axis can be one round, several rounds, or all 24 rounds.

`layout.rs` and `witness.rs` own the committed `A`/`D` witness layout and construction.

`shift_operands.rs` owns virtual `B` lowering and `D` correctness operands.

`shift_claims.rs` owns the production Shift schemas for chi operand pushback and `D` correctness.

Future `segment.rs` should own the segment boundary object and the handoff between outer and Shift reductions.

Future `prove.rs` should own the high-level prover/verifier entrypoints and transcript order.

### Residual-zero default

For the Binius-native implementation, prefer the residual-zero form:

$$
H = P Q - R - A_{\mathrm{next}}-\mathrm{Iota}.
$$

On the 64 base bit positions, $H$ is zero.
That matches the production BitAnd invariant and lets the prover send only the upper half of the doubled-domain NTT message.
This is more invasive than the scalar carried-claim invariant because $A_{\mathrm{next}}$ participates in the outer relation, but it is the right v0 bias inside Binius64.

If this makes the all-24-round global form awkward, use segmented proving:

- commit or otherwise bind segment boundaries;
- run the residual-zero BitAnd-style outer inside each segment;
- use the committed `A`/`D` plus production-Shift path inside the segment.

### Reuse points in Binius64

The implementation should read and adapt these production paths:

- `crates/prover/src/and_reduction/`
  for BitAnd first-round messages, sumcheck flow, and NTT lookup shape;
- `crates/prover/src/fold_word.rs`
  for byte-table/Four-Russians word folding;
- `crates/prover/src/protocols/shift/`
  for the generic shifted-word inner architecture;
- `crates/ip-prover/src/sumcheck/`
  for public sumcheck prover interfaces;
- `crates/iop/src/channel.rs`
  for final `OracleLinearRelation` boundary integration;
- `crates/circuits/src/keccak/`
  for existing circuit behavior and benchmark comparison.

The goal is to specialize production code, not to create a parallel toy implementation.

### Correctness guardrails

Even though the v0 prover should be parallel and NTT-backed, keep the cheap algebraic tests from the serial plan:

- rho/pi transpose identity;
- theta transpose identity;
- chi/iota residual identity on concrete states;
- one-segment direct table sum versus prover message on tiny batches;
- 24-round native Keccak trace agreement.

These tests should run before performance benchmarks.
They are not the architecture, but they are the seatbelts.

### Near-term milestone

The first milestone should be:

```text
cargo test -p binius-keccak-prove keccak_bitand_residual_matches_native_round
```

where the test:

1. samples a small power-of-two batch;
2. builds native Keccak pre-chi and next-state words for one or more rounds;
3. constructs $P,Q,R,A_{\mathrm{next}},\mathrm{Iota}$;
4. verifies that the residual is zero on all 64 base bit positions;
5. verifies that the upper-half NTT message agrees with a direct extension-domain computation.

The second milestone should prove and verify one segment using the Binius transcript and sumcheck machinery.
The third milestone should benchmark against the generic Keccak circuit path.
