# How the centres and subdivision were found

This document records the numerical search used to propose the finite data.
It is supplementary to the article *Borsuk number in four dimensions is
at most seven*. `VERIFICATION.md` describes the exact acceptance test;
`README.md` gives the commands.

The finder and verifier have different roles. The finder proposes centres,
chooses subdivisions, and suggests radii. The verifier proves coverage for
each retained proposal. Numerical sampling, tolerances, optimizer status,
and an unsuccessful local search never certify or disprove coverage of the
continuous outer set.

The target radius is $\rho=9999/20000$, its square is $\sigma=\rho^2$,
and the computational range is $299/1000\le r\le r_0=38401/100000$.
All notation for the parameter cube, contacts, and rational outer set is
defined in `VERIFICATION.md`.

## 1. Sampling a rational outer set

For a current parameter box $Q$, the rational outer set is

$$
\Omega_Q=\{x:\|x\|^2\le r_Q^+,\quad
\|x\|^2-2\langle A_i,x\rangle\le b_i\ (i=1,2,3)\}.
$$

The numbers $b_i$ are positive for this family, so $0\in\Omega_Q$.
Along a unit direction $u$, its radial boundary is at $R(u)u$, where

$$
R(u)=\min\left\{\sqrt{r_Q^+},
\min_{1\le i\le3}
\left(\langle A_i,u\rangle+
\sqrt{\langle A_i,u\rangle^2+b_i}\right)\right\}.
$$

This follows by solving each quadratic inequality along the ray
$x=tu$, $t\ge0$. The finder generates directions from scrambled Sobol
points in $[0,1]^4$: it clips coordinates away from zero and one, applies
the inverse standard-normal distribution function coordinatewise, and
normalizes the resulting four-dimensional vectors.

The default optimization sample contains $2^{11}$ boundary points, the
origin, and half-radius copies of every sixteenth boundary point. A second
boundary sample of $2^{12}$ points locates poorly covered regions. Samples
are deterministic for their numerical seeds, but they are not covering
certificates. The exact test covers the interior as well as the boundary,
regardless of this sampling strategy.

## 2. Sampled minimax optimization

For a finite sample $X=\{x_i\}$, the objective is

$$
f_X(C_1,\ldots,C_7)=
\max_i\min_{1\le j\le7}\|x_i-C_j\|^2.
$$

All 28 coordinates of the centres are free. No final symmetry is imposed,
even though the outer set is invariant under reflection in the fourth
coordinate. The code alternates the following two procedures.

**Assignment and circumball updates.** Assign every sample point to
its nearest centre, then move each centre to an approximate centre of the
assigned cluster's circumball. The numerical circumball problem is solved
by SLSQP on an active pool. Initially the pool
contains up to twelve farthest points; further farthest points are added
when needed. Empty clusters retain their old centre. These solutions are
proposals, not exact circumball certificates.

**Smooth approximations.** For $d_{ij}=\|x_i-C_j\|^2$, temperature
$T>0$, and parameter $\gamma>0$, define

$$
g_i=-T\log\sum_{j=1}^7\exp(-d_{ij}/T),\qquad
f_{T,\gamma}=\frac1\gamma\log\sum_i\exp(\gamma g_i).
$$

The first expression approximates the minimum over centres and the second
approximates the maximum over points. The gradient is explicit. If

$$
v_{ij}=\frac{\exp(-d_{ij}/T)}{\sum_k\exp(-d_{ik}/T)},\qquad
w_i=\frac{\exp(\gamma g_i)}{\sum_k\exp(\gamma g_k)},
$$

then

$$
\nabla_{C_j}f_{T,\gamma}
=2\sum_iw_i v_{ij}(C_j-x_i).
$$

The implementation uses stable log-sum-exp evaluations and L-BFGS-B. The
three stages are:

| $T$ | $\gamma$ | Maximum L-BFGS-B iterations |
| --- | --- | --- |
| $0.003$ | $350$ | $65$ |
| $0.0006$ | $1600$ | $85$ |
| $0.00015$ | $6500$ | $90$ |

The first stage smooths the objective more strongly; later stages
approximate the sampled minimax objective more closely. In the code,
the function argument named `beta` denotes this numerical smoothing
parameter $\gamma$. It is unrelated to the obsolete fixed axial-centre
parameter from an earlier analytic proof.

Each refinement uses up to sixteen assignment-and-update steps, the three
smooth stages, and up to six further assignment-and-update steps. By
default there are at most three refinements for each of two starting
configurations. The local search time budget may interrupt this schedule.
The 64 worst points in the larger sample are added before another
refinement. Numerical vertices that appear troublesome in the lifted
test can supply additional radially projected sample points; they do not
by themselves prove that an uncovered point exists.

## 3. From a numerical proposal to an exact certificate

Each proposed coordinate is rounded to the nearest multiple of $2^{-20}$.
A numerical version of the lifted five-dimensional test first screens
the rounded centres. Its default acceptance threshold is a maximum
numerical vertex value below $-10^{-7}$. This threshold only decides
whether to attempt the exact test.

A promising proposal is passed to `exact_cover`, which checks the positive
centre simplex, every possible vertex basis, all feasibility inequalities,
and strict vertex negativity using rational arithmetic. Only this exact
success can mark a box as covered. Determinant thresholds and feasibility
tolerances used during numerical screening cannot weaken exact acceptance.

The fallback remains an unresolved box. A numerical failure can cause more
optimization or subdivision, but cannot create a successful certificate.

## 4. Reusing configurations and exploratory seeds

Neighbouring parameter boxes often admit the same centres. Before local
optimization, the finder tries the parent configuration and configurations
selected from the seed bank and a recent cache. Distance between models is
measured by the squared Euclidean distance between their three ordered
contact vectors. The default cache holds 128 configurations and up to four
nearby candidates are selected. Small common radial factors
$1$, $0.9998$, $0.9995$, and $1.0002$ are also tried, always followed by
rounding and the exact test when the numerical screen succeeds.

The supplied bank contains 28 exploratory configurations and 2,897
configurations from earlier searches. Only centre proposals are reused.
The current finder starts from the whole single parameter cube and imports
no earlier subdivision. The exploratory contact triangles are aligned with
the current coordinate convention by relabelling and an orthogonal
transformation of their first two coordinates.

The exploratory models use
$r=0.3335,0.35,0.37,0.3835$. For each $r$, the values of $\tau$ are
$5r/6$, the midpoint of $[5r/6,1/3]$, and $1/3$. For the first two
values, the initial triangle shapes are regular,
$A=B=\theta$, and $A=\theta$, $B=C=(\pi-\theta)/2$, where

$$
\theta=\arcsin\frac1{2\sqrt\tau}.
$$

Here $A,B,C$ are the three triangle angles. At $\tau=1/3$, only the
regular triangle is needed. This gives seven models at each of four
values of $r$.

Writing the triangle contacts as $P_i=(p_i',h,0)$ with
$p_i'\in\mathbb E^2$, the initial configurations used three centres
whose projections onto the first two coordinates were aligned with

$$
-\frac{p_j'+p_k'}{\|p_j'+p_k'\|},
\qquad \{i,j,k\}=\{1,2,3\},
$$

and four centres arranged on a circle in the complementary two-dimensional
coordinate plane. This describes an initialization only. It does not
assert that three of the final centres lie in the triangle's affine plane,
or impose a three-plus-four decomposition on the optimized centres. Small
perturbations and sampled minimax refinement produced the seed bank. These
isolated models did not certify the parameter gaps between them.

The earlier searches began with radius $1/2$ and later rechecked proposals
at $\rho$. They contributed the additional bank. A full verification of
the final certificate does not consult that search history.

## 5. Subdivision and default limits

When no certificate has been found, the box is bisected. For each possible
split coordinate, the finder computes the mean interval uncertainty of the
two children and chooses a coordinate with a small value. Here uncertainty
is the width of the enclosing $r$ interval plus the largest contact error
$E_i$. To avoid extreme imbalance, a coordinate is eligible only when its
dyadic level is at most eight above the smallest current level. Both
children are retained.

The default local optimization budget is eight seconds, with two starts
and three refinements. Optimization is attempted only when the outer-model
uncertainty is at most $0.035$; otherwise candidate reuse and subdivision
can continue. The default boundary sample exponent is 11 and the maximum
tree depth is 48. Run `find --help` for the command-line settings; the
remaining constants are collected in `SEARCH_DEFAULTS`.

The global `--hours` budget limits the search invocation, with current
worker jobs allowed to finish. Neither that limit nor the local budget is
a guarantee of completion. The finder saves a transactional SQLite
checkpoint. Its resume checks include the source-code hash, seed-file
hash, and search settings. A checkpoint must therefore be resumed with
the same script version that created it. The updated analytic checks in
this supplementary release change the source hash; they do not invalidate
the supplied stand-alone certificate, but an older checkpoint requires
its original script for resumption or export.

`--fresh-seeds` omits the earlier configuration bank. `--retry-blocked`
is available on resume to retry unresolved depth-limit leaves. A reported
`bottleneck_box` is already covered and records the largest retained
radius; completion is determined by the unresolved-leaf count and the
subsequent full replay.

## 6. Optional radius polishing

Once a box passes at squared radius $\sigma$, its rounded centres are
held fixed. Each candidate lifted vertex is affine in a new squared
radius $s$, because the constraint normals are fixed and the right-hand
sides are affine in $s$. Numerical linear solves can therefore be reused
for a 24-step bisection that suggests a smaller sufficient value.

The proposal is rounded upwards to a multiple of $10^{-7}$ and increased
by two more grid units. The exact verifier then checks that rational
value, including all vertex bases and boundedness. If it fails, the next
proposal is halfway towards the previously certified $\sigma$. At most
three exact proposals are tried; otherwise the original certificate is
retained. `--no-radius-polish` skips this step.

Every stored radius is an exactly proved upper bound. The sampled radius
is recorded separately with an explicit `heuristic` label. A retained
radius is not asserted to be optimal for those centres or for the outer
set.

## 7. What reproduction means

Checking the supplied data reproduces the finite proof. Searching again
is an experiment: scheduling, optimizer versions, seeds, and time limits
may change both the centres and subdivision. No byte-for-byte reproduction
of a new search is claimed. Every newly exported tree must be checked by
the separate exact `verify` command described in `README.md`.
