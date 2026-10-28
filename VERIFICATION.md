# Exact verification of the seven-ball covering certificate

This document contains the detailed interval and software arguments
supporting the article *Borsuk number in four dimensions is at most seven*.
It replaces the long technical appendices in the article. Commands and a
file inventory are in `README.md`.

**Repository URL to be replaced before publication:**
<https://github.com/USERNAME/REPOSITORY>

The notation below is independent of equation numbers in the article.
All balls are Euclidean. Write $\overline B(c,R)$ for a closed ball and
$B(c,R)$ for an open ball. The inner product and norm are denoted by
$\langle\cdot,\cdot\rangle$ and $\|\cdot\|$.

## 1. Constants and the parameter family

The circumball of a bounded set is its smallest closed enclosing ball.
For a set of diameter one in $\mathbb E^4$, the article places its
circumball centre at zero and denotes the **squared** circumradius by $r$.
The fixed constants are

$$
\rho=\frac{9999}{20000},\qquad
\sigma=\rho^2=\frac{99980001}{400000000},\qquad
r_-=\frac{299}{1000},\qquad r_0=\frac{38401}{100000}.
$$

Here $r_-$ is the fixed lower endpoint of the computational range; local
enclosure endpoints below are written $r_Q^-$ and $r_Q^+$. In the code,
`RADIUS`, `RADIUS_SQUARED`, `R_LO`, and `R_HI` denote $\rho$, $\sigma$,
$r_-$, and $r_0$, respectively. A proposed or stored squared covering
radius is denoted by $s$, with $0<s\le\sigma$.

The geometric reduction in the article produces three contact points
$P_1,P_2,P_3$ with common squared norm $r$, mutual distances at most one,
and squared triangle circumradius $\tau$ satisfying

$$
\frac56r\le\tau\le\min\{r,1/3\}.
$$

It is sufficient to cover the outer set

$$
H(r,P_1,P_2,P_3)
=\overline B(0,\sqrt r)\cap\bigcap_{i=1}^3\overline B(P_i,1).
$$

This outer set need not itself have diameter one. The certificate covers
the whole following family, including boundary and degenerate parameter
values, so it also covers every set required by the geometric reduction.

Put

$$
\delta_-=\frac13-\frac{299}{1000}=\frac{103}{3000},\qquad
\delta_+=r_0-\frac13=\frac{15203}{300000}.
$$

For $u=(u_1,u_2,u_3,u_4)\in[0,1]^4$, define

$$
\begin{aligned}
t&=u_1\sqrt{\delta_+}-(1-u_1)\sqrt{\delta_-},
&r&=\frac13+t|t|,\\
t_+&=\max\{t,0\},
&h&=(1-u_2)\sqrt{r/6}+u_2t_+,\\
\tau&=r-h^2,
&a_0&=\sqrt{\max\{0,1-1/(4\tau)\}},\\
a&=(1-u_3)/2+u_3a_0,
&b&=(1-u_4)\sqrt{(1-a)/2}+u_4a.
\end{aligned}
$$

The contact points are

$$
\begin{aligned}
P_1&=(\sqrt\tau,0,h,0),\\
P_2&=(\sqrt\tau(2a^2-1),\;2\sqrt\tau\,a\sqrt{1-a^2},\;h,0),\\
P_3&=(\sqrt\tau(2b^2-1),\;-2\sqrt\tau\,b\sqrt{1-b^2},\;h,0).
\end{aligned}
$$

The signed square $t|t|$ is used so that the height endpoint
$\sqrt{\max\{0,r-1/3\}}=t_+$ stays Lipschitz in $u_1$ at $r=1/3$.
The article proves that every required contact triangle occurs, up to an
isometry and relabelling, in this parameter cube. That geometric
surjectivity statement is not inferred from numerical sampling.

## 2. Outward interval arithmetic

Every interval is closed. Its endpoints are stored as integers on the grid
$2^{-44}\mathbb Z$. If $S=2^{44}$, the integer pair $(L,U)$ represents
$[L/S,U/S]$. Rational input constants are rounded outwards: the lower
integer is a floor and the upper integer is a ceiling.

- Addition and subtraction are exact on this grid.
- For multiplication, compute all four endpoint products, choose their
  minimum and maximum, and round the resulting bounds outwards to the grid.
- Division uses the four endpoint quotients and outward rounding; it is
  allowed only when the denominator interval excludes zero.
- For squaring, the lower endpoint is zero if the input interval contains
  zero; otherwise use the smaller squared endpoint. The upper endpoint is
  the larger squared endpoint. Round outwards.
- The positive-part operation encloses $\max\{0,x\}$ over the entire
  interval. It does not merely discard negative parameter values.
- Intersect an enclosure with another interval only when the latter is a
  proved bound for every actual value being enclosed. An empty intersection
  is rejected.

For a nonnegative integer $m$ representing $m/S$, let
$j=\operatorname{isqrt}(mS)$, where $\operatorname{isqrt}(N)$ is the
largest integer with square at most $N$. The lower square-root endpoint
is $j/S$. The upper endpoint is $j/S$ when $j^2=mS$, and $(j+1)/S$
otherwise. Thus square roots also use integer arithmetic only.

No floating-point tolerance is used for these operations. The
`Interval` class implements these rules; `require` raises an exception
when a condition fails. Certificate acceptance does not depend on Python
`assert` statements, which could be disabled by `python -O`.

## 3. Enclosing every member of a parameter box

A dyadic box is a product

$$
Q=\prod_{j=1}^4
\left[\frac{k_j}{2^{\ell_j}},\frac{k_j+1}{2^{\ell_j}}\right],
\qquad 0\le k_j<2^{\ell_j}.
$$

Its depth is $\ell_1+\cdots+\ell_4$. The function `enclose(Q)` uses
nonnegative integer levels and indices (the parser permits levels at most
128). It uses
the rules below to bound all parameter values in the box, not just its
centre or corners as independent samples.

**Squared radius and height.** Both $t$ and $r$ increase with $u_1$.
Evaluate the two endpoints of its interval by outward arithmetic, using

$$
t|t|=t_+^2-(-t)_+^2,
$$

and intersect the resulting radius enclosure with $[299/1000,r_0]$.
Both height endpoints $\sqrt{r/6}$ and $t_+$ increase with $u_1$.
Moreover, $t_+\le\sqrt{r/6}$, so $h$ decreases with $u_2$. Its
minimum is therefore bounded at the lower $u_1$ endpoint and upper
$u_2$ endpoint, and its maximum at the other two endpoints. This avoids
treating the two summands in the formula for $h$ as independent.

**Triangle radius.** Apply interval subtraction and squaring to
$\tau=r-h^2$, then intersect with the proved bounds

$$
\frac{299}{1200}\le\tau\le\frac13,
\qquad \frac56r\le\tau\le r.
$$

In the second pair, interval endpoints are used in the safe direction:
the lower endpoint of $(5/6)r$ bounds $\tau$ from below, and the upper
endpoint of $r$ bounds it from above.

**The variable $a$.** The quantity $a_0$ increases with $\tau$ and
belongs to $[0,1/2]$. Since $a_0\le1/2$, the interpolation defining
$a$ decreases with $u_3$. Bound $a$ using the appropriate endpoints
and intersect both $a_0$ and $a$ with $[0,1/2]$.

**The variable $b$.** Write $d_a=\sqrt{(1-a)/2}$. For
$0\le a\le1/2$, one has $d_a\ge a$, so $b$ decreases with $u_4$
for each actual $a$. If

$$
a\in[a_-,a_+],\quad d_a\in[d_-,d_+],\quad
u_4\in[w_-,w_+],
$$

then

$$
(1-w_+)d_-+w_+a_-\le b
\le(1-w_-)d_++w_-a_+.
$$

These are independent bounds for $a$ and $d_a$; no monotonicity of
$b$ in $a$ is assumed. Intersect the result with $[0,1]$.

**Contact coordinates.** The actual values $a,b$ lie in
$[0,1/\sqrt2]$. For $c=a$ or $b$, enclose $v=c^2$ and intersect
with $[0,1/2]$. The two functions

$$
2v-1,\qquad 2\sqrt{v(1-v)}
$$

increase on $[0,1/2]$. Evaluate their interval endpoints, multiply by
the enclosure for $\sqrt\tau$, and include the height enclosure.
This gives intervals $[L_{ij},U_{ij}]$ for every coordinate of every
$P_i$, and an interval $[r_Q^-,r_Q^+]$ for $r$.

All intersections in this procedure use analytic bounds valid on the
whole family. No bound is inferred from observed numerical values.

## 4. A common rational outer set

Replace each uncertain contact point by its rational midpoint

$$
A_i=\left(\frac{L_{ij}+U_{ij}}2\right)_{j=1}^4,
\qquad e_{ij}=\frac{U_{ij}-L_{ij}}2,
\qquad E_i=\frac23\sum_{j=1}^4e_{ij}.
$$

Every point $x$ of a true outer set satisfies
$|x_j|\le\|x\|\le\sqrt r\le\sqrt{2/5}<2/3$. Hence

$$
|\langle P_i-A_i,x\rangle|\le E_i.
$$

Expanding $\|x-P_i\|^2\le1$ and using $\|P_i\|^2=r$ gives
the common rational outer set

$$
\Omega_Q=\left\{x\in\mathbb E^4:
\|x\|^2\le r_Q^+,\quad
\|x\|^2-2\langle A_i,x\rangle\le b_i\ (i=1,2,3)\right\},
\qquad b_i=1-r_Q^-+2E_i.
$$

Thus every family member represented by $Q$ is contained in
$\Omega_Q$. Covering this outer set suffices for all parameters in $Q$.

## 5. The exact finite covering test

Let $C_1,\ldots,C_7\in\mathbb E^4$ be the proposed rational centres,
and let $s\in(0,\sigma]$ be their proposed squared covering radius.
Introduce a real variable $q$ and the polyhedron

$$
\mathcal P_Q=\left\{(x,q)\in\mathbb E^5:
\begin{aligned}
0&\le q\le r_Q^+,\\
q-2\langle A_i,x\rangle&\le b_i &&(i=1,2,3),\\
2\langle C_j,x\rangle-q&\le\|C_j\|^2-s &&(j=1,\ldots,7)
\end{aligned}\right\}.
$$

These are twelve linear inequalities. If a point $x\in\Omega_Q$ is
missed by all seven open balls $B(C_j,\sqrt s)$, then
$(x,\|x\|^2)\in\mathcal P_Q$.

The verifier first finds five affinely independent proposed centres
whose convex hull contains zero in its interior. It solves for their
strictly positive barycentric weights using exact rational arithmetic.
This proves boundedness of $\mathcal P_Q$: a recession direction
$(y,t)$ has $t=0$ from the bounds on $q$, and
$\langle C_j,y\rangle\le0$. The positive weighted sum of five such
inner products is zero, so all five vanish. The five selected centres
span $\mathbb E^4$, and therefore $y=0$.

Every vertex has five linearly independent active constraint normals.
The verifier enumerates all $\binom{12}{5}=792$ subsets of five
inequalities, solves each nonsingular linear system, and checks all
twelve inequalities at the resulting point. This finds every vertex,
even if more than five constraints are active or the polyhedron has
smaller dimension. A point may occur more than once; this is harmless.

For each feasible vertex it tests the strict rational inequality

$$
F(x,q):=\|x\|^2-q<0.
$$

A nonempty bounded polyhedron is the convex hull of its vertices, and
$F$ is convex. If it is strictly negative at every vertex, it is strictly
negative throughout the polyhedron. Therefore no point
$(x,\|x\|^2)$ arising from an uncovered point can lie in
$\mathcal P_Q$. If there are no feasible vertices, boundedness implies
that the polyhedron is empty, which also proves coverage.

Singularity, feasibility, positivity of barycentric weights, and strict
negativity are tested exactly. There are no floating-point tolerances in
this acceptance condition. Failure of this sufficient condition does not
prove that the proposed balls fail to cover the set.

## 6. Certificate format and completeness

The proof data are gzip-compressed JSON, with format identifier
`seven-balls-r4-single-v1` and parameterization identifier
`signed-square-cosine-v1`. The verifier checks these identifiers, the
fixed computational endpoints, the target squared radius, and the centre
denominator before reading the tree.

- `centers` is a list of seven-by-four integer arrays. Divide every
  coordinate by $2^{20}$ to obtain a centre configuration.
- `certificates` is a list of records. Each contains a centre-configuration
  index and an exact rational `radius_squared` string in $(0,\sigma]$.
- `tree` is a preorder encoding of the subdivision. A token in
  $\{0,1,2,3\}$ splits the indicated coordinate interval into its two
  equal halves. The left child is encoded before the right child. A
  negative token $-(j+1)$ is a terminal box using certificate $j$, with
  indexing beginning at zero. A `null` token is an unresolved box and is
  rejected in an ordinary complete verification run.

The verifier starts with the entire cube on a stack. Each split replaces
one pending box by both children. Each leaf removes one pending box.
The input must end exactly when the stack becomes empty. Missing children,
invalid tokens, and extra trailing tokens are rejected. This proves that
the entire cube is represented; the exact terminal-volume sum of one is
an additional check. Closed boxes and closed interval enclosures include
shared boundaries.

For every leaf, the verifier rebuilds the box, all interval enclosures,
the rational outer set, and the exact covering test. Stored numerical
scores and previous verification reports play no role.

The complete data have 23,494 leaves, 23,493 splits, maximum depth 29,
6,913 centre configurations, and 19,656 centre/radius pairs. The largest
stored squared radius is $\sigma$; it is used on 28 leaves. All other
leaves use a strictly smaller squared radius.

## 7. Analytic checks at the interfaces

The exact replay also checks the scalar constants used in the article's
two analytic ranges. These calculations guard the interfaces with the
computational interval; they do not replace the geometric arguments.

For the small range $r\le299/1000$, seven balls cover the entire
circumball with squared radius at most $5r/6\le299/1200$.

For the large range $r_0\le r\le2/5$, the article uses the five
contracted contact points $\alpha p_i$, where $\alpha=9/20$, and the
two centres $\pm\sqrt{r/2}\,e$. Here $e$ minimizes the weighted sum
of squared contact projections. The two truncated cones of half-angle
$\pi/4$ are covered with squared radius at most $r/2\le1/5$.
Outside those cones, the projection estimate reduces coverage to

$$
Q_s(q,r)=-q^2+B_s(r)q+C_s(r)<0,
$$

where

$$
B_s(r)=1+s+2\alpha-(1+\alpha^2+6\alpha)r,\qquad
C_s(r)=(\alpha^2r-s)(1-r).
$$

The maximum over real $q$ is
$\Delta_s(r)=B_s(r)^2/4+C_s(r)$. Since

$$
\Delta_s''(r)=\frac{2307121}{320000}>0,
$$

it suffices to check

$$
s>\frac{299}{1200},\qquad s>\frac15,\qquad
\Delta_s(r_0)<0,\qquad \Delta_s'(2/5)<0.
$$

The second inequality follows from the first but is retained in the
code as an explicit check of the cone covering. At $s=\sigma$, the
checks include

$$
\begin{aligned}
\Delta_\sigma(r_0)
&=-\frac{990904111119}{640000000000000000},\\
\Delta_\sigma'(2/5)
&=-\frac{274796380761}{320000000000},\\
\sigma-\frac15&=\frac{19980001}{400000000}.
\end{aligned}
$$

The revised verifier uses these projection-and-cone checks directly.
It contains no additional condition involving the older fixed axial
centres $\pm(2/5)e$.

The routine `analytic_radius_bound` performs 32 rational bisection steps
between $249/1000$ and $\sigma$, retaining only a squared radius that
passes all checks. It returns

$$
s_{\rm an}=\frac{107350418082054063}{429496729600000000},
\qquad \sqrt{s_{\rm an}}<0.4999446664.
$$

The uniform all-$r$ squared-radius bound is the larger of $s_{\rm an}$
and the largest verified leaf value $s_Q$. For this certificate it is
still $\sigma$. The smaller analytic bound alone therefore does not
improve the theorem's uniform numerical constant.

## 8. Reports, strictness, and numerical discovery

`verify_report.json` records the number of checked leaves, exact covered
and unresolved volumes, failures, largest feasible vertex value, analytic
constant checks, exact squared-radius bounds, upward-rounded decimal
radius bounds, code and data hashes, backend, worker count, and elapsed
times. The decimal bounds are computed using integer square roots.

`complete_computational_certificate` is true only after a full successful
replay with no failures and complete coverage of the root cube. Passing
`--sample` never yields a complete certificate. `--allow-partial` permits
diagnosis of a tree with unresolved leaves, but such a tree cannot yield
a complete certificate. An exception or failed check is not success.
Exit code zero alone is insufficient: a successful sampled or partial
diagnostic run also exits with code zero. The completeness field and the
reported coverage and failure fields must all be checked.

The test is strict at every feasible vertex. Thus a verified leaf gives
coverage by open balls of its stored radius. Together with the analytic
arguments, the theorem covers each compact normalized set by seven open
balls of radius $\rho$. The continuous function
$x\mapsto\min_j\|x-C_j\|$ then has maximum strictly below $\rho$ on
that compact set. Taking the closure first handles arbitrary sets of
diameter one. This proves the common-diameter bound strictly below
$2\rho=0.9999$. The report does not claim an explicitly quantified
uniform improvement below its displayed upper bound.

The largest vertex value measures slack in $F(x,q)$ at the individual
stored radii. It is not directly a permissible decrease of the common
covering radius.

The finder uses floating-point samples and optimization to propose
centres, choose subdivisions, and suggest smaller radii. Every accepted
proposal subsequently passes the exact test. After a box passes at
$\sigma$, optional radius polishing proposes a smaller squared radius,
rounds it upwards to a multiple of $10^{-7}$ and adds two grid units,
then reruns the exact test. Unsuccessful proposals are moved towards the
already verified radius, with at most three exact attempts. If none
succeeds, the original certificate is retained.

The full numerical objectives, schedules, seed construction, and search
limits are recorded in `SEARCH.md`.

The seed bank, checkpoint, random samples, optimization logs, and
`radius_by_box.csv.gz` help discovery or inspection; they are not proof
inputs for replay. Numerical search completion must be followed by a
separate full verification. The standard-library and GMP backends use the
same verifier algorithm, so agreement between them is not an independent
geometric proof.

## 9. Changes in this supplementary version

- The technical interval rules and software instructions are collected
  here instead of in long article appendices.
- The mathematical notation is $\rho$ for the target radius,
  $\sigma=\rho^2$ for its square, and $r_0$ for the upper computational
  endpoint. The existing code constant names are mapped explicitly above.
- The analytic checks now match the scalar-projection and $\pi/4$-cone
  proof. The obsolete fixed-axis inequalities were removed and the cone
  squared-radius gap was added.
- The parameterization, exact interval arithmetic, finite covering test,
  tree parser, centre coordinates, leaf radii, certificate bytes, and seed
  bytes are unchanged. The source-code hash changes because its analytic
  checks and documentation were updated.
- `verify_report.json` is a new full replay of the revised source. Consult
  its `wall_seconds`, `workers`, and `backend` fields for the recorded run;
  these are not a runtime guarantee for other machines.

Run `sha256sum -c SHA256SUMS` on systems providing `sha256sum` to check the
packaged file digests. The repository URL above is intentionally a
placeholder and must be replaced before publication.
