#!/usr/bin/env python3
"""Seven balls in E^4: one parameter cube, with exact radius certificates.

Python 3.10+. Exact verification needs only the standard library; optional
gmpy2 speeds up rational arithmetic. Only the finder needs NumPy and SciPy.

  python seven_balls_single.py verify --data search_certificate.json.gz --workers 12 --report verify_report.json
  python -m pip install numpy scipy gmpy2  # optional finder dependencies
  python seven_balls_single.py find --out single_run --hours 24 --workers 16 --backend gmp
  python seven_balls_single.py find --out single_run --hours 24 --workers 16 --backend gmp --resume
  python seven_balls_single.py verify --data single_run/search_certificate.json.gz --workers 16 --backend gmp --report single_run/verify_report.json

The four-dimensional root cube represents ALL .299 <= r <= .38401. Here r
is SQUARED enclosing radius. No old subdivision is imported. By default
the old successful centre configurations are reused as numerical starting
points only; --fresh-seeds restricts this to the 28 exploratory seeds.

Each accepted leaf first passes at radius .49995. We then numerically
estimate a smaller sufficient radius for those centres and check the
proposal EXACTLY. The stored rational squared radius is always certified.
It is an upper bound on the true covering radius, not a claim of optimality.
--no-radius-polish skips this optional improvement. Sampled radii are
explicitly labelled heuristic and are never used to accept a certificate.

The SQLite checkpoint is resumable with identical settings. Ctrl+C or a
time limit preserves unresolved boxes. A partial search is NOT a proof.
status.json and progress.jsonl give covered parameter volume, unresolved
counts, and the largest radius among covered boxes. The final export also
contains radius_by_box.csv.gz. Send search_certificate.json.gz and the
full verify_report.json to update the manuscript after a successful run.

PROOF DATA: gzip JSON, format seven-balls-r4-single-v1. Constants fixing
the geometry, centre grid, and target are checked. `centers` lists 7-by-4
integer matrices on 2**(-20) Z. `certificates` lists pairs of a centre index
and an EXACT rational radius_squared. `tree` is a preorder list: 0..3 split
that coordinate, -(j+1) is a leaf using certificates[j], null is unresolved.
Both children of every split must appear. All boxes and interval bounds
are rebuilt by the verifier, and every leaf radius is rechecked. Stored
scores, sample values, timestamps, and previous reports are not trusted.

The target radius is rho=9999/20000 and its square is sigma=rho**2.
Strict rational vertex inequalities prove coverage by OPEN balls. Together
with the manuscript's geometric reductions and compactness, a successful
complete replay gives seven balls of common diameter strictly below .9999
for every set of diameter one. The displayed decimal radius is an upper
bound, not a numerical measurement of that strict improvement.

Parts I--IV are proof-critical. Floating-point search is confined to Part V.
README.md gives the commands; VERIFICATION.md gives all interval rules,
the certificate format, and the precise scope of verification.
"""

from __future__ import annotations
import argparse
import csv
from concurrent.futures import ProcessPoolExecutor, wait, FIRST_COMPLETED
from contextlib import contextmanager
from fractions import Fraction as F
import gzip
import hashlib
from itertools import combinations, permutations
import json
from math import isqrt, sqrt
import multiprocessing as mp
import os
from pathlib import Path
import random
import sqlite3
import sys
import time

FORMAT = 'seven-balls-r4-single-v1'
PARAMETERIZATION = 'signed-square-cosine-v1'
SEED_FORMAT = 'seven-balls-r4-single-seeds-v1'
RADIUS = F(9999, 20000)
RADIUS_SQUARED = RADIUS * RADIUS
# Manuscript notation: RADIUS = rho; RADIUS_SQUARED = sigma; R_HI = r_0.
CENTER_BITS = 20
R_LO, R_HI = F(299, 1000), F(38401, 100000)
DELTA_LO, DELTA_HI = F(1, 3)-R_LO, R_HI-F(1, 3)
ROOT = ((0, 0),) * 4
BASES = tuple(combinations(range(12), 5))
DEFAULT_SEEDS = Path(__file__).with_name('seven_balls_seeds.json.gz')


# ======================================================================
# PART I. INTEGER INTERVAL ARITHMETIC AND THE SINGLE PARAMETERIZATION
# ======================================================================

def require(condition, message):
    """Never use Python assert for certificate acceptance (python -O removes it)."""
    if not condition:
        raise ValueError(message)


def ceil_div(a, b):
    """Integer ceiling of a/b; used below only with b > 0."""
    return -((-a) // b)


class Interval:
    """Closed interval with integer endpoints on the grid 2**(-44) Z.

    All intermediate operations are rounded OUTWARD. For multiplication and
    division, the extrema occur among the four endpoint pairs. Division is
    allowed only when the denominator interval excludes zero. The method
    positive() encloses max(0,x); it does not merely discard negative values.
    """
    SCALE = 1 << 44
    __slots__ = ('lo', 'hi')

    def __init__(self, lo, hi=None):
        a = F(lo)
        b = a if hi is None else F(hi)
        self.lo = a.numerator * self.SCALE // a.denominator
        self.hi = ceil_div(b.numerator * self.SCALE, b.denominator)
        require(self.lo <= self.hi, 'Reversed interval')

    @classmethod
    def raw(cls, lo, hi):
        require(type(lo) is int and type(hi) is int and lo <= hi,
                'Invalid fixed-point interval')
        obj = object.__new__(cls)
        obj.lo, obj.hi = lo, hi
        return obj

    @classmethod
    def coerce(cls, x):
        return x if isinstance(x, cls) else cls(x)

    def __add__(self, other):
        other = self.coerce(other)
        return self.raw(self.lo + other.lo, self.hi + other.hi)
    __radd__ = __add__

    def __neg__(self):
        return self.raw(-self.hi, -self.lo)

    def __sub__(self, other):
        return self + (-self.coerce(other))

    def __rsub__(self, other):
        return self.coerce(other) + (-self)

    def __mul__(self, other):
        other = self.coerce(other)
        values = [a*b for a in (self.lo, self.hi)
                  for b in (other.lo, other.hi)]
        return self.raw(min(values)//self.SCALE,
                        ceil_div(max(values), self.SCALE))
    __rmul__ = __mul__

    def __truediv__(self, other):
        other = self.coerce(other)
        if other.lo <= 0 <= other.hi:
            raise ZeroDivisionError('Interval denominator contains zero')
        if other.hi < 0:
            return (-self)/(-other)
        lows = [a*self.SCALE//b for a in (self.lo, self.hi)
                for b in (other.lo, other.hi)]
        highs = [ceil_div(a*self.SCALE, b) for a in (self.lo, self.hi)
                 for b in (other.lo, other.hi)]
        return self.raw(min(lows), max(highs))

    def __rtruediv__(self, other):
        return self.coerce(other)/self

    def square(self):
        lo = 0 if self.lo <= 0 <= self.hi else min(self.lo**2, self.hi**2)
        return self.raw(lo//self.SCALE,
                        ceil_div(max(self.lo**2, self.hi**2), self.SCALE))

    def sqrt(self):
        require(self.lo >= 0, 'Negative square-root argument')
        lo = isqrt(self.lo*self.SCALE)
        square = self.hi*self.SCALE
        hi = isqrt(square)
        return self.raw(lo, hi if hi*hi == square else hi+1)

    def positive(self):
        return self.raw(max(0, self.lo), max(0, self.hi))

    def intersect(self, other):
        other = self.coerce(other)
        return self.raw(max(self.lo, other.lo), min(self.hi, other.hi))

    def left(self):
        return self.raw(self.lo, self.lo)

    def right(self):
        return self.raw(self.hi, self.hi)

    def endpoints(self):
        return F(self.lo, self.SCALE), F(self.hi, self.SCALE)


def axis_bounds(axis):
    k, level = axis
    require(type(k) is int and type(level) is int
            and 0 <= level <= 128 and 0 <= k < 2**level,
            'Invalid dyadic coordinate')
    return F(k, 2**level), F(k+1, 2**level)


def children(box, dim):
    require(type(dim) is int and 0 <= dim < 4, 'Invalid split dimension')
    left, right = list(box), list(box)
    k, level = box[dim]
    left[dim], right[dim] = (2*k, level+1), (2*k+1, level+1)
    return tuple(left), tuple(right)


def depth(box):
    return sum(level for _, level in box)


def radius_decimal(s, digits=10):
    """Decimal UPPER bound on sqrt(s), computed with integers, never floats."""
    s = F(s)
    require(s >= 0, 'Negative squared radius')
    scale = 10**digits
    m = isqrt(s.numerator*scale*scale//s.denominator)
    if m*m*s.denominator < s.numerator*scale*scale:
        m += 1
    return f'{m//scale}.{m%scale:0{digits}d}'


def parameter_r_interval(axis):
    """r=1/3+t|t|, t=x sqrt(DELTA_HI)-(1-x)sqrt(DELTA_LO).

    t and r increase with x. Evaluating the two endpoints with outward
    intervals therefore bounds the entire interval, including t=0.
    The identity t|t|=t_+**2-(-t)_+**2 needs no division or sign decision.
    """
    I = Interval
    low, high = map(I, axis_bounds(axis))
    dl, du = I(DELTA_LO).sqrt(), I(DELTA_HI).sqrt()

    def at(x):
        t = x*du-(1-x)*dl
        r = (I(F(1, 3))+t.positive().square()-(-t).positive().square())
        return t, r.intersect(I(R_LO, R_HI))

    tl, rl = at(low)
    tu, ru = at(high)
    return I.raw(rl.lo, ru.hi), (tl, rl), (tu, ru)


def rational_outer_model(r_interval, contact_intervals):
    """Replace uncertain contacts by exact rational centers and error terms.

    For each true contact P_i, A_i is the coordinatewise midpoint and
    e_ij the half-width. Since |X_j| < 2/3 throughout every true body,
    |(P_i-A_i).X| <= E_i=(2/3)*sum_j e_ij. Therefore
        |X|^2 <= r_upper,  |X|^2 - 2 A_i.X <= 1-r_lower+2E_i.
    This outer approximation is common to ALL parameters in the box.
    """
    contacts = [[(lo+hi)/2 for lo, hi in row] for row in contact_intervals]
    errors = [F(2, 3)*sum((hi-lo)/2 for lo, hi in row)
              for row in contact_intervals]
    rl, ru = r_interval
    return {'r_upper': ru, 'contacts': contacts,
            'intercepts': [1-rl+2*e for e in errors],
            'uncertainty': ru-rl+max(errors)}


def enclose(box):
    """One exact enclosure for the WHOLE computational interval.

    With (x,y,z,w) in [0,1]^4, set
      t=x sqrt(r_hi-1/3)-(1-x)sqrt(1/3-r_lo), r=1/3+t|t|,
      h=(1-y)sqrt(r/6)+y max(t,0), tau=r-h*h,
      a0=sqrt(max(0,1-1/(4tau))), a=(1-z)/2+z*a0,
      b=(1-w)sqrt((1-a)/2)+w*a.
    Then 5r/6 <= tau <= min(r,1/3). The contact coordinates are
      (sqrt(tau),0,h,0),
      (sqrt(tau)(2a*a-1), 2sqrt(tau)a sqrt(1-a*a),h,0),
      (sqrt(tau)(2b*b-1),-2sqrt(tau)b sqrt(1-b*b),h,0).

    Height increases with x and decreases with y: both height endpoints
    increase with x, and max(t,0) <= sqrt(r/6) since r <= 2/5. Thus corner
    bounds avoid unnecessary dependency loss in the height interpolation.
    All remaining intersections below use proved analytic bounds, not
    numerical observations. Fixed-point roundoff is always OUTWARD.
    """
    require(len(box) == 4, 'A parameter box must have four coordinates')
    I = Interval
    r, (tl, rl), (tu, ru) = parameter_r_interval(box[0])
    y, z, w = [I(*axis_bounds(axis)) for axis in box[1:]]
    hl = (1-y.right())*(rl/6).sqrt()+y.right()*tl.positive()
    hu = (1-y.left())*(ru/6).sqrt()+y.left()*tu.positive()
    h = I.raw(hl.lo, hu.hi).positive()
    tau = (r-h.square()).intersect(I(F(299, 1200), F(1, 3)))
    tau = I.raw(max(tau.lo, (r*F(5, 6)).lo), min(tau.hi, r.hi))
    a0 = (1-1/(4*tau)).positive().sqrt().intersect(I(0, F(1, 2)))
    # For each actual tau, a decreases with z. a0 increases with tau.
    a = I.raw(((1-z.right())*F(1, 2)+z.right()*a0.left()).lo,
              ((1-z.left())*F(1, 2)+z.left()*a0.right()).hi)
    a = a.intersect(I(0, F(1, 2)))
    base = ((1-a)/2).sqrt()
    # base >= a for 0 <= a <= 1/2, so b decreases with w. Do NOT assume
    # b is monotone in a: independent endpoint bounds are used here.
    b = I.raw(((1-w.right())*base.left()+w.right()*a.left()).lo,
              ((1-w.left())*base.right()+w.left()*a.right()).hi)
    b = b.intersect(I(0, 1))

    def direction(c):
        # Both a and b lie in [0,1/sqrt(2)]. On that interval, 2c*c-1
        # and 2c sqrt(1-c*c) increase. Evaluate endpoint squares to avoid
        # treating c and sqrt(1-c*c) as unrelated interval variables.
        z = c.square().intersect(I(0, F(1, 2)))
        sl = (4*z.left()*(1-z.left())).positive().sqrt()
        su = (4*z.right()*(1-z.right())).positive().sqrt()
        return 2*z-1, I.raw(sl.lo, su.hi)

    ca, sa = direction(a)
    cb, sb = direction(b)
    v, zero = tau.sqrt(), I(0)
    raw = [[v, zero, h, zero], [v*ca, v*sa, h, zero],
           [v*cb, -v*sb, h, zero]]
    return rational_outer_model(r.endpoints(),
                               [[value.endpoints() for value in row] for row in raw])


# ======================================================================
# PART II. EXACT FINITE COVERING TEST IN FIVE VARIABLES
# ======================================================================

def rational_type(backend):
    require(backend in ('fraction', 'gmp', 'auto'), 'Unknown exact backend')
    if backend != 'fraction':
        try:
            from gmpy2 import mpq
            return mpq, 'gmp'
        except ImportError:
            if backend == 'gmp':
                raise RuntimeError('Install gmpy2, or use --backend fraction') from None
    return F, 'fraction'


def dot(a, b):
    return sum(x*y for x, y in zip(a, b))


def solve(matrix, rhs):
    """Exact Gauss--Jordan elimination. Return None iff the matrix is singular.

    No tolerance or pivot threshold occurs here. Every nonzero rational is
    a valid pivot; equality, feasibility, and strict negativity are exact.
    """
    n = len(rhs)
    a = [list(row)+[value] for row, value in zip(matrix, rhs)]
    for column in range(n):
        pivot = next((j for j in range(column, n) if a[j][column]), None)
        if pivot is None:
            return None
        a[column], a[pivot] = a[pivot], a[column]
        divisor = a[column][column]
        a[column] = [x/divisor for x in a[column]]
        for j in range(n):
            if j != column and a[j][column]:
                multiplier = a[j][column]
                a[j] = [x-multiplier*y for x, y in zip(a[j], a[column])]
    return tuple(row[-1] for row in a)


def check_centers(integers):
    require(isinstance(integers, list) and len(integers) == 7
            and all(isinstance(row, list) and len(row) == 4 for row in integers)
            and all(type(v) is int for row in integers for v in row),
            'Centers must be a 7 by 4 integer array')


def exact_cover(model, integers, backend='fraction', early_failure=False,
                radius_squared=RADIUS_SQUARED):
    """Certify coverage of a whole rational outer body by seven OPEN balls.

    A missed point X produces a point (X,q), q=|X|^2, in a rational
    polyhedron P with 12 inequalities. The test proves P lies strictly
    above the paraboloid q=|X|^2, by maximizing F(X,q)=|X|^2-q at vertices.

    Boundedness is checked FIRST. Five centers containing zero strictly
    inside their full-dimensional simplex rule out a nonzero recession
    direction. This is needed before treating an empty vertex list as an
    empty polyhedron. A nonempty bounded polyhedron always has a vertex.
    The simplex test is sufficient, not necessary: failure is inconclusive.

    At any vertex, active constraint normals span R^5, even if P itself
    has lower dimension. Enumerating all 792 five-row subsets therefore
    finds every vertex, with harmless duplicates in degenerate cases.
    A convex function on a polytope is bounded above by its largest vertex
    value. Strict negativity excludes every lifted point q=|X|^2.
    """
    start = time.process_time()
    check_centers(integers)
    Q, used_backend = rational_type(backend)
    convert = lambda x: Q(x.numerator, x.denominator) if isinstance(x, F) else Q(x)
    r = convert(model['r_upper'])
    radius_squared = F(radius_squared)
    require(0 < radius_squared <= RADIUS_SQUARED, 'Invalid certified squared radius')
    s = Q(radius_squared.numerator, radius_squared.denominator)
    p = [[convert(v) for v in row] for row in model['contacts']]
    rhs_contacts = [convert(v) for v in model['intercepts']]
    c = [[Q(v, 1 << CENTER_BITS) for v in row] for row in integers]
    simplex = None
    for indices in combinations(range(7), 5):
        matrix = [[c[j][i] for j in indices] for i in range(4)] + [[Q(1)]*5]
        weights = solve(matrix, [Q(0)]*4+[Q(1)])
        if weights is not None and min(weights) > 0:
            simplex = list(indices)
            break
    if simplex is None:
        return {'ok': False, 'reason': 'no_positive_simplex',
                'backend': used_backend, 'cpu_seconds': time.process_time()-start}
    a = [[-2*v for v in row]+[Q(1)] for row in p]
    a += [[2*v for v in row]+[Q(-1)] for row in c]
    a += [[Q(0)]*4+[Q(1)], [Q(0)]*4+[Q(-1)]]
    b = rhs_contacts+[dot(row, row)-s for row in c]+[r, Q(0)]
    maximum = None
    feasible_bases = 0
    for indices in BASES:
        vertex = solve([a[i] for i in indices], [b[i] for i in indices])
        if vertex is None or any(dot(row, vertex) > value for row, value in zip(a, b)):
            continue
        feasible_bases += 1
        value = dot(vertex[:4], vertex[:4])-vertex[4]
        if maximum is None or value > maximum:
            maximum = value
        if early_failure and value >= 0:
            break
    ok = maximum is None or maximum < 0
    return {'ok': ok, 'reason': 'empty_polytope' if maximum is None else
            ('negative_vertex_maximum' if ok else 'nonnegative_vertex'),
            'maximum': None if maximum is None else str(maximum),
            'feasible_bases': feasible_bases, 'positive_simplex': simplex,
            'backend': used_backend, 'cpu_seconds': time.process_time()-start}


# ======================================================================
# PART III. DATA AND COMPLETE PARAMETER-CUBE COVERAGE
# ======================================================================

def read_gzip_json(path):
    def invalid_constant(value):
        raise ValueError('Non-finite JSON value: '+value)
    with gzip.open(path, 'rt', encoding='utf-8') as stream:
        return json.load(stream, parse_constant=invalid_constant)


def load_data(path):
    data = read_gzip_json(path)
    require(data.get('format') == FORMAT, 'Unknown certificate format')
    require(data.get('parameterization') == PARAMETERIZATION,
            'Unknown parameterization')
    require(data.get('center_bits') == CENTER_BITS, 'Unexpected center denominator')
    require(data.get('target_radius_squared') == str(RADIUS_SQUARED),
            'Unexpected target radius')
    require(data.get('r_interval') == [str(R_LO), str(R_HI)],
            'Unexpected computational range')
    require(isinstance(data.get('centers'), list), 'Missing center list')
    for centers in data['centers']:
        check_centers(centers)
    require(isinstance(data.get('certificates'), list), 'Missing certificate list')
    for row in data['certificates']:
        require(isinstance(row, dict) and type(row.get('center')) is int
                and 0 <= row['center'] < len(data['centers']), 'Invalid center index')
        require(type(row.get('radius_squared')) is str, 'Radius must be a rational string')
        require(0 < F(row['radius_squared']) <= RADIUS_SQUARED, 'Invalid leaf radius')
    return data


def tree_leaves(tree, certificate_count, allow_pending=False):
    """Parse a complete preorder tree, reconstructing every dyadic box.

    Each pending stack item is one part of the root that still must appear.
    A split replaces exactly one part by BOTH of its exact halves. Reaching
    an empty stack exactly at the end proves no region was omitted. The
    rational volume check is supplementary, not the coverage argument.
    The returned leaf indices are preorder token positions (zero based).
    """
    require(isinstance(tree, list), 'Missing preorder tree')
    stack = [ROOT]
    leaves = []
    covered_volume = F(0)
    pending_volume = F(0)
    maximum_depth = 0
    splits = 0
    for position, token in enumerate(tree):
        require(stack, 'Trailing tokens after the root subtree')
        box = stack.pop()
        for axis in box:
            axis_bounds(axis)
        maximum_depth = max(maximum_depth, depth(box))
        if token is None:
            require(allow_pending, 'Unresolved leaf in certificate')
            pending_volume += F(1, 2**depth(box))
        elif type(token) is int and 0 <= token < 4:
            left, right = children(box, token)
            stack.extend([right, left])  # left is parsed first
            splits += 1
        else:
            require(type(token) is int and token < 0
                    and 0 <= -token-1 < certificate_count, 'Invalid leaf certificate index')
            leaves.append((position, box, -token-1))
            covered_volume += F(1, 2**depth(box))
    require(not stack, 'Truncated tree: one or more children are missing')
    require(covered_volume+pending_volume == 1, 'Root volume is not one')
    return leaves, {'leaves': len(leaves), 'splits': splits,
                    'max_depth': maximum_depth, 'covered_volume': str(covered_volume),
                    'unresolved_volume': str(pending_volume)}


def analytic_constants(radius_squared=RADIUS_SQUARED):
    """Check the rational scalar inequalities in the two analytic ranges.

    This does not replace their geometric proofs in the manuscript. It
    guards the exact constants used at the interfaces with the computation.
    The large-radius proof uses contracted contacts and two 45-degree cones.
    The latter have squared covering radius at most 1/5; no fixed axial
    centre parameter or eigenvalue calculation is used here.
    """
    alpha, r0 = F(9, 20), R_HI
    s = F(radius_squared)
    slope = 1+alpha*alpha+6*alpha
    B = lambda r: 1+s+2*alpha-slope*r
    C = lambda r: (alpha*alpha*r-s)*(1-r)
    delta = lambda r: B(r)*B(r)/4+C(r)
    derivative_right = -slope*B(F(2, 5))/2 + alpha*alpha*(1-2*F(2, 5))+s
    second_derivative = slope*slope/2-2*alpha*alpha
    values = {'small_range_gap': s-F(299, 1200), 'Delta_at_r0': delta(r0),
              'Delta_prime_at_2_over_5': derivative_right,
              'Delta_second_derivative': second_derivative,
              'cone_squared_radius_gap': s-F(1, 5)}
    require(values['small_range_gap'] > 0 and values['Delta_at_r0'] < 0
            and derivative_right < 0 and second_derivative > 0
            and values['cone_squared_radius_gap'] > 0, 'Analytic constant check failed')
    return {key: str(value) for key, value in values.items()}


# ======================================================================
# PART IV. FULL EXACT REPLAY (spawn-safe on Windows)
# ======================================================================

def analytic_radius_bound():
    """A rigorous bound for BOTH analytic ranges, with their fixed endpoints.

    The manuscript's target is sigma=rho**2. A candidate squared radius s
    remains valid when the small-ball, 45-degree-cone, and quadratic sign
    checks in analytic_constants(s) pass. Only a passing rational value
    is retained. This is a
    bound for those explicit constructions, not an optimal covering radius.
    """
    low, high = F(249, 1000), RADIUS_SQUARED
    analytic_constants(high)
    for _ in range(32):
        middle = (low+high)/2
        try:
            analytic_constants(middle)
        except ValueError:
            low = middle
        else:
            high = middle
    return high


def verify_job(job):
    position, box, integers, radius_squared, backend = job
    result = exact_cover(enclose(box), integers, backend,
                         radius_squared=F(radius_squared))
    return {'leaf': position, 'radius_squared': radius_squared, **result}


def write_json(path, obj):
    path = Path(path)
    temp = path.with_name(path.name+'.tmp')
    temp.write_text(json.dumps(obj, indent=2, allow_nan=False)+'\n', encoding='utf-8')
    temp.replace(path)


def verify_command(args):
    start = time.monotonic()
    data = load_data(args.data)
    _, backend = rational_type(args.backend)
    leaves, structure = tree_leaves(data['tree'], len(data['certificates']),
                                   allow_pending=args.allow_partial)
    if args.sample is not None and args.sample < len(leaves):
        leaves = sorted(random.Random(20261014).sample(leaves, args.sample))
    jobs = []
    for pos, box, index in leaves:
        certificate = data['certificates'][index]
        jobs.append((pos, box, data['centers'][certificate['center']],
                     certificate['radius_squared'], backend))
    failures, checked, cpu, largest, radius_max = [], 0, 0., None, None

    def record(result):
        nonlocal checked, cpu, largest, radius_max
        checked += 1
        cpu += result['cpu_seconds']
        if not result['ok']:
            failures.append(result)
        else:
            s = F(result['radius_squared'])
            radius_max = s if radius_max is None else max(radius_max, s)
        if result.get('maximum') is not None:
            value = F(result['maximum'])
            largest = value if largest is None else max(largest, value)

    if args.workers == 1:
        next_report = time.monotonic()+30
        for job in jobs:
            record(verify_job(job))
            if time.monotonic() >= next_report:
                print(f'Checked {checked}/{len(jobs)} leaves', flush=True)
                next_report = time.monotonic()+30
    else:
        with ProcessPoolExecutor(args.workers, mp_context=mp.get_context('spawn')) as pool:
            batch = max(32, args.workers*16)
            for first in range(0, len(jobs), batch):
                for result in pool.map(verify_job, jobs[first:first+batch], chunksize=4):
                    record(result)
                print(f'Checked {checked}/{len(jobs)} leaves', flush=True)
    complete = (not failures and args.sample is None
                and structure['covered_volume'] == '1')
    analytic_s = analytic_radius_bound()
    uniform_s = max(radius_max, analytic_s) if complete else None
    report = {
        'format': FORMAT, 'parameterization': PARAMETERIZATION,
        'target_radius': str(RADIUS), 'target_radius_squared': str(RADIUS_SQUARED),
        'r_interval': [str(R_LO), str(R_HI)], 'backend': backend,
        'workers': args.workers, 'sample_requested': args.sample,
        'allow_partial': args.allow_partial, **structure,
        'checked': checked, 'failures': failures,
        'complete_computational_certificate': complete,
        'data_sha256': hashlib.sha256(Path(args.data).read_bytes()).hexdigest(),
        'code_sha256': hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
        'largest_vertex_value': None if largest is None else str(largest),
        'largest_checked_radius_squared': None if radius_max is None else str(radius_max),
        'largest_checked_radius_upper_decimal': None if radius_max is None else radius_decimal(radius_max),
        'analytic_radius_squared_bound': str(analytic_s),
        'analytic_radius_upper_decimal': radius_decimal(analytic_s),
        'analytic_constant_checks': analytic_constants(analytic_s),
        'uniform_all_r_squared_radius_bound': None if uniform_s is None else str(uniform_s),
        'uniform_all_r_radius_upper_decimal': None if uniform_s is None else radius_decimal(uniform_s),
        'cpu_seconds': cpu, 'wall_seconds': time.monotonic()-start,
        'scope': 'Every accepted leaf is checked at its own exact rational squared radius. '
                 'The uniform all-r bound is emitted only after a complete successful replay. '
                 'It also uses the manuscript geometric reductions and analytic proofs, '
                 'with the scalar checks repeated at the reported analytic radius. '
                 'A sample or unresolved tree never certifies the whole parameter cube.'}
    if args.report:
        write_json(args.report, report)
    print(json.dumps(report, indent=2), flush=True)
    return 0 if not failures else 1


# ======================================================================
# PART V. NUMERICAL SEARCH -- NOTHING BELOW IS USED TO ACCEPT A PROOF
# (The finder calls exact_cover above before recording a covered box.)
# ======================================================================

SEARCH_DEFAULTS = {'search_seconds': 8.0, 'search_uncertainty': 0.035,
                   'starts': 2, 'refinements': 3, 'sample_power': 11,
                   'seed_candidates': 4, 'cache_size': 128,
                   'float_accept_margin': 1e-7, 'seed': 20261011}
SEARCH_DEADLINE = float('inf')
SEARCH_ENGINE = None


class SearchTimeout(Exception):
    pass


def search_budget(*unused):
    if time.monotonic() >= SEARCH_DEADLINE:
        raise SearchTimeout()


def numerical_imports():
    """Imported only by find; constrain BLAS to one thread in EACH worker."""
    for name in ('OPENBLAS_NUM_THREADS', 'MKL_NUM_THREADS', 'OMP_NUM_THREADS',
                 'NUMEXPR_NUM_THREADS', 'VECLIB_MAXIMUM_THREADS'):
        os.environ[name] = '1'
    global np, minimize, ndtri, logsumexp, qmc, FLOAT_BASES
    import numpy as np
    from scipy.optimize import minimize
    from scipy.special import ndtri, logsumexp
    from scipy.stats import qmc
    FLOAT_BASES = np.asarray(BASES, dtype=int)


def numerical_point(parameters):
    x, y, z, w = map(float, parameters)
    t = x*sqrt(float(DELTA_HI))-(1-x)*sqrt(float(DELTA_LO))
    r = 1/3+t*abs(t)
    h = (1-y)*sqrt(r/6)+y*max(t, 0.)
    tau = r-h*h
    a = (1-z)/2+z*sqrt(max(0, 1-1/(4*tau)))
    a = min(.5, max(0., a))
    b = (1-w)*sqrt((1-a)/2)+w*a
    v = sqrt(tau)
    return np.array([[v, 0, h, 0],
                     [v*(2*a*a-1), 2*v*a*sqrt(max(0, 1-a*a)), h, 0],
                     [v*(2*b*b-1), -2*v*b*sqrt(max(0, 1-b*b)), h, 0]])


def midpoint_contacts(box):
    return numerical_point([(a+b)/2 for a, b in map(axis_bounds, box)])


def distances(points, centers):
    return np.maximum(0, np.sum(points*points, axis=1)[:, None]
                      +np.sum(centers*centers, axis=1)[None, :]-2*points@centers.T)


def minimum_ball(points, old_center):
    """Approximate minimum enclosing ball of one assigned sample cluster.

    SLSQP sees only a small active pool, enlarged by farthest sample points.
    Its answer is a proposal, so numerical tolerances do not affect validity.
    """
    search_budget()
    if not len(points):
        return old_center
    squared = np.sum((points-old_center)**2, axis=1)
    pool = np.argsort(squared)[-min(12, len(points)):]
    center = old_center.copy()
    for _ in range(6):
        selected = points[pool]
        initial = np.r_[center, np.max(np.sum((selected-center)**2, axis=1))+1e-8]
        constraint = lambda z: z[4]-np.sum((selected-z[:4])**2, axis=1)
        jacobian = lambda z: np.c_[2*(selected-z[:4]), np.ones(len(selected))]
        answer = minimize(lambda z: z[4], initial,
                          jac=lambda z: np.array([0., 0., 0., 0., 1.]),
                          constraints={'type': 'ineq', 'fun': constraint, 'jac': jacobian},
                          method='SLSQP', callback=search_budget,
                          options={'maxiter': 70, 'ftol': 1e-11})
        center = answer.x[:4]
        squared = np.sum((points-center)**2, axis=1)
        if squared.max() <= answer.x[4]+2e-8:
            break
        pool = np.unique(np.r_[pool, np.argsort(squared)[-min(8, len(points)):]])
    return center


def lloyd_minimax(points, centers, steps):
    """Alternate nearest-center assignment and cluster enclosing-ball updates."""
    best_error, best = float('inf'), centers.copy()
    for _ in range(steps):
        search_budget()
        squared = distances(points, centers)
        error = float(squared.min(axis=1).max())
        if error < best_error:
            best_error, best = error, centers.copy()
        labels = squared.argmin(axis=1)
        new = np.array([minimum_ball(points[labels == j], centers[j]) for j in range(7)])
        if np.max(np.linalg.norm(new-centers, axis=1)) < 2e-6:
            centers = new
            break
        centers = new
    return centers if distances(points, centers).min(axis=1).max() < best_error else best


def smooth_minimax(points, centers, temperature, beta, iterations):
    """Smooth min over centers and max over sample points, with exact gradient.

    g_i = -T log sum_j exp(-|x_i-c_j|^2/T),
    objective = beta^(-1) log sum_i exp(beta*g_i).
    Decreasing T and increasing beta approximate the sampled minimax problem.
    """
    def objective(vector):
        c = vector.reshape(7, 4)
        squared = distances(points, c)
        logits = -squared/temperature
        normalizer = logsumexp(logits, axis=1)
        soft_min = -temperature*normalizer
        point_weights = np.exp(beta*(soft_min-soft_min.max()))
        point_weights /= point_weights.sum()
        center_weights = np.exp(logits-normalizer[:, None])
        weights = point_weights[:, None]*center_weights
        gradient = 2*(weights.sum(axis=0)[:, None]*c-weights.T@points)
        value = soft_min.max()+logsumexp(beta*(soft_min-soft_min.max()))/beta
        return value, gradient.ravel()
    answer = minimize(objective, centers.ravel(), jac=True, method='L-BFGS-B',
                      callback=search_budget, options={'maxiter': iterations,
                      'ftol': 1e-11, 'gtol': 1e-7, 'maxls': 35})
    return answer.x.reshape(7, 4)


class FloatModel:
    """Fast numerical version of the lifted test, used ONLY to select proposals."""
    def __init__(self, exact_model):
        self.r = float(exact_model['r_upper'])
        self.p = np.array([[float(v) for v in row] for row in exact_model['contacts']])
        self.b = np.array([float(v) for v in exact_model['intercepts']])

    def boundary(self, power, seed):
        # Gaussian inverse-CDF images of scrambled Sobol points give a useful
        # deterministic-seed collection of directions on S^3. Along direction
        # u, each contact inequality gives t <= A.u+sqrt((A.u)^2+b).
        directions = ndtri(np.clip(qmc.Sobol(4, scramble=True, seed=seed)
                                   .random_base2(power), 1e-12, 1-1e-12))
        directions /= np.linalg.norm(directions, axis=1)[:, None]
        products = directions@self.p.T
        radii = np.minimum(sqrt(self.r),
                           np.min(products+np.sqrt(products*products+self.b), axis=1))
        return directions*radii[:, None]

    def evaluate(self, centers, radius_squared=RADIUS_SQUARED):
        a = np.concatenate((np.c_[-2*self.p, np.ones(3)],
                            np.c_[2*centers, -np.ones(7)],
                            np.array([[0, 0, 0, 0, 1], [0, 0, 0, 0, -1]])))
        b = np.r_[self.b, np.sum(centers*centers, axis=1)-float(radius_squared), self.r, 0.]
        matrices = a[FLOAT_BASES]
        use = np.abs(np.linalg.det(matrices)) > 1e-12
        if not use.any():
            return float('inf'), np.empty((0, 4))
        try:
            vertices = np.linalg.solve(matrices[use], b[FLOAT_BASES[use]][:, :, None])[:, :, 0]
        except np.linalg.LinAlgError:
            return float('inf'), np.empty((0, 4))
        vertices = vertices[np.isfinite(vertices).all(axis=1)]
        vertices = vertices[(vertices@a.T-b).max(axis=1) <= 2e-8]
        if not len(vertices):
            return float('inf'), np.empty((0, 4))
        norm_squared = np.sum(vertices[:, :4]**2, axis=1)
        values = norm_squared-vertices[:, 4]
        use = ((values > 0) & (vertices[:, 4] >= 0)
               & (vertices[:, 4] <= min(self.b)) & (norm_squared > 0))
        # Project troublesome vertices radially to q=|X|^2. They enrich the
        # sample but do NOT prove that an actual uncovered point exists.
        witnesses = vertices[use, :4]*np.sqrt(vertices[use, 4]/norm_squared[use])[:, None]
        return float(values.max()), witnesses


def quantize(centers):
    c = np.asarray(centers, dtype=float)
    require(c.shape == (7, 4) and np.isfinite(c).all() and np.abs(c).max() < 10,
            'Invalid numerical center proposal')
    return np.rint(c*(1 << CENTER_BITS)).astype(np.int64).tolist()


def seed_configurations(data):
    """Align the 28 exploratory seed triangles with our angle convention."""
    answer = []
    for row in data.get('search_seeds', []):
        r, tau = float(F(row['r'])), float(F(row['tau']))
        h = sqrt(r-tau)
        gamma = np.arcsin(1/(2*sqrt(tau)))
        if row['shape'] == 'regular':
            A = B = np.pi/3
        elif row['shape'] == 'two_long_angles':
            A = B = gamma
        else:
            A, B = gamma, (np.pi-gamma)/2
        target = np.array([[sqrt(tau), 0],
                           [sqrt(tau)*np.cos(2*A), sqrt(tau)*np.sin(2*A)],
                           [sqrt(tau)*np.cos(2*B), -sqrt(tau)*np.sin(2*B)]])
        old = np.array([[float(F(v)) for v in p[:2]] for p in row['contacts']])
        centers = np.array([[float(F(v)) for v in c] for c in row['centers']])
        best = None
        for perm in permutations(range(3)):
            source = old[list(perm)]
            u, _, vh = np.linalg.svd(source.T@target)
            rotation = u@vh
            error = np.linalg.norm(source@rotation-target)
            if best is None or error < best[0]:
                best = error, rotation
        centers[:, :2] = centers[:, :2]@best[1]
        answer.append((np.c_[target, np.full(3, h), np.zeros(3)], centers))
    return answer


class SearchEngine:
    def __init__(self, config, data, use_bank):
        self.config = config
        self.seeds = seed_configurations(data)
        self.cache = []
        if use_bank:
            for row in data.get('bank', []):
                self.seeds.append((np.array(row['contacts'], dtype=float),
                                   np.array(row['centers'], dtype=float)/(1 << CENTER_BITS)))
        self.seed_contacts = np.array([p for p, _ in self.seeds])


    def process(self, task):
        global SEARCH_DEADLINE
        node_id, box, inherited, attempt = task
        cfg = self.config
        model = enclose(box)
        floating = FloatModel(model)
        target = midpoint_contacts(box)
        indices = np.argsort(np.sum((self.seed_contacts-target)**2, axis=(1, 2)))[:cfg['seed_candidates']]
        nearby = [self.seeds[i] for i in indices]+self.cache
        nearest = sorted(nearby, key=lambda pair: np.sum((pair[0]-target)**2))[:cfg['seed_candidates']]
        candidates = [] if inherited is None else [np.array(inherited, dtype=float)/(1 << CENTER_BITS)]
        candidates += [c.copy() for _, c in nearest]
        best, best_centers = float('inf'), None
        witnesses = np.empty((0, 4))
        certificate = None

        def evaluate(c):
            nonlocal best, best_centers, witnesses, certificate
            try:
                integers = quantize(c)
            except (ValueError, OverflowError):
                return False
            rounded = np.array(integers, dtype=float)/(1 << CENTER_BITS)
            score, new_witnesses = floating.evaluate(rounded)
            if score < best:
                best, best_centers, witnesses = score, integers, new_witnesses
            if score < -cfg['float_accept_margin']:
                exact = exact_cover(model, integers, cfg['backend'], early_failure=True)
                if exact['ok']:
                    certificate, best_centers, best = exact, integers, score
                    return True
            return False

        for c in candidates:
            # Slight radial rescalings were useful in passing from radius
            # 1/2 to 0.49995. All are rounded and subjected to the exact test.
            for scale in (1., .9998, .9995, 1.0002):
                if evaluate(scale*c):
                    break
            if certificate is not None:
                break
        if certificate is None and float(model['uncertainty']) <= cfg['search_uncertainty']:
            SEARCH_DEADLINE = time.monotonic()+cfg['search_seconds']
            seed = (cfg['seed']+node_id*1009+attempt*1000003) % (2**32-1)
            rng = np.random.default_rng(seed)
            try:
                points = floating.boundary(cfg['sample_power'], seed)
                hold = floating.boundary(cfg['sample_power']+1, (seed+1) % (2**32-1))
                # Interior samples are helpful but not logically necessary:
                # the final exact check covers the entire body in every case.
                points = np.r_[points, np.zeros((1, 4)), .5*points[::16]]
                for trial in range(cfg['starts']):
                    search_budget()
                    if best_centers is not None:
                        c = np.array(best_centers, dtype=float)/(1 << CENTER_BITS)
                        if trial:
                            c += rng.normal(scale=.015, size=(7, 4))
                    else:
                        c = .55*points[rng.choice(len(points), 7, replace=False)]
                    for _ in range(cfg['refinements']):
                        search_budget()
                        if len(witnesses):
                            points = np.r_[points, witnesses]
                        c = lloyd_minimax(points, c, 16)
                        for temperature, beta, iterations in ((.003, 350, 65),
                                (.0006, 1600, 85), (.00015, 6500, 90)):
                            search_budget()
                            c = smooth_minimax(points, c, temperature, beta, iterations)
                        c = lloyd_minimax(points, c, 6)
                        if evaluate(c):
                            break
                        farthest = np.argsort(distances(hold, c).min(axis=1))[-64:]
                        points = np.r_[points, hold[farthest]]
                    if certificate is not None:
                        break
            except SearchTimeout:
                pass
            finally:
                SEARCH_DEADLINE = float('inf')
        if best_centers is not None:
            self.cache.append((target, np.array(best_centers, dtype=float)/(1 << CENTER_BITS)))
            self.cache = self.cache[-cfg['cache_size']:]
        radius_info = {}
        if certificate is not None:
            radius_info = polish_radius(model, floating, best_centers, certificate, cfg, node_id)
        result = {'id': node_id, **radius_info, 'centers': best_centers,
                  'covered': certificate is not None, 'certificate': certificate}
        if certificate is None:
            scores = [sum(enclose(child)['uncertainty']
                          for child in children(box, d))/2 for d in range(4)]
            minimum = min(level for _, level in box)
            allowed = [d for d in range(4) if box[d][1] <= minimum+8]
            result['split'] = min(allowed, key=lambda d: (scores[d], box[d][1], d))
        return result


def polish_radius(model, floating, integers, original_certificate, config, node_id):
    """Lower a leaf's radius WITHOUT ever weakening exact acceptance.

    Floating-point bisection is only a proposal generator. Each vertex is
    affine in s because the constraint normals are fixed and the right-hand
    sides are affine in s. Precomputing the affine vertices makes bisection
    cheap. Determinant thresholds and feasibility tolerances in this step
    can only spoil a proposal: the final exact test enumerates ALL bases.

    The original target-radius certificate is kept as a fallback. After a
    failed exact proposal we move halfway back towards that known radius,
    at most three times. The result is a certified UPPER bound, not the
    exact covering radius, nor the smallest radius this test could prove.
    """
    centers = np.array(integers, dtype=float)/(1 << CENTER_BITS)
    seed = (config['seed']+node_id*9176) % (2**32-1)
    points = floating.boundary(10, seed)
    points = np.r_[points, np.zeros((1, 4)), .5*points[::16]]
    sampled = float(distances(points, centers).min(axis=1).max())
    certified_s, proof = RADIUS_SQUARED, original_certificate
    attempts = 0
    if config['radius_polish']:
        a = np.concatenate((np.c_[-2*floating.p, np.ones(3)],
                            np.c_[2*centers, -np.ones(7)],
                            np.array([[0, 0, 0, 0, 1], [0, 0, 0, 0, -1]])))
        base = np.r_[floating.b, np.sum(centers*centers, axis=1), floating.r, 0.]
        slope = np.r_[np.zeros(3), -np.ones(7), np.zeros(2)]
        matrices = a[FLOAT_BASES]
        use = np.abs(np.linalg.det(matrices)) > 1e-12
        if use.any():
            try:
                rhs = np.stack((base[FLOAT_BASES[use]], slope[FLOAT_BASES[use]]), axis=2)
                coefficients = np.linalg.solve(matrices[use], rhs)
            except np.linalg.LinAlgError:
                coefficients = None
            if coefficients is not None and np.isfinite(coefficients).all():
                def sufficient(s):
                    vertices = coefficients[:, :, 0]+s*coefficients[:, :, 1]
                    feasible = (vertices@a.T-(base+s*slope)).max(axis=1) <= 2e-8
                    if not feasible.any():
                        return True  # heuristic only; exact empty-polytope check follows
                    vertices = vertices[feasible]
                    return float((np.sum(vertices[:, :4]**2, axis=1)-vertices[:, 4]).max()) < -1e-9

                low, high = 0., float(RADIUS_SQUARED)
                if sufficient(high):
                    for _ in range(24):
                        middle = (low+high)/2
                        if sufficient(middle):
                            high = middle
                        else:
                            low = middle
                    # Round the proposed squared radius upwards on a rational grid.
                    # The rounding here is convenient, not proof-critical: the
                    # exact checker sees this final rational number, not the float.
                    grid = 10**7
                    proposal = F(int(np.ceil(high*grid))+2, grid)
                    for _ in range(3):
                        if not 0 < proposal < RADIUS_SQUARED:
                            break
                        attempts += 1
                        candidate = exact_cover(model, integers, config['backend'],
                                                early_failure=True, radius_squared=proposal)
                        if candidate['ok']:
                            certified_s, proof = proposal, candidate
                            break
                        proposal = (proposal+RADIUS_SQUARED)/2
    return {'radius_squared': str(certified_s),
            'radius_upper_decimal': radius_decimal(certified_s),
            'sampled_radius_squared_heuristic': sampled,
            'radius_exact_proposals': attempts,
            'radius_vertex_maximum': proof.get('maximum')}


def initialize_search(config, data, use_bank):
    global SEARCH_ENGINE
    numerical_imports()
    SEARCH_ENGINE = SearchEngine(config, data, use_bank)


def search_job(task):
    return SEARCH_ENGINE.process(task)


@contextmanager
def checkpoint_lock(directory):
    """OS lock is released even after a crash; a stale filename is harmless."""
    stream = (directory/'search.lock').open('a+b')
    stream.seek(0, 2)
    if stream.tell() == 0:
        stream.write(b'0')
        stream.flush()
    stream.seek(0)
    try:
        if os.name == 'nt':
            import msvcrt
            msvcrt.locking(stream.fileno(), msvcrt.LK_NBLCK, 1)
        else:
            import fcntl
            fcntl.flock(stream.fileno(), fcntl.LOCK_EX | fcntl.LOCK_NB)
        yield
    finally:
        stream.close()


def open_search(path, config, resume, data_hash, retry_blocked=False):
    exists = path.exists()
    require(exists == resume, 'Use --resume for an existing checkpoint, or choose a new --out')
    conn = sqlite3.connect(path)
    conn.row_factory = sqlite3.Row
    conn.execute('PRAGMA journal_mode=DELETE')
    conn.execute('PRAGMA synchronous=FULL')
    code_hash = hashlib.sha256(Path(__file__).read_bytes()).hexdigest()
    if exists:
        meta = json.loads(conn.execute('SELECT value FROM meta').fetchone()[0])
        require(meta['format'] == FORMAT, 'Incompatible checkpoint')
        require(meta['code_hash'] == code_hash, 'Checkpoint was made by a different script version')
        require(meta['data_hash'] == data_hash, 'Checkpoint seed data changed')
        require(meta['config'] == config, 'Search settings changed: resume with the same settings')
        with conn:
            conn.execute("UPDATE nodes SET status='pending' WHERE status='running'")
            if retry_blocked:
                conn.execute("UPDATE nodes SET status='pending' WHERE status='blocked'")
        return conn
    conn.executescript('''
        CREATE TABLE meta(value TEXT NOT NULL);
        CREATE TABLE nodes(id INTEGER PRIMARY KEY, parent INTEGER,
          side INTEGER, box TEXT NOT NULL, depth INTEGER NOT NULL, status TEXT NOT NULL,
          split INTEGER, centers TEXT, attempts INTEGER NOT NULL DEFAULT 0,
          radius_squared TEXT, sampled_radius_squared REAL, radius_vertex_maximum TEXT,
          radius_exact_proposals INTEGER, UNIQUE(parent,side));
        CREATE INDEX pending ON nodes(status,depth,id);
        CREATE INDEX child ON nodes(parent,side);
    ''')
    with conn:
        conn.execute('INSERT INTO meta VALUES(?)', (json.dumps({'format': FORMAT,
                     'config': config, 'code_hash': code_hash, 'data_hash': data_hash}),))
        conn.execute("INSERT INTO nodes(box,depth,status) VALUES(?,0,'pending')", (json.dumps(ROOT),))
    return conn


def checkpoint_status(conn, elapsed=0., processed=0):
    counts = dict(conn.execute('SELECT status,count(*) FROM nodes GROUP BY status'))
    by_depth = list(conn.execute('SELECT depth,count(*) FROM nodes WHERE status=\'covered\' GROUP BY depth'))
    covered_volume = sum((F(count, 2**d) for d, count in by_depth), F(0))
    worst = None
    for row in conn.execute("SELECT id,box,depth,radius_squared FROM nodes WHERE status='covered'"):
        s = F(row['radius_squared'])
        if worst is None or s > worst[0]:
            worst = s, row['id'], row['box'], row['depth']
    unresolved = sum(counts.get(key, 0) for key in ('pending', 'running', 'blocked'))
    complete = unresolved == 0 and covered_volume == 1
    result = {'utc': time.strftime('%Y-%m-%dT%H:%M:%SZ', time.gmtime()),
              'elapsed_seconds_this_invocation': elapsed,
              'nodes_processed_this_invocation': processed,
              'counts': counts, 'covered_parameter_volume': str(covered_volume),
              'covered_parameter_percent': float(100*covered_volume),
              'unresolved_leaves': unresolved, 'search_complete': complete,
              'largest_covered_radius_squared': None if worst is None else str(worst[0]),
              'largest_covered_radius_upper_decimal': None if worst is None else radius_decimal(worst[0]),
              'bottleneck_node': None if worst is None else worst[1],
              'bottleneck_box': None if worst is None else json.loads(worst[2]),
              'bottleneck_depth': None if worst is None else worst[3],
              'radius_scope': 'Largest among accepted boxes only. Unresolved boxes have no radius certificate.',
              'next_step': 'Run a full independent verify after search_complete becomes true.'}
    return result


def export_checkpoint(conn, path, csv_path=None):
    """Export a full tree, including explicit nulls for every unresolved region.

    A leaf index refers to a pair (centre configuration, exact squared
    radius). Thus identical centres can have different certified radii on
    different boxes. Neither the database's status nor the stored vertex
    maximum is trusted by the independent verify command.
    """
    rows = {row['id']: dict(row) for row in conn.execute('SELECT * FROM nodes')}
    child_map = {}
    for row in rows.values():
        if row['parent'] is not None:
            child_map.setdefault(row['parent'], {})[row['side']] = row['id']
    roots = [row for row in rows.values() if row['parent'] is None]
    require(len(roots) == 1, 'Expected exactly one parameter cube')
    centers, center_indices, certificates, certificate_indices = [], {}, [], {}
    tree, stack = [], [roots[0]['id']]
    leaf_rows = []
    while stack:
        row = rows[stack.pop()]
        if row['status'] == 'split':
            tree.append(row['split'])
            require(set(child_map[row['id']]) == {0, 1}, 'Checkpoint split missing child')
            stack.extend([child_map[row['id']][1], child_map[row['id']][0]])
        elif row['status'] == 'covered':
            c = json.loads(row['centers'])
            check_centers(c)
            key = tuple(v for p in c for v in p)
            if key not in center_indices:
                center_indices[key] = len(centers)
                centers.append(c)
            s = str(F(row['radius_squared']))
            require(0 < F(s) <= RADIUS_SQUARED, 'Invalid checkpoint radius')
            pair = center_indices[key], s
            if pair not in certificate_indices:
                certificate_indices[pair] = len(certificates)
                certificates.append({'center': pair[0], 'radius_squared': s})
            position = len(tree)
            tree.append(-certificate_indices[pair]-1)
            leaf_rows.append((position, row))
        else:
            tree.append(None)
    _, statistics = tree_leaves(tree, len(certificates), allow_pending=True)
    data = {'format': FORMAT, 'parameterization': PARAMETERIZATION,
            'center_bits': CENTER_BITS, 'target_radius_squared': str(RADIUS_SQUARED),
            'r_interval': [str(R_LO), str(R_HI)],
            'centers': centers, 'certificates': certificates, 'tree': tree,
            'search_statistics': statistics}
    save_data(path, data)
    if csv_path is not None:
        csv_path = Path(csv_path)
        temporary = csv_path.with_name(csv_path.name+'.tmp')
        with gzip.open(temporary, 'wt', encoding='utf-8', newline='') as stream:
            writer = csv.writer(stream)
            writer.writerow(['preorder_position', 'node_id', 'box', 'depth',
                             'r_lower_exact', 'r_upper_exact',
                             'certified_radius_squared_exact', 'certified_radius_upper_decimal',
                             'sampled_radius_squared_heuristic', 'radius_vertex_maximum_exact'])
            for pos, row in leaf_rows:
                box = tuple(map(tuple, json.loads(row['box'])))
                r, _, _ = parameter_r_interval(box[0])
                rl, ru = r.endpoints()
                s = F(row['radius_squared'])
                writer.writerow([pos, row['id'], row['box'], row['depth'], str(rl), str(ru),
                                 str(s), radius_decimal(s), row['sampled_radius_squared'],
                                 row['radius_vertex_maximum']])
        temporary.replace(csv_path)
    return statistics


def save_data(path, data):
    path = Path(path)
    temporary = path.with_name(path.name+'.tmp')
    with temporary.open('wb') as output:
        with gzip.GzipFile(filename='', mode='wb', fileobj=output, mtime=0) as zipped:
            zipped.write(json.dumps(data, separators=(',', ':'), allow_nan=False).encode('utf-8'))
    temporary.replace(path)


def find_command(args):
    data = read_gzip_json(args.seeds)
    require(data.get('format') == SEED_FORMAT and data.get('center_bits') == CENTER_BITS,
            'Incompatible seed bank')
    numerical_imports()  # fail before checkpoint creation if dependencies are missing
    _, backend = rational_type(args.backend)
    use_bank = not args.fresh_seeds
    config = dict(SEARCH_DEFAULTS, backend=backend, search_seconds=args.search_seconds,
                  starts=args.starts, sample_power=args.sample_power, use_bank=use_bank,
                  max_depth=args.max_depth, radius_polish=not args.no_radius_polish)
    args.out.mkdir(parents=True, exist_ok=True)
    data_hash = hashlib.sha256(args.seeds.read_bytes()).hexdigest()
    with checkpoint_lock(args.out):
        conn = open_search(args.out/'search.sqlite', config, args.resume, data_hash, args.retry_blocked)
        start = time.monotonic()
        deadline = start+3600*args.hours
        completed = 0
        next_report = start
        inflight = {}
        pool = ProcessPoolExecutor(args.workers, mp_context=mp.get_context('spawn'),
                   initializer=initialize_search, initargs=(config, data, use_bank))

        def report():
            status = checkpoint_status(conn, time.monotonic()-start, completed)
            write_json(args.out/'status.json', status)
            with (args.out/'progress.jsonl').open('a', encoding='utf-8') as stream:
                stream.write(json.dumps(status)+'\n')
            print(f"Elapsed {status['elapsed_seconds_this_invocation']/3600:.3f} h; "
                  f"covered {status['covered_parameter_percent']:.6f}% of parameter cube; "
                  f"{status['counts']}; max covered rho <= "
                  f"{status['largest_covered_radius_upper_decimal']}", flush=True)
            return status

        try:
            while True:
                can_submit = time.monotonic() < deadline
                if args.max_nodes:
                    can_submit &= completed+len(inflight) < args.max_nodes
                if can_submit:
                    slots = args.workers-len(inflight)
                    if args.max_nodes:
                        slots = min(slots, args.max_nodes-completed-len(inflight))
                    rows = conn.execute("SELECT * FROM nodes WHERE status='pending' ORDER BY depth,id LIMIT ?",
                                        (slots,)).fetchall()
                    for row in rows:
                        box = tuple(map(tuple, json.loads(row['box'])))
                        inherited = json.loads(row['centers']) if row['centers'] else None
                        with conn:
                            conn.execute("UPDATE nodes SET status='running',attempts=attempts+1 WHERE id=?", (row['id'],))
                        task = row['id'], box, inherited, row['attempts']
                        inflight[pool.submit(search_job, task)] = row['id']
                if not inflight:
                    break
                done, _ = wait(inflight, timeout=1, return_when=FIRST_COMPLETED)
                for future in done:
                    node_id = inflight.pop(future)
                    result = future.result()
                    require(result['id'] == node_id, 'Unexpected worker result')
                    row = conn.execute('SELECT * FROM nodes WHERE id=?', (node_id,)).fetchone()
                    centers = json.dumps(result['centers']) if result['centers'] is not None else row['centers']
                    with conn:
                        if result['covered']:
                            require(result['certificate']['ok'], 'Missing exact acceptance')
                            s = F(result['radius_squared'])
                            require(0 < s <= RADIUS_SQUARED, 'Invalid accepted radius')
                            conn.execute("UPDATE nodes SET status='covered',centers=?,radius_squared=?,"
                                         "sampled_radius_squared=?,radius_vertex_maximum=?,radius_exact_proposals=? WHERE id=?",
                                         (centers, str(s), result['sampled_radius_squared_heuristic'],
                                          result['radius_vertex_maximum'], result['radius_exact_proposals'], node_id))
                        elif row['depth'] >= args.max_depth:
                            conn.execute("UPDATE nodes SET status='blocked',centers=? WHERE id=?", (centers, node_id))
                        else:
                            dim = result['split']
                            conn.execute("UPDATE nodes SET status='split',split=?,centers=? WHERE id=?", (dim, centers, node_id))
                            box = tuple(map(tuple, json.loads(row['box'])))
                            for side, child in enumerate(children(box, dim)):
                                conn.execute("INSERT INTO nodes(parent,side,box,depth,status,centers) VALUES(?,?,?,?,'pending',?)",
                                             (node_id, side, json.dumps(child), row['depth']+1, centers))
                    completed += 1
                if time.monotonic() >= next_report:
                    report()
                    next_report = time.monotonic()+30
        except KeyboardInterrupt:
            print('Interrupted; finishing current workers, then preserving unresolved boxes.', flush=True)
        finally:
            pool.shutdown(wait=True, cancel_futures=True)
            with conn:
                conn.execute("UPDATE nodes SET status='pending' WHERE status='running'")
            export_checkpoint(conn, args.out/'search_certificate.json.gz', args.out/'radius_by_box.csv.gz')
            status = report()
            conn.close()
        print(json.dumps(status, indent=2), flush=True)
    return 0


def export_command(args):
    with checkpoint_lock(args.out):
        require((args.out/'search.sqlite').exists(), 'Checkpoint not found')
        conn = sqlite3.connect(args.out/'search.sqlite')
        conn.row_factory = sqlite3.Row
        try:
            meta = json.loads(conn.execute('SELECT value FROM meta').fetchone()[0])
            require(meta['format'] == FORMAT, 'Incompatible checkpoint')
            require(meta['code_hash'] == hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
                    'Use the script version that created this checkpoint')
            path = args.data or args.out/'search_certificate.json.gz'
            print(json.dumps(export_checkpoint(conn, path, args.out/'radius_by_box.csv.gz'), indent=2))
        finally:
            conn.close()
    return 0


def main():
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    sub = parser.add_subparsers(dest='command', required=True)
    verify = sub.add_parser('verify', help='Reconstruct and check each box at its own exact radius')
    verify.add_argument('--data', type=Path, required=True)
    verify.add_argument('--workers', type=int, default=1)
    verify.add_argument('--backend', choices=('fraction', 'gmp', 'auto'), default='fraction')
    verify.add_argument('--sample', type=int, help='Diagnostic only; never declares a complete certificate')
    verify.add_argument('--allow-partial', action='store_true', help='Check covered leaves but allow unresolved regions; NOT a full proof')
    verify.add_argument('--report', type=Path)
    verify.set_defaults(function=verify_command)
    find = sub.add_parser('find', help='One-cube adaptive search with exact radius certificates')
    find.add_argument('--seeds', type=Path, default=DEFAULT_SEEDS)
    find.add_argument('--workers', type=int, default=min(16, max(1, (os.cpu_count() or 2)-2)))
    find.add_argument('--backend', choices=('fraction', 'gmp', 'auto'), default='auto')
    find.add_argument('--out', type=Path, default=Path('single_run'))
    find.add_argument('--hours', type=float, default=24)
    find.add_argument('--resume', action='store_true')
    find.add_argument('--fresh-seeds', action='store_true', help='Use only the 28 exploratory seeds, omitting the old centre bank')
    find.add_argument('--no-radius-polish', action='store_true', help='Certify every accepted leaf at the fixed target .49995')
    find.add_argument('--retry-blocked', action='store_true', help='On resume, retry unresolved leaves at the depth limit')
    find.add_argument('--search-seconds', type=float, default=8)
    find.add_argument('--starts', type=int, default=2)
    find.add_argument('--sample-power', type=int, default=11)
    find.add_argument('--max-depth', type=int, default=48)
    find.add_argument('--max-nodes', type=int, default=0, help='Processed-node limit for this invocation; 0 means unlimited')
    find.set_defaults(function=find_command)
    export = sub.add_parser('export', help='Export a stopped checkpoint, including radius_by_box.csv.gz')
    export.add_argument('--out', type=Path, required=True)
    export.add_argument('--data', type=Path)
    export.set_defaults(function=export_command)
    args = parser.parse_args()
    if args.command in ('verify', 'find'):
        require(args.workers >= 1, 'Worker count must be positive')
    if args.command == 'verify':
        require(args.sample is None or args.sample > 0, 'Sample count must be positive')
    if args.command == 'find':
        require(args.hours > 0 and args.search_seconds > 0 and args.starts > 0
                and 8 <= args.sample_power <= 15 and 1 <= args.max_depth <= 128
                and args.max_nodes >= 0, 'Invalid search settings')
        require(not args.retry_blocked or args.resume, '--retry-blocked requires --resume')
    return args.function(args)


if __name__ == '__main__':
    mp.freeze_support()
    raise SystemExit(main())
