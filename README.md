# Borsuk number in four dimensions is at most seven

Computational supplement to the article by A. Arman, A. Bondarenko,
A. Prymak and D. Radchenko.

**Repository URL:**
<https://github.com/andriyprm/4d7balls>

The article proves that every subset of four-dimensional Euclidean space
of diameter one is covered by seven balls of a common diameter strictly
smaller than $0.9999$. Thus the Borsuk number satisfies
$b_4\le7$, and the equal-diameter ball-covering number satisfies $g_4\le7$.

The circumball of a bounded set is its smallest closed enclosing ball.
Write $r$ for its **squared** radius. The target covering radius and its square are

$$
\rho=\frac{9999}{20000}=0.49995,
\qquad \sigma=\rho^2=\frac{99980001}{400000000},
\qquad r_0=\frac{38401}{100000}.
$$

One four-parameter cube covers the computational range
$299/1000\le r\le r_0$. The small and large ranges are handled by the
analytic arguments in the article. No subdivision from a previous search
is needed to verify this cube.

## Files and their roles

| File | Role |
| --- | --- |
| `seven_balls.tex`, `seven_balls.pdf` | Article: geometric reductions, analytic covering arguments, and finite covering condition. |
| `seven_balls_single.py` | Exact verifier and optional numerical finder. Parts I–IV contain the verification routines; Part V contains the numerical search. |
| `search_certificate.json.gz` | Complete proof data: the subdivision tree, integer centres, and exact rational squared radii. |
| `VERIFICATION.md` | Parameter formulas, exact interval rules, finite test, data format, analytic constant checks, and scope of the computation. |
| `SEARCH.md` | Numerical sampling, smooth objectives, optimization schedules, seeds, subdivision, and radius polishing. These are discovery details, not proof inputs. |
| `verify_report.json` | Result of a full exact replay of the supplied code and data. This records evidence of a completed run; it is not an input trusted by the verifier. |
| `seven_balls_seeds.json.gz` | Optional starting configurations for the finder. Not used by the verifier. |
| `SHA256SUMS` | Digests of the computational files and documentation. |

Only the Python program and `search_certificate.json.gz` are needed to run
the exact replay. The geometric reductions in the article and the interval
arguments in `VERIFICATION.md` are also part of the mathematical proof.

The LaTeX source contains its bibliography and both TikZ figures. With a
standard TeX Live or MiKTeX installation, compile it twice:

```sh
pdflatex seven_balls.tex
pdflatex seven_balls.tex
```

## Verify the complete certificate

Use Python 3.10 or later. From this directory, the standard-library route is:

```sh
python seven_balls_single.py verify --data search_certificate.json.gz --workers 12 --backend fraction --report verify_report.json
```

For faster exact rational arithmetic, install the optional GMP interface:

```sh
python -m pip install gmpy2
python seven_balls_single.py verify --data search_certificate.json.gz --workers 12 --backend gmp --report verify_report.json
```

On Windows, `py` can replace `python`. Each command above is a single line.
Reduce `--workers` if fewer processor cores are available. Both backends
perform exact rational arithmetic; neither NumPy nor SciPy is imported by
the verifier. The seed file is not read during verification.

A successful complete replay has all the following fields:

```json
{
  "complete_computational_certificate": true,
  "checked": 23494,
  "covered_volume": "1",
  "unresolved_volume": "0",
  "failures": [],
  "uniform_all_r_radius_upper_decimal": "0.4999500000"
}
```

The tree has 23,494 leaves, 23,493 splits, and maximum depth 29. The data
contain 6,913 distinct centre configurations and 19,656 centre/radius
pairs. Of the leaves, 23,466 have a certified squared radius below
$\sigma$; 28 retain $\sigma$. These are sufficient radii, not claims of
optimality.

The decimal `0.4999500000` is an upward-rounded upper bound. The proof
establishes coverage by **open** balls of the stored radii. Compactness
then permits a strictly smaller common radius, giving the strict diameter
bound $2\rho=0.9999$. A decimal value equal to the target does not mean
that strictness has been lost.

The verifier reconstructs every parameter box, interval enclosure, and
linear system. It does not trust previously computed interval bounds,
sampled radii, finder logs, or the supplied replay report. A sampled run
(`--sample`) or a tree with unresolved boxes does not prove the theorem.
Check the `complete_computational_certificate` field, not just the process
exit code: successful diagnostic runs can also exit with code zero.

## Find a new certificate (optional)

The numerical method and its default schedules are documented in
`SEARCH.md`; the commands below start or resume an experiment.

Finding centres requires NumPy and SciPy; GMP remains optional but is
recommended for the exact checks made during the search:

```sh
python -m pip install numpy scipy gmpy2
python seven_balls_single.py find --out single_run --hours 24 --workers 12 --backend gmp
```

The default seed file contains 28 exploratory configurations and 2,897
configurations from earlier searches. These are numerical proposals only.
The finder starts from the entire parameter cube and imports no earlier
subdivision. Add `--fresh-seeds` to use only the exploratory configurations.
Add `--no-radius-polish` to retain the fixed target radius instead of trying
to improve each accepted leaf's radius.

To resume, keep the same settings and append `--resume`:

```sh
python seven_balls_single.py find --out single_run --hours 24 --workers 12 --backend gmp --resume
```

The 24-hour limit is a search budget, not a completion guarantee. Current
worker jobs may finish after the budget expires. The SQLite checkpoint
preserves unresolved boxes. A timeout or unsuccessful numerical
optimization never makes a pending box disappear.

Verify a newly exported certificate with a separate full replay:

```sh
python seven_balls_single.py verify --data single_run/search_certificate.json.gz --workers 12 --backend gmp --report single_run/verify_report.json
```

The `export` command can export a checkpoint without continuing the search;
run `python seven_balls_single.py export --help` for its options. The export
also includes `radius_by_box.csv.gz` for analysis. This table is not used
to accept the proof.

**A reported `bottleneck_box` is already covered.** It identifies a leaf
with the largest retained radius, and remains in the progress output after
completion. Search completion requires `search_complete: true` and
`unresolved_leaves: 0`; the proof check is the subsequent successful full
verifier run.

## Reproducibility and scope

The supplied certificate gives a deterministic exact replay. Finding a
new certificate is a different task: processor scheduling, random seeds,
optimization libraries, and time limits can change the discovered tree.
No claim is made that a fresh search will reproduce the same centres.

The verifier checks a finite mathematical certificate. It is not a
formalization in a proof assistant, and its two rational backends are not
independent implementations of the geometric method. The article's
geometric reductions and the enclosure rules still require mathematical
review. See `VERIFICATION.md` for those rules and for the exact point at
which each numerical quantity becomes a certified one.
