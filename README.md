# qsi-eem-paper

[![license](https://img.shields.io/github/license/stk-kriging/qsi-eem-paper)](LICENSE.md)
[![last commit](https://img.shields.io/github/last-commit/stk-kriging/qsi-eem-paper/main)](https://github.com/stk-kriging/qsi-eem-paper/commits/main)

This repository contains the Matlab implementation and the numerical
experiments of the article:

Romain Ait Abdelmalek-Lomenech, Julien Bect and Emmanuel Vazquez (†)  
__Bayesian Active Learning of (small) Quantile Sets through Expected
Estimator Modification__
([hal-05111516](https://hal.science/hal-05111516),
[arXiv:2506.13211](https://arxiv.org/abs/2506.13211))

(†) Université Paris-Saclay, CNRS, CentraleSupélec,
[Laboratoire des signaux et systèmes](https://l2s.centralesupelec.fr/),
Gif-sur-Yvette, France.

The article deals with quantile set inversion (QSI): estimating the set
of designs x for which some outputs of an expensive simulator f(x, S),
subject to uncertain inputs S, stay within prescribed bounds with high
probability. Evaluation points are selected sequentially by the
Expected Estimator Modification (EEM) criterion — the expected change
of the estimated set brought by a candidate batch, under Gaussian
process models of the outputs — and, when the target set is too small
to be found directly, the algorithm tracks it through a sequence of
progressively tightened constraint sets, sampled by sequential Monte
Carlo (`QSI_EEM_SMC`).

The QSI-SUR criterion and the other Bayesian set inversion strategies
used for comparison come from the
[`contrib-qsi`](https://github.com/stk-kriging/contrib-qsi) repository,
as do the moderately difficult test problems (modified Branin-Hoo,
double camel, Hartmann-4).


## Setup

### Matlab/Octave

The code in this repository should work with any reasonably recent
release of Matlab (e.g., R2024a).

Reproducing the paper campaigns requires the Parallel Computing
Toolbox: the `run_*` and `extract_*` scripts of `experiments/` drive
their runs through 50-worker pools. The algorithms themselves
(`QSI_EEM`, `QSI_EEM_SMC`, the subset-simulation reference path) and
the demo are serial — a single run needs no toolbox beyond those used
by STK.

The core of the code (everything except the campaign scripts) should
run under Octave (with the statistics package): `startup`, the
algorithms (`QSI_EEM`, `QSI_EEM_SMC`), the test cases and the
subset-simulation reference path. Still Matlab-only: the campaign
drivers (Parallel Computing Toolbox), the figure and table
scripts (`make_*`, which rely on script-local functions), and the
comparison algorithms of `contrib-qsi` (plain `readmatrix` calls).

A C compiler is required to compile some MEX-files located in the
dependencies (STK, BSS).

### Dependencies

This repository depends on a few other mandatory Matlab/Octave
packages/repositories: STK, BSS (contrib-bss), QSI (contrib-qsi).

Under Octave, the
[statistics package](https://gnu-octave.github.io/packages/statistics/)
is also required (install it with `pkg install -forge statistics`);
`startup.m` loads it automatically.

The script `startup.m`, located at the root of this repository,
fetches them at pinned commits and sets the path (assuming git is
installed).

NB: `startup.m` is run automatically if you start Matlab/Octave from
the root of this repository. Otherwise, you have to run it manually.
All the other scripts assume that `startup.m` has been run.

NB: In Matlab, STK uses MathWorks' Optimization and Statistics &
Machine Learning toolboxes, if they are available. The results
presented in the paper have been obtained using these toolboxes.


## Repository layout

| Directory | Content |
|---|---|
| [`lib/algorithms/qsi_eem`](./lib/algorithms/qsi_eem) | the QSI-EEM criterion, its SMC variant (`QSI_EEM_SMC.m`) and the subset-simulation reference path |
| [`lib/utils`](./lib/utils) | campaign utilities (initial designs, extraction, figures) |
| [`test_cases`](./test_cases) | the four small-quantile-set cases; see [its README](./test_cases/README.md) for a catalogue of every problem of the paper |
| [`experiments`](./experiments) | the scripts reproducing the numerical experiments |
| [`demo`](./demo) | demonstration script |

The directories `lib/algorithms/qsi` and `lib/algorithms/bss` are the
dependency checkouts created by `startup.m`; they are not part of this
repository.


## Quick start

[`demo/qsieempaper_demo.m`](./demo/qsieempaper_demo.m) runs QSI-EEM on
the modified Branin-Hoo example and plots the selected evaluation points;
optionally, it also runs the comparison algorithms and compares the
designs they produce. It is the natural entry point for a first contact
with the method.

A minimal QSI-EEM-SMC run on an existing test case:

```matlab
[problem, config] = paper_get_test_case ('piston_smc', true);  % true: debug configuration

data_dir = fullfile (tempdir (), 'qsi_eem_quickstart');
paper_prepare_data_dir (data_dir);

paper_generate_initial_designs (problem, config, 1, data_dir);  % run id 1
QSI_EEM_SMC (problem, config, 1, 1, data_dir);                  % batches of size 1
```

With the debug configuration the run takes a few seconds; drop the
`true` for the configuration used in the paper.

Every algorithm reads and writes files named after the case, the run id
and the batch size:

| File | Content |
|---|---|
| `doe_init/doe_init_<name>_<it>_init.csv` | initial design (from `paper_generate_initial_designs`) |
| `results/design/doe_QSI_EEM_smc_<name>_<it>_<batch>.csv` | current design, updated after each step |
| `results/param/param_..._<m>_<name>_<it>_<batch>.csv` | GP parameters of output m, one row per step |
| `results/param/cov_..._<m>_<name>_<it>_<batch>.csv` | selected covariance family |
| `results/time/time_..._<name>_<it>_<batch>.csv` | criterion wall times |
| `results/smc_logs/log_..._<name>_<it>_<batch>.csv` | SMC level history |

The extraction scripts locate the runs through the same convention.


## Describing a test case

A test case defines a quantile set inversion problem: estimate the set
of designs x such that P(f(x, S) in C) <= alpha, where S is the vector
of uncertain variables and C the critical region. Each case consists of
three or four files in a subdirectory of [`test_cases`](./test_cases)
(see `piston_smc` for a minimal example, `rotor37` for a full-featured
one):

- `<name>.m` — the function itself: takes an n x (dim_x + dim_s) matrix,
  design variables first, and returns an n x M matrix of outputs;
- `<name>_struct.m` — the problem description (see below);
- `<name>_config.m` — the algorithm configuration, built on the shared
  defaults of `eem_smc_config_base` with case-specific overrides;
- `<name>_s_trnsf.m` — maps a uniform sample on `BOXs` to the actual
  distribution of S (identity if S is uniform; see `otl_smc_s_trnsf`
  for a truncated normal).

### The problem structure

`<name>_struct.m` returns `[prm, f, s_trnsf]`, with `f = @<name>`. Its
skeleton, on the piston example:

```matlab
function [prm, f, s_trnsf] = piston_smc_struct ()

f = @piston_smc;
s_trnsf = @piston_smc_s_trnsf;

prm = struct ();
prm.M = 1;                            % one output
prm.dim_x = 4;  prm.dim_s = 2;        % 4 design + 2 uncertain variables
prm.BOX  = stk_hrect (...);           % input domain, BOXx and BOXs its factors
prm.const = [-inf; 1.148];            % critical region (interval form)
prm.alpha = 0.05;
prm.name = "piston_smc";
prm.list_cov = {'stk_expcov_aniso'; ...};

end
```

The fields of `prm`:

| Field | Meaning |
|---|---|
| `M` | number of outputs |
| `dim_x`, `dim_s` | numbers of design and uncertain variables |
| `BOX`, `BOXx`, `BOXs` | input domain (`stk_hrect`) and its two factors |
| `const` | critical region, one column per output (see below) |
| `alpha` | probability threshold of the quantile set |
| `tol` | half-width of the band, for band-type regions only |
| `name` | case name, used in every result file name |
| `list_cov` | candidate covariance families for the GP models |

Two forms of critical region are supported (see `threshold_to_const`):

- **interval**: `C = {z : const(1,m) < z_m < const(2,m) for some m}`,
  with `const(1,m) < const(2,m)`. The SMC levels move the upper bound
  of the first output. Example: `prm.const = [-inf; 1.148]` (piston).
- **band complement**: `C = {z : |z_m| > tol for some m}`, declared by
  setting `prm.tol`; then `prm.const = prm.tol * [1; -1]` per output,
  and the SMC levels tighten the tolerance itself down to `prm.tol`.
  Example: ROTOR37, with `tol = 0.18` on two outputs expressed as
  relative deviations from their baselines.


### Registration and checks

Add the case to the switch of
[`test_cases/paper_get_test_case.m`](./test_cases/paper_get_test_case.m),
then:

```matlab
% quick end-to-end run (debug configuration)
[problem, config] = paper_get_test_case ('<name>', true);

% reference volume of the quantile set, by subset simulation
vol = deterministic_subset ('<name>');
```

For a very small quantile set, expect `deterministic_subset` to reach
its final level in a number of stages roughly proportional to
-log(vol). The deviation of an estimated set with respect to the true
one is computed by `extract_precision_small` from the design files of a
completed run.


## Paper experiments

Each campaign lives in its own subdirectory of
[`experiments`](./experiments), with `run_*` (runs the algorithms),
`extract_*` (post-processing) and `make_*` (figures/tables) scripts.
The header of each script states what it produces; the table below maps
the campaigns to the paper.

| Directory | Paper part |
|---|---|
| `small_synthetic` | scatter figures and summary table of the "Synthetic functions" section |
| `small_rotor37` | "Application to ROTOR37" section (requires the proprietary ROTOR37 code) |
| `qsi_comparison` | metric figures of the "Results on moderately difficult problems" section |
| `eem_sur_doe_crit_comparison` | criterion maps and design comparison of the "Approximation of the QSI-EEM criterion" section |
| `time_eem_sur` | criterion cost figure and runtime tables |
| `quadrature_size_comparison` | quadrature-size study of the supplement |
| `eem_smc_parameters_comparison` | SMC-parameters sensitivity of the supplement (incl. Table SM1) |

The campaigns are sized for a machine with many cores: the scripts
start a pool of `POOL` workers (all physical cores but one by default,
a single worker in debug mode), and dispatch the parfor iterations one
at a time (`SubrangeSize`, 1) so that a long run never blocks a queued
one behind it.


## The ROTOR37 test case

The ROTOR37 test case used in the "Application to ROTOR37" section of
the article relies on metamodels provided by Safran Tech, which are
proprietary and not publicly available. Only the
definition of the QSI problem is provided, in
[`test_cases/rotor37`](./test_cases/rotor37); it can be used together
with the proprietary code if it is available to you.


## Acknowledgements

This work has been funded by the French National Research Agency
(ANR), in the context of the project SAMOURAI (ANR-20-CE46-0013).

The authors are grateful to Sébastien Da Veiga (ENSAI, formerly Safran
Tech) and Safran Tech for providing the original metamodel of the
ROTOR37 model, which was adapted for this work.


## Copyright & license

Copyright 2025, 2026 CentraleSupélec

These computer programs are free software: you can redistribute them
and/or modify them under the terms of the GNU General Public License
as published by the Free Software Foundation, either version 3 of the
license, or (at your option) any later version.

They are distributed in the hope that they will be useful, but WITHOUT
ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
FITNESS FOR A PARTICULAR PURPOSE.  See the license for more details.

You should have received a copy of the license along with the software
(see [LICENSE.md](./LICENSE.md)).  If not, see <http://www.gnu.org/licenses/>.
