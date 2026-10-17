# Test cases

This directory contains the four test cases whose quantile sets are
small enough to require the SMC machinery (`QSI_EEM_SMC`). The three
moderately difficult problems of the paper live in the
[`contrib-qsi`](https://github.com/stk-kriging/contrib-qsi)
dependency (`lib/algorithms/qsi/test_functions` after `startup`) and
are catalogued here as well, so that this page lists every problem used
in the article. See the [root README](../README.md#describing-a-test-case)
for how to build a new case, and `paper_get_test_case.m` for the
registered names.

In every case, the problem is to estimate the set of designs x such
that P(f(x, S) in C) <= alpha.

## Overview

| Case | dim x + s | Critical region C | alpha | S distribution | vol(Γ)/vol(X) | Used in |
|---|---|---|---|---|---|---|
| `branin_mod` | 1 + 1 | z < 7.5 | 0.05 | 15 · Beta(7.5, 1.9) | — | "Results on moderately difficult problems", demo |
| `double_camel` | 2 + 2 | z < 1.2 | 0.15 | uniform | — | "Results on moderately difficult problems" |
| `hart4` | 2 + 2 | z > −1.1 | 0.6 | uniform | — | "Results on moderately difficult problems" |
| `piston_smc` | 4 + 2 | z < 1.148 | 0.05 | Beta(1/2, 1/2) margins | 4.7e-6 | "Synthetic functions" |
| `otl_smc` | 5 + 1 | z < 97.35 | 0.05 | N(175, 50) truncated to [50, 300] | 8.2e-7 | "Synthetic functions" |
| `trid_smc` | 4 + 3 | z < 5.3 | 0.1 | uniform | 1.9e-6 | "Synthetic functions" |
| `rotor37` | 13 + 5 | \|z1\| > 0.18 or \|z2\| > 0.18 | 0.05 | uniform | 1.7e-8 | "Application to ROTOR37" |

The volumes of the small quantile sets are measured by subset
simulation (`deterministic_subset('<name>')`, full configuration, one
seed — five for ROTOR37). The
sets of the first three cases are not small: their estimation error is
measured directly on a product grid (`extract_deviation`), and no
reference volume is needed.

## The moderately difficult cases (contrib-qsi)

**`branin_mod`** — a modified Branin-Hoo function on [0, 10] x [0, 15],
one design and one uncertain variable, introduced with the QSI-SUR
criterion (Ait Abdelmalek-Lomenech et al., Technometrics, 2025).
Illustrated in Figure 1 of the paper and used by the demo.

**`double_camel`** — a six-hump-camel-type function with two design and
two uncertain variables.

**`hart4`** — the four-dimensional Hartmann function, two design and
two uncertain variables, with a critical region of the form
{z > threshold} and a large alpha.

## The small-quantile-set cases

**`piston_smc`** — the piston simulation function (Kenett & Zacks,
1998): cycle time of a piston as a function of four design variables
(M, S, V0, T0) and two uncertain ones (P0, Ta), with arcsine
(Beta(1/2, 1/2)) margins on the uncertain variables.

**`otl_smc`** — the output transformerless (OTL) circuit function
(Ben-Ari & Steinberg, 2007): mid-point voltage as a function of five
design resistances and the uncertain current gain beta. The output is
returned reflected (100 − Vm) so that the constraint reads as an upper
bound.

**`trid_smc`** — the Trid function on [−49, 49]^7, four design and
three uncertain variables, rescaled as (10000 − Trid(u)) / 1000.

**`rotor37`** — the industrial case: Gaussian process metamodels of the
ROTOR37 transonic compressor (Safran Tech), 13 design and 5 uncertain
variables, two constrained outputs (mass flow and pressure ratio,
expressed as relative deviations from their baselines, tolerance 0.18).
The metamodels are proprietary and not publicly available at the
moment: only the definition of the QSI problem is provided. If it is
available to you, the proprietary code must be placed in
[`rotor37/code`](./rotor37/code).
