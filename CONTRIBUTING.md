# Contributing / reviewing

Thanks for taking a look. This package accompanies the method paper
(doi:[10.5281/zenodo.20033072](https://doi.org/10.5281/zenodo.20033072)) and is
still v0.1 — review of both the statistics and the code is genuinely useful at
this stage, not a formality.

## Where to start

- **[METHOD.md](METHOD.md)** is the substantive document: the estimand, the
  four validity conditions, known limitations ranked by severity, open
  questions, and (§9) a full treatment of the anchored difference-in-differences
  fallback including a novelty assessment against the DTA meta-analysis
  literature. §9's novelty claims are an explicit request for scrutiny — if
  you know of prior art that anticipates the raw-scale, Bucher-framed
  sensitivity DiD described there, that is exactly the kind of finding worth
  opening an issue about.
- **[README.md](README.md)** covers usage and architecture.
- `tests/run_tests.R` (61 checks) is the correctness contract — anything you
  change should keep it green (`PKG_ROOT="$(pwd)" Rscript tests/run_tests.R`).

## How to review

- **Statistical review** (no R needed): read METHOD.md, comment via a GitHub
  Issue. Most useful on §4 (validity conditions), §5 (limitations), and §9
  (the DiD estimator and its literature position).
- **Code review**: open a PR. `R CMD check` runs automatically via GitHub
  Actions on every PR (the `stan`/`stan_joint` engines are skipped in CI since
  they need a C++ toolchain — run those locally if you touch `R/meta.R`'s Stan
  paths).
- **Reproducing a result**: `inst/examples/` has three scripts —
  standardization, the S_k diagnostics (forest plots, prior-sensitivity sweep),
  and the DiD fallback. Each is runnable standalone once the package loads.

## Local setup

```r
# from the package root
for (f in list.files("R", pattern = "\\.R$", full.names = TRUE)) source(f)
for (d in list.files("data", full.names = TRUE)) load(d)
```

or `R CMD INSTALL .` after cloning. `Suggests: brms, metafor, lme4, dplyr,
tidyr, patchwork, knitr` — only `ggplot2` and base `stats` are hard
dependencies.

## Scope note

v0.1 implements spectrum-standardized **sensitivity** only. Specificity,
predictive values, and ROC-based measures are explicitly out of scope for this
release (METHOD.md §3) — proposals for those are welcome as issues, but are
larger asks than a PR.
