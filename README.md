# SSDTA — Spectrum-Standardized Diagnostic Test Accuracy

Diagnostic sensitivity is not a property of a test alone — it depends on the
bacillary-burden spectrum of the population. **SSDTA** estimates stratum-specific
sensitivity across studies and reweights it to a chosen target population's
burden distribution, giving a standardized, context-specific sensitivity. It
propagates posterior uncertainty in stratum-specific sensitivity and sampling
uncertainty in observed case-mix counts, conditional on the chosen `pi` scenario
and the independence assumptions described below:

```
S_std = Σ_k  w_k · S_k
```

where `S_k` is the sensitivity in burden stratum `k` (e.g. Xpert Ultra
semi-quant grade) and `w_k` is that stratum's share in the target population.

📄 **Method paper:** Schumacher SG (2026), *Spectrum-Standardised Diagnostic
Accuracy*, Zenodo preprint,
doi:[10.5281/zenodo.20033072](https://doi.org/10.5281/zenodo.20033072)
(CC BY-SA 4.0) · `citation("SSDTA")`
📘 **[METHOD.md](METHOD.md)** — validity conditions, limitations, open questions,
planned extensions, and the data release gate. Read it before reporting an
estimate.

Two things it does:

1. **Standardize one test** — report its expected sensitivity in a research,
   routine, or screening population (useful in its own right).
2. **Compare two tests without head-to-head data** — standardize each to a
   common target and difference. When one test lacks a stratum breakdown, an
   **anchored difference-in-differences** fallback partially deconfounds instead.

The stratifier is generic (Xpert Ultra semi-quant for TB; CD4 bands, Ct values,
etc. for other settings). Despite the package name, the current implemented
estimand is **sensitivity**; specificity, predictive values, and ROC curves need
different target distributions and are not yet implemented.

## Install / load

```r
for (f in list.files("SSDTA/R", pattern="\\.R$", full.names=TRUE)) source(f)
library(ggplot2); for (d in list.files("SSDTA/data", full.names=TRUE)) load(d)
```
or `install.packages("SSDTA", repos = NULL, type = "source")` then `library(SSDTA)`.

## Quick start — standardization (reproduces the paper)

```r
# S_k: shipped posterior draws from the paper's Bayesian meta-analysis
ssdta_standardize(npoc_posterior_draws, burden_profiles$research,          pi = "observed")
#   Research  80.9%  (68.1, 86.1)
ssdta_standardize(npoc_posterior_draws, burden_profiles$routine_pakistan,  pi = 0.20)
#   Routine   66.8%  (57.5, 73.7)
ssdta_standardize(npoc_posterior_draws, burden_profiles$screening_kendall, pi = 0.30)
#   Screening 53.3%  (43.6, 63.7)

# new Pakistan screening data + crude pool with Kendall
scr <- pool_profiles(burden_profiles$screening_kendall, burden_profiles$screening_pakistan)
ssdta_standardize(npoc_posterior_draws, scr, pi = 0.30)                     # 51.8%

# pi sensitivity grid + figures
grid <- ssdta_standardize(npoc_posterior_draws, burden_profiles$routine_pakistan, pi = seq(0.05, 0.4, 0.05))
plot_casemix(burden_profiles)
plot_pi_sensitivity(list(routine = grid))

# mechanism figure: w_k -> w_k * S_k -> S_std
plot_spectrum_mechanism(
  npoc_posterior_draws,
  burden_profiles[c("research", "routine_pakistan", "screening_kendall")],
  pi = list("observed", 0.20, 0.30)
)
```

Full scripts: `inst/examples/example_standardization.R`,
`inst/examples/example_meta_diagnostics.R`, and
`inst/examples/example_did_fallback.R`.

## Architecture (five layers)

| Layer | Function(s) | What it does |
|---|---|---|
| Estimation | `ssdta_meta()` | per-stratum sensitivity meta-analysis; engines `stan` (brms, the paper's method — default), `stan_joint` (shared tau), `glmm` (metafor RE), `conjugate` (fast, no deps) |
| Diagnostics | `plot_forest_stratum()`, `ssdta_meta_sensitivity()` | study-level forest within each stratum; prior / model / leave-one-out sweep |
| Weights | `burden_profile()`, `pool_profiles()` | target-population burden distributions; crude count-pooling; Dirichlet case-mix draws |
| Standardization | `ssdta_standardize()` | reweight S_k draws by w_k with posterior × Dirichlet uncertainty; `pi` scalar, grid, or `"observed"` |
| Comparison | `ssdta_compare()`, `anchored_did()` + `did_*` | difference of two standardized sensitivities; anchored-DiD fallback |

## The estimation engines

`ssdta_meta(data, engine = ...)` produces posterior sensitivity draws per
stratum (drop-in for `ssdta_standardize`). Default `"stan"` reproduces the
paper (Bayesian binomial-logit RE, Normal(0,1.5) / Exponential(1) priors; needs
`brms`). `"glmm"` (metafor logit-normal RE) and `"conjugate"` (complete-pooling
Beta-Binomial) are dependency-light alternatives that track the Stan fit closely
(e.g. glmm gives High 94% / Low 88% / Trace 38% / Not-detected 20% vs the paper's
94 / 89 / 38 / 18). Strata are fitted separately (no cross-stratum partial
pooling by default); `"stan_joint"` is the alternative — one model,
`~ 0 + stratum + (1 | study)`, sharing tau across strata, which matters because
k = 3–5 studies leaves tau weakly identified within any single stratum. Its
draws share one posterior index across strata, so the cross-stratum correlation
it induces survives into the weighted sum.

## Diagnostics: is S_k trustworthy before you reweight it?

Standardization inherits every weakness of the S_k it reweights, so two checks
belong with any reported S_std.

**Do studies agree within a stratum?** That is the conditional-transportability
condition stated as a picture. `plot_forest_stratum()` draws each study's
`tp/n` with a Wilson interval, faceted by stratum, with the pooled estimate as
a diamond — the sparse strata (Very low, Trace, Not-detected) are where one
study can carry the pooled value.

```r
plot_forest_stratum(npoc_by_ultra_stratum, pooled = npoc_posterior_draws)
```

**Does the answer survive the prior and the study set?** `ssdta_meta_sensitivity()`
refits across specifications and tabulates them side by side. On the `stan`
engine the defaults are the paper's analyses — N(0,5) intercept, joint model
with shared tau, Logistic(0,1) intercept (= Uniform(0,1) on the probability
scale), Half-t(3,0,2.5) tau — and `loo = TRUE` adds one leave-one-study-out arm
per study.

```r
sweep <- ssdta_meta_sensitivity(npoc_by_ultra_stratum, loo = TRUE)  # engine = "stan"
sweep                                   # wide "median (lo–hi)" comparison table
plot_meta_sensitivity(sweep)
as.data.frame(sweep)                    # the same table, for write.csv()
```

Build arbitrary arms with `meta_spec(label, engine=, prior_intercept=, prior_sd=,
exclude=)`. Priors only bite on the stan engines; setting them under `glmm` /
`conjugate` warns rather than silently producing a duplicate arm.

The question is never "did S_k move" but "did S_std move" — push each arm's
draws back through `ssdta_standardize()` with the same target profile before
concluding anything. In the shipped NPOC data the sparse strata swing widely
under leave-one-out (Very low: 44–72%) while the routine standardized estimate
barely moves (66.7–69.4%), because those strata carry little weight in that
target distribution.

## Two ways to compare two tests — and when each applies

- **Standardized difference** (`ssdta_compare`) — both tests have a stratum
  breakdown. Standardize each to a common `w*`, difference. This removes
  differences in the **measured** case mix at the chosen stratum resolution,
  conditional on transportability of each test's stratum-specific sensitivity.
  The function enforces identical profiles/`pi` and applies one shared
  realization of the target weights to both tests. Use
  `draw_pairing = "paired"` for aligned joint posterior draws and
  `"independent"` for separately fitted models.
- **Anchored DiD** (`anchored_did`, the folded-in module) — one test lacks a
  breakdown (e.g. **pooled Ultra** has no per-specimen semi-quant grade). Anchor
  on a common test and take a difference-in-differences: the diagnostic-accuracy
  analogue of the Bucher method for indirect treatment comparisons. This
  **partially** deconfounds — it removes the anchor level-shift but leaves the
  burden-mix × gap-slope interaction, and its reference population is implicit.
  See **[METHOD.md §9](METHOD.md#9-the-anchored-difference-in-differences-fallback--full-method)**
  for the estimator, its position in the literature, and open questions —
  including a novelty assessment we'd specifically like reviewers to
  stress-test.

The two worked examples are exactly this pairing: `example_standardization.R`
(the principled method) and `example_did_fallback.R` (NPOC vs pooled, where no
semi-quant breakdown exists).

> Note: `ssdta_compare` and the stratified/standardized comparison are the same
> estimator — "difference of standardized sensitivities" equals "standardized
> difference of within-stratum sensitivities" by linearity, provided the **same**
> weights `w*` are used for both tests.

All scalar standardization and comparison objects support `print()`, `summary()`,
`as.data.frame()`, and `plot()`.

## Data

| Object | What |
|---|---|
| `npoc_posterior_draws` | NPOC-NAAT S_k posterior draws (8000 × 6 strata) from the paper's Stan fit — authoritative S_k input |
| `npoc_stratum_summary` | per-stratum posterior medians + CrIs |
| `npoc_by_ultra_stratum` | real per-study TP/FN by Ultra stratum (GDG PICO1 Table 3) — for `ssdta_meta` |
| `npoc_studies` | characteristics of those 5 studies (country, case-finding, specimen, reference standard, n, publication status, DOI) — joins to `npoc_by_ultra_stratum` on `study` |
| `burden_profiles` | `research`, `routine_pakistan` (GxAlert 2019+2021, n=88,560), `screening_kendall` (Uganda, n=89), `screening_pakistan` (new, n=118) |

**pi (the Not-detected weight)** — the Ultra false-negative stratum is not
observable in routine data, so its weight `pi` enters as a scalar/grid (report a
central value + sensitivity analysis). Research data observed it directly, so
`pi = "observed"` uses a full 6-category Dirichlet there.

A fixed/grid `pi` is a scenario analysis, not probabilistic propagation of
uncertainty in `pi`. When defensible external information supports a probability
distribution for `pi`, a future API should accept `pi_draws` (or a Beta/logit-
normal prior) and integrate it jointly with the other draws.

## Aggregating burden distributions across settings

`pool_profiles()` sums counts — the crude aggregate — which flows through the
Dirichlet and propagates the combined sampling uncertainty. Pool only **within**
a setting type (e.g. two screening surveys); never across screening vs
programmatic (those are distinct reference targets). A random-effects /
Dirichlet-multinomial meta-analysis of case-mix distributions (carrying
between-dataset heterogeneity) is the principled next step and is not yet
implemented.

## Limitations

- Validity rests on conditional transportability of stratum-specific sensitivity
  (within a stratum, sensitivity is constant across settings) and on the
  stratifier adequately capturing burden variation — see METHOD.md §4.
- Index-to-index pairing of sensitivity and Dirichlet draws samples the product
  distribution and is coherent when those uncertainty sources are independent.
  If they share data or parameters, a joint model is required.
- `pool_profiles()` estimates an individual-weighted aggregate. It does not carry
  between-dataset heterogeneity or estimate the distribution for a new setting.
- A `pi` grid is structural sensitivity analysis; it does not average over a
  posterior for the unobserved Not-detected share.
- The stan engine needs `brms`; the shipped posterior draws let standardization
  run without re-fitting.
- Non-burden between-study confounders (reference standard, operators) are not
  removed by stratification. `plot_forest_stratum()` and
  `ssdta_meta_sensitivity(loo = TRUE)` surface within-stratum disagreement and
  single-study influence, but neither adjusts for them; quantitative bias
  analysis over residual per-stratum offsets is not yet implemented.
- The specification sweep is a robustness display, not model averaging — it
  shows what each arm gives, and does not combine them into one estimate.
- Specificity needs a separate spectrum model among people without the target
  condition; predictive values additionally require target prevalence.

Each of these is expanded, severity-ranked and paired with the planned fix in
**[METHOD.md](METHOD.md)** (§5 limitations, §6 open questions, §7 extensions).

> The shipped per-study data include studies that were unpublished when the
> package was assembled; publication status was re-verified and the release
> gate cleared 2026-08-17 — see [METHOD.md §8](METHOD.md#8-shipped-data--provenance-and-release-gate).

## Tests

```sh
PKG_ROOT="$(pwd)" Rscript tests/run_tests.R    # 53 checks incl. paper reproduction
```

Note the two harness paths: run from the package root it sources `R/` directly;
under `R CMD check` it runs against the *installed* package. Internal
(dot-prefixed) helpers used in a check must be added to the `getFromNamespace`
list at the top of the file, or the test passes locally and fails only on build.

## Building a distributable tarball

```sh
R CMD build SSDTA                                    # -> SSDTA_0.1.0.tar.gz (1.0 MB)
_R_CHECK_FORCE_SUGGESTS_=false R CMD check SSDTA_0.1.0.tar.gz
R CMD INSTALL SSDTA_0.1.0.tar.gz                     # what a downloader runs
```

`_R_CHECK_FORCE_SUGGESTS_=false` is needed only where `brms` is absent. Current
status: **check passes clean — 0 errors, 0 warnings, 0 notes.**

`.Rbuildignore` keeps working notes and `data-raw/` out of the CRAN-style
tarball (both still live in the git repo). Man pages are generated from the
roxygen comments in `R/*.R`:

```r
roxygen2::roxygenise(roclets = c("rd"))   # NOT "namespace"
```

`NAMESPACE` is hand-written and carries `importFrom` directives that no roxygen
tag in the source declares — regenerating it would silently drop the imports.
roxygen2 refuses to overwrite it on its own; keep it that way.

> The data release gate in
> [METHOD.md §8](METHOD.md#8-shipped-data--provenance-and-release-gate) was
> cleared 2026-08-17.

## Acknowledgments

The spectrum-standardization method, the anchored difference-in-differences
idea, and the overall statistical approach are Samuel Schumacher's own work.
Implementation — turning that analysis into this documented, tested R package
— was built with substantial AI assistance (Claude Code, Anthropic). See
[METHOD.md's Acknowledgments](METHOD.md#acknowledgments) for the fuller note,
including which review passes used AI and what that means for their
verification status.
