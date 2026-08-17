# SSDTA — method, assumptions, limitations and roadmap

Companion to [`README.md`](README.md), which covers *how to run* the package.
This file covers *what it assumes, what it does not do, and what is missing* —
the material that belongs beside any reported estimate.

**Method paper.** Schumacher SG (2026). *Spectrum-Standardised Diagnostic
Accuracy: A Framework for Adjusting Test Sensitivity for Case-Mix Variation
Across Clinical Settings.* World Health Organization, Geneva. Zenodo preprint,
CC BY-SA 4.0. **doi:[10.5281/zenodo.20033072](https://doi.org/10.5281/zenodo.20033072)**
· record: <https://zenodo.org/records/20033072>

`citation("SSDTA")` returns this plus the software entry.

---

## 1. The problem

Diagnostic sensitivity is not a property of a test. It is a joint function of
test performance and the *disease spectrum* of the population evaluated —
the spectrum effect (Ransohoff & Feinstein 1978; Mulherin & Miller 2002).
Research studies systematically over-represent patients with high mycobacterial
burden, who are easier to recruit and confirm, which inflates observed
sensitivity relative to routine diagnostic and screening contexts.

Two consequences follow: accuracy estimates from diagnostic trials predict field
performance poorly, and raw pooled sensitivities from non-comparative studies
cannot be fairly compared across tests evaluated in differently composed
populations.

## 2. The estimand

Estimate sensitivity *conditional on burden stratum*, then reweight by a target
population's burden distribution:

```
S_std = Σ_k  w_k · S_k
```

- `S_k` — sensitivity of the index test in burden stratum `k`, estimated from
  studies that applied a reference standard.
- `w_k` — the share of *true* cases in stratum `k` in the target population,
  estimated independently of the index test.

This is ordinary direct standardization applied to sensitivity. The algebra is
not where the difficulty lies; identification and transportability are (§4).

### The six-stratum structure

With Xpert Ultra semi-quantitative categories as the stratifier, five strata are
observable in routine data and one is not:

| Stratum | Observable in routine data? | Note |
|---|---|---|
| High, Medium, Low, Very low, Trace | yes | Ultra-positive |
| **Not-detected** | **no** | Ultra-negative but culture-positive — Ultra's own false negatives |

### Weight derivation and `pi`

Routine programmatic data identify the distribution *among detected cases*, so
the observed five-category distribution must be rescaled:

```
w_k = (1 - pi) · v_k     for k = 1..5
w_Not-detected = pi
```

where `v` is the observed detected-stratum composition and `pi` is the Ultra
false-negative rate in the target setting. `pi` is **not** observable in routine
data. In the package it enters as a scalar, a grid, or `"observed"` where a
study applied culture to Ultra-negatives (the research profile).

## 3. What the package implements — and what it does not

| | Status |
|---|---|
| Spectrum-standardized **sensitivity** | implemented |
| Specificity | **not implemented** — needs a separate spectrum model among people without the target condition |
| Predictive values | **not implemented** — additionally require target prevalence |
| ROC / threshold-varying measures | **not implemented** |
| Probabilistic `pi` | **not implemented** — `pi` is a scenario grid, not an integrated posterior (§5.2) |
| Hierarchical case-mix pooling | **not implemented** — `pool_profiles()` is crude count-summing (§5.5) |
| Cross-stratum partial pooling | optional via `engine = "stan_joint"`; per-stratum models are the default |
| Residual-transportability bias analysis (`delta_k`) | **not implemented** (§5.1) |
| Anchored difference-in-differences (indirect comparison, no stratum breakdown needed) | implemented — see §9 |

Despite the package name, v0.1 implements spectrum-standardized **sensitivity
only**. Do not apply the diseased-population weights to non-diseased
participants.

## 4. Validity conditions

Three conditions must hold for `S_std` to mean what it claims, plus one design
requirement on how the weights are obtained. Wording follows the method paper
(v1.1).

**C1 — Conditional transportability of stratum-specific sensitivity.** Within a
given burden stratum, sensitivity is constant across settings; differences in
test performance between populations are fully explained by differences in
burden distribution. Violated by differences in specimen quality, processing
delays, operator practice, or concurrent treatment.

**C2 — Adequacy of the burden measure as a stratifier.** The stratifier captures
biologically relevant variation in test performance at sufficient resolution.
Too coarse a measure leaves residual spectrum variation *within* strata. Where
the stratifier is itself a test output (Ultra semi-quant), it must reflect true
mycobacterial burden rather than test-specific measurement properties.

**C3 — Estimability of `pi` in the target setting.** The Not-detected weight
cannot be observed in routine data and must be modelled or borrowed, with
uncertainty quantified.

**D1 — Stratifier independence from the index test (design condition).** The
weights must come from data independent of the index test under evaluation.
Using the index test's own output to characterise the burden distribution is
circular. This is a design requirement, not a statistical assumption.

## 5. Known limitations

Ranked by severity. These should lead any write-up of a result from this
package; they are not footnotes.

**5.1 Conditional transportability is the identifying assumption (high).**
Standardization removes case-mix imbalance *at the resolution of the measured
stratifier* — it does not fully deconfound. Ultra semi-quant grade may leave
residual within-grade variation, and it captures nothing about specimen quality,
specimen type, treatment exposure, storage and transport, operator effects,
platform version, reference-standard differences, HIV status, or cavitation.
Diagnose with `plot_forest_stratum()` and `ssdta_meta_sensitivity(loo = TRUE)`;
neither adjusts for a violation it finds.

**5.2 `pi` is an unidentified target parameter, not sampling noise (high).**
A fixed or gridded `pi` is a *structural scenario analysis* and must not be
described as full uncertainty propagation. Report one explicitly justified
central scenario, a prespecified plausible range, the complete `S_std(pi)`
curve, and the source and population each anchor applies to.

**5.3 Verification and selection mechanisms can bias `S_k` (high).** Complete
verification, or a defensible correction model, is an input-quality condition
the package cannot check. Distinct from `pi`: `pi` concerns target composition,
verification bias concerns estimation of the conditional sensitivities.

**5.4 Independent per-stratum models are robust but inefficient (moderate).**
Sparse Very-low, Trace and Not-detected strata have unstable intercept and
heterogeneity posteriors — with k = 3–5 studies, tau is weakly identified, so
the between-study SD prior does real work there.

**5.5 Crude case-mix pooling answers one estimand only (moderate).**
`pool_profiles()` sums counts, which is correct for an *aggregate* distribution
across all sampled individuals. It is misleading for an *average setting* or a
*new setting*: a large programmatic dataset dominates smaller ones and
between-dataset heterogeneity disappears. Pool only within a setting type;
never across screening and programmatic.

**5.6 The anchored DiD fallback assumes more than its name suggests
(moderate).** It removes the anchor level-shift but leaves the burden-mix ×
gap-slope interaction, and its reference population is implicit. It is a
labelled fallback, not an equivalent of standardization.

**5.7 Independence of the two uncertainty sources.** Pairing one sensitivity
draw with one independently generated Dirichlet weight draw is a valid Monte
Carlo sample from the product posterior **only when those sources are
independent**. If they share data or parameters, a joint model is required. A
joint model is needed specifically when study-level composition, conditional
sensitivity, and shared random effects are estimated from the same studies and
modelled as correlated — plausible here, since case mix and sensitivity are
both study-level quantities, though the research application partially
factorizes into a stratum-frequency model and conditional binomial outcome
models. The `Dirichlet(counts + 1)` weight prior corresponds to a uniform
prior on composition; it is negligible for the large routine dataset and
material for small screening samples, so should be reported and included in
sensitivity analysis wherever a profile is sparse.

## 6. Open questions

Unresolved in the method itself:

- **Estimating `pi` in a target setting** — the central open problem. Routine
  data cannot identify it; borrowing across settings needs a defensible model.
  (`oq-spectrum-standardized-accuracy-pi-estimation`)
- **Partial verification and screening** — whether the `S_k` inputs are
  estimable and unbiased when source studies use CAD pre-screening or apply the
  reference standard to only some participants.
  (`oq-spectrum-standardization-partial-verification-screening`)

## 7. Planned extensions

Recommended order:

1. **`pi_draws` interface** — accept a Beta/logit-normal posterior for `pi` and
   integrate it jointly, with reporting that keeps scenario-grid and
   integrated-posterior results visibly separate. Starting model: with `y` Ultra
   false negatives among `n` true cases, `pi ~ Beta(a+y, b+n-y)`. A hierarchical
   logistic model is preferable when borrowing across settings. Do not assign a
   generic default prior without defining the target setting and evidence base.
2. **Leave-one-study-out predictive validation** — partly done:
   `ssdta_meta_sensitivity(loo = TRUE)` refits and tabulates. Still missing:
   predicting a *held-out study's marginal sensitivity from its own case mix*,
   which is the actual test of C1.
3. **Residual-transportability harness** — shift target stratum sensitivities by
   plausible offsets `delta_k` (logit scale) and report a tornado plot.
4. **Hierarchical case-mix model** — Dirichlet-multinomial or logistic-normal
   multinomial over datasets, with an explicit aggregate / mean-setting /
   new-setting estimand, replacing crude pooling (§5.5).
5. **Ordinal partial pooling as a v2 sensitivity analysis** — random-walk prior
   on stratum logits, or a monotone mean curve with regularized deviations.
   Compare against independent models; avoid a hard isotonic constraint as the
   sole analysis.
6. **A non-TB worked example** before claiming broad generic performance.
   Published stratum-specific data already support application in principle for
   COVID-19 rapid antigen testing, colorectal cancer screening, and malaria
   diagnostics.

Also outstanding as *software* rather than method work: `man/` is empty, so
`?ssdta_standardize` returns nothing — the roxygen comments in `R/*.R` have
never been rendered to `.Rd`. `NAMESPACE` is deliberately hand-written, so
running `roxygen2::roxygenise()` would need care.

## 8. Shipped data — provenance and release gate

See [`README.md`](README.md#data) for the object-by-object listing and `R/data.R`
for per-object documentation. Provenance in brief:

| Object | Origin |
|---|---|
| `npoc_by_ultra_stratum` | Per-study TP/FN by Ultra stratum, 5 studies, WHO TB Dx GDG PICO 1 systematic review (Horne et al. 2025), Table 3 |
| `npoc_studies` | Characteristics of those same 5 studies — Horne et al. 2025 Table 1 and the "Characteristics of included studies" narrative; `n_culture_pos`/`n_ultra_pos` derived from Table 3, citations/DOIs verified against Crossref |
| `npoc_posterior_draws`, `npoc_stratum_summary` | The paper's brms/Stan fit to the above; Union 2026 abstract UNIONCONF2026:3197 |
| `burden_profiles$research` | Pooled Ultra-positive culture-confirmed cases from the same 5 GDG PICO1 studies (n=507, plus 40 observed Not-detected) |
| `burden_profiles$routine_pakistan` | GxAlert National TB Reference Laboratory reports 2019 + 2021, pooled (n=88,560) |
| `burden_profiles$screening_kendall` | Kendall et al., community screening, Uganda (n=89) |
| `burden_profiles$screening_pakistan` | Pluslife/MiniDock screening evaluation (n=118) |

> **Publication status.** Two of the five source studies are peer-reviewed
> journal publications — Steadman (eBioMedicine 2025;121:105991) and Yerlikaya
> (NEJM 2026;394:1710–1722). The other three are preprints, not yet
> peer-reviewed — Hartati (medRxiv), Mbuli (medRxiv), and Rockman (Research
> Square). `npoc_studies$pub_status` and `npoc_studies$doi` record this per
> study.
>
> The shipped counts predate publication and differ slightly from the
> now-published papers — Yerlikaya's NEJM paper reports 1380 enrolled and 226
> culture-confirmed, against 1300 enrolled / 223-with-Ultra-grade used here.
> Expect small discrepancies when cross-checking against the published papers
> directly; that reconciliation has not yet been done.
>
> South African NHLS programmatic data (~11.4M Ultra tests) is larger than the
> Pakistan routine source but was not used: the full five-category breakdown is
> unpublished and data permissions did not allow it.

## 9. The anchored difference-in-differences fallback — full method

`anchored_did()` and the `did_*` functions cover the case standardization
cannot: comparing two tests when only one has a stratum breakdown. This is a
distinct estimator from `S_std = Σ w_k S_k`, with its own idea, assumption, and
literature position, laid out in full here rather than compressed into the
README's usage note.

### 9.1 The idea

When two tests, A and B, have never been compared head-to-head in the same
population against the same reference standard, but each has separately been
compared against a common third "anchor" test C — in different studies,
different populations, possibly different reference standards — estimate the
A-vs-B difference indirectly using a **difference-in-differences (DiD)**
contrast, borrowed from the econometric method of the same name.

This is structurally the diagnostic-accuracy analogue of the **Bucher method**
for anchored indirect treatment comparisons in network meta-analysis (Bucher et
al. 1997), applied to sensitivity rather than a treatment effect.

### 9.2 The estimator

Let `Se_P(T)` denote the sensitivity of test `T` as estimated in study or
population `P`, using that study's own reference standard.

**Naive (confounded) indirect comparison** — subtracting two separately-pooled
estimates:

```
Δ_naive = Se_PB(B) - Se_PA(A)
```

Biased whenever `P_A ≠ P_B` in spectrum, prevalence, or reference standard.

**Anchored DiD** — using common comparator `C`, measured within each study
alongside the test of interest:

```
Δ_A = Se_PA(A) - Se_PA(C)
Δ_B = Se_PB(B) - Se_PB(C)
Δ_DiD = Δ_B - Δ_A = [Se_PB(B) - Se_PA(A)] - [Se_PB(C) - Se_PA(C)]
```

The second bracket is the between-study drift in the anchor's own sensitivity.
The DiD estimator subtracts it out, on the assumption that this drift also
captures the population/reference-standard effect on the A-vs-C and B-vs-C
gaps — the diagnostic-accuracy analogue of the "parallel trends" assumption in
econometric DiD.

The package's worked example (`pooling_vs_npoc_example()`,
[`inst/examples/example_did_fallback.R`](inst/examples/example_did_fallback.R))
anchors NPOC and pooled Xpert Ultra on individual Xpert Ultra, giving
`Δ_DiD = -4.4` percentage points (NPOC sensitivity below pooled Ultra), with
four variance methods (Wald, meta-analytic, bootstrap, Bayesian conjugate
Beta-Binomial with prior-sensitivity analysis) agreeing closely.

### 9.3 Why a separate estimator, when standardization exists

Naive indirect comparison — subtracting two separately-pooled global estimates
of two tests — is confounded whenever the studies behind them differ in case
mix or reference standard. Anchored DiD is a candidate correction: instead of
anchoring to a *global* pooled estimate of the comparator, it anchors to the
comparator as measured *within each specific study*, which should remove
exactly that confounding, provided the parallel-trends-type assumption holds.

Its practical role in this package is for the case standardization cannot
handle: when one of the two tests has no stratum breakdown at all (e.g. pooled
Ultra is reported only as an aggregate sensitivity, with no per-specimen
semi-quant grade), so there is no `S_k` to reweight. See §5.6 for what this
estimator does and does not deconfound.

### 9.4 Position in the literature

No exact precedent was found for a raw-scale (not logit/odds-ratio) sensitivity
difference-in-differences through a common anchor test, explicitly framed as
Bucher-analogous with a parallel-trends identifying assumption, in the
diagnostic test accuracy literature.

**Not novel** — well-established, and not claimed as new here: the network
topology A–C, B–C ⇒ A–B and the broader direct/indirect-evidence framework
(standard in diagnostic network meta-analysis); the recognition that naive
between-study test comparison is confounded; using a common anchor test to
compare two tests indirectly; the algebra itself, which is Bucher-on-risk-
differences applied to sensitivity.

**Novel or underexploited:**

1. **Estimand and framing.** The raw percentage-point sensitivity difference as
   the target, with an explicit Bucher / difference-in-differences framing and
   a parallel-trends-style identifying assumption stated on the
   sensitivity-vs-anchor gap. The formal DTA-network literature works on the
   logit / odds-ratio / relative-diagnostic-odds-ratio scale, or full
   sensitivity-specificity (SROC) surfaces — a raw-difference,
   parallel-trends-labelled estimator for DTA is essentially absent from it.
2. **A "when is the simple estimator legitimate?" argument.** The methods
   literature's arm-based, bivariate machinery is motivated by threshold
   variation across tests, an SROC curve defined by accuracy levels rather than
   differences, and both sensitivity and specificity mattering jointly. A
   fixed-threshold molecular assay comparison removes each of those premises:
   fixed manufacturer thresholds (no SROC to reconstruct), near-perfect
   homogeneous specificity with sensitivity as the outcome of interest (the
   bivariate object collapses to one operating point), and — in the TB
   application — a single dominant, measurable driver of sensitivity
   heterogeneity (bacillary burden). In that collapsed regime, a transparent
   raw-scale anchored difference is arguably the natural estimator rather than
   a crude shortcut.
3. **A checkable identifying assumption.** General DTA-network methods treat
   case-mix as an unmodelled nuisance and can only assume the anchor's gap is
   exchangeable across studies. When the dominant driver of that gap is known
   and measurable (bacillary burden), the assumption becomes checkable and
   adjustable: both index-minus-anchor gaps can be standardized to a common
   burden distribution via §1–§4's method before subtracting.

**Caveats:**

- This is a special case, not a new estimator class: it is Bucher-on-risk-
  differences under stated conditions, and the contribution is the conditions
  and the framing, not the algebra. General diagnostic network meta-analysis
  can already produce this comparison as a byproduct; the case for the simpler
  estimator is transparency and communicability, not capability.
- A published warning is relevant here: an individual-patient-data study
  found that adjusting for indirectness did **not** successfully remove bias
  in a comparative test-accuracy case, even with full IPD (Wang et al. 2015,
  *J Clin Epidemiol*, ovarian-reserve testing). Its failure modes — variable
  thresholds on continuous tests, heterogeneous prognostic reference standards,
  multi-dimensional unmeasured confounding, a genuinely bivariate target — are
  largely the ones a fixed-threshold, high-specificity, sensitivity-only
  molecular TB comparison removes, but it is a real cautionary data point that
  adjustment can fail even with good data.
- Variance estimation for this two-step estimator under mixed paired/unpaired
  anchor reporting is genuinely open, not just a write-up detail — see §5.6 and
  §5.7.

**Bottom line:** not "a new way to compare diagnostic tests indirectly," but
the conditions under which a transparent anchored sensitivity
difference-in-differences is a legitimate and preferable alternative to full
diagnostic network meta-analysis — conditions met by fixed-threshold,
high-specificity molecular TB assays — with the identifying assumption made
checkable through spectrum standardization.

### 9.5 Open questions specific to this estimator

- **Assumption strength.** Requires that the population/reference-standard
  effect on the test-vs-anchor gap is constant across the two studies (no
  interaction). Untestable with only two studies and no direct A-vs-B data.
- **Variance under mixed paired/unpaired anchor reporting.** The package
  currently assumes within-pair independence for the bootstrap and Bayesian
  engines unless raw paired 2×2 tables are supplied; the general case with
  correlated anchors needs further work.
- **Relationship to full diagnostic network meta-analysis.** DiD with a single
  anchor is essentially a two-study, three-node diagnostic network; whether
  existing DNMA methodology already subsumes this as a special case has not
  been formally checked.

## 10. Status of this document

Sections 1–4 and 8 draw directly on the Zenodo preprint and the WHO GDG PICO1
reports. §9.1–9.2 follow the standard DiD/Bucher formalism. Sections 5–7 and
§9.3–9.5 are argued methodological positions from an internal review pass
(2026-07) and two independent literature searches (2026-07-08) — reasoned, but
**not yet independently checked against the primary methods literature by a
third-party reviewer**. Treat them as a starting point for scrutiny, not a
settled verdict. Exact wording of the validity conditions in §4 should be
checked against manuscript v1.1 before quotation. Corrections and disagreement
on the §9 novelty claims are exactly the kind of review this document is
published to invite.

## Acknowledgments

The spectrum-standardization method, the anchored difference-in-differences
idea, and the overall statistical approach in this package are Samuel
Schumacher's own work. Turning that analysis into a documented, tested R
package — the estimation engines, diagnostics, standardization and comparison
code, visualizations, self-tests, and this documentation — was built with
substantial AI assistance (Claude Code, Anthropic). Statistical review passes
that informed §5 and §9.3–9.5 also used AI tools; their conclusions carry the
verification caveat above, independent of who or what produced the code.

Last updated 2026-08-17.
