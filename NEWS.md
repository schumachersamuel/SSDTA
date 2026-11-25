# SSDTA 0.1.1 (2026-11-23)

- `ssdta_meta(engine = "stan_joint")` (and `ssdta_meta_sensitivity()` with a
  joint-model spec) no longer stops when a stratum name contains a non-syntactic
  character such as "Not-detected". The model is now fitted on internal labels and
  mapped back by position. Before, brms renamed the coefficient to
  `b_stratumNotMdetected` and the fit stopped with "Could not map strata".
- A `control` list passed to the joint model is now merged over the default
  (`adapt_delta = 0.95`) instead of failing with "formal argument 'control'
  matched by multiple actual arguments".
- New optional regression test (set `SSDTA_TEST_STAN=1`).

# SSDTA 0.1.0

- First release.
