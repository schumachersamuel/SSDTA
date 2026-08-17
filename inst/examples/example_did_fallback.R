# =============================================================================
# Example 2 -- Anchored difference-in-differences fallback
#
# When you WANT the spectrum-standardized comparison but one test lacks a
# semi-quant breakdown, use the anchored DiD instead. Here: NPOC vs pooled
# Ultra. Pooled Ultra has NO Ultra semi-quant breakdown (you cannot stratify a
# pooled result by an individual-specimen grade), so ssdta_standardize /
# ssdta_compare are not available for this pair -- the anchored DiD, anchored on
# individual Ultra, partially deconfounds instead.
#
# See the anchored-DiD concept page for why this only PARTIALLY deconfounds
# (it removes the anchor level-shift but leaves the burden-mix x gap-slope
# interaction), whereas the standardized comparison in Example 1 removes it fully.
# =============================================================================

library(SSDTA)

ex <- pooling_vs_npoc_example()

## Headline: contrast form (reported GDG difference CIs used directly)
xc <- did_input(delta_A = ex$contrasts$delta_A, ci_A = ex$contrasts$ci_A,
                delta_B = ex$contrasts$delta_B, ci_B = ex$contrasts$ci_B,
                labels  = ex$contrasts$labels)
did_wald(xc)          # Delta_DiD = -4.44 pp, 95% CI (-10.0, +1.2)

## All four interval methods on the component form
xs  <- did_input(sens = ex$components$sens, labels = ex$components$labels)
fit <- anchored_did(xs, methods = c("wald", "bootstrap", "bayes"), seed = 1)
print(fit)

## Meta-analytic pooling of the real per-study paired contrasts
xm <- did_input(delta_A = ex$contrasts$delta_A, delta_B = ex$contrasts$delta_B,
                delta_A_i = ex$per_study$delta_A_i, se_A_i = ex$per_study$se_A_i,
                delta_B_i = ex$per_study$delta_B_i, se_B_i = ex$per_study$se_B_i)
did_meta(xm)

## Plots
plot_parallel(xs); plot_methods(fit); plot_posterior(fit)
