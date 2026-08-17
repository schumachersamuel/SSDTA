# =============================================================================
# Example 3 -- Diagnostics for the stratum-specific meta-analysis (S_k)
#
# Standardization inherits every weakness of the S_k it reweights. Two checks
# belong with every reported S_std:
#
#   (a) Do the studies agree *within* a stratum? That is the conditional-
#       transportability condition, and a study-level forest is where a
#       violation is visible.
#   (b) Does the answer survive the prior and the study set? With k = 3-5
#       studies per stratum, tau is weakly identified, so the between-study SD
#       prior does real work in the sparse strata.
#
# Reproduces the diagnostics from the paper's Bayesian meta-analysis script.
# =============================================================================

library(SSDTA)
library(ggplot2)

## ---- 1. Study-level forest, with the paper's pooled posterior ---------------
## Circles: each study's tp/n with a Wilson interval. Diamond: the pooled S_k.
## Watch the sparse strata -- Very low, Trace, Not-detected -- where a single
## study can carry the pooled estimate.
p_forest <- plot_forest_stratum(npoc_by_ultra_stratum,
                                pooled = npoc_posterior_draws)

## ---- 2. Prior- and specification-sensitivity sweep --------------------------
## Defaults on the "stan" engine reproduce the paper's analyses:
##   Primary  N(0,1.5) intercept, Exponential(1) tau
##   8a       N(0,5) intercept
##   8b       joint model, shared tau across strata
##   8d       Logistic(0,1) intercept  [= Uniform(0,1) on the probability scale]
##   8e       Half-t(3,0,2.5) tau      [Gelman's hierarchical-scale prior]
## `loo = TRUE` adds one leave-one-study-out arm per study (the paper's 8c
## generalised), which is the influence check the package review asked for.
if (requireNamespace("brms", quietly = TRUE)) {
  sweep <- ssdta_meta_sensitivity(npoc_by_ultra_stratum, engine = "stan", loo = TRUE)
} else {
  # dependency-light stand-in: priors do nothing here, so sweep the study set
  # and the pooling assumption instead of the priors.
  sweep <- ssdta_meta_sensitivity(
    npoc_by_ultra_stratum, engine = "glmm", loo = TRUE,
    specs = list(meta_spec("Primary (RE)"),
                 meta_spec("Complete pooling", engine = "conjugate")))
}
print(sweep)                     # wide "median (lo-hi)" comparison table
p_sweep <- plot_meta_sensitivity(sweep)

## The wide table is the paper's table_sensitivity_analyses.csv:
write.csv(as.data.frame(sweep),
          file.path(tempdir(), "table_sensitivity_analyses.csv"), row.names = FALSE)

## ---- 3. Does the spread in S_k matter downstream? ---------------------------
## The question is never "did S_k move" but "did S_std move". Push each arm's
## draws through the same target distribution and compare.
std_by_arm <- lapply(names(sweep$fits), function(nm)
  ssdta_standardize(sweep$fits[[nm]]$draws, burden_profiles$routine_pakistan,
                    pi = 0.20, seed = 442))
names(std_by_arm) <- names(sweep$fits)
cat("\nStandardized sensitivity (routine, pi = 0.20) by specification:\n")
for (nm in names(std_by_arm))
  cat(sprintf("  %-26s %.1f%%  (%.1f, %.1f)\n", nm,
              100 * std_by_arm[[nm]]$summary$median,
              100 * std_by_arm[[nm]]$summary$lo,
              100 * std_by_arm[[nm]]$summary$hi))
p_std_by_arm <- plot_standardized(std_by_arm)

## ---- 4. Save ----------------------------------------------------------------
outdir <- file.path(tempdir(), "SSDTA_diagnostics"); dir.create(outdir, showWarnings = FALSE)
ggsave(file.path(outdir, "fig_forest.png"),            p_forest,     width = 10, height = 9, dpi = 130)
ggsave(file.path(outdir, "fig_meta_sensitivity.png"),  p_sweep,      width = 11, height = 6, dpi = 130)
ggsave(file.path(outdir, "fig_std_by_spec.png"),       p_std_by_arm, width = 8,  height = 5, dpi = 130)
cat("\nFigures written to:", outdir, "\n")
