# =============================================================================
# Example 1 -- Spectrum standardization (the main SSDTA use case)
#
# NPOC-NAAT sensitivity standardized to research, routine (Pakistan GxAlert),
# and screening (Kendall / new Pakistan / pooled) populations, using the
# stratum-specific posterior sensitivities from the paper's Bayesian
# meta-analysis. Reproduces the paper's numbers and adds the new screening data.
# =============================================================================

library(SSDTA)
library(ggplot2)

## ---- 1. Standardize across settings ---------------------------------------
## Research: pi observed (research studies applied culture to Ultra-negatives).
## Routine/Screening: pi unknown -> report at a central value + sensitivity grid.
research  <- ssdta_standardize(npoc_posterior_draws, burden_profiles$research,
                               pi = "observed")
routine   <- ssdta_standardize(npoc_posterior_draws, burden_profiles$routine_pakistan,
                               pi = 0.20)
screening <- ssdta_standardize(npoc_posterior_draws, burden_profiles$screening_kendall,
                               pi = 0.30)
print(research)    # 80.9% (68.1, 86.1)
print(routine)     # 66.8% (57.5, 73.7)
print(screening)   # 53.3% (43.6, 63.7)

## ---- 2. New Pakistan screening data + crude pool with Kendall --------------
scr_pak <- ssdta_standardize(npoc_posterior_draws, burden_profiles$screening_pakistan,
                             pi = 0.30)
print(scr_pak)     # 50.9% -- more low-burden-skewed than Kendall

scr_pool_profile <- pool_profiles(burden_profiles$screening_kendall,
                                  burden_profiles$screening_pakistan,
                                  setting = "Screening",
                                  source  = "Kendall + Pakistan (crude pool, n=207)")
scr_pooled <- ssdta_standardize(npoc_posterior_draws, scr_pool_profile, pi = 0.30)
print(scr_pooled)  # 51.8%

## ---- 3. pi sensitivity grids ----------------------------------------------
grids <- list(
  routine   = ssdta_standardize(npoc_posterior_draws, burden_profiles$routine_pakistan,
                                pi = seq(0.05, 0.40, 0.05)),
  screening = ssdta_standardize(npoc_posterior_draws, scr_pool_profile,
                                pi = seq(0.05, 0.55, 0.05))
)

## ---- 4. (optional) re-fit the stratum meta-analysis from raw counts --------
## Default engine "stan" reproduces the paper (needs brms). "glmm" and
## "conjugate" are dependency-light alternatives that track it closely.
meta_glmm <- ssdta_meta(npoc_by_ultra_stratum, engine = "glmm")
print(meta_glmm$summary)

## ---- 5. Figures ------------------------------------------------------------
p_casemix   <- plot_casemix(burden_profiles)                    # spectrum across settings
p_stratum   <- plot_stratum_sensitivity(npoc_posterior_draws)   # S_k forest
p_pi        <- plot_pi_sensitivity(grids)                       # S_std vs pi
p_std       <- plot_standardized(list(research, routine, scr_pooled))
p_mechanism <- plot_spectrum_mechanism(
  npoc_posterior_draws,
  burden_profiles[c("research", "routine_pakistan", "screening_kendall")],
  pi = list("observed", 0.20, 0.30))                            # w_k x S_k -> S_std

outdir <- file.path(tempdir(), "SSDTA_plots"); dir.create(outdir, showWarnings = FALSE)
for (nm in c("p_casemix","p_stratum","p_pi","p_std","p_mechanism"))
  ggsave(file.path(outdir, paste0(nm, ".png")), get(nm), width = 7.5, height = 4.6, dpi = 110)
cat("Figures written to:", outdir, "\n")
