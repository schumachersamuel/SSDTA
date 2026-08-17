# SSDTA self-tests. Run from package root:  PKG_ROOT="$(pwd)" Rscript tests/run_tests.R
root <- Sys.getenv("PKG_ROOT", unset = getwd())
suppressMessages(library(ggplot2))
if (dir.exists(file.path(root, "R")) && dir.exists(file.path(root, "data"))) {
  for (f in list.files(file.path(root, "R"), pattern = "\\.R$", full.names = TRUE)) source(f)
  for (d in list.files(file.path(root, "data"), full.names = TRUE)) load(d, envir = globalenv())
} else {
  suppressPackageStartupMessages(library(SSDTA))
  # Internal helpers exercised directly below. Keep this list in step with the
  # dot-prefixed functions used in the checks -- omitting one fails only under
  # R CMD check (installed package), not when sourcing R/ during development.
  for (fn in c(".weight_draws", ".wilson_ci"))
    assign(fn, getFromNamespace(fn, "SSDTA"))
}

ok <- 0L; fail <- 0L
check <- function(cond, msg) {
  if (isTRUE(cond)) { ok <<- ok + 1L; cat("  ok  ", msg, "\n") }
  else { fail <<- fail + 1L; cat("  FAIL", msg, "\n") }
}
near <- function(a, b, tol = 0.005) abs(a - b) < tol
errors <- function(expr, pattern = NULL) {
  msg <- tryCatch({ force(expr); NULL }, error = conditionMessage)
  !is.null(msg) && (is.null(pattern) || grepl(pattern, msg, fixed = TRUE))
}

# ---- data integrity -------------------------------------------------------
check(nrow(npoc_posterior_draws) == 48000 &&
      length(levels(npoc_posterior_draws$category)) == 6, "posterior draws: 48000 rows, 6 strata")
check(sum(subset(npoc_by_ultra_stratum, stratum != "Not-detected")[, c("tp","fn")]) == 507,
      "per-study positive strata total = 507 (matches research burden counts)")
check(all(c("research","routine_pakistan","screening_kendall","screening_pakistan") %in%
          names(burden_profiles)), "4 burden profiles incl. new Pakistan screening")
check(sum(burden_profiles$routine_pakistan$counts_detected) == 88560, "Pakistan routine n = 88,560")
check(sum(burden_profiles$screening_pakistan$counts_detected) == 118, "new Pakistan screening n = 118")

# ---- study table reconciles with the count table ---------------------------
check(nrow(npoc_studies) == 5 &&
      identical(sort(npoc_studies$study),
                sort(unique(as.character(npoc_by_ultra_stratum$study)))),
      "npoc_studies: 5 rows, study names match npoc_by_ultra_stratum")
.per_study <- function(d) {
  v <- tapply(d$tp + d$fn, factor(d$study, levels = npoc_studies$study), sum)
  as.integer(ifelse(is.na(v), 0L, v))
}
check(identical(npoc_studies$n_culture_pos, .per_study(npoc_by_ultra_stratum)) &&
      sum(npoc_studies$n_culture_pos) == 547L,
      "npoc_studies n_culture_pos = per-study TP+FN over all strata (total 547)")
check(identical(npoc_studies$n_ultra_pos,
                .per_study(subset(npoc_by_ultra_stratum, stratum != "Not-detected"))) &&
      sum(npoc_studies$n_ultra_pos) == 507L,
      "npoc_studies n_ultra_pos = positive-strata TP+FN (total 507)")
check(all(npoc_studies$n_enrolled >= npoc_studies$n_culture_pos) &&
      all(npoc_studies$n_culture_pos >= npoc_studies$n_ultra_pos),
      "npoc_studies: n_enrolled >= n_culture_pos >= n_ultra_pos for every study")
check(sum(npoc_studies$n_enrolled) == 3349L &&
      min(npoc_studies$n_enrolled) == 80L && max(npoc_studies$n_enrolled) == 1300L,
      "npoc_studies enrolment: 3349 total, range 80-1300 (PICO1 Table 1)")
check(all(!is.na(npoc_studies$doi)) &&
      all(grepl("^10\\.[0-9]{4,5}/", npoc_studies$doi)) &&
      !anyDuplicated(npoc_studies$doi),
      "npoc_studies: every study has a distinct, well-formed DOI")
check(all(npoc_studies$pub_status %in% c("Published", "Preprint")) &&
      sum(npoc_studies$pub_status == "Published") == 2L,
      "npoc_studies: publication status set; 2 of 5 fully published (release gate)")
check(length(unique(npoc_studies$index_test)) == 1L &&
      length(unique(npoc_studies$comparator)) == 1L,
      "npoc_studies: one index test and one comparator across all 5 studies")

# ---- standardization reproduces the paper ---------------------------------
r  <- ssdta_standardize(npoc_posterior_draws, burden_profiles$research, pi = "observed", seed = 42)
rt <- ssdta_standardize(npoc_posterior_draws, burden_profiles$routine_pakistan, pi = 0.20, seed = 442)
sc <- ssdta_standardize(npoc_posterior_draws, burden_profiles$screening_kendall, pi = 0.30, seed = 342)
check(near(r$summary$median,  0.8095), "Research standardized = 80.9% (paper: 0.809)")
check(near(rt$summary$median, 0.6680), "Routine pi=0.20 = 66.8% (paper: 0.668)")
check(near(sc$summary$median, 0.5325), "Screening pi=0.30 = 53.3% (paper: 0.533)")
check(near(r$summary$lo, 0.681, 0.002) && near(r$summary$hi, 0.861, 0.002) &&
      near(rt$summary$lo, 0.575, 0.002) && near(rt$summary$hi, 0.737, 0.002) &&
      near(sc$summary$lo, 0.436, 0.002) && near(sc$summary$hi, 0.637, 0.002),
      "paper headline 95% intervals reproduced")
check(r$summary$lo < r$summary$median && r$summary$median < r$summary$hi, "CrI brackets the median")

# local seeding must not alter the caller's random stream
set.seed(818); expected_rng <- runif(3)
set.seed(818); invisible(ssdta_standardize(npoc_posterior_draws,
  burden_profiles$routine_pakistan, pi = 0.2, seed = 99)); observed_rng <- runif(3)
check(identical(expected_rng, observed_rng), "standardization preserves caller RNG state")

# ---- monotone in pi, and lower for more low-burden populations ------------
g <- ssdta_standardize(npoc_posterior_draws, burden_profiles$routine_pakistan, pi = seq(0.05, 0.4, 0.05))
check(all(diff(g$median) < 0), "standardized sensitivity decreases monotonically in pi")
sp <- ssdta_standardize(npoc_posterior_draws, burden_profiles$screening_pakistan, pi = 0.30, seed = 342)
check(sp$summary$median < sc$summary$median,
      "more low-burden Pakistan screening < Kendall screening at same pi")

# ---- crude pooling sums counts --------------------------------------------
pooled <- pool_profiles(burden_profiles$screening_kendall, burden_profiles$screening_pakistan,
                        setting = "Screening")
check(pooled$counts_detected[["Trace"]] == 23 + 44, "pooled Trace count = 23 + 44")
check(pooled$n == 89 + 118, "pooled n = 207")
sp_pool <- ssdta_standardize(npoc_posterior_draws, pooled, pi = 0.30, seed = 342)
check(sp_pool$summary$median > sp$summary$median &&
      sp_pool$summary$median < sc$summary$median, "pooled screening lies between its two inputs")

# ---- weight draws sum to 1, respect pi ------------------------------------
w <- .weight_draws(burden_profiles$screening_kendall$counts_detected, pi = 0.3, n_draws = 500)
check(all(abs(rowSums(w) - 1) < 1e-9), "6-stratum weight draws sum to 1")
check(all(abs(w[, "Not-detected"] - 0.3) < 1e-9), "Not-detected weight fixed at pi")
at_one <- ssdta_standardize(npoc_posterior_draws, burden_profiles$screening_kendall,
                            pi = 1, seed = 4)
nd_draws <- subset(npoc_posterior_draws, category == "Not-detected")$sens_pooled
check(max(abs(at_one$draws - nd_draws)) < 1e-12, "pi=1 reduces exactly to Not-detected sensitivity")
check(errors(burden_profile(c(High = 0, Low = 0)), "positive total"),
      "all-zero burden profile fails clearly")
bad_mat <- matrix(0.5, nrow = 10, ncol = 5,
                  dimnames = list(NULL, c("High","Medium","Low","Very low","Trace")))
check(errors(ssdta_standardize(bad_mat, burden_profiles$screening_kendall, 0.3),
             "missing strata"), "missing sensitivity stratum fails clearly")

# ---- meta engines track the paper's Stan S_k ------------------------------
mg <- ssdta_meta(npoc_by_ultra_stratum, engine = "glmm", n_draws = 3000, seed = 1)
mc <- ssdta_meta(npoc_by_ultra_stratum, engine = "conjugate", n_draws = 3000, seed = 1)
paper_high <- 0.94; paper_low <- 0.89
check(near(mg$summary$Sensitivity[mg$summary$Category=="High"], paper_high, 0.03),
      "glmm High ~ 94% (paper Stan)")
check(near(mg$summary$Sensitivity[mg$summary$Category=="Low"], paper_low, 0.03),
      "glmm Low ~ 89% (paper Stan)")
check(nrow(mc$draws) == 6 * 3000 && all(levels(mc$draws$category) ==
      c("High","Medium","Low","Very low","Trace","Not-detected")), "conjugate draws shape/strata ok")
# re-fitted S_k feed standardization
sc2 <- ssdta_standardize(mg$draws, burden_profiles$screening_kendall, pi = 0.30, seed = 342)
check(is.finite(sc2$summary$median) && sc2$summary$median > 0.3 && sc2$summary$median < 0.7,
      "re-fitted (glmm) S_k feed ssdta_standardize to a plausible value")

check(inherits(mg, "ssdta_meta") && grepl("engine: glmm", paste(capture.output(print(mg)), collapse = " ")),
      "ssdta_meta returns a classed object with a readable print method")
check(errors(ssdta_meta(npoc_by_ultra_stratum, engine = "stan_joint", n_draws = 10), "brms") ||
      requireNamespace("brms", quietly = TRUE),
      "stan_joint engine names its missing dependency (or brms is installed)")

# ---- Wilson intervals ------------------------------------------------------
w1 <- .wilson_ci(1, 2)
check(near(w1[, "lower"], 0.0945, 1e-3) && near(w1[, "upper"], 0.9055, 1e-3),
      "Wilson CI for 1/2 = (9.5%, 90.6%)")
wb <- .wilson_ci(c(0, 4), c(4, 4))
check(near(wb[1, "lower"], 0, 1e-9) && near(wb[2, "upper"], 1, 1e-9) &&
      wb[1, "upper"] > 0.1 && wb[2, "lower"] < 0.9 &&
      all(wb >= 0) && all(wb <= 1),
      "Wilson CI stays in [0,1] and keeps width at the boundaries")

# ---- study-level forest plot ----------------------------------------------
fp <- plot_forest_stratum(npoc_by_ultra_stratum, pooled = npoc_posterior_draws)
n_cells <- sum(npoc_by_ultra_stratum$tp + npoc_by_ultra_stratum$fn > 0)
check(nrow(fp$data) == n_cells + 6, "forest rows = non-empty study cells + 6 pooled")
check(sum(fp$data$type == "Pooled") == 6 && all(is.na(fp$data$count_label[fp$data$type == "Pooled"])),
      "forest carries one pooled row per stratum, without a count label")
# first level of a discrete axis renders at the bottom: pooled must lead each facet
lv <- levels(fp$data$y_key)
check(all(grepl("^01\\|Pooled$|^0", lv[1])) && grepl("Pooled", lv[1]) &&
      grepl("Pooled", lv[grep("^02\\|", lv)][1]),
      "forest puts the pooled diamond at the bottom of every facet")
studies_high <- fp$data[fp$data$type == "Study" & fp$data$stratum == "High", ]
check(identical(order(studies_high$n),
                order(match(as.character(studies_high$y_key), lv))),
      "forest orders studies by n within a stratum")
check(nrow(plot_forest_stratum(npoc_by_ultra_stratum)$data) == n_cells,
      "forest without `pooled` shows study rows only")
check(errors(plot_forest_stratum(npoc_by_ultra_stratum, pooled = list(nope = 1)),
             "must be an ssdta_meta"), "forest rejects an unusable `pooled`")
check(nrow(plot_forest_stratum(npoc_by_ultra_stratum, pooled = mg)$data) == n_cells + 6 &&
      nrow(plot_forest_stratum(npoc_by_ultra_stratum, pooled = mg$summary)$data) == n_cells + 6,
      "forest accepts ssdta_meta results and plain summary tables")

# ---- specification sweep ---------------------------------------------------
sw <- suppressWarnings(ssdta_meta_sensitivity(npoc_by_ultra_stratum, engine = "glmm",
                                              n_draws = 800, loo = TRUE, seed = 1))
n_studies <- length(unique(npoc_by_ultra_stratum$study))
check(length(sw$fits) == n_studies + 1 && nlevels(sw$table$analysis) == n_studies + 1,
      "loo sweep = Primary + one arm per study")
check(nrow(sw$table) == (n_studies + 1) * 6, "sweep table has one row per arm x stratum")
wide <- as.data.frame(sw)
check(ncol(wide) == n_studies + 2 && nrow(wide) == 6 && !anyNA(wide),
      "wide comparison table is complete (6 strata x arms + Stratum column)")
check(all(sw$table$lo <= sw$table$sensitivity & sw$table$sensitivity <= sw$table$hi),
      "sweep intervals bracket every point estimate")
# leaving out a study must actually change the fit
prim <- sw$table[sw$table$analysis == "Primary" & sw$table$stratum == "Very low", "sensitivity"]
drop <- sw$table[sw$table$analysis == "Excl. Hartati 2025" & sw$table$stratum == "Very low", "sensitivity"]
check(abs(prim - drop) > 0.05, "leave-one-out changes a sparse stratum materially")
sw2 <- suppressWarnings(ssdta_meta_sensitivity(npoc_by_ultra_stratum, n_draws = 500, seed = 1,
  specs = list(meta_spec("RE", engine = "glmm"), meta_spec("Pooled", engine = "conjugate"))))
check(nlevels(sw2$table$analysis) == 2 &&
      sw2$table$sensitivity[sw2$table$analysis == "RE"][1] !=
      sw2$table$sensitivity[sw2$table$analysis == "Pooled"][1],
      "explicit specs override the engine per arm")
check(errors(ssdta_meta_sensitivity(npoc_by_ultra_stratum, engine = "glmm",
        specs = list(meta_spec("A"), meta_spec("A"))), "labels must be unique"),
      "sweep rejects duplicate spec labels")
check(errors(ssdta_meta_sensitivity(npoc_by_ultra_stratum, engine = "glmm", n_draws = 50,
        specs = meta_spec("X", exclude = unique(as.character(npoc_by_ultra_stratum$study)))),
        "excludes every study"), "sweep rejects a spec that drops all studies")
check(errors(meta_spec("X", engine = "banana"), "must be one of"),
      "meta_spec validates its engine")

# ---- ssdta_compare (two independent tests, same target) -------------------
set.seed(7)
worse <- npoc_posterior_draws
worse$sens_pooled <- stats::plogis(stats::qlogis(pmin(pmax(worse$sens_pooled,1e-4),1-1e-4)) - 0.5)
cA <- ssdta_standardize(npoc_posterior_draws, burden_profiles$screening_kendall, pi = 0.30, seed = 1)
cB <- ssdta_standardize(worse, burden_profiles$screening_kendall, pi = 0.30, seed = 2)
cmp <- ssdta_compare(cA, cB, labels = c("NPOC","Worse test"))
check(cmp$summary$median < 0, "compare: worse test has lower standardized sensitivity")
check(cmp$prob_B_lower > 0.9, "compare: P(worse < NPOC) high")
# A common target realization must cancel for identical paired sensitivity draws,
# even if the two standardizations were originally called with different seeds.
same_A <- ssdta_standardize(npoc_posterior_draws, burden_profiles$screening_kendall, 0.3, seed = 11)
same_B <- ssdta_standardize(npoc_posterior_draws, burden_profiles$screening_kendall, 0.3, seed = 99)
same_cmp <- ssdta_compare(same_A, same_B)
check(max(abs(same_cmp$draws)) < 1e-12, "compare enforces shared target weights")
check(errors(ssdta_compare(same_A, rt), "different burden profiles"),
      "compare rejects different target profiles")
check(is.data.frame(summary(r)) && is.data.frame(as.data.frame(cmp)),
      "summary/as.data.frame methods return tidy results")

# ---- DiD fallback still works (folded-in module) --------------------------
ex <- pooling_vs_npoc_example()
xc <- did_input(delta_A = ex$contrasts$delta_A, ci_A = ex$contrasts$ci_A,
                delta_B = ex$contrasts$delta_B, ci_B = ex$contrasts$ci_B, labels = ex$contrasts$labels)
check(near(xc$did, -0.0444, 1e-6), "DiD fallback point estimate = -4.44 pp")
check({w <- did_wald(xc); w$lower < w$estimate && w$estimate < w$upper}, "DiD Wald CI brackets estimate")

# ---- all figures build ----------------------------------------------------
p_ok <- tryCatch({
  ggplot_build(plot_casemix(burden_profiles))
  ggplot_build(plot_stratum_sensitivity(npoc_posterior_draws))
  ggplot_build(plot_forest_stratum(npoc_by_ultra_stratum, pooled = npoc_posterior_draws))
  ggplot_build(plot_forest_stratum(npoc_by_ultra_stratum, counts = FALSE))
  # the sweep plot must dodge: one visible point per (stratum x arm)
  sb <- ggplot_build(plot_meta_sensitivity(sw))
  stopifnot(length(unique(round(sb$data[[2]]$y, 4))) == 6 * (n_studies + 1))
  ggplot_build(plot(sw))
  ggplot_build(plot_pi_sensitivity(list(routine = g)))
  ggplot_build(plot_standardized(list(r, rt, sc)))
  ggplot_build(plot_spectrum_mechanism(
    npoc_posterior_draws,
    burden_profiles[c("research", "routine_pakistan", "screening_kendall")],
    list("observed", 0.20, 0.30)))
  ggplot_build(plot(r))
  ggplot_build(plot(cmp))
  ggplot_build(plot_parallel(did_input(sens = ex$components$sens, labels = ex$components$labels)))
  TRUE
}, error = function(e) { cat("   plot error:", conditionMessage(e), "\n"); FALSE })
check(p_ok, "all spectrum + DiD figures build")

cat(sprintf("\n%d passed, %d failed\n", ok, fail))
if (fail > 0) quit(status = 1)
