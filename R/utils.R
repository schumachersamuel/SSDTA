# ---------------------------------------------------------------------------
# Shared internal helpers.
# (`.pp` and `%||%` are defined in the folded-in DiD module files
#  did_plots.R / did_variance.R and are reused across the package.)
# ---------------------------------------------------------------------------

# Canonical Xpert Ultra semi-quant strata (TB default). The stratifier is
# generic; these are only the built-in default labels for the TB use case.
.ultra_strata_detected <- c("High", "Medium", "Low", "Very low", "Trace")
.ultra_strata_all      <- c(.ultra_strata_detected, "Not-detected")

# Colour-blind-safe palettes reused across spectrum plots. The stratum palette
# follows the ordered burden scale but remains distinguishable in print.
.stratum_palette <- c(
  "High"         = "#0072B2",
  "Medium"       = "#56B4E9",
  "Low"          = "#009E73",
  "Very low"     = "#E69F00",
  "Trace"        = "#D55E00",
  "Not-detected" = "#4d4d4d"
)
.setting_palette <- c(
  "Research"           = "#333333",
  "Routine diagnostic" = "#0072B2",
  "Screening"          = "#D55E00"
)

.ssdta_theme <- function(base_size = 11) {
  ggplot2::theme_minimal(base_size = base_size) +
    ggplot2::theme(
      panel.grid.minor = ggplot2::element_blank(),
      panel.grid.major.y = ggplot2::element_blank(),
      plot.title.position = "plot",
      plot.title = ggplot2::element_text(face = "bold"),
      plot.subtitle = ggplot2::element_text(colour = "grey35"),
      axis.title = ggplot2::element_text(colour = "grey20"),
      legend.position = "bottom"
    )
}

.validate_level <- function(level) {
  if (!is.numeric(level) || length(level) != 1L || !is.finite(level) ||
      level <= 0 || level >= 1)
    stop("`level` must be one finite number strictly between 0 and 1.")
  invisible(level)
}

.validate_counts <- function(counts, name = "counts", require_positive_total = TRUE) {
  if (!is.numeric(counts) || !length(counts) || any(!is.finite(counts)) || any(counts < 0))
    stop("`", name, "` must be a non-empty, finite, non-negative numeric vector.")
  if (require_positive_total && sum(counts) <= 0)
    stop("`", name, "` must have a positive total.")
  invisible(counts)
}

# Evaluate code under a reproducible seed without changing the caller's RNG
# stream. This matters in analysis pipelines that interleave SSDTA with other
# simulations.
.with_seed <- function(seed, code) {
  if (!is.numeric(seed) || length(seed) != 1L || !is.finite(seed))
    stop("`seed` must be one finite number.")
  had_seed <- exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)
  if (had_seed) old_seed <- get(".Random.seed", envir = .GlobalEnv, inherits = FALSE)
  on.exit({
    if (had_seed) assign(".Random.seed", old_seed, envir = .GlobalEnv)
    else if (exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE))
      rm(".Random.seed", envir = .GlobalEnv)
  }, add = TRUE)
  set.seed(as.integer(seed))
  force(code)
}

# Wilson score interval for a binomial proportion, vectorised over x and n.
# Used for the study-level points in forest plots: several study x stratum cells
# are at the boundary (tp = n, or tp = 0), where Wald intervals collapse to zero
# width and exact intervals are needlessly wide.
.wilson_ci <- function(x, n, level = 0.95) {
  .validate_level(level)
  if (any(n <= 0)) stop("Wilson interval needs n > 0.")
  z <- stats::qnorm(1 - (1 - level) / 2)
  p <- x / n
  denom <- 1 + z^2 / n
  centre <- (p + z^2 / (2 * n)) / denom
  half <- z * sqrt(p * (1 - p) / n + z^2 / (4 * n^2)) / denom
  cbind(lower = pmax(0, centre - half), upper = pmin(1, centre + half))
}

# equal-tailed interval summary of a numeric vector of draws
.summ_draws <- function(x, level = 0.95) {
  .validate_level(level)
  if (!is.numeric(x) || !length(x) || any(!is.finite(x)))
    stop("draws must be a non-empty, finite numeric vector.")
  a <- (1 - level) / 2
  c(median = stats::median(x),
    lo = as.numeric(stats::quantile(x, a)),
    hi = as.numeric(stats::quantile(x, 1 - a)),
    lo80 = as.numeric(stats::quantile(x, 0.10)),
    hi80 = as.numeric(stats::quantile(x, 0.90)))
}
