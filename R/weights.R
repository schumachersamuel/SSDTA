# ---------------------------------------------------------------------------
# Burden (case-mix) distribution profiles and the case-mix weight draws.
#
# A burden profile holds the observed counts of true (culture-positive) cases
# in each *detected* stratum of the target population. The weight of the
# Not-detected stratum (Ultra false-negatives) is NOT observable in routine
# data and enters separately as pi (see ssdta_standardize).
# ---------------------------------------------------------------------------

#' Construct a burden (case-mix) profile
#'
#' @param counts named numeric vector of case counts in the detected strata
#'   (e.g. \code{c(High=, Medium=, Low=, "Very low"=, Trace=)}). Order is
#'   preserved; names are the stratum labels (generic -- not TB-specific).
#' @param setting free-text setting label (e.g. "Screening").
#' @param count_notdetected optional observed count of culture-positive,
#'   test-negative cases (the Not-detected / false-negative stratum). Available
#'   only when the reference standard was applied to test-negatives too (e.g.
#'   research studies). When present, \code{ssdta_standardize(pi = "observed")}
#'   uses a full 6-category Dirichlet instead of a fixed pi.
#' @param country,source,note optional metadata.
#' @return object of class \code{burden_profile}.
#' @export
burden_profile <- function(counts, setting = NA, country = NA,
                           count_notdetected = NA, source = NA, note = NA) {
  if (is.null(names(counts)) || any(names(counts) == ""))
    stop("`counts` must be a fully named numeric vector.")
  if (anyDuplicated(names(counts))) stop("`counts` stratum names must be unique.")
  .validate_counts(counts)
  if (length(count_notdetected) != 1L ||
      (!is.na(count_notdetected) && (!is.numeric(count_notdetected) ||
       !is.finite(count_notdetected) || count_notdetected < 0)))
    stop("`count_notdetected` must be NA or one finite, non-negative number.")
  structure(list(setting = setting, country = country, source = source,
                 note = note, n = sum(counts), counts_detected = counts,
                 count_notdetected = count_notdetected),
            class = "burden_profile")
}

#' @export
print.burden_profile <- function(x, ...) {
  cat(sprintf("<burden profile> %s%s  (n = %s)\n",
              x$setting %||% "", if (!is.na(x$country)) paste0(" / ", x$country) else "",
              format(x$n, big.mark = ",")))
  cd <- x$counts_detected
  props <- cd / sum(cd)
  for (k in names(cd))
    cat(sprintf("  %-14s %8s  (%5.1f%%)\n", k, format(cd[[k]], big.mark=","), 100*props[[k]]))
  if (!is.na(x$source)) cat("  source:", x$source, "\n")
  if (!is.na(x$note))   cat("  note:  ", x$note, "\n")
  invisible(x)
}

#' Pool several burden profiles into one (crude count aggregation)
#'
#' Sums the per-stratum counts across profiles -- the simplest coherent
#' aggregation, which propagates the combined multinomial sampling uncertainty
#' through the Dirichlet in \code{ssdta_standardize}. Intended for combining
#' multiple datasets from the SAME setting type (e.g. two screening surveys).
#'
#' @param ... two or more \code{burden_profile} objects with identical stratum names.
#' @param setting,country,source,note metadata for the pooled profile.
#' @return a pooled \code{burden_profile}.
#' @note This is the "crude aggregate". A random-effects / Dirichlet-multinomial
#'   meta-analysis of case-mix distributions (carrying between-dataset
#'   heterogeneity) is the principled next step and is not yet implemented.
#'   Do NOT pool across different setting types (screening vs programmatic) --
#'   those are distinct reference targets meant to be reported separately.
#' @export
pool_profiles <- function(..., setting = "Pooled", country = NA,
                          source = NA, note = "crude count aggregation") {
  profs <- list(...)
  if (length(profs) < 2) stop("Provide at least two profiles.")
  if (!all(vapply(profs, inherits, logical(1), what = "burden_profile")))
    stop("Every object in `...` must be a burden_profile.")
  nm <- names(profs[[1]]$counts_detected)
  for (p in profs)
    if (!identical(names(p$counts_detected), nm))
      stop("All profiles must have identical stratum names in the same order.")
  settings <- unique(vapply(profs, function(p) as.character(p$setting), character(1)))
  if (length(settings) > 1)
    warning("Pooling profiles from different settings (", paste(settings, collapse=", "),
            "). Standardization targets should usually be a single setting type.")
  total <- Reduce(`+`, lapply(profs, function(p) p$counts_detected))
  nd <- vapply(profs, function(p) p$count_notdetected, numeric(1))
  pooled_nd <- if (all(!is.na(nd))) sum(nd) else NA_real_
  if (any(!is.na(nd)) && any(is.na(nd)))
    warning("Some profiles have observed Not-detected counts and others do not; ",
            "the pooled profile leaves count_notdetected unknown.")
  burden_profile(total, setting = setting, country = country,
                 count_notdetected = pooled_nd,
                 source = source %||% paste("pool of", length(profs), "profiles"), note = note)
}

# Dirichlet(counts + 1) draws given a vector of alpha counts. [n_draws x K]
.dirichlet_draws <- function(counts, n_draws) {
  .validate_counts(counts)
  if (!is.numeric(n_draws) || length(n_draws) != 1L || !is.finite(n_draws) ||
      n_draws < 1 || n_draws != as.integer(n_draws))
    stop("`n_draws` must be a positive integer.")
  K <- length(counts); alpha <- counts + 1
  g <- matrix(stats::rgamma(n_draws * K, shape = rep(alpha, each = n_draws), rate = 1),
              nrow = n_draws, ncol = K)
  g / rowSums(g)
}

# 6-stratum weight draws. Two modes, matching spectrum_standardization_pipeline:
#   pi scalar   -> 5 detected via Dirichlet(counts+1), scaled by (1-pi); ND = pi.
#   pi="observed" -> full 6-category Dirichlet on (detected, ND) observed counts.
.weight_draws <- function(counts_detected, pi, n_draws, count_notdetected = NA,
                          not_detected_label = "Not-detected") {
  .validate_counts(counts_detected, "counts_detected")
  if (identical(pi, "observed")) {
    if (length(count_notdetected) != 1L || is.na(count_notdetected))
      stop('pi = "observed" needs profile$count_notdetected (observed FN count).')
    .validate_counts(count_notdetected, "count_notdetected", require_positive_total = FALSE)
    w6 <- .dirichlet_draws(c(counts_detected, count_notdetected), n_draws)
    colnames(w6) <- c(names(counts_detected), not_detected_label)
    return(w6)
  }
  if (!is.numeric(pi) || length(pi) != 1L || !is.finite(pi) || pi < 0 || pi > 1)
    stop("`pi` must be one finite number in [0, 1], or \"observed\".")
  w5 <- .dirichlet_draws(counts_detected, n_draws)
  w6 <- cbind(w5 * (1 - pi), pi)
  colnames(w6) <- c(names(counts_detected), not_detected_label)
  w6
}
