# ---------------------------------------------------------------------------
# Comparing two tests.
#
# PRINCIPLED (when both tests have stratum-specific data): standardize each to a
# common target population, then difference. By linearity this is equivalent to
# a standardized within-stratum difference. It removes measured case-mix
# imbalance at the chosen stratum resolution under conditional transportability.
#
# FALLBACK (when one test lacks a stratum breakdown, e.g. pooled Ultra): the
# anchored difference-in-differences in this package's did_* functions, which
# partially deconfounds via a common anchor. Use ssdta_compare() when you can;
# fall back to anchored_did() when you must.
# ---------------------------------------------------------------------------

#' Compare two tests by differencing their standardized sensitivities
#'
#' Both tests must be standardized to the \strong{same} target population and pi
#' (use identical \code{profile}/\code{pi} in the two \code{ssdta_standardize}
#' calls). Returns the standardized sensitivity difference (test B minus test A)
#' with a credible interval.
#'
#' @param std_A,std_B \code{ssdta_std} objects (scalar-pi results) for the two
#'   tests, standardized to the same target.
#' @param labels length-2 character vector naming tests A and B.
#' @param level credible level.
#' @param draw_pairing whether sensitivity posterior rows correspond between
#'   tests ("paired", the default) or come from independent fits. Independent
#'   draws are deterministically permuted before differencing.
#' @param seed seed used only for independent draw pairing.
#' @return an \code{ssdta_compare} object (difference draws + summary).
#' @export
ssdta_compare <- function(std_A, std_B, labels = c("A", "B"), level = 0.95,
                          draw_pairing = c("paired", "independent"), seed = 1) {
  if (!inherits(std_A, "ssdta_std") || !inherits(std_B, "ssdta_std"))
    stop("`std_A` and `std_B` must both be ssdta_std objects.")
  .validate_level(level)
  draw_pairing <- match.arg(draw_pairing)
  if (length(labels) != 2L) stop("`labels` must contain exactly two test names.")
  if (!isTRUE(all.equal(std_A$profile$counts_detected, std_B$profile$counts_detected)))
    stop("std_A and std_B use different burden profiles; a common-reference comparison is required.")
  if (!identical(std_A$pi, std_B$pi))
    stop("std_A and std_B use different pi; a common-reference comparison is required.")

  n <- min(length(std_A$draws), length(std_B$draws))
  ia <- seq_len(n); ib <- seq_len(n)
  if (draw_pairing == "independent") ib <- .with_seed(seed, sample.int(n))

  # Reapply a single realization of the common target weights to both tests.
  # Older ssdta_std objects lack the stored components and fall back to their
  # precomputed draws with a warning.
  if (!is.null(std_A$weight_draws) && !is.null(std_A$sensitivity_draws) &&
      !is.null(std_B$sensitivity_draws)) {
    w <- std_A$weight_draws[ia, , drop = FALSE]
    A <- rowSums(std_A$sensitivity_draws[ia, colnames(w), drop = FALSE] * w)
    B <- rowSums(std_B$sensitivity_draws[ib, colnames(w), drop = FALSE] * w)
  } else {
    warning("Legacy ssdta_std object: common target-weight draws could not be enforced.")
    A <- std_A$draws[ia]; B <- std_B$draws[ib]
    w <- NULL
  }
  diff <- B - A
  s <- .summ_draws(diff, level)
  structure(list(
    draws = diff, labels = labels, level = level,
    setting = std_A$profile$setting, pi = std_A$pi,
    draw_pairing = draw_pairing, weight_draws = w,
    prob_B_lower = mean(diff < 0),
    summary = data.frame(comparison = sprintf("%s - %s", labels[2], labels[1]),
                         setting = std_A$profile$setting,
                         median = s["median"], lo = s["lo"], hi = s["hi"],
                         row.names = NULL)
  ), class = "ssdta_compare")
}

#' @export
summary.ssdta_compare <- function(object, ...) object$summary

#' @export
as.data.frame.ssdta_compare <- function(x, ...) x$summary

#' @export
plot.ssdta_compare <- function(x, ...) {
  d <- data.frame(difference = 100 * x$draws)
  ggplot2::ggplot(d, ggplot2::aes(.data$difference)) +
    ggplot2::geom_density(fill = "#D9ECF4", colour = "#0072B2", linewidth = 0.8) +
    ggplot2::geom_vline(xintercept = 0, linetype = 2, colour = "grey45") +
    ggplot2::labs(x = paste0(x$summary$comparison, " (percentage points)"),
                  y = "Posterior density", title = "Standardized sensitivity difference") +
    .ssdta_theme()
}

#' @export
print.ssdta_compare <- function(x, ...) {
  cat(sprintf("Standardized sensitivity difference: %s (target: %s)\n",
              x$summary$comparison, x$setting))
  s <- x$summary
  cat(sprintf("  %+.1f pp   %d%% CrI (%+.1f, %+.1f)   P(%s < %s) = %.3f\n",
              100*s$median, round(100*x$level), 100*s$lo, 100*s$hi,
              x$labels[2], x$labels[1], x$prob_B_lower))
  invisible(x)
}
