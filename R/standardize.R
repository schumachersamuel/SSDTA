# ---------------------------------------------------------------------------
# The core method: reweight stratum-specific sensitivity posterior draws to a
# target population's burden distribution.
#
#   S_std = sum_k  w_k * S_k
#
# Uncertainty is propagated on both sides: S_k as posterior draws (from the
# meta-analysis), w_k as Dirichlet case-mix draws. One weight draw is paired
# with each sensitivity draw.
# ---------------------------------------------------------------------------

# Reshape a long draws data.frame (category, draw, sens_pooled) into a matrix
# [n_draws x n_strata] with columns in `strata` order. `strata` includes the
# Not-detected stratum last.
.sens_matrix <- function(sens_draws, strata) {
  if (is.matrix(sens_draws)) {
    if (is.null(colnames(sens_draws)) || anyDuplicated(colnames(sens_draws)))
      stop("sens_draws matrix must have unique stratum column names.")
    if (!all(strata %in% colnames(sens_draws)))
      stop("sens_draws matrix is missing strata: ",
           paste(setdiff(strata, colnames(sens_draws)), collapse = ", "))
    m <- sens_draws[, strata, drop = FALSE]
    if (!nrow(m) || !is.numeric(m) || any(!is.finite(m)) || any(m < 0 | m > 1))
      stop("sens_draws matrix values must be finite probabilities in [0, 1].")
    return(m)
  }
  need <- c("category", "draw", "sens_pooled")
  if (!all(need %in% names(sens_draws)))
    stop("sens_draws must be a matrix or a data.frame with columns ",
         paste(need, collapse = ", "), ".")
  sd <- sens_draws
  sd$category <- as.character(sd$category)
  if (!all(strata %in% unique(sd$category)))
    stop("sens_draws is missing strata: ",
         paste(setdiff(strata, unique(sd$category)), collapse = ", "))
  if (anyNA(sd$draw) || anyDuplicated(sd[c("category", "draw")]))
    stop("sens_draws must contain exactly one row per (category, draw) pair.")
  draw_ids <- sort(unique(sd$draw))
  if (!length(draw_ids)) stop("sens_draws contains no draws.")
  m <- matrix(NA_real_, nrow = length(draw_ids), ncol = length(strata),
              dimnames = list(NULL, strata))
  for (k in strata) {
    rows <- sd[sd$category == k, c("draw", "sens_pooled")]
    if (!setequal(rows$draw, draw_ids))
      stop("sens_draws has inconsistent draw IDs across strata.")
    m[match(rows$draw, draw_ids), k] <- rows$sens_pooled
  }
  if (anyNA(m)) stop("sens_draws has missing (draw x stratum) cells.")
  if (!is.numeric(m) || any(!is.finite(m)) || any(m < 0 | m > 1))
    stop("sens_draws values must be finite probabilities in [0, 1].")
  m
}

#' Spectrum-standardize a test's sensitivity to a target population
#'
#' Reweights stratum-specific sensitivity posterior draws by a target
#' population's burden distribution, returning the standardized sensitivity
#' with a credible interval. If \code{pi} is a vector, runs the standardization
#' across the grid (a pi sensitivity analysis) and returns a tidy table.
#'
#' @param sens_draws stratum-specific sensitivity draws: either the long
#'   data.frame form (columns \code{category}, \code{draw}, \code{sens_pooled};
#'   e.g. \code{npoc_posterior_draws}) or a matrix \code{[n_draws x n_strata]}
#'   with stratum names as columns. Must include the Not-detected stratum.
#' @param profile a \code{burden_profile} (or a name in \code{burden_profiles}).
#' @param pi weight of the Not-detected stratum (Ultra false-negative rate).
#'   Scalar for a single estimate; a numeric vector for a sensitivity grid; or
#'   the string \code{"observed"} to use the profile's observed Not-detected
#'   count in a full 6-category Dirichlet (only when
#'   \code{profile$count_notdetected} is available, e.g. research data).
#' @param level credible level.
#' @param not_detected_label name of the Not-detected stratum in \code{sens_draws}.
#' @param seed RNG seed for the Dirichlet weight draws.
#' @return For scalar/observed \code{pi}: an \code{ssdta_std} object (draws +
#'   summary). For a vector \code{pi}: a data.frame with one row per pi.
#' @export
ssdta_standardize <- function(sens_draws, profile, pi,
                              level = 0.95, not_detected_label = "Not-detected",
                              seed = 42) {
  if (is.character(profile) && length(profile) == 1 && !identical(profile, "observed") &&
      is.null(dim(profile)) && !inherits(profile, "burden_profile")) {
    builtins <- get0("burden_profiles", inherits = TRUE)
    if (is.null(builtins) || is.null(builtins[[profile]]))
      stop("Unknown burden profile name: ", profile)
    profile <- builtins[[profile]]
  }
  if (!inherits(profile, "burden_profile")) stop("`profile` must be a burden_profile.")
  .validate_level(level)
  if (!(identical(pi, "observed") ||
        (is.numeric(pi) && length(pi) >= 1L && all(is.finite(pi)) && all(pi >= 0 & pi <= 1))))
    stop("`pi` must be \"observed\" or a finite numeric scalar/grid in [0, 1].")
  det <- profile$counts_detected
  strata_all <- c(names(det), not_detected_label)
  sens_mat <- .sens_matrix(sens_draws, strata_all)
  nd <- nrow(sens_mat)

  # Draw detected-case proportions once and reuse them across a pi grid. These
  # common random numbers remove avoidable Monte-Carlo wiggle and make each
  # draw's sensitivity curve exactly coherent across pi.
  base_w <- .with_seed(seed, .dirichlet_draws(det, nd))
  numeric_weights <- function(p) {
    w <- cbind(base_w * (1 - p), p)
    colnames(w) <- strata_all
    w
  }

  if (identical(pi, "observed") || length(pi) == 1L) {
    w6 <- if (identical(pi, "observed"))
      .with_seed(seed, .weight_draws(det, pi, nd, profile$count_notdetected,
                                    not_detected_label)) else numeric_weights(pi)
    S <- rowSums(sens_mat * w6)
    s <- .summ_draws(S, level)
    structure(list(
      draws = S, sensitivity_draws = sens_mat, weight_draws = w6,
      pi = pi, level = level, profile = profile,
      summary = data.frame(setting = profile$setting, pi = pi,
                           median = s["median"], lo = s["lo"], hi = s["hi"],
                           lo80 = s["lo80"], hi80 = s["hi80"], row.names = NULL)
    ), class = "ssdta_std")
  } else {
    rows <- lapply(seq_along(pi), function(j) {
      S <- rowSums(sens_mat * numeric_weights(pi[j]))
      s <- .summ_draws(S, level)
      data.frame(setting = profile$setting, pi = pi[j],
                 median = s["median"], lo = s["lo"], hi = s["hi"],
                 lo80 = s["lo80"], hi80 = s["hi80"], row.names = NULL)
    })
    out <- do.call(rbind, rows)
    attr(out, "profile") <- profile
    out
  }
}

#' @export
summary.ssdta_std <- function(object, ...) object$summary

#' @export
as.data.frame.ssdta_std <- function(x, ...) x$summary

#' @export
plot.ssdta_std <- function(x, ...) plot_standardized(list(x), ...)

#' @export
print.ssdta_std <- function(x, ...) {
  pilab <- if (is.numeric(x$pi)) sprintf("pi = %.2f", x$pi) else sprintf("pi = %s", x$pi)
  cat(sprintf("Spectrum-standardized sensitivity -- %s (%s)\n", x$profile$setting, pilab))
  s <- x$summary
  cat(sprintf("  S_std = %.1f%%   %d%% CrI (%.1f%%, %.1f%%)\n",
              100*s$median, round(100*x$level), 100*s$lo, 100*s$hi))
  invisible(x)
}
