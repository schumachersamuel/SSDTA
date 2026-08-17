# ---------------------------------------------------------------------------
# Handling correlation between the two sensitivities being differenced ("rho").
#
# Distinguish two correlations (easy to conflate):
#   (1) PATIENT-LEVEL, within a single study: the same patients are tested by
#       both the index test and the anchor. Positive (shared bacillary burden),
#       and it lowers the within-study variance of that study's difference. It
#       is already captured inside a reported paired contrast CI.
#   (2) MARGINAL, across studies: correlation between the two RANDOM-EFFECTS
#       POOLED marginal sensitivities when you difference them. THIS is what the
#       package's `rho` controls, in the component form only.
#
# These are not the same number. For the TB worked example, backing rho out of
# the SRs' own reported contrast CIs gives rho ~ 0.02-0.03 -- effectively
# independent -- because large between-study heterogeneity in the difference
# swamps the patient-level pairing benefit. So the package default is rho = 0,
# and a "moderate" rho would wrongly narrow the interval below what the SR
# reported. Use implied_rho() to recover the data-driven value; prefer the
# contrast form (reported difference CI) when it is available.
# ---------------------------------------------------------------------------

#' Back-calculate the marginal correlation implied by a reported contrast CI
#'
#' Given a reported difference (test minus anchor) with its 95\% CI, and the two
#' marginal sensitivity CIs, returns the correlation \code{rho} between the two
#' pooled marginals that reconciles them:
#' rho = (SE_test^2 + SE_anchor^2 - SE_diff^2) / (2 * SE_test * SE_anchor).
#' Use it to set \code{rho} for the component form so the bootstrap/Bayesian
#' intervals are consistent with the SR's own reported contrast.
#'
#' @param diff_ci length-2 c(lower, upper) 95\% CI of the reported difference.
#' @param test_ci,anchor_ci length-2 95\% CIs of the two marginal sensitivities.
#' @return the implied correlation (scalar). Values near 0 mean the reported
#'   contrast CI is essentially the independence combination of the marginals.
#' @export
implied_rho <- function(diff_ci, test_ci, anchor_ci) {
  z <- stats::qnorm(0.975)
  se <- function(ci) (ci[2] - ci[1]) / (2 * z)
  sD <- se(diff_ci); sT <- se(test_ci); sA <- se(anchor_ci)
  (sT^2 + sA^2 - sD^2) / (2 * sT * sA)
}

#' Per-study paired sensitivity contrast and its SE, from 2x2 counts
#'
#' delta_i = Se(test1)_i - Se(test2)_i, with
#' Var(delta_i) = Var(Se1) + Var(Se2) - 2*rho*SD(Se1)*SD(Se2)
#' (marginal binomial SEs; rho supplies the within-pair correlation that
#' the two separate 2x2 tables cannot identify on their own).
#'
#' @param tp1,fn1 true-positive / false-negative counts for test 1, per study.
#' @param tp2,fn2 true-positive / false-negative counts for test 2, per study
#'   (same length, same study order).
#' @param rho assumed within-study correlation (scalar, or vector recycled to
#'   study length).
#' @param correct apply a Haldane-Anscombe continuity correction
#'   (p_c = (tp+0.5)/(n+1)) to the SE -- but not the point estimate -- whenever
#'   a study has perfect (100%) sensitivity on either test. Without this, a
#'   perfect-sensitivity study gets SE = 0, which assigns it infinite weight in
#'   any downstream meta-analysis (breaks tau^2 estimation). Default TRUE.
#' @return data.frame with delta, se, se1, se2, n1, n2 per study.
#' @export
paired_contrasts <- function(tp1, fn1, tp2, fn2, rho = 0, correct = TRUE) {
  lens <- lengths(list(tp1, fn1, tp2, fn2))
  if (!length(tp1) || length(unique(lens)) != 1L)
    stop("tp1, fn1, tp2, and fn2 must be non-empty vectors of equal length.")
  vals <- c(tp1, fn1, tp2, fn2)
  if (!is.numeric(vals) || any(!is.finite(vals)) || any(vals < 0))
    stop("tp/fn inputs must be finite, non-negative counts.")
  if (any(tp1 + fn1 <= 0) || any(tp2 + fn2 <= 0))
    stop("Each test-by-study denominator must be positive.")
  if (!is.numeric(rho) || any(!is.finite(rho)) || any(abs(rho) > 1))
    stop("`rho` must contain finite correlations in [-1, 1].")
  n1 <- tp1 + fn1; p1 <- tp1 / n1
  n2 <- tp2 + fn2; p2 <- tp2 / n2
  if (correct) {
    p1c <- ifelse(p1 %in% c(0, 1), (tp1 + 0.5) / (n1 + 1), p1)
    p2c <- ifelse(p2 %in% c(0, 1), (tp2 + 0.5) / (n2 + 1), p2)
  } else { p1c <- p1; p2c <- p2 }
  s1 <- sqrt(p1c * (1 - p1c) / n1)
  s2 <- sqrt(p2c * (1 - p2c) / n2)
  delta <- p1 - p2                                    # point estimate: uncorrected
  se <- sqrt(s1^2 + s2^2 - 2 * rho * s1 * s2)          # SE: continuity-corrected
  data.frame(delta = delta, se = se, se1 = s1, se2 = s2, n1 = n1, n2 = n2)
}

# NORTA (Normal-to-Anything): draw a correlation-rho bivariate standard normal,
# then map each margin through its own inverse-CDF. Exact marginals, an
# approximately-rho correlation on the target scale (approximation is due to
# discreteness for the binomial case; adequate for interval estimation here).
.rmvn2 <- function(B, rho) {
  if (!is.numeric(rho) || length(rho) != 1L || !is.finite(rho) || abs(rho) > 1)
    stop("`rho` must be one finite correlation in [-1, 1].")
  z1 <- stats::rnorm(B)
  z2 <- rho * z1 + sqrt(1 - rho^2) * stats::rnorm(B)
  cbind(z1, z2)
}

#' Correlated Binomial pair via a Gaussian copula
#' @keywords internal
.rcorr_binom_pair <- function(B, n1, p1, n2, p2, rho) {
  if (is.null(rho) || rho == 0) {
    return(list(x1 = stats::rbinom(B, n1, p1), x2 = stats::rbinom(B, n2, p2)))
  }
  z <- .rmvn2(B, rho)
  list(x1 = stats::qbinom(stats::pnorm(z[, 1]), n1, p1),
       x2 = stats::qbinom(stats::pnorm(z[, 2]), n2, p2))
}

#' Correlated Beta pair via a Gaussian copula (exact marginals)
#' @keywords internal
.rcorr_beta_pair <- function(B, a1, b1, a2, b2, rho) {
  if (is.null(rho) || rho == 0) {
    return(list(x1 = stats::rbeta(B, a1, b1), x2 = stats::rbeta(B, a2, b2)))
  }
  z <- .rmvn2(B, rho)
  list(x1 = stats::qbeta(stats::pnorm(z[, 1]), a1, b1),
       x2 = stats::qbeta(stats::pnorm(z[, 2]), a2, b2))
}
