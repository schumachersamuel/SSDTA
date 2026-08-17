# ---------------------------------------------------------------------------
# Input harmonisation for the anchored difference-in-differences estimator.
#
# The estimand is:
#   Delta_A  = Se(A) - Se(C)   estimated in the study set comparing A to anchor C
#   Delta_B  = Se(B) - Se(C)   estimated in the study set comparing B to anchor C
#   Delta_DiD = Delta_B - Delta_A
#
# Anchor C is the common comparator (in the TB application: individual Xpert
# Ultra, which is also the spectrum stratifier via its semiquantitative grade).
# ---------------------------------------------------------------------------

#' Build a canonical anchored-DiD input object
#'
#' Accepts either (i) the two paired contrasts directly (Delta_A, Delta_B with
#' SEs or CIs -- the form WHO GDG evidence tables report), or (ii) the four
#' component sensitivities Se(A), Se(C|a), Se(B), Se(C|b) with per-sensitivity
#' uncertainty supplied as a standard error, a 95\% CI, or a denominator n
#' (number of culture/reference positives). Component input additionally enables
#' the bootstrap and Bayesian methods.
#'
#' All sensitivities and differences are on the 0-1 probability scale
#' internally; pass proportions (0.848), not percentages (84.8).
#'
#' @param delta_A,delta_B point estimates of the two paired contrasts (B/A minus C).
#' @param se_A,se_B standard errors of delta_A / delta_B.
#' @param ci_A,ci_B length-2 numeric c(lower, upper) 95\% CIs; alternative to se_*.
#' @param se_A_i,se_B_i optional numeric vectors of per-study contrast SEs (for did_meta).
#' @param delta_A_i,delta_B_i optional numeric vectors of per-study contrasts (for did_meta).
#' @param sens named list/vector for the component form, with elements
#'   \code{A}, \code{Ca}, \code{B}, \code{Cb} (point sensitivities), plus one of
#'   \code{n} (named vector of denominators) or \code{se} (named vector of SEs)
#'   or \code{tp}/\code{fn} (named vectors of true-positive / false-negative counts).
#' @param rho correlation between the two POOLED marginal sensitivity estimates
#'   being differenced -- Se(test) and Se(anchor) -- used by the component form
#'   only, to compute contrast SEs
#'   (Var(delta) = SE_test^2 + SE_anchor^2 - 2*rho*SE_test*SE_anchor) and, via
#'   \code{did_boot}/\code{did_bayes}, to induce the matching correlation in the
#'   resampled/posterior draws (Gaussian-copula NORTA; see \code{paired.R}).
#'
#'   \strong{This is NOT the patient-level pairing.} Two distinct correlations
#'   are easy to conflate: (1) within a single study, the same patients are
#'   tested by both the index test and the anchor -- positive, and it lowers the
#'   within-study variance of that study's difference; (2) across studies, the
#'   two RANDOM-EFFECTS POOLED marginal sensitivities co-vary -- this is what
#'   \code{rho} controls. When you feed the reported paired contrast CI (the
#'   contrast form) these are both already baked in and \code{rho} is not used.
#'   \code{rho} matters only when you rebuild the contrast from the two separate
#'   marginal sensitivities.
#'
#'   Default \strong{0} (independence of the pooled marginals). This looks
#'   surprising -- surely the tests are correlated? -- but it is what the data
#'   demand: for the TB worked example, backing rho out of the reports' own
#'   contrast CIs gives rho ~ 0.02-0.03 on both sides (see \code{implied_rho}).
#'   The reason is that between-study heterogeneity in the difference is large
#'   (e.g. one study's NPOC-Ultra gap is -33pp vs -4 to -10pp elsewhere) and
#'   swamps the patient-level pairing benefit, so the pooled marginals behave as
#'   near-independent. Using a "moderate" rho = 0.5 here would produce a DiD
#'   interval NARROWER than the SR's own reported contrast CIs -- indefensible.
#'   If you have a single low-heterogeneity study, or the true joint
#'   (test x anchor) concordance table, a higher rho may be right; back it out
#'   with \code{implied_rho()} rather than guessing. Best of all: use the
#'   contrast form when the reported difference CI is available.
#' @param labels character vector of length 3 naming tests A, B, C.
#' @return An object of class \code{did_input}.
#' @export
did_input <- function(delta_A = NULL, delta_B = NULL,
                      se_A = NULL, se_B = NULL,
                      ci_A = NULL, ci_B = NULL,
                      delta_A_i = NULL, delta_B_i = NULL,
                      se_A_i = NULL, se_B_i = NULL,
                      sens = NULL, rho = 0,
                      labels = c("A", "B", "C")) {

  z <- stats::qnorm(0.975)
  comp <- NULL
  if (!is.numeric(rho) || length(rho) != 1L || !is.finite(rho) || abs(rho) > 1)
    stop("`rho` must be one finite correlation in [-1, 1].")
  if (length(labels) != 3L || anyNA(labels))
    stop("`labels` must contain exactly three non-missing test names.")

  # ---- component-sensitivity form --------------------------------------
  if (!is.null(sens)) {
    need <- c("A", "Ca", "B", "Cb")
    if (!all(need %in% names(sens)))
      stop("`sens` must contain elements A, Ca, B, Cb.")
    se_pt <- vapply(need, function(k) as.numeric(sens[[k]]), numeric(1))
    if (any(!is.finite(se_pt)) || any(se_pt < 0 | se_pt > 1))
      stop("Component sensitivities must be proportions in [0, 1].")

    # per-component counts / SE
    tp <- fn <- n <- se_c <- stats::setNames(rep(NA_real_, 4), need)
    if (!is.null(sens$tp) && !is.null(sens$fn)) {
      tp[need] <- as.numeric(sens$tp[need]); fn[need] <- as.numeric(sens$fn[need])
      n[need]  <- tp[need] + fn[need]
      se_pt    <- tp[need] / n[need]                       # override with exact
      se_c[need] <- sqrt(se_pt * (1 - se_pt) / n[need])
    } else if (!is.null(sens$n)) {
      n[need]  <- as.numeric(sens$n[need])
      tp[need] <- round(se_pt * n[need]); fn[need] <- n[need] - tp[need]
      se_c[need] <- sqrt(se_pt * (1 - se_pt) / n[need])
    } else if (!is.null(sens$se)) {
      se_c[need] <- as.numeric(sens$se[need])
      n[need]  <- se_pt * (1 - se_pt) / se_c[need]^2       # effective n
      tp[need] <- round(se_pt * n[need]); fn[need] <- round(n[need] - tp[need])
    } else {
      stop("Component form needs one of sens$n, sens$se, or sens$tp+sens$fn.")
    }

    comp <- data.frame(
      part = need, sens = se_pt, se = se_c[need],
      n = n[need], tp = tp[need], fn = fn[need],
      row.names = NULL, stringsAsFactors = FALSE
    )

    delta_A <- se_pt["A"]  - se_pt["Ca"]
    delta_B <- se_pt["B"]  - se_pt["Cb"]
    # paired within a study set: rho defaults to 0 (conservative independence).
    se_A <- sqrt(se_c["A"]^2  + se_c["Ca"]^2 - 2 * rho * se_c["A"]  * se_c["Ca"])
    se_B <- sqrt(se_c["B"]^2  + se_c["Cb"]^2 - 2 * rho * se_c["B"]  * se_c["Cb"])
  }

  # ---- contrast form ---------------------------------------------------
  if (is.null(delta_A) || is.null(delta_B))
    stop("Provide the two contrasts (delta_A, delta_B) or the `sens` component form.")

  if (is.null(se_A)) {
    if (!is.null(ci_A)) se_A <- (ci_A[2] - ci_A[1]) / (2 * z)
    else if (!is.null(se_A_i) && !is.null(delta_A_i)) se_A <- NA_real_  # filled by meta
    else stop("Need se_A, ci_A, or per-study delta_A_i + se_A_i.")
  }
  if (is.null(se_B)) {
    if (!is.null(ci_B)) se_B <- (ci_B[2] - ci_B[1]) / (2 * z)
    else if (!is.null(se_B_i) && !is.null(delta_B_i)) se_B <- NA_real_
    else stop("Need se_B, ci_B, or per-study delta_B_i + se_B_i.")
  }

  structure(list(
    delta_A = as.numeric(delta_A), delta_B = as.numeric(delta_B),
    se_A = as.numeric(se_A), se_B = as.numeric(se_B),
    did = as.numeric(delta_B - delta_A),
    delta_A_i = delta_A_i, delta_B_i = delta_B_i,
    se_A_i = se_A_i, se_B_i = se_B_i,
    components = comp, labels = labels, rho = rho
  ), class = "did_input")
}

#' @export
print.did_input <- function(x, ...) {
  lab <- x$labels
  cat("<anchored-DiD input>\n")
  cat(sprintf("  Tests: A = %s, B = %s ; anchor C = %s\n", lab[1], lab[2], lab[3]))
  cat(sprintf("  Delta_A = Se(%s) - Se(%s) = %+.1f pp%s\n",
              lab[1], lab[3], 100 * x$delta_A,
              if (!is.na(x$se_A)) sprintf("  (SE %.1f)", 100 * x$se_A) else ""))
  cat(sprintf("  Delta_B = Se(%s) - Se(%s) = %+.1f pp%s\n",
              lab[2], lab[3], 100 * x$delta_B,
              if (!is.na(x$se_B)) sprintf("  (SE %.1f)", 100 * x$se_B) else ""))
  cat(sprintf("  Delta_DiD (%s vs %s) = %+.1f pp\n", lab[2], lab[1], 100 * x$did))
  if (!is.null(x$components)) {
    cat(sprintf("  Component sensitivities available (bootstrap + Bayesian enabled); rho = %.2f\n",
                x$rho))
  }
  invisible(x)
}

#' Reconstruct approximate true-positive / false-negative counts
#'
#' From a sensitivity plus one of an SE, a 95\% CI, or a denominator n. The SE
#' route uses the binomial effective sample size n_eff = Se(1-Se)/SE^2.
#'
#' @param sens sensitivity (proportion).
#' @param se optional standard error of the sensitivity.
#' @param ci optional c(lower, upper) 95\% CI.
#' @param n optional denominator (number of reference positives).
#' @return list(tp, fn, n, se).
#' @export
reconstruct_counts <- function(sens, se = NULL, ci = NULL, n = NULL) {
  z <- stats::qnorm(0.975)
  if (!is.numeric(sens) || length(sens) != 1L || !is.finite(sens) || sens < 0 || sens > 1)
    stop("`sens` must be one finite proportion in [0, 1].")
  if (is.null(se) && !is.null(ci)) se <- (ci[2] - ci[1]) / (2 * z)
  if (is.null(n)) {
    if (is.null(se)) stop("Provide se, ci, or n.")
    n <- sens * (1 - sens) / se^2
  }
  if (!is.numeric(n) || length(n) != 1L || !is.finite(n) || n <= 0)
    stop("`n` must be one positive finite number.")
  tp <- round(sens * n); fn <- round(n - tp)
  if (is.null(se)) se <- sqrt(sens * (1 - sens) / n)
  list(tp = tp, fn = fn, n = n, se = se)
}
