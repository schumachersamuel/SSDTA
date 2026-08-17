# ---------------------------------------------------------------------------
# Four approaches to interval estimation for Delta_DiD, matching the tabs in the
# methods reference. Each returns a tidy one-row-per-method data frame via the
# internal .fit_row() helper and is wrapped by anchored_did().
# ---------------------------------------------------------------------------

.fit_row <- function(method, est, se, lower, upper, extra = list()) {
  data.frame(method = method, estimate = est, se = se,
             lower = lower, upper = upper, stringsAsFactors = FALSE)
}

#' Method 1 -- Wald / delta (Bucher sum of variances)
#'
#' Var(Delta_DiD) = Var(Delta_A) + Var(Delta_B); the two contrasts are
#' independent because they come from separate study sets. Fast and transparent;
#' the normal approximation degrades when a component sensitivity is near 0 or 1
#' with small counts (rule of thumb: worry when Se >= ~0.90 or <= ~0.10 and the
#' per-arm reference-positive count is under ~30-50).
#' @param x a did_input object.
#' @param level confidence level.
#' @export
did_wald <- function(x, level = 0.95) {
  stopifnot(inherits(x, "did_input"))
  .validate_level(level)
  if (is.na(x$se_A) || is.na(x$se_B))
    stop("did_wald needs se_A and se_B (or a CI); use did_meta for per-study data.")
  z <- stats::qnorm(1 - (1 - level) / 2)
  se <- sqrt(x$se_A^2 + x$se_B^2)
  .fit_row("Wald / delta", x$did, se, x$did - z * se, x$did + z * se)
}

# --- small base-R random-effects pooler (DL and Paule-Mandel), with optional
#     Hartung-Knapp small-sample adjustment. Avoids a hard metafor dependency.
.meta_pool <- function(yi, sei, tau2_method = c("PM", "DL"), hakn = TRUE,
                       level = 0.95) {
  tau2_method <- match.arg(tau2_method)
  vi <- sei^2; k <- length(yi)
  wi_fe <- 1 / vi
  y_fe <- sum(wi_fe * yi) / sum(wi_fe)
  Q <- sum(wi_fe * (yi - y_fe)^2)

  if (tau2_method == "DL") {
    C <- sum(wi_fe) - sum(wi_fe^2) / sum(wi_fe)
    tau2 <- max(0, (Q - (k - 1)) / C)
  } else { # Paule-Mandel: solve sum w_i(tau2) (yi - ybar)^2 = k - 1
    f <- function(t2) {
      w <- 1 / (vi + t2)
      yb <- sum(w * yi) / sum(w)
      sum(w * (yi - yb)^2) - (k - 1)
    }
    tau2 <- if (f(0) <= 0) 0 else {
      hi <- 10 * max(vi) + 1
      while (f(hi) > 0 && hi < 1e6) hi <- hi * 2
      stats::uniroot(f, c(0, hi))$root
    }
  }

  w <- 1 / (vi + tau2)
  est <- sum(w * yi) / sum(w)
  var_est <- 1 / sum(w)

  if (hakn && k > 1) {
    q <- sum(w * (yi - est)^2) / (k - 1)
    se <- sqrt(q * var_est)
    crit <- stats::qt(1 - (1 - level) / 2, df = k - 1)
  } else {
    se <- sqrt(var_est)
    crit <- stats::qnorm(1 - (1 - level) / 2)
  }
  list(est = est, se = se, tau2 = tau2, k = k,
       lower = est - crit * se, upper = est + crit * se, hakn = hakn)
}

#' Method 2 -- Meta-analytic propagation
#'
#' Random-effects meta-analysis of the per-study within-pair contrasts on each
#' side, then subtract the two pooled estimates and add their variances. Carries
#' between-study heterogeneity (tau^2) into each side. tau^2 is poorly estimated
#' with few studies: with < 5 studies per side treat it cautiously (or use
#' \code{hakn = TRUE}, the default, which widens the interval via a t-reference);
#' >= 10 studies gives a stable tau^2.
#' @param x a did_input object carrying delta_*_i and se_*_i per-study vectors.
#' @param tau2_method "PM" (Paule-Mandel, default) or "DL".
#' @param hakn logical; Hartung-Knapp-Sidik-Jonkman small-sample adjustment.
#' @param level confidence level.
#' @export
did_meta <- function(x, tau2_method = "PM", hakn = TRUE, level = 0.95) {
  stopifnot(inherits(x, "did_input"))
  .validate_level(level)
  if (is.null(x$delta_A_i) || is.null(x$delta_B_i))
    stop("did_meta needs per-study contrasts: delta_A_i/se_A_i and delta_B_i/se_B_i.")
  pa <- .meta_pool(x$delta_A_i, x$se_A_i, tau2_method, hakn, level)
  pb <- .meta_pool(x$delta_B_i, x$se_B_i, tau2_method, hakn, level)
  est <- pb$est - pa$est
  se  <- sqrt(pa$se^2 + pb$se^2)
  crit <- stats::qnorm(1 - (1 - level) / 2)
  out <- .fit_row("Meta-analytic", est, se, est - crit * se, est + crit * se)
  attr(out, "sides") <- list(A = pa, B = pb)
  out
}

#' Method 3 -- Parametric bootstrap
#'
#' Resamples true-positive counts for each of the four component sensitivities
#' from their Binomial(n, Se) sampling model and recomputes Delta_DiD. Respects
#' the [0,1] bound and small-sample skew. Requires the component form (counts or
#' a reconstructable n).
#'
#' The (A, Ca) pair and the (B, Cb) pair are each resampled with the within-pair
#' correlation \code{rho} (see \code{did_input}), via a Gaussian-copula (NORTA)
#' construction that preserves each margin's exact Binomial distribution. A and
#' B are independent of each other (separate study sets/populations), matching
#' Bucher's independent-arm assumption. Set \code{rho = 0} to recover plain
#' independent resampling.
#' @param x a did_input object with components.
#' @param B number of bootstrap replicates.
#' @param level confidence level.
#' @param rho within-pair correlation; defaults to \code{x$rho} (set in
#'   \code{did_input()}, default 0).
#' @param seed optional RNG seed.
#' @export
did_boot <- function(x, B = 5000, level = 0.95, rho = NULL, seed = NULL) {
  stopifnot(inherits(x, "did_input"))
  .validate_level(level)
  if (!is.numeric(B) || length(B) != 1L || !is.finite(B) || B < 1 || B != as.integer(B))
    stop("`B` must be a positive integer.")
  if (is.null(x$components)) stop("did_boot needs the component sensitivity form.")
  if (!is.null(seed)) set.seed(seed)
  if (is.null(rho)) rho <- x$rho %||% 0
  cp <- x$components; rownames(cp) <- cp$part

  pairAC <- .rcorr_binom_pair(B, round(cp["A","n"]),  cp["A","sens"],
                              round(cp["Ca","n"]), cp["Ca","sens"], rho)
  pairBC <- .rcorr_binom_pair(B, round(cp["B","n"]),  cp["B","sens"],
                              round(cp["Cb","n"]), cp["Cb","sens"], rho)
  sA  <- pairAC$x1 / round(cp["A","n"]);  sCa <- pairAC$x2 / round(cp["Ca","n"])
  sB  <- pairBC$x1 / round(cp["B","n"]);  sCb <- pairBC$x2 / round(cp["Cb","n"])
  did <- (sB - sCb) - (sA - sCa)

  a <- (1 - level) / 2
  qs <- stats::quantile(did, c(a, 1 - a), names = FALSE)
  out <- .fit_row(sprintf("Bootstrap (rho=%.2f)", rho), x$did, stats::sd(did), qs[1], qs[2])
  attr(out, "draws") <- did
  out
}

`%||%` <- function(a, b) if (is.null(a)) b else a

#' Construct a Beta prior for a sensitivity
#'
#' @param name label for the prior.
#' @param a,b Beta shape parameters. Provide these OR center+concentration.
#' @param center prior mean (with \code{concentration}).
#' @param concentration prior "sample size" a+b (with \code{center}).
#' @return list(name, a, b).
#' @export
prior_beta <- function(name, a = NULL, b = NULL, center = NULL, concentration = NULL) {
  if (!is.null(center) && !is.null(concentration)) {
    a <- center * concentration; b <- (1 - center) * concentration
  }
  if (is.null(a) || is.null(b)) stop("Give a,b or center+concentration.")
  if (!is.numeric(a) || !is.numeric(b) || length(a) != 1L || length(b) != 1L ||
      !is.finite(a) || !is.finite(b) || a <= 0 || b <= 0)
    stop("Beta prior shape parameters must be finite and strictly positive.")
  list(name = name, a = a, b = b)
}

# built-in reasonable priors for a diagnostic sensitivity
.default_priors <- function() list(
  prior_beta("Jeffreys Beta(0.5,0.5)", a = 0.5, b = 0.5),
  prior_beta("Uniform Beta(1,1)",     a = 1,   b = 1),
  prior_beta("Weakly-informative (mean 0.85, k=10)", center = 0.85, concentration = 10),
  prior_beta("Weakly-informative (mean 0.85, k=40)", center = 0.85, concentration = 40)
)

#' Method 4 -- Conjugate Beta-Binomial Bayesian
#'
#' Independent Beta priors on each component sensitivity, updated by the
#' reconstructed TP/FN counts to Beta posteriors, then Monte-Carlo the posterior
#' of Delta_DiD. Works on aggregate data (counts reconstructed from Se + SE/CI/n
#' via the binomial effective sample size). Returns a posterior median and
#' equal-tailed credible interval.
#'
#' As in \code{did_boot}, the (A, Ca) and (B, Cb) posterior draws are linked by
#' the within-pair correlation \code{rho} via a Gaussian copula, so each
#' component's posterior is an exact Beta marginal while the pair carries the
#' assumed correlation. Set \code{rho = 0} for independent posteriors.
#' @param x a did_input object with components.
#' @param prior a \code{prior_beta()} object applied to all four sensitivities
#'   (default Jeffreys), or a named list with elements A, Ca, B, Cb.
#' @param draws number of posterior Monte-Carlo draws.
#' @param level credible level.
#' @param rho marginal-pair correlation; defaults to \code{x$rho} (default 0).
#' @param seed optional RNG seed.
#' @export
did_bayes <- function(x, prior = prior_beta("Jeffreys Beta(0.5,0.5)", a = 0.5, b = 0.5),
                      draws = 20000, level = 0.95, rho = NULL, seed = NULL) {
  stopifnot(inherits(x, "did_input"))
  .validate_level(level)
  if (!is.numeric(draws) || length(draws) != 1L || !is.finite(draws) ||
      draws < 1 || draws != as.integer(draws))
    stop("`draws` must be a positive integer.")
  if (is.null(x$components)) stop("did_bayes needs the component sensitivity form.")
  if (!is.null(seed)) set.seed(seed)
  if (is.null(rho)) rho <- x$rho %||% 0
  cp <- x$components; rownames(cp) <- cp$part
  get_prior <- function(part) if (!is.null(prior$name)) prior else prior[[part]]
  shp <- function(part) { pr <- get_prior(part); c(a = pr$a + cp[part,"tp"], b = pr$b + cp[part,"fn"]) }
  shA <- shp("A"); shCa <- shp("Ca"); shB <- shp("B"); shCb <- shp("Cb")

  pairAC <- .rcorr_beta_pair(draws, shA["a"], shA["b"], shCa["a"], shCa["b"], rho)
  pairBC <- .rcorr_beta_pair(draws, shB["a"], shB["b"], shCb["a"], shCb["b"], rho)
  sA <- pairAC$x1; sCa <- pairAC$x2; sB <- pairBC$x1; sCb <- pairBC$x2
  did <- (sB - sCb) - (sA - sCa)

  a <- (1 - level) / 2
  qs <- stats::quantile(did, c(a, 0.5, 1 - a), names = FALSE)
  pr_name <- if (!is.null(prior$name)) prior$name else "custom per-component"
  out <- .fit_row(sprintf("Bayesian [%s, rho=%.2f]", pr_name, rho),
                  qs[2], stats::sd(did), qs[1], qs[3])
  attr(out, "draws") <- did
  attr(out, "prob_negative") <- mean(did < 0)
  out
}

#' Bayesian prior-sensitivity analysis
#'
#' Refits the Bayesian estimator across a set of priors and returns a combined
#' table plus the posterior draws for each, for overlaying densities.
#' @param x a did_input object with components.
#' @param priors list of \code{prior_beta()} objects (defaults to four sensible ones).
#' @param draws,level,rho,seed passed to \code{did_bayes}.
#' @return list(table, draws) with class "did_prior_sensitivity".
#' @export
did_bayes_prior_sensitivity <- function(x, priors = NULL, draws = 20000,
                                        level = 0.95, rho = NULL, seed = 1) {
  if (is.null(priors)) priors <- .default_priors()
  rows <- list(); dl <- list()
  for (p in priors) {
    f <- did_bayes(x, prior = p, draws = draws, level = level, rho = rho, seed = seed)
    f$prior <- p$name
    f$prob_negative <- attr(f, "prob_negative")
    rows[[p$name]] <- f[, c("prior", "estimate", "lower", "upper", "prob_negative")]
    dl[[p$name]] <- attr(f, "draws")
  }
  structure(list(table = do.call(rbind, rows), draws = dl),
            class = "did_prior_sensitivity")
}
