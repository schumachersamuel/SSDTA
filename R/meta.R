# ---------------------------------------------------------------------------
# Stratum-specific meta-analysis of sensitivity.
#
# Produces posterior (or approximate posterior) draws of sensitivity for each
# stratum, in the long form (category, draw, sens_pooled[, sens_predicted]) so
# the output feeds ssdta_standardize() directly. One separate model per stratum
# (default) -- inter-stratum heterogeneity is treated as structural, and cross-
# stratum partial pooling is off by default (per design decision).
#
# Engines:
#   "stan"       -- Bayesian binomial-logit random-effects (brms/Stan). The
#                   paper's method. Requires the 'brms' package.
#   "stan_joint" -- one Bayesian model across all strata with a shared between-
#                   study SD (tau): TP | trials(n) ~ 0 + stratum + (1 | study).
#                   Partially pools heterogeneity across strata, which matters
#                   when k = 3-5 studies leaves tau weakly identified per
#                   stratum. Also needs brms.
#   "glmm"       -- frequentist logit-normal random-effects via metafor; draws
#                   are a normal approximation on the logit scale. Fast.
#   "conjugate"  -- complete-pooling Beta-Binomial (sum counts, Beta posterior).
#                   Fastest, no dependencies; ignores between-study heterogeneity.
# ---------------------------------------------------------------------------

#' Stratum-specific sensitivity meta-analysis
#'
#' @param data data.frame with columns \code{study}, \code{stratum}, \code{tp},
#'   \code{fn} (true-positive / false-negative counts per study x stratum).
#'   \code{stratum} may be a factor; its level order is preserved.
#' @param engine "stan" (default; Bayesian RE per stratum, needs brms),
#'   "stan_joint" (one model across strata with a shared tau, needs brms),
#'   "glmm" (frequentist RE via metafor), or "conjugate" (complete-pooling
#'   Beta-Binomial).
#' @param n_draws number of posterior/approximate-posterior draws per stratum.
#' @param prior_intercept,prior_sd (stan engines only) prior specs. Defaults
#'   reproduce the paper: intercept Normal(0, 1.5), between-study SD
#'   Exponential(1). Pass alternatives for prior-sensitivity analysis, e.g.
#'   \code{prior_intercept = "student_t(3, 0, 1)"}, or sweep a set of them with
#'   \code{\link{ssdta_meta_sensitivity}}.
#' @param seed RNG seed.
#' @param ... passed to \code{brms::brm} (stan engines).
#' @return an \code{ssdta_meta} object: a list with \code{draws} (long
#'   data.frame: category, draw, sens_pooled, sens_predicted) and \code{summary}
#'   (per-stratum median + 95\% CrI). The \code{draws} element is a drop-in for
#'   \code{ssdta_standardize()}.
#' @export
ssdta_meta <- function(data, engine = c("stan", "glmm", "conjugate", "stan_joint"),
                       n_draws = 4000,
                       prior_intercept = "normal(0, 1.5)",
                       prior_sd = "exponential(1)",
                       seed = 1, ...) {
  engine <- match.arg(engine)
  need <- c("study", "stratum", "tp", "fn")
  if (!is.data.frame(data) || !all(need %in% names(data)))
    stop("`data` must be a data.frame with columns: ", paste(need, collapse = ", "), ".")
  if (!nrow(data)) stop("`data` must contain at least one row.")
  if (!is.numeric(n_draws) || length(n_draws) != 1L || !is.finite(n_draws) ||
      n_draws < 1 || n_draws != as.integer(n_draws))
    stop("`n_draws` must be a positive integer.")
  if (anyNA(data[, need]) || any(!is.finite(data$tp)) || any(!is.finite(data$fn)) ||
      any(data$tp < 0 | data$fn < 0))
    stop("study/stratum cannot be missing; tp/fn must be finite and non-negative.")
  if (any(data$tp != floor(data$tp) | data$fn != floor(data$fn)))
    stop("tp and fn must be whole-number counts.")
  if (any((data$tp + data$fn) == 0))
    warning("Dropping study-by-stratum rows with tp + fn = 0.")
  data <- data[(data$tp + data$fn) > 0, , drop = FALSE]
  if (!nrow(data)) stop("No informative rows remain after dropping empty cells.")
  strata <- if (is.factor(data$stratum)) levels(droplevels(data$stratum)) else unique(data$stratum)

  if (identical(engine, "stan_joint")) {
    jf <- .with_seed(seed, .fit_joint_stan(data, strata, n_draws = n_draws,
            prior_intercept = prior_intercept, prior_sd = prior_sd, seed = seed, ...))
    per <- lapply(strata, function(k)
      data.frame(category = k, draw = seq_len(n_draws),
                 sens_pooled = jf$pooled[[k]], sens_predicted = jf$predicted[[k]]))
  } else {
    fit_one <- switch(engine,
      conjugate = .fit_stratum_conjugate,
      glmm      = .fit_stratum_glmm,
      stan      = .fit_stratum_stan)

    per <- lapply(strata, function(k) {
      d <- data[as.character(data$stratum) == k, , drop = FALSE]
      dr <- .with_seed(seed + match(k, strata) - 1L,
        fit_one(d, n_draws = n_draws, prior_intercept = prior_intercept,
                prior_sd = prior_sd, seed = seed + match(k, strata) - 1L, ...))
      data.frame(category = k, draw = seq_len(n_draws),
                 sens_pooled = dr$pooled, sens_predicted = dr$predicted)
    })
  }
  draws <- do.call(rbind, per)
  draws$category <- factor(draws$category, levels = strata)

  summ <- do.call(rbind, lapply(split(draws, draws$category), function(g) {
    s <- .summ_draws(g$sens_pooled)
    data.frame(Category = g$category[1], Sensitivity = s["median"],
               Lower_95_CrI = s["lo"], Upper_95_CrI = s["hi"], row.names = NULL)
  }))
  structure(list(draws = draws, summary = summ, engine = engine,
                 n_studies = length(unique(as.character(data$study)))),
            class = "ssdta_meta")
}

#' @export
print.ssdta_meta <- function(x, ...) {
  cat(sprintf("Stratum-specific sensitivity meta-analysis (engine: %s, %d studies)\n",
              x$engine, x$n_studies))
  s <- x$summary
  for (i in seq_len(nrow(s)))
    cat(sprintf("  %-13s %5.1f%%   95%% CrI (%.1f, %.1f)\n",
                as.character(s$Category[i]), 100 * s$Sensitivity[i],
                100 * s$Lower_95_CrI[i], 100 * s$Upper_95_CrI[i]))
  invisible(x)
}

# --- engines ---------------------------------------------------------------

.fit_stratum_conjugate <- function(d, n_draws, ...) {
  d <- d[(d$tp + d$fn) > 0, , drop = FALSE]
  if (!nrow(d)) stop("stratum has no informative data.")
  TP <- sum(d$tp); FN <- sum(d$fn)
  s <- stats::rbeta(n_draws, 0.5 + TP, 0.5 + FN)   # Jeffreys prior, complete pooling
  list(pooled = s, predicted = s)                  # no heterogeneity model
}

.fit_stratum_glmm <- function(d, n_draws, ...) {
  if (!requireNamespace("metafor", quietly = TRUE))
    stop("engine = 'glmm' needs the 'metafor' package.")
  d <- d[(d$tp + d$fn) > 0, , drop = FALSE]        # drop empty (n = 0) study cells
  if (nrow(d) == 0) stop("stratum has no data after dropping empty cells.")
  if (nrow(d) == 1L)                               # single study -> conjugate fallback
    return(.fit_stratum_conjugate(d, n_draws))
  es <- metafor::escalc(measure = "PLO", xi = d$tp, ni = d$tp + d$fn)  # logit props
  fit <- tryCatch(metafor::rma(yi = es$yi, vi = es$vi, method = "REML"),
                  error = function(e) metafor::rma(yi = es$yi, vi = es$vi, method = "DL"))
  mu <- as.numeric(fit$b); se <- as.numeric(fit$se); tau <- sqrt(as.numeric(fit$tau2))
  logit_pooled    <- stats::rnorm(n_draws, mu, se)                 # uncertainty in the mean
  logit_predicted <- stats::rnorm(n_draws, mu, sqrt(se^2 + tau^2)) # predictive (new study)
  list(pooled = stats::plogis(logit_pooled), predicted = stats::plogis(logit_predicted))
}

.fit_stratum_stan <- function(d, n_draws, prior_intercept, prior_sd, seed = 1, ...) {
  if (!requireNamespace("brms", quietly = TRUE))
    stop("engine = 'stan' needs the 'brms' package (Stan backend). ",
         "Use engine = 'glmm' or 'conjugate' for a dependency-light alternative.")
  d$n <- d$tp + d$fn
  # single-study strata cannot fit a (1|study) RE term -> fall back to no-RE
  has_re <- length(unique(d$study)) > 1
  form <- if (has_re) brms::bf(tp | trials(n) ~ 1 + (1 | study)) else brms::bf(tp | trials(n) ~ 1)
  pr <- brms::prior_string(prior_intercept, class = "Intercept")
  if (has_re) pr <- pr + brms::prior_string(prior_sd, class = "sd")
  fit <- brms::brm(form, data = d, family = brms::brmsfamily("binomial", link = "logit"),
                   prior = pr, chains = 4, iter = 2000, refresh = 0, silent = 2,
                   seed = seed, ...)
  post <- brms::as_draws_df(fit)
  b <- post$b_Intercept
  idx <- sample(length(b), n_draws, replace = length(b) < n_draws)
  b <- b[idx]
  sdre <- if (has_re && !is.null(post$sd_study__Intercept)) post$sd_study__Intercept[idx] else 0
  list(pooled = stats::plogis(b),
       predicted = stats::plogis(b + stats::rnorm(n_draws, 0, sdre)))
}

# Joint model across strata with a shared between-study SD. Returns per-stratum
# draw vectors that share ONE posterior draw index: the joint fit induces cross-
# stratum correlation (shared tau, shared study effects) and dropping it would
# misstate uncertainty in the weighted sum computed downstream by
# ssdta_standardize(), which pairs strata by draw id.
.fit_joint_stan <- function(d, strata, n_draws, prior_intercept, prior_sd,
                            seed = 1, ...) {
  if (!requireNamespace("brms", quietly = TRUE))
    stop("engine = 'stan_joint' needs the 'brms' package (Stan backend). ",
         "Use engine = 'glmm' or 'conjugate' for a dependency-light alternative.")
  d$n <- d$tp + d$fn
  d$stratum <- factor(as.character(d$stratum), levels = strata)
  has_re <- length(unique(d$study)) > 1
  form <- if (has_re) brms::bf(tp | trials(n) ~ 0 + stratum + (1 | study))
          else brms::bf(tp | trials(n) ~ 0 + stratum)
  pr <- brms::prior_string(prior_intercept, class = "b")
  if (has_re) pr <- pr + brms::prior_string(prior_sd, class = "sd")
  fit <- brms::brm(form, data = d, family = brms::brmsfamily("binomial", link = "logit"),
                   prior = pr, chains = 4, iter = 4000, refresh = 0, silent = 2,
                   seed = seed, control = list(adapt_delta = 0.95), ...)
  post <- brms::as_draws_df(fit)

  # brms strips non-alphanumeric characters from coefficient names, so match on
  # a normalised key rather than assuming the mangling rule. Any stratum that
  # fails to map is an error -- silently dropping one would quietly change the
  # estimand of every downstream standardization.
  cols <- grep("^b_stratum", colnames(post), value = TRUE)
  norm <- function(z) tolower(gsub("[^A-Za-z0-9]", "", z))
  idx <- match(norm(paste0("b_stratum", strata)), norm(cols))
  if (anyNA(idx))
    stop("Could not map strata to joint-model coefficients: ",
         paste(strata[is.na(idx)], collapse = ", "),
         ". Coefficients found: ", paste(cols, collapse = ", "), ".")

  n_post <- nrow(post)
  sel <- sample(n_post, n_draws, replace = n_post < n_draws)
  sdre <- if (has_re && "sd_study__Intercept" %in% colnames(post))
    post$sd_study__Intercept[sel] else rep(0, n_draws)
  pooled <- predicted <- stats::setNames(vector("list", length(strata)), strata)
  for (i in seq_along(strata)) {
    b <- post[[cols[idx[i]]]][sel]
    pooled[[i]] <- stats::plogis(b)
    predicted[[i]] <- stats::plogis(b + stats::rnorm(n_draws, 0, sdre))
  }
  list(pooled = pooled, predicted = predicted)
}


# ---------------------------------------------------------------------------
# Specification sweep: refit the stratum meta-analysis under alternative priors,
# model structures, and study exclusions, and tabulate them side by side. This
# is the transportability / robustness half of the framework -- the standardized
# estimate is only as defensible as the S_k that feed it.
# ---------------------------------------------------------------------------

#' One arm of a meta-analysis specification sweep
#'
#' @param label unique display name for this specification.
#' @param engine optional engine override (defaults to the sweep's engine).
#' @param prior_intercept,prior_sd optional prior overrides (stan engines only).
#' @param exclude optional character vector of \code{study} values to drop --
#'   for leave-one-study-out and influence checks.
#' @export
meta_spec <- function(label, engine = NULL, prior_intercept = NULL,
                      prior_sd = NULL, exclude = NULL) {
  if (!is.character(label) || length(label) != 1L || is.na(label) || !nzchar(label))
    stop("`label` must be a single non-empty string.")
  if (!is.null(engine) && !engine %in% c("stan", "stan_joint", "glmm", "conjugate"))
    stop("`engine` must be one of stan, stan_joint, glmm, conjugate.")
  if (!is.null(exclude) && (!is.character(exclude) || !length(exclude) || anyNA(exclude)))
    stop("`exclude` must be a non-empty character vector of study names.")
  structure(list(label = label, engine = engine, prior_intercept = prior_intercept,
                 prior_sd = prior_sd, exclude = exclude), class = "meta_spec")
}

# The paper's sensitivity analyses 8a, 8b, 8d, 8e. 8c (excluding one study) is
# the leave-one-out family, reachable generically via `loo = TRUE`.
.default_meta_specs <- function(engine) {
  if (!engine %in% c("stan", "stan_joint")) return(list(meta_spec("Primary")))
  list(
    meta_spec("Primary"),
    meta_spec("N(0,5) intercept",         prior_intercept = "normal(0, 5)"),
    meta_spec("Joint model (shared tau)", engine = "stan_joint"),
    meta_spec("Logistic(0,1) intercept",  prior_intercept = "logistic(0, 1)"),
    meta_spec("Half-t(3,0,2.5) tau",      prior_sd = "student_t(3, 0, 2.5)")
  )
}

#' Prior- and specification-sensitivity sweep of the stratum meta-analysis
#'
#' Refits \code{\link{ssdta_meta}} under a set of specifications and returns
#' their stratum-specific sensitivities side by side. With \code{engine = "stan"}
#' the defaults reproduce the paper's sensitivity analyses: an N(0,5) intercept
#' prior, a joint model with shared tau, a Logistic(0,1) intercept prior, and a
#' Half-t(3,0,2.5) prior on tau. \code{loo = TRUE} adds one leave-one-study-out
#' arm per study.
#'
#' Priors matter here more than they usually would: with k = 3-5 studies per
#' stratum, tau is weakly identified, so the between-study SD prior does real
#' work in the sparse strata (Very low, Trace, Not-detected).
#'
#' @param data per-study counts, as for \code{\link{ssdta_meta}}.
#' @param specs a \code{\link{meta_spec}} or list of them; \code{NULL} uses the
#'   defaults for \code{engine}.
#' @param engine base engine for specs that do not override it.
#' @param n_draws,seed,... passed to \code{\link{ssdta_meta}}.
#' @param level credible level for the reported intervals.
#' @param loo append one leave-one-study-out specification per study.
#' @param prior_intercept,prior_sd baseline priors for specs that do not
#'   override them.
#' @return an \code{ssdta_meta_sensitivity} object: \code{table} (long: analysis,
#'   stratum, sensitivity, lo, hi) and \code{fits} (the \code{ssdta_meta}
#'   results). \code{as.data.frame()} gives the wide comparison table.
#' @export
ssdta_meta_sensitivity <- function(data, specs = NULL,
                                   engine = c("stan", "glmm", "conjugate", "stan_joint"),
                                   n_draws = 4000, level = 0.95, loo = FALSE,
                                   prior_intercept = "normal(0, 1.5)",
                                   prior_sd = "exponential(1)",
                                   seed = 1, ...) {
  engine <- match.arg(engine)
  .validate_level(level)
  need <- c("study", "stratum", "tp", "fn")
  if (!is.data.frame(data) || !all(need %in% names(data)))
    stop("`data` must be a data.frame with columns: ", paste(need, collapse = ", "), ".")
  if (is.null(specs)) specs <- .default_meta_specs(engine)
  if (inherits(specs, "meta_spec")) specs <- list(specs)
  if (!is.list(specs) || !length(specs) ||
      !all(vapply(specs, inherits, logical(1), what = "meta_spec")))
    stop("`specs` must be a meta_spec, or a non-empty list of meta_spec objects.")

  if (isTRUE(loo)) {
    studies <- unique(as.character(data$study))
    if (length(studies) < 3)
      stop("Leave-one-study-out needs at least 3 studies; `data` has ", length(studies), ".")
    specs <- c(specs, lapply(studies, function(s) meta_spec(paste0("Excl. ", s), exclude = s)))
  }
  labs <- vapply(specs, function(s) s$label, character(1))
  if (anyDuplicated(labs))
    stop("`specs` labels must be unique; duplicated: ",
         paste(unique(labs[duplicated(labs)]), collapse = ", "), ".")

  fits <- stats::setNames(vector("list", length(specs)), labs)
  rows <- vector("list", length(specs))
  for (i in seq_along(specs)) {
    sp <- specs[[i]]
    eng <- if (is.null(sp$engine)) engine else sp$engine
    if (!eng %in% c("stan", "stan_joint") &&
        (!is.null(sp$prior_intercept) || !is.null(sp$prior_sd)))
      warning("Spec '", sp$label, "' sets priors, but engine '", eng,
              "' ignores them -- this arm will duplicate the primary fit.")
    d <- data
    if (!is.null(sp$exclude)) {
      keep <- !as.character(d$study) %in% sp$exclude
      if (!any(keep)) stop("Spec '", sp$label, "' excludes every study.")
      if (all(keep))
        warning("Spec '", sp$label, "' excludes no study present in `data`.")
      d <- d[keep, , drop = FALSE]
    }
    f <- ssdta_meta(d, engine = eng, n_draws = n_draws,
                    prior_intercept = if (is.null(sp$prior_intercept)) prior_intercept
                                      else sp$prior_intercept,
                    prior_sd = if (is.null(sp$prior_sd)) prior_sd else sp$prior_sd,
                    seed = seed, ...)
    fits[[i]] <- f
    rows[[i]] <- do.call(rbind, lapply(split(f$draws, f$draws$category), function(g) {
      if (!nrow(g)) return(NULL)
      s <- .summ_draws(g$sens_pooled, level)
      data.frame(analysis = sp$label, stratum = as.character(g$category[1]),
                 sensitivity = s[["median"]], lo = s[["lo"]], hi = s[["hi"]],
                 row.names = NULL, stringsAsFactors = FALSE)
    }))
  }
  tab <- do.call(rbind, rows)
  tab$analysis <- factor(tab$analysis, levels = labs)
  tab$stratum <- .ordered_strata(tab$stratum)
  structure(list(table = tab, fits = fits, level = level, engine = engine),
            class = "ssdta_meta_sensitivity")
}

# wide "median (lo-hi)" comparison table -- the paper's table_sensitivity_analyses
.meta_sensitivity_wide <- function(x) {
  t <- x$table
  t$cell <- sprintf("%.1f (%.1f-%.1f)", 100 * t$sensitivity, 100 * t$lo, 100 * t$hi)
  strata <- levels(droplevels(t$stratum))
  out <- data.frame(Stratum = strata, stringsAsFactors = FALSE)
  for (a in levels(t$analysis)) {
    sub <- t[t$analysis == a, , drop = FALSE]
    out[[a]] <- sub$cell[match(strata, as.character(sub$stratum))]
  }
  out
}

#' @export
as.data.frame.ssdta_meta_sensitivity <- function(x, ..., wide = TRUE) {
  if (isTRUE(wide)) .meta_sensitivity_wide(x) else x$table
}

#' @export
summary.ssdta_meta_sensitivity <- function(object, ...) .meta_sensitivity_wide(object)

#' @export
plot.ssdta_meta_sensitivity <- function(x, ...) plot_meta_sensitivity(x, ...)

#' @export
print.ssdta_meta_sensitivity <- function(x, ...) {
  cat(sprintf("Meta-analysis specification sweep -- %d arms, %d%% CrI (base engine: %s)\n",
              length(x$fits), round(100 * x$level), x$engine))
  print(.meta_sensitivity_wide(x), row.names = FALSE)
  invisible(x)
}
