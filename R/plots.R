# ---------------------------------------------------------------------------
# Spectrum-standardization visualisations (ggplot2). All return ggplot objects
# with colour-blind-safe, report-size defaults and remain fully themeable.
# ---------------------------------------------------------------------------

.stratum_levels <- c("High", "Medium", "Low", "Very low", "Trace", "Not-detected")

.plot_palette <- function(values, known = .stratum_palette) {
  values <- unique(as.character(values))
  out <- known[intersect(names(known), values)]
  missing <- setdiff(values, names(out))
  if (length(missing)) {
    extra <- grDevices::hcl.colors(length(missing), palette = "Dark 3")
    out <- c(out, stats::setNames(extra, missing))
  }
  out[values]
}

.ordered_strata <- function(x) {
  x <- as.character(x)
  known <- intersect(.stratum_levels, unique(x))
  factor(x, levels = c(known, setdiff(unique(x), known)))
}

.pi_for_profiles <- function(pi, profiles) {
  n <- length(profiles); nms <- names(profiles)
  if (is.null(pi)) return(rep(list(NULL), n))
  if (is.list(pi)) {
    if (!is.null(names(pi)) && all(nms %in% names(pi))) return(pi[nms])
    if (length(pi) == n) return(pi)
  }
  if (length(pi) == 1L) return(rep(list(pi), n))
  if (is.numeric(pi) && length(pi) == n) return(as.list(pi))
  if (!is.null(names(pi)) && all(nms %in% names(pi))) return(as.list(pi[nms]))
  stop("`pi` must be NULL, a scalar, or one value per profile (optionally named).")
}

#' Case-mix (burden distribution) across target populations
#'
#' Stacked bars of stratum proportions. By default only detected strata are
#' shown. Supply \code{pi} to add the Not-detected stratum and display the full
#' target distribution; use \code{"observed"} for a profile with observed
#' Not-detected counts.
#' @param profiles a named list of \code{burden_profile} objects.
#' @param pi optional scalar or one value per profile; numeric values must be in
#'   [0,1], and \code{"observed"} uses an observed Not-detected count.
#' @export
plot_casemix <- function(profiles, pi = NULL) {
  if (!is.list(profiles) || !length(profiles) ||
      !all(vapply(profiles, inherits, logical(1), what = "burden_profile")))
    stop("`profiles` must be a non-empty list of burden_profile objects.")
  if (is.null(names(profiles))) names(profiles) <- paste0("profile_", seq_along(profiles))
  pis <- .pi_for_profiles(pi, profiles)
  df <- do.call(rbind, lapply(seq_along(profiles), function(i) {
    p <- profiles[[i]]; cd <- p$counts_detected
    if (is.null(pis[[i]])) {
      prop <- cd / sum(cd)
    } else if (identical(pis[[i]], "observed")) {
      if (is.na(p$count_notdetected))
        stop('pi = "observed" needs count_notdetected for profile ', names(profiles)[i], ".")
      prop <- c(cd, `Not-detected` = p$count_notdetected)
      prop <- prop / sum(prop)
    } else {
      if (!is.numeric(pis[[i]]) || length(pis[[i]]) != 1L ||
          !is.finite(pis[[i]]) || pis[[i]] < 0 || pis[[i]] > 1)
        stop("numeric pi values must be finite scalars in [0, 1].")
      prop <- c(cd / sum(cd) * (1 - pis[[i]]), `Not-detected` = pis[[i]])
    }
    data.frame(
      label = sprintf("%s%s\n(n=%s)", p$setting,
                      if (!is.na(p$country)) paste0(", ", p$country) else "",
                      format(p$n, big.mark = ",")),
      stratum = names(prop), prop = as.numeric(prop), stringsAsFactors = FALSE)
  }))
  df$stratum <- .ordered_strata(df$stratum)
  df$label <- factor(df$label, levels = unique(df$label))
  pal <- .plot_palette(levels(df$stratum))

  ggplot2::ggplot(df, ggplot2::aes(.data$label, .data$prop, fill = .data$stratum)) +
    ggplot2::geom_col(width = 0.68, colour = "white", linewidth = 0.25) +
    ggplot2::scale_fill_manual(values = pal, drop = FALSE, name = "Burden stratum") +
    ggplot2::scale_y_continuous(labels = function(v) paste0(round(100 * v), "%"),
                               breaks = seq(0, 1, 0.2), expand = c(0, 0)) +
    ggplot2::labs(x = NULL, y = if (is.null(pi)) "Share of detected cases" else "Share of true cases",
                  title = "Case mix across target populations",
                  subtitle = if (is.null(pi))
                    "Detected strata only; supply pi to show the full target distribution" else
                    "Full target distribution, including the Not-detected stratum") +
    .ssdta_theme() +
    ggplot2::theme(panel.grid.major.x = ggplot2::element_blank())
}

#' Stratum-specific sensitivity with nested posterior intervals
#'
#' @param sens_draws long draws data.frame (category, draw, sens_pooled), e.g.
#'   \code{npoc_posterior_draws} or \code{ssdta_meta(...)$draws}.
#' @param level outer credible level; the inner interval is 80\%.
#' @export
plot_stratum_sensitivity <- function(sens_draws, level = 0.95) {
  .validate_level(level)
  need <- c("category", "sens_pooled")
  if (!is.data.frame(sens_draws) || !all(need %in% names(sens_draws)))
    stop("`sens_draws` must contain category and sens_pooled columns.")
  a <- (1 - level) / 2
  df <- do.call(rbind, lapply(split(sens_draws, sens_draws$category), function(g) {
    if (!nrow(g)) return(NULL)
    q <- stats::quantile(g$sens_pooled, c(a, 0.10, 0.5, 0.90, 1 - a), names = FALSE)
    data.frame(stratum = as.character(g$category[1]), median = q[3],
               lo = q[1], lo80 = q[2], hi80 = q[4], hi = q[5])
  }))
  ord <- rev(levels(.ordered_strata(df$stratum)))
  df$stratum <- factor(df$stratum, levels = ord)
  pal <- .plot_palette(as.character(df$stratum))

  ggplot2::ggplot(df, ggplot2::aes(100 * .data$median, .data$stratum,
                                   colour = .data$stratum)) +
    ggplot2::geom_errorbar(ggplot2::aes(xmin = 100 * .data$lo, xmax = 100 * .data$hi),
                           orientation = "y", width = 0, linewidth = 0.55) +
    ggplot2::geom_errorbar(ggplot2::aes(xmin = 100 * .data$lo80, xmax = 100 * .data$hi80),
                           orientation = "y", width = 0, linewidth = 2.5) +
    ggplot2::geom_point(size = 3.1, shape = 21, fill = "white", stroke = 1.1) +
    ggplot2::scale_colour_manual(values = pal, guide = "none") +
    ggplot2::scale_x_continuous(breaks = seq(0, 100, 20)) +
    ggplot2::coord_cartesian(xlim = c(0, 100)) +
    ggplot2::labs(x = "Sensitivity (%)", y = NULL, title = "Stratum-specific sensitivity",
                  subtitle = sprintf("Posterior median; thick 80%% and thin %d%% credible intervals",
                                     round(100 * level))) +
    .ssdta_theme()
}

# Coerce whatever the caller passes as `pooled` into a common summary frame.
.pooled_summary <- function(pooled, level, label = "Pooled") {
  if (is.list(pooled) && !is.data.frame(pooled) && !is.null(pooled$draws))
    pooled <- pooled$draws
  if (is.data.frame(pooled) && all(c("category", "sens_pooled") %in% names(pooled))) {
    out <- do.call(rbind, lapply(split(pooled, pooled$category), function(g) {
      if (!nrow(g)) return(NULL)
      s <- .summ_draws(g$sens_pooled, level)
      data.frame(stratum = as.character(g$category[1]), label = label,
                 est = s[["median"]], lo = s[["lo"]], hi = s[["hi"]],
                 n = NA_real_, tp = NA_real_, count_label = NA_character_,
                 type = "Pooled", stringsAsFactors = FALSE)
    }))
    return(out)
  }
  if (is.data.frame(pooled) &&
      all(c("Category", "Sensitivity", "Lower_95_CrI", "Upper_95_CrI") %in% names(pooled))) {
    if (!isTRUE(all.equal(level, 0.95)))
      warning("`pooled` is a 95% summary table; ignoring level = ", level,
              ". Pass posterior draws to honour a different level.")
    return(data.frame(stratum = as.character(pooled$Category), label = label,
                      est = pooled$Sensitivity, lo = pooled$Lower_95_CrI,
                      hi = pooled$Upper_95_CrI, n = NA_real_, tp = NA_real_,
                      count_label = NA_character_, type = "Pooled",
                      stringsAsFactors = FALSE))
  }
  stop("`pooled` must be an ssdta_meta() result, a draws data.frame (category, ",
       "sens_pooled), or a summary table (Category, Sensitivity, Lower_95_CrI, ",
       "Upper_95_CrI).")
}

#' Study-level forest plot of sensitivity within each stratum
#'
#' Per-study sensitivity with Wilson intervals, faceted by burden stratum, with
#' the pooled estimate as a diamond. This is the transportability diagnostic
#' behind the whole method: standardization assumes stratum-specific sensitivity
#' is constant across settings, and a stratum where the study points disagree
#' more than their intervals allow is exactly where that assumption fails.
#'
#' @param data data.frame with columns \code{study}, \code{stratum}, \code{tp},
#'   \code{fn} (as accepted by \code{\link{ssdta_meta}}). Cells with
#'   \code{tp + fn == 0} are dropped and reported in the caption.
#' @param pooled optional pooled estimates to overlay as diamonds: an
#'   \code{\link{ssdta_meta}} result, a draws data.frame (\code{category},
#'   \code{sens_pooled}) such as \code{npoc_posterior_draws}, or a summary table.
#'   \code{NULL} plots study estimates only.
#' @param level interval level for the study Wilson intervals, and for the
#'   pooled credible intervals when \code{pooled} carries draws.
#' @param ncol facet columns.
#' @param counts annotate each study row with its \code{tp/n}.
#' @param pooled_label row label for the pooled estimate.
#' @export
plot_forest_stratum <- function(data, pooled = NULL, level = 0.95, ncol = 2,
                                counts = TRUE, pooled_label = "Pooled") {
  .validate_level(level)
  need <- c("study", "stratum", "tp", "fn")
  if (!is.data.frame(data) || !all(need %in% names(data)))
    stop("`data` must be a data.frame with columns: ", paste(need, collapse = ", "), ".")
  if (anyNA(data[, need])) stop("`data` must not contain missing study/stratum/tp/fn.")
  if (any(data$tp < 0 | data$fn < 0) || any(!is.finite(data$tp)) || any(!is.finite(data$fn)))
    stop("tp and fn must be finite and non-negative.")

  d <- data
  d$n <- d$tp + d$fn
  n_dropped <- sum(d$n == 0)
  d <- d[d$n > 0, , drop = FALSE]
  if (!nrow(d)) stop("No study-by-stratum cells with n > 0.")
  ci <- .wilson_ci(d$tp, d$n, level)
  df <- data.frame(
    stratum = as.character(d$stratum), label = as.character(d$study),
    est = d$tp / d$n, lo = ci[, "lower"], hi = ci[, "upper"],
    n = as.double(d$n), tp = as.double(d$tp),
    count_label = sprintf("%d/%d", as.integer(d$tp), as.integer(d$n)),
    type = "Study", stringsAsFactors = FALSE)

  if (!is.null(pooled)) {
    ps <- .pooled_summary(pooled, level, pooled_label)
    ps <- ps[ps$stratum %in% df$stratum, , drop = FALSE]
    if (!nrow(ps))
      warning("`pooled` shares no stratum with `data`; plotting study estimates only.")
    df <- rbind(df, ps)
  }
  df$stratum <- .ordered_strata(df$stratum)
  # First level sits at the bottom of a discrete axis: pooled first, then
  # studies ascending in n, so the largest study is at the top of each facet.
  df$sort_val <- ifelse(df$type == "Pooled", -Inf, df$n)
  df$y_key <- sprintf("%02d|%s", as.integer(df$stratum), df$label)
  df$y_key <- factor(df$y_key,
                     levels = unique(df$y_key[order(df$stratum, df$sort_val, df$label)]))
  x_max <- if (isTRUE(counts)) 1.17 else 1

  p <- ggplot2::ggplot(df, ggplot2::aes(.data$est, .data$y_key)) +
    ggplot2::geom_errorbar(ggplot2::aes(xmin = .data$lo, xmax = .data$hi),
                           orientation = "y", width = 0.22, linewidth = 0.45,
                           colour = "grey45") +
    ggplot2::geom_point(ggplot2::aes(shape = .data$type, size = .data$type,
                                     colour = .data$type))
  if (isTRUE(counts))
    p <- p + ggplot2::geom_text(
      data = df[df$type == "Study", , drop = FALSE],
      ggplot2::aes(x = 1.02, label = .data$count_label),
      hjust = 0, size = 2.7, colour = "grey40")

  p +
    ggplot2::scale_shape_manual(values = c(Study = 16, Pooled = 18), name = NULL) +
    ggplot2::scale_size_manual(values = c(Study = 2, Pooled = 3.6), name = NULL) +
    ggplot2::scale_colour_manual(values = c(Study = "grey25", Pooled = "#0072B2"),
                                 name = NULL) +
    ggplot2::scale_x_continuous(limits = c(0, x_max), breaks = seq(0, 1, 0.2),
                                labels = function(v) paste0(round(100 * v), "%")) +
    ggplot2::scale_y_discrete(labels = function(z) sub("^[^|]+\\|", "", z)) +
    ggplot2::facet_wrap(~ stratum, scales = "free_y", ncol = ncol) +
    ggplot2::labs(
      x = "Sensitivity", y = NULL,
      title = "Study-level sensitivity by burden stratum",
      subtitle = sprintf("Circles: study estimate (Wilson %d%% CI)%s",
                         round(100 * level),
                         if (is.null(pooled)) "" else
                           sprintf("  |  Diamond: %s", tolower(pooled_label))),
      caption = if (n_dropped > 0)
        sprintf("%d study-by-stratum cell(s) with no observations omitted.", n_dropped)
        else NULL) +
    .ssdta_theme() +
    ggplot2::theme(strip.text = ggplot2::element_text(face = "bold"),
                   panel.grid.major.x = ggplot2::element_blank())
}

#' Compare stratum sensitivities across meta-analysis specifications
#'
#' Dot-and-whisker plot of the specification sweep from
#' \code{\link{ssdta_meta_sensitivity}} -- one row per stratum, one dodged point
#' per specification.
#'
#' @param sens an \code{ssdta_meta_sensitivity} object.
#' @param dodge vertical dodge width between specifications.
#' @param ... ignored.
#' @export
plot_meta_sensitivity <- function(sens, dodge = 0.72, ...) {
  if (!inherits(sens, "ssdta_meta_sensitivity"))
    stop("`sens` must be an ssdta_meta_sensitivity object.")
  df <- sens$table
  df$stratum <- factor(as.character(df$stratum),
                       levels = rev(levels(droplevels(df$stratum))))
  arms <- levels(df$analysis)
  pal <- .plot_palette(arms, known = c(Primary = "#2c3e50"))
  shapes <- stats::setNames(rep(c(16, 17, 15, 18, 8, 7, 3, 4),
                               length.out = length(arms)), arms)
  pd <- ggplot2::position_dodge(width = dodge)

  ggplot2::ggplot(df, ggplot2::aes(.data$sensitivity, .data$stratum,
                                   colour = .data$analysis, shape = .data$analysis)) +
    ggplot2::geom_errorbar(ggplot2::aes(xmin = .data$lo, xmax = .data$hi),
                           orientation = "y", width = 0.3, linewidth = 0.45,
                           position = pd) +
    ggplot2::geom_point(size = 2.4, position = pd) +
    ggplot2::scale_colour_manual(values = pal, name = NULL) +
    ggplot2::scale_shape_manual(values = shapes, name = NULL) +
    ggplot2::scale_x_continuous(breaks = seq(0, 1, 0.2),
                                labels = function(v) paste0(round(100 * v), "%")) +
    ggplot2::coord_cartesian(xlim = c(0, 1)) +
    ggplot2::labs(x = "Sensitivity", y = NULL,
                  title = "Stratum sensitivity across specifications",
                  subtitle = sprintf("Point estimate and %d%% CrI under each prior, model and exclusion",
                                     round(100 * sens$level))) +
    .ssdta_theme()
}

#' pi-sensitivity curves of standardized sensitivity across settings
#'
#' @param grids named list of pi-grid data.frames returned by
#'   \code{ssdta_standardize(..., pi = <vector>)}.
#' @export
plot_pi_sensitivity <- function(grids) {
  if (!is.list(grids) || !length(grids)) stop("`grids` must be a non-empty list.")
  if (is.null(names(grids))) names(grids) <- paste0("series_", seq_along(grids))
  df <- do.call(rbind, lapply(names(grids), function(nm) {
    g <- grids[[nm]]
    need <- c("setting", "pi", "median", "lo", "hi", "lo80", "hi80")
    if (!is.data.frame(g) || !all(need %in% names(g)))
      stop("Every grid must be a pi-grid result from ssdta_standardize().")
    g$series <- nm; g
  }))
  # Prefer human-readable setting names, but retain list names when settings collide.
  if (!anyDuplicated(unique(df[c("series", "setting")] )$setting)) df$series <- df$setting
  df$series <- factor(df$series, levels = unique(df$series))
  pal <- .plot_palette(levels(df$series), .setting_palette)

  ggplot2::ggplot(df, ggplot2::aes(.data$pi, 100 * .data$median,
                                   colour = .data$series, fill = .data$series)) +
    ggplot2::geom_ribbon(ggplot2::aes(ymin = 100 * .data$lo, ymax = 100 * .data$hi),
                         alpha = 0.10, colour = NA) +
    ggplot2::geom_ribbon(ggplot2::aes(ymin = 100 * .data$lo80, ymax = 100 * .data$hi80),
                         alpha = 0.19, colour = NA) +
    ggplot2::geom_line(linewidth = 0.95) +
    ggplot2::scale_colour_manual(values = pal, aesthetics = c("colour", "fill"), name = NULL) +
    ggplot2::scale_x_continuous(labels = function(v) paste0(round(100 * v), "%")) +
    ggplot2::labs(x = expression(pi ~ "(Not-detected share among true cases)"),
                  y = "Standardized sensitivity (%)",
                  title = "Sensitivity to the unobserved Not-detected weight",
                  subtitle = "Line: median; darker 80% and lighter 95% credible bands") +
    .ssdta_theme()
}

#' Forest of standardized sensitivity across settings
#'
#' @param stds list of \code{ssdta_std} objects (scalar-pi results).
#' @export
plot_standardized <- function(stds) {
  if (!is.list(stds) || !length(stds) ||
      !all(vapply(stds, inherits, logical(1), what = "ssdta_std")))
    stop("`stds` must be a non-empty list of ssdta_std objects.")
  df <- do.call(rbind, lapply(stds, function(s) {
    r <- s$summary
    ctry <- if (!is.null(s$profile$country) && !is.na(s$profile$country))
      paste0(", ", s$profile$country) else ""
    data.frame(label = sprintf("%s%s  (pi=%s)", r$setting, ctry,
                               if (is.numeric(s$pi)) sprintf("%.2f", s$pi) else s$pi),
               setting = r$setting, median = r$median, lo = r$lo, hi = r$hi,
               lo80 = r$lo80, hi80 = r$hi80)
  }))
  df$label <- make.unique(df$label)
  df$label <- factor(df$label, levels = rev(df$label))
  pal <- .plot_palette(unique(df$setting), .setting_palette)

  ggplot2::ggplot(df, ggplot2::aes(100 * .data$median, .data$label,
                                   colour = .data$setting)) +
    ggplot2::geom_errorbar(ggplot2::aes(xmin = 100 * .data$lo, xmax = 100 * .data$hi),
                           orientation = "y", width = 0, linewidth = 0.55) +
    ggplot2::geom_errorbar(ggplot2::aes(xmin = 100 * .data$lo80, xmax = 100 * .data$hi80),
                           orientation = "y", width = 0, linewidth = 2.6) +
    ggplot2::geom_point(size = 3.2, shape = 21, fill = "white", stroke = 1.1) +
    ggplot2::scale_colour_manual(values = pal, guide = "none") +
    ggplot2::scale_x_continuous(breaks = seq(0, 100, 20)) +
    ggplot2::coord_cartesian(xlim = c(0, 100)) +
    ggplot2::labs(x = "Standardized sensitivity (%)", y = NULL,
                  title = "Spectrum-standardized sensitivity by target population",
                  subtitle = "Posterior median; thick 80% and thin 95% credible intervals") +
    .ssdta_theme()
}

#' Show how case-mix weights produce standardized sensitivity
#'
#' A two-panel mechanism figure. The upper panel shows the full target case mix
#' (w[k]); the lower panel shows each stratum's contribution (w[k] * S[k]) to
#' standardized sensitivity. Points and intervals in the lower panel show the
#' resulting posterior S_std.
#'
#' @param sens_draws stratum-specific posterior draws accepted by
#'   \code{ssdta_standardize()}.
#' @param profiles named list of burden profiles.
#' @param pi scalar or one value per profile; may include \code{"observed"}.
#' @param level credible level.
#' @param seed reproducibility seed.
#' @export
plot_spectrum_mechanism <- function(sens_draws, profiles, pi, level = 0.95, seed = 42) {
  if (!is.list(profiles) || !length(profiles) ||
      !all(vapply(profiles, inherits, logical(1), what = "burden_profile")))
    stop("`profiles` must be a non-empty list of burden_profile objects.")
  if (is.null(names(profiles))) names(profiles) <- paste0("profile_", seq_along(profiles))
  pis <- .pi_for_profiles(pi, profiles)
  if (any(vapply(pis, is.null, logical(1)))) stop("`pi` is required for every profile.")

  pieces <- totals <- vector("list", length(profiles))
  for (i in seq_along(profiles)) {
    z <- ssdta_standardize(sens_draws, profiles[[i]], pis[[i]], level = level, seed = seed)
    label <- sprintf("%s%s", profiles[[i]]$setting,
                     if (!is.na(profiles[[i]]$country)) paste0(", ", profiles[[i]]$country) else "")
    w <- apply(z$weight_draws, 2, stats::median)
    contribution <- apply(z$weight_draws * z$sensitivity_draws, 2, stats::median)
    pieces[[i]] <- rbind(
      data.frame(label, panel = "A  Target case mix:  w[k]", stratum = names(w), value = w),
      data.frame(label, panel = "B  Sensitivity contribution:  w[k] x S[k]",
                 stratum = names(contribution), value = contribution))
    totals[[i]] <- data.frame(label, panel = "B  Sensitivity contribution:  w[k] x S[k]",
                              median = z$summary$median, lo = z$summary$lo, hi = z$summary$hi)
  }
  df <- do.call(rbind, pieces); total <- do.call(rbind, totals)
  df$stratum <- .ordered_strata(df$stratum)
  df$label <- factor(df$label, levels = unique(df$label))
  total$label <- factor(total$label, levels = levels(df$label))
  pal <- .plot_palette(levels(df$stratum))

  ggplot2::ggplot(df, ggplot2::aes(.data$label, .data$value, fill = .data$stratum)) +
    ggplot2::geom_col(width = 0.68, colour = "white", linewidth = 0.25) +
    ggplot2::geom_errorbar(data = total,
      ggplot2::aes(x = .data$label, y = .data$median, ymin = .data$lo, ymax = .data$hi),
      inherit.aes = FALSE, width = 0.12, linewidth = 0.65, colour = "black") +
    ggplot2::geom_point(data = total,
      ggplot2::aes(x = .data$label, y = .data$median), inherit.aes = FALSE,
      shape = 21, size = 2.7, fill = "white", colour = "black") +
    ggplot2::facet_wrap(~ panel, ncol = 1, scales = "free_y", strip.position = "top") +
    ggplot2::scale_fill_manual(values = pal, drop = FALSE, name = "Burden stratum") +
    ggplot2::scale_y_continuous(labels = function(v) paste0(round(100 * v), "%"),
                               expand = ggplot2::expansion(mult = c(0, 0.08))) +
    ggplot2::labs(x = NULL, y = "Population share / sensitivity contribution",
                  title = "How a case-mix shift changes expected sensitivity",
                  subtitle = "Reweighting moves mass toward strata with lower S[k]; point and whisker show S_std median and 95% CrI") +
    .ssdta_theme() +
    ggplot2::theme(strip.text.x = ggplot2::element_text(face = "bold", hjust = 0),
                   strip.background = ggplot2::element_blank(),
                   panel.grid.major.x = ggplot2::element_blank())
}
