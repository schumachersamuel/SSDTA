# ---------------------------------------------------------------------------
# ggplot2 visualisations. All return a ggplot object the caller can further
# theme or save with ggsave().
# ---------------------------------------------------------------------------

.pp <- function(v) 100 * v

#' Parallel-trends diagram for the anchored contrasts
#'
#' Draws Se(C) drifting between the two study populations and the counterfactual
#' parallel projection of Se(A); the gap to Se(B) is Delta_DiD. Needs the
#' component form.
#' @param x a did_input object with components.
#' @export
plot_parallel <- function(x) {
  stopifnot(inherits(x, "did_input"), !is.null(x$components))
  cp <- x$components; rownames(cp) <- cp$part; lab <- x$labels
  seA <- .pp(cp["A","sens"]);  seCa <- .pp(cp["Ca","sens"])
  seB <- .pp(cp["B","sens"]);  seCb <- .pp(cp["Cb","sens"])
  cf  <- seA + (seCb - seCa)                      # counterfactual A in pop b
  df_line <- rbind(
    data.frame(pop = c(1, 2), y = c(seCa, seCb), series = paste0("Anchor Se(", lab[3], ")")),
    data.frame(pop = c(1, 2), y = c(seA, cf),   series = paste0("Se(", lab[1], ") (obs + counterfactual)")),
    data.frame(pop = c(1, 2), y = c(NA, seB),   series = paste0("Se(", lab[2], ") (observed)"))
  )
  df_pt <- data.frame(
    pop = c(1, 1, 2, 2), y = c(seCa, seA, seCb, seB),
    series = c(paste0("Anchor Se(", lab[3], ")"),
               paste0("Se(", lab[1], ") (obs + counterfactual)"),
               paste0("Anchor Se(", lab[3], ")"),
               paste0("Se(", lab[2], ") (observed)"))
  )
  series <- unique(c(df_line$series, df_pt$series))
  pal <- stats::setNames(c("#333333", "#0072B2", "#D55E00")[seq_along(series)], series)
  ggplot() +
    geom_line(data = df_line[df_line$series != paste0("Se(", lab[2], ") (observed)"), ],
              aes(.data$pop, .data$y, colour = .data$series, linetype = .data$series),
              linewidth = 0.9, na.rm = TRUE) +
    geom_point(data = df_pt, aes(.data$pop, .data$y, colour = .data$series),
               size = 3, na.rm = TRUE) +
    annotate("segment", x = 2, xend = 2, y = cf, yend = seB,
             linewidth = 1.1, colour = "grey30",
             arrow = grid::arrow(ends = "both", length = grid::unit(0.15, "cm"))) +
    annotate("text", x = 2.04, y = (cf + seB) / 2,
             label = sprintf("Delta[DiD] == %+.1f", seB - cf),
             parse = TRUE, hjust = 0, size = 3.5) +
    scale_x_continuous(breaks = c(1, 2),
                       labels = c("Population a\n(A vs C studies)",
                                  "Population b\n(B vs C studies)"),
                       limits = c(0.9, 2.4)) +
    scale_colour_manual(values = pal) +
    labs(x = NULL, y = "Sensitivity (%)",
         title = "Parallel-trends view of the anchored difference-in-differences",
         colour = NULL, linetype = NULL) +
    .ssdta_theme()
}

#' Forest-style comparison of the interval methods
#'
#' @param fit an \code{anchored_did_fit}.
#' @export
plot_methods <- function(fit) {
  stopifnot(inherits(fit, "anchored_did_fit"))
  t <- fit$table
  t$method <- factor(t$method, levels = rev(t$method))
  ggplot(t, aes(x = .pp(.data$estimate), y = .data$method)) +
    geom_vline(xintercept = 0, linetype = 2, colour = "grey60") +
    geom_errorbar(aes(xmin = .pp(.data$lower), xmax = .pp(.data$upper)),
                  orientation = "y", width = 0.18, linewidth = 0.7) +
    geom_point(size = 3, colour = "#0072B2") +
    labs(x = expression(Delta[DiD] ~ "(percentage points)"), y = NULL,
         title = sprintf("%s vs %s, anchored on %s",
                         fit$input$labels[2], fit$input$labels[1], fit$input$labels[3]),
         subtitle = sprintf("%d%% intervals by method", round(100 * fit$level))) +
    .ssdta_theme()
}

#' The two component contrasts and the resulting DiD
#'
#' @param x a did_input object.
#' @param level confidence level for the component / DiD Wald bars.
#' @export
plot_components <- function(x, level = 0.95) {
  stopifnot(inherits(x, "did_input")); lab <- x$labels
  z <- stats::qnorm(1 - (1 - level) / 2)
  se_did <- if (!is.na(x$se_A) && !is.na(x$se_B)) sqrt(x$se_A^2 + x$se_B^2) else NA
  df <- data.frame(
    part = c(sprintf("Delta_A: Se(%s)-Se(%s)", lab[1], lab[3]),
             sprintf("Delta_B: Se(%s)-Se(%s)", lab[2], lab[3]),
             sprintf("Delta_DiD: %s vs %s", lab[2], lab[1])),
    est = c(x$delta_A, x$delta_B, x$did),
    se  = c(x$se_A, x$se_B, se_did),
    kind = c("component", "component", "difference-in-differences")
  )
  df$part <- factor(df$part, levels = rev(df$part))
  ggplot(df, aes(.pp(.data$est), .data$part, colour = .data$kind)) +
    geom_vline(xintercept = 0, linetype = 2, colour = "grey60") +
    geom_errorbar(aes(xmin = .pp(.data$est - z * .data$se),
                      xmax = .pp(.data$est + z * .data$se)),
                  orientation = "y", width = 0.16, linewidth = 0.7, na.rm = TRUE) +
    geom_point(size = 3, na.rm = TRUE) +
    scale_colour_manual(values = c("component" = "#565656",
                                   "difference-in-differences" = "#D55E00")) +
    labs(x = "Percentage points", y = NULL, colour = NULL,
         title = "From two anchored contrasts to the difference-in-differences") +
    .ssdta_theme()
}

#' Posterior density of Delta_DiD (Bayesian method)
#'
#' @param fit an \code{anchored_did_fit} that included the bayes method, or a
#'   \code{did_boot}/\code{did_bayes} result carrying a "draws" attribute.
#' @param level shaded credible level.
#' @export
plot_posterior <- function(fit, level = 0.95) {
  draws <- if (inherits(fit, "anchored_did_fit")) fit$draws$bayes else attr(fit, "draws")
  if (is.null(draws)) stop("No posterior/bootstrap draws found in this object.")
  a <- (1 - level) / 2
  qs <- stats::quantile(draws, c(a, 0.5, 1 - a), names = FALSE)
  d <- data.frame(x = .pp(draws))
  ggplot(d, aes(.data$x)) +
    geom_density(fill = "#D9ECF4", colour = "#0072B2", linewidth = 0.8) +
    geom_vline(xintercept = 0, linetype = 2, colour = "grey55") +
    geom_vline(xintercept = .pp(qs[2]), colour = "#D55E00", linewidth = 0.8) +
    annotate("segment", x = .pp(qs[1]), xend = .pp(qs[3]),
             y = 0, yend = 0, linewidth = 1.4, colour = "#0072B2") +
    labs(x = expression(Delta[DiD] ~ "(percentage points)"), y = "Posterior density",
         title = "Posterior of the anchored difference-in-differences",
         subtitle = sprintf("median %+.1f pp; %d%% CrI (%+.1f, %+.1f); P(<0) = %.3f",
                            .pp(qs[2]), round(100 * level), .pp(qs[1]), .pp(qs[3]),
                            mean(draws < 0))) +
    .ssdta_theme()
}

#' Overlay posterior densities across priors (sensitivity analysis)
#'
#' @param sens a \code{did_prior_sensitivity} object from
#'   \code{did_bayes_prior_sensitivity()}.
#' @export
plot_prior_sensitivity <- function(sens) {
  stopifnot(inherits(sens, "did_prior_sensitivity"))
  d <- do.call(rbind, lapply(names(sens$draws), function(nm)
    data.frame(prior = nm, x = .pp(sens$draws[[nm]]))))
  pal <- .plot_palette(unique(d$prior), known = character())
  ggplot(d, aes(.data$x, colour = .data$prior)) +
    geom_density(linewidth = 0.9) +
    geom_vline(xintercept = 0, linetype = 2, colour = "grey55") +
    scale_colour_manual(values = pal) +
    labs(x = expression(Delta[DiD] ~ "(percentage points)"), y = "Posterior density",
         colour = "Prior",
         title = "Prior-sensitivity analysis",
         subtitle = "Beta priors on each component sensitivity") +
    .ssdta_theme()
}
