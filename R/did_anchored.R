# ---------------------------------------------------------------------------
# Unified entry point: run one or several interval methods on a did_input.
# ---------------------------------------------------------------------------

#' Estimate the anchored difference-in-differences with chosen method(s)
#'
#' @param x a \code{did_input} object.
#' @param methods character vector, any of "wald", "meta", "bootstrap", "bayes".
#'   Methods needing data the input lacks are skipped with a message.
#' @param level confidence / credible level.
#' @param ... passed to the individual method functions (e.g. B, draws, prior,
#'   tau2_method, hakn, seed).
#' @return object of class \code{anchored_did_fit}: a table of estimates by
#'   method plus the underlying input and any stored draws.
#' @export
anchored_did <- function(x, methods = c("wald", "meta", "bootstrap", "bayes"),
                         level = 0.95, ...) {
  stopifnot(inherits(x, "did_input"))
  methods <- match.arg(methods, several.ok = TRUE)
  dots <- list(...)
  call1 <- function(f, ...) {
    args <- c(list(x = x, level = level), list(...))
    keep <- args[names(args) %in% names(formals(f)) | names(args) == "x"]
    do.call(f, keep)
  }
  rows <- list(); store <- list()
  for (m in methods) {
    res <- tryCatch({
      switch(m,
        wald      = do.call(did_wald,  c(list(x = x, level = level))),
        meta      = do.call(did_meta,  c(list(x = x, level = level),
                             dots[intersect(names(dots), c("tau2_method", "hakn"))])),
        bootstrap = do.call(did_boot,  c(list(x = x, level = level),
                             dots[intersect(names(dots), c("B", "rho", "seed"))])),
        bayes     = do.call(did_bayes, c(list(x = x, level = level),
                             dots[intersect(names(dots), c("prior", "draws", "rho", "seed"))]))
      )
    }, error = function(e) {
      message(sprintf("  [skipped %s] %s", m, conditionMessage(e))); NULL
    })
    if (!is.null(res)) {
      rows[[m]] <- res[, c("method", "estimate", "se", "lower", "upper")]
      if (!is.null(attr(res, "draws"))) store[[m]] <- attr(res, "draws")
    }
  }
  if (!length(rows)) stop("No method could be run on this input.")
  tab <- do.call(rbind, rows); rownames(tab) <- NULL
  structure(list(table = tab, input = x, level = level, draws = store),
            class = "anchored_did_fit")
}

#' @export
print.anchored_did_fit <- function(x, ...) {
  lab <- x$input$labels
  cat(sprintf("Anchored difference-in-differences: Se(%s) - Se(%s), anchored on %s\n",
              lab[2], lab[1], lab[3]))
  cat(sprintf("Point estimate: %+.1f pp   (%d%% intervals below)\n\n",
              100 * x$input$did, round(100 * x$level)))
  t <- x$table
  disp <- data.frame(
    Method = t$method,
    Est_pp = sprintf("%+.1f", 100 * t$estimate),
    SE_pp  = ifelse(is.na(t$se), "-", sprintf("%.1f", 100 * t$se)),
    CI_pp  = sprintf("(%+.1f, %+.1f)", 100 * t$lower, 100 * t$upper),
    stringsAsFactors = FALSE
  )
  print(disp, row.names = FALSE)
  invisible(x)
}

#' @export
summary.anchored_did_fit <- function(object, ...) {
  print(object)
  if (!is.null(object$draws$bayes)) {
    cat(sprintf("\nBayesian posterior P(Delta_DiD < 0) = %.3f\n",
                mean(object$draws$bayes < 0)))
  }
  invisible(object$table)
}
