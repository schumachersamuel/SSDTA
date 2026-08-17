# ---------------------------------------------------------------------------
# Real study-level data: Pooled Ultra (A) vs NPOC/MiniDock (B), anchored on
# individual Xpert Ultra (C), from the November 2025 WHO TB Dx GDG reviews.
#
#   Side A -- PICO 4, Horne pooled review (Start4All), Figure 7 (1:4 ratio):
#       7 fully-paired sub-studies, pooled Ultra vs individual Ultra, same
#       specimens. TP/FP/FN/TN read directly off the published forest plot.
#   Side B -- PICO 1, Horne respiratory review, Figure 6:
#       5 paired studies, MiniDock (NPOC) vs Xpert Ultra, same specimens.
#       TP/FP/FN/TN read directly off the published forest plot.
#
# Aggregated sensitivities reproduce the reports' summary figures:
#   Pooled Ultra 84.5% (report: 84.8%), Individual Ultra 87.9% (report: 87.9%),
#   NPOC 85.6% (report: 85.6%), Ultra head-to-head 92.7% raw-sum vs 93.1%
#   random-effects (small, expected RE-vs-raw-sum gap).
#   => Delta_A = -3.1 pp (report), Delta_B = -7.5 pp (report), Delta_DiD = -4.4 pp.
#
# Source: Horne et al. 2025, "Near point-of-care nucleic acid amplification
# tests (NPOC-NAATs) using respiratory specimens for detection of pulmonary
# tuberculosis" (PICO 1) and "LC-aNAAT on pooled versus individual sputum"
# (PICO 4), WHO TB Dx GDG systematic reviews, November 2025 (pre-publication).
# 01_Sources/systematic-reviews/Horne-2025-NPOC-NAAT-respiratory-SR-GDG-PICO1.docx
# 01_Sources/systematic-reviews/Horne-2025-LC-aNAAT-pooled-respiratory-SR-GDG-PICO4.docx
# ---------------------------------------------------------------------------

# Internal study-level data for the DiD fallback. These are NOT shipped as
# lazy-loaded datasets -- they are package-internal objects reached through
# pooling_vs_npoc_example(), which is the exported entry point.
#
# .pico1_npoc_vs_ultra
#   Five studies, same specimens tested by both index tests against culture
#   (PICO 1, Figure 6). Columns: TP/FP/FN/TN for NPOC (MiniDock) and Ultra.
#   Source: WHO TB Dx GDG PICO 1 (Horne NPOC respiratory review), Nov 2025.
#
# .pico4_pooled_vs_individual
#   Seven Start4All sub-studies, same specimens pooled (1:4) and tested
#   individually against culture (PICO 4, Figure 7). Columns: TP/FP/FN/TN
#   for each arm.
#   Source: WHO TB Dx GDG PICO 4 (Horne pooled respiratory review), Nov 2025.

.pico1_npoc_vs_ultra <- data.frame(
  study    = c("Hartati 2025", "Mbuli 2025", "Steadman 2025", "Yerlikaya 2025", "Rockman 2025"),
  npoc_tp  = c(98, 114, 62, 191, 5),  npoc_fp  = c(12, 20, 4, 25, 2),
  npoc_fn  = c(18, 18, 7, 32, 4),     npoc_tn  = c(442, 844, 217, 1031, 53),
  ultra_tp = c(110, 126, 65, 200, 8), ultra_fp = c(27, 24, 3, 18, 2),
  ultra_fn = c(6, 6, 4, 23, 1),       ultra_tn = c(427, 840, 218, 1038, 53),
  stringsAsFactors = FALSE
)

.pico4_pooled_vs_individual <- data.frame(
  study         = c("Start4All Malawi", "Start4All Brazil", "Start4All Bangladesh",
                    "Start4All Kenya", "Start4All Cameroon", "Start4All Vietnam",
                    "Start4All Nigeria"),
  pooled_tp     = c(12, 12, 148, 37, 181, 287, 86),
  pooled_fp     = c(7, 6, 31, 23, 23, 33, 70),
  pooled_fn     = c(2, 0, 24, 4, 25, 60, 25),
  pooled_tn     = c(246, 192, 2744, 614, 2317, 1904, 3174),
  individual_tp = c(12, 12, 153, 38, 186, 301, 92),
  individual_fp = c(7, 8, 51, 25, 28, 45, 96),
  individual_fn = c(2, 0, 19, 3, 20, 46, 19),
  individual_tn = c(246, 190, 2724, 612, 2312, 1892, 3148),
  stringsAsFactors = FALSE
)

#' Build the worked-example anchored-DiD input from the real GDG study data
#'
#' Pooled Ultra (A) vs NPOC (B), anchored on individual Xpert Ultra (C).
#'
#' \strong{Three input routes are exposed; on this data they all give
#' Delta_DiD approx -4.4 pp:}
#' \itemize{
#'   \item \code{contrasts} uses the GDG reports' own random-effects paired
#'     contrasts and 95\% CIs (the numbers in the evidence tables):
#'     Delta_A = -3.10pp (-6.95, +0.76), Delta_B = -7.54pp (-11.6, -3.47),
#'     Delta_DiD = -4.44pp (-10.0, +1.2). This is the headline route.
#'   \item \code{components} uses the SAME reported random-effects summary
#'     sensitivities and their 95\% CIs (NPOC 85.6\% [82.4-88.3], Ultra 93.1\%
#'     [89.7-95.5], pooled 84.8\% [81.3-87.8], individual 87.9\% [85.6-89.9]),
#'     reconstructing an effective n per sensitivity from each CI. This
#'     reproduces Delta_DiD = -4.4pp and, at the default rho = 0, a Wald
#'     interval essentially identical to the contrast form -- because the
#'     reported contrast CIs themselves imply a near-zero marginal correlation
#'     (see the rho note in \code{?did_input}).
#' }
#' The three routes agree closely: \code{contrasts} -> \code{did_wald} gives
#' -4.44pp (-10.0, +1.2); \code{components} -> \code{did_boot}/\code{did_bayes}
#' reproduces it with boundary-robust intervals; \code{per_study} ->
#' \code{did_meta} random-effects pools the real per-study contrasts and lands
#' near -4.4 to -4.5pp (it depends on \code{rho_patient} -- see below -- and is
#' approximate because the primary studies' joint concordance tables are not in
#' the reports).
#'
#' \strong{Two different correlations, at two levels -- do not conflate:}
#' \itemize{
#'   \item \emph{Patient-level} (within a study, same people tested twice):
#'     moderate-to-high; it lowers each study's per-study contrast SE and so
#'     feeds \code{per_study} / \code{did_meta}. Controlled by \code{rho_patient}
#'     (default 0.5).
#'   \item \emph{Marginal-level} (between the two RE-pooled sensitivities you
#'     difference): ~0.02-0.03 here (see \code{implied_rho}), because
#'     between-study heterogeneity swamps the pairing; it feeds the
#'     \code{components} route via \code{did_input(rho=)} (default 0).
#' }
#' (An earlier version derived the component sensitivities by simple summation
#' of the raw per-study counts instead of using the reported RE figures; that
#' gave a spurious -3.7pp because raw-sum pooling != RE pooling, most visibly for
#' the Ultra anchor -- 92.7\% raw-sum vs 93.1\% RE. Fixed.)
#'
#' @param rho_patient PATIENT-LEVEL within-study correlation used to compute the
#'   per-study contrast SEs in \code{per_study} (for \code{did_meta}). This is a
#'   DIFFERENT quantity from the marginal-level rho of the component form: at the
#'   patient level the same people are tested by both the index test and the
#'   anchor, so the correlation is genuinely moderate-to-high (default 0.5).
#'   Because the primary studies' joint (index x anchor) concordance tables are
#'   not in the GDG reports, this is an approximation -- the exact per-study
#'   paired SEs would come from those tables. The \code{components} route is
#'   independent of this argument (it uses the reported pooled marginals, whose
#'   own implied marginal correlation is ~0; see \code{implied_rho}).
#' @return list(contrasts, components, per_study, raw).
#' @export
pooling_vs_npoc_example <- function(rho_patient = 0.5) {
  p1 <- .pico1_npoc_vs_ultra
  p4 <- .pico4_pooled_vs_individual
  z  <- stats::qnorm(0.975)
  se_ci <- function(lo, hi) (hi - lo) / (2 * z)

  # Reported random-effects summary sensitivities (proportions) and 95% CIs
  sA  <- 0.848; ciA_s  <- c(0.813, 0.878)   # Pooled Ultra          (PICO4)
  sCa <- 0.879; ciCa_s <- c(0.856, 0.899)   # Individual Ultra      (PICO4)
  sB  <- 0.856; ciB_s  <- c(0.824, 0.883)   # NPOC / MiniDock       (PICO1)
  sCb <- 0.931; ciCb_s <- c(0.897, 0.955)   # Xpert Ultra (anchor)  (PICO1)

  # real per-study paired contrasts (for did_meta); patient-level rho
  ct_A <- paired_contrasts(p4$pooled_tp, p4$pooled_fn, p4$individual_tp, p4$individual_fn, rho = rho_patient)
  ct_B <- paired_contrasts(p1$npoc_tp,   p1$npoc_fn,   p1$ultra_tp,      p1$ultra_fn,      rho = rho_patient)

  list(
    contrasts = list(
      # exactly as reported in the SRs (RE-pooled paired contrasts)
      delta_A = -0.0310, ci_A = c(-0.0695,  0.0076),  # pooled vs individual Ultra
      delta_B = -0.0754, ci_B = c(-0.1160, -0.0347),  # NPOC vs Ultra
      labels  = c("Pooled Ultra", "NPOC", "Individual Ultra")
    ),
    components = list(
      # reported RE summary sensitivities + SEs from their CIs (NOT raw sums)
      sens = list(A = sA, Ca = sCa, B = sB, Cb = sCb,
                  se = c(A  = se_ci(ciA_s[1],  ciA_s[2]),
                         Ca = se_ci(ciCa_s[1], ciCa_s[2]),
                         B  = se_ci(ciB_s[1],  ciB_s[2]),
                         Cb = se_ci(ciCb_s[1], ciCb_s[2]))),
      labels = c("Pooled Ultra", "NPOC", "Individual Ultra")
    ),
    per_study = list(
      delta_A_i = ct_A$delta, se_A_i = ct_A$se, study_A = p4$study,
      delta_B_i = ct_B$delta, se_B_i = ct_B$se, study_B = p1$study
    ),
    raw = list(pico1 = p1, pico4 = p4)
  )
}

# retained for backward compatibility with earlier scripts
.make_pooling_npoc <- function(rho = 0.5) pooling_vs_npoc_example(rho_patient = rho)
