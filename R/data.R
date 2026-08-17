# ---------------------------------------------------------------------------
# Documentation for the shipped datasets. All from the WHO TB Dx GDG November
# 2025 reviews and Samuel Schumacher's spectrum-standardized DTA project.
# ---------------------------------------------------------------------------

#' NPOC-NAAT stratum-specific sensitivity posterior draws
#'
#' Posterior draws of NPOC-NAAT (Pluslife MiniDock) sensitivity within each
#' Xpert Ultra semi-quant stratum, from the paper's Bayesian binomial-logit
#' random-effects meta-analysis (brms/Stan; Normal(0,1.5) intercept,
#' Exponential(1) between-study SD). 8000 draws x 6 strata. This is the
#' authoritative S_k input for \code{\link{ssdta_standardize}}.
#'
#' @format data.frame, 48000 rows: \code{category} (factor, 6 strata),
#'   \code{draw}, \code{sens_pooled} (posterior sensitivity), \code{tau}
#'   (between-study SD draw), \code{sens_predicted} (predictive, new study).
#' @source Schumacher et al., Union 2026 abstract UNIONCONF2026:3197; GDG PICO1.
"npoc_posterior_draws"

#' NPOC-NAAT stratum-specific sensitivity summary
#'
#' Posterior medians, 95\%/80\% credible intervals, predictive intervals and
#' convergence diagnostics per stratum (the summary of
#' \code{npoc_posterior_draws}).
#' @format data.frame, 6 rows (one per stratum).
#' @source As \code{npoc_posterior_draws}.
"npoc_stratum_summary"

#' Bacillary-burden distribution profiles (case-mix weights)
#'
#' A named list of \code{burden_profile} objects giving the Xpert Ultra
#' semi-quant distribution of true (culture-positive) cases in different target
#' populations: \code{research} (GDG PICO1 studies, with observed Not-detected
#' count), \code{routine_pakistan} (GxAlert National reports 2019+2021,
#' n=88,560), \code{screening_kendall} (Kendall community screening, Uganda,
#' n=89), \code{screening_pakistan} (Pluslife/MiniDock screening evaluation,
#' n=118). Used as the target distribution in \code{\link{ssdta_standardize}}.
#' @format named list of \code{burden_profile} objects.
#' @source GDG PICO1/PICO4 reports; GxAlert Pakistan national reports; Kendall
#'   et al.; 'Data for Samuel.pptx' slide 1.
"burden_profiles"

#' NPOC-NAAT counts by Xpert Ultra semi-quant stratum, per study
#'
#' True-positive / false-negative counts (against culture) of NPOC-NAAT within
#' each Ultra semi-quant stratum, for the 5 GDG PICO1 studies. Real data behind
#' the S_k meta-analysis; use with \code{\link{ssdta_meta}} to re-fit. Positive
#' strata sum to 507 (matches the research burden counts).
#' @format data.frame, 30 rows: \code{study}, \code{stratum} (factor, 6 levels),
#'   \code{tp}, \code{fn}.
#' @source Horne et al. 2025, GDG PICO1 systematic review, Table 3.
"npoc_by_ultra_stratum"

#' Characteristics of the 5 studies behind \code{npoc_by_ultra_stratum}
#'
#' One row per study: design, setting, specimen handling, reference standard,
#' contributed sample sizes, and publication status. All five evaluated the same
#' index test (Pluslife MiniDock MTB) against culture, with Xpert Ultra as the
#' comparator that defines the semi-quant stratifier. \code{study} matches
#' \code{npoc_by_ultra_stratum$study} exactly, so the two join directly.
#'
#' \code{n_culture_pos} and \code{n_ultra_pos} are derived from
#' \code{npoc_by_ultra_stratum} rather than transcribed (547 and 507 in total),
#' so the two objects cannot drift apart. \code{n_enrolled} is the Table 1 study
#' size and is larger: it counts all participants, not just culture-positives
#' with an Ultra grade.
#'
#' \code{pub_status} drives the package release gate (see METHOD.md section 8):
#' four of the five were unpublished when the GDG review was written, so the
#' per-study counts were internal. Status was re-verified against Crossref on
#' 2026-08-17 -- Steadman (eBioMedicine) and Yerlikaya (NEJM) are now published
#' and the other three have preprinted. The gate decision itself has not been
#' changed here; check with the maintainer before sharing externally.
#'
#' @format data.frame, 5 rows x 22 columns. Identity and design:
#'   \code{study}, \code{first_author}, \code{year}, \code{cohort} (EVIDENT /
#'   R2D2 / NA), \code{countries}, \code{case_finding}, \code{population}.
#'   Tests and specimens: \code{index_test}, \code{comparator},
#'   \code{specimen}, \code{specimen_cond} (fresh / frozen / stored),
#'   \code{testing_site}, \code{ref_standard}. Sizes: \code{n_enrolled},
#'   \code{n_culture_pos} (all 6 strata), \code{n_ultra_pos} (5 Ultra-detected
#'   strata), \code{hiv_pct}, \code{smear_neg_pct} (among culture-positives;
#'   \code{NA} = not reported), \code{invalid_pct}. Provenance:
#'   \code{pub_status} ("Published" / "Preprint"), \code{citation}, \code{doi}.
#' @seealso \code{\link{npoc_by_ultra_stratum}} for the counts themselves.
#' @source Horne et al. 2025, GDG PICO1 systematic review, Table 1 and the
#'   "Characteristics of included studies" narrative; sizes derived from
#'   Table 3. Citations and DOIs verified against Crossref, 2026-08-17.
"npoc_studies"
