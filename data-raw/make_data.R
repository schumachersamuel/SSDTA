# Build SSDTA package data objects (run from data-raw/).  Rscript data-raw/make_data.R
raw   <- Sys.getenv("RAW", "data-raw")
outd  <- Sys.getenv("OUTD", "data")
relab <- function(x) ifelse(x == "Ultra-FN", "Not-detected", as.character(x))

# ---- 1. NPOC stratum-specific posterior draws (the fitted S_k, 6 strata x 8000)
pd <- readRDS(file.path(raw, "NPOC_sensitivity_posterior_draws.rds"))
pd$category <- factor(relab(pd$category),
                      levels = c("High","Medium","Low","Very low","Trace","Not-detected"))
npoc_posterior_draws <- pd[, c("category","draw","sens_pooled","tau","sens_predicted")]

# ---- 2. Stratum summary (posterior medians + CrIs, from the paper's Stan fit)
sm <- readRDS(file.path(raw, "NPOC_sensitivity_summary.rds"))
sm$Category <- factor(relab(sm$Category),
                      levels = c("High","Medium","Low","Very low","Trace","Not-detected"))
npoc_stratum_summary <- sm

# ---- 3. Burden (bacillary-load) distribution profiles ----------------------
# Each: counts among culture-positive cases in the 5 Ultra-detected strata.
# The Not-detected (Ultra false-negative) stratum weight is handled separately
# via pi (see ssdta_standardize), because it is not observable in routine data.
mkprof <- function(setting, country, source, n, High, Medium, Low, VeryLow, Trace,
                   count_notdetected = NA, note = NA) {
  structure(list(
    setting = setting, country = country, source = source, n = n, note = note,
    counts_detected = c(High = High, Medium = Medium, Low = Low,
                        `Very low` = VeryLow, Trace = Trace),
    count_notdetected = count_notdetected),
    class = "burden_profile")
}

burden_profiles <- list(
  research = mkprof(
    "Research", NA,
    "Pooled Ultra-positive culture-confirmed cases, 5 NPOC-NAAT DTA studies (GDG PICO1)",
    507, 212, 108, 128, 35, 24, count_notdetected = 40,
    note = "40 Ultra-negative culture-positive cases observed; use pi='observed'"),
  routine_pakistan = mkprof(
    "Routine diagnostic", "Pakistan",
    "GxAlert National TB Reference Laboratory reports 2019 + 2021 (pooled)",
    88560, 28039, 15868, 22277, 9544, 12832,
    note = "Largest published programmatic 5-category breakdown; pi unknown -> sensitivity grid"),
  screening_kendall = mkprof(
    "Screening", "Uganda",
    "Kendall et al. community screening (observed Ultra semi-quant counts)",
    89, 14, 11, 21, 20, 23,
    note = "pi unknown -> sensitivity grid"),
  screening_pakistan = mkprof(
    "Screening", "Pakistan",
    "Pluslife/MiniDock evaluation, screening context (slide 1, 'Data for Samuel.pptx'); country per S. Schumacher",
    118, 12, 9, 34, 19, 44,
    note = "Second screening Sq distribution; poolable with Kendall (both screening)")
)

# ---- 4. Real per-study NPOC-NAAT counts by Ultra semi-quant stratum --------
# From Horne 2025 PICO1 Table 3 ("Sensitivity of NPOC-NAAT by Xpert Ultra
# semiquantitative status"), 5 studies. TP = NPOC-positive among culture-positives
# in that Ultra stratum; FN = NPOC-negative among culture-positives. The
# "Not-detected" stratum = culture-positive, Ultra-negative specimens.
# Verified: positive-stratum totals sum to 507 (= research burden counts).
npoc_by_ultra_stratum <- do.call(rbind, list(
  data.frame(study="Hartati 2025",  stratum="High",         tp=39, fn=2),
  data.frame(study="Mbuli 2025",    stratum="High",         tp=56, fn=4),
  data.frame(study="Rockman 2025",  stratum="High",         tp=1,  fn=1),
  data.frame(study="Steadman 2025", stratum="High",         tp=24, fn=0),
  data.frame(study="Yerlikaya 2025",stratum="High",         tp=85, fn=0),
  data.frame(study="Hartati 2025",  stratum="Medium",       tp=13, fn=4),
  data.frame(study="Mbuli 2025",    stratum="Medium",       tp=22, fn=0),
  data.frame(study="Rockman 2025",  stratum="Medium",       tp=1,  fn=0),
  data.frame(study="Steadman 2025", stratum="Medium",       tp=17, fn=0),
  data.frame(study="Yerlikaya 2025",stratum="Medium",       tp=50, fn=1),
  data.frame(study="Hartati 2025",  stratum="Low",          tp=28, fn=6),
  data.frame(study="Mbuli 2025",    stratum="Low",          tp=26, fn=2),
  data.frame(study="Rockman 2025",  stratum="Low",          tp=2,  fn=0),
  data.frame(study="Steadman 2025", stratum="Low",          tp=18, fn=1),
  data.frame(study="Yerlikaya 2025",stratum="Low",          tp=40, fn=5),
  data.frame(study="Hartati 2025",  stratum="Very low",     tp=0,  fn=8),
  data.frame(study="Mbuli 2025",    stratum="Very low",     tp=5,  fn=2),
  data.frame(study="Rockman 2025",  stratum="Very low",     tp=1,  fn=2),
  data.frame(study="Steadman 2025", stratum="Very low",     tp=2,  fn=1),
  data.frame(study="Yerlikaya 2025",stratum="Very low",     tp=12, fn=2),
  data.frame(study="Hartati 2025",  stratum="Trace",        tp=2,  fn=6),
  data.frame(study="Mbuli 2025",    stratum="Trace",        tp=4,  fn=5),
  data.frame(study="Rockman 2025",  stratum="Trace",        tp=0,  fn=0),
  data.frame(study="Steadman 2025", stratum="Trace",        tp=1,  fn=1),
  data.frame(study="Yerlikaya 2025",stratum="Trace",        tp=2,  fn=3),
  data.frame(study="Hartati 2025",  stratum="Not-detected", tp=3,  fn=3),
  data.frame(study="Mbuli 2025",    stratum="Not-detected", tp=1,  fn=5),
  data.frame(study="Rockman 2025",  stratum="Not-detected", tp=0,  fn=1),
  data.frame(study="Steadman 2025", stratum="Not-detected", tp=0,  fn=4),
  data.frame(study="Yerlikaya 2025",stratum="Not-detected", tp=2,  fn=21)
))
npoc_by_ultra_stratum$stratum <- factor(
  npoc_by_ultra_stratum$stratum,
  levels = c("High","Medium","Low","Very low","Trace","Not-detected"))
rownames(npoc_by_ultra_stratum) <- NULL

# ---- 5. Study characteristics for the 5 studies in npoc_by_ultra_stratum ---
# From Horne 2025 PICO1 Table 1 ("Characteristics of included studies") and the
# narrative in "Characteristics of included studies" / "Respiratory specimen
# testing". `n_enrolled` is the Table 1 study size; `n_culture_pos` and
# `n_ultra_pos` are derived from Table 3 below, not hand-entered, so the study
# table cannot drift away from the counts.
#
# Publication status verified 2026-08-17 against Crossref. At the time the GDG
# review was written (Nov 2025) four of the five were unpublished; three have
# since preprinted and one (Steadman) has appeared in eBioMedicine. See
# METHOD.md section 8 -- the release gate keys off this column.
npoc_studies <- data.frame(
  study          = c("Hartati 2025", "Mbuli 2025", "Rockman 2025",
                     "Steadman 2025", "Yerlikaya 2025"),
  first_author   = c("Hartati", "Mbuli", "Rockman", "Steadman", "Yerlikaya"),
  year           = c(2025L, 2025L, 2025L, 2025L, 2025L),
  cohort         = c("EVIDENT", NA, NA, "R2D2", "R2D2"),
  countries      = c("Indonesia",
                     "Cameroon",
                     "South Africa",
                     "Uganda; Vietnam",
                     "Nigeria; Philippines; South Africa; Uganda; Zambia; India; Vietnam"),
  case_finding   = c("Passive",
                     "Passive + active (32.6% via community ACF)",
                     "Passive",
                     "Passive/intensified/active",
                     "Passive/intensified/active"),
  population     = c(">28 days old, cough >=2 wks + >=1 other symptom; 13.2% aged <10 yrs",
                     ">=15 yrs, symptoms or clinical risk factors; ACF arm screened by CXR-AI (qXR v4)",
                     ">=18 yrs, presumptive TB, unable to produce sputum ('sputum scarce')",
                     ">=12 yrs, outpatient presumptive TB (cough >=2 wks, or risk factor + abnormal CXR/CRP)",
                     ">=12 yrs, outpatient presumptive TB"),
  index_test     = rep("Pluslife MiniDock MTB", 5),
  comparator     = rep("Xpert MTB/RIF Ultra", 5),
  specimen       = c("Expectorated sputum (100%)",
                     "Expectorated sputum (100%)",
                     "Induced sputum (100%)",
                     "Expectorated (90%) / induced (10%) sputum",
                     "Expectorated (89.8%) / induced (10.2%) sputum"),
  specimen_cond  = c("Fresh",
                     "Fresh; and swabs stored at 2-8 C, tested <24 h",
                     "Frozen",
                     "Fresh",
                     "Fresh"),
  testing_site   = c("Research laboratory",
                     "Facility: Xpert & reference labs; ACF: on-site testing",
                     "Research laboratory",
                     "Research laboratory (Uganda) / hospital laboratory (Vietnam)",
                     "Intermediate laboratory"),
  ref_standard   = c("Liquid culture x1",
                     "Solid & liquid culture x1",
                     "Liquid culture x1",
                     "Liquid culture x2",
                     "Liquid culture x2"),
  n_enrolled     = c(651L, 996L, 80L, 322L, 1300L),
  hiv_pct        = c(0.5, 18.8, 35.0, 20.0, 18.6),
  smear_neg_pct  = c(NA, 32.0, NA, 33.3, 39.0),
  invalid_pct    = c(0.2, 0.2, 0.0, 1.2, 1.1),
  pub_status     = c("Preprint", "Preprint", "Preprint", "Published", "Published"),
  citation       = c(
    "Hartati S, Koesoemadinata RC, Sharples K, McAllister S, et al. Assessment of Minidock MTB for the diagnosis of tuberculosis from sputum in patients presenting to health facilities in Indonesia. medRxiv 2026-04-30 (preprint).",
    "Mbuli C, Jean Bosco TF, Nsamenang R, Nestor B, et al. Diagnostic performance of the Pluslife MiniDock MTB and Molbio MTB Ultima assays to detect tuberculosis from tongue and sputum swabs among outpatients and in active case finding in Cameroon. medRxiv 2025-09-12 (preprint).",
    "Rockman L, Abdulgader S, Minnies S, Naidoo T, et al. Diagnostic accuracy of Pluslife MiniDock MTB on tongue swabs in sputum scarce people with presumptive TB: a retrospective analysis. Research Square 2026-03-23 (preprint).",
    "Steadman A, Kumar KM, Asege L, et al. Diagnostic accuracy of swab-based molecular tests for tuberculosis using near-point-of-care platforms: a multi-country evaluation. eBioMedicine 2025;121:105991.",
    "Yerlikaya S, Chirwa M, Ajide B, Castro MM, et al. Pulmonary Tuberculosis Detection with MiniDock MTB Using Swab Samples. N Engl J Med 2026;394:1710-1722."),
  doi            = c("10.64898/2026.04.22.26351245",
                     "10.1101/2025.09.11.25335435",
                     "10.21203/rs.3.rs-6225528/v1",
                     "10.1016/j.ebiom.2025.105991",
                     "10.1056/NEJMoa2509761"),
  stringsAsFactors = FALSE
)

# Derive the contributed-n columns from the count table so the two objects
# cannot disagree. n_culture_pos = all 6 strata (547); n_ultra_pos = the 5
# Ultra-detected strata (507, = burden_profiles$research counts).
.tot <- function(d) {
  v <- tapply(d$tp + d$fn, factor(d$study, levels = npoc_studies$study), sum)
  as.integer(ifelse(is.na(v), 0L, v))
}
npoc_studies$n_culture_pos <- .tot(npoc_by_ultra_stratum)
npoc_studies$n_ultra_pos   <- .tot(subset(npoc_by_ultra_stratum, stratum != "Not-detected"))
npoc_studies <- npoc_studies[, c(
  "study", "first_author", "year", "cohort", "countries", "case_finding",
  "population", "index_test", "comparator", "specimen", "specimen_cond",
  "testing_site", "ref_standard", "n_enrolled", "n_culture_pos", "n_ultra_pos",
  "hiv_pct", "smear_neg_pct", "invalid_pct", "pub_status", "citation", "doi")]
rownames(npoc_studies) <- NULL
stopifnot(
  identical(sort(npoc_studies$study), sort(unique(as.character(npoc_by_ultra_stratum$study)))),
  sum(npoc_studies$n_ultra_pos)   == 507L,
  sum(npoc_studies$n_culture_pos) == 547L)

dir.create(outd, showWarnings = FALSE)
save(npoc_posterior_draws,  file = file.path(outd, "npoc_posterior_draws.rda"),  compress = "xz")
save(npoc_stratum_summary,  file = file.path(outd, "npoc_stratum_summary.rda"),  compress = "xz")
save(burden_profiles,       file = file.path(outd, "burden_profiles.rda"),       compress = "xz")
save(npoc_by_ultra_stratum, file = file.path(outd, "npoc_by_ultra_stratum.rda"), compress = "xz")
save(npoc_studies,          file = file.path(outd, "npoc_studies.rda"),          compress = "xz")

cat("data objects written to", outd, "\n")
cat("posterior draws:", nrow(npoc_posterior_draws), "rows,",
    length(unique(npoc_posterior_draws$category)), "strata\n")
cat("burden profiles:", paste(names(burden_profiles), collapse=", "), "\n")
cat("per-study rows:", nrow(npoc_by_ultra_stratum),
    "| positive-stratum TP+FN total:",
    sum(subset(npoc_by_ultra_stratum, stratum != "Not-detected")[,c("tp","fn")]), "\n")
cat("studies:", nrow(npoc_studies), "| enrolled:", sum(npoc_studies$n_enrolled),
    "| culture-pos with Ultra grade:", sum(npoc_studies$n_culture_pos),
    "| Ultra-positive:", sum(npoc_studies$n_ultra_pos),
    "| published:", sum(npoc_studies$pub_status == "Published"), "of", nrow(npoc_studies), "\n")
