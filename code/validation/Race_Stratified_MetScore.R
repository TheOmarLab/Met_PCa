################################################################################
# Race_Stratified_MetScore.R
#
# Durham race-stratified Met-Score analysis. Race codes follow the clin_codebook
# sheet of Durham_cohort_clinical_data_022526.xlsx (row 6): 1=Black, 3=White,
# 10=Asian/Pacific Islander, 13=American Indian/Alaska Native. The frozen locked
# Met-Score probabilities and High/Low classes are used as provided; nothing is
# refit, rescaled for prediction, recalibrated, or rethresholded.
#
# Formal cause-specific survival inference is limited to Black and White patients,
# which together hold all 40 metastatic events (Black n=305/21, White n=238/19);
# the other race groups have no events and enter only the count ledger. A single
# robust cause-specific Cox model carries a Met-Score-by-race interaction, and the
# race-specific Met-Score hazard ratios per SD are read off as model contrasts.
# The frozen High-vs-Low interaction is kept only as a secondary sensitivity.
# Locked-score distributions are compared between Black and White patients.
################################################################################

suppressWarnings(suppressMessages(library(survival)))

ROOT <- local({ e <- Sys.getenv("MET_PCA_ROOT", ""); if (nzchar(e) && dir.exists(e)) e else normalizePath(".") })
OUT <- file.path(ROOT, "outs"); dir.create(OUT, showWarnings = FALSE, recursive = TRUE)

.warn <- new.env(); .warn$log <- list()
cap <- function(ctx, expr) withCallingHandlers(expr, warning = function(w) {
  .warn$log[[length(.warn$log) + 1L]] <- data.frame(context = ctx, message = conditionMessage(w),
                                                     stringsAsFactors = FALSE)
  invokeRestart("muffleWarning") })

# competing clock: metastasis = event; death without prior metastasis is censored
# at its own death time; analysis time is the metastasis time otherwise.
clock <- function(met, mt, dead, dt) {
  st <- ifelse(met == 1L, 1L, ifelse(dead == 1L & met == 0L & is.finite(dt) & dt <= mt, 2L, 0L))
  list(ev = as.integer(st == 1L), time = ifelse(st == 2L, dt, mt)) }
gg_ord <- function(gs, pr) { o <- rep(NA_real_, length(gs)); o[gs <= 6] <- 1; o[gs == 7 & pr == 3] <- 2
  o[gs == 7 & pr == 4] <- 3; o[gs == 8] <- 4; o[gs >= 9] <- 5; o }
pt_collapse <- function(x) { x <- toupper(trimws(as.character(x))); o <- rep(NA_character_, length(x))
  o[grepl("^T2", x)] <- "T2"; o[grepl("^T3", x)] <- "T3"; o[grepl("^T4", x)] <- "T4"; o }

RACE_MAP <- c("1" = "Black", "3" = "White", "10" = "Asian/Pacific Islander",
              "13" = "American Indian/Alaska Native")
ETH_MAP  <- c("1" = "Hispanic/Latino", "2" = "Not Hispanic/Latino")

## ---- assemble the Durham analysis frame ---------------------------------
e <- new.env(); load(file.path(ROOT, "output", "Durham", "durham_metscore_batchcorrected.rda"), envir = e)
D0 <- e$clin_valid
cl <- clock(as.integer(D0$mets), as.numeric(D0$surgmets), as.integer(D0$dead), as.numeric(D0$limbo))
D <- data.frame(
  time = cl$time, ev = cl$ev,
  ms_prob = as.numeric(D0$MetScore_prob),
  ms_z = as.numeric(scale(as.numeric(D0$MetScore_prob))),        # per cohort SD (association scale)
  cls = ifelse(as.character(D0$MetScoreClass) == "High risk", "High", "Low"),
  gg = gg_ord(as.numeric(as.character(D0$PathGleason)), as.numeric(as.character(D0$pogl1))),
  log2PSA = log2(as.numeric(D0$psapresurg) + 1),
  pT34 = { p <- pt_collapse(D0$stg); ifelse(p %in% c("T3", "T4"), 1L, ifelse(p == "T2", 0L, NA_integer_)) },
  race = unname(RACE_MAP[as.character(D0$race)]),
  eth = unname(ETH_MAP[as.character(D0$ethnicity)]),
  stringsAsFactors = FALSE)

## ---- race / ethnicity count + event ledger ------------------------------
race_levels <- c("Black", "White", "Asian/Pacific Islander", "American Indian/Alaska Native")
ledger <- do.call(rbind, lapply(race_levels, function(r) { s <- D[which(D$race == r), ]
  data.frame(group_type = "race", group = r, n = nrow(s), events = sum(s$ev),
             inference = ifelse(sum(s$ev) > 0, "formal cause-specific Cox", "counts only (no events)"),
             stringsAsFactors = FALSE) }))
n_unmapped <- sum(is.na(D$race))
if (n_unmapped > 0) ledger <- rbind(ledger, data.frame(group_type = "race", group = "Other/undocumented code",
  n = n_unmapped, events = sum(D$ev[is.na(D$race)]), inference = "counts only", stringsAsFactors = FALSE))
eth_levels <- c("Hispanic/Latino", "Not Hispanic/Latino")
ledger <- rbind(ledger, do.call(rbind, lapply(eth_levels, function(x) { s <- D[which(D$eth == x), ]
  data.frame(group_type = "ethnicity", group = x, n = nrow(s), events = sum(s$ev),
             inference = "counts only", stringsAsFactors = FALSE) })))
n_eth_na <- sum(is.na(D$eth))
if (n_eth_na > 0) ledger <- rbind(ledger, data.frame(group_type = "ethnicity", group = "Unknown/undocumented code",
  n = n_eth_na, events = sum(D$ev[is.na(D$eth)]), inference = "counts only", stringsAsFactors = FALSE))
write.csv(ledger, file.path(OUT, "Race_Stratified_MetScore_Durham_ledger.csv"), row.names = FALSE)

## ---- one robust cause-specific Cox with a Met-Score x race interaction ---
BW <- D[D$race %in% c("Black", "White") & is.finite(D$time) & D$time > 0 &
        is.finite(D$ms_z) & is.finite(D$gg) & is.finite(D$log2PSA) & !is.na(D$pT34), ]
BW$race <- relevel(factor(BW$race), ref = "White")
fit <- cap("interaction_cox", survival::coxph(Surv(time, ev) ~ ms_z * race + gg + log2PSA + pT34,
                                              data = BW, robust = TRUE))
b <- coef(fit); V <- vcov(fit)
im <- which(names(b) == "ms_z"); ix <- grep("ms_z:raceBlack", names(b))
white_b <- unname(b[im]); white_se <- unname(sqrt(V[im, im]))
black_b <- unname(b[im] + b[ix]); black_se <- unname(sqrt(V[im, im] + V[ix, ix] + 2 * V[im, ix]))
int_b <- unname(b[ix]); int_se <- unname(sqrt(V[ix, ix]))
converged <- all(is.finite(b)) && all(is.finite(diag(V)))
zph <- cap("cox_zph", survival::cox.zph(fit)); ph_global <- unname(zph$table["GLOBAL", "p"])
contrast <- rbind(
  data.frame(quantity = "White Met-Score HR per SD", estimate = round(exp(white_b), 5),
             ci_lo = round(exp(white_b - 1.96 * white_se), 5), ci_hi = round(exp(white_b + 1.96 * white_se), 5),
             p = signif(2 * pnorm(-abs(white_b / white_se)), 5), stringsAsFactors = FALSE),
  data.frame(quantity = "Black Met-Score HR per SD (contrast)", estimate = round(exp(black_b), 5),
             ci_lo = round(exp(black_b - 1.96 * black_se), 5), ci_hi = round(exp(black_b + 1.96 * black_se), 5),
             p = signif(2 * pnorm(-abs(black_b / black_se)), 5), stringsAsFactors = FALSE),
  data.frame(quantity = "Met-Score x race interaction ratio (Black vs White)", estimate = round(exp(int_b), 5),
             ci_lo = round(exp(int_b - 1.96 * int_se), 5), ci_hi = round(exp(int_b + 1.96 * int_se), 5),
             p = signif(2 * pnorm(-abs(int_b / int_se)), 5), stringsAsFactors = FALSE))
contrast$n_white <- sum(BW$race == "White"); contrast$n_black <- sum(BW$race == "Black")
contrast$events_white <- sum(BW$ev[BW$race == "White"]); contrast$events_black <- sum(BW$ev[BW$race == "Black"])
contrast$global_ph_p <- signif(ph_global, 5); contrast$converged <- converged
contrast$estimator <- "robust cause-specific Cox (Met-Score x race + ordinal GG + log2(PSA+1) + pT3/4; White ref.)"
write.csv(contrast, file.path(OUT, "Race_Stratified_MetScore_Durham_interaction.csv"), row.names = FALSE)

## ---- proportional-hazards diagnostics (per term + global) ----------------
ph <- as.data.frame(zph$table); ph$term <- rownames(zph$table)
ph <- ph[, c("term", "chisq", "df", "p")]; ph[c("chisq", "p")] <- lapply(ph[c("chisq", "p")], function(x) signif(x, 5))
write.csv(ph, file.path(OUT, "Race_Stratified_MetScore_Durham_PH.csv"), row.names = FALSE)

## ---- Grade Group-stratified sensitivity (same Black+White set) ----------
# The primary global PH is non-significant, but the Grade Group term departs from PH
# (p=0.0215); stratifying by Grade Group removes it from the hazard. This gives the
# same null-interaction conclusion.
fitg <- cap("gg_strat_cox", survival::coxph(Surv(time, ev) ~ ms_z * race + log2PSA + pT34 + strata(gg),
                                            data = BW, robust = TRUE))
bg <- coef(fitg); Vg <- vcov(fitg)
img <- which(names(bg) == "ms_z"); ixg <- grep("ms_z:raceBlack", names(bg))
gw_b <- unname(bg[img]); gw_se <- unname(sqrt(Vg[img, img]))
gb_b <- unname(bg[img] + bg[ixg]); gb_se <- unname(sqrt(Vg[img, img] + Vg[ixg, ixg] + 2 * Vg[img, ixg]))
gi_b <- unname(bg[ixg]); gi_se <- unname(sqrt(Vg[ixg, ixg]))
converged_g <- all(is.finite(bg)) && all(is.finite(diag(Vg)))
zphg <- cap("gg_strat_zph", survival::cox.zph(fitg)); ph_global_g <- unname(zphg$table["GLOBAL", "p"])
gg_row <- function(section, quantity, estimate = NA_real_, ci_lo = NA_real_, ci_hi = NA_real_,
                   p = NA_real_, chisq = NA_real_, df = NA_real_)
  data.frame(section = section, quantity = quantity, estimate = estimate, ci_lo = ci_lo, ci_hi = ci_hi,
             p = p, chisq = chisq, df = df, stringsAsFactors = FALSE)
gg_contrasts <- rbind(
  gg_row("contrast", "White Met-Score HR per SD", round(exp(gw_b), 5),
         round(exp(gw_b - 1.96 * gw_se), 5), round(exp(gw_b + 1.96 * gw_se), 5), signif(2 * pnorm(-abs(gw_b / gw_se)), 5)),
  gg_row("contrast", "Black Met-Score HR per SD (contrast)", round(exp(gb_b), 5),
         round(exp(gb_b - 1.96 * gb_se), 5), round(exp(gb_b + 1.96 * gb_se), 5), signif(2 * pnorm(-abs(gb_b / gb_se)), 5)),
  gg_row("contrast", "Met-Score x race interaction ratio (Black vs White)", round(exp(gi_b), 5),
         round(exp(gi_b - 1.96 * gi_se), 5), round(exp(gi_b + 1.96 * gi_se), 5), signif(2 * pnorm(-abs(gi_b / gi_se)), 5)))
gg_ph <- do.call(rbind, lapply(rownames(zphg$table), function(t)
  gg_row("schoenfeld", t, chisq = signif(unname(zphg$table[t, "chisq"]), 5),
         df = unname(zphg$table[t, "df"]), p = signif(unname(zphg$table[t, "p"]), 5))))
gg_out <- rbind(gg_contrasts, gg_ph)
gg_out$n_white <- sum(BW$race == "White"); gg_out$n_black <- sum(BW$race == "Black")
gg_out$events_white <- sum(BW$ev[BW$race == "White"]); gg_out$events_black <- sum(BW$ev[BW$race == "Black"])
gg_out$converged <- converged_g; gg_out$global_ph_p <- signif(ph_global_g, 5)
gg_out$model <- "robust cause-specific Cox: ms_z*race + log2(PSA+1) + pT3/4 + strata(Grade Group); White ref."
write.csv(gg_out, file.path(OUT, "Race_Stratified_MetScore_Durham_GG_stratified_sensitivity.csv"), row.names = FALSE)

## ---- secondary: frozen High-vs-Low class x race (unadjusted) ------------
BWb <- BW; BWb$grp <- relevel(factor(BWb$cls), ref = "Low")
fb <- cap("binary_cox", survival::coxph(Surv(time, ev) ~ grp * race, data = BWb, robust = TRUE))
bb <- coef(fb); Vb <- vcov(fb); jg <- which(names(bb) == "grpHigh"); jx <- grep("grpHigh:raceBlack", names(bb))
wb <- unname(bb[jg]); wse <- unname(sqrt(Vb[jg, jg]))
kb <- unname(bb[jg] + bb[jx]); kse <- unname(sqrt(Vb[jg, jg] + Vb[jx, jx] + 2 * Vb[jg, jx]))
ib <- unname(bb[jx]); ise <- unname(sqrt(Vb[jx, jx]))
binsens <- rbind(
  data.frame(quantity = "White High-vs-Low HR (frozen class, unadjusted)", estimate = round(exp(wb), 5),
             ci_lo = round(exp(wb - 1.96 * wse), 5), ci_hi = round(exp(wb + 1.96 * wse), 5),
             p = signif(2 * pnorm(-abs(wb / wse)), 5), stringsAsFactors = FALSE),
  data.frame(quantity = "Black High-vs-Low HR (frozen class, unadjusted, contrast)", estimate = round(exp(kb), 5),
             ci_lo = round(exp(kb - 1.96 * kse), 5), ci_hi = round(exp(kb + 1.96 * kse), 5),
             p = signif(2 * pnorm(-abs(kb / kse)), 5), stringsAsFactors = FALSE),
  data.frame(quantity = "High-vs-Low x race interaction ratio", estimate = round(exp(ib), 5),
             ci_lo = round(exp(ib - 1.96 * ise), 5), ci_hi = round(exp(ib + 1.96 * ise), 5),
             p = signif(2 * pnorm(-abs(ib / ise)), 5), stringsAsFactors = FALSE))
binsens$note <- "secondary sensitivity only; not evidence of racial heterogeneity"
write.csv(binsens, file.path(OUT, "Race_Stratified_MetScore_Durham_binary_sensitivity.csv"), row.names = FALSE)

## ---- locked-score distribution: Black vs White --------------------------
mb <- D$ms_prob[which(D$race == "Black")]; mw <- D$ms_prob[which(D$race == "White")]
mb <- mb[is.finite(mb)]; mw <- mw[is.finite(mw)]
wilx <- cap("wilcoxon", wilcox.test(mb, mw))
# adjusted difference in cohort-SD units via linear model with HC3 robust SE (z-based Wald)
LD <- D[D$race %in% c("Black", "White") & is.finite(D$ms_z) & is.finite(D$gg) &
        is.finite(D$log2PSA) & !is.na(D$pT34), ]
LD$race <- relevel(factor(LD$race), ref = "White")
lm1 <- lm(ms_z ~ race + gg + log2PSA + pT34, data = LD)
X <- model.matrix(lm1); hh <- hat(X); rr <- residuals(lm1); XtXi <- solve(crossprod(X))
Vhc3 <- XtXi %*% (t(X) %*% diag(rr^2 / (1 - hh)^2) %*% X) %*% XtXi
k <- which(colnames(X) == "raceBlack"); adj_est <- unname(coef(lm1)[k]); adj_se <- unname(sqrt(Vhc3[k, k]))
dist <- rbind(
  data.frame(metric = "Median locked probability (Black)", value = round(median(mb), 5),
             lo = round(unname(quantile(mb, .25)), 5), hi = round(unname(quantile(mb, .75)), 5),
             p = NA_real_, stringsAsFactors = FALSE),
  data.frame(metric = "Median locked probability (White)", value = round(median(mw), 5),
             lo = round(unname(quantile(mw, .25)), 5), hi = round(unname(quantile(mw, .75)), 5),
             p = NA_real_, stringsAsFactors = FALSE),
  data.frame(metric = "Wilcoxon rank-sum (Black vs White)", value = NA_real_, lo = NA_real_, hi = NA_real_,
             p = signif(wilx$p.value, 5), stringsAsFactors = FALSE),
  data.frame(metric = "Adjusted Black-vs-White difference (cohort SD; HC3)", value = round(adj_est, 5),
             lo = round(adj_est - 1.96 * adj_se, 5), hi = round(adj_est + 1.96 * adj_se, 5),
             p = signif(2 * pnorm(-abs(adj_est / adj_se)), 5), stringsAsFactors = FALSE))
dist$n_black <- length(mb); dist$n_white <- length(mw)
write.csv(dist, file.path(OUT, "Race_Stratified_MetScore_Durham_score_distribution.csv"), row.names = FALSE)

## ---- warning / convergence ledger ---------------------------------------
wl <- if (length(.warn$log)) do.call(rbind, .warn$log) else
  data.frame(context = "none", message = "no warnings captured", stringsAsFactors = FALSE)
write.csv(unique(wl), file.path(OUT, "Race_Stratified_MetScore_Durham_warnings.csv"), row.names = FALSE)

## ---- console summary -----------------------------------------------------
cat("\n== Durham race-stratified Met-Score analysis ==\n")
cat("ledger:\n"); print(ledger, row.names = FALSE)
cat(sprintf("\nBlack+White analysis set n=%d events=%d\n", nrow(BW), sum(BW$ev)))
cat("continuous interaction contrasts:\n"); print(contrast[, c("quantity", "estimate", "ci_lo", "ci_hi", "p")], row.names = FALSE)
cat(sprintf("primary global PH p = %.4f (Grade Group term p = %.4f) ; model converged = %s\n",
            ph_global, signif(unname(zph$table["gg", "p"]), 4), converged))
cat("Grade Group-stratified sensitivity contrasts:\n"); print(gg_contrasts[, c("quantity", "estimate", "ci_lo", "ci_hi", "p")], row.names = FALSE)
cat(sprintf("Grade Group-stratified global PH p = %.4f ; converged = %s\n", ph_global_g, converged_g))
cat("frozen High-vs-Low sensitivity:\n"); print(binsens[, c("quantity", "estimate", "ci_lo", "ci_hi", "p")], row.names = FALSE)
cat("score distribution:\n"); print(dist[, c("metric", "value", "lo", "hi", "p")], row.names = FALSE)
cat(sprintf("warnings captured: %d\n", length(.warn$log)))
cat("Done.\n")
