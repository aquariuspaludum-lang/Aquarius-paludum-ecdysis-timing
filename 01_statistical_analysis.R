#!/usr/bin/env Rscript
# =============================================================================
# 01_statistical_analysis.R
#
# Statistical analyses for:
# "Light-dependent effects of feeding timing on the daily timing of ecdysis
#  in the water strider Aquarius paludum (Hemiptera: Gerridae)"
# Manabu Kishi
#
# Analyses performed:
#   1. Hierarchical log-linear model (Light x Feeding x Time)            -> Table 1, Fig. 2
#   2. Pearson chi-squared tests by light regime (feeding groups)        -> Table 1, Fig. 2
#   3. Pearson chi-squared tests by light regime (nymphal instars)       -> Fig. 3
#   4. Synchrony index and paired Wilcoxon signed-rank test              -> Fig. 5
#   5. Within-instar temporal shift: G2 (likelihood-ratio) tests         -> Table 2, Fig. 4
#
# Revision notes (2026-05-23):
#   - Section 4: paired Wilcoxon now computed WITHOUT continuity correction.
#                z is derived analytically from V_plus / E(V) / SD(V) with
#                tie correction, so the test statistic and the reported
#                two-sided p-value are mutually consistent.
#   - Section 5: switched from Pearson chi-squared to G2 (likelihood-ratio)
#                test for all 2x3 (day x observation time) tables. The G2
#                statistic is well-defined for cells with zero observations
#                (using the convention 0*log(0) = 0) and provides a
#                consistent treatment across all instars.
#
# Required packages: MASS
# R version: 4.5.1 (2025-06-13 ucrt)
# =============================================================================

# -----------------------------------------------------------------------------
# 0. Setup
# -----------------------------------------------------------------------------
if (!requireNamespace("MASS", quietly = TRUE)) install.packages("MASS")
library(MASS)

dat <- read.csv("Aquarius_paludum_ecdysis_timing.csv")
names(dat)[1] <- "light_regime"   # remove UTF-8 BOM if present

dat$Light   <- factor(dat$light_regime,
                      levels = c(1, 2),
                      labels = c("15.5L:8.5D", "24L"))
dat$Feeding <- factor(dat$feeding_group,
                      levels = 1:4,
                      labels = c("A", "B", "C", "D"))
dat$Instar  <- factor(dat$instar,
                      levels = c(2, 3, 4, 5, 6),
                      labels = c("2nd", "3rd", "4th", "5th", "Adult"))
dat$Time    <- factor(dat$observation_time,
                      levels = 1:3,
                      labels = c("Morning", "Midday", "Evening"))

cat("=============================================================\n")
cat("  Statistical analysis: Aquarius_paludum_ecdysis_timing.csv\n")
cat("=============================================================\n\n")
cat(sprintf("Total ecdysis events: %d\n\n", sum(dat$n_ecdysis)))

# -----------------------------------------------------------------------------
# Helper: G2 (likelihood-ratio) test for a contingency table
# -----------------------------------------------------------------------------
# Returns a list with statistic (G2), parameter (df), and p.value
# (from asymptotic chi-squared distribution). Cells with O = 0 contribute
# 0 to G2 (using the convention 0 * log(0) = 0).
G2_test <- function(tab) {
  O <- as.matrix(tab)
  row_tot <- rowSums(O)
  col_tot <- colSums(O)
  grand   <- sum(O)
  E       <- outer(row_tot, col_tot) / grand
  contrib <- ifelse(O > 0, O * log(O / E), 0)
  G2 <- 2 * sum(contrib)
  df <- (nrow(O) - 1) * (ncol(O) - 1)
  p  <- pchisq(G2, df, lower.tail = FALSE)
  list(statistic = G2, parameter = df, p.value = p)
}

# =============================================================================
# 1. Hierarchical log-linear model: Light x Feeding x Time (Table 1, Fig. 2)
# =============================================================================
cat("-------------------------------------------------------------\n")
cat("1. Hierarchical log-linear model: Light x Feeding x Time\n")
cat("-------------------------------------------------------------\n")

agg3 <- aggregate(n_ecdysis ~ Light + Feeding + Time, data = dat, FUN = sum)
tab3 <- xtabs(n_ecdysis ~ Light + Feeding + Time, data = agg3)

# Saturated model [LFT]
fit_sat <- loglm(~ Light * Feeding * Time, data = tab3)

# Model without 3-way interaction [LF][LT][FT]
fit_no3 <- loglm(~ Light * Feeding + Light * Time + Feeding * Time,
                 data = tab3)

# Test of 3-way interaction (difference in deviance)
G2_3way <- fit_no3$lrt - fit_sat$lrt
df_3way <- fit_no3$df  - fit_sat$df
p_3way  <- pchisq(G2_3way, df_3way, lower.tail = FALSE)

cat(sprintf("  Three-way interaction (Light x Feeding x Time):\n"))
cat(sprintf("  G\u00b2 = %.2f, df = %d, P = %.6f\n\n", G2_3way, df_3way, p_3way))

# Two-way interactions
fit_noLT <- loglm(~ Light * Feeding + Feeding * Time, data = tab3)
fit_noFT <- loglm(~ Light * Feeding + Light  * Time,  data = tab3)
fit_noLF <- loglm(~ Light * Time    + Feeding * Time,  data = tab3)

for (info in list(
  list(fit_noLT, "Light x Time"),
  list(fit_noFT, "Feeding x Time"),
  list(fit_noLF, "Light x Feeding")
)) {
  G2 <- info[[1]]$lrt - fit_no3$lrt
  df <- info[[1]]$df  - fit_no3$df
  p  <- pchisq(G2, df, lower.tail = FALSE)
  cat(sprintf("  %-20s G\u00b2 = %6.2f, df = %d, P = %.4f\n",
              info[[2]], G2, df, p))
}
cat("\n")

# =============================================================================
# 2. Pearson chi-squared tests: feeding groups within each light regime
#    (Table 1, Fig. 2)
# =============================================================================
cat("-------------------------------------------------------------\n")
cat("2. Pearson chi-squared tests: feeding groups within each light regime\n")
cat("   (pooled across instars; Table 1, Fig. 2)\n")
cat("-------------------------------------------------------------\n")

for (lv in c("15.5L:8.5D", "24L")) {
  t <- chisq.test(tab3[lv, , ])
  cat(sprintf("  %-10s  \u03c7\u00b2 = %.2f, df = %d, P = %.4f\n",
              lv, t$statistic, t$parameter, t$p.value))
}
cat("\n")

# =============================================================================
# 3. Pearson chi-squared tests: nymphal instars within each light regime
#    (Fig. 3)
# =============================================================================
cat("-------------------------------------------------------------\n")
cat("3. Pearson chi-squared tests: nymphal instars within each light regime\n")
cat("   (pooled across feeding groups; Fig. 3)\n")
cat("-------------------------------------------------------------\n")

agg3b <- aggregate(n_ecdysis ~ Light + Instar + Time, data = dat, FUN = sum)
tab3b <- xtabs(n_ecdysis ~ Light + Instar + Time, data = agg3b)

for (lv in c("15.5L:8.5D", "24L")) {
  t <- chisq.test(tab3b[lv, , ])
  cat(sprintf("  %-10s  \u03c7\u00b2 = %.2f, df = %d, P = %.6f\n",
              lv, t$statistic, t$parameter, t$p.value))
}
cat("\n")

# Print ecdysis counts per instar per light regime
cat("  Ecdysis counts per instar:\n")
for (lv in c("15.5L:8.5D", "24L")) {
  counts <- apply(tab3b[lv, , ], 1, sum)
  cat(sprintf("  %s: %s\n", lv,
              paste(names(counts), counts, sep = " = ", collapse = ";  ")))
}
cat("\n")

# =============================================================================
# 4. Synchrony index and paired Wilcoxon signed-rank test (Fig. 5)
#
# Pairing: each Feeding x Instar combination is paired between the two
# light regimes (n = 20 pairs). The test statistic V (sum of positive
# ranks) is taken from wilcox.test(); z is computed analytically from
# V_plus, its expected value E(V) = n(n+1)/4, and its variance with tie
# correction. No continuity correction is applied so that z and the
# two-sided p-value (2 * P(Z > |z|)) are mutually consistent.
# =============================================================================
cat("-------------------------------------------------------------\n")
cat("4. Synchrony index and paired Wilcoxon signed-rank test (Fig. 5)\n")
cat("   Index = 1 - H_norm, where H_norm = H / log(3)\n")
cat("   H = Shannon entropy of ecdysis proportions at 3 observation times\n")
cat("-------------------------------------------------------------\n")

rows <- list()
for (lr in 1:2) {
  for (fg in 1:4) {
    for (ins in c(2, 3, 4, 5, 6)) {
      sub <- dat[dat$light_regime == lr &
                 dat$feeding_group == fg &
                 dat$instar == ins, ]
      tot <- sum(sub$n_ecdysis)
      if (tot == 0) next
      pv <- tapply(sub$n_ecdysis, sub$Time, sum)
      pv[is.na(pv)] <- 0
      pv <- pv / sum(pv)
      H  <- -sum(ifelse(pv > 0, pv * log(pv), 0))
      rows[[length(rows) + 1]] <- data.frame(
        light   = lr,
        feeding = fg,
        instar  = ins,
        sync    = 1 - H / log(3)
      )
    }
  }
}
sync_df <- do.call(rbind, rows)

# Sort to guarantee one-to-one pairing by (Feeding, Instar)
sync_df <- sync_df[order(sync_df$feeding, sync_df$instar, sync_df$light), ]
s1 <- sync_df$sync[sync_df$light == 1]   # 15.5L:8.5D
s2 <- sync_df$sync[sync_df$light == 2]   # 24L

cat(sprintf("  n per group: 15.5L:8.5D = %d,  24L = %d\n",
            length(s1), length(s2)))
cat(sprintf("  Median synchrony index: 15.5L:8.5D = %.3f,  24L = %.3f\n",
            median(s1), median(s2)))

# Paired Wilcoxon signed-rank test WITHOUT continuity correction
wt <- wilcox.test(s1, s2, paired = TRUE,
                  exact = FALSE, correct = FALSE)

# Analytical z computed from V_plus (tie correction included)
d     <- s1 - s2
abs_d <- abs(d[d != 0])
n_nz  <- length(abs_d)
ranks_abs <- rank(abs_d)
V_plus <- sum(ranks_abs[d[d != 0] > 0])

E_V   <- n_nz * (n_nz + 1) / 4
ties  <- table(abs_d)
tie_correction <- sum(ties^3 - ties) / 48
SD_V  <- sqrt(n_nz * (n_nz + 1) * (2 * n_nz + 1) / 24 - tie_correction)
z_val <- (V_plus - E_V) / SD_V
p_val <- 2 * (1 - pnorm(abs(z_val)))

cat(sprintf("  Paired Wilcoxon signed-rank test (no continuity correction):\n"))
cat(sprintf("    V = %.0f, z = %.2f, P = %.4f\n",
            wt$statistic, z_val, p_val))
cat(sprintf("    (n = %d pairs)\n\n", length(s1)))

# =============================================================================
# 5. Within-instar temporal shift: G2 (likelihood-ratio) tests (Table 2, Fig. 4)
#
# For each instar within each light regime, the two days with the highest
# ecdysis counts are selected and a 2 x 3 contingency table (day x
# observation time) is constructed. The G2 (likelihood-ratio) statistic
# is computed and referred to the asymptotic chi-squared distribution
# with 2 degrees of freedom. G2 is used uniformly across all rows because
# it remains well defined when some observed cells contain zero counts
# (the cell contribution 0 * log(0) is treated as 0), and it provides a
# consistent treatment of all 2 x 3 tables irrespective of the magnitude
# of individual cells.
# =============================================================================
cat("-------------------------------------------------------------\n")
cat("5. Within-instar temporal shift: G2 (likelihood-ratio) tests\n")
cat("   (Table 2, Fig. 4)\n")
cat("-------------------------------------------------------------\n")

cat(sprintf("  %-12s %-7s %-13s %8s %4s %8s\n",
            "Light", "Instar", "Days", "G2", "df", "P"))
cat(paste(rep("-", 58), collapse = ""), "\n")

for (lv in c("15.5L:8.5D", "24L")) {
  for (ins in c("2nd", "3rd", "4th", "5th", "Adult")) {
    sub      <- dat[dat$Light == lv & dat$Instar == ins, ]
    day_tot  <- tapply(sub$n_ecdysis, sub$observation_day, sum)
    top2     <- sort(as.integer(names(sort(day_tot, decreasing = TRUE)[1:2])))

    tab <- matrix(0, nrow = 2, ncol = 3)
    for (i in 1:2) {
      d <- sub[sub$observation_day == top2[i], ]
      for (j in 1:3) {
        tab[i, j] <- sum(d$n_ecdysis[d$observation_time == j])
      }
    }
    t <- G2_test(tab)
    cat(sprintf("  %-12s %-7s Day%d vs Day%-3d %7.2f %4d %8.4f\n",
                lv, ins, top2[1], top2[2],
                t$statistic, t$parameter, t$p.value))
  }
}

cat("\n=============================================================\n")
cat("  Analysis complete.\n")
cat("=============================================================\n")
