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
#   1. Hierarchical log-linear model (Light x Feeding x Time)       -> Fig. 2
#   2. Chi-squared tests by light regime (feeding groups)           -> Fig. 2
#   3. Chi-squared tests by light regime (larval instars)           -> Fig. 3
#   4. Synchrony index and Wilcoxon signed-ranks test               -> Fig. 4
#   5. Within-instar temporal shift: chi-squared tests              -> Fig. 5
#
# Required packages: MASS
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

# =============================================================================
# 1. Hierarchical log-linear model: Light x Feeding x Time (Fig. 2)
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
# 2. Chi-squared tests: feeding groups within each light regime (Fig. 2)
# =============================================================================
cat("-------------------------------------------------------------\n")
cat("2. Chi-squared tests: feeding groups within each light regime\n")
cat("   (pooled across instars; Fig. 2)\n")
cat("-------------------------------------------------------------\n")

for (lv in c("15.5L:8.5D", "24L")) {
  t <- chisq.test(tab3[lv, , ])
  cat(sprintf("  %-10s  \u03c7\u00b2 = %.2f, df = %d, P = %.4f\n",
              lv, t$statistic, t$parameter, t$p.value))
}
cat("\n")

# =============================================================================
# 3. Chi-squared tests: larval instars within each light regime (Fig. 3)
# =============================================================================
cat("-------------------------------------------------------------\n")
cat("3. Chi-squared tests: larval instars within each light regime\n")
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
# 4. Synchrony index and Wilcoxon signed-ranks test (Fig. 4)
# =============================================================================
cat("-------------------------------------------------------------\n")
cat("4. Synchrony index (Fig. 4)\n")
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

s1 <- sync_df$sync[sync_df$light == 1]   # 15.5L:8.5D
s2 <- sync_df$sync[sync_df$light == 2]   # 24L

cat(sprintf("  n per group: 15.5L:8.5D = %d,  24L = %d\n",
            length(s1), length(s2)))
cat(sprintf("  Median synchrony index: 15.5L:8.5D = %.3f,  24L = %.3f\n",
            median(s1), median(s2)))

wt <- wilcox.test(s1, s2, paired = TRUE, exact = FALSE)
z  <- qnorm(wt$p.value / 2, lower.tail = FALSE)
cat(sprintf("  Wilcoxon signed-ranks test: W = %.0f, z = %.2f, P = %.4f\n\n",
            wt$statistic, z, wt$p.value))

# =============================================================================
# 5. Within-instar temporal shift: chi-squared tests (Fig. 5)
# =============================================================================
cat("-------------------------------------------------------------\n")
cat("5. Within-instar temporal shift: chi-squared tests (Fig. 5)\n")
cat("   Comparison of ecdysis timing between the two days with the\n")
cat("   highest ecdysis counts within each instar x light regime\n")
cat("-------------------------------------------------------------\n")

cat(sprintf("  %-12s %-7s %-10s %8s %4s %8s\n",
            "Light", "Instar", "Days", "chi2", "df", "P"))
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
    t <- suppressWarnings(chisq.test(tab))
    cat(sprintf("  %-12s %-7s Day%d vs Day%d  %7.2f %4d %8.4f\n",
                lv, ins, top2[1], top2[2],
                t$statistic, t$parameter, t$p.value))
  }
}

cat("\n=============================================================\n")
cat("  Analysis complete.\n")
cat("=============================================================\n")
