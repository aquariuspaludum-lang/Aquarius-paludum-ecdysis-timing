#!/usr/bin/env Rscript
# =============================================================================
# 02_gating_analysis.R
#
# Statistical analyses for:
# "Light-dependent effects of feeding timing on the daily timing of ecdysis
#  in the water strider Aquarius paludum (Hemiptera: Gerridae)"
# Manabu Kishi
#
# This script quantifies the within-instar temporal shift of ecdysis timing
# (circadian gating evidence) by summarising daily ecdysis distributions
# within each instar under each light regime.
#
# Required packages: none (base R only)
# =============================================================================

# -----------------------------------------------------------------------------
# 0. Setup
# -----------------------------------------------------------------------------
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
cat("  Within-instar gating analysis\n")
cat("=============================================================\n\n")

# =============================================================================
# 1. Daily ecdysis distribution within each instar (pooled across feeding groups)
# =============================================================================
cat("-------------------------------------------------------------\n")
cat("1. Daily ecdysis proportions within each instar\n")
cat("   (pooled across feeding groups; days with >= 10 events shown)\n")
cat("-------------------------------------------------------------\n\n")

for (lv in c("15.5L:8.5D", "24L")) {
  cat(sprintf("=== %s ===\n", lv))
  for (ins in c("2nd", "3rd", "4th", "5th", "Adult")) {
    sub     <- dat[dat$Light == lv & dat$Instar == ins, ]
    days    <- sort(unique(sub$observation_day))
    cat(sprintf("  Instar %s:\n", ins))
    cat(sprintf("    %-8s %5s  %7s  %7s  %7s\n",
                "Day", "n", "Morning", "Midday", "Evening"))
    for (d in days) {
      dd  <- sub[sub$observation_day == d, ]
      tot <- sum(dd$n_ecdysis)
      if (tot < 10) next
      m   <- sum(dd$n_ecdysis[dd$observation_time == 1])
      mid <- sum(dd$n_ecdysis[dd$observation_time == 2])
      e   <- sum(dd$n_ecdysis[dd$observation_time == 3])
      cat(sprintf("    Day %3d  %5d  %6.1f%%  %6.1f%%  %6.1f%%\n",
                  d, tot,
                  m   / tot * 100,
                  mid / tot * 100,
                  e   / tot * 100))
    }
    cat("\n")
  }
}

# =============================================================================
# 2. Evening ecdysis proportion shift (Day 1 -> Day 2 within instar)
# =============================================================================
cat("-------------------------------------------------------------\n")
cat("2. Evening ecdysis proportion: Day 1 vs Day 2 within instar\n")
cat("   (positive shift = evening-dominant on Day 1, decline on Day 2)\n")
cat("-------------------------------------------------------------\n\n")

cat(sprintf("  %-12s %-7s %8s %8s %8s\n",
            "Light", "Instar", "Day1(%)", "Day2(%)", "Shift"))
cat(paste(rep("-", 52), collapse = ""), "\n")

for (lv in c("15.5L:8.5D", "24L")) {
  for (ins in c("2nd", "3rd", "4th", "5th", "Adult")) {
    sub     <- dat[dat$Light == lv & dat$Instar == ins, ]
    day_tot <- tapply(sub$n_ecdysis, sub$observation_day, sum)
    top2    <- sort(as.integer(names(sort(day_tot, decreasing = TRUE)[1:2])))
    if (length(top2) < 2) next

    get_eve_pct <- function(day) {
      dd  <- sub[sub$observation_day == day, ]
      tot <- sum(dd$n_ecdysis)
      if (tot == 0) return(NA)
      sum(dd$n_ecdysis[dd$observation_time == 3]) / tot * 100
    }
    e1 <- get_eve_pct(top2[1])
    e2 <- get_eve_pct(top2[2])
    cat(sprintf("  %-12s %-7s %7.1f%%  %7.1f%%  %+7.1f\n",
                lv, ins, e1, e2, e1 - e2))
  }
}
cat("\n")

# =============================================================================
# 3. Summary: proportion of ecdysis at each time across all instars and days
# =============================================================================
cat("-------------------------------------------------------------\n")
cat("3. Overall summary by light regime\n")
cat("-------------------------------------------------------------\n\n")

for (lv in c("15.5L:8.5D", "24L")) {
  sub <- dat[dat$Light == lv, ]
  tot <- sum(sub$n_ecdysis)
  m   <- sum(sub$n_ecdysis[sub$observation_time == 1])
  mid <- sum(sub$n_ecdysis[sub$observation_time == 2])
  e   <- sum(sub$n_ecdysis[sub$observation_time == 3])
  cat(sprintf("  %s (n = %d):\n", lv, tot))
  cat(sprintf("    Morning: %5.1f%%  Midday: %5.1f%%  Evening: %5.1f%%\n\n",
              m / tot * 100, mid / tot * 100, e / tot * 100))
}

cat("=============================================================\n")
cat("  Analysis complete.\n")
cat("=============================================================\n")
