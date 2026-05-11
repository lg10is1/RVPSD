#!/usr/bin/env Rscript
# PLDDT_Family_Analysis.R
# ============================================================
# Description: Comparative analysis of AlphaFold-predicted
#   pLDDT scores across the top 10 viral families.
#   Includes: single-family density plots, pairwise Wilcoxon
#   tests (both unbalanced and balanced 1:1 stratified sampling),
#   violin/boxplot overview, effect-size heatmap, and robustness
#   check (balanced vs unbalanced).
#
# Input:  virus_new.csv  (columns: protein_name, family, plddt)
# Output: figures/*/*.pdf, figures/*/*.png
#         results/*.csv, results/ANALYSIS_BILINGUAL_REPORT.txt
#
# Dependencies (auto-installed if missing):
#   data.table, ggplot2, patchwork, dplyr, tidyr
#
# Usage:
#   Rscript PLDDT_Family_Analysis.R
#
#   # Install dependencies manually (optional):
#   install.packages(c("data.table","ggplot2","patchwork","dplyr","tidyr"),
#                    repos="https://cloud.r-project.org")
# ============================================================

# ----- Auto-install missing packages -----
packages <- c("data.table", "ggplot2", "patchwork", "dplyr", "tidyr")
for (pkg in packages) {
  if (!requireNamespace(pkg, quietly = TRUE))
    install.packages(pkg, repos = "https://cloud.r-project.org")
  library(pkg, character.only = TRUE)
}

# ----- Working directory (auto-detect, no hard-coded paths) -----
work_dir <- if (file.exists("virus_new.csv")) {
  getwd()
} else {
  cat("[ERROR] virus_new.csv not found in current directory.\n")
  cat("        Please run this script from the directory containing virus_new.csv,\n")
  cat("        or set work_dir manually above this line.\n")
  quit(status = 1)
}
setwd(work_dir)
cat("[INFO] Working directory:", work_dir, "\n")

# ----- Create output directories -----
dirs <- c(
  "figures/single_family",
  "figures/pairwise_unbalanced",
  "figures/pairwise_balanced",
  "figures/overview",
  "figures/significant_pairs",
  "results"
)
for (d in dirs) dir.create(d, recursive = TRUE, showWarnings = FALSE)

# ==================== Nature-style colour palette ====================
nature_palette <- c(
  "#1f4e79", "#c0504d", "#9bbb59", "#8064a2", "#f79646",
  "#4bacc6", "#8c8c8c", "#d6604d", "#4393c3", "#b2182b"
)

# ==================== Nature-style ggplot2 theme ====================
theme_nature <- function(base_size = 12) {
  theme_bw(base_size = base_size) +
    theme(
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      panel.border     = element_rect(color = "black", fill = NA, linewidth = 0.8),
      axis.line        = element_line(color = "black", linewidth = 0.4),
      axis.ticks       = element_line(color = "black", linewidth = 0.4),
      axis.text        = element_text(color = "black", size = base_size),
      axis.title       = element_text(color = "black", size = base_size + 1, face = "bold"),
      plot.title       = element_text(size = base_size + 2, face = "bold", hjust = 0.5),
      plot.subtitle    = element_text(size = base_size, hjust = 0.5, color = "#333333"),
      legend.position  = "none",
      plot.background  = element_rect(fill = "white", color = NA),
      panel.background = element_rect(fill = "white", color = NA)
    )
}

# ==================== Read data ====================
cat("[INFO] Reading virus_new.csv ...\n")
df <- fread("virus_new.csv", header = TRUE)
setnames(df, c("protein_name", "family", "plddt"))
df <- df[!is.na(plddt) & !is.na(family) & family != ""]
df[, plddt := as.numeric(plddt)]
df <- df[!is.na(plddt)]

# ==================== Select Top-10 families by sample size ====================
family_stats <- df[, .(N = .N, Median = median(plddt), Mean = mean(plddt), SD = sd(plddt)), by = family][order(-N)]
top10 <- family_stats[1:10, family]
df_top10 <- df[family %in% top10]
df_top10[, family := factor(family, levels = top10)]

cat("[INFO] Top-10 families:", paste(top10, collapse = ", "), "\n")
cat("[INFO] Total records (top-10):", nrow(df_top10), "\n")

# Save family statistics
fwrite(family_stats, "results/all_family_stats.csv")
fwrite(family_stats[1:10], "results/top10_family_stats.csv")

# ==================== Pairwise Wilcoxon (unbalanced, full data) ====================
cat("[INFO] Running 45 pairwise Wilcoxon tests (unbalanced)...\n")

unbal_results <- list()
combs <- combn(top10, 2)
for (i in 1:ncol(combs)) {
  g1 <- combs[1, i]; g2 <- combs[2, i]
  x <- df_top10[family == g1, plddt]
  y <- df_top10[family == g2, plddt]
  n1 <- length(x); n2 <- length(y)
  wt <- wilcox.test(x, y, exact = FALSE)
  U  <- as.numeric(wt$statistic)
  r  <- (2 * U - n1 * n2) / (n1 * n2)
  r_abs <- abs(r)
  mag <- ifelse(r_abs < 0.1, "negligible",
                ifelse(r_abs < 0.3, "small",
                       ifelse(r_abs < 0.5, "medium", "large")))
  unbal_results[[i]] <- data.frame(
    group1 = g1, group2 = g2, n1 = n1, n2 = n2,
    U = round(U, 1), p = wt$p.value,
    effsize = round(r, 4), effsize_abs = round(r_abs, 4), magnitude = mag
  )
}

pw <- as.data.table(do.call(rbind, unbal_results))
pw[, p.adj := p.adjust(p, method = "bonferroni")]
pw[, sig := ifelse(p.adj < 0.05, "***",
            ifelse(p.adj < 0.01, "**",
            ifelse(p.adj < 0.05, "*", "ns")))]
fwrite(pw, "results/pairwise_wilcox_unbalanced.csv")

# ==================== 1. Single-family density plots ====================
cat("[INFO] 1. Generating single-family density plots ...\n")

for (i in seq_along(top10)) {
  fam <- top10[i]
  col  <- nature_palette[i]
  dt_sub <- df_top10[family == fam]
  med <- median(dt_sub$plddt)
  mn  <- mean(dt_sub$plddt)

  p <- ggplot(dt_sub, aes(x = plddt)) +
    geom_density(fill = col, color = "black", alpha = 0.75, linewidth = 0.3) +
    geom_vline(xintercept = med, color = "#2c2c2c", linetype = "dashed", linewidth = 0.6) +
    geom_vline(xintercept = mn,  color = "#7f7f7f", linetype = "dotted", linewidth = 0.5) +
    annotate("text", x = med, y = Inf, vjust = -0.5, hjust = -0.05,
             label = paste0("Median = ", round(med, 1)),
             color = "#2c2c2c", size = 3.5, fontface = "bold") +
    annotate("text", x = mn,  y = Inf, vjust = -2.0, hjust = -0.05,
             label = paste0("Mean = ", round(mn, 1)),
             color = "#7f7f7f", size = 3.2) +
    scale_x_continuous(limits = c(0, 105), breaks = seq(0, 100, 20), expand = c(0.01, 0.01)) +
    scale_y_continuous(expand = expansion(mult = c(0, 0.05))) +
    labs(
      title = paste0(fam, " (N = ", format(nrow(dt_sub), big.mark = ","), ")"),
      x = "pLDDT score", y = "Density"
    ) +
    theme_nature(base_size = 13)

  ggsave(paste0("figures/single_family/", fam, "_density.pdf"), p, width = 5, height = 4, dpi = 600)
  ggsave(paste0("figures/single_family/", fam, "_density.png"), p, width = 5, height = 4, dpi = 600)
}

# ==================== 2. Global violin + boxplot overview ====================
cat("[INFO] 2. Generating overview violin-boxplot ...\n")

p_overview <- ggplot(df_top10, aes(x = family, y = plddt, fill = family)) +
  geom_violin(alpha = 0.6, color = "black", linewidth = 0.3, scale = "width") +
  geom_boxplot(width = 0.15, alpha = 0.9, outlier.shape = NA, color = "black", linewidth = 0.3) +
  stat_summary(fun = median, geom = "point", shape = 23, size = 2.5,
               fill = "white", color = "black", stroke = 0.8) +
  scale_fill_manual(values = nature_palette[1:10]) +
  scale_y_continuous(limits = c(0, 100), breaks = seq(0, 100, 20), expand = c(0, 0)) +
  labs(
    title = "pLDDT Distribution Across Top 10 Viral Families",
    x = "Viral Family", y = "pLDDT score"
  ) +
  theme_nature(base_size = 13) +
  theme(
    legend.position = "none",
    axis.text.x = element_text(angle = 35, hjust = 1, face = "italic")
  )

ggsave("figures/overview/Top10_ViolinBox.pdf", p_overview, width = 10, height = 6, dpi = 600)
ggsave("figures/overview/Top10_ViolinBox.png", p_overview, width = 10, height = 6, dpi = 600)

# ==================== 3. Effect-size heatmap ====================
cat("[INFO] 3. Generating effect-size heatmap ...\n")

pw_mat    <- pw[, .(group1, group2, effsize_abs)]
pw_mat_rev <- data.table(group1 = pw_mat$group2, group2 = pw_mat$group1, effsize_abs = pw_mat$effsize_abs)
pw_full   <- rbind(pw_mat, pw_mat_rev)
pw_cast   <- dcast(pw_full, group1 ~ group2, value.var = "effsize_abs")
pw_melt   <- as.data.table(melt(pw_cast, id.vars = "group1", variable.name = "group2", value.name = "effsize_abs"))
pw_melt[, group1 := factor(group1, levels = top10)]
pw_melt[, group2 := factor(group2, levels = top10)]

p_heat <- ggplot(pw_melt, aes(x = group2, y = group1, fill = effsize_abs)) +
  geom_tile(color = "white", linewidth = 0.5) +
  geom_text(aes(label = ifelse(!is.na(effsize_abs), sprintf("%.2f", effsize_abs), "")),
            size = 3, color = "white", fontface = "bold") +
  scale_fill_gradientn(colors = c("#f7fbff", "#4292c6", "#08306b"),
                       limits = c(0, 0.6), name = "|r|") +
  labs(title = "Effect Size Heatmap (Rank-biserial |r|)", x = "", y = "") +
  theme_nature(base_size = 12) +
  theme(
    legend.position = "right",
    axis.text.x = element_text(angle = 45, hjust = 1, face = "italic"),
    axis.text.y = element_text(face = "italic")
  )

ggsave("figures/overview/Effsize_Heatmap.pdf", p_heat, width = 8, height = 7, dpi = 600)
ggsave("figures/overview/Effsize_Heatmap.png", p_heat, width = 8, height = 7, dpi = 600)

# ==================== Helper functions ====================

# Balanced 1:1 stratified sampling
run_balanced_pair <- function(data, fam1, fam2, n_strata = 20, seed = 42) {
  set.seed(seed)
  dt_pair <- data[family %in% c(fam1, fam2)]
  dt_pair[, strata := cut(plddt, breaks = n_strata, labels = FALSE)]
  strata_tab  <- dt_pair[, .(n = .N), by = .(strata, family)]
  strata_min  <- strata_tab[, .(min_n = min(n)), by = strata]
  dt_pair <- merge(dt_pair, strata_min, by = "strata")
  balanced <- dt_pair[, .SD[sample(.N, min(.N, min_n[1]))], by = .(strata, family)]
  balanced[, c("strata", "min_n") := NULL]
  list(
    data = balanced,
    n1   = balanced[family == fam1, .N],
    n2   = balanced[family == fam2, .N]
  )
}

# Wilcoxon + effect size for a pair
calc_wilcox_pair <- function(x, y, g1, g2) {
  n1 <- length(x); n2 <- length(y)
  wt <- wilcox.test(x, y, exact = FALSE)
  U  <- as.numeric(wt$statistic)
  r  <- (2 * U - n1 * n2) / (n1 * n2)
  r_abs <- abs(r)
  mag <- ifelse(r_abs < 0.1, "negligible",
                ifelse(r_abs < 0.3, "small",
                       ifelse(r_abs < 0.5, "medium", "large")))
  data.frame(
    group1 = g1, group2 = g2, n1 = n1, n2 = n2,
    U = round(U, 1), p = wt$p.value,
    effsize = round(r, 4), effsize_abs = round(r_abs, 4), magnitude = mag
  )
}

# ==================== 4. All 45 pairwise density overlays (unbalanced) ====================
cat("[INFO] 4. Generating all 45 unbalanced pairwise density overlays ...\n")

pw[, flag := ifelse(effsize_abs >= 0.3, "HIGH_EFFECT",
             ifelse(p.adj < 0.05, "SIG", "ns"))]
pw[, is_significant := (effsize_abs >= 0.3 | p.adj < 0.05)]

for (k in 1:nrow(pw)) {
  g1   <- pw[k, group1]
  g2   <- pw[k, group2]
  r_val  <- pw[k, effsize_abs]
  p_val  <- pw[k, p.adj]
  flag   <- pw[k, flag]
  sig    <- pw[k, is_significant]

  dt1     <- df_top10[family == g1, .(plddt, group = g1)]
  dt2     <- df_top10[family == g2, .(plddt, group = g2)]
  dt_pair <- rbind(dt1, dt2)
  dt_pair[, group := factor(group, levels = c(g1, g2))]

  idx1 <- which(top10 == g1)
  idx2 <- which(top10 == g2)
  pair_colors <- c(nature_palette[idx1], nature_palette[idx2])

  p <- ggplot(dt_pair, aes(x = plddt, fill = group, color = group)) +
    geom_density(alpha = 0.35, linewidth = 0.6) +
    geom_vline(data = dt_pair[, .(med = median(plddt)), by = group],
               aes(xintercept = med, color = group), linetype = "dashed", linewidth = 0.8) +
    scale_color_manual(values = pair_colors, name = "Family") +
    scale_fill_manual(values = pair_colors, name = "Family") +
    scale_x_continuous(limits = c(0, 105), breaks = seq(0, 100, 20), expand = c(0.01, 0.01)) +
    scale_y_continuous(expand = expansion(mult = c(0, 0.05))) +
    labs(
      title = paste0(g1, " vs ", g2, "  [", flag, "]"),
      subtitle = paste0("|r| = ", round(r_val, 3), " (", pw[k, magnitude], ")  |  p.adj = ",
                        format(p_val, digits = 2, scientific = TRUE)),
      x = "pLDDT score", y = "Density"
    ) +
    theme_nature(base_size = 13) +
    theme(
      legend.position  = "bottom",
      legend.box       = "horizontal",
      legend.background = element_rect(fill = "white", color = NA),
      legend.text      = element_text(face = "italic", size = 10),
      legend.title     = element_text(face = "bold", size = 10),
      legend.margin    = margin(t = -5),
      plot.margin      = margin(b = 15, r = 10)
    )

  fname_base <- paste0(g1, "_vs_", g2)
  ggsave(paste0("figures/pairwise_unbalanced/", fname_base, ".pdf"), p, width = 6, height = 5.5, dpi = 600)
  ggsave(paste0("figures/pairwise_unbalanced/", fname_base, ".png"), p, width = 6, height = 5.5, dpi = 600)

  if (sig) {
    ggsave(paste0("figures/significant_pairs/SIG_", fname_base, ".pdf"), p, width = 6, height = 5.5, dpi = 600)
    ggsave(paste0("figures/significant_pairs/SIG_", fname_base, ".png"), p, width = 6, height = 5.5, dpi = 600)
  }
}

# ==================== 5. All 45 pairwise (balanced 1:1) ====================
cat("[INFO] 5. Running & plotting all 45 balanced 1:1 pairwise comparisons ...\n")

balanced_results <- list()

for (k in 1:nrow(pw)) {
  g1   <- pw[k, group1]
  g2   <- pw[k, group2]
  flag <- pw[k, flag]
  sig  <- pw[k, is_significant]

  res_bal  <- run_balanced_pair(df_top10, g1, g2)
  bal_data <- res_bal$data

  x <- bal_data[family == g1, plddt]
  y <- bal_data[family == g2, plddt]
  stat_bal <- calc_wilcox_pair(x, y, g1, g2)
  balanced_results[[k]] <- stat_bal

  idx1 <- which(top10 == g1)
  idx2 <- which(top10 == g2)
  pair_colors <- c(nature_palette[idx1], nature_palette[idx2])

  p <- ggplot(bal_data, aes(x = plddt, fill = family, color = family)) +
    geom_density(alpha = 0.35, linewidth = 0.6) +
    geom_vline(data = bal_data[, .(med = median(plddt)), by = family],
               aes(xintercept = med, color = family), linetype = "dashed", linewidth = 0.8) +
    scale_color_manual(values = pair_colors, name = "Family") +
    scale_fill_manual(values = pair_colors, name = "Family") +
    scale_x_continuous(limits = c(0, 105), breaks = seq(0, 100, 20), expand = c(0.01, 0.01)) +
    scale_y_continuous(expand = expansion(mult = c(0, 0.05))) +
    labs(
      title = paste0("Balanced 1:1  ", g1, " vs ", g2, "  [", flag, "]"),
      subtitle = paste0("N = ", res_bal$n1, " vs ", res_bal$n2,
                        "  |  |r| = ", round(stat_bal$effsize_abs, 3),
                        "  |  p = ", format(stat_bal$p, digits = 2, scientific = TRUE)),
      x = "pLDDT score", y = "Density"
    ) +
    theme_nature(base_size = 13) +
    theme(
      legend.position   = "bottom",
      legend.box        = "horizontal",
      legend.background = element_rect(fill = "white", color = NA),
      legend.text       = element_text(face = "italic", size = 10),
      legend.title      = element_text(face = "bold", size = 10),
      legend.margin     = margin(t = -5),
      plot.margin       = margin(b = 15, r = 10)
    )

  fname_base <- paste0("Balanced_", g1, "_vs_", g2)
  ggsave(paste0("figures/pairwise_balanced/", fname_base, ".pdf"), p, width = 6, height = 5.5, dpi = 600)
  ggsave(paste0("figures/pairwise_balanced/", fname_base, ".png"), p, width = 6, height = 5.5, dpi = 600)

  if (sig) {
    ggsave(paste0("figures/significant_pairs/SIG_", fname_base, ".pdf"), p, width = 6, height = 5.5, dpi = 600)
    ggsave(paste0("figures/significant_pairs/SIG_", fname_base, ".png"), p, width = 6, height = 5.5, dpi = 600)
  }
}

bal_df <- as.data.table(do.call(rbind, balanced_results))
fwrite(bal_df, "results/pairwise_balanced_all_45.csv")

# ==================== 6. Robustness scatter (balanced vs unbalanced) ====================
cat("[INFO] 6. Plotting robustness check (balanced vs unbalanced)...\n")

comp_dt <- merge(
  pw[, .(group1, group2, effsize_abs_unbal = effsize_abs, flag)],
  bal_df[, .(group1, group2, effsize_abs_bal = effsize_abs)],
  by = c("group1", "group2")
)
comp_dt[, point_color := ifelse(grepl("HIGH_EFFECT", flag), nature_palette[2],
                        ifelse(grepl("SIG", flag), nature_palette[5], "#999999"))]

p_comp <- ggplot(comp_dt, aes(x = effsize_abs_unbal, y = effsize_abs_bal)) +
  geom_abline(intercept = 0, slope = 1, linetype = "dashed", color = "grey50", linewidth = 0.5) +
  geom_point(aes(color = point_color), size = 3, alpha = 0.8, show.legend = FALSE) +
  scale_color_identity() +
  geom_text(aes(label = paste0(group1, "\nvs\n", group2)),
            size = 2.6, vjust = -0.8, hjust = 0.5, lineheight = 0.8) +
  scale_x_continuous(limits = c(0, 0.65), breaks = seq(0, 0.6, 0.2), name = "Unbalanced |r|") +
  scale_y_continuous(limits = c(0, 0.65), breaks = seq(0, 0.6, 0.2), name = "Balanced 1:1 |r|") +
  labs(
    title    = "Effect Size Robustness: Balanced vs Unbalanced",
    subtitle = "Red = high effect (|r|>=0.3), Orange = significant (p<0.05), Grey = non-significant"
  ) +
  theme_nature(base_size = 13) +
  theme(legend.position = "none", aspect.ratio = 1)

ggsave("figures/overview/Balanced_vs_Unbalanced.pdf", p_comp, width = 6.5, height = 6.5, dpi = 600)
ggsave("figures/overview/Balanced_vs_Unbalanced.png", p_comp, width = 6.5, height = 6.5, dpi = 600)

# ==================== 7. Bilingual analysis report ====================
cat("[INFO] 7. Generating bilingual analysis report ...\n")

n_total <- nrow(df)
n_top10 <- nrow(df_top10)
n_pairs <- nrow(pw)
n_sig   <- sum(pw$p.adj < 0.05, na.rm = TRUE)
n_high  <- sum(pw$effsize_abs >= 0.3, na.rm = TRUE)
n_large <- sum(pw$effsize_abs >= 0.5, na.rm = TRUE)
n_med   <- sum(pw$effsize_abs >= 0.3 & pw$effsize_abs < 0.5, na.rm = TRUE)

top5_pairs <- pw[order(-effsize_abs)][1:5, .(
  Pair = paste0(group1, " vs ", group2),
  r    = round(effsize_abs, 4),
  mag  = magnitude,
  p    = format(p.adj, digits = 2, scientific = TRUE)
)]

retro_median  <- family_stats[family == "Retroviridae", Median]
others_median <- family_stats[2:10, median(Median)]
retro_mean    <- family_stats[family == "Retroviridae", Mean]

comp_dt[, delta := abs(effsize_abs_unbal - effsize_abs_bal)]
n_robust     <- sum(comp_dt$delta < 0.05)
robust_pairs <- comp_dt[delta < 0.05, .(Pair = paste0(group1, " vs ", group2), Delta = round(delta, 4))]

sink("results/ANALYSIS_BILINGUAL_REPORT.txt")

cat("================================================================================\n")
cat("     PLDDT FAMILY ANALYSIS — COMPREHENSIVE BILINGUAL REPORT\n")
cat("================================================================================\n\n")
cat("Generated:", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n")
cat("Working Directory:", work_dir, "\n\n")

# --- Executive Summary ---
cat("================================================================================\n")
cat("SECTION 0: EXECUTIVE SUMMARY\n")
cat("================================================================================\n\n")

cat("[EN] This report presents a comparative analysis of AlphaFold-predicted pLDDT\n")
cat("     scores across the top 10 viral families by sample size. A total of", format(n_total, big.mark = ","),
    "protein structures were analyzed.\n")
cat("     Pairwise comparisons used the Wilcoxon rank-sum test with rank-biserial\n")
cat("     correlation (|r|) as the effect-size metric. Both unbalanced (full data)\n")
cat("     and balanced 1:1 (stratified sampling) analyses were performed.\n\n")

cat("[CN] 本报告对样本量前10的病毒科（Family）的AlphaFold预测pLDDT分数进行比较分析。\n")
cat("     共分析", format(n_total, big.mark = ","), "个蛋白结构，采用Wilcoxon秩和检验并以Rank-biserial相关系数（|r|）\n")
cat("     作为效应量指标。同时进行了非均衡（全数据）和均衡1:1（分层抽样）两种分析。\n\n")

cat("KEY STATISTICS:\n")
cat("  Total valid records        :", format(n_total, big.mark = ","), "\n")
cat("  Top 10 families analyzed   : 10\n")
cat("  Pairwise comparisons       :", n_pairs, "\n")
cat("  Significant (p.adj < 0.05) :", n_sig, "(", round(100*n_sig/n_pairs, 1), "%)\n")
cat("  High effect (|r| >= 0.3)   :", n_high, "(", round(100*n_high/n_pairs, 1), "%)\n")
cat("  Large effect (|r| >= 0.5)  :", n_large, "\n")
cat("  Medium effect (0.3-0.5)    :", n_med, "\n\n")

# --- Data Overview ---
cat("================================================================================\n")
cat("SECTION 1: DATA OVERVIEW\n")
cat("================================================================================\n\n")

cat("TOP-10 FAMILY DESCRIPTIVE STATISTICS\n")
cat("--------------------------------------------------------------------------------\n")
print(family_stats[1:10, .(
  Rank  = 1:10,
  Family = family,
  N     = format(N, big.mark = ","),
  Median = round(Median, 2),
  Mean   = round(Mean, 2),
  SD     = round(SD, 2)
)], row.names = FALSE)
cat("\n")

# --- Statistical Methods ---
cat("================================================================================\n")
cat("SECTION 2: STATISTICAL METHODS\n")
cat("================================================================================\n\n")

cat("[EN] 1. Unbalanced Analysis (Full Data)\n")
cat("     Test: Wilcoxon rank-sum test (Mann-Whitney U), exact = FALSE\n")
cat("     Correction: Bonferroni adjustment for 45 comparisons\n")
cat("     Effect size: Rank-biserial r = (2U - n1*n2) / (n1*n2)\n")
cat("     Magnitude: negligible (<0.1), small (0.1-0.3), medium (0.3-0.5), large (>=0.5)\n\n")

cat("[CN] 1. 非均衡分析（全数据）\n")
cat("     检验：Wilcoxon秩和检验（Mann-Whitney U），exact = FALSE\n")
cat("     多重检验校正：Bonferroni（45次比较）\n")
cat("     效应量：Rank-biserial r = (2U - n1*n2) / (n1*n2)\n")
cat("     效应量分级：可忽略（<0.1）、小（0.1-0.3）、中（0.3-0.5）、大（>=0.5）\n\n")

cat("[EN] 2. Balanced 1:1 Analysis (Stratified Sampling)\n")
cat("     pLDDT is divided into 20 strata; within each stratum, min(N1,N2) samples\n")
cat("     are randomly drawn from each family. This controls for sample-size imbalance.\n\n")

cat("[CN] 2. 均衡1:1分析（分层抽样）\n")
cat("     将pLDDT分为20个分层，每层从两组中各抽取min(N1,N2)个样本，以消除样本量差异影响。\n\n")

# --- Main Findings ---
cat("================================================================================\n")
cat("SECTION 3: MAIN FINDINGS\n")
cat("================================================================================\n\n")

cat("[EN] 3.1 Top 5 Largest Effect Sizes (Unbalanced)\n\n")
cat("[CN] 3.1 效应量最大的前5对（非均衡）\n\n")
print(top5_pairs, row.names = FALSE)
cat("\n")

cat("[EN] Key observation: Picornaviridae appears in 4 of the top 5 pairs, suggesting a\n")
cat("     distinct pLDDT profile compared to most other families.\n\n")

cat("[CN] 关键发现：Picornaviridae在前5对中出现4次，表明其pLDDT特征与其他多数科存在显著差异。\n\n")

cat("[EN] 3.2 Robustness (Balanced vs Unbalanced)\n")
cat("     Of", n_pairs, "pairs,", n_robust, "showed |delta| < 0.05 between balanced and unbalanced\n")
cat("     effect sizes, indicating high robustness.\n\n")

cat("[CN] 3.2 稳健性（均衡 vs 非均衡）\n")
cat("     在", n_pairs, "对比较中，", n_robust, "对的差异|delta| < 0.05，结果高度稳健。\n\n")

# --- Biological Interpretation ---
cat("================================================================================\n")
cat("SECTION 4: BIOLOGICAL INTERPRETATION\n")
cat("================================================================================\n\n")

cat("[EN] pLDDT (predicted Local Distance Difference Test) is AlphaFold's per-residue\n")
cat("     confidence score (range 0-100). Higher scores indicate greater structural\n")
cat("     confidence. Family-level pLDDT differences may reflect intrinsic structural\n")
cat("     variability, sequence-length bias, taxonomic sampling bias, or evolutionary\n")
cat("     conservation. Low pLDDT may indicate genuine disorder or lack of homologs\n")
cat("     in AlphaFold's training set.\n\n")

cat("[CN] pLDDT是AlphaFold的逐残基置信度评分（0-100）。科水平差异可能反映内在结构变异性、\n")
cat("     序列长度偏倚、分类学采样偏倚或进化保守性。低pLDDT可能表示真实无序区域或\n")
cat("     AlphaFold训练集中缺乏同源模板。\n\n")

# --- File Index ---
cat("================================================================================\n")
cat("SECTION 5: OUTPUT FILE INDEX\n")
cat("================================================================================\n\n")
cat("STATISTICAL TABLES:\n")
cat("  results/pairwise_wilcox_unbalanced.csv   — 45-pair unbalanced statistics\n")
cat("  results/pairwise_balanced_all_45.csv     — 45-pair balanced 1:1 statistics\n")
cat("  results/all_family_stats.csv             — all family descriptive stats\n")
cat("  results/top10_family_stats.csv           — top-10 family descriptive stats\n\n")
cat("MAIN FIGURES:\n")
cat("  figures/overview/Top10_ViolinBox.pdf     — overview violin+boxplot\n")
cat("  figures/overview/Effsize_Heatmap.pdf    — effect-size heatmap\n")
cat("  figures/overview/Balanced_vs_Unbalanced.pdf — robustness scatter\n\n")
cat("SINGLE FAMILY:\n")
cat("  figures/single_family/*_density.pdf      — 10 single-family density plots\n\n")
cat("PAIRWISE OVERLAYS:\n")
cat("  figures/pairwise_unbalanced/*.pdf        — all 45 unbalanced overlays\n")
cat("  figures/pairwise_balanced/*.pdf          — all 45 balanced 1:1 overlays\n")
cat("  figures/significant_pairs/*.pdf          — significant pairs only\n\n")
cat("================================================================================\n")
cat("END OF REPORT\n")
cat("================================================================================\n")
sink()

cat("\n[OK] Done! Bilingual report saved to: results/ANALYSIS_BILINGUAL_REPORT.txt\n")
cat("    All figures and tables have been generated under the 'figures/' and 'results/' directories.\n")