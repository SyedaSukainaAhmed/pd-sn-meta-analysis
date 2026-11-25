# 07_figures.R - regenerates Figures 1-7 as PNG files in figures/.
# NOTE: rebuilt from the saved result files; styling may differ slightly from the figures in the preprint.
suppressPackageStartupMessages(library(ggplot2))
source("scripts/00_helpers.R")
dir.create("figures", showWarnings = FALSE)
th <- theme_bw(base_size = 11)
v2 <- readRDS("data/v2_results.rds"); oc <- readRDS("data/overadjustment_checks.rds")
wg <- readRDS("data/within_group_cor.rds"); pw <- readRDS("data/pathway_unsimplified.rds")
rn <- readRDS("data/RNAseq_GSE136666_results.rds")
cohort_res <- lapply(cohort_names, function(n) readRDS(file.path("data", paste0(n, "_results.rds"))))
names(cohort_res) <- cohort_names

# Fig 1: dopaminergic score by group, per cohort
d1 <- do.call(rbind, lapply(cohort_names, function(n) data.frame(cohort = n, score = v2$scores[[n]],
        group = factor(cohort_res[[n]]$group, levels = c("Control", "PD")))))
p1 <- ggplot(d1, aes(group, score, colour = group)) + geom_boxplot(outlier.shape = NA) + geom_jitter(width = 0.15) +
  facet_wrap(~ cohort, nrow = 1) + labs(y = "Dopaminergic marker score (mean z)", x = NULL) + th + theme(legend.position = "none")
ggsave("figures/Fig1_DA_score.png", p1, width = 10, height = 3.5, dpi = 300)

# Fig 2: genes significant in the meta-analysis, without and with adjustment
d2 <- data.frame(analysis = rep(c("Unadjusted", "DA-adjusted"), each = 2), direction = rep(c("Down", "Up"), 2),
  n = c(length(sigset(v2$res$unadjusted, "dn")), length(sigset(v2$res$unadjusted, "up")),
        length(sigset(v2$res$DA_adjusted, "dn")), length(sigset(v2$res$DA_adjusted, "up"))))
p2 <- ggplot(d2, aes(analysis, n, fill = direction)) + geom_col(position = "dodge") + geom_text(aes(label = n), position = position_dodge(0.9), vjust = -0.3) +
  labs(y = "Genes significant in meta-analysis (BH < 0.05)", x = NULL) + th
ggsave("figures/Fig2_significant_genes.png", p2, width = 5, height = 4, dpi = 300)

# Fig 3: stability (top-200 overlap, mean of up and down lists)
st <- rbind(cbind(v2$stab_raw, analysis = "Unadjusted"), cbind(v2$stab_adj, analysis = "DA-adjusted"))
st$overlap <- (st$top200_up + st$top200_dn) / 2
st$left_out <- factor(st$n_dropped, labels = c("1 cohort left out", "2 cohorts left out"))
p3 <- ggplot(st, aes(left_out, overlap, colour = analysis)) + geom_boxplot(outlier.shape = NA) + geom_jitter(width = 0.15) +
  labs(y = "Top-200 overlap with full-data result", x = NULL) + th
ggsave("figures/Fig3_stability.png", p3, width = 6, height = 4, dpi = 300)

# Fig 4: negative controls
d4 <- rbind(data.frame(set = "No covariate", n = 451),
  data.frame(set = "Control A: random genes", n = rowSums(oc$ctlA)),
  data.frame(set = "Control B: random signature genes", n = rowSums(oc$ctlB)),
  data.frame(set = "DA score", n = 0))
p4 <- ggplot(d4, aes(set, n)) + geom_jitter(width = 0.15) + stat_summary(fun = median, geom = "crossbar", colour = "red", width = 0.5) +
  labs(y = "Genes significant in meta-analysis", x = NULL) + th + theme(axis.text.x = element_text(angle = 20, hjust = 1))
ggsave("figures/Fig4_negative_controls.png", p4, width = 6, height = 4.5, dpi = 300)

# Fig 5: within-group correlation with the score
d5 <- do.call(rbind, lapply(cohort_names, function(n) do.call(rbind, lapply(names(wg$gene_cors[[n]]), function(s)
  data.frame(cohort = n, set = s, r = unname(wg$gene_cors[[n]][[s]]))))))
p5 <- ggplot(d5, aes(cohort, r, fill = set)) + geom_boxplot(outlier.size = 0.3) + geom_hline(yintercept = 0, linetype = 2) +
  labs(y = "Within-group correlation with DA score", x = NULL, fill = NULL) + th
ggsave("figures/Fig5_within_group_correlation.png", p5, width = 9, height = 4, dpi = 300)

# Fig 6: microarray vs RNA-seq
d6 <- data.frame(mean_t = rn$mean_t[rn$shared], wald = rn$rs[rn$shared, "stat"])
p6 <- ggplot(d6, aes(mean_t, wald)) + geom_point(alpha = 0.2, size = 0.5) + geom_smooth(method = "lm", colour = "red", se = FALSE) +
  labs(x = "Mean t across 5 microarray cohorts", y = "DESeq2 Wald statistic (GSE136666)") + th
ggsave("figures/Fig6_RNAseq_concordance.png", p6, width = 5, height = 4.5, dpi = 300)

# Fig 7: gene vs pathway reproducibility
cm <- pw$cmp; d7 <- do.call(rbind, lapply(1:2, function(i) data.frame(level = cm$level[i],
  threshold = c(">= 2 cohorts", ">= 3 cohorts", "all 5"), fraction = unlist(cm[i, c("in_ge2_cohorts", "in_ge3_cohorts", "in_all_5")]))))
p7 <- ggplot(d7, aes(threshold, fraction, fill = level)) + geom_col(position = "dodge") +
  labs(y = "Fraction of meta-analysis items found in cohorts", x = NULL, fill = NULL) + th
ggsave("figures/Fig7_gene_vs_pathway.png", p7, width = 6, height = 4, dpi = 300)

writeLines(capture.output(sessionInfo()), "sessionInfo.txt")
