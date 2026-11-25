# 04_overadjustment_checks.R
# (1) Effect-size attenuation, (2) negative-control covariates (random 7-gene scores),
# (3) within-group correlation of signature genes with the DA score.
# Random draws use fixed seeds; exact counts can differ slightly between R versions.
source("scripts/00_helpers.R")
v2 <- readRDS("data/v2_results.rds"); sens <- readRDS("data/sensitivity_results.rds")
de_raw <- v2$de_raw; de_adj <- v2$de_adj; res <- v2$res; scores <- v2$scores
common <- v2$common; core <- sens$core
cohort_res <- lapply(cohort_names, function(n) readRDS(file.path("data", paste0(n, "_results.rds"))))
names(cohort_res) <- cohort_names
grp  <- lapply(cohort_res, function(x) factor(x$group, levels = c("Control", "PD")))
mats <- lapply(cohort_names, function(n) collapse_probes(cohort_res[[n]]$rma, cohort_db[[n]]))
names(mats) <- cohort_names

set.seed(1)
run_one <- function(n, covar = NULL) {
  des <- if (is.null(covar)) model.matrix(~ grp[[n]]) else model.matrix(~ grp[[n]] + covar)
  topTable(eBayes(lmFit(mats[[n]], des)), coef = 2, number = Inf, sort.by = "none")
}
meta_count <- function(de) {
  m <- meta_rra(de)
  c(up = length(sigset(m, "up")), down = length(sigset(m, "dn")))
}

## 1. Attenuation for the 451 meta-analysis genes excluding the 7 score genes
sig <- setdiff(c(sigset(res$unadjusted, "up"), sigset(res$unadjusted, "dn")), da_markers)
att <- sapply(cohort_names, function(n) {
  tr <- de_raw[[n]][sig, "t"]; ta <- de_adj[[n]][sig, "t"]
  c(median_abs_t_ratio = median(abs(ta) / abs(tr)), same_sign = mean(sign(ta) == sign(tr)),
    still_p05 = mean(de_adj[[n]][sig, "P.Value"] < 0.05))
})
cat("\n== Attenuation (", length(sig), "genes) ==\n"); print(round(att, 2))

## 2. Negative controls
mean_expr <- rowMeans(sapply(mats, function(m) rowMeans(m[common, ])))
pool_mid  <- names(mean_expr)[mean_expr > median(mean_expr)]
meta_dn_sig <- setdiff(sigset(res$unadjusted, "dn"), c(da_markers, core))
run_control <- function(genes) {
  de <- lapply(cohort_names, function(n) run_one(n, score_of(mats[[n]], genes))); names(de) <- cohort_names
  meta_count(de)
}
nrep <- 20
ctlA <- t(replicate(nrep, run_control(sample(pool_mid, 7))))      # random expressed genes
ctlB <- t(replicate(nrep, run_control(sample(meta_dn_sig, 7))))   # random non-core, non-score signature genes
cat("\n== Control A =="); print(summary(ctlA)); cat("\n== Control B =="); print(summary(ctlB))
saveRDS(list(att = att, ctlA = ctlA, ctlB = ctlB), "data/overadjustment_checks.rds")

## 3. Within-group correlation (group means removed from score and from each gene)
set.seed(2)
sig_dn <- setdiff(sigset(res$unadjusted, "dn"), da_markers)
sig_up <- setdiff(sigset(res$unadjusted, "up"), da_markers)
rand   <- sample(setdiff(common, c(sig_dn, sig_up, da_markers)), 400)
gene_cors <- function(n, genes) {
  m  <- mats[[n]][intersect(genes, rownames(mats[[n]])), , drop = FALSE]
  rs <- resid(lm(scores[[n]] ~ grp[[n]]))
  rg <- t(apply(m, 1, function(y) resid(lm(y ~ grp[[n]]))))
  apply(rg, 1, function(y) cor(y, rs))
}
gc <- lapply(cohort_names, function(n)
  list(down_sig = gene_cors(n, sig_dn), up_sig = gene_cors(n, sig_up), random400 = gene_cors(n, rand)))
names(gc) <- cohort_names
out <- sapply(gc, function(x) sapply(x, median))
cat("\nMedian within-group correlation with the DA score\n"); print(round(out, 2))
saveRDS(list(median = out, gene_cors = gc), "data/within_group_cor.rds")
