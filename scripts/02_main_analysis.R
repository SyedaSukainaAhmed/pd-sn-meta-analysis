# 02_main_analysis.R
# Per-cohort limma (with and without the dopaminergic score), robust rank aggregation,
# leave-one/two-cohort-out stability. Writes results/*.csv and data/v2_results.rds.
source("scripts/00_helpers.R")
dir.create("results", showWarnings = FALSE)

cohort_res <- lapply(cohort_names, function(n) readRDS(file.path("data", paste0(n, "_results.rds"))))
names(cohort_res) <- cohort_names
grp  <- lapply(cohort_res, function(x) factor(x$group, levels = c("Control", "PD")))
mats <- lapply(cohort_names, function(n) collapse_probes(cohort_res[[n]]$rma, cohort_db[[n]]))
names(mats) <- cohort_names
print(sapply(mats, dim))

scores <- lapply(mats, da_score)
for (n in cohort_names) {
  tt <- t.test(scores[[n]] ~ grp[[n]])
  cat(n, "DA-score t-test p =", signif(tt$p.value, 3), "\n")   # expect 0.002, 6e-5, 0.006, 0.004, 1.2e-4
}

de_raw <- lapply(cohort_names, function(n) run_limma(mats[[n]], grp[[n]]))
de_adj <- lapply(cohort_names, function(n) run_limma(mats[[n]], grp[[n]], scores[[n]]))
names(de_raw) <- names(de_adj) <- cohort_names
print(data.frame(cohort = cohort_names,
                 unadjusted  = sapply(de_raw, function(d) sum(d$adj.P.Val < 0.05)),
                 DA_adjusted = sapply(de_adj, function(d) sum(d$adj.P.Val < 0.05))))

common <- Reduce(intersect, lapply(de_raw, rownames))
cat("Genes common to all 5 cohorts:", length(common), "\n")          # expect 12154
res <- list(unadjusted = meta_rra(de_raw), DA_adjusted = meta_rra(de_adj))
print(sapply(res, function(m) c(up = length(sigset(m, "up")), down = length(sigset(m, "dn")))))
# expect unadjusted 45 up / 406 down; DA_adjusted 0 / 0

stab_raw <- stability(de_raw, res$unadjusted)
stab_adj <- stability(de_adj, res$DA_adjusted)

write.csv(res$unadjusted$up,   "results/RRA_up_unadjusted.csv",   row.names = FALSE)
write.csv(res$unadjusted$dn,   "results/RRA_down_unadjusted.csv", row.names = FALSE)
write.csv(res$DA_adjusted$up,  "results/RRA_up_DAadjusted.csv",   row.names = FALSE)
write.csv(res$DA_adjusted$dn,  "results/RRA_down_DAadjusted.csv", row.names = FALSE)
write.csv(stab_raw, "results/stability_unadjusted.csv", row.names = FALSE)
write.csv(stab_adj, "results/stability_DAadjusted.csv", row.names = FALSE)
saveRDS(list(de_raw = de_raw, de_adj = de_adj, res = res, scores = scores, common = common,
             stab_raw = stab_raw, stab_adj = stab_adj), "data/v2_results.rds")
