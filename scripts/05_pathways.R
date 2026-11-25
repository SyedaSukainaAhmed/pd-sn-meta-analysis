# 05_pathways.R
# GO biological-process enrichment (clusterProfiler, unsimplified) and gene- vs pathway-level
# reproducibility across cohorts (each cohort = its 400 most down-regulated genes).
suppressPackageStartupMessages({ library(clusterProfiler); library(org.Hs.eg.db) })
source("scripts/00_helpers.R")
v2 <- readRDS("data/v2_results.rds"); de_raw <- v2$de_raw; res <- v2$res; common <- v2$common
meta_dn <- sigset(res$unadjusted, "dn")

ego_run <- function(genes) enrichGO(genes, OrgDb = org.Hs.eg.db, keyType = "SYMBOL", ont = "BP",
                                    universe = common, pAdjustMethod = "BH", pvalueCutoff = 0.05,
                                    qvalueCutoff = 0.05, minGSSize = 10, maxGSSize = 500, readable = FALSE)
top_dn <- lapply(de_raw, function(d) { t <- setNames(d[common, "t"], common); names(sort(t))[1:400] })

ego_meta   <- ego_run(meta_dn)
meta_terms <- ego_meta@result$ID[ego_meta@result$p.adjust < 0.05]
cat("Enriched terms among meta-down genes:", length(meta_terms), "\n")      # expect 105
ego_coh    <- lapply(top_dn, ego_run)
coh_terms  <- lapply(ego_coh, function(e) if (is.null(e)) character(0) else e@result$ID[e@result$p.adjust < 0.05])
print(sapply(coh_terms, length))                                              # expect 0, 92, 58, 105, 29

gene_n <- sapply(meta_dn,   function(g)  sum(sapply(top_dn,    function(v) g  %in% v)))
term_n <- sapply(meta_terms, function(tm) sum(sapply(coh_terms, function(v) tm %in% v)))
cmp <- data.frame(level = c("genes", "pathways"), n_items = c(length(gene_n), length(term_n)),
  in_ge2_cohorts = c(mean(gene_n >= 2), mean(term_n >= 2)),
  in_ge3_cohorts = c(mean(gene_n >= 3), mean(term_n >= 3)),
  in_all_5       = c(mean(gene_n == 5), mean(term_n == 5)))
cmp[, 3:5] <- round(cmp[, 3:5], 3); print(cmp)
write.csv(cmp, "results/pathway_reproducibility.csv", row.names = FALSE)
write.csv(ego_meta@result[ego_meta@result$p.adjust < 0.05, c("ID", "Description", "Count", "p.adjust")],
          "results/GO_meta_down_terms.csv", row.names = FALSE)
saveRDS(list(ego_meta = ego_meta, ego_coh = ego_coh, cmp = cmp, gene_n = gene_n, term_n = term_n),
        "data/pathway_unsimplified.rds")
