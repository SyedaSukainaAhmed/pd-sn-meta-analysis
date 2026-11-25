# 00_helpers.R - shared functions. Sourced by scripts 02-07; nothing to run here.
suppressPackageStartupMessages({
  library(limma); library(AnnotationDbi); library(RobustRankAggreg); library(Biobase)
  library(hgu133a.db); library(hgu133plus2.db)
})

cohort_names <- c("GSE7621", "GSE20163", "GSE20292", "GSE49036", "GSE8397")
cohort_db <- list(GSE7621 = hgu133plus2.db, GSE20163 = hgu133a.db, GSE20292 = hgu133a.db,
                  GSE49036 = hgu133plus2.db, GSE8397 = hgu133a.db)
da_markers <- c("TH", "SLC6A3", "SLC18A2", "DDC", "NR4A2", "KCNJ6", "DRD2")

# One probe per gene, chosen without using the PD/control labels
collapse_probes <- function(eset, db) {
  ex  <- exprs(eset)
  map <- AnnotationDbi::select(db, keys = rownames(ex), columns = "SYMBOL", keytype = "PROBEID")
  map <- map[!is.na(map$SYMBOL), ]
  nsym <- tapply(map$SYMBOL, map$PROBEID, function(x) length(unique(x)))
  keep <- names(nsym)[nsym == 1]
  map  <- unique(map[map$PROBEID %in% keep & !grepl("_x_at$", map$PROBEID), ])
  map  <- map[map$PROBEID %in% rownames(ex), ]
  map$mean <- rowMeans(ex[map$PROBEID, , drop = FALSE])
  map  <- map[order(map$mean, decreasing = TRUE), ]
  map  <- map[!duplicated(map$SYMBOL), ]
  out  <- ex[map$PROBEID, , drop = FALSE]
  rownames(out) <- map$SYMBOL
  out
}

# Mean z-score of the marker genes, per sample
score_of <- function(mat, genes = da_markers) {
  m <- mat[intersect(genes, rownames(mat)), , drop = FALSE]
  colMeans(t(scale(t(m))))
}
da_score <- function(mat) score_of(mat, da_markers)

run_limma <- function(mat, grp, covar = NULL) {
  design <- if (is.null(covar)) model.matrix(~ grp) else model.matrix(~ grp + covar)
  topTable(eBayes(lmFit(mat, design)), coef = 2, number = Inf, sort.by = "none")
}

# Robust rank aggregation of per-cohort t-statistics (needs `common` in the calling environment)
meta_rra <- function(de_list, keys = names(de_list)) {
  up_lists <- lapply(de_list[keys], function(d) {
    t <- d[common, "t"]; names(t) <- common
    names(sort(t, decreasing = TRUE))
  })
  up <- aggregateRanks(up_lists, N = length(common))
  dn <- aggregateRanks(lapply(up_lists, rev), N = length(common))
  up$adj <- p.adjust(up$Score, "BH"); dn$adj <- p.adjust(dn$Score, "BH")
  list(up = up, dn = dn)
}
sigset <- function(m, dir) { x <- m[[dir]]; x$Name[x$adj < 0.05] }
topN   <- function(m, dir, n = 200) m[[dir]]$Name[order(m[[dir]]$Score)][1:n]

stability <- function(de_list, full, n = 200) {
  nm <- names(de_list); rows <- list()
  for (k in 1:2) for (drop in combn(nm, k, simplify = FALSE)) {
    m <- meta_rra(de_list, setdiff(nm, drop))
    rows[[length(rows) + 1]] <- data.frame(
      dropped = paste(drop, collapse = "+"), n_dropped = k,
      sig_up = length(sigset(m, "up")), sig_dn = length(sigset(m, "dn")),
      top200_up = length(intersect(topN(m, "up", n), topN(full, "up", n))),
      top200_dn = length(intersect(topN(m, "dn", n), topN(full, "dn", n))))
  }
  do.call(rbind, rows)
}
