# 03_sensitivity.R
# Scenarios s0-s4: outlier removal and scan-month batch, each with/without the DA score.
# Also defines the 69-gene "core" (down genes significant in both s0 and s3).
source("scripts/00_helpers.R")
v2 <- readRDS("data/v2_results.rds")
de_raw <- v2$de_raw; common <- v2$common; res <- v2$res
cohort_res <- lapply(cohort_names, function(n) readRDS(file.path("data", paste0(n, "_results.rds"))))
names(cohort_res) <- cohort_names
grp  <- lapply(cohort_res, function(x) factor(x$group, levels = c("Control", "PD")))
mats <- lapply(cohort_names, function(n) collapse_probes(cohort_res[[n]]$rma, cohort_db[[n]]))
names(mats) <- cohort_names

# Scan month per sample (from CEL headers, saved in script 01)
batch <- lapply(cohort_names, function(n) {
  b <- factor(format(cohort_res[[n]]$scan_date, "%Y-%m")); names(b) <- names(cohort_res[[n]]$scan_date)
  b[colnames(mats[[n]])]
})
names(batch) <- cohort_names

# Four samples flagged as outliers (robust z > 2.5 on RLE IQR or PCA distance, see QC below)
drop_samples <- c("GSM184358.CEL.gz", "GSM508735.CEL.gz", "GSM1192714_BR56_4SN.CEL.gz", "GSM208635.cel.gz")

# QC flags (RLE IQR and PCA distance, robust z-scores)
for (n in cohort_names) {
  ex  <- mats[[n]]
  iqr <- apply(sweep(ex, 1, apply(ex, 1, median)), 2, IQR)
  z_iqr <- (iqr - median(iqr)) / mad(iqr)
  pc  <- prcomp(t(ex))$x[, 1:2]
  d   <- sqrt(rowSums(sweep(pc, 2, apply(pc, 2, median))^2))
  z_d <- (d - median(d)) / mad(d)
  f <- data.frame(sample = colnames(ex), group = grp[[n]], RLE_IQR_z = round(z_iqr, 1), PCA_dist_z = round(z_d, 1))
  cat("\n==", n, "==\n"); print(f[f$RLE_IQR_z > 2.5 | f$PCA_dist_z > 2.5, ])
}

# s0, s1, s2: all samples / outliers removed / outliers removed + GSE20292 batch term
run_scenario <- function(drop = character(0), batch_20292 = FALSE) {
  de_r <- list(); de_a <- list()
  for (n in cohort_names) {
    keep <- !(colnames(mats[[n]]) %in% drop)
    m  <- mats[[n]][, keep, drop = FALSE]
    df <- data.frame(g = droplevels(grp[[n]][keep]), s = da_score(m))
    blk <- batch_20292 && n == "GSE20292"
    if (blk) df$blk <- factor(ifelse(as.numeric(sub("GSM([0-9]+).*", "\\1", colnames(m))) <= 508722, "early", "late"))
    f_r <- if (blk) ~ g + blk else ~ g
    f_a <- if (blk) ~ g + blk + s else ~ g + s
    getde <- function(f) topTable(eBayes(lmFit(m, model.matrix(f, df))), coef = "gPD", number = Inf, sort.by = "none")
    de_r[[n]] <- getde(f_r); de_a[[n]] <- getde(f_a)
  }
  mr <- meta_rra(de_r); ma <- meta_rra(de_a)
  list(mr = mr, ma = ma, summary = data.frame(analysis = c("unadjusted", "DA_adjusted"),
       up = c(length(sigset(mr, "up")), length(sigset(ma, "up"))),
       down = c(length(sigset(mr, "dn")), length(sigset(ma, "dn")))))
}
s0 <- run_scenario()                                   # expect 45/406 and 0/0
s1 <- run_scenario(drop = drop_samples)                # expect 53/486 and 0/0
s2 <- run_scenario(drop = drop_samples, batch_20292 = TRUE)   # expect 75/522 and 0/0

# s3, s4: scan month as batch covariate (dropped if inseparable from group)
run_scenario2 <- function(drop = character(0)) {
  de_r <- list(); de_a <- list()
  for (n in cohort_names) {
    keep <- !(colnames(mats[[n]]) %in% drop)
    m  <- mats[[n]][, keep, drop = FALSE]
    df <- data.frame(g = droplevels(grp[[n]][keep]), s = da_score(m), b = droplevels(batch[[n]][keep]))
    mk <- function(f) {
      if (nlevels(df$b) < 2) f <- update(f, . ~ . - b)
      d <- model.matrix(f, df)
      if (!is.null(nonEstimable(d))) { message(n, ": batch not separable from group; batch term dropped"); d <- model.matrix(update(f, . ~ . - b), df) }
      d
    }
    getde <- function(f) topTable(eBayes(lmFit(m, mk(f))), coef = "gPD", number = Inf, sort.by = "none")
    de_r[[n]] <- getde(~ g + b); de_a[[n]] <- getde(~ g + b + s)
  }
  mr <- meta_rra(de_r); ma <- meta_rra(de_a)
  list(mr = mr, ma = ma, de_r = de_r, de_a = de_a, summary = data.frame(analysis = c("unadjusted", "DA_adjusted"),
       up = c(length(sigset(mr, "up")), length(sigset(ma, "up"))),
       down = c(length(sigset(mr, "dn")), length(sigset(ma, "dn")))))
}
s3 <- run_scenario2()                                  # expect 0 up / 72 down; adjusted 0/0
s4 <- run_scenario2(drop = drop_samples)               # expect 0 up / 81 down; adjusted 0/0
print(rbind(cbind(scenario = "s0", s0$summary), cbind(scenario = "s1", s1$summary), cbind(scenario = "s2", s2$summary),
            cbind(scenario = "s3", s3$summary), cbind(scenario = "s4", s4$summary)))

core <- intersect(sigset(s0$mr, "dn"), sigset(s3$mr, "dn"))
cat("Core genes:", length(core), "\n")                  # expect 69
write.csv(data.frame(gene = core), "results/core_genes.csv", row.names = FALSE)
saveRDS(list(s0 = s0, s1 = s1, s2 = s2, s3 = s3, s4 = s4, batch = batch, core = core,
             drop_samples = drop_samples), "data/sensitivity_results.rds")
