# 06_rnaseq_check.R
# Independent RNA-seq check in GSE136666 (substantia nigra, 5 control vs 5 PD; each sample is a pool of 2 donors).
# NOTE: reconstructed from the described methods; compare its output with the paper (165 significant genes,
# 22,150 tested, Spearman rho 0.54 over 10,305 shared genes, 40/41 up and 361/382 down in the same direction).
suppressPackageStartupMessages({ library(GEOquery); library(DESeq2); library(EnsDb.Hsapiens.v86); library(AnnotationDbi) })
source("scripts/00_helpers.R")
options(timeout = 1200)
getGEOSuppFiles("GSE136666", baseDir = "data")
untar("data/GSE136666/GSE136666_RAW.tar", exdir = "data/GSE136666")
ph <- pData(getGEO("GSE136666", GSEMatrix = TRUE)[[1]])
ph <- ph[ph$`tissue:ch1` == "Substantia nigra", ]
ph$group <- factor(ifelse(grepl("control", ph$`disease state:ch1`, ignore.case = TRUE), "Control", "PD"), levels = c("Control", "PD"))
stopifnot(identical(as.integer(table(ph$group)), c(5L, 5L)))

files <- list.files("data/GSE136666", pattern = "Counts_SN_.*\\.txt\\.gz$", full.names = TRUE)
names(files) <- sub("_Counts.*", "", basename(files))
files <- files[rownames(ph)]
cnt <- lapply(files, function(f) { x <- read.delim(f, header = TRUE, row.names = 1, check.names = FALSE); x[, 1, drop = FALSE] })
ids <- Reduce(intersect, lapply(cnt, rownames))
counts <- round(sapply(cnt, function(x) x[ids, 1])); rownames(counts) <- ids
rownames(counts) <- sub("\\..*$", "", rownames(counts))     # strip Ensembl version suffix if present

sym <- mapIds(EnsDb.Hsapiens.v86, keys = rownames(counts), column = "SYMBOL", keytype = "GENEID", multiVals = "first")
counts <- counts[!is.na(sym), ]; rownames(counts) <- NULL
sym <- sym[!is.na(sym)]
counts <- rowsum(counts, sym)                                # sum duplicate symbols
keep <- rowSums(counts >= 10) >= 3
cat("Genes tested:", sum(keep), "\n")                        # expect 22150

dds <- DESeqDataSetFromMatrix(counts[keep, ], data.frame(group = ph$group, row.names = colnames(counts)), ~ group)
rs  <- as.data.frame(results(DESeq(dds), name = "group_PD_vs_Control"))
cat("Significant genes (padj < 0.05):", sum(rs$padj < 0.05, na.rm = TRUE), "\n")   # expect 165

v2 <- readRDS("data/v2_results.rds"); de_raw <- v2$de_raw; res <- v2$res; common <- v2$common
mean_t <- rowMeans(sapply(de_raw, function(d) setNames(d[common, "t"], common)))
sh <- intersect(rownames(rs), names(mean_t))
cat("Shared genes:", length(sh), "\n")                       # expect 10305
rho <- cor.test(rs[sh, "stat"], mean_t[sh], method = "spearman")
cat("Spearman rho:", round(rho$estimate, 2), "\n")           # expect 0.54
up_g <- intersect(sigset(res$unadjusted, "up"), sh); dn_g <- intersect(sigset(res$unadjusted, "dn"), sh)
cat("Up genes same direction:", sum(rs[up_g, "stat"] > 0, na.rm = TRUE), "/", length(up_g), "\n")
cat("Down genes same direction:", sum(rs[dn_g, "stat"] < 0, na.rm = TRUE), "/", length(dn_g), "\n")
write.csv(cbind(gene = rownames(rs), rs), "results/RNAseq_GSE136666_DESeq2.csv", row.names = FALSE)
saveRDS(list(rs = rs, mean_t = mean_t, shared = sh), "data/RNAseq_GSE136666_results.rds")
