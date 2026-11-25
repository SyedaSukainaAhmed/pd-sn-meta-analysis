# 01_download_preprocess.R
# Downloads the five microarray cohorts from GEO, RMA-normalises each one separately,
# and saves normalised expression, group labels and scan dates to data/<GSE>_results.rds.
# Needs ~3 GB disk for the raw CEL files and ~30-60 min. Raw files are not in the repo.
# Run from the project folder:  setwd("<path>/pd-sn-meta-analysis")
suppressPackageStartupMessages({ library(GEOquery); library(oligo); library(affy); library(Biobase) })
options(timeout = 1200)
dir.create("data", showWarnings = FALSE)
gsm <- function(x) sub("\\.(CEL|cel)(\\.gz)?$", "", basename(x))
fetch <- function(acc) {
  getGEOSuppFiles(acc, baseDir = "data")
  untar(file.path("data", acc, paste0(acc, "_RAW.tar")), exdir = file.path("data", acc))
  pData(getGEO(acc, GSEMatrix = TRUE)[[1]])
}
check <- function(g, nC, nPD, acc) {
  stopifnot(identical(as.integer(table(g)), c(as.integer(nC), as.integer(nPD))))
  message(acc, ": ", nC, " control / ", nPD, " PD - OK")
}
save_cohort <- function(acc, rma, group, dates) {
  names(dates) <- sampleNames(rma)
  saveRDS(list(rma = rma, group = group, scan_date = dates), file.path("data", paste0(acc, "_results.rds")))
}

## GSE7621 (GPL570, oligo): 9 control / 16 PD
ph <- fetch("GSE7621")
cel <- list.files("data/GSE7621", pattern = "\\.CEL\\.gz$", full.names = TRUE)
raw <- read.celfiles(cel); rma <- oligo::rma(raw)
lab <- ph$characteristics_ch1[match(gsm(sampleNames(rma)), ph$geo_accession)]
g <- factor(ifelse(grepl("control", lab, ignore.case = TRUE), "Control", "PD"), levels = c("Control", "PD"))
check(g, 9, 16, "GSE7621")
save_cohort("GSE7621", rma, g, as.Date(substr(pData(protocolData(raw))$dates, 1, 8), "%m/%d/%y"))

## GSE20163 (GPL96, affy): 9 control / 8 PD (controls have no diagnosis entry)
ph <- fetch("GSE20163")
raw <- ReadAffy(celfile.path = "data/GSE20163"); rma <- affy::rma(raw)
d  <- ph$`diagnosis:ch1`[match(gsm(sampleNames(rma)), ph$geo_accession)]
g  <- factor(ifelse(is.na(d), "Control", "PD"), levels = c("Control", "PD"))
check(g, 9, 8, "GSE20163")
save_cohort("GSE20163", rma, g, as.Date(substr(pData(protocolData(raw))$ScanDate, 1, 8), "%m/%d/%y"))

## GSE20292 (GPL96, affy): 15 control / 11 PD; GSM606624-26 excluded (not in the original dataset)
ph <- fetch("GSE20292")
excl <- c("GSM606624", "GSM606625", "GSM606626")
cel <- list.files("data/GSE20292", pattern = "\\.CEL\\.gz$", full.names = TRUE)
cel <- cel[!gsm(cel) %in% excl]
raw <- ReadAffy(filenames = cel); rma <- affy::rma(raw)
lab <- ph$`disease state:ch1`[match(gsm(sampleNames(rma)), ph$geo_accession)]
g <- factor(ifelse(grepl("control", lab, ignore.case = TRUE), "Control", "PD"), levels = c("Control", "PD"))
check(g, 15, 11, "GSE20292")
save_cohort("GSE20292", rma, g, as.Date(substr(pData(protocolData(raw))$ScanDate, 1, 8), "%m/%d/%y"))

## GSE49036 (GPL570, oligo): 8 Braak-0 controls vs 8 Braak 5-6 samples (stages 1-4 not used)
fetch("GSE49036")
cel_all <- list.files("data/GSE49036", pattern = "\\.CEL\\.gz$", full.names = TRUE)
fn  <- basename(cel_all)
grp <- ifelse(grepl("_C[0-9]", fn), "Control", ifelse(grepl("_BR56_", fn), "PD", "intermediate"))
cel <- cel_all[grp != "intermediate"]
g   <- factor(grp[grp != "intermediate"], levels = c("Control", "PD"))
check(g, 8, 8, "GSE49036")
raw <- read.celfiles(cel); rma <- oligo::rma(raw)
save_cohort("GSE49036", rma, g, as.Date(substr(pData(protocolData(raw))$dates, 1, 10)))

## GSE8397 (GPL96, affy): lateral substantia nigra only, 7 control / 9 PD
ph <- fetch("GSE8397")
ph <- ph[grepl("Lateral substantia nigra", ph$title), ]
cel_all <- list.files("data/GSE8397", pattern = "cel", ignore.case = TRUE, full.names = TRUE)
cel <- cel_all[gsm(cel_all) %in% ph$geo_accession]
raw <- ReadAffy(filenames = cel); rma <- affy::rma(raw)
t_  <- ph$title[match(gsm(sampleNames(rma)), ph$geo_accession)]
g   <- factor(ifelse(grepl("control", t_, ignore.case = TRUE), "Control", "PD"), levels = c("Control", "PD"))
check(g, 7, 9, "GSE8397")
save_cohort("GSE8397", rma, g, as.Date(substr(pData(protocolData(raw))$ScanDate, 1, 8), "%m/%d/%y"))
message("Done. Next: 02_main_analysis.R")
