library(yolo)
library(RSQLite)
library(SummarizedExperiment)
library(ggplot2)
library(reshape2)
library(irlba)

if (basename(getwd()) != "code") setwd("code")
sqlite <- "../data/tenX.sqlite"

# Import rowspace
rD <- read.table("../data/genes.tsv")
names(rD) <- c("ensembl", "geneid")

yhMake <- function(sample){
  con <- RSQLite::dbConnect(drv=RSQLite::SQLite(), dbname=sqlite)
  rowmax <- as.numeric(RSQLite::dbGetQuery(conn=con, statement= paste0('select max(row) from ', sample)))
  columnmax <- as.numeric(RSQLite::dbGetQuery(conn=con, statement= paste0('select max(column) from ', sample)))
  dbDisconnect(con)
  cD <- DataFrame(name = paste0(sample, seq(1, columnmax)))
  
  d <- yoloHandleMake(rD, cD, lookupFileName = sqlite, lookupTableName = sample)
  return(d)
}

#samples <- c("aml027_post_transplant","aml027_pre_transplant",
#  "aml035_post_transplant","aml035_pre_transplant","b_cells","cd4_t_helper",
#  "cd14_monocytes","cd34","cd56_nk","cytotoxic_t","fresh_68k_pbmc_donor_a",
#  "frozen_bmmc_healthy_donor1","frozen_bmmc_healthy_donor2","frozen_pbmc_b_c_50_50",
#  "frozen_pbmc_b_c_90_10","frozen_pbmc_b_c_99_1","frozen_pbmc_donor_a","frozen_pbmc_donor_b",
#  "frozen_pbmc_donor_c", "jurkat_293t_99_1", "jurkat_293t_50_50", "jurkat","memory_t",
#  "naive_cytotoxic","naive_t","regulatory_t","t293", "pbmc3k", "pbmc6k", "pbmc33k")

samples <- c("aml027_pre_transplant")

x <- assay(getvalues(Reduce('+', lapply(samples, yhMake))), 1)

rowvar <-  rowSums((x - Matrix::rowMeans(x))^2)/(dim(x)[2] - 1)

saveRDS(cd, file = "fullVariance.rds")