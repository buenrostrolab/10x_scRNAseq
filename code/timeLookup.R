#detach("package:yolo", unload=TRUE)
library(yolo)
library(RSQLite)
library(SummarizedExperiment)
library(ggplot2)
library(reshape2)

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

time <- function(i, ns, d){
  print(i)
  bb <- sapply(ns, function(n){
    s <- sample(1:dim(d)[2],min(dim(d)[2],n))
    m <- d[,s]
    ptm <- proc.time()
    o1 <- getvalues(m)
    as.numeric(proc.time() - ptm)[2]
  })
  return(bb)
}

aml <- sapply(1:100, time, seq(800,4000,800), yhMake("aml027_post_transplant"))
saveRDS(aml,"../data/timeLookup/aml027_post_transplant_times.rds")

#fro <- sapply(1:100, time, seq(400,2000,400), yhMake("frozen_bmmc_healthy_donor1"))
#saveRDS(fro,"../data/timeLookup/frozen_bmmc_healthy_donor1_times.rds")

#bcells <- sapply(1:100, time, seq(2100,10100,2000), yhMake("b_cells"))
#saveRDS(bcells,"../data/timeLookup/b_cells_times.rds")

#froB <- sapply(1:100, time, seq(1550,8000,1560), yhMake("frozen_pbmc_donor_b"))
#saveRDS(froB, "../data/timeLookup/frozen_pbmc_donor_b_times.rds")

#a35po <- sapply(1:100, time, seq(180,1000,185), yhMake("aml035_post_transplant"))
#saveRDS(a35po ,"../data/timeLookup/aml035_post_transplant_times.rds")




