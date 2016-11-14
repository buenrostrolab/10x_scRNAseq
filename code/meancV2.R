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

rowVar <- function(x) rowSums((data.matrix(x) - Matrix::rowMeans(x))^2)/(dim(x)[2] - 1)

a <- assay(getvalues(yhMake("aml035_post_transplant")), 1)
true_mean <- Matrix::rowMeans(a)
true_cv <- sqrt(rowVar(a))/true_mean

mean_cv_error <- function(i, n){
    s <- sample(1:dim(d)[2],min(dim(d)[2],n))
    m <- a[,s]
    mean <- Matrix::rowMeans(m)
    cv <- sqrt(rowVar(m))/true_mean
    return(c(mean(mean-true_mean), sum(cv[complete.cases(cv)] - true_cv[complete.cases(cv)])/length(complete.cases(cv))))
}

ntotal <- 909
pdf <- setNames(data.frame(t(sapply(seq(101,ntotal,101), function(n){
  return(round(c(n/ntotal, rowMeans(abs(sapply(1:3, mean_cv_error, n)))),3))
}))), c("n", "mean", "cv"))


