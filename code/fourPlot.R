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

#samples <- c("aml027_post_transplant","aml027_pre_transplant")
samples <- c("cd14_monocytes",
  "jurkat_293t_99_1", "jurkat_293t_50_50", "jurkat",
  "t293")

oall <- Reduce('+', lapply(samples, yhMake))
class <- as.character(colData(oall)$lookupTableName)
aa <- assay(getvalues(oall), 1)
print(dim(aa))
ab <- sweep(aa, 2, Matrix::colSums(aa), FUN="/") * 1000000
ac <- sparseMatrix(i=summary(ab)$i, j=summary(ab)$j, x=log2(summary(ab)$x + 1))

pr <- prcomp_irlba(ac, n = 10)
X <- abs(pr$rotation)

library(largeVis)

system.time({
  neighbors <- randomProjectionTreeSearch(t(X))
  edges <- buildEdgeMatrix(data = t(X), neighbors = neighbors)
  rm(neighbors)
  gc()
  wij <- buildWijMatrix(edges)
  rm(edges)
  gc()
})
system.time({
  coords <- projectKNNs(wij)
})

cd <- data.frame(scale(t(coords)), class)
names(cd) <- c("x", "y", "class")

saveRDS(cd, file = "../data/fourPlot.rds")

cd <- readRDS("../data/fourPlot.rds")

ggplot(cd[sample(nrow(cd)),], aes( x = x,  y = y, color = class)) +
  geom_point() + theme_bw()


