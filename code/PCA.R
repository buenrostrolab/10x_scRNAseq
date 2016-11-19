library(Matrix)
library(irlba)
library(yolo)
library(RSQLite)
library(SummarizedExperiment)
library(Rtsne)

if (basename(getwd()) != "code") setwd("code")
sqlite <- "../data/tenX.sqlite"

# Import rowspace
rD <- read.table("../data/genes.tsv")
names(rD) <- c("ensembl", "geneid")

# Make handle object without any particular column annotation
sample <- "aml035_post_transplant"
con <- RSQLite::dbConnect(drv=RSQLite::SQLite(), dbname=sqlite)
rowmax <- as.numeric(RSQLite::dbGetQuery(conn=con, statement= paste0('select max(row) from ', sample)))
columnmax <- as.numeric(RSQLite::dbGetQuery(conn=con, statement= paste0('select max(column) from ', sample)))
dbDisconnect(con)
cD <- DataFrame(name = paste0(sample, seq(1, columnmax)))

d <- yoloHandleMake(rD, cD, lookupFileName = sqlite, lookupTableName = sample)
a <- assay(getvalues(d),1)

n_comp <- 3
system.time({
  xt.x <- crossprod(a)
  x.means <- Matrix::colMeans(a)
  m <- dim(a)[1]
  xt.x <- (xt.x - m * tcrossprod(x.means)) / (m-1)
  svd.0 <- irlba(xt.x, nu=0, nv=n_comp, tol=1e-10)
})

pr <- prcomp_irlba(a, n = 10)
tsne <- Rtsne(pr$rotation, pca = FALSE)
plot(tsne$Y)

e <- predict(pr)


# Try large vis

library(largeVis)

neighbors <- randomProjectionTreeSearch(a)
edges <- buildEdgeMatrix(data = a, neighbors = neighbors)
rm(neighbors)
gc()
wij <- buildWijMatrix(edges)
rm(edges)
gc()
coords <- projectKNNs(wij)

cd <- data.frame(scale(t(coords)))
names(cd) <- c("x", "y")

ggplot(cd, aes( x = x,  y = y)) +
  geom_point() + theme_bw()
