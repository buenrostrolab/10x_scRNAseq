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

a <- getvalues(yhMake("fresh_68k_pbmc_donor_a"))
ft <- list(row="INTEGER", column="INTEGER", value="INTEGER")
f1name <- "dat.sqlite"
f2name <- "dat.csv"
n <- 10000

sizes <- sapply(1:10, function(i){
  print(i)
  s <- data.frame(summary(assay(a,1)[,sample(1:dim(a)[2],min(dim(a)[2],n))]))
  
  #SQLite
  db <- dbConnect(SQLite(), dbname=f1name)
  dbWriteTable(conn=db, name=paste0("data", i), value=s, field.types=ft)
  dbDisconnect(db)
  o1 <- file.info(f1name)$size
  
  # csv
  write.table(s, file = "dat.csv", append = TRUE, row.names = FALSE, quote = FALSE, sep = ",", col.names = FALSE)
  
  # gzip
  system("cp dat.csv dat2.csv")
  if(file.exists("dat2.csv.gz")) system("rm dat2.csv.gz")
  system("gzip dat2.csv")
  o2 <- file.info("dat2.csv.gz")$size
  o3 <- file.info("dat.csv")$size
  return(c(o1,o2,o3))
})

df <- data.frame(t(rbind(1:10*n, sizes/1048576)))
names(df) <- c("ncells", "sqlite", "gzip", "uncompressed")

saveRDS(df,"../data/fileSizes.rds")
