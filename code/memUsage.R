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

# Make handle object without any particular column annotation
sample <- "fresh_68k_pbmc_donor_a"
con <- RSQLite::dbConnect(drv=RSQLite::SQLite(), dbname=sqlite)
rowmax <- as.numeric(RSQLite::dbGetQuery(conn=con, statement= paste0('select max(row) from ', sample)))
columnmax <- as.numeric(RSQLite::dbGetQuery(conn=con, statement= paste0('select max(column) from ', sample)))
dbDisconnect(con)
cD <- DataFrame(name = paste0(sample, seq(1, columnmax)))

d <- yoloHandleMake(rD, cD, lookupFileName = sqlite, lookupTableName = sample)
a <- getvalues(d)

ns <- seq(1000,60000,1000)
mb_memusage <- sapply(ns, function(n){
  s <- sample(1:dim(a)[2],n)
  o1 <- as.numeric(strsplit(format(object.size(d[,s]), units = "Mb"), split = " ")[[1]][1])
  o2 <- as.numeric(strsplit(format(object.size(a[,s]), units = "Mb"), split = " ")[[1]][1]) + as.numeric(strsplit(format(object.size(assay(a[,s],1)), units = "Mb"), split = " ")[[1]][1])
  c(o1,o2)
})

df <- data.frame(t(rbind(mb_memusage, ns)))
names(df) <- c("yolo", "regular", "ncells")
mdf <- melt(df, id.vars = 3)
ggplot(mdf, aes(x=ncells, y=value, color = variable)) + geom_line() + theme_bw()+ labs(title = "scRNA-Seq Data Size", 
    x = "Number of Cells", y = "Memory Use (Mb)") + geom_point(data=mdf,aes(x=ncells,y=value, colour = variable))+labs(colour = "Data Type")

