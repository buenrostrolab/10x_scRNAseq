#detach("package:yolo", unload=TRUE)
library(yolo)
library(RSQLite)
library(SummarizedExperiment)
library(ggplot2)
library(reshape2)

if (basename(getwd()) != "code") setwd("code")
sqlite <- "../data/tenX.sqlite"

if(!file.exists("../data/memUsage.rds")){
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
  saveRDS(df, "../data/memUsage.rds")
} else {
  df <- readRDS("../data/memUsage.rds")
}

mdf <- melt(df, id.vars = 3)

# Project Data
nv <- c("ncells", "variable", "value")
my <- lm(yolo ~ ncells, df)
mr <- lm(regular ~ ncells, df)
vals <- c(100000,500000,1000000,5000000)
#yolopred <- predict(my, data.frame(ncells = vals))
regpred <- predict(mr, data.frame(ncells = vals))
#dfappend <- rbind(setNames(data.frame(vals, "regular", regpred), nv), setNames(data.frame(vals, "yolo", yolopred),nv))
dfappend <- setNames(data.frame(vals, "regular", regpred), nv)
                  
ggplot(rbind(mdf,dfappend), aes(x=ncells, y=value, color = variable)) + geom_line() + theme_bw()+ labs(title = "scRNA-Seq Data Size", 
    x = "Number of Cells", y = "RAM Requirements (Megabytes)") +
  scale_x_continuous(trans = "log2", breaks = c(1000, 10000, 100000, 1000000) ) + 
  scale_y_continuous(trans = "log2", breaks =  c(100, 1000, 10000, 10000)) +
  geom_point(data=mdf,aes(x=ncells,y=value, colour = variable))+labs(colour = "Data Type") +
  #geom_segment(aes(x = 100000, y = 7000, xend = 1000000, yend = 6312), arrow = arrow(length = unit(0.05, "npc")), colour = "black") +
  annotate("text", x = 300000, y = 6000, label = "~ 6.3 Gb", size = 7) +
  geom_point(data=data.frame(x = 1000000, y = 6312),aes(x,y),colour="black",size=2)

