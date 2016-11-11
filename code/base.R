#detach("package:yolo", unload=TRUE)
library(yolo)
library(RSQLite)
library(SummarizedExperiment)

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
str(d)

a <- getvalues(d)@assays[[1]]


## Do stuff with the mean
dim()

# Get Random samples
nSamples <- 100
randomSamples <- sort(sample(1:columnmax, 100))

truemean <- Matrix::rowMeans(getvalues(d[,])@assays[[1]])
