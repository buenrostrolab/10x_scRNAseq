library(cellrangerRkit)
library(RSQLite)

d <- "/Users/lareauc"
f1name <- "tenX.sqlite"
ft <- list(row="INTEGER", column="INTEGER", value="INTEGER")

#samples <- c("aml027_post_transplant","aml027_pre_transplant",
#  "aml035_post_transplant","aml035_pre_transplant","b_cells","cd4_t_helper",
#  "cd14_monocytes","cd34","cd56_nk","cytotoxic_t","ercc","fresh_68k_pbmc_donor_a",
#  "frozen_bmmc_healthy_donor1","frozen_bmmc_healthy_donor2","frozen_pbmc_b_c_50_50",
#  "frozen_pbmc_b_c_90_10","frozen_pbmc_b_c_99_1","frozen_pbmc_donor_a","frozen_pbmc_donor_b",
#  "frozen_pbmc_donor_c")
#samples <- c("jurkat_293t_99_1","jurkat","memory_t",
#  "naive_cytotoxic","naive_t","regulatory_t","293t")

samples <- c("jurkat_293t_50")
for(s in samples){
  download_sample(sample_name=s, sample_dir=d, host="https://s3-us-west-2.amazonaws.com/10x.files/samples/cell/")
  print(s)
  mat <- load_cellranger_matrix(d)
  if(dim(mat)[1] == 32738){ #removed ercc sample
    print(dim(mat))
    if(s == "293t") s <- "t293"
    db <- dbConnect(SQLite(), dbname=f1name)
    df <- data.frame(summary(exprs(mat)))
    if(max(df[,1]) < 32738){ # makes correct dimension in case of lots of zeros
      df <- rbind(df, c(32738, 1, 0))
    }
    names(df) <- c("row", "column", "value")
    dbWriteTable(conn=db, name=s, value=df, field.types=ft)
    dbDisconnect(db)
  }
}


# Check to see which samples are available
db <- dbConnect(SQLite(), dbname=f1name)
dbListTables(db)
length(dbListTables(db))
dbDisconnect(db)
