library(cellrangerRkit)

samples <- c("293t","aml027_post_transplant","aml027_pre_transplant",
  "aml035_post_transplant","aml035_pre_transplant","b_cells","cd4_t_helper",
  "cd14_monocytes","cd34","cd56_nk","cytotoxic_t","ercc","fresh_68k_pbmc_donor_a",
  "frozen_bmmc_healthy_donor1","frozen_bmmc_healthy_donor2","frozen_pbmc_b_c_50_50",
  "frozen_pbmc_b_c_90_10","frozen_pbmc_b_c_99_1","frozen_pbmc_donor_a","frozen_pbmc_donor_b",
  "frozen_pbmc_donor_c","jurkat_293t_50_50","jurkat_293t_99_1","jurkat","memory_t",
  "naive_cytotoxic","naive_t","regulatory_t")

d <- "/Users/lareauc/Desktop/"
for(s in samples){
  download_sample(sample_name=s,sample_dir=d, host="https://s3-us-west-2.amazonaws.com/10x.files/samples/cell/")
  print(s)
  assign(s, load_cellranger_matrix(d))
}
