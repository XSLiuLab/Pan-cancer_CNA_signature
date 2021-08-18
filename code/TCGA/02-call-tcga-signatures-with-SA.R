library(sigminer)

tally_X <- readRDS("data/TCGA/tcga_cn_tally_X.rds")
tally_X_noLOH <- readRDS("data/TCGA/tcga_cn_tally_X_noLOH.rds")

dim(tally_X$nmf_matrix)
dim(tally_X_noLOH$simplified_matrix)

sigs <- sig_auto_extract(
  tally_X$nmf_matrix,
  result_prefix = "TCGA_CN176",
  destdir = "SA",
  K0 = 40,
  nrun = 100,
  strategy = "stable",
  cores = 16,
  optimize = TRUE,
  skip = TRUE)
saveRDS(sigs, file = "data/TCGA/tcga_cn_sigs_CN176_SA.rds")

sigs2 <- sig_auto_extract(
  tally_X_noLOH$simplified_matrix,
  result_prefix = "TCGA_CN136",
  destdir = "SA",
  K0 = 40,
  nrun = 100,
  strategy = "stable",
  cores = 16,
  optimize = TRUE,
  skip = TRUE)
saveRDS(sigs2, file = "data/TCGA/tcga_cn_sigs_CN136_SA.rds")
