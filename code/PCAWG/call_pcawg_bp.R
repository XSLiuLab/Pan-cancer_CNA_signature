library(sigminer)

tally_X <- readRDS("pcawg_cn_tally_X.rds")
tally_X_noLOH <- readRDS("pcawg_cn_tally_X_noLOH.rds")

sigs <- bp_extract_signatures(
  tally_X$nmf_matrix,
  range = 2:30,
  n_bootstrap = 10,
  n_nmf_run = 10,
  cores = 50,
  cache_dir = "pcawg_bp_pynmf",
  keep_cache = TRUE,
  cores_solution = 30,
  pynmf = TRUE,
  use_conda = TRUE
)

saveRDS(sigs, file = "pcawg_cn_sigs_CN176_BP.rds")

sigs2 <- bp_extract_signatures(
  tally_X_noLOH$simplified_matrix,
  range = 2:30,
  n_bootstrap = 10,
  n_nmf_run = 10,
  cores = 50,
  cache_dir = "pcawg_bp_pynmf",
  keep_cache = TRUE,
  cores_solution = 30,
  pynmf = TRUE,
  use_conda = TRUE
)
saveRDS(sigs2, file = "pcawg_cn_sigs_CN136_BP.rds")

