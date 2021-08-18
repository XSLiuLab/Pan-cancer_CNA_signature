library(sigminer)

tally_X <- readRDS("pcawg_cn_tally_X.rds")
tally_X_noLOH <- readRDS("pcawg_cn_tally_X_noLOH.rds")

#sample_info <- readRDS("~/projects/CNSigs/data/pcawg_tidy_tumor_sample_info.rds")
#pcawg_type_info <- sample_info %>% dplyr::select(sample, cancer_type)
#saveRDS(pcawg_type_info, file = "data/pcawg_type_info.rds")
pcawg_types <- readRDS("pcawg_type_info.rds")

for (i in unique(pcawg_types$cancer_type)) {
  message("handling type: ", i)

  fn <- file.path("BP", paste0("PCAWG_CN176_", gsub("/", "-", i), ".rds"))

  if (!file.exists(fn)) {
    type_samples <- pcawg_types %>%
      dplyr::filter(cancer_type == i) %>%
      dplyr::pull(sample)
    message("Sample number: ", length(type_samples))

    sigs <- bp_extract_signatures(
      tally_X$nmf_matrix[type_samples, , drop = FALSE],
      range = 2:min(20, length(type_samples)-1),
      cores = 30,
      cores_solution = 20,
      cache_dir = "BP/BP_PCAWG_types",
      keep_cache = TRUE,
      only_core_stats = TRUE
    )
    saveRDS(sigs, file = fn)
  }
}

for (i in unique(pcawg_types$cancer_type)) {

  message("handling type: ", i)

  fn <- file.path("BP", paste0("PCAWG_CN136_", gsub("/", "-", i), ".rds"))

  if (!file.exists(fn)) {
    type_samples <- pcawg_types %>%
      dplyr::filter(cancer_type == i) %>%
      dplyr::pull(sample)

    message("Sample number: ", length(type_samples))

    sigs <- bp_extract_signatures(
      tally_X_noLOH$simplified_matrix[type_samples, , drop = FALSE],
      range = 2:min(20, length(type_samples)-1),
      cores = 30,
      cores_solution = 20,
      cache_dir = "BP/BP_PCAWG_types",
      keep_cache = TRUE,
      only_core_stats = TRUE
    )
    saveRDS(sigs, file = fn)
  }
}
