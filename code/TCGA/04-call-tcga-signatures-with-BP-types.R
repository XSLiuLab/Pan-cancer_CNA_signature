library(sigminer)

tally_X <- readRDS("tcga_cn_tally_X.rds")
tally_X_noLOH <- readRDS("tcga_cn_tally_X_noLOH.rds")

# sample_info <- readRDS("data/TCGA/tcga_cli.rds")
# TCGA_type_info <- sample_info %>% dplyr::select(sample, `cancer type abbreviation`) %>%
#   dplyr::filter(sample %in% rownames(tally_X$nmf_matrix))
# colnames(TCGA_type_info)[2] <- "cancer_type"
# saveRDS(TCGA_type_info, file = "data/TCGA/TCGA_type_info.rds")
tcga_types <- readRDS("TCGA_type_info.rds")

for (i in unique(tcga_types$cancer_type)) {

  message("handling type: ", i)

  fn <- file.path("BP", paste0("TCGA_CN176_", i, ".rds"))

  if (!file.exists(fn)) {
    type_samples <- tcga_types %>%
      dplyr::filter(cancer_type == i) %>%
      dplyr::pull(sample)

    message("Sample number: ", length(type_samples))

    sigs <- bp_extract_signatures(
      tally_X$nmf_matrix[type_samples, , drop = FALSE],
      range = 2:20,
      cores = 30,
      cores_solution = 4,
      cache_dir = "BP/BP_tcga_types",
      keep_cache = TRUE,
      only_core_stats = TRUE
    )
    saveRDS(sigs, file = fn)
  }
}

for (i in unique(tcga_types$cancer_type)) {

  message("handling type: ", i)

  fn <- file.path("BP", paste0("TCGA_CN136_", i, ".rds"))
  if (!file.exists(fn)) {
    type_samples <- tcga_types %>%
      dplyr::filter(cancer_type == i) %>%
      dplyr::pull(sample)

    message("Sample number: ", length(type_samples))

    sigs <- bp_extract_signatures(
      tally_X_noLOH$simplified_matrix[type_samples, , drop = FALSE],
      range = 2:20,
      cores = 30,
      cores_solution = 4,
      cache_dir = "BP/BP_tcga_types",
      keep_cache = TRUE,
      only_core_stats = TRUE
    )
    saveRDS(sigs, file = fn)
  }
}
