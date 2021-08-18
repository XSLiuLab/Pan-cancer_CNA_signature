dir.create("Xena")

# Phenotype (Specimen Centric) -----------------------------------------------

library(UCSCXenaTools)
pcawg_phenotype <- XenaData %>% dplyr::filter(XenaHostNames == "pcawgHub") %>%
  XenaScan("phenotype") %>% XenaScan("specimen centric") %>%
  XenaGenerate() %>%
  XenaQuery()
pcawg_phenotype <- pcawg_phenotype %>%
  XenaDownload(destdir = "Xena", trans_slash = TRUE)

phenotype_list <- XenaPrepare(pcawg_phenotype)
pcawg_samp_info_sp <- phenotype_list[1:4]

saveRDS(pcawg_samp_info_sp, file = "data/pcawg_samp_info_sp.rds")

# PCAWG (Specimen Centric) ------------------------------------------------

download.file("https://pcawg.xenahubs.net/download/20170119_final_consensus_copynumber_sp.gz",
              "Xena/pcawg_copynumber_sp.gz")

pcawg_cn <- data.table::fread("Xena/pcawg_copynumber_sp.gz")
pcawg_samp_info_sp <- readRDS("data/pcawg_samp_info_sp.rds")

sex_dt = pcawg_samp_info_sp$pcawg_donor_clinical_August2016_v9_sp %>%
  dplyr::select(xena_sample, donor_sex) %>%
  purrr::set_names(c("sample", "sex")) %>%
  data.table::as.data.table()

saveRDS(sex_dt, file = "data/pcawg_sex_sp.rds")

pcawg_cn <- pcawg_cn[!is.na(total_cn)]
pcawg_cn$value <- NULL
pcawg_cn <- pcawg_cn[, c(1:5, 7)]
colnames(pcawg_cn)[1:5] = c("sample", "Chromosome", "Start.bp", "End.bp", "modal_cn")

saveRDS(pcawg_cn, file = "data/pcawg_copynumber_sp.rds")

# LOH ---------------------------------------------------------------------

pcawg_cn <- readRDS("data/pcawg_copynumber_sp.rds")

pcawg_loh <- pcawg_cn %>%
  dplyr::filter(Chromosome %in% as.character(1:22)) %>%
  dplyr::mutate(len = End.bp - Start.bp + 1) %>%
  dplyr::group_by(sample) %>%
  dplyr::summarise(
    n_LOH = sum(minor_cn == 0 & modal_cn > 0 & len >= 1e4)
  ) %>%
  setNames(c("sample", "n_LOH")) %>%
  data.table::as.data.table()

saveRDS(pcawg_loh, file = "data/pcawg_loh.rds")
