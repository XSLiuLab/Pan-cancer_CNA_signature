library(tidyverse)

driver_info <- read_tsv("data/Compendium_Cancer_Genes.tsv")
driver_info %>%
  dplyr::select(SYMBOL, ROLE) %>%
  unique() %>%
  dplyr::filter(! ROLE %in% c(NA, "ambiguous")) %>%
  dplyr::pull(ROLE) -> x2

table(x2)

driver_info <- read_tsv("data/Compendium_Cancer_Genes.tsv")
driver_info %>%
  dplyr::select(SYMBOL, ROLE) %>%
  unique() %>%
  dplyr::filter(! ROLE %in% c(NA, "ambiguous")) -> driver_info

download.file("https://xena-pcawg.hiplot.com.cn/download/pcawg_whitelist_coding_drivers_v1_sep302016.sp.xena",
              "Xena/pcawg_coding_drivers.txt")

pcawg_drivers <- read_tsv("Xena/pcawg_coding_drivers.txt")

pcawg_driver_number <- pcawg_drivers %>%
  dplyr::count(sample)

pcawg_driver_gene <- pcawg_drivers %>%
  dplyr::left_join(driver_info, by = c("gene" = "SYMBOL")) %>%
  dplyr::filter(!is.na(ROLE)) %>%
  dplyr::select(sample, gene, ROLE) %>%
  unique()

pcawg_driver_gene

save(pcawg_driver_number, pcawg_driver_gene, file = "data/pcawg_driver_info.RData")
