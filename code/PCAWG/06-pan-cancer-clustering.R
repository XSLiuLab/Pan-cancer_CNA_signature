library(cola)
library(tidyverse)

act <- readRDS("data/pcawg_cn_sigs_CN176_activity.rds")

df <- purrr::reduce(
  list(
    act$absolute,
    act$similarity
  ),
  dplyr::left_join,
  by = "sample"
)

mat <- df %>%
  filter(similarity > 0.75) %>%
  select(sample, starts_with("CNS")) %>%
  column_to_rownames("sample") %>%
  t()
mat_adj <- adjust_matrix(mat)
#1 rows have been removed with too low variance (sd < 0.05 quantile)

rownames(mat_adj)

# Select suitable parameters ----------------------------------------------

ds <- colSums(mat_adj)
boxplot(ds)

set.seed(123)
select_samps <- sample(names(sort(ds)), 500)

boxplot(ds[select_samps])

rl_samp <- run_all_consensus_partition_methods(mat_adj[, select_samps], top_n = 13, mc.cores = 8, max_k = 10)
cola_report(rl_samp, output_dir = "output/cola_report/pcawg_sigs_500_sampls", mc.cores = 8)

rm(rl_samp)

# Run clustering ----------------------------------------------------------

final <- run_all_consensus_partition_methods(
  mat_adj,
  top_value_method = "ATC",
  partition_method = "skmeans",
  top_n = 13, mc.cores = 8, max_k = 10
)
cola_report(final, output_dir = "output/cola_report/pcawg_sigs_all_sampls", mc.cores = 8)

saveRDS(final, file = "data/pcawg_cola_result.rds")

# res = final["ATC", "skmeans"]
# select_partition_number(res)
# collect_classes(res)
