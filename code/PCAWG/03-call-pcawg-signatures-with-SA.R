library(sigminer)

tally_X <- readRDS("data/pcawg_cn_tally_X.rds")
tally_X_noLOH <- readRDS("data/pcawg_cn_tally_X_noLOH.rds")

dim(tally_X$nmf_matrix)

tcga_cn_sigs_CN176_SA <- readRDS("~/projects/CNX/data/TCGA/tcga_cn_sigs_CN176_SA.rds")
tcga_cn_sigs_CN136_SA <- readRDS("~/projects/CNX/data/TCGA/tcga_cn_sigs_CN136_SA.rds")

sigs <- sig_auto_extract(
  tally_X$nmf_matrix,
  result_prefix = "CN176",
  destdir = "SA",
  K0 = 40,
  nrun = 100,
  ref_sigs = tcga_cn_sigs_CN176_SA,
  strategy = "ms",
  cores = 20,
  optimize = TRUE,
  skip = TRUE)
saveRDS(sigs, file = "data/pcawg_cn_sigs_CN176_SA.rds")

sim_v1 <- get_sig_similarity(sigs, tcga_cn_sigs_CN176_SA)
pheatmap::pheatmap(sim_v1$similarity, cluster_cols = F, cluster_rows = F, display_numbers = T)

sigs2 <- sig_auto_extract(
  tally_X_noLOH$simplified_matrix,
  result_prefix = "CN136",
  destdir = "SA",
  K0 = 40,
  nrun = 100,
  ref_sigs = tcga_cn_sigs_CN136_SA,
  strategy = "ms",
  cores = 16,
  optimize = TRUE,
  skip = TRUE)
saveRDS(sigs2, file = "data/pcawg_cn_sigs_CN136_SA.rds")

sim_v2 <- get_sig_similarity(sigs2, tcga_cn_sigs_CN136_SA)
pheatmap::pheatmap(sim_v2$similarity, cluster_cols = F, cluster_rows = F, display_numbers = T)

show_sig_profile(sigs, style = "cosmic", mode = "copynumber", method = "X")
show_sig_profile(sigs2, style = "cosmic", mode = "copynumber", method = "X")

# sp = sigprofiler_import("SP/PCAWG_CN176X/", order_by_expo = TRUE, type = "all")
# sim <- purrr::map(sp$solution_list, function(x) {
#   sim <- get_sig_similarity(x, sigs)
#   sim <- sim$similarity
#   y <- diag(sim)
#   names(y) <- colnames(sim)[seq_along(y)]
#   y
# })
# sort(sapply(sim, mean))
#
# sp2 = sigprofiler_import("SP/PCAWG_CN136//", order_by_expo = TRUE, type = "all")
# sim2 <- purrr::map(sp2$solution_list, function(x) {
#   sim <- get_sig_similarity(x, sigs2)
#   sim <- sim$similarity
#   y <- diag(sim)
#   names(y) <- colnames(sim)[seq_along(y)]
#   y
# })
# sort(sapply(sim2, mean))
#
