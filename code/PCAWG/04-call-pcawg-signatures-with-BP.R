library(sigminer)

tally_X <- readRDS("data/pcawg_cn_tally_X.rds")
tally_X_noLOH <- readRDS("data/pcawg_cn_tally_X_noLOH.rds")

# library(profvis)
# debug(bp_extract_signatures)
#
# sigs <- bp_extract_signatures(
#   tally_X$nmf_matrix,
#   range = 2:30,
#   n_bootstrap = 5,
#   n_nmf_run = 5,
#   cores = 16,
#   cores_solution = 4,
#   one_batch = TRUE,
#   cache_dir = "BP/BP_PCAWG_test",
#   keep_cache = TRUE,
#   only_core_stats = TRUE
# )
#
# bp_show_survey(sigs, fixed_ratio = F, add_score = FALSE)
# saveRDS(sigs, file = "data/pcawg_cn_sigs_CN176_BP_test.rds")

library(NMF)
sigs <- sig_extract(tally_X$nmf_matrix, 10, nrun = 100, cores = 16)
saveRDS(sigs, file = "data/pcawg_cn_sigs_CN176_10sigs.rds")

sigs <- readRDS("data/pcawg_cn_sigs_CN176_10sigs.rds")
show_sig_profile(sigs, mode = "copynumber", method = "X", style = "cosmic")

sim <- get_sig_similarity(sigs, sigs)
pheatmap::pheatmap(sim$similarity, display_numbers = TRUE)
expo <- get_sig_exposure(sigs)
expo_rel <- get_sig_exposure(sigs, type = "relative")

show_cor(expo[, -1], cor_method = "pearson")
show_cor(expo_rel[, -1], cor_method = "pearson")
show_cor(expo[, -1])
show_cor(expo_rel[, -1])

plot(density(expo$Sig1))
plot(density(expo_rel$Sig1))

shapiro.test(expo$Sig1)
shapiro.test(expo_rel$Sig1)
# sigs <- bp_extract_signatures(
#   tally_X_noLOH$nmf_matrix,
#   range = 2:30,
#   n_nmf_run = 5,
#   cores = 16
# )
# saveRDS(sigs2, file = "data/pcawg_cn_sigs_CN136_BP_test.rds")

sigs <- sig_extract(tally_X$nmf_matrix, 11, nrun = 100, cores = 10)
saveRDS(sigs, file = "data/pcawg_cn_sigs_CN176_11sigs.rds")

# Compare BP result and Standard NMF --------------------------------------

library(sigminer)
solution1000 <- readRDS("BP/BP_PCAWG_1000_Extraction_Result.rds")
bp_show_survey(solution1000, add_score = F, fixed_ratio = F)

solution100 <- readRDS("data/pcawg_cn_sigs_CN176_BP.rds")
bp_show_survey(solution100, add_score = F, fixed_ratio = F)

sigs <- readRDS("data/pcawg_cn_sigs_CN176_11sigs.rds")

s1 <- get_sig_similarity(solution1000$object$K11, sigs)
s2 <- get_sig_similarity(solution100$object$K11, sigs)
s3 <- get_sig_similarity(solution1000$object$K11, solution100$object$K11)

ComplexHeatmap::pheatmap(s1$similarity, cluster_cols = F, cluster_rows = F)
ComplexHeatmap::pheatmap(s2$similarity, cluster_cols = F, cluster_rows = F)
ComplexHeatmap::pheatmap(s3$similarity, cluster_cols = F, cluster_rows = F)

mean(subset(solution1000$stats_sample, signature_number == 11)$cosine_distance_mean)
mean(subset(solution100$stats_sample, signature_number == 11)$cosine_distance_mean)
# Optimizing exposures ----------------------------------------------------

library(sigminer)
tally_X <- readRDS("data/pcawg_cn_tally_X.rds")
pcawg_types <- readRDS("data/pcawg_type_info.rds")

sc <- pcawg_types$cancer_type
names(sc) <- pcawg_types$sample

pcawg_expo <- bp_attribute_activity(
  solution1000$object$K11,
  sample_class = sc,
  nmf_matrix = tally_X$nmf_matrix,
  method = "bt",
  bt_use_prop = FALSE,
  return_class = "data.table",
  use_parallel = 12
)

saveRDS(pcawg_expo, file = "data/pcawg_cn_sigs_CN176_activity.rds")

# pcawg_expo <- bp_attribute_activity(
#   solution1000$object$K11,
#   sample_class = sc,
#   nmf_matrix = tally_X$nmf_matrix,
#   method = "bt",
#   bt_use_prop = TRUE,
#   return_class = "data.table",
#   use_parallel = 12
# )

pcawg_expo2 <- bp_attribute_activity(
  solution1000$object$K11,
  sample_class = sc,
  nmf_matrix = tally_X$nmf_matrix,
  method = "stepwise",
  return_class = "data.table",
  use_parallel = 12
)

saveRDS(pcawg_expo2, file = "data/pcawg_cn_sigs_CN176_activity_stepwise.rds")

get_sig_rec_similarity(solution1000$object$K11, tally_X$nmf_matrix) -> x
summary(x$similarity)

summary(pcawg_expo$similarity)
summary(pcawg_expo2$similarity)

hist(pcawg_expo$similarity, breaks = 100)

scount <- apply(as.matrix(pcawg_expo$rel_activity[, -1]), 1, function(x) {
  sum(x > 0.001)
})

hist(scount, breaks = 100)

saveRDS(solution1000$object$K11, file = "data/pcawg_cn_sigs_CN176_signature.rds")

# Signature profile -------------------------------------------------------

p <- show_sig_profile_loop(solution1000$object$K11,
                           mode = "copynumber",
                           method = "X", style = "cosmic", font_scale = 0.8
)
ggplot2::ggsave("output/pcawg_cn_sigs_loop_c176_BP.pdf", plot = p, width = 14, height = 18)

show_cor(pcawg_expo$abs_activity[, -1])
show_cor(pcawg_expo$rel_activity[, -1])
