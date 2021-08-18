library(sigminer)
# NOTE: the silhouette reported by BP seems strange, need to take a check the clustering with match algo
# sig number 大于 5 左右之后轮廓系数下降的很严重

# test <- sig_estimate(
#   t(pcawg[[5]]$catalog_matrix),
#   range = 2:20, cores = 10, nrun = 100,
#   verbose = TRUE
# )
bp_show_survey2(pcawg[[5]]) # 14

bp_extract_signatures(t(pcawg[[5]]$catalog_matrix), range = 7, n_bootstrap = 5, n_nmf_run = 2,
                      cores = 10, only_core_stats = TRUE)

cancer_type_files <- list.files("data/cancer_types/", pattern = "PCAWG_CN176", full.names = TRUE)
pcawg <- lapply(cancer_type_files, readRDS)

bp_show_survey2(pcawg[[1]]) # 11
bp_show_survey2(pcawg[[2]]) # 14
bp_show_survey2(pcawg[[3]]) # 13
bp_show_survey2(pcawg[[4]]) # 8
bp_show_survey2(pcawg[[5]]) # 14
bp_show_survey2(pcawg[[6]]) # 11
bp_show_survey2(pcawg[[7]]) # 12
bp_show_survey2(pcawg[[8]]) # 8
bp_show_survey2(pcawg[[9]]) # 7
bp_show_survey2(pcawg[[10]]) # 5
bp_show_survey2(pcawg[[11]]) # 13
bp_show_survey2(pcawg[[12]]) # 9
bp_show_survey2(pcawg[[13]]) # 15
bp_show_survey2(pcawg[[14]]) # 7
bp_show_survey2(pcawg[[15]]) # 11
bp_show_survey2(pcawg[[16]]) # 12
bp_show_survey2(pcawg[[17]]) # 12
bp_show_survey2(pcawg[[18]]) # 17
bp_show_survey2(pcawg[[19]]) # 12
bp_show_survey2(pcawg[[20]]) # 6
bp_show_survey2(pcawg[[21]]) # 3
bp_show_survey2(pcawg[[22]]) # 7
bp_show_survey2(pcawg[[23]]) # 11
bp_show_survey2(pcawg[[24]]) # 19
bp_show_survey2(pcawg[[25]]) # 8
bp_show_survey2(pcawg[[26]]) # 10
bp_show_survey2(pcawg[[27]]) # 11
bp_show_survey2(pcawg[[28]]) # 9
bp_show_survey2(pcawg[[29]]) # 9
bp_show_survey2(pcawg[[30]]) # 15
bp_show_survey2(pcawg[[31]]) # 6
bp_show_survey2(pcawg[[32]]) # 17

type_sig <- data.frame(
  cancer_type = sub(".*_CN176_([^_]+).rds", "\\1", basename(cancer_type_files)),
  signum = c(11, 14, 13, 8, 14, 11, 12, 8,
             7, 5, 13, 9, 15, 7, 11, 12,
             12, 17, 12, 6, 3, 7, 11, 19,
             8, 10, 11, 9, 9, 15, 6, 17)
)

pcawg_sig_list <- lapply(
  seq_along(pcawg),
  function(i) {
    bp_get_sig_obj(pcawg[[i]], signum = type_sig$signum[i])
  }
)

names(pcawg_sig_list) <- type_sig$cancer_type
saveRDS(pcawg_sig_list, file = "data/pcawg_type_CNS.rds")


##
pcawg_sig_list <- readRDS("data/pcawg_type_CNS.rds")
pcawg_sig_list$`Biliary-AdenoCA`

debug(bp_cluster_iter_list)
undebug(bp_cluster_iter_list)
cls <- bp_cluster_iter_list(pcawg_sig_list, k = 2:30)

cls$sil_summary
plot(cls$sil_summary$k, cls$sil_summary$mean, type = "b", xlab = "Total signatures", ylab = "Mean silhouette")

# Similar Nature Cancer (extend figure 7,8), we use 0.6 similarity as a cutoff (i.e. distance 0.4)
# dend <- cls$cluster %>% as.dendrogram()
# library(dendextend)
# par(cex=0.1,font=3)

plot(cls$cluster)
factoextra::fviz_dend(cls$cluster, h = 0.4, cex = 0.15,
                      main = "",
                      xlab = "",
                      ylab = "1 - cosine similarity, average linkage")


cls_tree <- cutree(cls$cluster, h = 0.4)
table(cls_tree)

grp_sigs <- bp_get_clustered_sigs(cls, cls_tree)

# Compare similarity
pcawg_sigs <- readRDS("data/pcawg_cn_sigs_CN176_signature.rds")

sim1 <- get_sig_similarity(pcawg_sigs, grp_sigs$grp_sigs)
pheatmap::pheatmap(sim1$similarity, cluster_rows = FALSE, cluster_cols = FALSE, display_numbers = TRUE)

select_sigs <- paste0("Sig", names(table(cls_tree)[table(cls_tree) > 1]))
pheatmap::pheatmap(sim1$similarity[, select_sigs], cluster_rows = FALSE, cluster_cols = FALSE, display_numbers = TRUE)


mc <- apply(sim1$similarity, 2, which.max)
mc_sim <- apply(sim1$similarity, 2, max)
mc_df <- data.frame(
  grp_sig = names(mc),
  global_sig = paste0("CNS", mc),
  sim = as.numeric(mc_sim)
) %>%
  dplyr::mutate(
    label = paste0(grp_sig, "(", global_sig, ", ", sim, ")"),
    grp_sig = factor(grp_sig, levels = paste0("Sig", seq_len(nrow(.))))
  ) %>%
  dplyr::arrange(grp_sig)

ggpubr::ggboxplot(grp_sigs$sim_sig_to_grp_mean, #%>%
                    #dplyr::filter(grp_sig %in% select_sigs),
                  x = "grp_sig", y = "sim", width = 0.3,
                  xlab = FALSE, ylab = "cosine similarity",
                  title = "Cosine similarity between signatures in each group") +
  ggpubr::rotate_x_text()

ggpubr::ggboxplot(grp_sigs$sim_sig_to_grp_mean, #%>%
                  #dplyr::filter(grp_sig %in% select_sigs),
                  x = "grp_sig", y = "sim", width = 0.3,
                  xlab = FALSE, ylab = "Cosine similarity",
                  title = "Cosine similarity between signatures in each group") +
  ggpubr::rotate_x_text() +
  scale_x_discrete(labels = mc_df$label)

# Signature network -------------------------------------------------------

# library(igraph)
# g  <- graph.adjacency(1 - cls$distance, weighted=TRUE)
# #g_mst <- mst(g)
# #plot(g_mst, vertex.size = 3, edge.arrow.size = 0.01)
#
# library(RColorBrewer)
# n <- 60
# qual_col_pals = brewer.pal.info[brewer.pal.info$category == 'qual',]
# col_vector = unlist(mapply(brewer.pal, qual_col_pals$maxcolors, rownames(qual_col_pals)))
# # Make label colors
# set.seed(1)
# col_df <- dplyr::left_join(
#   grp_sigs$sim_sig_to_grp_mean,
#   data.frame(
#     grp_sig = paste0("Sig", seq_len(length(unique(cls_tree)))),
#     col = sample(col_vector, length(unique(cls_tree)))
#   )
# )
# my_color <- col_df$col
# names(my_color) <- col_df$sig
#
# pdf("test.pdf", width = 8, height = 8)
# set.seed(1)
# g_mst <- mst(g)
# plot(g_mst,
#      vertex.size = 0.5,
#      #vertex.frame.color = "white",
#      vertex.label.color = as.character(my_color[colnames(cls$distance)]),
#      vertex.label.cex = 0.3,
#      vertex.label.dist = 0,
#      edge.arrow.size = 0.01,
#      edge.color="lightgray"
#      )
# dev.off()

# res.sim <- 1 - cls$distance
# diag(res.sim) <- NA
# res.sim <- res.sim %>%
#   as.data.frame() %>%
#   tibble::rownames_to_column("term") %>%
#   dplyr::as_tibble()
# library(corrr)
# p <- res.sim[1:9, 1:10] %>%
#   network_plot(min_cor = 0.7,
#                colours = rev(c("indianred2", "white", "skyblue1")))
# p + scale_size(range=c(0.3, 0.5))
#
# p <- res.sim %>%
#   network_plot(min_cor = 0.7,
#                legend = FALSE,
#                colours = rev(c("indianred2", "white", "skyblue1")))
# ggsave("corrr_network.png",
#        plot = p + scale_size(range=c(0.3, 0.5)),
#        height = 20, width = 20, dpi = 300, limitsize = FALSE)
