library(sigminer)

# Only focus autosome, suggested by Prof. Liu
pcawg_cn <- readRDS("data/pcawg_copynumber_sp.rds")
table(pcawg_cn$Chromosome)
pcawg_cn <- pcawg_cn[!Chromosome %in% c("X", "Y")]
cn_obj <- read_copynumber(pcawg_cn,
                          add_loh = TRUE,
                          loh_min_len = 1e4,
                          loh_min_frac = 0.05,
                          max_copynumber = 1000L,
                          genome_build = "hg19",
                          complement = FALSE,
                          genome_measure = "called")
saveRDS(cn_obj, file = "data/pcawg_cn_obj.rds")

# Check the LOH length fraction in a segment in merged situation
loh_frac = na.omit(cn_obj@data$.loh_frac)
hist(loh_frac, breaks = 1000)
hist(loh_frac[loh_frac < 0.9], breaks = 1000)
hist(loh_frac[loh_frac < 0.1], breaks = 1000)
abline(v = 0.05, col = "red")
#cn_obj@data[order(cn_obj@data$.loh_frac)][1000:1100]

# Tally -------------------------------------------------------------------
library(sigminer)
cn_obj <- readRDS("data/pcawg_cn_obj.rds")
table(cn_obj@data$chromosome)

tally_X <- sig_tally(cn_obj,
                     method = "X",
                     add_loh = TRUE,
                     cores = 10)
saveRDS(tally_X, file = "data/pcawg_cn_tally_X.rds")
str(tally_X$all_matrices, max.level = 1)

p <- show_catalogue(tally_X, mode = "copynumber", method = "X",
                    style = "cosmic", y_tr = function(x) log10(x + 1),
                    y_lab = "log10(count +1)")
p
ggplot2::ggsave("output/pcawg_catalogs_tally_X.pdf", plot = p, width = 16, height = 2.5)

p <- show_catalogue(tally_X, mode = "copynumber", method = "X",
                    style = "cosmic", y_tr = function(x) log10(x + 1),
                    y_lab = "log10(count +1)", by_context = TRUE)
p
ggplot2::ggsave("output/pcawg_catalogs_tally_X_by_context.pdf", plot = p, width = 16, height = 2.5)

head(sort(colSums(tally_X$nmf_matrix)), n = 50)

## classes without LOH labels
tally_X_noLOH <- sig_tally(cn_obj,
                           method = "X",
                           add_loh = FALSE,
                           cores = 10)
saveRDS(tally_X_noLOH$all_matrices, file = "data/pcawg_cn_tally_X_noLOH.rds")
str(tally_X_noLOH$all_matrices, max.level = 1)

p <- show_catalogue(tally_X_noLOH, mode = "copynumber", method = "X",
                    style = "cosmic", y_tr = function(x) log10(x + 1),
                    y_lab = "log10(count +1)")
p
ggplot2::ggsave("output/pcawg_catalogs_tally_X_noLOH.pdf", plot = p, width = 16, height = 2.5)
