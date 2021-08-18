# Huimin collected TCGA data from GDC portal
# and generate file with format needed by sigminer
# Here I will firstly further clean up the data

# library(data.table)
# x = readRDS("data/TCGA/datamicopy.rds")
# setDT(x)

# colnames(x)[6] = "minor_cn"
# head(x)
# x[, sample := substr(sample, 1, 15)]
# saveRDS(x, file = "data/TCGA/tcga_cn.rds")

# rm(list = ls())

# Generate segment counting matrices
library(sigminer)

# Only focus autosomes
tcga_cn = readRDS("data/TCGA/tcga_cn.rds")
table(tcga_cn$chromosome)
tcga_cn = tcga_cn[!chromosome %in% c("chrX", "chrY")]
cn_obj = read_copynumber(
  tcga_cn,
  seg_cols = c("chromosome", "start", "end", "segVal"),
  add_loh = TRUE,
  loh_min_len = 1e4,
  loh_min_frac = 0.05,
  max_copynumber = 1000L,
  genome_build = "hg38",
  complement = FALSE,
  genome_measure = "called"
)
saveRDS(cn_obj, file = "data/TCGA/tcga_cn_obj.rds")

# Tally step
# Generate the counting matrices
# Same as what have done for PCAWG dataset
library(sigminer)
cn_obj = readRDS("data/tcga_cn_obj.rds")
table(cn_obj@data$chromosome)

tally_X <- sig_tally(cn_obj,
                     method = "X",
                     add_loh = TRUE,
                     cores = 10)
saveRDS(tally_X, file = "data/TCGA/tcga_cn_tally_X.rds")


p <- show_catalogue(tally_X, mode = "copynumber", method = "X",
                    style = "cosmic", y_tr = function(x) log10(x + 1),
                    y_lab = "log10(count +1)")
p
ggplot2::ggsave("output/tcga_catalogs_tally_X.pdf", plot = p, width = 16, height = 2.5)

p <- show_catalogue(tally_X, mode = "copynumber", method = "X",
                    style = "cosmic", y_tr = function(x) log10(x + 1),
                    y_lab = "log10(count +1)", by_context = TRUE)
p
ggplot2::ggsave("output/tcga_catalogs_tally_X_by_context.pdf", plot = p, width = 16, height = 2.5)

head(sort(colSums(tally_X$nmf_matrix)), n = 50)

## classes without LOH labels
tally_X_noLOH <- sig_tally(cn_obj,
                           method = "X",
                           add_loh = FALSE,
                           cores = 10)
saveRDS(tally_X_noLOH$all_matrices, file = "data/TCGA/tcga_cn_tally_X_noLOH.rds")
str(tally_X_noLOH$all_matrices, max.level = 1)

p <- show_catalogue(tally_X_noLOH, mode = "copynumber", method = "X",
                    style = "cosmic", y_tr = function(x) log10(x + 1),
                    y_lab = "log10(count +1)")
p
ggplot2::ggsave("output/tcga_catalogs_tally_X_noLOH.pdf", plot = p, width = 16, height = 2.5)


# Phenotype ---------------------------------------------------------------

download.file("https://pancanatlas.xenahubs.net/download/Survival_SupplementalTable_S1_20171025_xena_sp.gz",
              "Xena/Survival_SupplementalTable_S1_20171025_xena_sp.gz")

tcga_cli <- data.table::fread("Xena/Survival_SupplementalTable_S1_20171025_xena_sp.gz")
table(tcga_cli$`cancer type abbreviation`)
saveRDS(tcga_cli, file = "data/TCGA/tcga_cli.rds")
