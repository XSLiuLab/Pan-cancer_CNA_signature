library(sigminer)

tally_X <- readRDS("tcga_cn_tally_X.rds")
tally_X_noLOH <- readRDS("tcga_cn_tally_X_noLOH.rds")

sigprofiler_extract(tally_X$nmf_matrix, 
                    output = "TCGA_CN176X", 
                    range = 2:30, 
                    nrun = 100,
                    cores = 36,
                    init_method = "random",
                    is_exome = FALSE,
                    use_conda = TRUE)

sigprofiler_extract(tally_X_noLOH$simplified_matrix, 
                    output = "TCGA_CN136", 
                    range = 2:30, 
                    nrun = 100,
                    cores = 36,
                    init_method = "random",
                    is_exome = FALSE,
                    use_conda = TRUE)


