# Repair already processed NMF runs based on previous sigminer version

run_files <- list.files("BP/pcawg_bp/", full.names = TRUE)
run_files[1]
run_file_basenames = basename(run_files)

run_info <- dplyr::tibble(
  run_file = run_files,
  basename = run_file_basenames
)

library(dplyr)

run_info <- run_info %>%
  tidyr::separate("basename", into = c("data_tag", "nmf", "K", "seed", "seed_number"), sep = "_") %>%
  dplyr::select(-c(nmf, seed)) %>%
  tidyr::separate("seed_number", into = c("seed", "ext")) %>%
  dplyr::select(-ext)

run_info$seed <- as.integer(run_info$seed)
run_info$data_tag <- NULL

run_info <- run_info %>%
  dplyr::mutate(grp = floor((seed - 123456) / 50))

run_info
table(run_info$grp)

run_nest <- run_info %>%
  group_nest(K)

nrow(run_nest$data[[1]])
table(run_nest$data[[1]]$grp)

run_nest$data[[1]] %>%
  group_split(grp) -> test

test[[1]]

# Source code from sigminer
filter_nmf <- function(s, bt_flag, RTOL) {
  KLD_list <- purrr::map_dbl(s, "KLD")
  if (bt_flag) {
    ki <- KLD_list <= min(KLD_list) * (1 + RTOL)
    s <- s[ki]
    if (length(s) > 10) {
      # Limits 10 best runs
      KLD_list <- KLD_list[ki]
      s <- s[order(KLD_list)[1:10]]
    }
  } else if (length(s) > 100 & !bt_flag) {
    s <- s[order(KLD_list)[1:100]]
  }
  s
}

solutions <- list()

K = run_nest$K
for (i in seq_along(K)) {
  message("Handling signature number #", K[i])
  solutions[[K[i]]] <- purrr::map(
    run_nest$data[[i]] %>% group_split(grp), function(x) {
    nmf_list <- purrr::map(x$run_file, readRDS)
    nmf_list <- filter_nmf(nmf_list, TRUE, 0.001)
    nmf_list
  }) %>% purrr::flatten()
}


tally_X <- readRDS("data/pcawg_cn_tally_X.rds")

# future::plan("multisession", workers = 4, gc = TRUE, .skip = TRUE)
#
# solutions <- furrr::future_map(
#   solutions,
#   .f = sigminer:::process_solution,
#   catalogue_matrix = t(tally_X$nmf_matrix),
#   report_integer_exposure = FALSE,
#   only_core_stats = TRUE,
#   .progress = TRUE,
#   .options = furrr::furrr_options(seed = 123456)
# )

solutions_bk <- solutions
n <- names(solutions)
solutions2 <- list()
for (i in seq_along(solutions)) {
  message("Processing solution: ", n[i])
  solutions2[[n[i]]] <- sigminer:::process_solution(
    solutions[[i]],
    catalogue_matrix = t(handle_hyper_mutation(tally_X$nmf_matrix)),
    report_integer_exposure = FALSE,
    only_core_stats = TRUE
  )
}

saveRDS(solutions2, file = "BP/BP_PCAWG_1000_NMFRuns.rds")

solutions <- solutions2
solutions <- purrr::transpose(solutions)
solutions$stats <- purrr::reduce(solutions$stats, rbind)
solutions$stats_signature <- purrr::reduce(solutions$stats_signature, rbind)
solutions$stats_sample <- purrr::reduce(solutions$stats_sample, rbind)
# 追加属性
solutions$object <- purrr::map(solutions$object, .f = function(obj) {
  attr(obj, "nrun") <- 1000
  attr(obj, "seed") <- 123456
  obj
})


rank_score <- sigminer:::rank_solutions(solutions$stats)
suggested <- rank_score[order(
  rank_score$aggregated_score,
  rank_score$signature_number,
  decreasing = TRUE
), ][1, ]

solutions$rank_score <- rank_score
solutions$suggested <- as.integer(suggested$signature_number)
class(solutions) <- "ExtractionResult"

saveRDS(solutions, file = "BP/BP_PCAWG_1000_Extraction_Result.rds")

bp_show_survey(solutions, add_score = F, fixed_ratio = F)

# 10 or 11 signatures
# Check another 100 NMF runs for make sure
