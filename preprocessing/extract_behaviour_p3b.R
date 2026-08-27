#!/usr/bin/env Rscript
# =====================================================================
# Extract P3b task behaviour (accuracy + RT) from the BIDS events files
# =====================================================================

setwd("/Users/miftahfaizah/Library/CloudStorage/OneDrive-UniversityofLeeds/PHD JOURNEY/YBMAP/CN_DATASET/CN_analysis/YBMAP_P3b_review/2_Preprocessing/ND_26")

beh_dir <- "bids_data_behavior"
summary_file <- "preprocessing_output/p3b_erp_summary_autoreject_nd_26_allch.csv"
qc_file <- "manual_channel_qc_nd26.csv"
out_file <- "behaviour_p3b_nd26.csv"
invalid_threshold <- 0.30   # subjects with > 30% multi-press / invalid trials are excluded (accuracy unreliable)

is_true <- function(x) {
  tolower(trimws(as.character(x))) %in% c("true", "1", "yes", "t")
}

read_csv_base <- function(path, ...) {
  read.csv(path, stringsAsFactors = FALSE, check.names = FALSE, ...)
}

summ <- read_csv_base(summary_file)
qc   <- read_csv_base(qc_file)

# Included analysis subjects = AutoReject pass + not manually excluded.
manual_excluded <- unique(qc$subject[is_true(qc$exclude_subject)])
included <- sort(setdiff(
  unique(summ$subject[is_true(summ$analysis_include)]),
  manual_excluded
))

calc_accuracy <- function(response_correct, mask, valid) {
  m <- mask & valid
  if (!sum(m)) return(NA_real_)
  round(mean(response_correct[m] == "True"), 3)
}

empty_row <- function(subject, source) {
  data.frame(
    subject = subject,
    acc_overall = NA_real_,
    acc_target = NA_real_,
    acc_nontarget = NA_real_,
    rt_target_mean = NA_real_,
    rt_overall_mean = NA_real_,
    n_valid = NA_real_,
    invalid_rate = NA_real_,
    source = source,
    stringsAsFactors = FALSE,
    check.names = FALSE
  )
}

rows <- lapply(included, function(sub) {
  event_file <- file.path(beh_dir, paste0(sub, "_task-p3b_events.tsv"))

  if (!file.exists(event_file)) {
    return(empty_row(sub, "edf_only_no_behaviour"))
  }

  header <- names(read.delim(event_file, nrows = 1, stringsAsFactors = FALSE, check.names = FALSE))
  if (!"response_correct" %in% header) {
    return(empty_row(sub, "edf_only_no_behaviour"))
  }

  d <- read.delim(
    event_file,
    stringsAsFactors = FALSE,
    check.names = FALSE,
    colClasses = "character"
  )
  d <- d[d$trial_type %in% c("target", "non-target"), , drop = FALSE]

  rc <- as.character(d$response_correct)
  valid <- rc %in% c("True", "False")
  rt <- suppressWarnings(as.numeric(d$response_rt))

  # exclude subjects with too many invalid (multi-press) trials - accuracy unreliable
  inv_rate <- round(mean(rc == "invalid"), 3)
  if (!is.na(inv_rate) && inv_rate > invalid_threshold) {
    row <- empty_row(sub, "high_invalid_excluded")
    row$invalid_rate <- inv_rate
    row$n_valid <- as.numeric(sum(valid))
    return(row)
  }

  target <- d$trial_type == "target"
  non_target <- d$trial_type == "non-target"
  correct <- rc == "True"

  data.frame(
    subject = sub,
    acc_overall = calc_accuracy(rc, rep(TRUE, length(rc)), valid),
    acc_target = calc_accuracy(rc, target, valid),
    acc_nontarget = calc_accuracy(rc, non_target, valid),
    rt_target_mean = round(mean(rt[target & correct], na.rm = TRUE), 3),
    rt_overall_mean = round(mean(rt[correct], na.rm = TRUE), 3),
    n_valid = as.numeric(sum(valid)),
    invalid_rate = round(mean(rc == "invalid"), 3),
    source = "psychopy",
    stringsAsFactors = FALSE,
    check.names = FALSE
  )
})

beh <- do.call(rbind, rows)
write.csv(beh, out_file, row.names = FALSE, quote = FALSE, na = "")

message(
  "wrote ", basename(out_file), ": ", nrow(beh), " subjects (",
  sum(!is.na(beh$acc_overall)), " with behaviour, ",
  sum(is.na(beh$acc_overall)), " EDF-only -> NA)"
)
