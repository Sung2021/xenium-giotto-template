source("scripts/helpers.R")
config <- load_project_config()
ensure_directories(config)

library(Giotto)

xenium_path <- file.path(config$objects_dir, "02_clustered.rds")
reference_path <- Sys.getenv("GIOTTO_REFERENCE_OBJECT", "")
label_column <- Sys.getenv("GIOTTO_REFERENCE_LABEL", "cell_type")

if (!file.exists(xenium_path)) stop("Run scripts/02_process_cluster.R first.")
if (!nzchar(reference_path) || !file.exists(reference_path)) {
  stop(
    "Set GIOTTO_REFERENCE_OBJECT to a processed, annotated Giotto scRNA-seq RDS. ",
    "The object must contain the label column selected by GIOTTO_REFERENCE_LABEL."
  )
}

xenium_gobject <- readRDS(xenium_path)
reference_gobject <- readRDS(reference_path)

message(
  "Reference label transfer is dataset-specific. Inspect feature overlap, QC, ",
  "batch effects, and reference annotations before accepting the result."
)

transfer_result <- labelTransfer(
  xenium_gobject,
  reference_gobject,
  labels = label_column,
  integration_method = "harmony",
  k = 10,
  dimensions_to_use = seq_len(config$pca_dims)
)

output_path <- file.path(config$objects_dir, "04_label_transfer.rds")
saveRDS(transfer_result, output_path)
message("Saved label-transfer result: ", output_path)
