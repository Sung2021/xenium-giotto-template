source("scripts/helpers.R")
config <- load_project_config()
ensure_directories(config)

library(Giotto)
library(data.table)

xenium_gobject <- load_stage(
  config, "06_markers", "Run the core Xenium pipeline through step 06."
)
joined_gobject <- load_stage(
  config, "08_harmony", "Run scripts/08_integrate_harmony.R first."
)

label_column <- config$reference_label_column
require_metadata_column(joined_gobject, label_column, "joined Harmony object")

source_cell_ids <- spatIDs(joined_gobject, subset = list_ID == "sc")
joined_gobject <- labelTransfer(
  joined_gobject,
  source_cell_ids = source_cell_ids,
  k = config$label_transfer_k,
  labels = label_column,
  reduction_method = "harmony",
  reduction_name = "harmony",
  dimensions_to_use = seq_len(config$integration_dims)
)

transfer_column <- paste0("trnsfr_", label_column)
probability_column <- paste0(transfer_column, "_prob")
require_metadata_column(joined_gobject, transfer_column, "label-transfer result")
require_metadata_column(joined_gobject, probability_column, "label-transfer result")

plotUMAP(
  joined_gobject,
  cell_color = transfer_column,
  dim_reduction_name = "umap_harmony",
  point_size = 0.3,
  point_border_stroke = 0,
  save_param = list(save_name = "joined_umap_transferred_labels")
)

transfer_metadata <- as.data.table(pDataDT(joined_gobject))
transfer_metadata <- transfer_metadata[
  list_ID == "xen",
  c("original_cell_ID", transfer_column, probability_column),
  with = FALSE
]
if (anyNA(transfer_metadata$original_cell_ID)) {
  stop("Joined Xenium metadata contains missing original cell identifiers.")
}
if (anyDuplicated(transfer_metadata$original_cell_ID)) {
  stop("Joined Xenium metadata contains duplicated original cell identifiers.")
}
data.table::setnames(transfer_metadata, "original_cell_ID", "cell_ID")

xenium_gobject <- addCellMetadata(
  xenium_gobject,
  new_metadata = transfer_metadata,
  spat_unit = config$aggregation_unit,
  by_column = TRUE,
  column_cell_ID = "cell_ID"
)

spatInSituPlotPoints(
  xenium_gobject,
  polygon_feat_type = config$aggregation_unit,
  polygon_fill = transfer_column,
  polygon_fill_as_factor = TRUE,
  polygon_line_size = 0,
  save_param = list(save_name = "xenium_spatial_transferred_labels")
)

write_table(transfer_metadata, config, "xenium_transferred_labels.csv")
save_stage(joined_gobject, config, "09_joined_labeled")
save_stage(xenium_gobject, config, "09_xenium_labeled")
