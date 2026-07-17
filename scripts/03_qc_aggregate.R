source("scripts/helpers.R")
config <- load_project_config()
ensure_directories(config)

library(Giotto)
xenium_gobject <- load_stage(
  config, "02_images", "Run scripts/02_attach_images.R first."
)

if (config$load_transcripts) {
  message(
    "Aggregating RNA detections inside ", config$aggregation_unit,
    " polygons at QV >= ", config$qv_threshold, "."
  )
  xenium_gobject <- calculateOverlapRaster(
    xenium_gobject,
    spatial_info = config$aggregation_unit,
    feat_info = "rna"
  )
  xenium_gobject <- overlapToMatrix(
    xenium_gobject,
    poly_info = config$aggregation_unit,
    feat_info = "rna",
    name = "raw"
  )
} else {
  message(
    "Low-memory mode: using the 10x expression matrix and cell metadata. ",
    "Transcript-level QC and de novo polygon aggregation are unavailable."
  )
}

showGiottoExpression(xenium_gobject)

cell_metadata <- pDataDT(xenium_gobject, spat_unit = config$aggregation_unit)
qc_summary <- data.table::data.table(
  metric = c(
    "spatial_unit", "load_transcripts", "qv_threshold",
    "cells_before_filtering", "features_before_filtering"
  ),
  value = c(
    config$aggregation_unit,
    as.character(config$load_transcripts),
    as.character(config$qv_threshold),
    as.character(nrow(cell_metadata)),
    as.character(length(featIDs(xenium_gobject, feat_type = "rna")))
  )
)
write_table(qc_summary, config, "xenium_import_qc_summary.csv")

save_stage(xenium_gobject, config, "03_aggregated")
