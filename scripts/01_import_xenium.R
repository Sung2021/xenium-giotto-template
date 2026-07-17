source("scripts/helpers.R")
config <- load_project_config()
ensure_directories(config)
require_xenium_directory(config$xenium_dir)

library(Giotto)
set_giotto_python_path(config$python_path)

instructions <- createGiottoInstructions(
  save_dir = config$figures_dir,
  save_plot = TRUE,
  show_plot = FALSE,
  return_plot = FALSE,
  python_path = config$python_path
)

import_args <- list(
  xenium_dir = config$xenium_dir,
  qv_threshold = config$qv_threshold,
  load_transcripts = config$load_transcripts,
  load_expression = !config$load_transcripts,
  load_cellmeta = !config$load_transcripts,
  instructions = instructions
)

if (config$legacy_probe_splitting) {
  import_args$feat_type <- c(
    "rna", "UnassignedCodeword", "NegControlCodeword", "NegControlProbe"
  )
  import_args$split_keyword <- list(
    c("BLANK", "UnassignedCodeword"),
    "NegControlCodeword",
    c("NegControlProbe", "antisense")
  )
}

message("Importing Xenium data from: ", config$xenium_dir)
xenium_gobject <- do.call(createGiottoXeniumObject, import_args)

if (config$load_transcripts) {
  message("Aggregating RNA transcript detections within cell polygons.")
  xenium_gobject <- calculateOverlapRaster(
    xenium_gobject,
    spatial_info = "cell",
    feat_info = "rna"
  )
  xenium_gobject <- overlapToMatrix(
    xenium_gobject,
    poly_info = "cell",
    feat_info = "rna",
    name = "raw"
  )
}

output_path <- file.path(config$objects_dir, "01_imported.rds")
saveRDS(xenium_gobject, output_path)
message("Saved imported object: ", output_path)
