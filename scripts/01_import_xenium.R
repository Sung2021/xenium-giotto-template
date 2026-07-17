source("scripts/helpers.R")
config <- load_project_config()
ensure_directories(config)
require_directory(config$xenium_dir, "Xenium output directory")

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

message(
  "Importing Xenium data in ",
  if (config$load_transcripts) "full transcript" else "low-memory matrix",
  " mode from: ", config$xenium_dir
)
xenium_gobject <- do.call(createGiottoXeniumObject, import_args)

save_stage(xenium_gobject, config, "01_imported")
