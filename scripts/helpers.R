`%||%` <- function(value, fallback) {
  if (is.null(value) || !length(value)) fallback else value
}

load_project_config <- function() {
  config_path <- file.path(Sys.getenv("PROJECT_ROOT", "."), "config", "project_config.R")
  if (!file.exists(config_path)) {
    stop("Run this script from the repository root, or set PROJECT_ROOT.")
  }
  source(config_path, local = TRUE)$value
}

ensure_directories <- function(config) {
  dirs <- c(
    config$results_dir, config$objects_dir, config$figures_dir,
    config$tables_dir, config$image_export_dir
  )
  invisible(vapply(dirs, dir.create, logical(1), recursive = TRUE, showWarnings = FALSE))
}

require_directory <- function(path, label = "directory") {
  if (!dir.exists(path)) stop(label, " not found: ", path)
  invisible(path)
}

require_file <- function(path, label = "file") {
  if (is.null(path) || !nzchar(path) || !file.exists(path)) {
    stop(label, " not found: ", path %||% "<not configured>")
  }
  invisible(path)
}

stage_path <- function(config, stage) {
  file.path(config$objects_dir, paste0(stage, ".rds"))
}

load_stage <- function(config, stage, hint = NULL) {
  path <- stage_path(config, stage)
  if (!file.exists(path)) {
    stop(
      "Required stage object not found: ", path,
      if (!is.null(hint)) paste0("\n", hint) else ""
    )
  }
  readRDS(path)
}

save_stage <- function(object, config, stage) {
  path <- stage_path(config, stage)
  saveRDS(object, path)
  message("Saved stage object: ", path)
  invisible(path)
}

existing_markers <- function(gobject, markers, feat_type = "rna") {
  available <- Giotto::featIDs(gobject, feat_type = feat_type)
  keep <- intersect(markers, available)
  missing <- setdiff(markers, keep)
  if (length(missing)) {
    message("Markers absent from this panel: ", paste(missing, collapse = ", "))
  }
  keep
}

read_annotation_map <- function(path, label_column) {
  require_file(path, "annotation mapping CSV")
  mapping <- data.table::fread(path)
  required <- c("cluster", label_column)
  missing <- setdiff(required, names(mapping))
  if (length(missing)) {
    stop("Annotation CSV is missing columns: ", paste(missing, collapse = ", "))
  }
  mapping <- unique(mapping[, ..required])
  if (anyNA(mapping[["cluster"]]) || anyNA(mapping[[label_column]])) {
    stop("Annotation CSV contains missing cluster or label values.")
  }
  mapping
}

annotation_vector <- function(mapping, label_column) {
  stats::setNames(as.character(mapping[[label_column]]), as.character(mapping$cluster))
}

require_metadata_column <- function(gobject, column, label = "Giotto object") {
  metadata <- Giotto::pDataDT(gobject)
  if (!column %in% names(metadata)) {
    stop(label, " does not contain metadata column '", column, "'.")
  }
  invisible(column)
}

write_table <- function(value, config, filename) {
  path <- file.path(config$tables_dir, filename)
  data.table::fwrite(value, path)
  message("Saved table: ", path)
  invisible(path)
}
