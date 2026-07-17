load_project_config <- function() {
  config_path <- file.path(Sys.getenv("PROJECT_ROOT", "."), "config", "project_config.R")
  if (!file.exists(config_path)) {
    stop("Run this script from the repository root, or set PROJECT_ROOT.")
  }
  source(config_path, local = TRUE)$value
}

ensure_directories <- function(config) {
  dirs <- c(config$results_dir, config$objects_dir, config$figures_dir)
  invisible(vapply(dirs, dir.create, logical(1), recursive = TRUE, showWarnings = FALSE))
}

require_xenium_directory <- function(path) {
  if (!dir.exists(path)) {
    stop(
      "Xenium directory not found: ", path,
      "\nSet XENIUM_DATA_DIR to the standard Xenium output directory."
    )
  }
}

existing_markers <- function(gobject, markers) {
  available <- Giotto::featIDs(gobject, feat_type = "rna")
  keep <- intersect(markers, available)
  missing <- setdiff(markers, keep)
  if (length(missing)) {
    message("Markers absent from this panel: ", paste(missing, collapse = ", "))
  }
  keep
}
