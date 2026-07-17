env_flag <- function(name, default = FALSE) {
  value <- tolower(Sys.getenv(name, if (default) "true" else "false"))
  value %in% c("1", "true", "yes", "y")
}

project_root <- normalizePath(Sys.getenv("PROJECT_ROOT", "."), mustWork = TRUE)

config <- list(
  project_root = project_root,
  xenium_dir = Sys.getenv("XENIUM_DATA_DIR", file.path(project_root, "data", "outs")),
  results_dir = Sys.getenv("XENIUM_RESULTS_DIR", file.path(project_root, "results")),
  python_path = {
    value <- Sys.getenv("GIOTTO_PYTHON_PATH", "")
    if (nzchar(value)) value else NULL
  },
  load_transcripts = env_flag("XENIUM_LOAD_TRANSCRIPTS", FALSE),
  qv_threshold = as.numeric(Sys.getenv("XENIUM_QV_THRESHOLD", "20")),
  expression_threshold = 1,
  feature_min_cells = 3,
  cell_min_features = 5,
  pca_dims = as.integer(Sys.getenv("XENIUM_PCA_DIMS", "10")),
  neighbor_k = 10,
  leiden_resolution = as.numeric(Sys.getenv("XENIUM_LEIDEN_RESOLUTION", "0.25")),
  leiden_iterations = 100,
  marker_genes = c(
    "EPCAM", "KRT8", "KRT18",
    "PTPRC", "CD3D", "CD8A", "MS4A1",
    "LYZ", "COL1A1", "PECAM1", "MKI67"
  ),
  # Enable only for legacy/pre-release Xenium feature naming.
  legacy_probe_splitting = env_flag("XENIUM_LEGACY_PROBES", FALSE)
)

config$objects_dir <- file.path(config$results_dir, "objects")
config$figures_dir <- file.path(config$results_dir, "figures")
