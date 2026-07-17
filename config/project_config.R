env_flag <- function(name, default = FALSE) {
  value <- tolower(Sys.getenv(name, if (default) "true" else "false"))
  value %in% c("1", "true", "yes", "y")
}

env_number_pair <- function(name) {
  value <- Sys.getenv(name, "")
  if (!nzchar(value)) return(NULL)
  parsed <- as.numeric(strsplit(value, ",", fixed = TRUE)[[1]])
  if (length(parsed) != 2L || anyNA(parsed)) {
    stop(name, " must contain two comma-separated numbers.")
  }
  parsed
}

env_character_vector <- function(name, default = character()) {
  value <- Sys.getenv(name, "")
  if (!nzchar(value)) return(default)
  trimws(strsplit(value, ",", fixed = TRUE)[[1]])
}

env_integer_vector <- function(name, default = integer()) {
  values <- env_character_vector(name)
  if (!length(values)) return(default)
  parsed <- as.integer(values)
  if (anyNA(parsed)) stop(name, " must contain comma-separated integers.")
  parsed
}

env_optional_path <- function(name) {
  value <- Sys.getenv(name, "")
  if (nzchar(value)) value else NULL
}

project_root <- normalizePath(Sys.getenv("PROJECT_ROOT", "."), mustWork = TRUE)

config <- list(
  project_root = project_root,
  xenium_dir = Sys.getenv("XENIUM_DATA_DIR", file.path(project_root, "data", "outs")),
  results_dir = Sys.getenv("XENIUM_RESULTS_DIR", file.path(project_root, "results")),
  python_path = env_optional_path("GIOTTO_PYTHON_PATH"),

  # Xenium import and aggregation
  load_transcripts = env_flag("XENIUM_LOAD_TRANSCRIPTS", FALSE),
  qv_threshold = as.numeric(Sys.getenv("XENIUM_QV_THRESHOLD", "20")),
  aggregation_unit = Sys.getenv("XENIUM_AGGREGATION_UNIT", "cell"),
  legacy_probe_splitting = env_flag("XENIUM_LEGACY_PROBES", FALSE),

  # Optional aligned post-Xenium images
  he_ometif = env_optional_path("XENIUM_HE_OMETIF"),
  he_alignment = env_optional_path("XENIUM_HE_ALIGNMENT"),
  if_ometif = env_optional_path("XENIUM_IF_OMETIF"),
  if_alignment = env_optional_path("XENIUM_IF_ALIGNMENT"),
  if_pages = env_integer_vector("XENIUM_IF_PAGES", 1L),
  if_channel_names = env_character_vector("XENIUM_IF_CHANNEL_NAMES"),
  zoom_xlim = env_number_pair("XENIUM_ZOOM_XLIM"),
  zoom_ylim = env_number_pair("XENIUM_ZOOM_YLIM"),
  zoom_image_name = Sys.getenv("XENIUM_ZOOM_IMAGE", "dapi"),

  # Xenium filtering, reduction, and clustering
  expression_threshold = as.numeric(Sys.getenv("XENIUM_EXPRESSION_THRESHOLD", "1")),
  feature_min_cells = as.integer(Sys.getenv("XENIUM_FEATURE_MIN_CELLS", "3")),
  cell_min_features = as.integer(Sys.getenv("XENIUM_CELL_MIN_FEATURES", "5")),
  pca_dims = as.integer(Sys.getenv("XENIUM_PCA_DIMS", "10")),
  run_tsne = env_flag("XENIUM_RUN_TSNE", TRUE),
  neighbor_k = as.integer(Sys.getenv("XENIUM_NEIGHBOR_K", "10")),
  leiden_resolution = as.numeric(Sys.getenv("XENIUM_LEIDEN_RESOLUTION", "0.25")),
  leiden_iterations = as.integer(Sys.getenv("XENIUM_LEIDEN_ITERATIONS", "100")),

  # Marker exploration. Replace these with tissue- and panel-specific markers.
  marker_genes = env_character_vector(
    "XENIUM_MARKER_GENES",
    c(
      "EPCAM", "KRT8", "MKI67", "PTPRC", "CD3D", "CD8A", "CD4",
      "MS4A1", "JCHAIN", "LYZ", "TYROBP", "COL1A1", "DCN", "PECAM1"
    )
  ),
  marker_methods = env_character_vector("XENIUM_MARKER_METHODS", c("scran", "gini")),
  marker_top_n = as.integer(Sys.getenv("XENIUM_MARKER_TOP_N", "5")),
  xenium_annotation_csv = env_optional_path("XENIUM_CLUSTER_ANNOTATIONS"),

  # Optional scRNA-seq reference processing
  reference_h5 = env_optional_path("GIOTTO_REFERENCE_H5"),
  reference_annotation_csv = env_optional_path("GIOTTO_REFERENCE_ANNOTATIONS"),
  reference_label_column = Sys.getenv("GIOTTO_REFERENCE_LABEL", "cell_type"),
  reference_mito_max = as.numeric(Sys.getenv("GIOTTO_REFERENCE_MITO_MAX", "50")),
  reference_run_scrublet = env_flag("GIOTTO_REFERENCE_RUN_SCRUBLET", TRUE),
  reference_feature_min_cells = as.integer(Sys.getenv("GIOTTO_REFERENCE_FEATURE_MIN_CELLS", "20")),
  reference_cell_min_features = as.integer(Sys.getenv("GIOTTO_REFERENCE_CELL_MIN_FEATURES", "200")),
  reference_pca_dims = as.integer(Sys.getenv("GIOTTO_REFERENCE_PCA_DIMS", "30")),
  reference_leiden_resolution = as.numeric(Sys.getenv("GIOTTO_REFERENCE_LEIDEN_RESOLUTION", "1.5")),

  # Joining, Harmony, and label transfer
  integration_dims = as.integer(Sys.getenv("GIOTTO_INTEGRATION_DIMS", "10")),
  integration_min_shared_features = as.integer(Sys.getenv("GIOTTO_MIN_SHARED_FEATURES", "20")),
  label_transfer_k = as.integer(Sys.getenv("GIOTTO_LABEL_TRANSFER_K", "10"))
)

if (!config$aggregation_unit %in% c("cell", "nucleus")) {
  stop("XENIUM_AGGREGATION_UNIT must be either 'cell' or 'nucleus'.")
}
if (config$aggregation_unit == "nucleus" && !config$load_transcripts) {
  stop(
    "Nuclear aggregation requires XENIUM_LOAD_TRANSCRIPTS=true because the ",
    "10x supplied low-memory expression matrix is cell based."
  )
}

positive_integer_fields <- c(
  "pca_dims", "neighbor_k", "leiden_iterations", "marker_top_n",
  "reference_feature_min_cells", "reference_cell_min_features",
  "reference_pca_dims", "integration_dims",
  "integration_min_shared_features", "label_transfer_k"
)
for (field in positive_integer_fields) {
  value <- config[[field]]
  if (length(value) != 1L || is.na(value) || value < 1L) {
    stop("Configuration value '", field, "' must be a positive integer.")
  }
}

config$objects_dir <- file.path(config$results_dir, "objects")
config$figures_dir <- file.path(config$results_dir, "figures")
config$tables_dir <- file.path(config$results_dir, "tables")
config$image_export_dir <- file.path(config$results_dir, "converted_images")
