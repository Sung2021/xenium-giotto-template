source("scripts/helpers.R")
config <- load_project_config()
ensure_directories(config)

library(Giotto)
set_giotto_python_path(config$python_path)
require_file(config$reference_h5, "scRNA-seq 10x H5 reference")

instructions <- createGiottoInstructions(
  save_dir = config$figures_dir,
  save_plot = TRUE,
  show_plot = FALSE,
  return_plot = FALSE,
  python_path = config$python_path
)

reference_matrix <- get10Xmatrix_h5(config$reference_h5)$`Gene Expression`
reference_gobject <- createGiottoObject(
  expression = reference_matrix,
  instructions = instructions
)

mitochondrial_features <- grep(
  "^MT-",
  featIDs(reference_gobject),
  value = TRUE
)
if (length(mitochondrial_features)) {
  reference_gobject <- addFeatsPerc(
    reference_gobject,
    expression_values = "raw",
    feats = mitochondrial_features,
    vector_name = "mito"
  )
  mito_max <- config$reference_mito_max
  reference_gobject <- subset(reference_gobject, subset = mito < mito_max)
} else {
  warning("No ^MT- mitochondrial features were found; mitochondrial QC was skipped.")
}

if (config$reference_run_scrublet) {
  reference_gobject <- doScrubletDetect(
    reference_gobject,
    return_gobject = TRUE,
    seed = 0
  )
  reference_gobject <- subset(reference_gobject, subset = doublet == FALSE)
}

reference_gobject <- reference_gobject |>
  filterGiotto(
    expression_threshold = 1,
    feat_det_in_min_cells = config$reference_feature_min_cells,
    min_det_feats_per_cell = config$reference_cell_min_features
  ) |>
  normalizeGiotto() |>
  calculateHVF() |>
  runPCA()

reference_dims <- seq_len(config$reference_pca_dims)
reference_gobject <- reference_gobject |>
  runUMAP(dimensions_to_use = reference_dims) |>
  createNearestNetwork(dimensions_to_use = reference_dims) |>
  doLeidenCluster(resolution = config$reference_leiden_resolution)

screePlot(
  reference_gobject,
  ncp = max(50L, config$reference_pca_dims),
  save_param = list(save_name = "reference_scree")
)
dimPlot2D(
  reference_gobject,
  cell_color = "leiden_clus",
  point_size = 0.4,
  point_border_stroke = 0,
  save_param = list(save_name = "reference_umap_leiden")
)

known_markers <- existing_markers(reference_gobject, config$marker_genes)
if (length(known_markers)) {
  dimFeatPlot2D(
    reference_gobject,
    feats = known_markers,
    point_size = 0.3,
    point_border_stroke = 0,
    background_color = "black",
    cow_n_col = 4,
    show_legend = FALSE,
    save_param = list(save_name = "reference_known_markers", base_height = 8)
  )
}

for (method in config$marker_methods) {
  reference_markers <- findMarkers_one_vs_all(
    reference_gobject,
    cluster_column = "leiden_clus",
    method = method
  )
  write_table(
    reference_markers,
    config,
    paste0("reference_markers_", method, ".csv")
  )
}

if (!is.null(config$reference_annotation_csv)) {
  label_column <- config$reference_label_column
  mapping <- read_annotation_map(config$reference_annotation_csv, label_column)
  reference_gobject <- annotateGiotto(
    reference_gobject,
    annotation_vector = annotation_vector(mapping, label_column),
    cluster_column = "leiden_clus",
    name = label_column
  )
  dimPlot2D(
    reference_gobject,
    cell_color = label_column,
    point_size = 0.4,
    point_border_stroke = 0,
    save_param = list(save_name = "reference_umap_cell_type")
  )
} else {
  message(
    "Reference processing is complete but no GIOTTO_REFERENCE_ANNOTATIONS CSV ",
    "was supplied. Review marker tables, create the mapping, and rerun this step ",
    "before Harmony/label transfer."
  )
}

save_stage(reference_gobject, config, "07_reference")
