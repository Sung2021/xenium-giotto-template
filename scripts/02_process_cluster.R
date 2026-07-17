source("scripts/helpers.R")
config <- load_project_config()
ensure_directories(config)

library(Giotto)
set_giotto_python_path(config$python_path)

input_path <- file.path(config$objects_dir, "01_imported.rds")
if (!file.exists(input_path)) stop("Run scripts/01_import_xenium.R first.")
xenium_gobject <- readRDS(input_path)

xenium_gobject <- xenium_gobject |>
  filterGiotto(
    spat_unit = "cell",
    expression_threshold = config$expression_threshold,
    feat_det_in_min_cells = config$feature_min_cells,
    min_det_feats_per_cell = config$cell_min_features
  ) |>
  normalizeGiotto() |>
  addStatistics() |>
  runPCA(feats_to_use = NULL)

dims <- seq_len(config$pca_dims)

xenium_gobject <- xenium_gobject |>
  runUMAP(dimensions_to_use = dims) |>
  createNearestNetwork(dimensions_to_use = dims, k = config$neighbor_k) |>
  doLeidenCluster(
    resolution = config$leiden_resolution,
    n_iterations = config$leiden_iterations
  )

output_path <- file.path(config$objects_dir, "02_clustered.rds")
saveRDS(xenium_gobject, output_path)
message("Saved clustered object: ", output_path)
