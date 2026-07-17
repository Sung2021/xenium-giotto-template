source("scripts/helpers.R")
config <- load_project_config()
ensure_directories(config)

library(Giotto)
set_giotto_python_path(config$python_path)
xenium_gobject <- load_stage(
  config, "03_aggregated", "Run scripts/03_qc_aggregate.R first."
)

xenium_gobject <- xenium_gobject |>
  filterGiotto(
    spat_unit = config$aggregation_unit,
    expression_threshold = config$expression_threshold,
    feat_det_in_min_cells = config$feature_min_cells,
    min_det_feats_per_cell = config$cell_min_features
  ) |>
  normalizeGiotto() |>
  addStatistics() |>
  runPCA(feats_to_use = NULL)

dims <- seq_len(config$pca_dims)

screePlot(
  xenium_gobject,
  ncp = max(20L, config$pca_dims),
  save_param = list(save_name = "xenium_scree")
)
plotPCA(
  xenium_gobject,
  point_size = 0.1,
  save_param = list(save_name = "xenium_pca")
)

if (config$run_tsne) {
  xenium_gobject <- runtSNE(xenium_gobject, dimensions_to_use = dims)
}
xenium_gobject <- xenium_gobject |>
  runUMAP(dimensions_to_use = dims) |>
  createNearestNetwork(dimensions_to_use = dims, k = config$neighbor_k) |>
  doLeidenCluster(
    resolution = config$leiden_resolution,
    n_iterations = config$leiden_iterations
  )

save_stage(xenium_gobject, config, "04_clustered")
