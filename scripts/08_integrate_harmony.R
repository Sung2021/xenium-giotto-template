source("scripts/helpers.R")
config <- load_project_config()
ensure_directories(config)

library(Giotto)
set_giotto_python_path(config$python_path)

xenium_gobject <- load_stage(
  config, "06_markers", "Run the core Xenium pipeline through step 06."
)
reference_gobject <- load_stage(
  config, "07_reference", "Run scripts/07_process_scrna_reference.R first."
)
require_metadata_column(
  reference_gobject,
  config$reference_label_column,
  "scRNA-seq reference object"
)

shared_features <- intersect(
  featIDs(xenium_gobject, feat_type = "rna"),
  featIDs(reference_gobject)
)
if (length(shared_features) < config$integration_min_shared_features) {
  stop(
    "Only ", length(shared_features), " shared features were found; at least ",
    config$integration_min_shared_features, " are required."
  )
}
write_table(
  data.table::data.table(feature = shared_features),
  config,
  "integration_shared_features.csv"
)

reference_join <- setGiotto(
  giotto(),
  reference_gobject[c("expression", "spatial_locs"), "raw"]
)
reference_metadata <- data.table::as.data.table(pDataDT(reference_gobject))[
  , c("cell_ID", config$reference_label_column),
  with = FALSE
]
reference_join <- addCellMetadata(
  reference_join,
  new_metadata = reference_metadata,
  by_column = TRUE,
  column_cell_ID = "cell_ID"
)
xenium_join <- setGiotto(
  giotto(),
  xenium_gobject[
    c("expression", "spatial_locs"),
    "raw",
    spat_unit = config$aggregation_unit
  ]
)
xenium_metadata <- data.table::as.data.table(
  pDataDT(xenium_gobject, spat_unit = config$aggregation_unit)
)[, .(cell_ID, original_cell_ID = cell_ID)]
xenium_join <- addCellMetadata(
  xenium_join,
  new_metadata = xenium_metadata,
  spat_unit = config$aggregation_unit,
  by_column = TRUE,
  column_cell_ID = "cell_ID"
)

joined_gobject <- joinGiottoObjects(
  list(reference_join[shared_features], xenium_join[shared_features]),
  gobject_names = c("sc", "xen")
)
instructions(joined_gobject) <- instructions(xenium_gobject)
require_metadata_column(
  joined_gobject,
  config$reference_label_column,
  "joined object"
)

joined_gobject <- joined_gobject |>
  filterGiotto(
    expression_threshold = 1,
    feat_det_in_min_cells = 1,
    min_det_feats_per_cell = 10
  ) |>
  normalizeGiotto() |>
  runPCA(feats_to_use = NULL)

integration_dims <- seq_len(config$integration_dims)
joined_gobject <- runUMAP(
  joined_gobject,
  dimensions_to_use = integration_dims,
  name = "umap_unintegrated"
)
plotUMAP(
  joined_gobject,
  cell_color = "list_ID",
  dim_reduction_name = "umap_unintegrated",
  point_size = 0.3,
  point_border_stroke = 0,
  save_param = list(save_name = "joined_umap_unintegrated")
)

joined_gobject <- runGiottoHarmony(
  joined_gobject,
  vars_use = "list_ID",
  dim_reduction_name = "pca",
  dimensions_to_use = integration_dims,
  name = "harmony"
)
joined_gobject <- runUMAP(
  joined_gobject,
  dim_reduction_to_use = "harmony",
  dim_reduction_name = "harmony",
  name = "umap_harmony",
  dimensions_to_use = integration_dims
)
plotUMAP(
  joined_gobject,
  cell_color = "list_ID",
  dim_reduction_name = "umap_harmony",
  point_size = 0.3,
  point_border_stroke = 0,
  save_param = list(save_name = "joined_umap_harmony_by_source")
)

save_stage(joined_gobject, config, "08_harmony")
