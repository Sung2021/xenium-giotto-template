source("scripts/helpers.R")
config <- load_project_config()
ensure_directories(config)

library(Giotto)
xenium_gobject <- load_stage(
  config, "05_visualized", "Run scripts/05_visualize_subcellular.R first."
)

known_markers <- existing_markers(xenium_gobject, config$marker_genes)
if (length(known_markers)) {
  dimFeatPlot2D(
    xenium_gobject,
    feats = known_markers,
    point_size = 0.2,
    point_border_stroke = 0,
    background_color = "black",
    cow_n_col = 4,
    show_legend = FALSE,
    gradient_style = "sequential",
    save_param = list(save_name = "xenium_known_markers", base_height = 8)
  )
}

for (method in config$marker_methods) {
  marker_table <- findMarkers_one_vs_all(
    xenium_gobject,
    cluster_column = "leiden_clus",
    method = method
  )
  write_table(marker_table, config, paste0("xenium_markers_", method, ".csv"))

  if (all(c("cluster", "feats") %in% names(marker_table))) {
    marker_table <- data.table::as.data.table(marker_table)
    top_markers <- marker_table[, utils::head(.SD, config$marker_top_n), by = "cluster"]
    violinPlot(
      xenium_gobject,
      feats = unique(top_markers$feats),
      cluster_column = "leiden_clus",
      save_param = list(save_name = paste0("xenium_top_markers_", method))
    )
  }
}

if (!is.null(config$xenium_annotation_csv)) {
  mapping <- read_annotation_map(config$xenium_annotation_csv, "cell_type")
  xenium_gobject <- annotateGiotto(
    xenium_gobject,
    annotation_vector = annotation_vector(mapping, "cell_type"),
    cluster_column = "leiden_clus",
    name = "cell_type"
  )
  plotUMAP(
    xenium_gobject,
    cell_color = "cell_type",
    point_size = 0.05,
    point_border_stroke = 0,
    save_param = list(save_name = "xenium_umap_cell_type")
  )
  spatInSituPlotPoints(
    xenium_gobject,
    polygon_feat_type = config$aggregation_unit,
    polygon_fill = "cell_type",
    polygon_fill_as_factor = TRUE,
    polygon_line_size = 0,
    save_param = list(save_name = "xenium_spatial_cell_type")
  )
} else {
  message(
    "No XENIUM_CLUSTER_ANNOTATIONS configured. Marker tables are ready for ",
    "manual cluster review; no Xenium cell-type labels were assigned."
  )
}

save_stage(xenium_gobject, config, "06_markers")
