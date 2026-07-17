source("scripts/helpers.R")
config <- load_project_config()
ensure_directories(config)

library(Giotto)

input_path <- file.path(config$objects_dir, "02_clustered.rds")
if (!file.exists(input_path)) stop("Run scripts/02_process_cluster.R first.")
xenium_gobject <- readRDS(input_path)

plotUMAP(
  xenium_gobject,
  cell_color = "leiden_clus",
  point_size = 0.1,
  point_shape = "no_border",
  show_legend = TRUE,
  save_param = list(save_name = "umap_leiden")
)

spatPlot2D(
  xenium_gobject,
  plot_method = "scattermore",
  cell_color = "leiden_clus",
  point_size = 0.1,
  point_shape = "no_border",
  background_color = "black",
  save_param = list(save_name = "spatial_leiden")
)

markers <- existing_markers(xenium_gobject, config$marker_genes)
if (length(markers)) {
  dimFeatPlot2D(
    xenium_gobject,
    feats = markers,
    point_size = 0.2,
    point_border_stroke = 0,
    background_color = "black",
    show_legend = FALSE,
    cow_n_col = 4,
    save_param = list(save_name = "marker_umap")
  )
}

if (config$load_transcripts && length(markers) >= 3) {
  spatInSituPlotPoints(
    xenium_gobject,
    feats = list(rna = markers[seq_len(3)]),
    point_size = 0.1,
    polygon_feat_type = "cell",
    polygon_fill = "leiden_clus",
    polygon_fill_as_factor = TRUE,
    polygon_alpha = 0.25,
    polygon_line_size = 0.01,
    save_param = list(save_name = "transcripts_spatial_example")
  )
}

message("Plots were written to: ", config$figures_dir)
