source("scripts/helpers.R")
config <- load_project_config()
ensure_directories(config)

library(Giotto)
xenium_gobject <- load_stage(
  config, "04_clustered", "Run scripts/04_process_cluster.R first."
)

if (config$run_tsne) {
  plotTSNE(
    xenium_gobject,
    cell_color = "leiden_clus",
    show_legend = FALSE,
    point_size = 0.05,
    point_shape = "no_border",
    save_param = list(save_name = "xenium_tsne_leiden")
  )
}

plotUMAP(
  xenium_gobject,
  cell_color = "leiden_clus",
  point_size = 0.05,
  point_shape = "no_border",
  save_param = list(save_name = "xenium_umap_leiden")
)

spatPlot2D(
  xenium_gobject,
  spat_unit = config$aggregation_unit,
  plot_method = "scattermore",
  cell_color = "leiden_clus",
  point_size = 0.1,
  point_shape = "no_border",
  background_color = "black",
  save_param = list(save_name = "xenium_spatial_leiden")
)

spatInSituPlotPoints(
  xenium_gobject,
  polygon_feat_type = config$aggregation_unit,
  polygon_alpha = 1,
  polygon_line_size = 0.01,
  polygon_color = "black",
  polygon_fill = "leiden_clus",
  polygon_fill_as_factor = TRUE,
  save_param = list(save_name = "xenium_polygons_leiden")
)

if (config$load_transcripts) {
  transcript_markers <- existing_markers(xenium_gobject, config$marker_genes)
  transcript_markers <- utils::head(transcript_markers, 3L)

  if (length(transcript_markers)) {
    plot_args <- list(
      gobject = xenium_gobject,
      feats = list(rna = transcript_markers),
      point_size = 0.1,
      plot_last = "polygons",
      polygon_feat_type = config$aggregation_unit,
      polygon_alpha = 0.3,
      polygon_line_size = 0.01,
      polygon_color = "black",
      polygon_fill = "leiden_clus",
      polygon_fill_as_factor = TRUE,
      show_image = TRUE,
      image_name = config$zoom_image_name,
      save_param = list(save_name = "xenium_transcripts_polygons")
    )
    if (!is.null(config$zoom_xlim) && !is.null(config$zoom_ylim)) {
      plot_args$xlim <- config$zoom_xlim
      plot_args$ylim <- config$zoom_ylim
    } else {
      message(
        "No XENIUM_ZOOM_XLIM/YLIM configured; transcript plot uses the full field. ",
        "Set both variables for a practical subcellular zoom."
      )
    }
    do.call(spatInSituPlotPoints, plot_args)
  }
} else {
  message("Skipping transcript-point plots because XENIUM_LOAD_TRANSCRIPTS=false.")
}

save_stage(xenium_gobject, config, "05_visualized")
