source("scripts/helpers.R")
config <- load_project_config()
ensure_directories(config)

library(Giotto)
xenium_gobject <- load_stage(
  config, "01_imported", "Run scripts/01_import_xenium.R first."
)

has_he <- !is.null(config$he_ometif) || !is.null(config$he_alignment)
has_if <- !is.null(config$if_ometif) || !is.null(config$if_alignment)

if (!has_he && !has_if) {
  message(
    "No post-Xenium aligned images configured. Keeping morphology images ",
    "loaded from the standard Xenium output."
  )
  save_stage(xenium_gobject, config, "02_images")
  quit(save = "no", status = 0)
}

set_giotto_python_path(config$python_path)
xenium_loader <- importXenium(config$xenium_dir)
image_objects <- list()

if (has_he) {
  require_file(config$he_ometif, "H&E OME-TIFF")
  require_file(config$he_alignment, "H&E alignment CSV")
  converted_he <- GiottoClass::ometif_to_tif(
    config$he_ometif,
    output_dir = config$image_export_dir
  )
  image_objects[["post_he"]] <- xenium_loader$load_aligned_image(
    name = "post_he",
    path = converted_he,
    imagealignment_path = config$he_alignment
  )
}

if (has_if) {
  require_file(config$if_ometif, "IF OME-TIFF")
  require_file(config$if_alignment, "IF alignment CSV")

  if_pages <- config$if_pages
  if_names <- config$if_channel_names
  if (!length(if_names)) {
    metadata_names <- GiottoClass::ometif_metadata(
      config$if_ometif,
      node = "Channel"
    )$Name
    if_names <- as.character(metadata_names[if_pages])
  }
  if (length(if_names) != length(if_pages) || any(!nzchar(if_names))) {
    stop("IF page count and XENIUM_IF_CHANNEL_NAMES must match.")
  }

  for (index in seq_along(if_pages)) {
    converted_if <- GiottoClass::ometif_to_tif(
      config$if_ometif,
      page = if_pages[[index]],
      output_dir = config$image_export_dir
    )
    channel_name <- if_names[[index]]
    image_objects[[channel_name]] <- xenium_loader$load_aligned_image(
      name = channel_name,
      path = converted_if,
      imagealignment_path = config$if_alignment
    )
  }
}

xenium_gobject <- setGiotto(xenium_gobject, image_objects)
save_stage(xenium_gobject, config, "02_images")
