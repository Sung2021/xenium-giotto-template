# Xenium Giotto analysis template

Reusable R workflow for importing, quality-checking, clustering, visualizing, and annotating 10x Genomics Xenium data with Giotto Suite.

The project follows the full analysis arc of Giotto's [Xenium breast cancer tutorial](https://giottosuite.com/articles/xenium_breast_cancer.html), but it is not a copy of the legacy breast-cancer analysis. Dataset paths, image channels, QC thresholds, marker genes, cluster annotations, and the scRNA-seq reference are configurable.

## Workflow coverage

| Step | Script | What it does | Required? |
|---:|---|---|---|
| 00 | `00_setup.R` | Installs Giotto, Arrow/ZSTD, Python image codecs, and Scrublet | Once per environment |
| 01 | `01_import_xenium.R` | Imports standard Xenium outputs in low-memory or transcript mode | Yes |
| 02 | `02_attach_images.R` | Converts OME-TIFF and attaches aligned H&E/IF images | Optional, pass-through if unset |
| 03 | `03_qc_aggregate.R` | Records import QC and aggregates transcripts within cell or nuclear polygons | Yes |
| 04 | `04_process_cluster.R` | Filtering, normalization, PCA, t-SNE, UMAP, neighbors, and Leiden clustering | Yes |
| 05 | `05_visualize_subcellular.R` | UMAP, spatial centroids, polygons, images, and transcript-point plots | Yes; transcript plots require transcript mode |
| 06 | `06_detect_markers.R` | Known-marker plots, Scran/Gini marker tables, and optional manual annotations | Yes |
| 07 | `07_process_scrna_reference.R` | scRNA QC, mitochondrial filtering, Scrublet, HVFs, clustering, markers, and manual annotation | Optional and dataset-specific |
| 08 | `08_integrate_harmony.R` | Feature matching, object joining, unintegrated comparison, and Harmony integration | Optional and dataset-specific |
| 09 | `09_transfer_labels.R` | kNN label transfer, confidence table, UMAP, and spatial label visualization | Optional and dataset-specific |

The first six analysis steps form the reusable Xenium pipeline. Steps 07-09 require a biologically suitable, manually reviewed scRNA-seq reference and are never treated as an automatic default.

## Requirements

- R 4.4 or newer is recommended
- Python and Conda accessible through `reticulate`
- A standard Xenium output directory with `transcripts.parquet`, `cells.parquet`, polygon files, and `gene_panel.json`
- A high-memory machine when individual transcript detections are loaded

The source tutorial reports approximately 70 GB peak RAM with transcripts and 5 GB without transcripts for its full legacy dataset. Treat these as rough planning figures, not guarantees for other Xenium runs.

## Quick start

Clone the template, configure paths, and run the one-time setup:

```bash
git clone https://github.com/Sung2021/xenium-giotto-template.git
cd xenium-giotto-template

export XENIUM_DATA_DIR=/absolute/path/to/xenium/outs
export XENIUM_RESULTS_DIR=/absolute/path/to/results

Rscript scripts/00_setup.R
```

Run the core Xenium workflow:

```bash
bash run_pipeline.sh
```

`scripts/00_setup.R` is intentionally not part of the runner because dependency installation should normally happen once. The runner executes steps 01-06 every time.

## Import modes

### Low-memory mode: default

```bash
export XENIUM_LOAD_TRANSCRIPTS=false
```

This uses the expression matrix and cell metadata supplied by 10x. It supports cell-level QC, clustering, markers, and spatial polygon plots, but it cannot display individual transcript detections or reaggregate with a different QV threshold or segmentation unit.

### Full transcript mode

```bash
export XENIUM_LOAD_TRANSCRIPTS=true
export XENIUM_QV_THRESHOLD=20
export XENIUM_AGGREGATION_UNIT=cell  # or nucleus
```

This loads transcript points and rebuilds the RNA matrix from transcript-polygon overlaps. Use it when transcript-level plots, a different QV threshold, or nuclear aggregation are required.

For old pre-release datasets with legacy control-probe names:

```bash
export XENIUM_LEGACY_PROBES=true
```

Negative controls are kept separate from the `rna` feature type so they do not enter RNA normalization.

## Optional aligned H&E and IF images

Standard Xenium morphology images are loaded during import. To add post-Xenium images aligned in Xenium Explorer, provide OME-TIFF and affine-alignment CSV pairs:

```bash
export XENIUM_HE_OMETIF=/path/to/post_xenium_he.ome.tif
export XENIUM_HE_ALIGNMENT=/path/to/he_alignment.csv

export XENIUM_IF_OMETIF=/path/to/post_xenium_if.ome.tif
export XENIUM_IF_ALIGNMENT=/path/to/if_alignment.csv
export XENIUM_IF_PAGES=1,2,3
export XENIUM_IF_CHANNEL_NAMES=CD20,HER2,DAPI
```

The image script converts requested pages to regular TIFF files before attaching them. If no paths are provided, step 02 simply forwards the imported object.

## Subcellular zoom

Set a tissue region and an attached image name for transcript/polygon/image overlays:

```bash
export XENIUM_ZOOM_XLIM=1000,2000
export XENIUM_ZOOM_YLIM=-3000,-2000
export XENIUM_ZOOM_IMAGE=dapi
```

These coordinates are dataset-specific. Without them, the script creates a full-field plot and prints a reminder.

## Filtering and clustering

Common environment variables are:

| Variable | Default |
|---|---:|
| `XENIUM_EXPRESSION_THRESHOLD` | `1` |
| `XENIUM_FEATURE_MIN_CELLS` | `3` |
| `XENIUM_CELL_MIN_FEATURES` | `5` |
| `XENIUM_PCA_DIMS` | `10` |
| `XENIUM_RUN_TSNE` | `true` |
| `XENIUM_NEIGHBOR_K` | `10` |
| `XENIUM_LEIDEN_RESOLUTION` | `0.25` |
| `XENIUM_LEIDEN_ITERATIONS` | `100` |

The defaults mirror the tutorial's targeted-panel workflow, including PCA on all retained genes rather than a highly-variable subset. Revisit them for larger panels or substantially different cell counts.

## Marker detection and manual Xenium annotation

Configure genes and marker methods:

```bash
export XENIUM_MARKER_GENES=EPCAM,KRT8,MKI67,PTPRC,CD3D,CD8A,MS4A1,LYZ,DCN,PECAM1
export XENIUM_MARKER_METHODS=scran,gini
export XENIUM_MARKER_TOP_N=5
```

Step 06 writes full marker tables to `results/tables/`. It does not infer biological cell types automatically. To apply reviewed cluster labels, copy and edit `config/xenium_cluster_annotations.example.csv`, then set:

```bash
export XENIUM_CLUSTER_ANNOTATIONS=/path/to/xenium_cluster_annotations.csv
```

The CSV must contain `cluster` and `cell_type` columns.

## Optional scRNA reference, Harmony, and label transfer

This section is dataset-specific and not part of the core runner by default. A reference is appropriate only when its tissue, biological condition, annotation quality, and gene overlap support the intended transfer.

Prepare an annotated reference mapping and configure:

```bash
export GIOTTO_REFERENCE_H5=/path/to/filtered_feature_bc_matrix.h5
export GIOTTO_REFERENCE_ANNOTATIONS=/path/to/reference_cluster_annotations.csv
export GIOTTO_REFERENCE_LABEL=cell_type

export GIOTTO_REFERENCE_MITO_MAX=50
export GIOTTO_REFERENCE_RUN_SCRUBLET=true
export GIOTTO_REFERENCE_PCA_DIMS=30
export GIOTTO_REFERENCE_LEIDEN_RESOLUTION=1.5
```

Run steps 07-09 with the core workflow:

```bash
export XENIUM_RUN_REFERENCE_PIPELINE=true
bash run_pipeline.sh
```

Or run them individually after the core pipeline:

```bash
Rscript scripts/07_process_scrna_reference.R
Rscript scripts/08_integrate_harmony.R
Rscript scripts/09_transfer_labels.R
```

Step 08 saves the shared-gene list and both unintegrated and Harmony UMAPs. Inspect these before interpreting transferred labels. Step 09 saves a per-cell label and probability table; low-confidence or biologically implausible assignments should not be used downstream without review.

## Outputs and restart points

```text
results/
├── converted_images/   # TIFF exports from optional OME-TIFF inputs
├── figures/            # plots saved by Giotto
├── objects/            # numbered RDS checkpoints after each stage
└── tables/             # QC, markers, shared features, and transferred labels
```

Each script loads the previous numbered object and writes its own checkpoint. A failed downstream step can therefore be rerun without repeating the full import.

## Important limitations

- The original Giotto vignette uses a legacy pre-release breast-cancer dataset. File names, probe conventions, image types, and segmentation behavior may differ in current Xenium outputs.
- Cluster annotations in the example CSV files are placeholders, not biological recommendations.
- Segmentation quality directly affects transcript assignment. Always inspect cell and nuclear polygons in relevant tissue regions.
- Harmony can reduce technical separation but can also obscure real biology. Compare integrated and unintegrated embeddings.
- Label transfer is a hypothesis-generating annotation step. Review shared features, reference labels, transfer probabilities, marker expression, and spatial plausibility.

## License

MIT
