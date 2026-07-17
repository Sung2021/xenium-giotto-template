# Xenium Giotto analysis template

Reusable R project for importing, processing, clustering, and visualizing a 10x Genomics Xenium dataset with Giotto Suite.

This repository follows the main workflow in Giotto's [Xenium breast cancer tutorial](https://giottosuite.com/articles/xenium_breast_cancer.html), but it is organized as a dataset-independent template. It does not include the tutorial data or copy its dataset-specific annotations.

## What is included

1. Giotto and Python environment setup
2. Xenium import with full or low-memory mode
3. Optional transcript-to-cell aggregation
4. Filtering, normalization, PCA, UMAP, and Leiden clustering
5. Spatial, cluster, and marker plots
6. An optional scaffold for single-cell reference label transfer

## Requirements

- R 4.4 or newer is recommended
- Python 3.10 or another Giotto-compatible Python environment
- A standard Xenium output directory containing files such as `transcripts.parquet`, `cells.parquet`, `cell_boundaries.parquet`, and `gene_panel.json`
- Substantially more memory when individual transcripts are loaded. The source tutorial reports roughly 70 GB peak RAM with transcripts and 5 GB without them for its full dataset.

## Quick start

```bash
git clone https://github.com/Sung2021/xenium-giotto-template.git
cd xenium-giotto-template

export XENIUM_DATA_DIR=/absolute/path/to/xenium/outs
export XENIUM_RESULTS_DIR=/absolute/path/to/results

Rscript scripts/00_setup.R
Rscript scripts/01_import_xenium.R
Rscript scripts/02_process_cluster.R
Rscript scripts/03_visualize.R
```

Run all analysis steps after setup with:

```bash
bash run_pipeline.sh
```

`scripts/00_setup.R` is a one-time manual environment setup step. `run_pipeline.sh` intentionally runs only the repeatable core analysis steps: import, processing/clustering, and visualization. It does not rerun setup or the optional reference label-transfer scaffold.

## Choose an import mode

The default is the lower-memory mode:

```bash
export XENIUM_LOAD_TRANSCRIPTS=false
```

It loads the 10x expression matrix and cell metadata, but skips individual transcript detections. This supports cell-level processing and clustering, but not transcript-point plots.

To load transcript detections and rebuild the expression matrix from transcript/polygon overlaps:

```bash
export XENIUM_LOAD_TRANSCRIPTS=true
```

This gives more control over the QV threshold and segmentation unit, but may require an HPC node.

## Configuration

Edit [config/project_config.R](config/project_config.R) or set environment variables:

| Variable | Default | Purpose |
|---|---:|---|
| `XENIUM_DATA_DIR` | `data/outs` | Xenium output directory |
| `XENIUM_RESULTS_DIR` | `results` | Analysis outputs |
| `XENIUM_LOAD_TRANSCRIPTS` | `false` | Load individual transcript detections |
| `XENIUM_QV_THRESHOLD` | `20` | Transcript quality threshold |
| `XENIUM_PCA_DIMS` | `10` | Dimensions used for UMAP and neighbors |
| `XENIUM_LEIDEN_RESOLUTION` | `0.25` | Leiden clustering resolution |
| `GIOTTO_PYTHON_PATH` | unset | Optional explicit Python executable |

Project data and results are ignored by Git. Keep raw Xenium data outside the repository when possible.

## Outputs

The pipeline writes intermediate Giotto objects to `results/objects/` and figures to `results/figures/`. Each step saves an object so a failed downstream step can be restarted without repeating import.

## Optional, dataset-specific reference label transfer

This is not part of the core pipeline. The source tutorial integrates Xenium with paired scRNA-seq data using Harmony and transfers cell-type labels, but that step depends on the reference dataset, its annotations, and feature overlap with the Xenium panel.

`scripts/04_label_transfer_template.R` is therefore a guarded scaffold, not a ready-to-run default. Supply a processed and annotated Giotto single-cell object, then verify reference quality, shared features, batch effects, transferred-label probabilities, and biological plausibility before using the labels downstream.

## Notes

- The referenced breast cancer vignette is explicitly marked as a legacy pre-release dataset. This template uses the standard Xenium output directory as its primary input and keeps old probe-splitting rules optional.
- Marker genes in `config/project_config.R` are examples. Replace them with markers appropriate for your tissue, panel, and biological question.
- Check plots and QC distributions before accepting clusters or transferred labels.

## License

MIT
