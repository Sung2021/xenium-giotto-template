#!/usr/bin/env bash
set -euo pipefail

Rscript scripts/01_import_xenium.R
Rscript scripts/02_process_cluster.R
Rscript scripts/03_visualize.R
