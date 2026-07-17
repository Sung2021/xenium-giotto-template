#!/usr/bin/env bash
set -euo pipefail

core_steps=(
  scripts/01_import_xenium.R
  scripts/02_attach_images.R
  scripts/03_qc_aggregate.R
  scripts/04_process_cluster.R
  scripts/05_visualize_subcellular.R
  scripts/06_detect_markers.R
)

for step in "${core_steps[@]}"; do
  Rscript "$step"
done

case "${XENIUM_RUN_REFERENCE_PIPELINE:-false}" in
  1|true|TRUE|yes|YES)
    reference_steps=(
      scripts/07_process_scrna_reference.R
      scripts/08_integrate_harmony.R
      scripts/09_transfer_labels.R
    )
    for step in "${reference_steps[@]}"; do
      Rscript "$step"
    done
    ;;
esac
