#!/usr/bin/env bash
set -euo pipefail

repo_root="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
cd "$repo_root"

if ! command -v Rscript >/dev/null 2>&1; then
  echo "Rscript was not found on PATH. Install R before running the pipeline." >&2
  exit 127
fi

core_steps=(
  scripts/01_import_xenium.R
  scripts/02_attach_images.R
  scripts/03_qc_aggregate.R
  scripts/04_process_cluster.R
  scripts/05_visualize_subcellular.R
  scripts/06_detect_markers.R
)

run_steps() {
  local step
  for step in "$@"; do
    if [[ ! -f "$step" ]]; then
      echo "Pipeline step not found: $repo_root/$step" >&2
      exit 1
    fi
  done
  for step in "$@"; do
    echo "==> Running $step"
    Rscript "$step"
  done
}

for required_file in config/project_config.R scripts/helpers.R; do
  if [[ ! -f "$required_file" ]]; then
    echo "Required pipeline file not found: $repo_root/$required_file" >&2
    exit 1
  fi
done

run_steps "${core_steps[@]}"

case "${XENIUM_RUN_REFERENCE_PIPELINE:-false}" in
  1|true|TRUE|yes|YES)
    reference_steps=(
      scripts/07_process_scrna_reference.R
      scripts/08_integrate_harmony.R
      scripts/09_transfer_labels.R
    )
    run_steps "${reference_steps[@]}"
    ;;
  0|false|FALSE|no|NO|"")
    ;;
  *)
    echo "XENIUM_RUN_REFERENCE_PIPELINE must be true or false." >&2
    exit 2
    ;;
esac
