#!/usr/bin/env bash
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO_ROOT="$(cd "${SCRIPT_DIR}/../.." && pwd)"

PROJECT_DIR="${PR141_THERMO_PROJECT_DIR:-${REPO_ROOT}}"
JULIA="${JULIA:-julia}"
EXECUTE="${RUN_PR141_EXTENDED_THERMO_EXECUTE:-1}"
RESULTS_DIR="${PR141_THERMO_RESULTS_DIR:-${SCRIPT_DIR}/results}"

echo "PROJECT_DIR=${PROJECT_DIR}"
echo "RESULTS_DIR=${RESULTS_DIR}"
echo "Spec file: ${SCRIPT_DIR}/thermodynamics_case_specs.toml"

if [[ ! -f "${PROJECT_DIR}/Project.toml" ]]; then
  echo "Missing ClimaSeaIce project at ${PROJECT_DIR}" >&2
  exit 2
fi

if [[ "${EXECUTE}" != "1" ]]; then
  echo "Dry run only. Set RUN_PR141_EXTENDED_THERMO_EXECUTE=1 to execute."
  exit 0
fi

PR141_THERMO_RESULTS_DIR="${RESULTS_DIR}" \
  "${JULIA}" --project="${PROJECT_DIR}" \
  "${SCRIPT_DIR}/analyze_thermodynamics_comparison.jl"

PR141_THERMO_RESULTS_DIR="${RESULTS_DIR}" \
  "${JULIA}" --project="${PROJECT_DIR}" \
  "${SCRIPT_DIR}/plot_thermodynamics_comparison.jl"
