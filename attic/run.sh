#!/bin/bash
# run.sh -- LRZ cluster launch script for Registration Benchmark v2
#
# Usage:
#   bash run.sh [study] [n_cores]
#
# Arguments:
#   study:   A, B, C, D, pilot, all (default: pilot)
#   n_cores: number of parallel cores (default: 6)

set -euo pipefail

STUDY="${1:-pilot}"
NCORES="${2:-6}"
SCRIPT_DIR="$(cd "$(dirname "$0")" && pwd)"

echo "=== Registration Benchmark v2 ==="
echo "Study: ${STUDY}"
echo "Cores: ${NCORES}"
echo "Dir:   ${SCRIPT_DIR}"
echo "Start: $(date)"
echo ""

# Change to package root for devtools::load_all()
cd "${SCRIPT_DIR}/../.."

# Run
Rscript "${SCRIPT_DIR}/sim-run.R" "${STUDY}" "${NCORES}" 2>&1 | \
  tee "${SCRIPT_DIR}/benchmark_${STUDY}.log"

echo ""
echo "=== Done: $(date) ==="
