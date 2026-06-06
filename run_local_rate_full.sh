#!/usr/bin/env bash
#
# Run the natural-gradient local-rate split pipeline end to end (CPU):
#
#   1. operator grid        -> Lambda_hat   (configs/.../smoke.yaml)
#   2. linearized-rate grid -> gamma_loc + Lambda_hat, saves slow eigenvectors
#                              (configs/.../smoke.yaml)
#   3. flow validation      -> reuses the eigenvectors from stage 2
#                              (configs/.../flow_validation.yaml)
#   4. figures              -> fig1..fig6 (PNG + PDF)
#
# This drives the CPU split-stage runners on the smoke grid as a quick, fully
# local end-to-end check. The FINAL production evidence is produced separately by
# the joint GPU runner on configs/.../gpu_lowdim_operator_full.yaml (see README).
#
# Usage:
#   ./run_local_rate_full.sh              # fresh run (overwrites existing results)
#   RESUME=1 ./run_local_rate_full.sh     # skip stages that already completed
#   PYTHON=python ./run_local_rate_full.sh
#
# Tip: run unattended with
#   nohup ./run_local_rate_full.sh >/dev/null 2>&1 &
# (output is also written to the log file printed below).

set -euo pipefail

# Always operate from the repository root (this script's directory).
cd "$(dirname "$0")"

PYTHON="${PYTHON:-python3}"
CFG_DIR="configs/natural_gradient_local_rate"
SCRIPT_DIR="scripts/natural_gradient_local_rate"
OUT_DIR="outputs/natural_gradient_local_rate"

# Fail fast (with guidance) if the chosen interpreter lacks the dependencies,
# rather than crashing several lines into stage 1.
if ! "$PYTHON" -c "import yaml, numpy, scipy, matplotlib, pandas" 2>/dev/null; then
  echo "ERROR: '$PYTHON' is missing required packages (pyyaml/numpy/scipy/matplotlib/pandas)." >&2
  echo "       ($("$PYTHON" --version 2>&1))" >&2
  echo "Fix it with either:" >&2
  echo "  1) point the script at an interpreter that has them, e.g.:" >&2
  echo "       PYTHON=/Users/ryanlai/miniconda3/bin/python3 $0" >&2
  echo "  2) install them into this interpreter:" >&2
  echo "       $PYTHON -m pip install -r requirements.txt" >&2
  exit 1
fi

# Fresh run overwrites prior results (e.g. leftover smoke data); RESUME=1 skips
# any stage whose results already exist.
OVERWRITE_FLAG="--overwrite"
if [[ "${RESUME:-0}" == "1" ]]; then
  OVERWRITE_FLAG=""
fi

mkdir -p outputs/logs
LOG="outputs/logs/local_rate_full_$(date +%Y%m%d_%H%M%S).log"

# Mirror all output to the console and the log file.
exec > >(tee -a "$LOG") 2>&1

echo "######################################################################"
echo "# Natural-gradient local-rate FULL run"
echo "#   python   : $($PYTHON --version 2>&1)"
echo "#   mode     : $([[ -n "$OVERWRITE_FLAG" ]] && echo 'fresh (overwrite)' || echo 'resume')"
echo "#   log      : $LOG"
echo "#   started  : $(date)"
echo "######################################################################"

run_stage () {
  local title="$1"; shift
  echo
  echo "======================================================================"
  echo "STAGE: $title"
  echo "  cmd: $*"
  echo "  at : $(date)"
  echo "======================================================================"
  local t0=$SECONDS
  "$@"
  echo "  [stage done in $((SECONDS - t0))s]"
}

START=$SECONDS

run_stage "1/4  operator grid (Lambda_hat)" \
  "$PYTHON" "$SCRIPT_DIR/run_operator_grid.py" \
  --config "$CFG_DIR/smoke.yaml" $OVERWRITE_FLAG

run_stage "2/4  linearized-rate grid (gamma_loc + eigenvectors)" \
  "$PYTHON" "$SCRIPT_DIR/run_linearized_rate_grid.py" \
  --config "$CFG_DIR/smoke.yaml" $OVERWRITE_FLAG

run_stage "3/4  flow validation" \
  "$PYTHON" "$SCRIPT_DIR/run_flow_validation.py" \
  --config "$CFG_DIR/smoke.yaml" $OVERWRITE_FLAG

run_stage "4/4  figures" \
  "$PYTHON" "$SCRIPT_DIR/plot_results.py" \
  --input "$OUT_DIR" --outdir "$OUT_DIR/figures"

echo
echo "######################################################################"
echo "# ALL DONE in $((SECONDS - START))s"
echo "#   operator grid : $OUT_DIR/operator_grid/{results_long,summary}.csv"
echo "#   linearized    : $OUT_DIR/linearized_rate_grid/{results_long,summary}.csv (+ eigenvectors/)"
echo "#   flow          : $OUT_DIR/flow_validation/{summary,trajectories}.csv"
echo "#   figures       : $OUT_DIR/figures/  (fig1..fig6 PNG+PDF)"
echo "#   log           : $LOG"
echo "######################################################################"
