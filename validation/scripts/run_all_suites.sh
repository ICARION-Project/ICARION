#!/bin/bash
# Run physics + instruments + performance suites sequentially into a single run directory.

set -euo pipefail

ORIGINAL_ARGS=("$@")

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
VALIDATION_DIR="$(cd "$SCRIPT_DIR/.." && pwd)"
REPO_ROOT="$(cd "$VALIDATION_DIR/.." && pwd)"

RUN_ID=${RUN_ID:-""}
RUN_DIR=${RUN_DIR:-""}

# Forwarded options
PYTHON_BIN=${PYTHON_BIN:-python3}
THERMALIZATION_MODE=${THERMALIZATION_MODE:-"full"}
CONFIG_JOBS=${CONFIG_JOBS:-8}

INSTRUMENT_SUITE_JOBS=${INSTRUMENT_SUITE_JOBS:-1}
INSTRUMENT_JOBS=${INSTRUMENT_JOBS:-""}
INSTRUMENT_THREADS=${INSTRUMENT_THREADS:-""}

print_usage() {
  cat <<'EOF'
Usage: ./run_all_suites.sh [options]

Options:
  --run-id ID                Run identifier (default: YYYYmmdd_HHMMSS)
  --run-dir PATH             Output directory for this run (default: validation/runs/<run-id>)

  --python PATH              Python interpreter for physics suite (default: $PYTHON_BIN)
  --thermalization-mode MODE Physics thermalization mode (quick|subset|full; default: full)
  --config-jobs N            Physics suite parallelism where applicable (default: 8)

  -J, --suite-jobs N         Instruments: number of instruments to run concurrently (default: 1)
  -j, --jobs N               Instruments: parallel jobs per instrument runner (optional)
  -t, --threads N            Instruments: threads per simulation (optional)

  -h, --help                 Show this help
EOF
}

while [[ $# -gt 0 ]]; do
  case "$1" in
    --run-id)
      RUN_ID="$2"; shift 2 ;;
    --run-dir)
      RUN_DIR="$2"; shift 2 ;;

    --python)
      PYTHON_BIN="$2"; shift 2 ;;
    --thermalization-mode)
      THERMALIZATION_MODE="$2"; shift 2 ;;
    --config-jobs)
      CONFIG_JOBS="$2"; shift 2 ;;

    -J|--suite-jobs)
      INSTRUMENT_SUITE_JOBS="$2"; shift 2 ;;
    -j|--jobs)
      INSTRUMENT_JOBS="$2"; shift 2 ;;
    -t|--threads)
      INSTRUMENT_THREADS="$2"; shift 2 ;;

    -h|--help)
      print_usage; exit 0 ;;
    *)
      echo "Unknown argument: $1" >&2
      print_usage >&2
      exit 1
      ;;
  esac
done

if [[ -z "$RUN_ID" ]]; then
  RUN_ID=$(date +%Y%m%d_%H%M%S)
fi
if [[ -z "$RUN_DIR" ]]; then
  RUN_DIR="$VALIDATION_DIR/runs/$RUN_ID"
fi

RUN_DIR_ABS=$(cd "$(dirname "$RUN_DIR")" && pwd)/"$(basename "$RUN_DIR")"
mkdir -p "$RUN_DIR_ABS/logs" "$RUN_DIR_ABS/figures" "$RUN_DIR_ABS/results"

# Stable top-level manifest for the overall run (suites write manifest.<suite>.json)
if command -v python3 >/dev/null 2>&1; then
  cmd=$(printf '%q ' "$0" "${ORIGINAL_ARGS[@]}")
  SUITE_NAME="full" \
  RUN_ID="$RUN_ID" \
  RUN_DIR="$RUN_DIR_ABS" \
  REPO_ROOT="$REPO_ROOT" \
  VALIDATION_DIR="$VALIDATION_DIR" \
  COMMAND_LINE="$cmd" \
  python3 - "$RUN_DIR_ABS/manifest.json" <<'PY'
import json
import os
import subprocess
import sys
from datetime import datetime, timezone

out_file = sys.argv[1]

def _cmd(argv):
  try:
    return subprocess.check_output(argv, stderr=subprocess.DEVNULL, text=True).strip()
  except Exception:
    return None

repo_root = os.environ.get("REPO_ROOT")

data = {
  "suite": os.environ.get("SUITE_NAME"),
  "run_id": os.environ.get("RUN_ID"),
  "run_dir": os.environ.get("RUN_DIR"),
  "timestamp_utc": datetime.now(timezone.utc).isoformat(),
  "repo_root": repo_root,
  "validation_dir": os.environ.get("VALIDATION_DIR"),
  "command": os.environ.get("COMMAND_LINE"),
}

data["git_commit"] = _cmd(["git", "-C", repo_root, "rev-parse", "HEAD"]) if repo_root else None
git_status = _cmd(["git", "-C", repo_root, "status", "--porcelain"]) if repo_root else None
data["git_dirty"] = bool(git_status)

os.makedirs(os.path.dirname(out_file), exist_ok=True)
with open(out_file, "w", encoding="utf-8") as handle:
  json.dump(data, handle, indent=2)
  handle.write("\n")
PY
fi

echo "============================================================"
echo "ICARION Full Validation Run"
echo "============================================================"
echo "Repo root : $REPO_ROOT"
echo "Run ID    : $RUN_ID"
echo "Run dir   : $RUN_DIR_ABS"
echo "Started   : $(date)"
echo "============================================================"
echo ""

failures=()

run_step() {
  local label="$1"; shift
  local log_file="$RUN_DIR_ABS/logs/${label}.log"

  echo "============================================================"
  echo "$label"
  echo "Log: $log_file"
  echo "============================================================"

  # Run from repo root to avoid any relative output surprises.
  (cd "$REPO_ROOT" && "$@") 2>&1 | tee "$log_file"
}

# Physics
if ! run_step "suite_physics" \
  bash "$VALIDATION_DIR/scripts/run_physics_suite.sh" \
    --python "$PYTHON_BIN" \
    --thermalization-mode "$THERMALIZATION_MODE" \
    --config-jobs "$CONFIG_JOBS" \
    --run-id "$RUN_ID" \
    --run-dir "$RUN_DIR_ABS"; then
  failures+=("suite_physics")
fi

# Instruments
instrument_args=(
  bash "$VALIDATION_DIR/scripts/run_instrument_suite.sh"
  --run-id "$RUN_ID"
  --run-dir "$RUN_DIR_ABS"
  --suite-jobs "$INSTRUMENT_SUITE_JOBS"
)
if [[ -n "$INSTRUMENT_JOBS" ]]; then
  instrument_args+=(--jobs "$INSTRUMENT_JOBS")
fi
if [[ -n "$INSTRUMENT_THREADS" ]]; then
  instrument_args+=(--threads "$INSTRUMENT_THREADS")
fi
if ! run_step "suite_instruments" "${instrument_args[@]}"; then
  failures+=("suite_instruments")
fi

# Performance
if ! run_step "suite_performance" \
  bash "$VALIDATION_DIR/scripts/performance/run_performance_suite.sh" \
    --run-id "$RUN_ID" \
    --run-dir "$RUN_DIR_ABS"; then
  failures+=("suite_performance")
fi

echo ""
echo "============================================================"
echo "All suites finished"
echo "Run dir: $RUN_DIR_ABS"
echo "Finished: $(date)"
echo "============================================================"

if [[ ${#failures[@]} -gt 0 ]]; then
  echo ""
  echo "Failures: ${failures[*]}" >&2
  exit 1
fi
