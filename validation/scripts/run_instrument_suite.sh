#!/bin/bash
# Sequentially (or concurrently) execute the instrument validation runners.

set -euo pipefail

ORIGINAL_ARGS=("$@")

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
VALIDATION_DIR="$(cd "$SCRIPT_DIR/.." && pwd)"
RUNNER="$VALIDATION_DIR/scripts/run_instrument_tests.sh"
DEFAULT_CONFIG_ROOT="$VALIDATION_DIR/configs/instruments"
DEFAULT_OUTPUT_ROOT="$VALIDATION_DIR/results/v1.0.0_test/instruments"
DEFAULT_INSTRUMENTS=(ims fticr lqit orbitrap tof quadrupole)

SELECTED_INSTRUMENTS=()
JOBS=""
THREADS=""
BINARY=""
CONFIG_ROOT=""
OUTPUT_ROOT=""
SUITE_JOBS=1
RUN_ID=${RUN_ID:-""}
RUN_DIR=${RUN_DIR:-""}
BASELINE_OUTPUT=false

print_usage() {
  cat <<'EOF'
Usage: ./run_instrument_suite.sh [options] [instrument ...]

Options:
  -j, --jobs N          Parallel jobs passed to each run_instrument_tests invocation
  -t, --threads N       Threads per simulation forwarded to the runner
  -b, --binary PATH     Override path to icarion_main
  -c, --config-root DIR Use DIR/<instrument> for config overrides
  -o, --output-root DIR Use DIR/<instrument> for output overrides
  --run-id ID         Run identifier (default: YYYYmmdd_HHMMSS)
  --run-dir PATH      Output directory for this run (default: validation/runs/<run-id>)
  --baseline-output   Write into validation/results/v1.0.0_test/instruments (legacy)
  -J, --suite-jobs N    Number of instruments to run concurrently (default: 1)
  --list                Show normalized instrument keys and exit
  -h, --help            Show this help text

If no instruments are specified, the suite runs: ims, fticr, lqit, orbitrap, tof, quadrupole.
EOF
}

list_instruments() {
  cat <<'EOF'
Available instruments:
  ims          Ion Mobility Spectrometry (delegates to IMS runner)
  fticr        Fourier Transform Ion Cyclotron Resonance
  lqit         Linear Quadrupole Ion Trap
  orbitrap     Orbitrap axial frequency suite
  tof          Time-of-flight mass spectrometer
  quadrupole   Quadrupole stability map
EOF
}

require_positive_integer() {
  local value="$1"
  local label="$2"
  if ! [[ "$value" =~ ^[1-9][0-9]*$ ]]; then
    echo "Error: $label must be a positive integer (got '$value')." >&2
    exit 1
  fi
}

abs_path() {
  local target="$1"
  if [[ "$target" = /* ]]; then
    printf '%s' "$target"
  else
    printf '%s/%s' "$(pwd)" "$target"
  fi
}

normalize_instrument() {
  local token="${1,,}"
  case "$token" in
    ims|ion-mobility|ionmobility) echo "ims" ;;
    fticr|ft-icr) echo "fticr" ;;
    lqit) echo "lqit" ;;
    orbitrap) echo "orbitrap" ;;
    tof|time-of-flight|timeofflight) echo "tof" ;;
    quadrupole|quad|quadrupole_stability) echo "quadrupole" ;;
    all) echo "all" ;;
    *) return 1 ;;
  esac
}

dedupe_array() {
  local -n arr=$1
  local -A seen=()
  local unique=()
  for item in "${arr[@]:-}"; do
    if [[ -n "$item" && -z "${seen[$item]:-}" ]]; then
      unique+=("$item")
      seen[$item]=1
    fi
  done
  arr=("${unique[@]:-}")
}

while [[ $# -gt 0 ]]; do
  case "$1" in
    -j|--jobs)
      [[ $# -lt 2 ]] && print_usage && exit 1
      JOBS="$2"
      shift 2
      ;;
    -t|--threads)
      [[ $# -lt 2 ]] && print_usage && exit 1
      THREADS="$2"
      shift 2
      ;;
    -b|--binary)
      [[ $# -lt 2 ]] && print_usage && exit 1
      BINARY="$2"
      shift 2
      ;;
    -c|--config-root)
      [[ $# -lt 2 ]] && print_usage && exit 1
      CONFIG_ROOT="$2"
      shift 2
      ;;
    -o|--output-root)
      [[ $# -lt 2 ]] && print_usage && exit 1
      OUTPUT_ROOT="$2"
      shift 2
      ;;
    --run-id)
      [[ $# -lt 2 ]] && print_usage && exit 1
      RUN_ID="$2"
      shift 2
      ;;
    --run-dir)
      [[ $# -lt 2 ]] && print_usage && exit 1
      RUN_DIR="$2"
      shift 2
      ;;
    --baseline-output)
      BASELINE_OUTPUT=true
      shift
      ;;
    -J|--suite-jobs)
      [[ $# -lt 2 ]] && print_usage && exit 1
      SUITE_JOBS="$2"
      shift 2
      ;;
    --list)
      list_instruments
      exit 0
      ;;
    -h|--help)
      print_usage
      exit 0
      ;;
    --)
      shift
      break
      ;;
    -* )
      echo "Unknown option: $1" >&2
      print_usage >&2
      exit 1
      ;;
    *)
      if ! normalized=$(normalize_instrument "$1"); then
        echo "Unknown instrument: $1" >&2
        exit 1
      fi
      if [[ "$normalized" == "all" ]]; then
        SELECTED_INSTRUMENTS+=("${DEFAULT_INSTRUMENTS[@]}")
      else
        SELECTED_INSTRUMENTS+=("$normalized")
      fi
      shift
      ;;
  esac
done

if [[ ${#SELECTED_INSTRUMENTS[@]} -eq 0 ]]; then
  SELECTED_INSTRUMENTS=("${DEFAULT_INSTRUMENTS[@]}")
fi

dedupe_array SELECTED_INSTRUMENTS

require_positive_integer "$SUITE_JOBS" "suite jobs"

if [[ ! -x "$RUNNER" ]]; then
  echo "Error: Expected runner not found at $RUNNER" >&2
  exit 1
fi

if [[ -n "$BINARY" ]]; then
  BINARY="$(abs_path "$BINARY")"
fi
if [[ -n "$CONFIG_ROOT" ]]; then
  CONFIG_ROOT="$(abs_path "$CONFIG_ROOT")"
fi
if [[ -n "$OUTPUT_ROOT" ]]; then
  OUTPUT_ROOT="$(abs_path "$OUTPUT_ROOT")"
fi

if [[ -z "$RUN_ID" ]]; then
  RUN_ID=$(date +%Y%m%d_%H%M%S)
fi

if [[ -z "$RUN_DIR" ]]; then
  RUN_DIR="$VALIDATION_DIR/runs/$RUN_ID"
fi

RUN_DIR_ABS=$(cd "$(dirname "$RUN_DIR")" && pwd)/"$(basename "$RUN_DIR")"
mkdir -p "$RUN_DIR_ABS/logs" "$RUN_DIR_ABS/figures" "$RUN_DIR_ABS/results"
export ICARION_VALIDATION_RUN_DIR="$RUN_DIR_ABS"

write_manifest() {
  local out_file="$1"
  if ! command -v python3 >/dev/null 2>&1; then
    return 0
  fi
  local cmd
  cmd=$(printf '%q ' "$0" "${ORIGINAL_ARGS[@]}")

  SUITE_NAME="instruments" \
  RUN_ID="$RUN_ID" \
  RUN_DIR="$RUN_DIR_ABS" \
  REPO_ROOT="$(cd "$VALIDATION_DIR/.." && pwd)" \
  VALIDATION_DIR="$VALIDATION_DIR" \
  ICARION_BIN="${BINARY:-}" \
  COMMAND_LINE="$cmd" \
  python3 - "$out_file" <<'PY'
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
    "icarion_bin": os.environ.get("ICARION_BIN") or None,
}

data["git_commit"] = _cmd(["git", "-C", repo_root, "rev-parse", "HEAD"]) if repo_root else None
git_status = _cmd(["git", "-C", repo_root, "status", "--porcelain"]) if repo_root else None
data["git_dirty"] = bool(git_status)

os.makedirs(os.path.dirname(out_file), exist_ok=True)
with open(out_file, "w", encoding="utf-8") as handle:
    json.dump(data, handle, indent=2)
PY
}

write_manifest "$RUN_DIR_ABS/manifest.instruments.json" || true

if [[ -z "$OUTPUT_ROOT" ]]; then
  if $BASELINE_OUTPUT; then
    OUTPUT_ROOT="$DEFAULT_OUTPUT_ROOT"
  else
    OUTPUT_ROOT="$RUN_DIR_ABS/results/instruments"
  fi
fi

TOTAL=0
PASSED=0
SKIPPED=0
FAILED=0
FAILED_LIST=()
SKIPPED_LIST=()

STATUS_LOG="$(mktemp "$VALIDATION_DIR/suite_status_XXXX.tsv")"
declare -a RUNNING_PIDS=()
ACTIVE_JOBS=0

cleanup() {
  if [[ ${#RUNNING_PIDS[@]} -gt 0 ]]; then
    kill "${RUNNING_PIDS[@]}" 2>/dev/null || true
  fi
  rm -f "$STATUS_LOG"
}
trap cleanup EXIT INT TERM

printf '==============================================\n'
printf 'Instrument Validation Suite\n'
printf '==============================================\n'
printf 'Runner      : %s\n' "$RUNNER"
printf 'Jobs        : %s\n' "${JOBS:-default}"
printf 'Threads     : %s\n' "${THREADS:-default}"
printf 'Suite jobs  : %s\n' "$SUITE_JOBS"
printf 'Binary      : %s\n' "${BINARY:-auto}"
printf 'Config root : %s\n' "${CONFIG_ROOT:-$DEFAULT_CONFIG_ROOT}" 
printf 'Run ID      : %s\n' "$RUN_ID"
printf 'Run dir     : %s\n' "$RUN_DIR_ABS"
printf 'Output root : %s\n' "$OUTPUT_ROOT"
printf 'Targets     : %s\n' "${SELECTED_INSTRUMENTS[*]}"
printf '==============================================\n\n'

run_single() {
  local instrument="$1"

  # Proactively skip instruments without config directories (or without any json configs)
  # so the suite remains reproducible across repos/branches where some instruments
  # are intentionally not shipped.
  local cfg_dir=""
  if [[ -n "${CONFIG_ROOT:-}" ]]; then
    cfg_dir="$CONFIG_ROOT/$instrument"
  else
    cfg_dir="$DEFAULT_CONFIG_ROOT/$instrument"
  fi
  if [[ ! -d "$cfg_dir" ]]; then
    echo "[${instrument}] ⏭ skipped (config dir not found: $cfg_dir)"
    return 2
  fi
  shopt -s nullglob
  local cfg_matches=("$cfg_dir"/*.json)
  shopt -u nullglob
  if [[ ${#cfg_matches[@]} -eq 0 ]]; then
    echo "[${instrument}] ⏭ skipped (no *.json configs in: $cfg_dir)"
    return 2
  fi

  local cmd=("$RUNNER")
  if [[ -n "$JOBS" ]]; then
    cmd+=("-j" "$JOBS")
  fi
  if [[ -n "$THREADS" ]]; then
    cmd+=("-t" "$THREADS")
  fi
  if [[ -n "$BINARY" ]]; then
    cmd+=("-b" "$BINARY")
  fi
  if [[ -n "$CONFIG_ROOT" ]]; then
    cmd+=("-c" "$CONFIG_ROOT/$instrument")
  fi
  if [[ -n "$OUTPUT_ROOT" ]]; then
    cmd+=("-o" "$OUTPUT_ROOT/$instrument")
  fi
  cmd+=("$instrument")

  echo "[${instrument}] starting..."
  if "${cmd[@]}"; then
    echo "[${instrument}] ✅ completed"
    return 0
  elif [[ $? -eq 2 ]]; then
    # (kept for future: runner-level skip)
    echo "[${instrument}] ⏭ skipped"
    return 2
  else
    echo "[${instrument}] ❌ failed" >&2
    return 1
  fi
}

for instrument in "${SELECTED_INSTRUMENTS[@]}"; do
  TOTAL=$((TOTAL + 1))
  (
    set +e
    run_single "$instrument"
    status=$?
    set -e
    printf '%s\t%d\n' "$instrument" "$status" >>"$STATUS_LOG"
    exit "$status"
  ) &
  pid=$!
  RUNNING_PIDS+=("$pid")
  ((ACTIVE_JOBS+=1))

  if [[ $ACTIVE_JOBS -ge $SUITE_JOBS ]]; then
    wait -n || true
    ((ACTIVE_JOBS-=1)) || true
  fi
done

while [[ $ACTIVE_JOBS -gt 0 ]]; do
  wait -n || true
  ((ACTIVE_JOBS-=1)) || true
done

while IFS=$'\t' read -r instrument status; do
  if [[ -z "$instrument" ]]; then
    continue
  fi
  if [[ "$status" -eq 0 ]]; then
    PASSED=$((PASSED + 1))
  elif [[ "$status" -eq 2 ]]; then
    SKIPPED=$((SKIPPED + 1))
    SKIPPED_LIST+=("$instrument")
  else
    FAILED=$((FAILED + 1))
    FAILED_LIST+=("$instrument")
  fi
done <"$STATUS_LOG" || true

rm -f "$STATUS_LOG"

printf '\n==============================================\n'
printf 'Suite summary\n'
printf '==============================================\n'
printf 'Total instruments : %d\n' "$TOTAL"
printf 'Succeeded         : %d\n' "$PASSED"
printf 'Skipped           : %d\n' "$SKIPPED"
printf 'Failed            : %d\n' "$FAILED"
if [[ $SKIPPED -gt 0 ]]; then
  printf 'Skipped list      : %s\n' "${SKIPPED_LIST[*]}"
fi
if [[ $FAILED -gt 0 ]]; then
  printf 'Failures          : %s\n' "${FAILED_LIST[*]}"
fi
printf 'Finished          : %s\n' "$(date)"
printf '==============================================\n'

if [[ $FAILED -gt 0 ]]; then
  exit 1
fi

exit 0
