#!/bin/bash
# Unified analyzer runner for instrument validation datasets.

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
VALIDATION_DIR="$(cd "$SCRIPT_DIR/.." && pwd)"
REPO_ROOT="$(cd "$VALIDATION_DIR/.." && pwd)"
INSTRUMENT_DIR="$VALIDATION_DIR/scripts/instrumentation"
CONFIG_DIR="$VALIDATION_DIR/configs/instruments"

PYTHON_BIN=${PYTHON_BIN:-python3}
RESULTS_ROOT_OVERRIDE=""

AVAILABLE_ANALYZERS=(
  ims-drift
  ims-en
  fticr
  lqit
  orbitrap
  tof
  quadrupole
)

declare -A ANALYZER_LABELS=(
  [ims-drift]="IMS drift velocity"
  [ims-en]="IMS E/N mapping"
  [fticr]="FTICR cyclotron frequencies"
  [lqit]="LQIT confinement and secular frequency"
  [orbitrap]="Orbitrap axial frequency"
  [tof]="TOF flight time scaling"
  [quadrupole]="Quadrupole stability map"
)

declare -A ANALYZER_SCRIPTS=(
  [ims-drift]="$INSTRUMENT_DIR/ims/analyze_ims_drift.py"
  [ims-en]="$INSTRUMENT_DIR/ims/analyze_ims_EN_mapping.py"
  [fticr]="$INSTRUMENT_DIR/fticr/analyze_fticr_frequencies.py"
  [lqit]="$INSTRUMENT_DIR/lqit/analyze_lqit_stability.py"
  [orbitrap]="$INSTRUMENT_DIR/orbitrap/analyze_orbitrap_frequency.py"
  [tof]="$INSTRUMENT_DIR/tof/analyze_tof_flight_time.py"
  [quadrupole]="$INSTRUMENT_DIR/quadrupole/analyze_quad_stability.py"
)

declare -A ANALYZER_INSTRUMENT_DIR=(
  [ims-drift]="ims"
  [ims-en]="ims"
  [fticr]="fticr"
  [lqit]="lqit"
  [orbitrap]="orbitrap"
  [tof]="tof"
  [quadrupole]="quadrupole"
)

QUAD_CONFIG_DIR="$CONFIG_DIR/quadrupole"
RESULT_ROOTS=()
SELECTED_ANALYZERS=()

print_usage() {
  cat <<EOF
Usage: $(basename "$0") [options] [analyzer ...]

Options:
  --python PATH        Python interpreter to use (default: $PYTHON_BIN)
  --results-root DIR   Prepend DIR to instrument results search path
  --list               Show available analyzers and exit
  -h, --help           Show this help message

Analyzers:
EOF
  list_analyzers "  "
}

list_analyzers() {
  local prefix="${1:-}"
  for key in "${AVAILABLE_ANALYZERS[@]}"; do
    printf '%s- %-14s %s\n' "$prefix" "$key" "${ANALYZER_LABELS[$key]}"
  done
  printf '%s- %-14s %s\n' "$prefix" "ims" "Shortcut for ims-drift + ims-en"
  printf '%s- %-14s %s\n' "$prefix" "all" "Run every analyzer (default)"
}

dedupe_array() {
  local -n _arr=$1
  local -A _seen=()
  local _unique=()
  for item in "${_arr[@]:-}"; do
    if [[ -n "$item" && -z "${_seen[$item]:-}" ]]; then
      _unique+=("$item")
      _seen[$item]=1
    fi
  done
  _arr=("${_unique[@]:-}")
}

abs_path() {
  local target="$1"
  if [[ -z "$target" ]]; then
    return 1
  fi
  if [[ "$target" = /* ]]; then
    printf '%s\n' "$target"
  else
    printf '%s/%s\n' "$(pwd)" "$target"
  fi
}

add_selection() {
  local token="${1,,}"
  case "$token" in
    all)
      SELECTED_ANALYZERS+=("${AVAILABLE_ANALYZERS[@]}")
      ;;
    ims|ion-mobility|ionmobility)
      SELECTED_ANALYZERS+=("ims-drift" "ims-en")
      ;;
    ims-drift|ims_drift|drift)
      SELECTED_ANALYZERS+=("ims-drift")
      ;;
    ims-en|ims_en|en|mapping)
      SELECTED_ANALYZERS+=("ims-en")
      ;;
    fticr|ft-icr)
      SELECTED_ANALYZERS+=("fticr")
      ;;
    lqit)
      SELECTED_ANALYZERS+=("lqit")
      ;;
    orbitrap)
      SELECTED_ANALYZERS+=("orbitrap")
      ;;
    tof|time-of-flight|timeofflight)
      SELECTED_ANALYZERS+=("tof")
      ;;
    quadrupole|quad|quadrupole-stability)
      SELECTED_ANALYZERS+=("quadrupole")
      ;;
    *)
      echo "Unknown analyzer: $1" >&2
      exit 1
      ;;
  esac
}

while [[ $# -gt 0 ]]; do
  case "$1" in
    --python)
      [[ $# -lt 2 ]] && print_usage && exit 1
      PYTHON_BIN="$2"
      shift 2
      ;;
    --results-root)
      [[ $# -lt 2 ]] && print_usage && exit 1
      RESULTS_ROOT_OVERRIDE="$2"
      shift 2
      ;;
    --list)
      list_analyzers "  "
      exit 0
      ;;
    -h|--help)
      print_usage
      exit 0
      ;;
    *)
      add_selection "$1"
      shift
      ;;
  esac
done

if [[ ${#SELECTED_ANALYZERS[@]} -eq 0 ]]; then
  SELECTED_ANALYZERS=("${AVAILABLE_ANALYZERS[@]}")
fi

dedupe_array SELECTED_ANALYZERS

if ! command -v "$PYTHON_BIN" >/dev/null 2>&1; then
  echo "Error: Python interpreter '$PYTHON_BIN' not found" >&2
  exit 1
fi

build_results_roots() {
  RESULT_ROOTS=()
  if [[ -n "$RESULTS_ROOT_OVERRIDE" ]]; then
    RESULT_ROOTS+=("$(abs_path "$RESULTS_ROOT_OVERRIDE")")
  fi
  RESULT_ROOTS+=(
    "$VALIDATION_DIR/results/instruments"
    "$VALIDATION_DIR/results/v1.0_test/instruments"
    "$REPO_ROOT/results/v1.0_test/instruments"
  )
  dedupe_array RESULT_ROOTS
}

build_results_roots

resolve_results_dir() {
  local instrument="$1"
  for root in "${RESULT_ROOTS[@]}"; do
    local candidate="$root/$instrument"
    if [[ -d "$candidate" ]]; then
      printf '%s' "$candidate"
      return 0
    fi
  done
  return 1
}

run_python() {
  local script="$1"
  shift
  if [[ ! -f "$script" ]]; then
    echo "Error: Analyzer script not found: $script" >&2
    return 1
  fi
  "$PYTHON_BIN" "$script" "$@"
}

run_analysis() {
  local key="$1"
  local label="${ANALYZER_LABELS[$key]}"
  local instrument_dir="${ANALYZER_INSTRUMENT_DIR[$key]}"

  echo ""
  echo "============================================================"
  echo "$label"
  echo "============================================================"

  local data_dir
  if ! data_dir=$(resolve_results_dir "$instrument_dir"); then
    echo "✗ No results found for instrument '$instrument_dir'." >&2
    return 1
  fi
  echo "Data directory: $data_dir"

  case "$key" in
    ims-drift|lqit|orbitrap|tof)
      mapfile -t files < <(find "$data_dir" -maxdepth 1 -type f -name '*.h5' | sort)
      if [[ ${#files[@]} -eq 0 ]]; then
        echo "✗ No .h5 files detected under $data_dir" >&2
        return 1
      fi
      if [[ "$key" == "lqit" ]]; then
        echo "Found ${#files[@]} trajectory files"
      fi
      if ! run_python "${ANALYZER_SCRIPTS[$key]}" "${files[@]}"; then
        return 1
      fi
      ;;
    ims-en)
      if ! run_python "${ANALYZER_SCRIPTS[$key]}" "$data_dir"; then
        return 1
      fi
      ;;
    fticr)
      if ! run_python "${ANALYZER_SCRIPTS[$key]}" "$data_dir"; then
        return 1
      fi
      ;;
    quadrupole)
      if [[ ! -d "$QUAD_CONFIG_DIR" ]]; then
        echo "✗ Config directory not found: $QUAD_CONFIG_DIR" >&2
        return 1
      fi
      echo "Config directory: $QUAD_CONFIG_DIR"
      if ! run_python "${ANALYZER_SCRIPTS[$key]}" \
        --results-dir "$data_dir" \
        --config-dir "$QUAD_CONFIG_DIR"; then
        return 1
      fi
      ;;
    *)
      echo "✗ Internal error: unsupported analyzer '$key'" >&2
      return 1
      ;;
  esac

  echo "✓ Completed $label"
  return 0
}

TOTAL=0
FAILED=0
PASSED=0
FAILED_KEYS=()

for key in "${SELECTED_ANALYZERS[@]}"; do
  TOTAL=$((TOTAL + 1))
  if run_analysis "$key"; then
    PASSED=$((PASSED + 1))
  else
    FAILED=$((FAILED + 1))
    FAILED_KEYS+=("$key")
  fi
 done

echo ""
echo "============================================================"
echo "Instrument analysis summary"
echo "============================================================"
echo "Requested analyzers : ${SELECTED_ANALYZERS[*]}"
echo "Results roots       : ${RESULT_ROOTS[*]}"
echo "Succeeded            : $PASSED"
echo "Failed               : $FAILED"
echo "============================================================"

if [[ $FAILED -gt 0 ]]; then
  echo "Failures: ${FAILED_KEYS[*]}" >&2
  exit 1
fi

exit 0
