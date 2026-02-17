#!/bin/bash
# Unified analyzer runner for physics validation datasets.

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}" )" && pwd)"
VALIDATION_DIR="$(cd "$SCRIPT_DIR/.." && pwd)"
REPO_ROOT="$(cd "$VALIDATION_DIR/.." && pwd)"
THERMALIZATION_ANALYZER_DIR="$VALIDATION_DIR/scripts/thermalization"
PHYSICS_ANALYZER_DIR="$VALIDATION_DIR/scripts/physics"

THERMALIZATION_MODULE="analyze_thermalization_complete"
TRANSPORT_ANALYZER="$PHYSICS_ANALYZER_DIR/analyze_transport_drift.py"
SPACECHARGE_ANALYZER="$PHYSICS_ANALYZER_DIR/analyze_spacecharge.py"
REACTIONS_ANALYZER="$PHYSICS_ANALYZER_DIR/analyze_reactions.py"

PYTHON_BIN=${PYTHON_BIN:-python3}
RESULTS_ROOT_OVERRIDE=""
THERM_DIR_OVERRIDE=""
TRANSPORT_DIR_OVERRIDE=""
SPACECHARGE_DIR_OVERRIDE=""
REACTIONS_DIR_OVERRIDE=""
SELECTED_ANALYSES=()

AVAILABLE_ANALYSES=(thermalization transport spacecharge reactions)

declare -A ANALYSIS_LABELS=(
  [thermalization]="Comprehensive thermalization analysis"
  [transport]="Transport & combined drift analysis"
  [spacecharge]="Space charge expansion / IMS analysis"
  [reactions]="Reaction kinetics scenario summaries"
)

print_usage() {
  cat <<'EOF'
Usage: ./run_physics_analysis.sh [options] [analysis ...]

Options:
  --python PATH              Python interpreter (default: $PYTHON_BIN)
  --results-root DIR         Prepend DIR to the auto-detection list
  --thermalization-dir DIR   Explicit directory for thermalization outputs
  --transport-dir DIR        Explicit directory for transport/drift outputs
  --spacecharge-dir DIR      Explicit directory for space charge outputs
  --reactions-dir DIR        Explicit directory for reaction outputs
  --list                     Show available analyses and exit
  -h, --help                 Show this help message

Analyses (default: all): thermalization, transport, spacecharge, reactions
EOF
}

list_analyses() {
  for key in "${AVAILABLE_ANALYSES[@]}"; do
    printf '  %-13s %s\n' "$key" "${ANALYSIS_LABELS[$key]}"
  done
}

normalize_analysis() {
  local token="${1,,}"
  case "$token" in
    all|suite)
      printf '%s\n' "all"
      ;;
    thermalization|therm)
      printf '%s\n' "thermalization"
      ;;
    transport|drift|combined)
      printf '%s\n' "transport"
      ;;
    spacecharge|space-charge|sc)
      printf '%s\n' "spacecharge"
      ;;
    reactions|reaction|kinetics)
      printf '%s\n' "reactions"
      ;;
    *)
      return 1
      ;;
  esac
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
    --thermalization-dir)
      [[ $# -lt 2 ]] && print_usage && exit 1
      THERM_DIR_OVERRIDE="$2"
      shift 2
      ;;
    --transport-dir)
      [[ $# -lt 2 ]] && print_usage && exit 1
      TRANSPORT_DIR_OVERRIDE="$2"
      shift 2
      ;;
    --spacecharge-dir)
      [[ $# -lt 2 ]] && print_usage && exit 1
      SPACECHARGE_DIR_OVERRIDE="$2"
      shift 2
      ;;
    --reactions-dir)
      [[ $# -lt 2 ]] && print_usage && exit 1
      REACTIONS_DIR_OVERRIDE="$2"
      shift 2
      ;;
    --list)
      list_analyses
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
    -*)
      echo "Unknown option: $1" >&2
      print_usage >&2
      exit 1
      ;;
    *)
      if ! normalized=$(normalize_analysis "$1"); then
        echo "Unknown analysis: $1" >&2
        exit 1
      fi
      if [[ "$normalized" == "all" ]]; then
        SELECTED_ANALYSES=("${AVAILABLE_ANALYSES[@]}")
      else
        SELECTED_ANALYSES+=("$normalized")
      fi
      shift
      ;;
  esac
done

if [[ ${#SELECTED_ANALYSES[@]} -eq 0 ]]; then
  SELECTED_ANALYSES=("${AVAILABLE_ANALYSES[@]}")
fi

dedupe_array SELECTED_ANALYSES

if ! command -v "$PYTHON_BIN" >/dev/null 2>&1; then
  echo "Error: Python interpreter '$PYTHON_BIN' not found" >&2
  exit 1
fi

abs_path() {
  local target="$1"
  if [[ "$target" = /* ]]; then
    printf '%s' "$target"
  else
    printf '%s/%s' "$(pwd)" "$target"
  fi
}

latest_match() {
  local pattern="$1"
  local match
  match=$(ls -1dt -- $pattern 2>/dev/null | head -n 1 || true)
  if [[ -n "$match" && -d "$match" ]]; then
    printf '%s' "$match"
    return 0
  fi
  return 1
}

candidate_patterns() {
  local key="$1"
  local -a entries=()
  if [[ -n "$RESULTS_ROOT_OVERRIDE" ]]; then
    case "$key" in
      thermalization)
        entries+=("$RESULTS_ROOT_OVERRIDE/thermalization")
        ;;
      transport)
        entries+=("$RESULTS_ROOT_OVERRIDE/transport" "$RESULTS_ROOT_OVERRIDE/combined_drift")
        ;;
      spacecharge)
        entries+=("$RESULTS_ROOT_OVERRIDE/spacecharge")
        ;;
      reactions)
        entries+=("$RESULTS_ROOT_OVERRIDE/reactions")
        ;;
    esac
  fi
  case "$key" in
    thermalization)
      entries+=(
        "$VALIDATION_DIR/results/v1.0.0_test/physics/thermalization"
        "$VALIDATION_DIR/results/physics/thermalization"
        "$VALIDATION_DIR/results/thermalization"
        "$REPO_ROOT/results/v1.0.0_test/physics/thermalization"
        "$VALIDATION_DIR/scripts/results/thermalization_session_*"
      )
      ;;
    transport)
      entries+=(
        "$VALIDATION_DIR/results/combined_drift"
        "$VALIDATION_DIR/results/physics/combined_drift"
        "$VALIDATION_DIR/results/physics/gas_flow_transport"
        "$VALIDATION_DIR/results/v1.0.0_test/physics/transport/drift"
        "$REPO_ROOT/results/v1.0.0_test/physics/transport/drift"
      )
      ;;
    spacecharge)
      entries+=(
        "$VALIDATION_DIR/results/v1.0.0_test/physics/spacecharge"
        "$VALIDATION_DIR/results/physics/spacecharge"
        "$VALIDATION_DIR/results/spacecharge"
        "$REPO_ROOT/results/v1.0.0_test/physics/spacecharge"
      )
      ;;
    reactions)
      entries+=(
        "$VALIDATION_DIR/results/v1.0.0_test/physics/reactions"
        "$VALIDATION_DIR/results/physics/reactions"
        "$VALIDATION_DIR/results/reactions"
        "$REPO_ROOT/results/v1.0.0_test/physics/reactions"
      )
      ;;
  esac
  printf '%s\n' "${entries[@]}"
}

resolve_data_dir() {
  local key="$1"
  local override=""
  case "$key" in
    thermalization) override="$THERM_DIR_OVERRIDE" ;;
    transport) override="$TRANSPORT_DIR_OVERRIDE" ;;
    spacecharge) override="$SPACECHARGE_DIR_OVERRIDE" ;;
    reactions) override="$REACTIONS_DIR_OVERRIDE" ;;
  esac
  if [[ -n "$override" ]]; then
    printf '%s' "$override"
    return 0
  fi
  local pattern
  while IFS= read -r pattern; do
    [[ -z "$pattern" ]] && continue
    if [[ "$pattern" == *[*?[]* ]]; then
      if match=$(latest_match "$pattern"); then
        printf '%s' "$match"
        return 0
      fi
    elif [[ -d "$pattern" ]]; then
      printf '%s' "$pattern"
      return 0
    fi
  done < <(candidate_patterns "$key")
  return 1
}

run_in_repo() {
  (cd "$REPO_ROOT" && "$@")
}

run_thermalization_analysis() {
  local data_dir="$1"
  local python_path="$THERMALIZATION_ANALYZER_DIR${PYTHONPATH:+:$PYTHONPATH}"
  run_in_repo env PYTHONPATH="$python_path" "$PYTHON_BIN" - "$data_dir" <<'PY'
import sys
from analyze_thermalization_complete import analyze_batch
if len(sys.argv) < 2:
    raise SystemExit("Missing results directory for thermalization analysis")
analyze_batch(sys.argv[1])
PY
}

run_transport_analysis() {
  local data_dir="$1"
  run_in_repo "$PYTHON_BIN" "$TRANSPORT_ANALYZER" "$data_dir"
}

run_spacecharge_analysis() {
  local data_dir="$1"
  run_in_repo "$PYTHON_BIN" "$SPACECHARGE_ANALYZER" "$data_dir"
}

run_reactions_analysis() {
  local data_dir="$1"
  run_in_repo "$PYTHON_BIN" "$REACTIONS_ANALYZER" "$data_dir"
}

run_analysis() {
  local key="$1"
  local label="${ANALYSIS_LABELS[$key]}"
  local data_dir

  echo "============================================================"
  echo "$label"
  echo "============================================================"

  if ! data_dir=$(resolve_data_dir "$key"); then
    echo "✗ Could not locate data directory for '$key'" >&2
    return 1
  fi
  echo "Data directory: $data_dir"

  case "$key" in
    thermalization)
      run_thermalization_analysis "$data_dir"
      ;;
    transport)
      run_transport_analysis "$data_dir"
      ;;
    spacecharge)
      run_spacecharge_analysis "$data_dir"
      ;;
    reactions)
      run_reactions_analysis "$data_dir"
      ;;
    *)
      echo "Internal error: unknown analysis '$key'" >&2
      return 1
      ;;
  esac

  echo "✓ Completed $label"
}

TOTAL=0
PASSED=0
FAILED=0
FAILED_KEYS=()

START_TIME=$(date)

echo "=============================================="
echo "Physics Analysis Suite"
echo "=============================================="
echo "Repo root      : $REPO_ROOT"
echo "Python         : $PYTHON_BIN"
echo "Results root   : ${RESULTS_ROOT_OVERRIDE:-auto}" 
echo "Overrides      : thermalization=${THERM_DIR_OVERRIDE:-auto}, transport=${TRANSPORT_DIR_OVERRIDE:-auto}, spacecharge=${SPACECHARGE_DIR_OVERRIDE:-auto}, reactions=${REACTIONS_DIR_OVERRIDE:-auto}"
echo "Analyses       : ${SELECTED_ANALYSES[*]}"
echo "Started        : $START_TIME"
echo "=============================================="
echo

for key in "${SELECTED_ANALYSES[@]}"; do
  TOTAL=$((TOTAL + 1))
  if run_analysis "$key"; then
    PASSED=$((PASSED + 1))
  else
    FAILED=$((FAILED + 1))
    FAILED_KEYS+=("$key")
  fi
  echo
 done

END_TIME=$(date)

echo "=============================================="
echo "Analysis summary"
echo "=============================================="
echo "Requested analyses : ${SELECTED_ANALYSES[*]}"
echo "Succeeded          : $PASSED"
echo "Failed             : $FAILED"
echo "Finished           : $END_TIME"
if [[ $FAILED -gt 0 ]]; then
  echo "Failures           : ${FAILED_KEYS[*]}" >&2
fi
echo "=============================================="

if [[ $FAILED -gt 0 ]]; then
  exit 1
fi

exit 0
