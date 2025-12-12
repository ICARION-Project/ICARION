#!/bin/bash
# Orchestrator for the core physics validation drivers.

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
VALIDATION_DIR="$(cd "$SCRIPT_DIR/.." && pwd)"
REPO_ROOT="$(cd "$VALIDATION_DIR/.." && pwd)"
THERMALIZATION_RUNNER="$VALIDATION_DIR/scripts/thermalization/run_thermalization_tests.sh"
PHYSICS_SCRIPT_DIR="$VALIDATION_DIR/scripts/physics"

PYTHON_BIN=${PYTHON_BIN:-python3}
THERMALIZATION_MODE="quick"
ICARION_BIN=""
DRY_RUN=false
SELECTED_TARGETS=()

AVAILABLE_TARGETS=(
  thermalization
  gas_flow
  combined_drift
  gas_mixture_mobility
  mixture_thermalization
  reactions
)

declare -A TARGET_LABELS=(
  [thermalization]="Session 2 thermalization sweep"
  [gas_flow]="Gas flow transport (E=0)"
  [combined_drift]="Combined drift (E + gas flow)"
  [gas_mixture_mobility]="Blanc's law mixture mobility"
  [mixture_thermalization]="Gas mixture thermalization"
  [reactions]="Reaction kinetics scenarios"
)

print_usage() {
  cat <<'EOF'
Usage: ./run_physics_suite.sh [options] [target ...]

Options:
  --python PATH              Python interpreter to use (default: $PYTHON_BIN)
  --thermalization-mode MODE Mode passed to run_thermalization_tests.sh (quick|subset|full)
  --icarion-bin PATH         Forward PATH to validate_reaction_kinetics.py
  --dry-run                  Show the commands without executing them
  --list                     Show available targets and exit
  -h, --help                 Show this help

Targets (default: all):
  thermalization             ${TARGET_LABELS[thermalization]}
  gas_flow                   ${TARGET_LABELS[gas_flow]}
  combined_drift             ${TARGET_LABELS[combined_drift]}
  gas_mixture_mobility       ${TARGET_LABELS[gas_mixture_mobility]}
  mixture_thermalization     ${TARGET_LABELS[mixture_thermalization]}
  reactions                  ${TARGET_LABELS[reactions]}
EOF
}

normalize_target() {
  local token="${1,,}"
  case "$token" in
    all|suite)
      printf '%s\n' "all"
      ;;
    thermalization|therm|session2)
      printf '%s\n' "thermalization"
      ;;
    gas-flow|gas_flow|gasflow)
      printf '%s\n' "gas_flow"
      ;;
    combined-drift|combined_drift|drift)
      printf '%s\n' "combined_drift"
      ;;
    gas-mixture|gas_mixture|mixture-mobility|blanc|mobility)
      printf '%s\n' "gas_mixture_mobility"
      ;;
    mixture-thermalization|mixture_thermalization|mix-therm)
      printf '%s\n' "mixture_thermalization"
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

list_targets() {
  cat <<EOF
Available targets:
  thermalization         ${TARGET_LABELS[thermalization]}
  gas_flow               ${TARGET_LABELS[gas_flow]}
  combined_drift         ${TARGET_LABELS[combined_drift]}
  gas_mixture_mobility   ${TARGET_LABELS[gas_mixture_mobility]}
  mixture_thermalization ${TARGET_LABELS[mixture_thermalization]}
  reactions              ${TARGET_LABELS[reactions]}
EOF
}

while [[ $# -gt 0 ]]; do
  case "$1" in
    --python)
      [[ $# -lt 2 ]] && print_usage && exit 1
      PYTHON_BIN="$2"
      shift 2
      ;;
    --thermalization-mode)
      [[ $# -lt 2 ]] && print_usage && exit 1
      THERMALIZATION_MODE="$2"
      shift 2
      ;;
    --icarion-bin)
      [[ $# -lt 2 ]] && print_usage && exit 1
      ICARION_BIN="$2"
      shift 2
      ;;
    --dry-run)
      DRY_RUN=true
      shift
      ;;
    --list)
      list_targets
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
      if ! normalized=$(normalize_target "$1"); then
        echo "Unknown target: $1" >&2
        exit 1
      fi
      if [[ "$normalized" == "all" ]]; then
        SELECTED_TARGETS=("${AVAILABLE_TARGETS[@]}")
      else
        SELECTED_TARGETS+=("$normalized")
      fi
      shift
      ;;
  esac
end

if [[ ${#SELECTED_TARGETS[@]} -eq 0 ]]; then
  SELECTED_TARGETS=("${AVAILABLE_TARGETS[@]}")
fi

dedupe_array SELECTED_TARGETS

if ! command -v "$PYTHON_BIN" >/dev/null 2>&1; then
  echo "Error: Python interpreter '$PYTHON_BIN' not found" >&2
  exit 1
fi

if [[ ! -x "$THERMALIZATION_RUNNER" ]]; then
  echo "Error: Thermalization runner not found: $THERMALIZATION_RUNNER" >&2
  exit 1
fi

run_in_repo() {
  (cd "$REPO_ROOT" && "$@")
}

build_command() {
  local target="$1"
  case "$target" in
    thermalization)
      printf '%q ' "$THERMALIZATION_RUNNER" "$THERMALIZATION_MODE"
      ;;
    gas_flow)
      printf '%q ' "$PYTHON_BIN" "$PHYSICS_SCRIPT_DIR/validate_gas_flow_transport.py"
      ;;
    combined_drift)
      printf '%q ' "$PYTHON_BIN" "$PHYSICS_SCRIPT_DIR/validate_combined_drift.py"
      ;;
    gas_mixture_mobility)
      printf '%q ' "$PYTHON_BIN" "$PHYSICS_SCRIPT_DIR/validate_gas_mixture_mobility.py"
      ;;
    mixture_thermalization)
      printf '%q ' "$PYTHON_BIN" "$PHYSICS_SCRIPT_DIR/validate_mixture_thermalization.py"
      ;;
    reactions)
      printf '%q ' "$PYTHON_BIN" "$PHYSICS_SCRIPT_DIR/validate_reaction_kinetics.py"
      if [[ -n "$ICARION_BIN" ]]; then
        printf '%q ' --icarion-bin "$ICARION_BIN"
      fi
      ;;
    *)
      return 1
      ;;
  esac
}

run_target() {
  local target="$1"
  local label="${TARGET_LABELS[$target]}"
  local cmd_str
  if ! cmd_str=$(build_command "$target"); then
    echo "Internal error: unsupported target '$target'" >&2
    return 1
  fi

  echo "============================================================"
  echo "$label"
  echo "============================================================"
  echo "Command: $cmd_str"

  if $DRY_RUN; then
    echo "[dry-run] Skipping execution"
    return 0
  fi

  if run_in_repo bash -c "$cmd_str"; then
    echo "[$target] ✅ completed"
    return 0
  else
    echo "[$target] ❌ failed" >&2
    return 1
  fi
}

TOTAL=0
PASSED=0
FAILED=0
FAILED_TARGETS=()

START_TIME=$(date)

echo "=============================================="
echo "Physics Validation Suite"
echo "=============================================="
echo "Repo root      : $REPO_ROOT"
echo "Python         : $PYTHON_BIN"
echo "Therm mode     : $THERMALIZATION_MODE"
echo "ICARION bin    : ${ICARION_BIN:-default from scripts}"
echo "Dry run        : $DRY_RUN"
echo "Targets        : ${SELECTED_TARGETS[*]}"
echo "Started        : $START_TIME"
echo "=============================================="
echo

for target in "${SELECTED_TARGETS[@]}"; do
  TOTAL=$((TOTAL + 1))
  if run_target "$target"; then
    PASSED=$((PASSED + 1))
  else
    FAILED=$((FAILED + 1))
    FAILED_TARGETS+=("$target")
  fi
  echo
 done

END_TIME=$(date)

echo "=============================================="
echo "Suite summary"
echo "=============================================="
echo "Requested targets : ${SELECTED_TARGETS[*]}"
echo "Succeeded         : $PASSED"
echo "Failed            : $FAILED"
echo "Finished          : $END_TIME"
if [[ $FAILED -gt 0 ]]; then
  echo "Failures          : ${FAILED_TARGETS[*]}" >&2
fi
echo "=============================================="

if [[ $FAILED -gt 0 ]]; then
  exit 1
fi

exit 0
