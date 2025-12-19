#!/bin/bash
# Orchestrator for the core physics validation drivers.

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
VALIDATION_DIR="$(cd "$SCRIPT_DIR/.." && pwd)"
REPO_ROOT="$(cd "$VALIDATION_DIR/.." && pwd)"
THERMALIZATION_RUNNER="$VALIDATION_DIR/scripts/thermalization/run_thermalization_tests.sh"
PHYSICS_SCRIPT_DIR="$VALIDATION_DIR/scripts/physics"
THERMALIZATION_RUNNER_DIR="$(dirname "$THERMALIZATION_RUNNER")"

PYTHON_BIN=${PYTHON_BIN:-python3}
# Default to full sweep to cover all thermalization configs
THERMALIZATION_MODE="full"
ICARION_BIN=""
DRY_RUN=false
CONFIG_JOBS=${CONFIG_JOBS:-8}
SELECTED_TARGETS=()

AVAILABLE_TARGETS=(
  thermalization
  gas_flow
  combined_drift
  gas_mixture_mobility
  mixture_thermalization
  spacecharge
  reactions
)

declare -A TARGET_LABELS=(
  [thermalization]="Session 2 thermalization sweep"
  [gas_flow]="Gas flow transport (E=0)"
  [combined_drift]="Combined drift (E + gas flow)"
  [gas_mixture_mobility]="Blanc's law mixture mobility"
  [mixture_thermalization]="Gas mixture thermalization"
  [spacecharge]="Space charge expansion / IMS space charge"
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
  --config-jobs N            Max configs to run in parallel (where applicable, default: $CONFIG_JOBS)
  --list                     Show available targets and exit
  -h, --help                 Show this help

Targets (default: all):
  thermalization             ${TARGET_LABELS[thermalization]}
  gas_flow                   ${TARGET_LABELS[gas_flow]}
  combined_drift             ${TARGET_LABELS[combined_drift]}
  gas_mixture_mobility       ${TARGET_LABELS[gas_mixture_mobility]}
  mixture_thermalization     ${TARGET_LABELS[mixture_thermalization]}
  spacecharge                ${TARGET_LABELS[spacecharge]}
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

list_targets() {
  cat <<EOF
Available targets:
  thermalization         ${TARGET_LABELS[thermalization]}
  gas_flow               ${TARGET_LABELS[gas_flow]}
  combined_drift         ${TARGET_LABELS[combined_drift]}
  gas_mixture_mobility   ${TARGET_LABELS[gas_mixture_mobility]}
  mixture_thermalization ${TARGET_LABELS[mixture_thermalization]}
  spacecharge            ${TARGET_LABELS[spacecharge]}
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
    --config-jobs)
      [[ $# -lt 2 ]] && print_usage && exit 1
      CONFIG_JOBS="$2"
      shift 2
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
done

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

NEED_ICARION=false
for tgt in "${SELECTED_TARGETS[@]}"; do
  if [[ "$tgt" == "spacecharge" || "$tgt" == "reactions" ]]; then
    NEED_ICARION=true
    break
  fi
done

if $NEED_ICARION; then
  if [[ -z "$ICARION_BIN" ]]; then
    ICARION_BIN="$REPO_ROOT/build/src/icarion_main"
  fi
  if [[ ! -x "$ICARION_BIN" ]]; then
    echo "Error: ICARION binary not found or not executable: $ICARION_BIN" >&2
    echo "Hint: set --icarion-bin or build the project first." >&2
    exit 1
  fi
fi

run_in_repo() {
  (cd "$REPO_ROOT" && "$@")
}

run_spacecharge_configs() {
  local bin="${ICARION_BIN:-$REPO_ROOT/build/src/icarion_main}"
  local cfg_dir="$REPO_ROOT/validation/configs/physics/spacecharge"
  local -a cfgs=("$cfg_dir"/*.json)

  if [[ ! -x "$bin" ]]; then
    echo "[spacecharge] icarion_main not found or not executable at '$bin'" >&2
    return 1
  fi
  if [[ ${#cfgs[@]} -eq 0 ]]; then
    echo "[spacecharge] No configs found in $cfg_dir" >&2
    return 1
  fi

  local -a pids=()
  local failures=0

  for cfg in "${cfgs[@]}"; do
    [[ ! -f "$cfg" ]] && continue
    echo "[spacecharge] running $cfg"
    ("$bin" "$cfg" --threads 1 >/dev/null 2>&1) &
    pids+=("$!")

    while [[ ${#pids[@]} -ge $CONFIG_JOBS && $CONFIG_JOBS -gt 0 ]]; do
      if ! wait -n; then
        failures=$((failures + 1))
      fi
      # prune finished
      local new_pids=()
      for pid in "${pids[@]}"; do
        if kill -0 "$pid" 2>/dev/null; then
          new_pids+=("$pid")
        fi
      done
      pids=("${new_pids[@]}")
    done
  done

  for pid in "${pids[@]:-}"; do
    wait "$pid" || failures=$((failures + 1))
  done

  if [[ $failures -gt 0 ]]; then
    echo "[spacecharge] $failures config(s) failed" >&2
    return 1
  fi
  return 0
}

build_command() {
  local target="$1"
  case "$target" in
    thermalization)
      printf '%q ' bash -c \
        "cd \"$THERMALIZATION_RUNNER_DIR\" && THERM_JOBS=\"$CONFIG_JOBS\" ./$(basename "$THERMALIZATION_RUNNER") \"$THERMALIZATION_MODE\""
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
    spacecharge)
      printf '%s' "__SPACECHARGE_INTERNAL__"
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

  if [[ "$cmd_str" == "__SPACECHARGE_INTERNAL__" ]]; then
    echo "============================================================"
    echo "$label"
    echo "============================================================"
    echo "Command: run_spacecharge_configs (jobs=$CONFIG_JOBS)"

    if $DRY_RUN; then
      echo "[dry-run] Skipping execution"
      return 0
    fi

    if run_spacecharge_configs; then
      echo "[$target] ✅ completed"
      return 0
    else
      echo "[$target] ❌ failed" >&2
      return 1
    fi
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
