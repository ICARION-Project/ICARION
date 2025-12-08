#!/bin/bash
# Run ICARION GPU/Hybrid performance benchmarks with timing capture

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PROJECT_ROOT="$(cd "$SCRIPT_DIR/../../.." && pwd)"
CONFIG_DIR="$PROJECT_ROOT/validation/configs/performance/gpu"
RESULT_DIR="$PROJECT_ROOT/validation/results/v1.0_test/performance/gpu_logs"
ICARION_BIN="$PROJECT_ROOT/build/src/icarion_main"

mkdir -p "$RESULT_DIR"

if [ ! -x "$ICARION_BIN" ]; then
  echo "Error: ICARION binary not found at $ICARION_BIN" >&2
  exit 1
fi

declare -A CATEGORY_GLOBS=(
  [cpu_scaling]="RK4_cpu_*.json"
  [gpu_scaling]="RK4_gpu_*.json"
  [integrators]="integrator_*.json"
  [threshold]="threshold_*.json"
  [long]="long_*.json"
)

print_usage() {
  cat <<EOF
Usage: $(basename "$0") [category ...]

Categories:
  cpu_scaling   CPU baselines that mirror GPU configs
  gpu_scaling   GPU-enabled counterparts of the CPU runs
  integrators   GPU integrator comparisons (RK4/RK45/Boris)
  threshold     Threshold validation near GPU handoff
  long          Long-duration GPU efficiency tests
  all           Run every category above (default)
EOF
}

if [ ${#} -eq 0 ]; then
  CATEGORIES=(cpu_scaling gpu_scaling integrators threshold long)
else
  if [[ "$1" == "-h" || "$1" == "--help" ]]; then
    print_usage
    exit 0
  fi
  if [[ "$1" == "all" ]]; then
    CATEGORIES=(cpu_scaling gpu_scaling integrators threshold long)
  else
    CATEGORIES=("$@")
  fi
fi

run_config() {
  local config_path="$1"
  local name="$(basename "$config_path" .json)"
  local log_file="$RESULT_DIR/$name.log"
  local time_file="$RESULT_DIR/$name.time"

  echo "→ Running $name"
  if /usr/bin/time -f "%e" -o "$time_file" "$ICARION_BIN" "$config_path" > "$log_file" 2>&1; then
    local elapsed="$(cat "$time_file")"
    echo "  ✓ Completed ($elapsed s)"
  else
    local status=$?
    echo "  ✗ Failed (exit $status)"
    echo "    Log: $log_file"
  fi
}

TOTAL_RUNS=0
shopt -s nullglob
for category in "${CATEGORIES[@]}"; do
  if [[ "$category" == "all" ]]; then
    category="cpu_scaling"
  fi
  pattern=${CATEGORY_GLOBS[$category]:-}
  if [ -z "$pattern" ]; then
    echo "Unknown category: $category" >&2
    continue
  fi
  echo ""
  echo "========================================"
  echo "Category: $category"
  echo "========================================"
  matches=0
  for config in "$CONFIG_DIR"/$pattern; do
    matches=$((matches + 1))
    run_config "$config"
  done
  if [ $matches -eq 0 ]; then
    echo "(no configs matched pattern '$pattern')"
  fi
  TOTAL_RUNS=$((TOTAL_RUNS + matches))
  echo ""
  echo "Completed $matches runs for '$category'"

done
shopt -u nullglob

if [ $TOTAL_RUNS -eq 0 ]; then
  echo "No GPU benchmarks were executed."
else
  echo ""
  echo "Logs: $RESULT_DIR"
fi
