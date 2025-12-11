#!/bin/bash
# Unified ICARION performance benchmark runner (CPU + GPU)

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PROJECT_ROOT="$(cd "$SCRIPT_DIR/../../.." && pwd)"
BUILD_DIR="$PROJECT_ROOT/build"
CPU_CONFIG_DIR="$PROJECT_ROOT/validation/configs/performance"
GPU_CONFIG_DIR="$PROJECT_ROOT/validation/configs/performance/gpu"
CPU_RESULT_DIR="$PROJECT_ROOT/validation/results/v1.0_test/performance/logs"
GPU_RESULT_DIR="$PROJECT_ROOT/validation/results/v1.0_test/performance/gpu_logs"
ICARION_BIN="$BUILD_DIR/src/icarion_main"
CMAKE_CACHE="$BUILD_DIR/CMakeCache.txt"

declare -A CPU_CATEGORY_GLOBS=(
  [baseline]="scaling_baseline_*.json"
  [collision]="collision_overhead_*.json"
  [spacecharge]="spacecharge_overhead_*.json"
  [longrun]="scaling_longrun_*.json"
  [mixed]="scaling_mixedphysics_*.json"
)
declare -a CPU_DEFAULT_CATEGORIES=(baseline collision spacecharge longrun mixed)

declare -A GPU_CATEGORY_GLOBS=(
  [cpu_scaling]="RK4_cpu_*.json"
  [gpu_scaling]="RK4_gpu_*.json"
  [integrators]="integrator_*.json"
  [threshold]="threshold_*.json"
  [long]="long_*.json"
)
declare -a GPU_DEFAULT_CATEGORIES=(cpu_scaling gpu_scaling integrators threshold long)

print_usage() {
  cat <<EOF
Usage: $(basename "$0") [options] [category ...]

Options:
  --cpu-only        Run only CPU benchmark categories
  --gpu-only        Run only GPU benchmark categories (requires USE_GPU_ACCEL=ON)
  -h, --help        Show this help and exit

Categories:
  CPU: ${CPU_DEFAULT_CATEGORIES[*]}
  GPU: ${GPU_DEFAULT_CATEGORIES[*]}
  Special aliases:
    all       Run all CPU categories and GPU categories (GPU section skipped if disabled)
    cpu_all   Run all CPU categories
    gpu_all   Run all GPU categories (requires GPU build)
EOF
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

MODE="all" # all|cpu|gpu
USER_CATEGORIES=()
while [[ $# -gt 0 ]]; do
  case "$1" in
    --cpu-only)
      MODE="cpu"
      shift
      ;;
    --gpu-only)
      MODE="gpu"
      shift
      ;;
    -h|--help)
      print_usage
      exit 0
      ;;
    *)
      USER_CATEGORIES+=("$1")
      shift
      ;;
  esac
done

declare -a CPU_SELECTED=()
declare -a GPU_SELECTED=()

if [ ${#USER_CATEGORIES[@]} -eq 0 ]; then
  CPU_SELECTED=("${CPU_DEFAULT_CATEGORIES[@]}")
  GPU_SELECTED=("${GPU_DEFAULT_CATEGORIES[@]}")
else
  for category in "${USER_CATEGORIES[@]}"; do
    case "$category" in
      all)
        CPU_SELECTED+=("${CPU_DEFAULT_CATEGORIES[@]}")
        GPU_SELECTED+=("${GPU_DEFAULT_CATEGORIES[@]}")
        ;;
      cpu_all)
        CPU_SELECTED+=("${CPU_DEFAULT_CATEGORIES[@]}")
        ;;
      gpu_all)
        GPU_SELECTED+=("${GPU_DEFAULT_CATEGORIES[@]}")
        ;;
      *)
        if [[ -n "${CPU_CATEGORY_GLOBS[$category]:-}" ]]; then
          CPU_SELECTED+=("$category")
        elif [[ -n "${GPU_CATEGORY_GLOBS[$category]:-}" ]]; then
          GPU_SELECTED+=("$category")
        else
          echo "Unknown category: $category" >&2
          print_usage
          exit 1
        fi
        ;;
    esac
  done
fi

if [[ "$MODE" == "cpu" ]]; then
  GPU_SELECTED=()
elif [[ "$MODE" == "gpu" ]]; then
  CPU_SELECTED=()
fi

dedupe_array CPU_SELECTED
 dedupe_array GPU_SELECTED

if [ ! -x "$ICARION_BIN" ]; then
  echo "Error: ICARION binary not found at $ICARION_BIN" >&2
  exit 1
fi

GPU_ENABLED=0
if [ -f "$CMAKE_CACHE" ] && grep -q '^USE_GPU_ACCEL:BOOL=ON' "$CMAKE_CACHE"; then
  GPU_ENABLED=1
fi

if [ ${#GPU_SELECTED[@]} -gt 0 ] && [ $GPU_ENABLED -ne 1 ]; then
  if [[ "$MODE" == "gpu" ]]; then
    echo "Error: GPU benchmarks requested but this build was configured without USE_GPU_ACCEL=ON" >&2
    exit 1
  fi
  echo "Skipping GPU benchmarks because USE_GPU_ACCEL=OFF in $CMAKE_CACHE"
  GPU_SELECTED=()
fi

if [ ${#CPU_SELECTED[@]} -gt 0 ] && [ ! -d "$CPU_CONFIG_DIR" ]; then
  echo "Error: CPU config directory not found at $CPU_CONFIG_DIR" >&2
  exit 1
fi

if [ ${#GPU_SELECTED[@]} -gt 0 ] && [ ! -d "$GPU_CONFIG_DIR" ]; then
  echo "Error: GPU config directory not found at $GPU_CONFIG_DIR" >&2
  exit 1
fi

mkdir -p "$CPU_RESULT_DIR"
mkdir -p "$GPU_RESULT_DIR"

echo ""
echo "============================================================"
echo "ICARION Performance Benchmark Suite"
echo "============================================================"
echo "Binary:  $ICARION_BIN"
echo "CPU configs: $CPU_CONFIG_DIR"
echo "GPU configs: $GPU_CONFIG_DIR"
echo "============================================================"
echo ""

TOTAL_RUNS=0
CPU_RUNS=0
GPU_RUNS=0

run_config() {
  local config_path="$1"
  local result_dir="$2"
  local name="$(basename "$config_path" .json)"
  local log_file="$result_dir/$name.log"
  local time_file="$result_dir/$name.time"

  echo "→ Running $name"
  if /usr/bin/time -f "%e" -o "$time_file" "$ICARION_BIN" "$config_path" > "$log_file" 2>&1; then
    local elapsed
    elapsed="$(cat "$time_file")"
    echo "  ✓ Completed (${elapsed}s)"
  else
    local status=$?
    echo "  ✗ Failed (exit $status)"
    echo "    Log: $log_file"
  fi
}

run_category() {
  local config_dir="$1"
  local result_dir="$2"
  local category="$3"
  local pattern="$4"
  local kind="$5"

  echo ""
  echo "========================================"
  echo "Category: $category"
  echo "========================================"
  echo ""

  local matches=0
  shopt -s nullglob
  for config in "$config_dir"/$pattern; do
    matches=$((matches + 1))
    run_config "$config" "$result_dir"
  done
  shopt -u nullglob

  if [ $matches -eq 0 ]; then
    echo "(no configs matched pattern '$pattern')"
  else
    echo ""
    echo "Completed $matches runs for '$category'"
  fi

  TOTAL_RUNS=$((TOTAL_RUNS + matches))
  if [[ "$kind" == "cpu" ]]; then
    CPU_RUNS=$((CPU_RUNS + matches))
  else
    GPU_RUNS=$((GPU_RUNS + matches))
  fi
}

for category in "${CPU_SELECTED[@]:-}"; do
  pattern=${CPU_CATEGORY_GLOBS[$category]:-}
  if [ -n "$pattern" ]; then
    run_category "$CPU_CONFIG_DIR" "$CPU_RESULT_DIR" "$category" "$pattern" "cpu"
  fi
 done

if [ ${#GPU_SELECTED[@]} -gt 0 ]; then
  echo ""
  echo "GPU benchmarks enabled (USE_GPU_ACCEL=ON)"
fi

for category in "${GPU_SELECTED[@]:-}"; do
  pattern=${GPU_CATEGORY_GLOBS[$category]:-}
  if [ -n "$pattern" ]; then
    run_category "$GPU_CONFIG_DIR" "$GPU_RESULT_DIR" "$category" "$pattern" "gpu"
  fi
 done

if [ $TOTAL_RUNS -eq 0 ]; then
  echo "No benchmarks were executed."
  exit 0
fi

echo ""
echo "============================================================"
echo "Summary"
echo "============================================================"
echo "CPU runs: $CPU_RUNS (logs: $CPU_RESULT_DIR)"
if [ $GPU_RUNS -gt 0 ]; then
  echo "GPU runs: $GPU_RUNS (logs: $GPU_RESULT_DIR)"
else
  echo "GPU runs: 0"
fi
echo "Total runs: $TOTAL_RUNS"
echo "============================================================"
echo ""
