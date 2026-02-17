#!/bin/bash
#
# run_quad_tests.sh
#
# Batch runner for quadrupole stability map validation tests
# Runs the first stability-region (a,q) grid with CaffeineH+ ions
#
# Usage: 
#   ./run_quad_tests.sh [parallel_jobs]
#
# Example:
#   ./run_quad_tests.sh 2    # Run 2 simulations in parallel
#   nohup ./run_quad_tests.sh 4 > quad_stability.log 2>&1 &  # Background with 4 parallel jobs

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
VALIDATION_DIR="$(cd "$SCRIPT_DIR/../../.." && pwd)"
REPO_ROOT="$(cd "$VALIDATION_DIR/.." && pwd)"

ICARION_BIN_DEFAULT="$REPO_ROOT/build/src/icarion_main"
CONFIG_DIR_DEFAULT="$VALIDATION_DIR/configs/instruments/quadrupole"
OUTPUT_DIR_DEFAULT="$VALIDATION_DIR/results/v1.0.0_test/instruments/quadrupole_first_region"

RUN_DIR="${ICARION_VALIDATION_RUN_DIR:-}"
if [[ -n "$RUN_DIR" ]]; then
    OUTPUT_DIR_DEFAULT="$RUN_DIR/results/instruments/quadrupole_first_region"
fi

ICARION_BIN="${ICARION_BIN_OVERRIDE:-$ICARION_BIN_DEFAULT}"
CONFIG_DIR="${INSTRUMENT_CONFIG_DIR_OVERRIDE:-$CONFIG_DIR_DEFAULT}"
OUTPUT_DIR="${INSTRUMENT_OUTPUT_ROOT_OVERRIDE:-$OUTPUT_DIR_DEFAULT}"

JOBS_DEFAULT=${INSTRUMENT_JOBS:-2}
THREADS_DEFAULT=${INSTRUMENT_THREADS:-1}

MAX_PARALLEL="$JOBS_DEFAULT"
THREADS="$THREADS_DEFAULT"

print_usage() {
    cat <<'EOF'
Usage: run_quad_tests.sh [options]

Options:
    -j, --jobs N         Parallel jobs (default: from INSTRUMENT_JOBS or 2)
    -t, --threads N      Threads passed to icarion_main (default: from INSTRUMENT_THREADS or 1)
    -b, --binary PATH    Path to icarion_main (default: build/src/icarion_main)
    -c, --config-dir DIR Config directory (default: validation/configs/instruments/quadrupole)
    -o, --output-dir DIR Output directory root (default: validation/results/... or validation/runs/<id>/results/...)
    -h, --help           Show this help

Back-compat:
    run_quad_tests.sh 4  # treated as --jobs 4
EOF
}

while [[ $# -gt 0 ]]; do
    case "$1" in
        -j|--jobs)
            MAX_PARALLEL="$2"; shift 2 ;;
        -t|--threads)
            THREADS="$2"; shift 2 ;;
        -b|--binary)
            ICARION_BIN="$2"; shift 2 ;;
        -c|--config-dir)
            CONFIG_DIR="$2"; shift 2 ;;
        -o|--output-dir)
            OUTPUT_DIR="$2"; shift 2 ;;
        -h|--help)
            print_usage; exit 0 ;;
        --)
            shift; break ;;
        -* )
            echo "Unknown option: $1" >&2
            print_usage >&2
            exit 1
            ;;
        *)
            # Back-compat: first bare arg is jobs
            if [[ "$1" =~ ^[0-9]+$ ]]; then
                MAX_PARALLEL="$1"; shift
            else
                echo "Unknown argument: $1" >&2
                print_usage >&2
                exit 1
            fi
            ;;
    esac
done

if [[ "$ICARION_BIN" != /* ]]; then
    ICARION_BIN="$REPO_ROOT/$ICARION_BIN"
fi
if [[ "$CONFIG_DIR" != /* ]]; then
    CONFIG_DIR="$REPO_ROOT/$CONFIG_DIR"
fi
if [[ "$OUTPUT_DIR" != /* ]]; then
    OUTPUT_DIR="$REPO_ROOT/$OUTPUT_DIR"
fi

ICARION_BIN="$(cd "$(dirname "$ICARION_BIN")" && pwd)/$(basename "$ICARION_BIN")"
CONFIG_DIR="$(cd "$CONFIG_DIR" && pwd)"
mkdir -p "$OUTPUT_DIR"
OUTPUT_DIR="$(cd "$OUTPUT_DIR" && pwd)"

LOG_DIR="$OUTPUT_DIR/logs"
CFG_SNAPSHOT_DIR="$OUTPUT_DIR/run_configs"

# Create output directories
mkdir -p "${OUTPUT_DIR}"
mkdir -p "${LOG_DIR}"
mkdir -p "${CFG_SNAPSHOT_DIR}"

# Check if executable exists
if [[ ! -f "${ICARION_BIN}" ]]; then
    echo "❌ ERROR: ICARION executable not found: ${ICARION_BIN}"
    exit 1
fi

# Check if config directory exists
if [[ ! -d "${CONFIG_DIR}" ]]; then
    echo "❌ ERROR: Config directory not found: ${CONFIG_DIR}"
    exit 1
fi

# Count config files
CONFIG_FILES=(${CONFIG_DIR}/quad_stability_*.json)
TOTAL_CONFIGS=${#CONFIG_FILES[@]}

if [[ ${TOTAL_CONFIGS} -eq 0 ]]; then
    echo "❌ ERROR: No config files found in ${CONFIG_DIR}"
    exit 1
fi

echo "=============================================="
echo "Quadrupole Stability Map Validation"
echo "=============================================="
echo "Executable:  ${ICARION_BIN}"
echo "Configs:     ${CONFIG_DIR}"
echo "Output:      ${OUTPUT_DIR}"
echo "Total tests: ${TOTAL_CONFIGS}"
echo "Parallel:    ${MAX_PARALLEL} jobs"
echo "Threads:     ${THREADS}"
echo "=============================================="
echo ""

# Function to run a single test
run_test() {
    local config_file=$1
    local config_name=$(basename "${config_file}" .json)
    local log_file="${LOG_DIR}/${config_name}.log"

    # Snapshot config with overridden output.folder so artifacts land in OUTPUT_DIR.
    local run_cfg="${CFG_SNAPSHOT_DIR}/${config_name}.run_config.json"
    python3 - "${config_file}" "${run_cfg}" "${OUTPUT_DIR}" <<'PY'
import json
import sys

src, dst, out_dir = sys.argv[1], sys.argv[2], sys.argv[3]
with open(src, 'r', encoding='utf-8') as handle:
    cfg = json.load(handle)
cfg.setdefault('output', {})['folder'] = out_dir
with open(dst, 'w', encoding='utf-8') as handle:
    json.dump(cfg, handle, indent=2)
PY
    
    echo "[$(date '+%H:%M:%S')] Running ${config_name}..."
    
    # Run simulation and capture output
    if (cd "$REPO_ROOT" && "${ICARION_BIN}" "${run_cfg}" --threads "${THREADS}" > "${log_file}" 2>&1); then
        echo "[$(date '+%H:%M:%S')] ✅ ${config_name} completed"
        return 0
    else
        echo "[$(date '+%H:%M:%S')] ❌ ${config_name} FAILED (see ${log_file})"
        return 1
    fi
}

# Export function for parallel execution
export -f run_test
export ICARION_BIN
export REPO_ROOT
export LOG_DIR
export OUTPUT_DIR
export CFG_SNAPSHOT_DIR
export THREADS

# Start time
START_TIME=$(date +%s)
echo "Started at: $(date)"
echo ""

# Run tests in parallel using GNU parallel if available, otherwise xargs
if command -v parallel &> /dev/null; then
    # Use GNU parallel
    printf '%s\n' "${CONFIG_FILES[@]}" | parallel -j ${MAX_PARALLEL} --bar run_test {}
else
    # Fallback to xargs
    printf '%s\n' "${CONFIG_FILES[@]}" | xargs -P ${MAX_PARALLEL} -I {} bash -c 'run_test "$@"' _ {}
fi

# End time
END_TIME=$(date +%s)
ELAPSED=$((END_TIME - START_TIME))
ELAPSED_MIN=$((ELAPSED / 60))
ELAPSED_SEC=$((ELAPSED % 60))

echo ""
echo "=============================================="
echo "Batch run completed!"
echo "Total time: ${ELAPSED_MIN}m ${ELAPSED_SEC}s"
echo "=============================================="

# Count successes and failures
SUCCESS_COUNT=$(grep -El "Simulation completed successfully|=== Simulation Complete ===" "${LOG_DIR}"/*.log 2>/dev/null | wc -l)
FAIL_COUNT=$((TOTAL_CONFIGS - SUCCESS_COUNT))

echo ""
echo "Results summary:"
echo "  ✅ Successful: ${SUCCESS_COUNT}/${TOTAL_CONFIGS}"
echo "  ❌ Failed:     ${FAIL_COUNT}/${TOTAL_CONFIGS}"
echo ""

# List HDF5 output files
HDF5_COUNT=$(find "${OUTPUT_DIR}" -name "*.h5" 2>/dev/null | wc -l)
echo "HDF5 outputs:  ${HDF5_COUNT} files"
echo "Output dir:    ${OUTPUT_DIR}"
echo "Logs dir:      ${LOG_DIR}"
echo ""

if [[ ${FAIL_COUNT} -gt 0 ]]; then
    echo "⚠️  Some tests failed. Check logs in ${LOG_DIR}"
    echo ""
    echo "Failed tests:"
    for log in "${LOG_DIR}"/*.log; do
        if ! grep -Eq "Simulation completed successfully|=== Simulation Complete ===" "${log}" 2>/dev/null; then
            echo "  - $(basename "${log}" .log)"
        fi
    done
else
    echo "✅ All tests completed successfully!"
fi

echo ""
echo "Next steps:"
echo "  1. Analyze results: python3 validation/scripts/instrumentation/quadrupole/analyze_quad_stability.py"
echo "  2. Plot stability map: python3 validation/scripts/instrumentation/quadrupole/plot_quad_stability.py"
echo ""
