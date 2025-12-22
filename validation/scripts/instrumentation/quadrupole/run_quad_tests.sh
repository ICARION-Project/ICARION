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

set -e

# Configuration
ICARION_BIN="/home/chsch95/ICARION/build/src/icarion_main"
CONFIG_DIR="/home/chsch95/ICARION/validation/configs/instruments/quadrupole"
OUTPUT_DIR="/home/chsch95/ICARION/validation/results/v1.0_test/instruments/quadrupole_first_region"
LOG_DIR="${OUTPUT_DIR}/logs"
MAX_PARALLEL=${1:-2}  # Default to 2 parallel jobs

# Create output directories
mkdir -p "${OUTPUT_DIR}"
mkdir -p "${LOG_DIR}"

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
echo "=============================================="
echo ""

# Function to run a single test
run_test() {
    local config_file=$1
    local config_name=$(basename "${config_file}" .json)
    local log_file="${LOG_DIR}/${config_name}.log"
    
    echo "[$(date '+%H:%M:%S')] Running ${config_name}..."
    
    # Run simulation and capture output
    if "${ICARION_BIN}" "${config_file}" > "${log_file}" 2>&1; then
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
export LOG_DIR

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
