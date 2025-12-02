#!/bin/bash
################################################################################
# GPU Performance Benchmark Suite Runner
#
# Runs all GPU performance benchmarks and validates GPU acceleration
#
# Usage:
#   ./run_gpu_performance_tests.sh [category]
#
# Categories:
#   scaling     - GPU vs CPU scaling comparison
#   integrators - Integrator comparison (RK4, RK45, Boris)
#   threshold   - Threshold validation tests
#   long        - Long simulation efficiency
#   all         - Run all benchmarks (default)
################################################################################

set -euo pipefail

# Configuration
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PROJECT_ROOT="$(cd "$SCRIPT_DIR/../.." && pwd)"
BUILD_DIR="$PROJECT_ROOT/build"
CONFIG_DIR="$PROJECT_ROOT/validation/configs/performance/gpu"
RESULTS_DIR="$PROJECT_ROOT/validation/results/performance/gpu"
ICARION_BIN="$BUILD_DIR/src/icarion_main"

# Colors
RED='\033[0;31m'
GREEN='\033[0;32m'
YELLOW='\033[1;33m'
BLUE='\033[0;34m'
NC='\033[0m' # No Color

# Counters
TOTAL_TESTS=0
PASSED_TESTS=0
FAILED_TESTS=0

# Check if binary exists
if [ ! -f "$ICARION_BIN" ]; then
    echo -e "${RED}ERROR: icarion binary not found at $ICARION_BIN${NC}"
    echo "Please build the project first: cd build && make -j\$(nproc)"
    exit 1
fi

# Check if configs exist
if [ ! -d "$CONFIG_DIR" ]; then
    echo -e "${RED}ERROR: Config directory not found at $CONFIG_DIR${NC}"
    echo "Please generate configs first: cd validation/scripts && python3 generate_gpu_performance_configs.py"
    exit 1
fi

# Create results directory
mkdir -p "$RESULTS_DIR"

echo ""
echo "================================================================================"
echo "GPU Performance Benchmark Suite"
echo "================================================================================"
echo "Build:       $ICARION_BIN"
echo "Configs:     $CONFIG_DIR"
echo "Results:     $RESULTS_DIR"
echo "================================================================================"
echo ""

# Function to run a single benchmark
run_benchmark() {
    local config_file="$1"
    local category="$2"
    
    TOTAL_TESTS=$((TOTAL_TESTS + 1))
    
    local config_name=$(basename "$config_file" .json)
    local log_file="$RESULTS_DIR/${config_name}.log"
    
    echo -ne "${BLUE}[$TOTAL_TESTS]${NC} Running $config_name ... "
    
    # Run simulation
    if "$ICARION_BIN" "$config_file" > "$log_file" 2>&1; then
        echo -e "${GREEN}PASS${NC}"
        PASSED_TESTS=$((PASSED_TESTS + 1))
        return 0
    else
        echo -e "${RED}FAIL${NC}"
        FAILED_TESTS=$((FAILED_TESTS + 1))
        echo "  Log: $log_file"
        return 1
    fi
}

# Function to run benchmark category
run_category() {
    local pattern="$1"
    local category_name="$2"
    
    echo ""
    echo "================================================================================"
    echo "Category: $category_name"
    echo "================================================================================"
    echo ""
    
    local count=0
    for config in "$CONFIG_DIR"/$pattern; do
        if [ -f "$config" ]; then
            run_benchmark "$config" "$category_name"
            count=$((count + 1))
        fi
    done
    
    if [ $count -eq 0 ]; then
        echo -e "${YELLOW}WARNING: No configs found matching pattern: $pattern${NC}"
    fi
}

# Parse command line arguments
CATEGORY="${1:-all}"

case "$CATEGORY" in
    scaling)
        run_category "RK4_cpu_*.json" "GPU vs CPU Scaling (CPU)"
        run_category "RK4_gpu_*.json" "GPU vs CPU Scaling (GPU)"
        ;;
    integrators)
        run_category "integrator_*.json" "Integrator Comparison"
        ;;
    threshold)
        run_category "threshold_*.json" "Threshold Validation"
        ;;
    long)
        run_category "long_*.json" "Long Simulation Efficiency"
        ;;
    all)
        run_category "RK4_cpu_*.json" "GPU vs CPU Scaling (CPU)"
        run_category "RK4_gpu_*.json" "GPU vs CPU Scaling (GPU)"
        run_category "integrator_*.json" "Integrator Comparison"
        run_category "threshold_*.json" "Threshold Validation"
        run_category "long_*.json" "Long Simulation Efficiency"
        ;;
    *)
        echo -e "${RED}ERROR: Unknown category: $CATEGORY${NC}"
        echo "Valid categories: scaling, integrators, threshold, long, all"
        exit 1
        ;;
esac

# Summary
echo ""
echo "================================================================================"
echo "SUMMARY"
echo "================================================================================"
echo -e "Total tests:  $TOTAL_TESTS"
echo -e "Passed:       ${GREEN}$PASSED_TESTS${NC}"
echo -e "Failed:       ${RED}$FAILED_TESTS${NC}"
echo "================================================================================"
echo ""

if [ $FAILED_TESTS -eq 0 ]; then
    echo -e "${GREEN}All tests passed!${NC}"
    exit 0
else
    echo -e "${RED}Some tests failed. Check logs in $RESULTS_DIR${NC}"
    exit 1
fi
