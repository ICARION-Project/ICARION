#!/bin/bash
# Run thermalization validation tests (Session 2)
# Usage: ./run_thermalization_tests.sh [quick|subset|full]

MODE="${1:-quick}"
SCRIPT_DIR="$(cd "$(dirname "$0")" && pwd)"
CONFIG_DIR="$SCRIPT_DIR/../configs/physics/thermalization"
RESULTS_DIR="$SCRIPT_DIR/../results/thermalization_session_$(date +%Y%m%d_%H%M%S)"

echo "=============================================="
echo "ICARION Validation Suite - Session 2"
echo "Thermalization Tests"
echo "=============================================="
echo "Mode: $MODE"
echo "Config dir: $CONFIG_DIR"
echo "Results dir: $RESULTS_DIR"
echo "Time: $(date)"
echo "=============================================="

# Create results directory
mkdir -p "$RESULTS_DIR"

# Count configs
TOTAL_CONFIGS=$(find "$CONFIG_DIR" -name "*.json" | wc -l)
echo "Found $TOTAL_CONFIGS thermalization configs"

# Select test subset based on mode
case "$MODE" in
    "quick")
        # Quick test: Only room temperature, medium pressure, one ion, both models
        PATTERN="*300K_20.0Pa.json"
        echo "Quick mode: Testing room temperature (300K) at 20 Pa only"
        ;;
    "subset")
        # Subset test: All temperatures, medium pressure, one ion, both models  
        PATTERN="*H3Op*20.0Pa.json"
        echo "Subset mode: Testing H3O+ at all temperatures at 20 Pa"
        ;;
    "full")
        # Full test: All 90 configs
        PATTERN="*.json"
        echo "Full mode: Testing all $TOTAL_CONFIGS configs"
        ;;
    *)
        echo "Error: Unknown mode '$MODE'. Use: quick, subset, or full"
        exit 1
        ;;
esac

# Find matching configs
CONFIGS=($(find "$CONFIG_DIR" -name "$PATTERN" | sort))
NUM_CONFIGS=${#CONFIGS[@]}

echo "Selected $NUM_CONFIGS configs for testing"
echo ""

if [ $NUM_CONFIGS -eq 0 ]; then
    echo "No configs found matching pattern: $PATTERN"
    exit 1
fi

# Initialize summary
echo "# Thermalization Test Results" > "$RESULTS_DIR/summary.md"
echo "" >> "$RESULTS_DIR/summary.md"
echo "**Session:** 2 - Thermalization Tests" >> "$RESULTS_DIR/summary.md"
echo "**Mode:** $MODE" >> "$RESULTS_DIR/summary.md"
echo "**Configs:** $NUM_CONFIGS / $TOTAL_CONFIGS" >> "$RESULTS_DIR/summary.md"
echo "**Started:** $(date)" >> "$RESULTS_DIR/summary.md"
echo "" >> "$RESULTS_DIR/summary.md"
echo "| Config | Status | Runtime | Temperature Error | Notes |" >> "$RESULTS_DIR/summary.md"
echo "|--------|--------|---------|-------------------|-------|" >> "$RESULTS_DIR/summary.md"

# Run tests
PASSED=0
FAILED=0
START_SESSION=$(date +%s)

for i in "${!CONFIGS[@]}"; do
    CONFIG="${CONFIGS[$i]}"
    CONFIG_NAME=$(basename "$CONFIG" .json)
    
    echo "----------------------------------------"
    echo "Test $((i+1))/$NUM_CONFIGS: $CONFIG_NAME"
    echo "----------------------------------------"
    
    # Run single test
    TEST_START=$(date +%s)
    
    if cd "$SCRIPT_DIR" && ./test_single_config.sh "$CONFIG"; then
        TEST_END=$(date +%s)
        RUNTIME=$((TEST_END - TEST_START))
        PASSED=$((PASSED + 1))
        
        # TODO: Extract temperature analysis from HDF5 output
        # For now, just mark as PASS
        echo "| $CONFIG_NAME | ✅ PASS | ${RUNTIME}s | TBD | - |" >> "$RESULTS_DIR/summary.md"
        
    else
        TEST_END=$(date +%s)
        RUNTIME=$((TEST_END - TEST_START))
        FAILED=$((FAILED + 1))
        
        echo "| $CONFIG_NAME | ❌ FAIL | ${RUNTIME}s | - | Check logs |" >> "$RESULTS_DIR/summary.md"
    fi
    
    echo ""
done

# Session summary
END_SESSION=$(date +%s)
SESSION_RUNTIME=$((END_SESSION - START_SESSION))

echo "=============================================="
echo "Session 2 - Thermalization Tests Complete"
echo "=============================================="
echo "Passed: $PASSED"
echo "Failed: $FAILED" 
echo "Total:  $NUM_CONFIGS"
echo "Runtime: ${SESSION_RUNTIME}s"
echo "Success rate: $(echo "scale=1; $PASSED * 100 / $NUM_CONFIGS" | bc -l)%"
echo ""
echo "Results saved to: $RESULTS_DIR"

# Update summary
echo "" >> "$RESULTS_DIR/summary.md"
echo "**Completed:** $(date)" >> "$RESULTS_DIR/summary.md"
echo "**Runtime:** ${SESSION_RUNTIME}s" >> "$RESULTS_DIR/summary.md"
echo "**Results:** $PASSED passed, $FAILED failed ($(echo "scale=1; $PASSED * 100 / $NUM_CONFIGS" | bc -l)% success rate)" >> "$RESULTS_DIR/summary.md"

if [ $FAILED -gt 0 ]; then
    echo "⚠️  Some tests failed. Check individual logs in $RESULTS_DIR"
    exit 1
else
    echo "✅ All tests passed!"
    exit 0
fi