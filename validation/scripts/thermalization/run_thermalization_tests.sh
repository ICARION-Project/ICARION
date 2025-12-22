#!/bin/bash
# Run thermalization validation tests (Session 2)
# Usage: ./run_thermalization_tests.sh [quick|subset|full]

MODE="${1:-quick}"
SCRIPT_DIR="$(cd "$(dirname "$0")" && pwd)"
VALIDATION_DIR="$(cd "$SCRIPT_DIR/../.." && pwd)"
CONFIG_DIR="$VALIDATION_DIR/configs/physics/thermalization"
RESULTS_DIR="$SCRIPT_DIR/../results/thermalization_session_$(date +%Y%m%d_%H%M%S)"
THERM_JOBS=${THERM_JOBS:-1}

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

TMP_DIR="$(mktemp -d "$RESULTS_DIR/tmp.XXXX")"
declare -a PIDS=()

for i in "${!CONFIGS[@]}"; do
    CONFIG="${CONFIGS[$i]}"
    CONFIG_NAME=$(basename "$CONFIG" .json)
    {
        TEST_START=$(date +%s)
        if cd "$SCRIPT_DIR" && ./test_single_config.sh "$CONFIG"; then
            TEST_END=$(date +%s)
            RUNTIME=$((TEST_END - TEST_START))
            STATUS=0
            LINE="| $CONFIG_NAME | ✅ PASS | ${RUNTIME}s | TBD | - |"
        else
            TEST_END=$(date +%s)
            RUNTIME=$((TEST_END - TEST_START))
            STATUS=1
            LINE="| $CONFIG_NAME | ❌ FAIL | ${RUNTIME}s | - | Check logs |"
        fi
        echo "$STATUS $RUNTIME" > "$TMP_DIR/${i}_${CONFIG_NAME}.status"
        echo "$LINE" > "$TMP_DIR/${i}_${CONFIG_NAME}.line"
    } &
    PIDS+=("$!")

    while [[ ${#PIDS[@]} -ge $THERM_JOBS && $THERM_JOBS -gt 0 ]]; do
        wait -n || true
        # prune finished
        NEW_PIDS=()
        for pid in "${PIDS[@]}"; do
            if kill -0 "$pid" 2>/dev/null; then
                NEW_PIDS+=("$pid")
            fi
        done
        PIDS=("${NEW_PIDS[@]}")
    done
done

# wait for remaining
for pid in "${PIDS[@]:-}"; do
    wait "$pid" || true
done

# Collect results
for status_file in "$TMP_DIR"/*.status; do
    [[ -f "$status_file" ]] || continue
    read -r STATUS RUNTIME < "$status_file"
    base=${status_file%.status}
    line_file="${base}.line"
    if [[ "$STATUS" -eq 0 ]]; then
        PASSED=$((PASSED + 1))
    else
        FAILED=$((FAILED + 1))
    fi
    [[ -f "$line_file" ]] && cat "$line_file" >> "$RESULTS_DIR/summary.md"
done

rm -rf "$TMP_DIR"

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
