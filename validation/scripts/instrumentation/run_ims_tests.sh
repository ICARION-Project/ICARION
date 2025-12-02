#!/bin/bash
# Run IMS validation tests (Session 2)
# Validates drift velocity vs E/N for 3 collision models × 3 E/N values × 3 species

SCRIPT_DIR="$(cd "$(dirname "$0")" && pwd)"
CONFIG_DIR="$SCRIPT_DIR/../../configs/instruments/ims"
RESULTS_DIR="$SCRIPT_DIR/../../results/ims_session_$(date +%Y%m%d_%H%M%S)"
ICARION_BIN="$SCRIPT_DIR/../../../build/src/icarion_main"

echo "=============================================="
echo "ICARION Validation Suite - Session 2"
echo "IMS Drift Velocity Tests"
echo "=============================================="
echo "Config dir: $CONFIG_DIR"
echo "Results dir: $RESULTS_DIR"
echo "Binary: $ICARION_BIN"
echo "Time: $(date)"
echo "=============================================="

# Check binary exists
if [ ! -f "$ICARION_BIN" ]; then
    echo "❌ ERROR: ICARION binary not found at $ICARION_BIN"
    echo "Please build first: cd build && make -j$(nproc)"
    exit 1
fi

# Check configs exist
if [ ! -d "$CONFIG_DIR" ]; then
    echo "❌ ERROR: Config directory not found: $CONFIG_DIR"
    echo "Please generate configs first:"
    echo "  python3 validation/scripts/instrumentation/generate_ims_configs.py"
    exit 1
fi

# Create results directory
mkdir -p "$RESULTS_DIR"

# Count configs
CONFIG_COUNT=$(ls "$CONFIG_DIR"/*.json 2>/dev/null | wc -l)
echo "Found $CONFIG_COUNT IMS configurations"
echo ""

# Run tests
PASSED=0
FAILED=0
TOTAL=0

for config in "$CONFIG_DIR"/*.json; do
    TOTAL=$((TOTAL + 1))
    basename=$(basename "$config" .json)
    
    echo "[$TOTAL/$CONFIG_COUNT] Running: $basename"
    
    # Run simulation with output capture
    output_dir="$RESULTS_DIR/$basename"
    mkdir -p "$output_dir"
    
    if "$ICARION_BIN" "$config" > "$output_dir/stdout.log" 2> "$output_dir/stderr.log"; then
        # Check if HDF5 output was created
        if [ -f "$output_dir/output.h5" ]; then
            echo "  ✅ PASS"
            PASSED=$((PASSED + 1))
        else
            echo "  ❌ FAIL (no HDF5 output)"
            FAILED=$((FAILED + 1))
        fi
    else
        echo "  ❌ FAIL (non-zero exit code)"
        FAILED=$((FAILED + 1))
        cat "$output_dir/stderr.log"
    fi
done

echo ""
echo "=============================================="
echo "IMS VALIDATION SUMMARY"
echo "=============================================="
echo "Total tests:  $TOTAL"
echo "✅ Passed:     $PASSED"
echo "❌ Failed:     $FAILED"
echo "Success rate: $(awk "BEGIN {printf \"%.1f\", ($PASSED/$TOTAL)*100}")%"
echo "Results dir:  $RESULTS_DIR"
echo "=============================================="

# Generate summary
if [ $PASSED -gt 0 ]; then
    echo ""
    echo "Next step: Analyze results"
    echo "  python3 validation/scripts/instrumentation/analyze_ims_drift.py $RESULTS_DIR"
fi

exit $FAILED
