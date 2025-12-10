#!/bin/bash
# Run IMS validation tests (Session 2)
# Validates drift velocity vs E/N for 3 collision models × 3 E/N values × 3 species

SCRIPT_DIR="$(cd "$(dirname "$0")" && pwd)"
VALIDATION_DIR="$(cd "$SCRIPT_DIR/../../.." && pwd)"
REPO_ROOT="$(cd "$VALIDATION_DIR/.." && pwd)"
CONFIG_DIR="$VALIDATION_DIR/configs/instruments/ims"
RESULTS_DIR="$VALIDATION_DIR/results/ims_session_$(date +%Y%m%d_%H%M%S)"
ICARION_BIN="$REPO_ROOT/build/src/icarion_main"

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

# Run tests in parallel (2 at a time, 4 threads each)
PASSED=0
FAILED=0
TOTAL=0

# Function to run a single test
run_test() {
    local config="$1"
    local total="$2"
    local basename=$(basename "$config" .json)
    
    echo "[$total/$CONFIG_COUNT] Running: $basename"
    
    # Run simulation with output capture
    output_dir="$RESULTS_DIR/$basename"
    mkdir -p "$output_dir"
    
    if "$ICARION_BIN" "$config" --threads 4 > "$output_dir/stdout.log" 2> "$output_dir/stderr.log"; then
        echo "  ✅ PASS: $basename"
        return 0
    else
        echo "  ❌ FAIL (non-zero exit code): $basename"
        cat "$output_dir/stderr.log" | head -20
        return 1
    fi
}

export -f run_test
export ICARION_BIN RESULTS_DIR CONFIG_COUNT

# Collect all configs into array
CONFIGS=("$CONFIG_DIR"/*.json)

# Run tests 2 at a time
for ((i=0; i<${#CONFIGS[@]}; i+=2)); do
    TOTAL=$((i + 1))
    
    # Start first job in background
    run_test "${CONFIGS[i]}" "$TOTAL" &
    pid1=$!
    
    # Start second job if available
    if [ $((i+1)) -lt ${#CONFIGS[@]} ]; then
        TOTAL=$((i + 2))
        run_test "${CONFIGS[$((i+1))]}" "$TOTAL" &
        pid2=$!
        
        # Wait for both jobs
        wait $pid1 && ((PASSED++)) || ((FAILED++))
        wait $pid2 && ((PASSED++)) || ((FAILED++))
    else
        # Only one job left
        wait $pid1 && ((PASSED++)) || ((FAILED++))
    fi
done

TOTAL=${#CONFIGS[@]}

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
