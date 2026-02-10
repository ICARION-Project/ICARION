#!/bin/bash
# Run IMS validation tests (Session 2)
# Validates drift velocity vs E/N for 3 collision models × 3 E/N values × 3 species

SCRIPT_DIR="$(cd "$(dirname "$0")" && pwd)"
VALIDATION_DIR="$(cd "$SCRIPT_DIR/../../.." && pwd)"
REPO_ROOT="$(cd "$VALIDATION_DIR/.." && pwd)"

CONFIG_DIR_DEFAULT="$VALIDATION_DIR/configs/instruments/ims"
RESULTS_ROOT_DEFAULT="$VALIDATION_DIR/results/instruments/ims"
ICARION_BIN_DEFAULT="$REPO_ROOT/build/src/icarion_main"

RUN_DIR="${ICARION_VALIDATION_RUN_DIR:-}"
if [[ -n "$RUN_DIR" ]]; then
    RESULTS_ROOT_DEFAULT="$RUN_DIR/results/instruments/ims"
fi

CONFIG_DIR="${INSTRUMENT_CONFIG_DIR_OVERRIDE:-$CONFIG_DIR_DEFAULT}"
RESULTS_ROOT="${INSTRUMENT_OUTPUT_ROOT_OVERRIDE:-$RESULTS_ROOT_DEFAULT}"
ICARION_BIN="${ICARION_BIN_OVERRIDE:-$ICARION_BIN_DEFAULT}"

CONFIG_DIR="$(cd "$CONFIG_DIR" && pwd)"
mkdir -p "$RESULTS_ROOT"
RESULTS_ROOT="$(cd "$RESULTS_ROOT" && pwd)"

if [[ -n "$RUN_DIR" ]]; then
    RESULTS_DIR="$RESULTS_ROOT"
else
    RESULTS_DIR="$RESULTS_ROOT/session_$(date +%Y%m%d_%H%M%S)"
fi

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

JOBS=${INSTRUMENT_JOBS:-${IMS_JOBS:-2}}
THREADS=${INSTRUMENT_THREADS:-${IMS_THREADS:-4}}
PASSED=0
FAILED=0
TOTAL=0
STATUS_LOG="$RESULTS_DIR/status.tsv"
: >"$STATUS_LOG"
declare -a RUNNING_PIDS=()
ACTIVE_JOBS=0

cleanup() {
    if [[ ${#RUNNING_PIDS[@]} -gt 0 ]]; then
        kill "${RUNNING_PIDS[@]}" 2>/dev/null || true
    fi
}
trap cleanup INT TERM

# Function to run a single test
run_test() {
    local config="$1"
    local total="$2"
    local basename=$(basename "$config" .json)
    
    echo "[$total/$CONFIG_COUNT] Running: $basename"
    
    # Run simulation with output capture
    output_dir="$RESULTS_DIR/$basename"
    mkdir -p "$output_dir"
    
    if (cd "$REPO_ROOT" && "$ICARION_BIN" "$config" --threads "$THREADS" > "$output_dir/stdout.log" 2> "$output_dir/stderr.log"); then
        echo "  ✅ PASS: $basename"
        return 0
    else
        echo "  ❌ FAIL (non-zero exit code): $basename"
        cat "$output_dir/stderr.log" | head -20
        return 1
    fi
}

export -f run_test
export ICARION_BIN RESULTS_DIR CONFIG_COUNT REPO_ROOT

# Collect all configs into array
CONFIGS=("$CONFIG_DIR"/*.json)

for idx in "${!CONFIGS[@]}"; do
    config="${CONFIGS[$idx]}"
    ordinal=$((idx + 1))
    (
        if run_test "$config" "$ordinal"; then
            printf '%s\t0\t%s\n' "$$" "$config" >>"$STATUS_LOG"
            exit 0
        else
            printf '%s\t1\t%s\n' "$$" "$config" >>"$STATUS_LOG"
            exit 1
        fi
    ) &
    pid=$!
    RUNNING_PIDS+=("$pid")
    ((ACTIVE_JOBS+=1))

    if [[ $ACTIVE_JOBS -ge $JOBS ]]; then
        wait -n || true
        ((ACTIVE_JOBS-=1)) || true
    fi
done

while [[ $ACTIVE_JOBS -gt 0 ]]; do
    wait -n || true
    ((ACTIVE_JOBS-=1)) || true
done

while IFS=$'\t' read -r pid status config; do
    if [[ -z "$pid" ]]; then
        continue
    fi
    TOTAL=$((TOTAL + 1))
    if [[ "$status" -eq 0 ]]; then
        ((PASSED++))
    else
        ((FAILED++))
    fi
done <"$STATUS_LOG" || true

rm -f "$STATUS_LOG"

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
