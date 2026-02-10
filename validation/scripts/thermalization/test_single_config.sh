#!/bin/bash
# Single thermalization config runner (invokes icarion_main)
# Usage: ./test_single_config.sh path/to/config.json

set -euo pipefail

if [[ $# -ne 1 ]]; then
    echo "Usage: $0 <config.json>" >&2
    exit 1
fi

CONFIG="$1"
if [[ ! -f "$CONFIG" ]]; then
    echo "Error: Config file '$CONFIG' not found" >&2
    exit 1
fi

ICARION_ROOT="$(cd "$(dirname "$0")/../../.." && pwd)"
ICARION_BIN="${ICARION_BIN:-$ICARION_ROOT/build/src/icarion_main}"
THERM_THREADS=${THERM_THREADS:-1}
export OMP_NUM_THREADS="$THERM_THREADS"
if [[ ! -x "$ICARION_BIN" ]]; then
    echo "Error: icarion_main not found or not executable at '$ICARION_BIN'" >&2
    exit 1
fi

TEST_NAME=$(basename "$CONFIG" .json)
if [[ -n "${ICARION_VALIDATION_RUN_DIR:-}" ]]; then
    OUTPUT_DIR="$ICARION_VALIDATION_RUN_DIR/results/physics/thermalization_single/${TEST_NAME}_$(date +%Y%m%d_%H%M%S)"
else
    OUTPUT_DIR="$ICARION_ROOT/validation/results/thermalization_single/${TEST_NAME}_$(date +%Y%m%d_%H%M%S)"
fi
mkdir -p "$OUTPUT_DIR"
cp "$CONFIG" "$OUTPUT_DIR/config.json"

echo "==========================================="
echo "Thermalization single test"
echo "Config   : $CONFIG"
echo "Output   : $OUTPUT_DIR"
echo "Binary   : $ICARION_BIN"
echo "Threads  : $THERM_THREADS"
echo "Timestamp: $(date)"
echo "==========================================="

START_TIME=$(date +%s)
if (cd "$ICARION_ROOT" && "$ICARION_BIN" "$CONFIG" --threads "$THERM_THREADS" >"$OUTPUT_DIR/simulation.log" 2>&1); then
    END_TIME=$(date +%s)
    RUNTIME=$((END_TIME - START_TIME))
    echo "✅ SUCCESS (${RUNTIME}s)"
    echo "Status: PASS" >"$OUTPUT_DIR/status.txt"
    echo "Runtime: ${RUNTIME}s" >>"$OUTPUT_DIR/status.txt"
    echo "Timestamp: $(date)" >>"$OUTPUT_DIR/status.txt"
    exit 0
else
    END_TIME=$(date +%s)
    RUNTIME=$((END_TIME - START_TIME))
    echo "❌ FAILED (${RUNTIME}s)"
    echo "Status: FAIL" >"$OUTPUT_DIR/status.txt"
    echo "Runtime: ${RUNTIME}s" >>"$OUTPUT_DIR/status.txt"
    echo "Timestamp: $(date)" >>"$OUTPUT_DIR/status.txt"
    tail -40 "$OUTPUT_DIR/simulation.log" >&2 || true
    exit 1
fi
