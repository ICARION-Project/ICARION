#!/bin/bash
# Single config test helper for validation suite
# Usage: ./test_single_config.sh path/to/config.json

if [ $# -ne 1 ]; then
    echo "Usage: $0 <config.json>"
    exit 1
fi

CONFIG="$1"
if [ ! -f "$CONFIG" ]; then
    echo "Error: Config file '$CONFIG' not found"
    exit 1
fi

# Extract test name from config path
TEST_NAME=$(basename "$CONFIG" .json)
ICARION_ROOT="$(cd "$(dirname "$0")/../.." && pwd)"
OUTPUT_DIR="results/validation_$(date +%Y%m%d_%H%M%S)_${TEST_NAME}"

echo "==========================================="
echo "ICARION Validation Suite - Single Test"
echo "==========================================="
echo "Config: $CONFIG"
echo "Output: $OUTPUT_DIR"
echo "Time:   $(date)"
echo "==========================================="

# Create output directory
mkdir -p "$OUTPUT_DIR"

# Copy config for reproducibility
cp "$CONFIG" "$OUTPUT_DIR/config.json"

# Run simulation
echo "Starting simulation..."
START_TIME=$(date +%s)

if "$ICARION_ROOT/build/src/icarion_main" "$CONFIG" --threads 4 > "$OUTPUT_DIR/simulation.log" 2>&1; then
    END_TIME=$(date +%s)
    RUNTIME=$((END_TIME - START_TIME))
    
    echo "✅ SUCCESS (${RUNTIME}s)"
    echo "Status: PASS" > "$OUTPUT_DIR/status.txt"
    echo "Runtime: ${RUNTIME}s" >> "$OUTPUT_DIR/status.txt"
    echo "Timestamp: $(date)" >> "$OUTPUT_DIR/status.txt"
    
    # List output files
    echo ""
    echo "Output files:"
    ls -lh "$OUTPUT_DIR/"
    
    # Show simulation summary
    if [ -f "$OUTPUT_DIR/simulation.log" ]; then
        echo ""
        echo "Simulation Summary:"
        tail -15 "$OUTPUT_DIR/simulation.log"
    fi
    
    echo ""
    echo "✅ Test completed successfully"
    exit 0
else
    END_TIME=$(date +%s)
    RUNTIME=$((END_TIME - START_TIME))
    
    echo "❌ FAILED (${RUNTIME}s)"
    echo "Status: FAIL" > "$OUTPUT_DIR/status.txt"
    echo "Runtime: ${RUNTIME}s" >> "$OUTPUT_DIR/status.txt"
    echo "Timestamp: $(date)" >> "$OUTPUT_DIR/status.txt"
    
    # Show error log
    echo ""
    echo "Error log:"
    tail -20 "$OUTPUT_DIR/simulation.log"
    
    echo ""
    echo "❌ Test failed"
    exit 1
fi
