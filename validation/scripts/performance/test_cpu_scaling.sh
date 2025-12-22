#!/bin/bash
# Test OpenMP CPU core scaling with different thread counts

CONFIG="$1"
RESULTS_DIR="results/cpu_scaling_$(date +%Y%m%d_%H%M%S)"

if [ -z "$CONFIG" ]; then
    echo "Usage: $0 <config.json>"
    exit 1
fi

SCRIPT_DIR="$(cd "$(dirname "$0")" && pwd)"
ICARION_ROOT="$(cd "$SCRIPT_DIR/../../.." && pwd)"
ICARION_BIN="$ICARION_ROOT/build/src/icarion_main"

mkdir -p "$RESULTS_DIR"

echo "=========================================="
echo "OpenMP CPU Core Scaling Test"
echo "=========================================="
echo "Config: $CONFIG"
echo "Executable: $ICARION_BIN"
echo "Results: $RESULTS_DIR"
echo ""

# Test with different thread counts
for THREADS in 1 2 4 8 16 32; do
    echo "----------------------------------------"
    echo "Testing with OMP_NUM_THREADS=$THREADS"
    echo "----------------------------------------"
    
    export OMP_NUM_THREADS=$THREADS
    
    LOG_FILE="$RESULTS_DIR/threads_${THREADS}.log"
    TIME_FILE="$RESULTS_DIR/threads_${THREADS}.time"
    
    # Create unique output directory for each thread count to avoid HDF5 file locking
    THREAD_OUTPUT_DIR="$RESULTS_DIR/output_threads_${THREADS}"
    mkdir -p "$THREAD_OUTPUT_DIR"
    
    # Create temporary config with unique output path
    TEMP_CONFIG="$RESULTS_DIR/config_threads_${THREADS}.json"
    sed "s|scaling_cpu_cores.h5|scaling_threads_${THREADS}.h5|g" "$CONFIG" > "$TEMP_CONFIG"
    
    # Run simulation and capture timing
    /usr/bin/time -v "$ICARION_BIN" "$TEMP_CONFIG" > "$LOG_FILE" 2>&1
    EXIT_CODE=$?
    
    # Extract CPU time from log
    CPU_TIME=$(sed -n 's/.*CPU time:[[:space:]]*\([0-9.eE+-]\+\).*/\1/p' "$LOG_FILE" | tail -1)
    
    # Extract real/user/sys time from time command
    REAL_TIME=$(grep "Elapsed" "$LOG_FILE" | tail -1 | awk '{print $NF}')
    
    echo "  Exit code: $EXIT_CODE"
    echo "  CPU time (ICARION): ${CPU_TIME}s"
    echo "  Real time (wall): $REAL_TIME"
    echo ""
    
    # Store result
    echo "$THREADS,$CPU_TIME,$REAL_TIME,$EXIT_CODE" >> "$RESULTS_DIR/scaling_summary.csv"
done

echo "=========================================="
echo "Scaling Test Complete"
echo "=========================================="
echo ""
echo "Results saved to: $RESULTS_DIR/scaling_summary.csv"
echo ""
cat "$RESULTS_DIR/scaling_summary.csv"
