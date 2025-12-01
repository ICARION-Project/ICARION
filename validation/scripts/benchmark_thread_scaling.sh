#!/bin/bash
# ICARION Thread Scaling Benchmark
# Tests CPU performance with different thread counts

set -e

CONFIG="configs/performance/thread_scaling_benchmark.json"
RESULTS_DIR="./results/thread_scaling"

echo "============================================"
echo "ICARION CPU Thread Scaling Benchmark"
echo "============================================"
echo "Config: $CONFIG"
echo "Date: $(date)"
echo ""

# Create results directory
mkdir -p "$RESULTS_DIR"

# Detect system info
if [ -f /proc/cpuinfo ]; then
    TOTAL_CORES=$(grep -c ^processor /proc/cpuinfo)
    echo "System: $TOTAL_CORES CPU cores detected"
    
    # Check for NUMA
    if command -v numactl &> /dev/null; then
        NUMA_NODES=$(numactl --hardware | grep "available:" | awk '{print $2}')
        echo "NUMA: $NUMA_NODES nodes"
    fi
fi
echo ""

# Benchmark function
benchmark_threads() {
    local threads=$1
    local label=$2
    
    export OMP_NUM_THREADS=$threads
    
    echo -n "Testing OMP_NUM_THREADS=$threads ($label)... "
    
    # Run simulation and capture timing
    START_TIME=$(date +%s.%N)
    ../build/src/icarion_main "$CONFIG" > /dev/null 2>&1
    END_TIME=$(date +%s.%N)
    
    # Calculate elapsed time
    ELAPSED=$(echo "$END_TIME - $START_TIME" | bc)
    
    # Calculate speedup (relative to 1 thread)
    if [ "$threads" -eq 1 ]; then
        BASELINE=$ELAPSED
        SPEEDUP="1.00"
        EFFICIENCY="100.0"
    else
        SPEEDUP=$(echo "scale=2; $BASELINE / $ELAPSED" | bc)
        EFFICIENCY=$(echo "scale=1; ($SPEEDUP / $threads) * 100" | bc)
    fi
    
    echo "${ELAPSED}s (${SPEEDUP}× speedup, ${EFFICIENCY}% efficiency)"
    
    # Log to results file
    echo "$threads,$ELAPSED,$SPEEDUP,$EFFICIENCY" >> "$RESULTS_DIR/benchmark_$(date +%Y%m%d_%H%M%S).csv"
}

# Run benchmarks
echo "Results:"
echo "--------"

# Initialize CSV
TIMESTAMP=$(date +%Y%m%d_%H%M%S)
CSV_FILE="$RESULTS_DIR/benchmark_${TIMESTAMP}.csv"
echo "threads,time_s,speedup,efficiency_percent" > "$CSV_FILE"

# Test different thread counts
benchmark_threads 1 "baseline"
benchmark_threads 2 "dual-core"
benchmark_threads 4 "quad-core"

if [ "$TOTAL_CORES" -ge 8 ]; then
    benchmark_threads 8 "8-core"
fi

if [ "$TOTAL_CORES" -ge 16 ]; then
    benchmark_threads 16 "16-core"
fi

echo ""
echo "============================================"
echo "Results saved to: $CSV_FILE"
echo ""
echo "Summary:"
echo "--------"
cat "$CSV_FILE" | column -t -s,

echo ""
echo "Expected CPU scaling (memory bandwidth limited):"
echo "  1 thread:  1.0× (baseline)"
echo "  2 threads: 1.8-1.9× (90-95% efficiency)"
echo "  4 threads: 3.0-3.5× (75-85% efficiency)"
echo "  8 threads: 4.0-5.0× (50-60% efficiency)"
echo "  16+ threads: <6× (memory bandwidth saturates)"
