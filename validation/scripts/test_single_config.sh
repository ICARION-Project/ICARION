#!/bin/bash
CONFIG=$1
ICARION_ROOT="$(cd "$(dirname "$0")/../.." && pwd)"
OUTPUT_DIR="results/test_$(date +%Y%m%d_%H%M%S)"
echo "Testing: $CONFIG"
mkdir -p "$OUTPUT_DIR"
"$ICARION_ROOT/build/src/icarion_main" "$CONFIG" > "$OUTPUT_DIR/simulation.log" 2>&1
[ $? -eq 0 ] && echo "✅ SUCCESS" || echo "❌ FAILED"
ls -lh "$OUTPUT_DIR/"
