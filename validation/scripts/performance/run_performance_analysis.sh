#!/bin/bash
# Unified performance analysis wrapper; GPU benchmarks are runtime-disabled in v1.0.0.

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
VALIDATION_DIR="$(cd "$SCRIPT_DIR/../.." && pwd)"
PERF_DIR="$VALIDATION_DIR/scripts/performance"
ANALYZER="$PERF_DIR/analyze_gpu_performance.py"
PYTHON_BIN=${PYTHON_BIN:-python3}

print_usage() {
  cat <<'EOF'
Usage: run_performance_analysis.sh [--python PATH]

Options:
  --python PATH   Explicit Python interpreter (defaults to python3 or $PYTHON_BIN)
  -h, --help      Show this help message

Note: GPU benchmarks are runtime-disabled in v1.0.0; this runner skips GPU analysis.
EOF
}

while [[ $# -gt 0 ]]; do
  case "$1" in
    --python)
      if [[ $# -lt 2 ]]; then
        echo "Error: --python requires a path" >&2
        exit 1
      fi
      PYTHON_BIN="$2"
      shift 2
      ;;
    -h|--help)
      print_usage
      exit 0
      ;;
    *)
      echo "Unknown argument: $1" >&2
      print_usage >&2
      exit 1
      ;;
  esac
 done

echo "GPU performance benchmarks are runtime-disabled for v1.0.0; skipping GPU performance analysis."
exit 0
