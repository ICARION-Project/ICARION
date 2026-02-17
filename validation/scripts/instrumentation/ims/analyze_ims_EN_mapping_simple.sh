#!/bin/bash
# Analyze IMS E/N mapping results using the analyze_ims_drift.py output
# Parse simulation logs directly

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "$0")" && pwd)"
VALIDATION_DIR="$(cd "$SCRIPT_DIR/../../.." && pwd)"
REPO_ROOT="$(cd "$VALIDATION_DIR/.." && pwd)"

RUN_DIR="${ICARION_VALIDATION_RUN_DIR:-}"

pick_results_dir() {
    local candidate

    if [[ -n "$RUN_DIR" ]]; then
        candidate="$RUN_DIR/results/instruments/ims"
        if [[ -d "$candidate" ]]; then
            echo "$candidate"
            return 0
        fi
    fi

    candidate="$VALIDATION_DIR/results/v1.0.0_test/instruments/ims"
    if [[ -d "$candidate" ]]; then
        echo "$candidate"
        return 0
    fi

    candidate="$REPO_ROOT/results/v1.0.0_test/instruments/ims"
    if [[ -d "$candidate" ]]; then
        echo "$candidate"
        return 0
    fi

    echo "" 
    return 1
}

RESULTS_DIR="$(pick_results_dir)"
if [[ -z "$RESULTS_DIR" ]]; then
    echo "ERROR: Could not find IMS results directory." >&2
    echo "Tried (in order):" >&2
    echo "  - \$ICARION_VALIDATION_RUN_DIR/results/instruments/ims" >&2
    echo "  - $VALIDATION_DIR/results/v1.0.0_test/instruments/ims" >&2
    echo "  - $REPO_ROOT/results/v1.0.0_test/instruments/ims" >&2
    exit 1
fi

echo "================================================================================"
echo "IMS E/N Mapping Analysis (10-50 Td × 10,20,100,1000 Pa)"
echo "================================================================================"
echo ""

# Create output file
if [[ -n "$RUN_DIR" ]]; then
    OUTPUT="$RUN_DIR/logs/ims_EN_mapping_summary.txt"
else
    OUTPUT="$VALIDATION_DIR/logs/ims_EN_mapping_summary.txt"
fi
mkdir -p "$(dirname "$OUTPUT")"
echo "IMS E/N Mapping Results" > "$OUTPUT"
echo "Generated: $(date)" >> "$OUTPUT"
echo "================================================================================" >> "$OUTPUT"

# Models and parameters
MODELS="hss ehss friction"
EN_VALUES="10 20 30 40 50"
PRESSURES="10 20 100 1000"

for MODEL in $MODELS; do
    echo "" | tee -a "$OUTPUT"
    echo "================================================================================" | tee -a "$OUTPUT"
    echo "$MODEL Model" | tr '[:lower:]' '[:upper:]' | tee -a "$OUTPUT"
    echo "================================================================================" | tee -a "$OUTPUT"
    printf "%-10s %-10s %-15s %-15s %-12s %-8s\n" "E/N (Td)" "P (Pa)" "v_meas (m/s)" "v_exp (m/s)" "Error (%)" "Status" | tee -a "$OUTPUT"
    echo "--------------------------------------------------------------------------------" | tee -a "$OUTPUT"
    
    for EN in $EN_VALUES; do
        for P in $PRESSURES; do
            # Find matching HDF5 file
            H5_FILE=$(ls -1 "$RESULTS_DIR"/ims_${MODEL}_*${EN}Td*${P}Pa.h5 2>/dev/null | head -1)
            
            if [ -z "$H5_FILE" ]; then
                continue
            fi
            
            # Extract drift velocity from simulation log
            LOG_FILE="${RESULTS_DIR}/simulation.log"
            
            # For now, just report that file exists
            printf "%-10s %-10s %-15s %-15s %-12s %-8s\n" "$EN" "$P" "$(basename $H5_FILE)" "-" "-" "EXISTS" | tee -a "$OUTPUT"
        done
    done
done

echo "" | tee -a "$OUTPUT"
echo "================================================================================" | tee -a "$OUTPUT"
echo "Analysis complete. Output saved to: $OUTPUT"
echo "================================================================================" | tee -a "$OUTPUT"
echo ""
echo "Note: Full analysis requires h5py. Install with:"
echo "  python3 -m venv analysis_venv && source analysis_venv/bin/activate && pip install h5py numpy"
