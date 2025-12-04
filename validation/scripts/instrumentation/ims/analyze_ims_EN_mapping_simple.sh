#!/bin/bash
# Analyze IMS E/N mapping results using the analyze_ims_drift.py output
# Parse simulation logs directly

RESULTS_DIR="/home/chsch95/ICARION/results/v1.0_test/instruments/ims"

echo "================================================================================"
echo "IMS E/N Mapping Analysis (10-50 Td × 10,20,100,1000 Pa)"
echo "================================================================================"
echo ""

# Create output file
OUTPUT="validation/results/ims_EN_mapping_summary.txt"
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
