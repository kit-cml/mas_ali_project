#!/bin/bash

# Evaluate performance from drug candidates and obtain calibration drugs
MAX_JOBS=5  # Reset counter for this section
TMS_FOLDER="results_tms_nanion"
SENSITIVITY_FOLDER="results_sensitivity_nanion"
OUTPUT_FOLDER="results_performance_evals_nanion"
LOG_DIR="logs_performance_evals"
mkdir -p "$OUTPUT_FOLDER"

# Function to wait for available processing slots
wait_for_available_slot() {
    while [[ $(jobs -p | wc -l) -ge $MAX_JOBS ]]; do
        sleep 2  # Check every 2 seconds
    done
}

for tms_file in "$TMS_FOLDER"/tms_metrics*; do
    base_filename=$(basename "$tms_file")
    drug_candidates_file=$(echo "$SENSITIVITY_FOLDER/$base_filename"/*candidates.csv)
    sens_file=$(echo "$SENSITIVITY_FOLDER/$base_filename"/*sens.csv)
    echo "Processing: $tms_file"
    echo "Drug Candidates: $drug_candidates_file"
    echo "Sensitivity File: $sens_file"
    for i in {6..12}; do
        tms_folder="$OUTPUT_FOLDER/$base_filename"
        mkdir -p "$tms_folder/$LOG_DIR"
        logfile="$tms_folder/$LOG_DIR/performance_evals_$i.log"
        echo "Output Folder: $tms_folder"
        echo "Log File: $logfile"
        wait_for_available_slot  # Ensure a slot is available
        Rscript.exe performance_evaluation.R --drug_candidate_file "$drug_candidates_file" \
            --metrics_file "$tms_file" \
            --drug_size $i \
            --simdir "$tms_folder/$base_filename" \
            --evaluation_mode 1 \
            --sensitivity_file "$sens_file" > "$logfile" 2>&1 &  # Run in background
            # --th1 2.529741346 \
            # --th2 -4.401176331 > "$logfile" 2>&1 &  # Run in background
        echo "Sens file: $sens_file"
    done
    echo "=========================================="
done

wait