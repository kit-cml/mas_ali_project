#!/bin/bash

# Perform sensitivity analysis on all accepted TMS
MAX_JOBS=5  # Maximum parallel jobs
INPUT_FOLDER="results_tms_nanion"
OUTPUT_FOLDER="results_sensitivity_nanion"
LOG_DIR="logs_sensitivity_nanion"
NUM_DRUGS=10
CPR=0.75
mkdir -p "$OUTPUT_FOLDER/$LOG_DIR"

# Function to wait for available processing slots
wait_for_available_slot() {
    while [[ $(jobs -p | wc -l) -ge $MAX_JOBS ]]; do
        sleep 2  # Check every 2 seconds
    done
}

# Sensitivity Analysis Loop
for input_file in "$INPUT_FOLDER"/tms_metrics*; do
    filename=$(basename "$input_file")
    log_file="$OUTPUT_FOLDER/$LOG_DIR/$filename.log"
    r_script="sensitivity_analysis.R"
    wait_for_available_slot  # Ensure a slot is available
    Rscript "$r_script" --datadir "$input_file" \
        --dataset_name "nanion" \
        --outdir "$OUTPUT_FOLDER/$filename" \
        --influential_drugs $NUM_DRUGS \
        --cpr $CPR > "$log_file" 2>&1 &  # Run in background
done

# Wait for all sensitivity analysis jobs to finish
wait