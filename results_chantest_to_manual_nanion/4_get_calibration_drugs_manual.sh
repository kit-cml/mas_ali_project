#!/bin/bash

# Get the calibration drugs
MAX_JOBS=5
CURRENT_JOBS=0
TMS_FOLDER="results_tms_manual"
INPUT_FOLDER="results_performance_evals_manual"
OUTPUT_FOLDER="results_get_calibration_drugs"
LOG_DIR="logs_calibration_drugs"
mkdir $OUTPUT_FOLDER
for tms_file in "$TMS_FOLDER"/tms_metrics*; do
    tms_folder=$(echo "$INPUT_FOLDER/$(basename "$tms_file")")
    logfile=$(echo "$OUTPUT_FOLDER/$LOG_DIR/calibration_drugs_"$(basename "$tms_file")".log")
    mkdir "$OUTPUT_FOLDER/$LOG_DIR"
    Rscript.exe get_calibration_drugs.R --input_dir "$tms_folder/$(basename "$tms_file")" \
        --output_dir "$OUTPUT_FOLDER/$(basename "$tms_file")" > "$logfile" &
    ((CURRENT_JOBS++))
    # If max parallel jobs reached, wait for processes to complete
    if [[ $CURRENT_JOBS -ge $MAX_JOBS ]]; then
        wait
        CURRENT_JOBS=0
    fi
    echo "=========================================="
done

# Final wait to ensure all processes complete
wait
