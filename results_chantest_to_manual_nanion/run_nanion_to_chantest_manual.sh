#!/bin/bash

# Compute all accepted TMS
Rscript compute_tms.R --original_metrics "data/metrics_nanion_scaled.csv" \
    --original_training "data/nanion_training_scaled.csv" \
    --original_testing "data/nanion_testing_scaled.csv" \
    --accepted_models_dir "accepted_models_nanion_latest" \
    --result_folder "results_tms_nanion" \
    --max_dimension 7 > log_compute_tms_nanion.txt &

wait

# Perform sensitivity analysis on all accepted TMS
MAX_JOBS=12  # Maximum parallel jobs
INPUT_FOLDER="results_tms_nanion_all"
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

# Get the best model from model training and validation if required
# Rscript model_selection.R --input_folder "results_testing_manual" \
#     --results_filename "best_model_manual.csv" \
#     --max_input 11 > log_model_selection_manual.txt &

wait  # Wait for model selection to finish

# Evaluate performance from drug candidates and obtain calibration drugs
MAX_JOBS=12  # Reset counter for this section
TMS_FOLDER="results_tms_nanion_all"
SENSITIVITY_FOLDER="results_sensitivity_nanion"
OUTPUT_FOLDER="results_performance_evals_nanion"
LOG_DIR="logs_performance_evals"
mkdir -p "$OUTPUT_FOLDER"
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
        echo "Sens file: $sens_file"
    done
    echo "=========================================="
done

wait

# Get the calibration drugs
MAX_JOBS=12
CURRENT_JOBS=0
TMS_FOLDER="results_tms_nanion_all"
INPUT_FOLDER="results_performance_evals_nanion"
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

# Perform lab-specific validation of calibration drugs
# in model implementation lab (mil) == Chantest
MAX_JOBS=12
ACCEPTED_MODELS_FOLDER="accepted_models_nanion_latest"
SENSITIVITY_FOLDER="results_sensitivity_nanion"
GET_CALIBRATION_FOLDER="results_get_calibration_drugs"
TMS_FOLDER="results_tms_nanion_all"
OUTPUT_FOLDER="results_validate_calibration_drugs_nanion_to_chantest"
LOG_DIR="logs_validate_calibration_drugs_nanion_to_chantest"
mkdir $OUTPUT_FOLDER
for tms_metrics_file in "$TMS_FOLDER"/tms_metrics_*_*.csv; do
  base_tms_name=$(basename "$tms_metrics_file")
  tms_name=${base_tms_name%%.*}
  dimension=$(echo "$tms_name" | cut -d'_' -f3)
  pairid=$(echo "$tms_name" | cut -d'_' -f4)
  dataset_name="manual"
  sensitivity_file=$(echo "$SENSITIVITY_FOLDER/$base_tms_name"/*sens.csv)
  calibration_file=$(echo "$GET_CALIBRATION_FOLDER/$(basename "$tms_metrics_file")/calibration_drugs.csv")
  logfile=$(echo "$OUTPUT_FOLDER/$LOG_DIR/"$(basename "$base_tms_name")".log")
  base_tms_file=${base_tms_name/metrics_/}
  tms_file_input="$TMS_FOLDER/$base_tms_file"
  echo "tms_file_input = $tms_file_input"
  echo "Extracted dimension = $dimension"
  echo "Extracted pairid = $pairid"
  echo "dataset = $dataset_name"
  echo "sensitifity file =$sensitivity_file"
  echo "calibration file =$calibration_file"
  echo "logfile =$logfile"
  mkdir "$OUTPUT_FOLDER/$LOG_DIR"
  wait_for_available_slot
  Rscript validate_calibration_drugs.R -a "$ACCEPTED_MODELS_FOLDER" \
    --tms_dimension $dimension \
    --tms_pairid $pairid \
    --mdl_tms_file "$tms_file_input" \
    --mdl_sens_file "$sensitivity_file" \
    --calibration_drugs_file "$calibration_file" \
    --mil_training_file "data/chantest_training_scaled.csv" \
    --mil_testing_file "data/chantest_testing_scaled.csv" \
    --mil_metrics_file "data/metrics_chantest_scaled.csv" \
    --results_folder "$OUTPUT_FOLDER/$base_tms_name" > "$logfile" 2>&1 &
  echo "=============================="
done

wait

# Performance evaluation for the calibration drugs from model implementation lab
for tms_metrics_file in "$TMS_FOLDER"/tms_metrics_*_*.csv; do
  base_tms_name=$(basename "$tms_metrics_file")
  milth_file="$OUTPUT_FOLDER/$base_tms_name/mil_thresholds.csv"
  read th1 th2 < <(tail -n +2 "$milth_file" | tr -d '"' | tr ',' ' ')
  mil_drug_candidates_file="$OUTPUT_FOLDER/$base_tms_name/calibration_drugs.csv"
  mil_drug_size=$(tail -n +2 "$mil_drug_candidates_file" | wc -l)
  echo "Number of drugs = $mil_drug_size"
  mil_tms_file="$OUTPUT_FOLDER/$base_tms_name/mil_tms_metrics.csv"
  echo "MIL metrics file = $mil_tms_file"
  mil_simdir="$OUTPUT_FOLDER/$base_tms_name/mil_performance_evals"
  echo "MIL simdir = $mil_simdir"
  mil_logfile=$(echo "$OUTPUT_FOLDER/$LOG_DIR/"$(basename "$base_tms_name")"_evals.log")
  # Read th1 and th2 and trim whitespace
  IFS=$'\n' read -r -d '' thresholds < <(tail -n +2 "$milth_file" | tr -d '"' | tr ',' ' ')
  IFS=' ' read -r th1_raw th2_raw <<< "$thresholds"
  th1=$(echo "$th1_raw" | xargs) # Trim leading/trailing whitespace
  th2=$(echo "$th2_raw" | xargs) # Trim leading/trailing whitespace
  echo "th1 = $th1"
  echo "th2 = $th2"
  wait_for_available_slot
  Rscript.exe performance_evaluation.R --drug_candidate_file "$mil_drug_candidates_file" \
    --metrics_file "$mil_tms_file" \
    --drug_size $mil_drug_size \
    --simdir "$mil_simdir" \
    --evaluation_mode 1 \
    --th1="$th1" \
    --th2="$th2" > "$mil_logfile" 2>&1 &  # Run in background
  echo "=============================="
done

wait

# Perform lab-specific validation of calibration drugs
# in model implementation lab (mil) == Manual
MAX_JOBS=12
ACCEPTED_MODELS_FOLDER="accepted_models_nanion_latest"
SENSITIVITY_FOLDER="results_sensitivity_nanion"
GET_CALIBRATION_FOLDER="results_get_calibration_drugs"
TMS_FOLDER="results_tms_nanion_all"
OUTPUT_FOLDER="results_validate_calibration_drugs_nanion_to_manual"
LOG_DIR="logs_validate_calibration_drugs_nanion_to_manual"
mkdir $OUTPUT_FOLDER
for tms_metrics_file in "$TMS_FOLDER"/tms_metrics_*_*.csv; do
  base_tms_name=$(basename "$tms_metrics_file")
  tms_name=${base_tms_name%%.*}
  dimension=$(echo "$tms_name" | cut -d'_' -f3)
  pairid=$(echo "$tms_name" | cut -d'_' -f4)
  dataset_name="manual"
  sensitivity_file=$(echo "$SENSITIVITY_FOLDER/$base_tms_name"/*sens.csv)
  calibration_file=$(echo "$GET_CALIBRATION_FOLDER/$(basename "$tms_metrics_file")/calibration_drugs.csv")
  logfile=$(echo "$OUTPUT_FOLDER/$LOG_DIR/"$(basename "$base_tms_name")".log")
  base_tms_file=${base_tms_name/metrics_/}
  tms_file_input="$TMS_FOLDER/$base_tms_file"
  echo "tms_file_input = $tms_file_input"
  echo "Extracted dimension = $dimension"
  echo "Extracted pairid = $pairid"
  echo "dataset = $dataset_name"
  echo "sensitifity file =$sensitivity_file"
  echo "calibration file =$calibration_file"
  echo "logfile =$logfile"
  mkdir "$OUTPUT_FOLDER/$LOG_DIR"
  wait_for_available_slot
  Rscript validate_calibration_drugs.R -a "$ACCEPTED_MODELS_FOLDER" \
    --tms_dimension $dimension \
    --tms_pairid $pairid \
    --mdl_tms_file "$tms_file_input" \
    --mdl_sens_file "$sensitivity_file" \
    --calibration_drugs_file "$calibration_file" \
    --mil_training_file "data/manual_training_scaled.csv" \
    --mil_testing_file "data/manual_testing_scaled.csv" \
    --mil_metrics_file "data/metrics_manual_scaled.csv" \
    --results_folder "$OUTPUT_FOLDER/$base_tms_name" > "$logfile" 2>&1 &
  echo "=============================="
done

wait

# Performance evaluation for the calibration drugs from model implementation lab
for tms_metrics_file in "$TMS_FOLDER"/tms_metrics_*_*.csv; do
  base_tms_name=$(basename "$tms_metrics_file")
  milth_file="$OUTPUT_FOLDER/$base_tms_name/mil_thresholds.csv"
  read th1 th2 < <(tail -n +2 "$milth_file" | tr -d '"' | tr ',' ' ')
  mil_drug_candidates_file="$OUTPUT_FOLDER/$base_tms_name/calibration_drugs.csv"
  mil_drug_size=$(tail -n +2 "$mil_drug_candidates_file" | wc -l)
  echo "Number of drugs = $mil_drug_size"
  mil_tms_file="$OUTPUT_FOLDER/$base_tms_name/mil_tms_metrics.csv"
  echo "MIL metrics file = $mil_tms_file"
  mil_simdir="$OUTPUT_FOLDER/$base_tms_name/mil_performance_evals"
  echo "MIL simdir = $mil_simdir"
  mil_logfile=$(echo "$OUTPUT_FOLDER/$LOG_DIR/"$(basename "$base_tms_name")"_evals.log")
  # Read th1 and th2 and trim whitespace
  IFS=$'\n' read -r -d '' thresholds < <(tail -n +2 "$milth_file" | tr -d '"' | tr ',' ' ')
  IFS=' ' read -r th1_raw th2_raw <<< "$thresholds"
  th1=$(echo "$th1_raw" | xargs) # Trim leading/trailing whitespace
  th2=$(echo "$th2_raw" | xargs) # Trim leading/trailing whitespace
  echo "th1 = $th1"
  echo "th2 = $th2"
  wait_for_available_slot
  Rscript.exe performance_evaluation.R --drug_candidate_file "$mil_drug_candidates_file" \
    --metrics_file "$mil_tms_file" \
    --drug_size $mil_drug_size \
    --simdir "$mil_simdir" \
    --evaluation_mode 1 \
    --th1="$th1" \
    --th2="$th2" > "$mil_logfile" 2>&1 &  # Run in background
  echo "=============================="
done

wait 

echo "All processes finished."
