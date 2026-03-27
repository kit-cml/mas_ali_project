#!/bin/bash

# Perform lab-specific validation of calibration drugs
# in model implementation lab (mil) == Nanion
MAX_JOBS=8
ACCEPTED_MODELS_FOLDER="accepted_models_chantest_latest"
SENSITIVITY_FOLDER="results_sensitivity_chantest"
GET_CALIBRATION_FOLDER="results_get_calibration_drugs"
TMS_FOLDER="results_tms_chantest"
OUTPUT_FOLDER="results_validate_calibration_drugs_chantest_to_nanion"
LOG_DIR="logs_validate_calibration_drugs_chantest_to_nanion"
mkdir $OUTPUT_FOLDER

# Function to wait for available processing slots
wait_for_available_slot() {
    while [[ $(jobs -p | wc -l) -ge $MAX_JOBS ]]; do
        sleep 2  # Check every 2 seconds
    done
}

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
    --mil_training_file "data/nanion_training_scaled.csv" \
    --mil_testing_file "data/nanion_testing_scaled.csv" \
    --mil_metrics_file "data/metrics_nanion_scaled.csv" \
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
MAX_JOBS=8
ACCEPTED_MODELS_FOLDER="accepted_models_chantest_latest"
SENSITIVITY_FOLDER="results_sensitivity_chantest"
GET_CALIBRATION_FOLDER="results_get_calibration_drugs"
TMS_FOLDER="results_tms_chantest"
OUTPUT_FOLDER="results_validate_calibration_drugs_chantest_to_manual"
LOG_DIR="logs_validate_calibration_drugs_chantest_to_manual"
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