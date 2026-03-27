#!/bin/bash

# Compute all accepted TMS for manual dataset
Rscript compute_tms.R --original_metrics "data/metrics_manual_scaled.csv" \
  --original_training "data/manual_training_scaled.csv" \
  --original_testing "data/manual_testing_scaled.csv" \
  --accepted_models_dir "accepted_models_manual_latest" \
  --result_folder "results_tms_manual" \
  --max_dimension 10 > log_compute_tms_manual.txt &

wait