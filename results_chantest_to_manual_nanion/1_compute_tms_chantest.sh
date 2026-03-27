#!/bin/bash

# Compute all accepted TMS
Rscript compute_tms.R --original_metrics "data/metrics_chantest_scaled.csv" \
    --original_training "data/chantest_training_scaled.csv" \
    --original_testing "data/chantest_testing_scaled.csv" \
    --accepted_models_dir "accepted_models_chantest_latest" \
    --result_folder "results_tms_chantest" \
    --max_dimension 10 > log_compute_tms_chantest.txt &

wait