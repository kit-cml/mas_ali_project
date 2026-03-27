#!/bin/bash

# Compute all accepted TMS
Rscript compute_tms.R --original_metrics "data/metrics_nanion_scaled.csv" \
    --original_training "data/nanion_training_scaled.csv" \
    --original_testing "data/nanion_testing_scaled.csv" \
    --accepted_models_dir "accepted_models_nanion_latest" \
    --result_folder "results_tms_nanion" \
    --max_dimension 7 > log_compute_tms_nanion.txt &

wait