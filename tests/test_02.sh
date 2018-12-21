#!/usr/bin/env bash
# Regression test 02.
# Uses simulation parameters and data from the original SOFI3D code,
# from the benchmark `fullspace`.
. tests/functions.sh

MODEL="src/model_elastic.c"
TEST_PATH="tests/fixtures/test_02"
TEST_ID="TEST_02"

# Setup function prepares environment for the test (creates directories).
setup

# Preserve old model.
mv $MODEL ${MODEL}.bak.${TEST_ID}

on_exit() {
    # Clean up function to be called at the end of the test.
    mv ${MODEL}.bak.${TEST_ID} $MODEL
}

# Execute function 'on_exit' when this script exits.
trap on_exit INT TERM EXIT

# Copy test model.
cp "${TEST_PATH}/model_elastic.c"                  src/
cp "${TEST_PATH}/in_and_out/fullspace.json"        tmp/in_and_out/sofi3D.json
cp "${TEST_PATH}/sources/fullspace_sources.dat"    tmp/sources/

compile_code

# Run code.
log "Running solver. Output is captured to tmp/ASOFI3D.log"
./run_ASOFI3D.sh 16 tmp/ > tmp/ASOFI3D.log &
task_id=$!
animate_progress $task_id "TEST_02: Running solver"

code=$?
if [ "$code" -ne "0" ]; then
    error "Running ASOFI3D solver failed"
fi

# Convert seismograms in SEG-Y format to the Madagascar RSF format.
convert_segy_to_rsf tmp/su/fullspace_vx.sgy
convert_segy_to_rsf ${TEST_PATH}/su/fullspace_vx.sgy

convert_segy_to_rsf tmp/su/fullspace_p.sgy
convert_segy_to_rsf ${TEST_PATH}/su/fullspace_p.sgy

# Read the files.
# Compare with the recorded output.
tests/compare_datasets.py tmp/su/fullspace_vx.rsf \
                          ${TEST_PATH}/su/fullspace_vx.rsf \
                          --rtol=1e-8 --atol=1e-10
result=$?
if [ "$result" -ne "0" ]; then
    error "Vx seismograms differ"
fi

tests/compare_datasets.py tmp/su/fullspace_p.rsf \
                          ${TEST_PATH}/su/fullspace_p.rsf \
                          --rtol=1e-1 --atol=1e-5
result=$?
if [ "$result" -ne "0" ]; then
    error "Pressure seismograms differ"
fi

log "PASS"
