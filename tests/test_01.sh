#!/usr/bin/env bash
# Regression test 01.
# Uses simulation parameters and data recorded from the previous run,
# of the ASOFI3D code.
. tests/functions.sh

MODEL="src/model_elastic.c"
TEST_PATH="tests/fixtures/test_01"
TEST_ID="TEST_01"

# Setup function prepares environment for the test (creates directories).
setup

# Preserve old model.
mv $MODEL ${MODEL}.bak.${TEST_ID}

on_exit() {
    mv ${MODEL}.bak.${TEST_ID} $MODEL
}

# Execute function 'on_exit' when this script exits to avoid resource leak.
trap on_exit INT TERM EXIT

# Copy test model.
cp "${TEST_PATH}/src/model_elastic.c"       src/
cp "${TEST_PATH}/in_and_out/sofi3D.json"    tmp/in_and_out
cp "${TEST_PATH}/sources/source.dat"        tmp/sources/

compile_code

# Run code.
echo "${TEST_ID}: Running solver. Output is captured to tmp/ASOFI3D.log"
./run_ASOFI3D.sh 16 tmp/ > tmp/ASOFI3D.log &
task_id=$!
animate_progress $task_id "${TEST_ID}: Running solver"

code=$?
if [ "$code" -ne "0" ]; then
    echo ${TEST_ID}: FAIL Running ASOFI3D failed > /dev/stderr
    exit 1
fi

# Convert seismograms in SEG-Y format to the Madagascar RSF format.
convert_segy_to_rsf tmp/su/test_vx.sgy
convert_segy_to_rsf ${TEST_PATH}/su/test_vx.sgy

# Read the files.
# Compare with the old output.
tests/compare_datasets.py tmp/su/test_vx.rsf ${TEST_PATH}/su/test_vx.rsf \
                          --rtol=1e-12 --atol=1e-14
result=$?
if [ "$result" -ne "0" ]; then
    echo "${TEST_ID}: FAIL Velocity x-component seismograms differ" > /dev/stderr
    exit 1
fi

echo "${TEST_ID}: PASS"
