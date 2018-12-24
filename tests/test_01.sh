#!/usr/bin/env bash
# Regression test 01.
# Uses simulation parameters and data recorded from the previous run,
# of the ASOFI3D code.
. tests/functions.sh

readonly MODEL="src/model_elastic.c"
readonly TEST_PATH="tests/fixtures/test_01"
readonly TEST_ID="TEST_01"

# Setup function prepares environment for the test (creates directories).
setup

backup_default_model

# Copy test model.
cp "${TEST_PATH}/src/model_elastic.c"       src/
cp "${TEST_PATH}/in_and_out/sofi3D.json"    tmp/in_and_out
cp "${TEST_PATH}/sources/source.dat"        tmp/sources/

compile_code

run_solver np=16 dir=tmp log=ASOFI3D.log

# Convert seismograms in SEG-Y format to the Madagascar RSF format.
convert_segy_to_rsf tmp/su/test_vx.sgy
convert_segy_to_rsf ${TEST_PATH}/su/test_vx.sgy

# Read the files.
# Compare with the old output.
tests/compare_datasets.py tmp/su/test_vx.rsf ${TEST_PATH}/su/test_vx.rsf \
                          --rtol=1e-12 --atol=1e-14
result=$?
if [ "$result" -ne "0" ]; then
    error "Velocity x-component seismograms differ"
fi

log "PASS"
