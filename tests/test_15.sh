#!/usr/bin/env bash
# Regression test 01.
# Uses simulation parameters and data recorded from the previous run,
# of the ASOFI3D code.
. tests/functions.sh

readonly MODEL="src/model_elastic.c"
readonly TEST_PATH="tests/fixtures/test_15"
readonly TEST_ID="TEST_15"

# Setup function prepares environment for the test (creates directories).
setup

# Copy test model.
cp "${TEST_PATH}/asofi3D.json"   tmp/in_and_out
cp "${TEST_PATH}/source.dat"     tmp/sources/

compile_code

run_solver np=16 dir=tmp log=ASOFI3D.log

# Read the files.
# Compare with the old output.
tests/compare_datasets.py \
    tmp/snap/test.bin.div ${TEST_PATH}/snap/test.bin.div \
    --rtol=1e-12 --atol=1e-14
result=$?
if [ "$result" -ne "0" ]; then
    error "Snapshots .div differ"
fi

log "PASS"
