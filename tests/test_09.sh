#!/usr/bin/env bash
# Regression test 09.
# Check option readmod=1 (reading model from binary files).
# The test is organized as follows:
# 1. Run simulation in which model is generated on-the-fly and write the model
#    on disk.
# 2. Run simulation in which model is read from files written in the previous
#    step.
# 3. Compare that the seismograms are close to each other.
. tests/functions.sh

readonly TEST_PATH="tests/fixtures/test_09"
readonly TEST_ID="TEST_09"

setup

# Copy test model.
cp "${TEST_PATH}/in_and_out/sofi3D-readmod=-1.json"  tmp/in_and_out/sofi3D.json
cp "${TEST_PATH}/sources/source.dat"                 tmp/sources/

compile_code

# Run code.
logfile="ASOFI3D-readmod=-1.log"
run_solver np=16 dir=tmp log="${logfile}"

# Copy second parameter file.
cp "${TEST_PATH}/in_and_out/sofi3D-readmod=1.json"  tmp/in_and_out/sofi3D.json

# Run code.
logfile="ASOFI3D-readmod=1.log"
run_solver np=16 dir=tmp log="${logfile}"

# Convert seismograms in SEG-Y format to the Madagascar RSF format.
convert_segy_to_rsf tmp/su/test-readmod=-1_vx.sgy
convert_segy_to_rsf tmp/su/test-readmod=1_vx.sgy

# Read the files.
# Compare with the recorded output.
tests/compare_datasets.py \
    tmp/su/test-readmod=-1_vx.rsf tmp/su/test-readmod=1_vx.rsf \
    --rtol=1e-15 --atol=1e-15
result=$?
if [ "$result" -ne "0" ]; then
    error "Traces differ"
fi

log "PASS"
