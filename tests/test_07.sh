#!/usr/bin/env bash
# Regression test 07.
# Check the reciprocity.
# Precisely, we do the following:
# Run simulation 1 with source type "force in x-direction" and measure v_y
# in receiver;
# Run simulation 2 with source type "force in y-direction" and measure v_x
# in the receiver. Here the source and the receiver are swapped.
# Check that the seismograms from these simulations are the same.
# We consider isotropic half-space (two-layer) medium.
. tests/functions.sh

MODEL="src/model_elastic.c"
TEST_PATH="tests/fixtures/test_07"
TEST_ID="TEST_07"

setup

# Preserve old model.
mv $MODEL ${MODEL}.bak

# Copy test model.
cp "${TEST_PATH}/model_elastic.c"            src/
cp "${TEST_PATH}/sofi3D_force_in_x.json"     tmp/in_and_out/sofi3D.json
cp "${TEST_PATH}/source_force_in_x.dat"      tmp/sources/
cp "${TEST_PATH}/receiver_force_in_x.dat"    tmp/receiver/

compile_code

# Run code.
log "Running solver. Output is captured to tmp/ASOFI3D.log"
./run_ASOFI3D.sh 16 tmp/ > tmp/ASOFI3D.log &
task_id=$!
animate_progress $task_id "${TEST_ID}: Running solver"

code=$?
if [ "$code" -ne "0" ]; then
    error "Running ASOFI3D failed"
fi

cp "${TEST_PATH}/sofi3D_force_in_y.json"      tmp/in_and_out/sofi3D.json
cp "${TEST_PATH}/source_force_in_y.dat"       tmp/sources/
cp "${TEST_PATH}/receiver_force_in_y.dat"     tmp/receiver/

# Run code.
echo "${TEST_ID}: Running solver. Output is captured to tmp/ASOFI3D.log"
./run_ASOFI3D.sh 16 tmp/ > tmp/ASOFI3D.log &
task_id=$!
animate_progress $task_id "${TEST_ID}: Running solver"

code=$?
if [ "$code" -ne "0" ]; then
    error "Running ASOFI3D failed"
fi

# Convert seismograms in SEG-Y format to the Madagascar RSF format.
convert_segy_to_rsf tmp/su/force_in_x_vy.sgy
convert_segy_to_rsf tmp/su/force_in_y_vx.sgy

# Read the files.
# Compare with the recorded output.
tests/compare_datasets.py tmp/su/force_in_x_vy.rsf tmp/su/force_in_y_vx.rsf \
    --rtol=1e-15 --atol=1e-15
result=$?
if [ "$result" -ne "0" ]; then
    error "Traces differ"
fi

log "PASS"
