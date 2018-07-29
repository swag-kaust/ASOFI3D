#!/usr/bin/env bash
# Regression test 04.
# Check that the seismograms are the same if we swap a source and receiver
# while the medium is anisotropic.
. tests/functions.sh

MODEL="src/hh_elastic.c"
TEST_PATH="tests/fixtures/test_04"
TEST_ID="TEST_04"

# Setup function prepares environment for the test (creates directories).
setup

# Preserve old model.
mv $MODEL ${MODEL}.bak

# Copy test model.
cp "${TEST_PATH}/hh_elastic.c"     src/hh_elastic.c
cp "${TEST_PATH}/sofi3D.json"      tmp/in_and_out/sofi3D.json
cp "${TEST_PATH}/source.dat"       tmp/sources/
cp "${TEST_PATH}/receiver.dat"     tmp/receiver/

# Compile code.
cd src
make sofi3D > /dev/null
if [ "$?" -ne "0" ]; then
    cd ..
    echo "${TEST_ID}: FAIL" > /dev/stderr
    exit 1
fi
cd ..

# Run code.
echo "${TEST_ID}: Running solver. Output is captured to tmp/ASOFI3D.log"
./run_ASOFI3D.sh 16 tmp/ > tmp/ASOFI3D.log &
task_id=$!
animate_progress $task_id "${TEST_ID}: Running solver"

wait $task_id
code=$?

if [ "$code" -ne "0" ]; then
    echo "${TEST_ID}: FAIL Running ASOFI3D failed" > /dev/stderr
    exit 1
fi

# Convert seismograms in SEG-Y format to the Madagascar RSF format.
sfsegyread tape=tmp/su/test_p.sgy.shot1 \
    tfile=tmp/su/test_p_trace.rsf.shot1 \
    > tmp/su/test_p.rsf.shot1
sfsegyread tape=tmp/su/test_p.sgy.shot2 \
    tfile=tmp/su/test_p_trace.rsf.shot2 \
    > tmp/su/test_p.rsf.shot2

# Extract traces.
# For source 1 we extract trace 2, while for source 2 we extract trace 1.
# Madagascar enumerates traces starting with 0, that's why min2 and max2
# have such values.
sfwindow < tmp/su/test_p.rsf.shot1 min2=1 > tmp/su/trace2.rsf
sfwindow < tmp/su/test_p.rsf.shot2 max2=0 > tmp/su/trace1.rsf

# Read the files.
# Compare with the recorded output.
tests/compare_datasets.py \
    tmp/su/trace1.rsf tmp/su/trace2.rsf \
    --rtol=1e-10 --atol=1e-8
result=$?
if [ "$result" -ne "0" ]; then
    echo "${TEST_ID}: Traces differ" > /dev/stderr
    exit 1
fi

# Teardown.
# Restore default model.
git checkout -- ${MODEL}

echo "${TEST_ID}: PASS"
