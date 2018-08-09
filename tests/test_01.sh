#!/usr/bin/env bash
# Regression test 01.
# Uses simulation parameters and data recorded from the previous run,
# of the ASOFI3D code.
. tests/functions.sh

MODEL="src/hh_elastic.c"
TEST_PATH="tests/fixtures/test_01"

# Setup function prepares environment for the test (creates directories).
setup

# Preserve old model.
mv $MODEL ${MODEL}.bak

# Copy test model.
cp "${TEST_PATH}/src/hh_elastic.c"       src/
cp "${TEST_PATH}/in_and_out/sofi3D.json" tmp/in_and_out
cp "${TEST_PATH}/sources/source.dat"     tmp/sources/

# Compile code.
cd src
make sofi3D > /dev/null
if [ "$?" -ne "0" ]; then
    cd ..
    echo TEST_01: FAIL > /dev/stderr
    exit 1
fi
cd ..

# Run code.
echo "TEST_01: Running solver. Output is captured to tmp/ASOFI3D.log"
./run_ASOFI3D.sh 16 tmp/ > tmp/ASOFI3D.log &
task_id=$!
animate_progress $task_id "TEST_01: Running solver"

code=$?
if [ "$code" -ne "0" ]; then
    echo TEST_01: FAIL Running ASOFI3D failed > /dev/stderr
    exit 1
fi

echo $PATH | tr ":" "\n"

# Convert seismograms in SEG-Y format to the Madagascar RSF format.
sfsegyread < tmp/su/test_vx.sgy --out=stdout \
    tfile=tmp/su/test_vx_trace.rsf > tmp/su/test_vx.rsf

sfsegyread < ${TEST_PATH}/su/test_vx.sgy \
    tfile=${TEST_PATH}/su/test_vx_trace.rsf \
    > ${TEST_PATH}/su/test_vx.rsf

# Read the files.
# Compare with the old output.
tests/compare_datasets.py tmp/su/test_vx.rsf ${TEST_PATH}/su/test_vx.rsf
result=$?
if [ "$result" -ne "0" ]; then
    echo "TEST_01: FAIL Velocity x-component seismograms differ" > /dev/stderr
    exit 1
fi

# Teardown
git checkout -- ${MODEL}

echo "TEST_01: PASS"
