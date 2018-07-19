#!/usr/bin/env bash
# Regression test 02.
# Uses simulation parameters and data from the original SOFI3D code,
# from the benchmark `fullspace`.
MODEL="src/hh_elastic.c"
TEST_PATH="tests/fixtures/test_02"

# Preserve old model.
mv $MODEL ${MODEL}.bak

# Copy test model.
cp "${TEST_PATH}/model/fullspace.c"             src/hh_elastic.c
cp "${TEST_PATH}/in_and_out/fullspace.json"     tmp/in_and_out/
cp "${TEST_PATH}/sources/fullspace_sources.dat" tmp/sources/

# Compile code.
cd src && make sofi3D > /dev/null && cd ..

# Run code.
./run_ASOFI3D.sh 16 tmp/
code=$?

if [ "$code" -ne "0" ]; then
    echo TEST_01: FAIL > /dev/stderr
    exit 1
fi

# Convert seismograms in SEG-Y format to the Madagascar RSF format.
sfsegyread < tmp/su/fullspace_vx.sgy --out=stdout \
    tfile=tmp/su/fullspace_vx_trace.rsf > tmp/su/fullspace_vx.rsf

# Read the files.
# Compare with the old output.
tests/compare_datasets.py tmp/su/test_vx.rsf ${TEST_PATH}/su/test_vx.rsf
result=$?
if [ "$result" -eq "0" ]; then
    echo TEST_01: PASS
else
    echo TEST_01: FAIL
fi

# Teardown
mv $MODEL.bak ${MODEL}
