#!/usr/bin/env bash
MODEL="src/hh_elastic.c"
TEST_PATH="tests/fixtures/test_01"

# Preserve old model.
mv $MODEL ${MODEL}.bak

# Copy test model.
cp "${TEST_PATH}/src/hh_elastic.c"       src/
cp "${TEST_PATH}/in_and_out/sofi3D.json" tmp/in_and_out
cp "${TEST_PATH}/sources/source.dat"     tmp/sources/

# Compile code.
cd src && make sofi3D > /dev/null && cd ..

# Run code.
./run_ASOFI3D.sh 16 tmp/ > /dev/null
code=$?

if [ "$code" -ne "0" ]; then
    echo TEST_01: FAIL > /dev/stderr
    exit 1
fi

# Convert seismograms in SEG-Y format to the Madagascar RSF format.
sfsegyread < tmp/su/test_vx.sgy --out=stdout \
    tfile=tmp/su/test_vx_trace.rsf > tmp/su/test_vx.rsf

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
