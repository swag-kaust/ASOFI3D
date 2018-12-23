#!/usr/bin/env bash
# Regression test 14.
# Checks that the code is self-adjoint via the dot-product test.
# Precisely, simulation can be treated as a linear operator acting
# on the source term.
# The source term is chosen such that it is random all over the grid
# at time step 1.
# The result of the operator action is the pressure field, which is recorded
# on time step 2.
# We run two of such simulations and then check that the dot products
# s_1 .dot. p_2   and   s_2 .dot. p_1
# are equal to each other, where s, p are source and pressure, respectively,
# while subscripts 1 and 2 are the numbers of corresponding simulations.
# See also http://sepwww.stanford.edu/sep/prof/pvi/conj/paper_html/node9.html
. tests/functions.sh

readonly TEST_PATH="tests/fixtures/test_14"
readonly TEST_ID="TEST_14"

setup

# Copy test data.
cp "${TEST_PATH}/sofi3D.json"    tmp/in_and_out/
cp "${TEST_PATH}/source.dat"     tmp/sources/

compile_code

cd tmp
rm -rf snap_1 snap_2
rm -rf source_field_1 source_field_2

mkdir -p snap
mkdir -p source_field
cd ..

# Run code.
for i in 1 2; do
    log "Running solver. Output is captured to tmp/ASOFI3D_$i.log"
    ./run_ASOFI3D.sh 16 tmp/ > tmp/ASOFI3D_$i.log &
    task_id=$!
    animate_progress $task_id "Running solver"

    code=$?
    if [ "$code" -ne "0" ]; then
        error "Running ASOFI3D failed"
    fi
    mv tmp/snap tmp/snap_$i
    mv tmp/source_field tmp/source_field_$i

    mkdir tmp/snap
    mkdir tmp/source_field
done

${TEST_PATH}/analyze.py --directory=tmp/ --rtol=1e-6 --atol=0
result=$?
if [ "$result" -ne "0" ]; then
    error "Dot product test failed"
fi

log "PASS"
