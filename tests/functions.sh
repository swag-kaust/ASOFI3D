# Functions related to testing.

setup () {
    # Setup environment for the execution of a test.
    rm -rf tmp/
    mkdir tmp
    mkdir tmp/in_and_out
    mkdir tmp/log
    mkdir tmp/model
    mkdir tmp/snap
    mkdir tmp/sources
    mkdir tmp/su
}
