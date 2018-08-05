#!/usr/bin/env bash
function exists ()
{
    # Check that the program $1 exists.
    command -v "$1" >/dev/null 2>&1
}

function success ()
{
    echo -e " \033[0;32m+\033[0m  $*";
}

function error ()
{
    echo -e " \033[1;31m-\033[0m  $*";
}

success=true

if exists mpicc
then
    success "Program mpicc is found"
else
    success=false
    error "Program mpicc is not found"
fi

if exists python
then
    success "Program python is found"
else
    success=false
    error "Program python is not found"
fi

if exists sfsegyread
then
    success "Program sfsegyread is found"
else
    success=false
    error "Program sfsegyread is not found"
fi

python -c 'import rsf.api as rsf' >/dev/null 2>&1
if [ $? -eq 0 ]
then
    success "Importing rsf.api in python works"
else
    success=false
    error "Importing rsf.api in python does not work"
fi

if [ "$success" = true ]
then
    success "Requirements for tests are satisfied"
else
    error "Requirements for tests are not satisfied. Cannot proceed"
    exit 1
fi
