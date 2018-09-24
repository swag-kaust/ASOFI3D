#!/usr/bin/env bash
# Rename test.SOFI3D.v to test.v for all variables.
# Run it from the model subdirectory.

for i in *; do
      mv "$i" "${i/.SOFI3D/}"
done

