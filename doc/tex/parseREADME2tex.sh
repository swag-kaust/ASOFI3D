#!/bin/bash

# parses README.md to tex

grep -A10000 "## Obtaining" ../../README.md | \
pandoc -t latex -o README.tex

