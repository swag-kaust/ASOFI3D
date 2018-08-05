#!/usr/bin/env bash
# Install dependencies necessary for building and running ASOFI3D.

# Install build tools (gcc, openmpi, make, and cmake).
apt-get -qq update > /dev/null &&
apt-get --yes --quiet install build-essential cmake > /dev/null
apt-get --yes --quiet install libopenmpi-dev

apt-get --yes --quiet install wget > /dev/null

# Create a nonroot user.
adduser --disabled-password --gecos "" asofi3duser

# Install Python 2.7 and Madagascar.
MINICONDA_URL=https://repo.continuum.io/miniconda
MINICONDA_FILE=Miniconda2-latest-Linux-x86_64.sh

wget "$MINICONDA_URL/$MINICONDA_FILE" --output-file=wget.log;
bash "$MINICONDA_FILE" -b -p /opt/miniconda
export PATH=/opt/miniconda/bin:$PATH

conda create --yes --name env python=2.7
source activate env
conda install --quiet --yes --channel swag-kaust madagascar

# Checking the installed versions
python -V
gcc --version

NP=$(nproc)
echo ${NP}
