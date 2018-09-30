FROM ubuntu:16.04

MAINTAINER Dmitry Kabanov (dmitry.kabanov@kaust.edu.sa) 

# USAGE:
# To build:
# $ docker build . -t swag-kaust/swag-image:<TAG>
# To upload to the Docker hub:
# $ docker push swag-kaust/swag-image:latest
# $ docker push swag-kaust/swag-image:<TAG>
#
# In all commands above, <TAG> should be replaced by information that lets
# identify image version.
# Currently, tag format is YYYY-MM.BB, where
# YYYY - four-digit year,
# MM - two-digit month with leading zeroes,
# BB - two-digit build numer with leading zeroes.
# For example: 2018-07.04 is the tag of the fourth version of the image
# built in July 2018.

# Install build tools (gcc, openmpi, make, and cmake).
RUN apt-get -qq update > /dev/null && \
    apt-get --yes --quiet install build-essential \
                                  cmake \
                                  libopenmpi-dev \
                                  wget \
                                  git \
                                  > /dev/null && \
    rm -rf /var/lib/apt/lists/*

# Install Python 2.7 and Madagascar.
RUN MINICONDA_URL=https://repo.continuum.io/miniconda && \
    MINICONDA_FILE=Miniconda2-latest-Linux-x86_64.sh && \
    wget --quiet "$MINICONDA_URL/$MINICONDA_FILE" --output-file=wget.log && \
    bash "$MINICONDA_FILE" -b -p /opt/miniconda && \
    rm "$MINICONDA_FILE" && \
    # Enable conda for the current user.
    . /opt/miniconda/etc/profile.d/conda.sh && \
    # Enable conda for all users.
    ln -s /opt/miniconda/etc/profile.d/conda.sh /etc/profile.d/conda.sh && \
    conda activate base && \
    conda update --quiet conda && \
    conda create --yes --quiet --name envpy27 python=2.7 numpy && \
    conda activate envpy27 && \
    conda install --yes --quiet --channel swag-kaust madagascar

RUN adduser --disabled-password --gecos '' user
RUN echo ". /opt/miniconda/etc/profile.d/conda.sh" >> /home/user/.bashrc
RUN echo "conda activate envpy27" >> /home/user/.bashrc

USER user
WORKDIR /home/user

CMD [ "/bin/bash" ]
