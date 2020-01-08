#!/bin/bash
#
# This script requires INSTALL_PATH (typically /opt/umi_qc),
# UMI_QC_PATH (typically /opt/umi_qc/source), and
# CONDA_DEFAULT_ENV (typically /opt/miniconda) to be set.
#
# A miniconda install must exist at $CONDA_DEFAULT_ENV
# and $CONDA_DEFAULT_ENV/bin must be in the PATH
#
# Otherwise, this only requires the existence of the following files:
#	requirements-conda.txt
#	requirements-other.txt

set -e -o pipefail

apt-get update
apt-get install -y -qq --no-install-recommends \
	build-essential cmake g++-5 gcc-5 libz-dev \
	bc

CONDA_CHANNEL_STRING="--override-channels -c bioconda -c conda-forge -c anaconda"

mkdir -p $INSTALL_PATH/umi_qc-etc
if [ ! -f $INSTALL_PATH/umi_qc-etc/umi_qc ]; then
	ln -s $UMI_QC_PATH $INSTALL_PATH/umi_qc-etc/umi_qc
fi
if [ ! -f $INSTALL_PATH/umi_qc-etc/conda-env ]; then
	ln -s $CONDA_DEFAULT_ENV $INSTALL_PATH/umi_qc-etc/conda-env
fi

# setup/install viral-ngs directory tree and conda dependencies
sync

conda install -y \
	-q $CONDA_CHANNEL_STRING \
	--file "$UMI_QC_PATH/requirements-conda.txt"

# clean up
conda clean -y --all

# Download piledriver from requirements-git.txt and compile binaries
cd /opt/
while read line; do
	git clone $line
done < $UMI_QC_PATH/requirements-git.txt

# Use sbt to compile .jar
cd $FGBIO_BASE
git checkout al_extract_umi_qs
sbt assembly

# Use gcc-5 to compile binaries
update-alternatives --install /usr/bin/gcc gcc /usr/bin/gcc-5 10
update-alternatives --install /usr/bin/gcc gcc /usr/bin/gcc-5 20
update-alternatives --install /usr/bin/g++ g++ /usr/bin/g++-5 10
update-alternatives --install /usr/bin/g++ g++ /usr/bin/g++-5 20
update-alternatives --install /usr/bin/cc cc /usr/bin/gcc 30
update-alternatives --set cc /usr/bin/gcc
update-alternatives --install /usr/bin/c++ c++ /usr/bin/g++ 30
update-alternatives --set c++ /usr/bin/g++
mkdir $PILEDRIVER_PATH/build
cd $PILEDRIVER_PATH/build
cmake ..
make