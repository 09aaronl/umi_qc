FROM quay.io/broadinstitute/viral-baseimage:0.1.14

# Largely borrowed from https://github.com/broadinstitute/viral-ngs/blob/master/Dockerfile

LABEL maintainer "alin@broadinstitute.org"

# to build:
#   docker build . --build-arg <...>
#
# to run:
#   docker run --rm <image_ID> "<command>.py subcommand"
#
# to run interactively:
#   docker run --rm -it <image_ID>

ENV \
	INSTALL_PATH="/opt/umi_qc" \
	UMI_QC_PATH="/opt/umi_qc/source" \
	MINICONDA_PATH="/opt/miniconda" \
	PILEDRIVER_PATH="/opt/piledriver" \
	FGBIO_BASE="/opt/fgbio" \
	FGBIO_PATH="$FGBIO_BASE/target/scala-2.13"

ENV \
	PATH="$FGBIO_PATH:$PILEDRIVER_PATH/bin:$UMI_QC_PATH/scripts:$MINICONDA_PATH/bin:/usr/local/sbin:/usr/local/bin:/usr/sbin:/usr/bin:/sbin:/bin" \
	CONDA_DEFAULT_ENV=$MINICONDA_PATH \
	CONDA_PREFIX=$MINICONDA_PATH \
	JAVA_HOME=$MINICONDA_PATH

WORKDIR $INSTALL_PATH
COPY Dockerfile $UMI_QC_PATH/
COPY docker/install-*.sh $UMI_QC_PATH/docker/
COPY requirements-*.txt $UMI_QC_PATH/
COPY scripts/* $UMI_QC_PATH/scripts/
RUN $UMI_QC_PATH/docker/install-umi_qc.sh

RUN /bin/bash -c "set -e; echo -n 'multiqc version: '; multiqc --version"

CMD ["/bin/bash"]
