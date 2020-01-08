[![Docker Repository on Quay](https://quay.io/repository/aelin/umi_qc/status "Docker Repository on Quay")](https://quay.io/repository/aelin/umi_qc)

# umi_qc
This repo contains the Docker files needed to build the `quay.io/aelin/umi_qc` image. This image is based on the ubuntu
image `quay.io/broadinstitute/viral-baseimage:0.1.14` and contains various tools for manipulating and checking quality of reads and unique molecular indices (UMIs).

 - file manipulation: bedtools, samtools, bamtools
 - read trimming and merging: trimmomatic, flash, flash2
 - base and QS pileup: piledriver
 - umi-specific utilities: umi_tools, fgbio
 - QC: fastqc, multiqc
