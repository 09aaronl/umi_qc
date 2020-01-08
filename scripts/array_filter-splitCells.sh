#!/bin/bash

if [[ "$#" -ne 5 ]]; then
	echo "Input: a bam file containing mapped NGS reads tagged with cell barcodes"
	echo "Output: one bam file per cell, containing NGS reads mapped to a specific chr or chr:region"
	echo "Usage: array_filter-splitCells.sh input.bam boolean_filter chr:region tag_cellBarcode output_base"
	echo "Example: array_filter-splitCells.sh sample1.array1.bam true KU182905.1 XC sample1.array1.EBOV"
	exit 1
fi

bam_raw=$1
boolean_filtered=$2
chr_region=$(basename $3)
chr=${chr_region%:*}
bam_header=${bam_raw%.*}.header.txt
bam_filtered=${bam_raw%.*}.$chr_region.bam

tag=$4
output_base=$5

samtools index $bam_raw

if [ "$boolean_filtered" = "true" ]; then
	chr_SQ=$(samtools view -H $bam_raw | grep SN:$chr | head -n 1)
	samtools view -H $bam_raw > $bam_header
	sed -i '/@SQ/d' $bam_header
	echo $chr_SQ | tr " " "\t" >> $bam_header
	cat $bam_header <(samtools view $bam_raw $chr_region) | samtools view -h -b -o $bam_filtered -
fi

bamtools split -tag $tag -in $bam_filtered -stub $output_base
