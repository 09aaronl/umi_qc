#!/bin/bash

if [[ "$#" -ne 6 ]]; then
	echo "Input: a bam file containing mapped NGS reads originating from a single cell, tagged with an UMI transcript barcode. User also provides a list of UMIs to dedup"
	echo "Output: a bam file containing one consensus per UMI"
	echo "Usage: cell_splitUMIs-dedup.sh input.bam reference.fasta tag_cellBarcode tag_UMI tagValues.txt output.bam"
	echo "Example: cell_splitUMIs-dedup.sh sample1.array1.EBOV.cell1.bam KU182905.1.fa XC XM UMIsToDedup.txt sample1.array1.EBOV.cell1.mrg01.bam"
	exit 1
fi

bam_cell=$1
ref=$2
tag_CB=$3
tag_UMI=$4
tag_UMI_file=$5
bam_output=$6

while read line; do
	echo Starting tag $tag_UMI:$line at $(date | awk '{print $4}')
	bamtools filter -tag $tag_UMI:$line -in $bam_cell -out ${bam_cell%.*}.$line.bam
	echo Starting to dedup at $(date | awk '{print $4}')
	UMIbam_callConsensus.sh ${bam_cell%.*}.$line.bam $ref $tag_CB $tag_UMI ${bam_cell%.*}.$line.consensus.bam
	echo Finished dedup at $(date | awk '{print $4}')
done < $tag_UMI_file

samtools merge $bam_output ${bam_cell%.*}.*.consensus.bam
