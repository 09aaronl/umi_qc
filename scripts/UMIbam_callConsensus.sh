#!/bin/bash

if [[ "$#" -ne 5 ]]; then
	echo "Input: bam file containing UMI PCR duplicates of an individual RNA molecule, reference seq"
	echo "Output: bam file containing a single read with consensus seq/QS of UMI PCR duplicates"
	echo "Usage: UMIbam_callConsensus.sh input.bam reference.fasta tag_cellBarcode tag_UMI output.bam"
	echo "Example: UMIbam_callConsensus.sh sample1.array1.EBOV.cell1.UMI1.bam KU182905.1.fa XC XM sample1.array1.EBOV.cell1.UMI1.consensus.bam"
	exit 1
fi

bam_raw=$1
bam_stranded=${bam_raw%.*}.stranded.bam
bam_region=${bam_raw%.*}.region.bam
bam_region_pileup=${bam_raw%.*}.region.pileup.txt
sam_consensus=${bam_raw%.*}.consensus.sam
ref=$2

tag_CB=$3
tag_UMI=$4
XC=$(samtools view $bam_raw | head -n 1 | grep -oh "\<$tag_CB:Z:\w*")
XM=$(samtools view $bam_raw | head -n 1 | grep -oh "\<$tag_UMI:Z:\w*")

bam_consensus=$5

# Determine RNA molecule orientation
reads_fwd=$(samtools view -F 16 $bam_raw | wc -l)
reads_rev=$(samtools view -f 16 $bam_raw | wc -l)

# Discard fwd/rev comp reads, whichever is fewer
if [ "$reads_fwd" -gt "$reads_rev" ]; then
	samtools view -h -b -F 16 -o $bam_stranded $bam_raw
	samtools index $bam_stranded
	echo "Picked reads_fwd, $reads_fwd fwd > $reads_rev rev"
elif [ "$reads_fwd" -lt "$reads_rev" ]; then
	samtools view -h -b -f 16 -o $bam_stranded $bam_raw
	samtools index $bam_stranded
	echo "Picked reads_rev, $reads_fwd fwd < $reads_rev rev"
else
	# Have some flag if calling RNA orientation is confusing
	echo "Ambiguous strand. $reads_fwd fwd = $reads_rev rev. Discarding UMI"
	exit 1
fi

# Determine general RNA region, pick region with most reads nearby (currently, hardcoded for 1500 bp window)
num_reads=$(samtools view $bam_stranded | wc -l)
starts=$(samtools view $bam_stranded | awk '{print $4}')
max_reads=1
max_read_index=0
first_start=$(echo $starts | awk '{print $1}')
last_added=$first_start
last_start=$first_start

for (( a=1; a<=num_reads; a++ )); do
	curr_reads=1
	a_start=$(echo $starts | awk '{print $'"$a"'}')
	((a_dist = $a_start + 1500))
	for (( b=a+1; b<=num_reads; b++ )); do
		b_start=$(echo $starts | awk '{print $'"$b"'}')
		if [ "$a_dist" -ge "$b_start" ]; then
			((curr_reads+=1))
			last_added=$(echo $starts | awk '{print $'"$b"'}')
		fi
	done
	if [ "$curr_reads" -gt "$max_reads" ]; then
		max_reads=$curr_reads
		max_read_index=$a
		first_start=$(echo $starts | awk '{print $'"$a"'}')
		last_start=$last_added
	fi
done

# Recalculate first_start and last_start and max_reads, since samtools view region is inclusive
samtools view -h -b -o $bam_region $bam_stranded KU182905.1:"$first_start"-"$(($last_start + 88))"
first_start=$(samtools view $bam_region | awk '{print $4}' | head -n 1)
last_start=$(samtools view $bam_region | awk '{print $4}' | tail -n 1)
max_reads=$(samtools view $bam_region | wc -l)
echo Best cluster was $max_reads, spanning $first_start to $last_start

# Call bamtools piledriver to generate pileup
bamtools piledriver -fasta $ref -in $bam_region > $bam_region_pileup

# Loop through pileup to call consensus base and qs
header=true
consensus_base=""
consensus_qs=""
consensus_CIGAR=""
CIGAR_last="M"
pos=$first_start
while read line; do
	if [ "$header" = "false" ]; then
		# Skip until pos catches up to line
		line_end=$(echo $line | awk '{print $3}')
		while [ "$pos" -ne "$line_end" ]; do
			if [ "$CIGAR_last" = "N" ]; then
				((CIGAR_count++))
			else
				consensus_CIGAR+=$CIGAR_count
				consensus_CIGAR+=$CIGAR_last
				CIGAR_last="N"
				CIGAR_count=1
			fi
			((pos++))
		done
		qs_other=0
		qs_adj=0
		phred33_adj=0
		bases_all=$(echo $line | awk '{print $8,$9,$10,$11,$14,$15,$16,$17}')
		num_A=$(echo $bases_all | awk '{print $1}')
		num_C=$(echo $bases_all | awk '{print $2}')
		num_G=$(echo $bases_all | awk '{print $3}')
		num_T=$(echo $bases_all | awk '{print $4}')
		((qs_total_A=$(echo $bases_all | awk '{print $5}')-7*$num_A))
		((qs_total_C=$(echo $bases_all | awk '{print $6}')-7*$num_C))
		((qs_total_G=$(echo $bases_all | awk '{print $7}')-7*$num_G))
		((qs_total_T=$(echo $bases_all | awk '{print $8}')-7*$num_T))
		if [ "$qs_total_A" -lt "0" ]; then
			qs_total_A=0
		fi
		if [ "$qs_total_C" -lt "0" ]; then
			qs_total_C=0
		fi
		if [ "$qs_total_G" -lt "0" ]; then
			qs_total_G=0
		fi
		if [ "$qs_total_T" -lt "0" ]; then
			qs_total_T=0
		fi

		qs_max=$(echo "$qs_total_A $qs_total_C $qs_total_G $qs_total_T" | tr " " "\n" | sort -n | tail -n 1)
		qs_2nd=$(echo "$qs_total_A $qs_total_C $qs_total_G $qs_total_T" | tr " " "\n" | sort -n | tail -n 2 | head -n 1)

		if [ "$qs_max" -le "0" ] || [ "$qs_max" -eq "$qs_2nd" ]; then
			consensus_base+="N"
			consensus_qs+="#"
		else
			if [ "$qs_total_A" -eq "$qs_max" ]; then
				consensus_base+="A"
				((qs_other = $qs_total_C + $qs_total_G + $qs_total_T))
				((qs_adj=$qs_total_A-$qs_other))
			elif [ "$qs_total_C" -eq "$qs_max" ]; then
				consensus_base+="C"
				((qs_other = $qs_total_A + $qs_total_G + $qs_total_T))
				((qs_adj=$qs_total_C-$qs_other))
			elif [ "$qs_total_G" -eq "$qs_max" ]; then
				consensus_base+="G"
				((qs_other = $qs_total_A + $qs_total_C + $qs_total_T))
				((qs_adj=$qs_total_G-$qs_other))
			elif [ "$qs_total_T" -eq "$qs_max" ]; then
				consensus_base+="T"
				((qs_other = $qs_total_A + $qs_total_C + $qs_total_G))
				((qs_adj=$qs_total_T-$qs_other))
			fi
			if [ "$qs_adj" -lt "2" ]; then
				qs_adj=2
			fi
			if [ "$qs_adj" -gt "60" ]; then
				qs_adj=60
			fi
			phred33_adj=$(printf "\x$(printf %x $((qs_adj + 33)))")
			consensus_qs+="$phred33_adj"
		fi

		if [ "$CIGAR_last" = "M" ]; then
			((CIGAR_count++))
		else
			consensus_CIGAR+=$CIGAR_count
			consensus_CIGAR+=$CIGAR_last
			CIGAR_last="M"
			CIGAR_count=1
		fi

		((pos++))
	else
		header=false
	fi

done < $bam_region_pileup

consensus_CIGAR+=$CIGAR_count
consensus_CIGAR+=$CIGAR_last
echo pos is $pos and first_start is $first_start, consensus_CIGAR is $consensus_CIGAR

echo "Called bases $consensus_base"
echo "Called quals $consensus_qs"

# Generate a bam file
ref_base=$(basename $ref)
bam_base=$(basename $bam_raw)
echo -e "${bam_base%.*}\t0\t${ref_base%.*}\t$first_start\t255\t$consensus_CIGAR\t*\t0\t0\t$consensus_base\t$consensus_qs\t$XC\t$XM" > $sam_consensus
cat <(samtools view -H $bam_raw) $sam_consensus | samtools view -h -b -o $bam_consensus -
