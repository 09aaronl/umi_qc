#!/bin/bash

fctn=$1
if [[ $fctn = "split" ]]
then
	if [[ "$#" -eq 5 || "$#" -eq 7 ]]
	then
		input_fastq=$2
		length=$3
		output_fastq1=$4
		if [[ "$#" -eq 5 ]]
		then
			output_fastq2=$5
		else
			length2=$5
			output_fastq2=$6
			output_fastq3=$7
		fi
		count=0
		while read line; do
			if [[ $count -eq 0 ]]
			then
				header=$line
				echo $header >> $output_fastq1
				echo $header >> $output_fastq2
				if [[ "$#" -eq 7 ]]
				then
					echo $header >> $output_fastq3
				fi
				(( count++ ))
			elif [[ $count -eq 1 ]]
			then
				seq=$line
				echo ${seq:0:length} >> $output_fastq1
				seq2=$(echo $seq | tail -c +$(( length + 1 )))
				if [[ "$#" -eq 5 ]]
				then
					echo $seq2 >> $output_fastq2
				else
					echo ${seq2:0:$(( ${#seq2} - length2 ))} >> $output_fastq2
					echo $seq2 | tail -c +$(( ${#seq2} - length2 + 1 )) >> $output_fastq3
				fi
				(( count++ ))
			elif [[ $count -eq 2 ]]
			then
				echo "+" >> $output_fastq1
				echo "+" >> $output_fastq2
				if [[ "$#" -eq 7 ]]
				then
					echo "+" >> $output_fastq3
				fi
				(( count++ ))
			elif [[ $count -eq 3 ]]
			then
				qs=$line
				echo ${qs:0:length} >> $output_fastq1
				qs2=$(echo $qs | tail -c +$(( length + 1 )))
				if [[ "$#" -eq 5 ]]
				then
					echo $qs2 >> $output_fastq2
				else
					echo ${qs2:0:$(( ${#qs2} - length2 ))} >> $output_fastq2
					echo $qs2 | tail -c +$(( ${#qs2} - length2 + 1 )) >> $output_fastq3
				fi
				count=0
			fi
		done < $input_fastq
	else
		echo "Some options are missing"
		echo "Usage for splitting at one end: fastq_split-merge.sh split input.fastq length output1.fastq output2.fastq"
		echo "Usage for splitting at both ends: fastq_split-merge.sh split input.fastq length output1.fastq length2 output2.fastq output3.fastq"
		echo "Output fastq files will be in the same order as the input fastq file"
		echo
		exit 1
	fi				
elif [[ $fctn = "merge" ]]
then
	input_fastq1=$2
	input_fastq2=$3
	output_fastq=$4
	count=0
	warn=0

	if [[ -z $input_fastq1 || -z $input_fastq2 || -z $output_fastq ]]
	then
		echo "Some options are missing"
		echo "Usage: fastq_split-merge.sh merge input1.fastq input2.fastq output.fastq"
		echo "Input fastq files must be in the same read order. Oputput fastq file will have the same order"
		echo
		exit 1
	fi

	while read -r line1 && read -r line2 <&3; do
		if [[ $count -eq 0 ]]
		then
			if [[ ( $line1 != $line2 ) && ( $warn < 3 ) ]]
			then
				echo "WARNING! Header of $input_fastq1 does not match $input_fastq2!"
				echo "$input_fastq1: $line1"
				echo "$input_fastq2: $line2"
				echo
				(( warn++ ))
			fi
			header=$line1
			echo $header >> $output_fastq
			(( count++ ))
		elif [[ $count -eq 1 ]]
		then
			echo ${line1}${line2} >> $output_fastq
			(( count++ ))
		elif [[ $count -eq 2 ]]
		then
			echo "+" >> $output_fastq
			(( count++ ))
		elif [[ $count -eq 3 ]]
		then
			echo ${line1}${line2} >> $output_fastq
			(( count++ ))
			count=0
		fi
	done < $input_fastq1 3<$input_fastq2
else
	echo "fastq_split-merge.sh is a tool that can either (a) split reads into 2 fastq files (splitting at the 5' or 3' end) or 3 files (splitting at both ends), or (b) merge reads from 2 fastq files"
	echo "Available functions: split merge"
	echo
fi
