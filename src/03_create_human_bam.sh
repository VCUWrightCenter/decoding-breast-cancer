#!/bin/bash

# Copyright (c) 2019
# Amy L. Olex, Virginia Commonwealth University
# alolex at vcu.edu

## Usage Function
usage()
{
cat << EOF
usage: $0 options
This script is used on the Cell Ranger count output and must be run from the outs directory.  It extracts all the human barcodes from teh GEM analysis, filters the BAM file, sorts the BAM file and then exports to 2 Fastq files.

REQUIRED ARGUMENTS:
	-f	Input metadata file with the following columns for each sample: SampleName, PathToGEM, PathToBAM
OPTIONAL ARGUMENTS:
	-d		Directory of where the resulting fastq files should be saved.  Default is the current directory.

-f Sample list:
Tab delimited list of samples to process. Each row includes a unique sample name followed by the path to the GEM file for that sample and the the path to the BAM file for that sample.

-d RESULTS Directory:
Directory where results should be saved to.

EXAMPLE USAGE:
   >> ~/03_create_human_bam.sh -f samples.list  -d /path/to/results/dir/

EOF
}

## Parsing input arguments

FILE=
DIR= 
while getopts “hf:d:” OPTION
do
     case $OPTION in
         h)
             usage
             exit 1
             ;;
         f)
             FILE=$OPTARG
             ;;
         d)
             DIR=$OPTARG
             ;;
         ?)
             usage
             exit
             ;;
     esac
done

## Handle empty arguments

if [[ -z $FILE ]]; then usage; exit 1; fi
if [[ -z $DIR ]]; then DIR="."; fi

cat $FILE | while read sampleid gempath bampath
do
	echo Filtering Barcodes for $sampleid using GEM at $gempath

	tail -n +2 $gempath/gem_classification.csv | grep hg19 | cut -d ',' -f1 > $gempath/human_barcodes.txt 
	tail -n +2 $gempath/gem_classification.csv | grep mm10 | cut -d ',' -f1 > $gempath/mouse_barcodes.txt
	tail -n +2 $gempath/gem_classification.csv | grep Multiplet | cut -d ',' -f1 > $gempath/multiplet_barcodes.txt

	echo Subsetting BAM

	subset-bam --bam $bampath/possorted_genome_bam.bam --cell-barcodes $gempath/human_barcodes.txt --out-bam $bampath/human_genome_bam.bam --cores 20

	echo Converting BAM to Fastq

	bamtofastq $bampath/human_genome_bam.bam $DIR/$sampleid

	echo Renaming and Moving Files

	cd $DIR/$sampleid
	newdir=`ls`
	cd $newdir
	rename bamtofastq $sampleid bamtofastq*

	##mv * $DIR/
	echo Completed $sampleid
done
