#!/bin/bash

# Copyright (c) 2018
# Amy L. Olex, Virginia Commonwealth University
# alolex at vcu.edu

## Usage Function
usage()
{
cat << EOF
usage: $0 options
This script is used to run the single cell sequencing analysis tool CellRanger from 10X for a list of samples.

REQUIRED ARGUMENTS:
   	-f  Tab delimited input file with list of samples to process and meta data.
	-o	Output directory to save results in.
	-F	Force the number of cells to the specified number in the input file
	-n

-f INPUT File:
A tab delimited file with the samples to process along with metadata information and the location of the input raw fastq files.  
Column information must be the following in this order: SampleName, SampleID, Reference, ExpectedCells, FastqDirectory
The Reference should be listed as one of the following: grch38, mm10, or hg19mm10

-o Directory where the results should be saved. Subfolders will be created using the SampleID given in the input file.

-F Force the number of cells to thoe specified in the input file.

-n Max number of local cores to request at one time.

EXAMPLE USAGE:
   >> ~/02_runCellRanger.sh -f /path/to/inventoryFile.list -o /path/to/output/dir/ -n 35 -F

EOF
}

## Parsing input arguments

FILE=
ODIR=
FORCE=
CORES=
while getopts “hf:o:n:F” OPTION
do
     case $OPTION in
         h)
             usage
             exit 1
             ;;
         f)
             FILE=$OPTARG
             ;;
	 o)
	     ODIR=$OPTARG
	     ;; 
	 n)
             CORES=$OPTARG
	     ;;
	 F)
	     FORCE=1
	     ;;
         ?)
             usage
             exit
             ;;
     esac
done

## Handle empty arguments

if [[ -z $FILE ]]; then usage; exit 1; fi
if [[ -z $ODIR ]]; then usage; exit 1; fi
if [[ -z $CORES ]]; then usage; exit 1; fi
if [[ -z $FORCE ]]; then FORCE=0; fi

cd $ODIR

INPUT=$FILE

cat $INPUT | while read sample sampleid ref cellc fastqdir
do 
	if [[ "$ref" == "mm10" ]]
	then
		reference="/data/refGenomes/CellRanger/refdata-cellranger-mm10-3.0.0"
	elif [[ "$ref" == "hg19mm10" ]]
	then
		reference="/data/refGenomes/CellRanger/refdata-cellranger-hg19-and-mm10-3.0.0"
	elif [[ "$ref" == "grch38" ]]
	then
		reference="/data/refGenomes/CellRanger/refdata-cellranger-GRCh38-3.0.0"	
	elif [[ "$ref" == "hg19" ]]
        then
                reference="/data/refGenomes/CellRanger/refdata-cellranger-hg19-3.0.0"
	elif [[ "$ref" == "mm10firefly" ]]
	then
		reference="/vcu_gpfs2/home/mccbnfolab/ref_genomes/mouse/mm10firefly"
	else
		echo "ERROR: Reference Genome Not Found"; usage; exit 1
	fi

	
	if [[ $FORCE == 1 ]]
	then
                echo Processing $sample with reference $reference and raw fastq files located in $fastqdir;
                echo cellranger count --id=$sampleid --transcriptome=$reference --fastqs=$fastqdir --sample=$sample --force-cells=$cellc --localcores=$CORES
	
                cellranger count --id=$sampleid --transcriptome=$reference --fastqs=$fastqdir --sample=$sample --force-cells=$cellc --localcores=$CORES
	else
		echo Processing $sample with reference $reference;
		echo cellranger count --id=$sampleid --transcriptome=$reference --fastqs=$fastqdir --sample=$sample --expect-cells=$cellc --localcores=$CORES
	
		cellranger count --id=$sampleid --transcriptome=$reference --fastqs=$fastqdir --sample=$sample --expect-cells=$cellc --localcores=$CORES
	fi
	
	
done
