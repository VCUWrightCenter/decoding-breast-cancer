#!/bin/bash

# Copyright (c) 2018
# Amy L. Olex, Virginia Commonwealth University
# alolex at vcu.edu

## Usage Function
usage()
{
cat << EOF
usage: $0 options
This script is used to run FastQC on a list of FastQ files.

REQUIRED ARGUMENTS:
   	-f      Input file with list of FastQ files to process.

OPTIONAL ARGUMENTS:
	-d		Directory of where results should be saved.  Sub-folder 00_fastqc will be created here.

-f INPUT File:
Input file with list of FastQ files to process.

-d RESULTS Directory:
Directory where results should be saved to.

EXAMPLE USAGE:
   >> ~/00_fastqc.sh -f /path/to/inventoryFile.list -d /path/to/results/dir/

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
if [[ -z $DIR ]]; then usage; exit 1; fi


cd $DIR

INPUT=$FILE

mkdir -p 01_fastqc

OUTDIR=$DIR/01_fastqc

cat $INPUT | while read file
do 
	echo $file;
	fastqc -t 4 -o $OUTDIR --noextract $file;
done
