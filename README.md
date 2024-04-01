# Decoding breast cancer: Unraveling subtype and model differences through multi-model single-cell RNA sequencing data integration

This repository contains a description and in-house scripts and code used to process the single cell data for Altman et.al manuscript.

## Bioinformatics Pipeline Steps

### Step 01: FastQC and MultiFastQC

Step | Details
---  | ---
Tool Name and Version | FastQC v0.11.9
Pipeline Script | 01_fastqc.sh
Input File Format | Text file with one FastQ file name and absolut path per line.
Example Usage | ~/01_fastqc.sh -f /path/to/fastqFiles.list -d /path/to/output/dir/

Step | Details
---  | ---
Tool Name and Version | MultiQC, version 1.17
Pipeline Script | None, manually run MultiQC in same directory as FastQC output.
Output Files used in downstream steps? | No, FastQC and MultiQC reports are manually reviewed to inform on data quality.


### Step 02: Single Sample Alignment

Step | Details
---  | ---
Tool Name and Version | CellRanger, cellranger-6.0.1
Pipeline Script |  02_runCellRanger.sh
Input File Format | No header, TAB delimited, has the following 5 columns: sample ID used in fastq file name, result directory name, reference genome (grch38, mm10, or grch38mm10), number of cells to target, absolute path to raw fastq data directory.
Example Usage | ~/02_runCellRanger.sh -f /path/to/sampleFile.list -o /path/to/output/dir/ -n 35
Output Files used in downstream steps? | Yes: PDX Samples have a GEM file used to identify human cells in the step 3, removing mouse cells.  Pure Human/Mouse or Organoid Samples utilize the filtered_bc_matrix data files as input into step 5, finding dead cells.
Output File formats and type | CellRanger “outs” directory with GEM file and "filtered_feature_bc_matrix" directory of files.

### Step 03: Removal of Mouse Cells and Conversion Back to FastQ

Step | Details
---  | ---
Tool Name and Version | subset-bam 1.1.0
Tool Name and Version | bamtofastq v1.3.2
Pipeline Script | 03_create_human_bam.sh
Input File Format | TAB delimited, has the following 3 columns: sample ID (same as Step 02), absolute path to samples “outs/analysis” directory, absolute path to samples “outs” directory.
Example Usage | ~/03_create_human_bam.sh -f /path/to/sampleFile.list  -d /path/to/output/dir/
Output Files used in downstream steps? | Yes, FastQ files used in step 4, realignment to human genome.
Output File formats and type | BAM and FastQ files that only contain human reads.

### Step 04: Realign to Human only Genome

Step | Details
---  | ---
Tool Name and Version | CellRanger cellranger-6.0.1
Pipeline Script | 02_runCellRanger.sh
Input File Format | See Step 02; must edit the number of cells to those that are human only instead of the default targeted number.  Also add the "-F" for force cells option to stop CellRanger from doing further filtering of cells.
Example Usage | ~/02_runCellRanger.sh -f /path/to/sampleFile.list -o /path/to/output/dir/ -n 35 -F
Output Files used in downstream steps? | Yes, Human only data: filtered_feature_bc_matrix data files are directly input into step 5, finding dead cells.
Output File formats and type | See Step 02


### Step 05: Finding Dead Cells

Step | Details
---  | ---
Tool Name and Version | R 4.1.3
Tool Name and Version | hdf5-1.12.2
Pipeline Script | 05_DeadCellAnalysis.R
Input File Format | TAB delimited, has the following 2 columns: sample ID, absolute path to human-aligned “outs” directory.
Example Usage (use --help option for all parameters) | Rscript 05_DeadCellAnalysis.R -r UniqRunID -c sample_tabdelim_file.list -f /path/to/mitogenes/MitoCodingGenes13_human.txt -s human -o /path/to/results/directory
Output Files used in downstream steps? | Yes, the CSV file “cells2keeps” listing the good quality barcodes is used in step 6, sample merging.
Output File formats and type | Tab delim .txt report, .png images of violin plots before and after filtering, and .csv “cells2keep” files listing barcodes to keep.

### Step 06: Sample Merging

Step | Details
---  | ---
Tool Name and Version | R 4.1.3
Tool Name and Version | hdf5-1.12.2
Pipeline script | 06_SeuratMerge.R
Input File Format | Comma delimited, has header row (SampleName,DataType,SamplePath,Source,Condition,Sex,TumorType,Cells2Keep), has the following 8 columns (but more can be added and will be used as cell annotations): Sample ID, data type (10X or seurat), absolute path to sample’s “filtered_feature_bc_matrix” directory, sample source (e.g. PDX, MGT, etc), condition or treatment, sample sex (M or F), tumor type (e.g. HCI011, UCD52, etc), absolute path to CSV file listing barcodes of cells to keep in merge.
Example Usage (use --help option for all parameters) | Rscript 06_SeuratMerge.R -r UniqRunID -c /path/to/input.csv -i LogNormalize -t simple --parallel -n 32 --filter
Output Files used in downstream steps? | Yes, .h5 file is used by CellRanger Reanalyze to generate a Loupe file
Output File formats and type | .h5 file, RData file contining the merged Seurat object, multiple .CSV annotation files (one for each column in the input CSV file)

### Step 07: Create Loupe File

Note: During the development of this manuscript the RLoupe package that converts Seurat objects to Loupe-compatible files had not been released and this was our work-around. 

Step | Details
---  | ---
Tool Name and Version | CellRanger cellranger-6.0.1
Pipeline script | N/A, calls the reanalyze function directly.
Input File Format | the .h5 file generated from the merging script in step 06.
Example Usage | cellranger reanalyze --id=UniquRunID --description=RunDescription --matrix=/path/to/merged_file.h5 --params=/path/to/reanalyze_params.txt
Output Files used in downstream steps? | Yes, Loupe file used for figure generation and analysis.
Output File formats and type | CellRanger “outs” directory with a Loupe File.

Once the Loupe file is generated, it is opened and the annotation CSV files from step 06 are manually imported. This includes the UMAP and tSNE layouts and the SNN clusters. Because this was a work-around and not a supported conversion method by 10X the layout and clustering of the default generated data in the Loupe file cannot be trusted as can be seen by the SampleID that ever only has 99 barcodes in it. The LoupeR package fixes most of this; however, at this point the annotations are being imported incorrectly.

### Step 08: Identifying Malignant Cells using InferCNV

This step was only done once on a 33% subsampled dataset as using the entire dataset exceeded out compute resources.  
This script only processes all samples in the input file in batches of "k" at a time. The same 5 manually selected samples are also run with each batch and are used to calculate the Zscores. 

Usage NOTE: The MetaData CSV for inferCNV has to have 3 columns named exactly “Sample_ID,TissueType,inferCNV.normal”. The TissueType and inferCNV.normal columns are identical except the inferCNV column has "REFERENCE" in replacement of the tissue type for those samples that are to be used as the static number of reference samples to calculate the inferCNV scroe Z-score. 
Usage NOTE: The MetaData CSV for inferCNV has to use the EXACT name sample names that are used in step 06 for sample merging (i.e. the first column in the list CSV file).


Step | Details
---  | ---
Tool Names and Versions | R 4.1.3
Tool Name and Version | JAGS 4.3.1
Tool Name and Version | InferCNV 1.10.1
Pipeline Script | 08_inferCNV_id_malignent_cells.R
Input File Format | The .RData file from the merging step.
Example Usage | Rscript 08_inferCNV_id_malignent_cells.R -r UniqRunID -o /path/to/existing/output/directory/ -i /path/to/merged_file.RData -k 5 -f /path/to/sample/phenotype_MetaData.csv -g /path/to/reference/genome/genes.gtf
Output Files used in downstream steps? | Yes, the CSV files were loaded into the Loupe File and used to identify clusters of cells that are normal-like.
Output File formats and type | CSV files with scores, RDA file with inferCNV object

The InferCNV scores were overlaid on the Loupe plot to identify if the normal samples correlated with low z-scores. This was done by discretizing the scores to round to the closest integer and then placing everything greater than 6 in a single group ">=6".  Any non-normal cells that clustered with the normal reference groups were selected for removal to generate the Cancer-Only merge (step described below).

### Step 09: Cancer-Only Merge.

To generate the Cancer-Only datset, the 3 subsampled master merges that included all samples was utilized to select the normal clusters in the UMAP plot to remove those cells and any non-normal cells that clustered with them. The chosen cancer-only barcodes were downloaded into a single file and input into the script "cells2keep_convertMerged2Unmerged_050823.R", which divided them up into new "cells2keep" csv files and saved them in the appropriate directories.  ---I NEED TO ANNOTATE THIS PROCeSS AND DESCRIBE THE SCRIPT---

Step | Details
---  | ---
Tool Names and Versions | R 4.1.3
Pipeline script | 06_SeuratMerge.R
Input File Format | Comma delimited, has header row (SampleName,DataType,SamplePath,Source,Condition,Sex,TumorType,Cells2Keep), has the following 8 columns (but more can be added and will be used as cell annotations): Sample ID, data type (10X or seurat), absolute path to sample’s “filtered_feature_bc_matrix” directory, sample source (e.g. PDX, MGT, etc), condition or treatment, sample sex (M or F), tumor type (e.g. HCI011, UCD52, etc), absolute path to CSV file listing barcodes of cells to keep in merge.


## Funding

[ADD FUNDING HERE]


