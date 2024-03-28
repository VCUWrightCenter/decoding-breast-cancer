## Step 07: Create Loupe File

Note: During the development of this manuscript the RLoupe package that converts Seurat objects to Loupe-compatible files had not been released and this was our work-around. 

Step | Details
---  | ---
Tool Name and Version | CellRanger cellranger-6.0.1
Pipeline script | N/A, calls the reanalyze function directly.
Input File Format | the .h5 file generated from the merging script in step 06.
Example Usage | cellranger reanalyze --id=UniquRunID --description=RunDescription --matrix=/path/to/merged_file.h5 --params=/path/to/reanalyze_params.txt
Output Files used in downstream steps? | Yes, Loupe file used for figure generation and analysis.
Output File formats and type | CellRanger “outs” directory with a Loupe File.

Once the Loupe file is generated, it is opened and the annotation CSV files from step 06 are manually imported.  This includes the UMAP and tSNE layouts and the SNN clusters.  Because this was a work-around and not a supported conversion method by 10X the layout and clustering of the default generated data in the Loupe file cannot be trusted as can be seen by the SampleID that ever only has 99 barcodes in it.  The LoupeR package fixes most of this; however, at this point the annotations are being imported incorrectly.
