## Step 04: Realign to Human only Genome

Step | Details
---  | ---
Tool Name and Version | CellRanger cellranger-6.0.1
Pipeline Script | 02_runCellRanger.sh
Input File Format | See Step 02; must edit the number of cells to those that are human only instead of the default targeted number.  Also add the "-F" for force cells option to stop CellRanger from doing further filtering of cells.
Example Usage | ~/02_runCellRanger.sh -f /path/to/sampleFile.list -o /path/to/output/dir/ -n 35 -F
Output Files used in downstream steps? | Yes, Human only data: filtered_feature_bc_matrix data files are directly input into step 5, finding dead cells.
Output File formats and type | See Step 02
