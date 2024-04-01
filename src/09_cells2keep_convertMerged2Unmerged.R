## Amy Olex
## 5/8/2023
## Script to convert a list of barcodes to keep from the merged Loupe file to 
## individual txt files with a list of barcodes that can be re-imported by the merge script.

## INPUTS:
##         - Original CSV file used to do the original merging.
##         - the full list of cells to keep obtained from the merged Loupe browser
##         - An output directory.

library(optparse)
library(foreach)

option_list = list(
  make_option(c("-r", "--runid"), type="character", default=NULL, 
              help="Required. A unique name for this analysis.", metavar="character"),
  make_option(c("-c", "--configfile"), type="character", 
              help="Required. Input config file with sample information.", metavar="character"),
  make_option(c("-o", "--outdir"), type="character", 
              help="output directory for results report (must already exist). Default is current directory.",
              default = "./", metavar="character"),
  make_option(c("-k", "--keep"), type="character", 
              help="Optional. A TSV file with a list of barcodes that should be KEPT after merging (no header).  This was added as a specialty feature for merging very large datasets in batches. It needs to already be formatted with the numbered extensions for each sample.", 
              default = "", metavar="character")
); 

opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);

# test for NULL required arguments
if (is.null(opt$runid)){
  print_help(opt_parser)
  stop("A unique analysis name must be provided.", call.=FALSE)
}

if (is.null(opt$configfile)){
  print_help(opt_parser)
  stop("A configuration file with sample information must be provided.", call.=FALSE)
}

if (is.null(opt$keep)){
  print_help(opt_parser)
  stop("A file with the barcodes to keep must be provided.", call.=FALSE)
}


runID <- opt$runid
inFile <- opt$configfile
outDir <- opt$outdir
savedir <- paste0(outDir,runID)
keep <- opt$keep

print("Summary of input options:\n")
print(paste("Run Name: ", runID))
print(paste("ConfigFile: ", inFile))
print(paste("Output Directory:" , outDir))
print(paste("Save Directory:" , savedir))
print(paste("Keeping cells in file: ", keep))


# Read in the provided config file and loop for each row.
toProcess = read.table(inFile, header=TRUE, sep=",", stringsAsFactors = FALSE)
print(paste(dim(toProcess)[1], " rows were found."))

## Import the barcodes to keep
barcodes_to_keep = ""
tmp <- read.delim(keep, header=FALSE)
barcodes_to_keep <- as.data.frame(t(as.data.frame(strsplit(tmp$V1,split = "-"))))

for(i in 1:dim(toProcess)[1]){
  to_keep_this_sample <- barcodes_to_keep[barcodes_to_keep$V2 == i,,drop=FALSE]
  if(length(to_keep_this_sample$V1)>0){
    print(paste("Processing sample:", toProcess$SampleName[i]))
    to_keep_this_sample$barcode <- paste0(to_keep_this_sample$V1,"-1")  ##paste a "-1" on the end of each barcode for this sample.
    to_keep <- to_keep_this_sample[,"barcode",drop=FALSE]
    
    keep_file_name <- paste0(dirname(toProcess$Cells2Keep[i]), "/", toProcess$SampleName[i], "_", runID, ".csv")
    toProcess[i,"Cells2Keep"] <- keep_file_name
    print("Writing CSV...")
    write.csv(to_keep, file=keep_file_name, quote = FALSE, row.names = FALSE)
  }
  else{
    print(paste("No barcodes, excluding sample: ",toProcess$SampleName[i]))
    toProcess[i,"Cells2Keep"] <- "NO-CELLS-DELETE"
  }
}

print("Completed. Writing config file...")
write.csv(toProcess, file=paste0(savedir, "_", basename(inFile)), quote = FALSE, row.names = FALSE)

print("All DONE")



