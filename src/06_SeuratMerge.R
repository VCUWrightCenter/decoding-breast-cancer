## Amy Olex
## 10/03/2022 (updated from version on 4/5/2021)
## File to remove identified dead cells and merge scRNASeq data from single cell 10X experiments
## Steps for Filtering:
## 1) Import cells to keep and store in a list in same order as list of experiments to merge.
## 2) subset each Seurat object
## 3) proceed to merging

## Steps for Merging:
## 1) create Seurat .h5Seurat files for each sample
## 2) Rename Barcodes
## 3) Normalize each data set using SCT or Log
## 4) Create a list of Seurat objects to merge
## 5) Select integration features and find anchors
## 6) Merge Data Sets
## 7) Write out 10X files
## 8) Create sample annotations
## 9) Run PCA, tSNE, and UMAP, and output feature annotations
## 10) Cluster cells and output feature annotations
## 11) Save Merged and annotated Seurat Object as an h5Seurat file.

## INPUTS
## Path to cells to keep file
## Path to unfiltered h5Seurat file OR path to 10X filtered_features_bc_matrix folder (or equivalent)
## Flag to delineate input types.
## 
## sample file format:
## SampleName, DataType (Seurat|10X), SamplePath, Cells2KeepPath, Condition

## UPDATED 6/4/21 to also provide a simple merge option
## BUG Fixes 9/2/21 to output the conditions correctly
## UPDATED 9/2/21 to run PCA, UMAP, tSNE, and Clustering using a subset of features
## UPDATED 10/14/21 to fix the feature subset issue and add in CellCycle Annotations.
## UPDATED 10/03/22 updated to include filtering out dead cells from a provided list 
##         as well as not requiring the exact bc matrix directory name.
## UPDATED 12/12/2022 to also save an RDS file to save the Seurat object.
## UPDATED 3/25/2023 to allow the exclusion of specific cells BEFORE running any dead cell analyses.
## UPDATED 4/6/2023 to allow user to specify if they want the scaled RNA counts or the raw RNA counts to be output.  Previously is was only exporting the scaled counts.
## UPDATED 8/2/2023 to allow user option to do ambient RNA adjustments of read counts BEFORE any other processing is done on the data.
## UPDATED 1/23/2024 to utilize the new LoupeR package to directly generate a Loupe file.


library(Seurat)
#library(SeuratDisk)
library(stringr)
library(DropletUtils)
library(optparse)
library(foreach)
library(doParallel)
library(future)
library("biomaRt")
library(dplyr)
library("SoupX")
library("loupeR")
library("hdf5r")


### DEFINE FUNCTIONS
## The following file is saved on Fenn in case this URL goes dark: /vcu_gpfs2/home/alolex/src/WCCTR_RNASeq_Pipeline/SingleCell
mouse_human_genes = read.csv("http://www.informatics.jax.org/downloads/reports/HOM_MouseHumanSequence.rpt",sep="\t")
convert_human_to_mouse <- function(gene_list){
  ## This function was copied from https://support.bioconductor.org/p/129636/
  
  output = c()
  
  for(gene in gene_list){
    class_key = (mouse_human_genes %>% filter(Symbol == gene & Common.Organism.Name=="human"))[['DB.Class.Key']]
    if(!identical(class_key, integer(0)) ){
      mouse_genes = (mouse_human_genes %>% filter(DB.Class.Key == class_key & Common.Organism.Name=="mouse, laboratory"))[,"Symbol"]
      for(mouse_gene in mouse_genes){
        output = append(output,mouse_gene)
      }
    }
  }
  
  return (output)
}



## define functions
get_soup_groups <- function (sobj){
  sobj <- NormalizeData(sobj, verbose = FALSE)
  sobj <- FindVariableFeatures(sobj, selection.method = "vst", nfeatures = 2000, verbose=FALSE)
  sobj <- ScaleData(sobj, verbose=FALSE)
  sobj <- RunPCA(sobj, npcs = 20, verbose = FALSE)
  sobj <- FindNeighbors(sobj, dims=1:20) #1:20
  sobj <- FindClusters(sobj, resolution = 0.4, verbose=FALSE)
  return(sobj@meta.data[['seurat_clusters']])
}

add_soup_groups <- function (sobj){
  sobj$soup_group <- get_soup_groups(sobj)
  return(sobj)
}




localtest = FALSE
###########################################
#### Local Testing Block
if(localtest){
  setwd("/Users/alolex/Desktop/Harrell_LocalTesting")
  runID <- "TEST"
  inFile <- "./configtest2.csv"
  outDir <- "./"
  features <- ""
  savedir <- paste0(outDir,runID)
  numCores <- 2
  numAnchors <- 2000
  normalization <- "LogNormalize"
  mergeType <- "integration"
  parallel <- FALSE
  saveH5 <- TRUE
  species <- "mouse"
  regressCC <- FALSE
  exclude <- ""
  downsample <- 100
  filtercells <- TRUE
  exportCounts <- TRUE
  keep <- ""
  options(future.globals.maxSize = 3000 * 1024^2)
}
###########################################

#numCores <- detectCores()
#registerDoParallel(numCores)

options(future.globals.maxSize = 100000 * 1024^2)

start_time <- Sys.time()


option_list = list(
  make_option(c("-r", "--runid"), type="character", default=NULL, 
              help="Required. A unique name for this analysis.", metavar="character"),
  make_option(c("-c", "--configfile"), type="character", 
              help="Required. Input config file with sample information.", metavar="character"),
  make_option(c("-o", "--outdir"), type="character", 
              help="output directory for results report (must already exist). Default is current directory.",
              default = "./", metavar="character"),
  make_option(c("-i", "--normalization"), type="character", 
              help="Type of normalization to run prior to merging/integrating (SCT or LogNormalize). Default is SCT.",
              default = "SCT", metavar="character"),
  make_option(c("-t", "--type"), type="character", 
              help="Type of merging to do (simple or integration). Default is simple.",
              default = "simple", metavar="character"),
  make_option(c("-f", "--features"), type="character", 
              help="Optional. A list of features to use for the PCA, UMAP, tSNE, and Clustering.", 
              default = "", metavar="character"),
  make_option(c("-s", "--downsample"), type="numeric", 
              help="Optional. The percentage of cells to KEEP from each sample entered as an integer (eg. 20 for 20%). Default is 100.", 
              default = 100, metavar="character"),
  make_option(c("-e", "--exclude"), type="character", 
              help="Optional. A TSV file with a list of barcodes that should be removed FIRST before downsampling is performed (no header).  This was added as a specialty feature for merging very large datasets in batches. It needs to already be formatted with the numbered extensions for each sample.", 
              default = "", metavar="character"),
  make_option(c("-k", "--keep"), type="character", 
              help="Optional. A TSV file with a list of barcodes that should be KEPT after merging (no header).  This was added as a specialty feature for merging very large datasets in batches. It needs to already be formatted with the numbered extensions for each sample.", 
              default = "", metavar="character"),
  make_option(c("--parallel"), type="logical", 
              help="Optional. Use paralellization for Merging and Integration (note, normalization of individual samples is always parallelized).",
              default = FALSE, action = "store_true", metavar="logical"),
  make_option(c("--filter"), type="logical", 
              help="Optional. Filter to only selected cells for each file before merging. Default is False. If True, config file must include path to and name of file listing the barcodes of the cells to KEEP in the merged data set in a column named 'Cells2Keep'.",
              default = FALSE, action = "store_true", metavar="logical"),
  make_option(c("-n", "--numCores"), type="numeric", 
              help="Optional. Specify number of cores to use. Default is 1.",
              default = 1, metavar="character"),
  make_option(c("-a", "--numAnchors"), type="numeric", 
              help="Optional. Specify number of anchor genes to use. Default is 2000.",
              default = 2000, metavar="character"),
  make_option(c("-p", "--species"), type="character", 
              help="species to use for cell cycle scroing. Options are human or mouse (default = human).",
              default = "human", metavar="character"),
  make_option(c("--saveH5"), type="logical", 
              help="Optional. Save merged and annotated Seruat object as a .h5Seruat file. Default saves as a .RData file.",
              default = FALSE, action = "store_true", metavar="logical"),
  make_option(c("--regressCellCycle"), type="logical", 
              help="Optional. Regress out the cell cycle difference between S and G2M scores during normalization (Default = FALSE).",
              default = FALSE, action = "store_true", metavar="logical"),
  make_option(c("--ambientRNAadjust"), type="logical", 
              help="Optional. Uses SoupX to adjust for ambient RNA contamination. MUST have a column in input CSV named 'raw10Xdata' that lists the full path to the raw_feature_bc_matrix file in 10X format. (Default = FALSE).",
              default = FALSE, action = "store_true", metavar="logical")
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


runID <- opt$runid
inFile <- opt$configfile
outDir <- opt$outdir
savedir <- paste0(outDir,runID)
numCores <- opt$numCores
numAnchors <- opt$numAnchors
normalization <- opt$normalization
parallel <- opt$parallel
saveH5 <- opt$saveH5
mergeType <- opt$type
features <- opt$features
filtercells <- opt$filter
regressCC <- opt$regressCellCycle
downsample <- opt$downsample
exclude <- opt$exclude
species <- opt$species
exportCounts <- opt$exportCounts
keep <- opt$keep
ambientRNAadjust <- opt$ambientRNAadjust


if(species == "mouse"){
  print("Converting human Cell Cycle Genes to MOUSE...")
  s.genes <- convert_human_to_mouse(cc.genes.updated.2019$s.genes)
  g2m.genes <- convert_human_to_mouse(cc.genes.updated.2019$g2m.genes)
  cc.genes <- union(s.genes, g2m.genes)
}else{
  print("Using HUMAN Cell Cycle Genes...")
  s.genes <- cc.genes.updated.2019$s.genes
  g2m.genes <- cc.genes.updated.2019$g2m.genes
  cc.genes <- union(s.genes, g2m.genes)
}

#### Ok, print out all the current running options as a summary.

print("Summary of input options:\n")
print(paste("Run Name: ", runID))
print(paste("ConfigFile: ", inFile))
print(paste("Merge Type: ", mergeType))
print(paste("Output Directory:" , outDir))
print(paste("Cores Specified: ", numCores))
print(paste("Number of Anchor Genes: ", numAnchors))
print(paste("Integration/Normalization Type: ", normalization))
print(paste("Filtering Cells: ", filtercells))
print(paste("Using Parallel: ", parallel))
print(paste("Saving h5Seurat file: ", saveH5))
print(paste("Regressing out S-G2M cell cycle score: ", regressCC))
print(paste("Percent of cells to keep in downsampling: ", downsample))
print(paste("Excluding cells in file: ", exclude))


if(!parallel){
  print("WARNING: Seurat merge and integration functions will NOT be parallized.  Use the --parallel flag to parallelize these functions.")
}

############## PROCESS CONFIG FILE###########################
# Read in the provided config file and loop for each row.
toProcess = read.table(inFile, header=TRUE, sep=",", stringsAsFactors = FALSE)
#toProcess = read.table("/Users/alolex/Desktop/CCTR_Git_Repos/PBos_scRNASeq/MergeScriptTesting/configtest.csv", header=TRUE, sep=",", stringsAsFactors = FALSE)

print(paste(names(toProcess)))

toProcess$SeuratObj <- ""

print(paste(names(toProcess)))

print(paste(dim(toProcess)[1], " rows were found."))

## Import the barcodes to exclude if it exists
barcodes_to_remove = ""
if(exclude != ""){
  tmp <- read.delim(exclude, header=FALSE)
  barcodes_to_remove <- as.data.frame(t(as.data.frame(strsplit(tmp$V1,split = "-"))))
  
  #sub <- df[df$V2 == 21, "V1"]
}

## Import the barcodes to keep if it exists
barcodes_to_keep = ""
if(keep != ""){
  tmp <- read.delim(keep, header=FALSE)
  barcodes_to_keep <- as.data.frame(t(as.data.frame(strsplit(tmp$V1,split = "-"))))
  print("Keeping:")
  print(dim(barcodes_to_keep))
  #sub <- df[df$V2 == 21, "V1"]
}

#Register num cores for doParallel loop
registerDoParallel(numCores)

sc_start <- Sys.time()

seurat_list <- foreach(i=1:dim(toProcess)[1]) %dopar% {
  print(paste("Importing data for row", i, "from sample", toProcess[i,"SampleName"]))
  
  if(toProcess[i,"DataType"] == "Seurat"){
    print(paste0(toProcess[i,"SampleName"], ": Loading Seurat h5 File..."))
    library(SeuratDisk)
    if(file.exists(trimws(toProcess[i,"SamplePath"]))){ 
      h5 <- LoadH5Seurat(trimws(toProcess[i,"SamplePath"])) 
    }
    else{ 
      print(paste("ERROR, file not found: ",toProcess[i,"SampleName"]) ) 
    }
  }
  else if(toProcess[i,"DataType"] == "10X"){
    print(paste0(toProcess[i,"SampleName"], ": Loading 10X feature matrix..."))
    # Load the data set and create the Seurat object
    sc.data <- Read10X(data.dir = toProcess[i,"SamplePath"])  ## must point to the filtered_feature_bc_matrix directory.
    # Initialize the Seurat object
    h5 <- CreateSeuratObject(counts = sc.data, project = toProcess[i,"SampleName"], min.cells = 0, min.features = 0)
    
  }
  else{
    print("Error unknown data type.")
    quit(1)
  }
  
  # Run ambient RNA adjustment FIRST, before excluding or doing any further processing
  ###I need to finish this up to point to the correct column of the toProcess object.
  if(ambientRNAadjust){
    print("Ambient RNA Adjusment...")
    
    ## add Soup Groups to filtered feature data
    h5 <- add_soup_groups(h5)
    
    ##load in unfiltered raw data (must be in 10X format)
    raw_10x <- Read10X(data.dir = toProcess[i,"raw10Xdata"])
    rownames(raw_10x) <- gsub("_", "-", rownames(raw_10x))
    
    sc1 <- SoupChannel(raw_10x, h5@assays$RNA@counts)
    sc1 <- setClusters(sc1, h5$soup_group)
    
    png(file = paste0(savedir, toProcess[i,"SampleName"], "_SoupXestimates.png"), width = 1000, height = 500, res = 100)
      sc1 <- autoEstCont(sc1, doPlot=TRUE)
    dev.off()
    
    out1 <- adjustCounts(sc1, roundToInt = TRUE)
    h5[['originalcounts']] <- CreateAssayObject(counts = h5@assays$RNA@counts)
    h5@assays$RNA@counts <- out1
    #recalculate nCounts, nFeature, and create the UMAP plot
    h5$nCount_RNA = colSums(x = h5, slot = "counts")
    h5$nFeature_RNA = colSums(x = GetAssayData(object = h5, slot = "counts") > 0)
    
    print(paste0("Original Total Read Count: ", sum(h5@assays$originalcounts@counts)))
    print(paste0("Adjusted Total Read Count: ", sum(h5@assays$RNA@counts)))
    #reads_before_ambiant <- sum(h5@assays$originalcounts@counts)
    #reads_after_ambiant <- sum(h5@assays$RNA@counts)
    
  }
  
  if(filtercells){
    print(paste0(toProcess[i,"SampleName"], ": Loading list of cell barcodes to keep..."))
    file_name <- trimws(toProcess[i,"Cells2Keep"])
    if(file.exists(file_name)){
      cells2keep <- read.delim(file=file_name, header=TRUE) ##not sure if there is a header, check the file.
    }
    else {
      print(paste("ERROR, file not found: ", toProcess[i,"Cells2Keep"]))
    }
    
    print(paste0(toProcess[i,"SampleName"], ": Keeping ", length(cells2keep$barcode), " cells."))
    print(paste0(toProcess[i,"SampleName"], ": Downsampling to ", downsample, " percent."))
    
    if(downsample < 100){
      print(paste0(toProcess[i,"SampleName"], ": Downsampling kept cells to ", downsample, "%..."))
    }
    
    if(barcodes_to_remove != ""){
      to_remove_this_sample <- paste0(barcodes_to_remove[barcodes_to_remove$V2 == i, "V1"],"-1")
      print(paste0(toProcess[i,"SampleName"], ": Removing ",length(to_remove_this_sample)," cells to exclude from cells2keep."))
      cells2keep <- cells2keep[!(cells2keep$barcode %in% to_remove_this_sample),,drop=FALSE]
      print(paste0(toProcess[i,"SampleName"], ": After EXCLUSION keeping ", length(cells2keep$barcode), " cells."))
    }
    
    if(barcodes_to_keep != ""){
      to_keep_this_sample <- paste0(barcodes_to_keep[barcodes_to_keep$V2 == i, "V1"],"-1")
      print(paste0(toProcess[i,"SampleName"], ": Keeping ",length(to_keep_this_sample)," cells from cells2keep."))
      cells2keep <- cells2keep[cells2keep$barcode %in% to_keep_this_sample,,drop=FALSE]
      print(paste0(toProcess[i,"SampleName"], ": After FILTERING keeping ", length(cells2keep$barcode), " cells."))
      #print(head(cells2keep))
    }
    
    if(length(cells2keep$barcode) > 0){
      samplesize <- floor((downsample/100)*length(cells2keep$barcode))
      set.seed(i)
      sampledcells <- sample(x = cells2keep$barcode, size = samplesize, replace = F)
    
      print(paste0(toProcess[i,"SampleName"], ": After DOWNSAMPLING keeping ", length(sampledcells), " cells."))
    
    #print(paste0("Filtering cells2keep using file: ", toProcess[i,"Cells2Keep"]))
      #print(head(sampledcells))
      h5 <- subset(h5, cells = sampledcells)
    }
    else{
      print(paste0(toProcess[i,"SampleName"], ": Skipping iteration.."))
      next 
    }
  }
  
  #if(exclude == "" and !filtercells){
  #  to_remove_this_sample <- barcodes_to_remove[barcodes_to_remove$V2 == i, "V1"]
  #  print(paste("Excluding cells : ", str(length(to_remove_this_sample)) ))
  #  
  #  if(downsample < 100){
  #    print("Downsampling...")
  #  }
  #  
  #  samplesize <- floor((downsample/100)*length(cells2keep$barcode))
  #  sampledcells <- sample(x = cells2keep$barcode, size = samplesize, replace = F)
  #  
  #  h5 <- subset(h5, cells = sampledcells)
  #}
  
  print(paste0(toProcess[i,"SampleName"], ": Renaming Cells..."))
  ## Rename the barcodes in each file with a count greater than 1
  if(i > 1){
    h5 <- RenameCells(h5,  new.names = str_replace(names(h5$orig.ident), "-1", paste0("-",i)))
    #print(head(names(h5$orig.ident)))
  }
  
  mid_time <- Sys.time()
  
  if(normalization == "SCT"){
    print(paste0(toProcess[i,"SampleName"], ": SCTransform..."))
    ## Normalize each dataset with SCT
    h5 <- SCTransform(h5, vst.flavor = "v2", verbose = FALSE)
  }
  else if(normalization == "LogNormalize"){
    print(paste0(toProcess[i,"SampleName"], ": LogNormalize..."))
    h5 <- NormalizeData(h5, normalization.method = "LogNormalize", verbose = FALSE)
    #print(paste0(toProcess[i,"SampleName"], ": Fariable Features"))
    h5 <- FindVariableFeatures(h5, selection.method = "vst", nfeatures = numAnchors, verbose = FALSE)
    #print(paste0(toProcess[i,"SampleName"], ": scaleData"))
    h5 <- ScaleData(h5, verbose = FALSE)
    print(h5)
    
  }
  
  
  ## Add seurat obj to list now
  #seurat_list <- c(seurat_list,h5)  
  
  print(Sys.time() - mid_time)
  
  print(h5)
}  ## end processing of each sample.

### Remove the empty string elements before moving on in order to remove any samples that were completly removed from the data set.
#print("Removing empty elements...")
print(paste0("Seurat list length: ", length(seurat_list)))
#seurat_list <- unlist(lapply(seurat_list, function(z){ z[!is.na(z)]}))
#print(paste0("Seurat list length AFTER: ", length(seurat_list)))


print("Total Time to run sample normalization:")
print(Sys.time() - sc_start)

print("Sample import and normalization completed...")
print("Finding integration features...")
## Find integration features

#seurat_list
if(parallel){
  ## Set Future Plan for asynchronous execution
  plan("multicore", workers = numCores)
}

plan()


###### Run Integration
if(mergeType == "integration"){
  print("Running Integration...")
  ## Select integration features
  if(normalization == "SCT"){
    mid_time <- Sys.time()
    hfile_features <- SelectIntegrationFeatures(object.list = seurat_list, nfeatures = numAnchors)
    
    seurat_list <- PrepSCTIntegration(object.list = seurat_list, anchor.features = hfile_features, verbose = FALSE)
    print(Sys.time() - mid_time)
  } else{
    hfile_features <- numAnchors
  }
  ## Find the anchors and then integrate the data sets.
  mid_time <- Sys.time()
  print("Finding Anchors...")
  #seurat.anchors <- FindIntegrationAnchors(object.list = seurat_list, normalization.method = normalization, anchor.features = hfile_features, reference = 5, verbose = FALSE)
  seurat.anchors <- FindIntegrationAnchors(object.list = seurat_list, normalization.method = normalization, anchor.features = hfile_features, verbose = FALSE)
  

  if(regressCC){
    ## Ensure cell cycle genes are included in the list
    genes.to.integrate <- union(seurat.anchors@anchor.features, cc.genes)
  }
  else{
    genes.to.integrate <- seurat.anchors@anchor.features
  }
  
  ## Save Anchor Feature File
  write.csv(genes.to.integrate, file=paste0(savedir, "_Seurat", normalization, "Merge_AnchorList.csv"), quote = FALSE, row.names = FALSE)
  
  
  
  print(Sys.time() - mid_time)
  
  mid_time <- Sys.time()
  print("Performing integration...")
  #seurat.merged <- IntegrateData(anchorset = seurat.anchors, features.to.integrate = genes.to.integrate, normalization.method = normalization, verbose = FALSE)
  seurat.merged <- IntegrateEmbeddings(anchorset = seurat.anchors, verbose = TRUE, reduction = "pcaproject")
  print(Sys.time() - mid_time)
  
} else if(mergeType == "simple"){
  print("Running simple merge...")
  
  seurat.merged <- merge(seurat_list[[1]], y = seurat_list[2:length(seurat_list)], merge.data=TRUE, project=runID)
  seurat.merged <- FindVariableFeatures(seurat.merged, selection.method = "vst", nfeatures = numAnchors, verbose = FALSE)
}
#### End Integration


####
## Cell Cycle Scoreing
#####

seurat.merged <- ScaleData(seurat.merged, features = rownames(seurat.merged))
seurat.merged <- CellCycleScoring(seurat.merged, g2m.features = g2m.genes, s.features = s.genes, set.ident = TRUE)

cc <- data.frame("barcode" = names(seurat.merged$Phase), "CellCyclePhase" = seurat.merged$Phase)
write.csv(cc, paste0(savedir, "_CellCyclePhase_",mergeType,"MergedData.csv"), quote=FALSE, row.names = FALSE)
#DimPlot(seurat.merged, group.by = "Phase", reduction = "umap", label=FALSE)

if(regressCC){
  print("Regressing out CC...")
  seurat.merged$CC.Difference <- seurat.merged$S.Score - seurat.merged$G2M.Score
  seurat.merged <- ScaleData(seurat.merged, vars.to.regress = c("S.Score", "G2M.Score"), features = rownames(seurat.merged))
}

mid_time <- Sys.time()
print("Saving to 10X...")
## Saving to 10X format:
if(normalization == "SCT"){
  write10xCounts(x=seurat.merged@assays$SCT@counts, path=paste0(savedir, "_seurat_",mergeType,"Merge_",normalization,"_SCTdata_rawCounts.h5"), version="3")
  write10xCounts(x=seurat.merged@assays$RNA@counts, path=paste0(savedir, "_seurat_",mergeType,"Merge_",normalization,"_RNAdata_rawCounts.h5"), version="3")
  write10xCounts(x=seurat.merged@assays$SCT@data, path=paste0(savedir, "_seurat_",mergeType,"Merge_",normalization,"_SCTdata_scaledCounts.h5"), version="3")
  write10xCounts(x=seurat.merged@assays$RNA@data, path=paste0(savedir, "_seurat_",mergeType,"Merge_",normalization,"_RNAdata_scaledCounts.h5"), version="3")

  print(Sys.time() - mid_time)
} else if(normalization == "LogNormalize"){
  write10xCounts(x=seurat.merged@assays$RNA@counts, path=paste0(savedir, "_seurat_",mergeType,"Merge_",normalization,"_RNAdata_rawCounts.h5"), version="3")
  write10xCounts(x=seurat.merged@assays$RNA@data, path=paste0(savedir, "_seurat_",mergeType,"Merge_",normalization,"_RNAdata_scaledCounts.h5"), version="3")

  print(Sys.time() - mid_time)
}

print("Adding Condition...")
Idents(seurat.merged) <- seurat.merged@meta.data$orig.ident
mid_time <- Sys.time()
## Add Sample annotations
idents <- data.frame(barcode = names(seurat.merged@active.ident), LibraryID = seurat.merged@active.ident)
#condition <- idents
#names(condition) <- c("barcode","Condition")
#condition$Condition <- as.character(condition$Condition)

#for(i in 1:dim(toProcess)[1]){
#  condition[condition$Condition == toProcess[i,1],"Condition"] <- toProcess[i,"Condition"]
#}

#seurat.merged$Condition <- condition[,"Condition",drop=FALSE]

#print("Adding TumorType...")
#mid_time <- Sys.time()
## Add Sample annotations
#idents <- data.frame(barcode = names(seurat.merged@active.ident), LibraryID = seurat.merged@active.ident)
#ttype <- idents
#names(ttype) <- c("barcode","TumorType")
#ttype$TumorType <- as.character(ttype$TumorType)

#for(i in 1:dim(toProcess)[1]){
#  ttype[ttype$TumorType == toProcess[i,1],"TumorType"] <- toProcess[i,"TumorType"]
#}

#seurat.merged$TumorType <- ttype[,"TumorType",drop=FALSE]


## Save Sample Annotations
write.csv(idents, file=paste0(savedir, "_Seurat_",mergeType,"Merge_", normalization, "_LibraryID.csv"), quote = FALSE, row.names = FALSE)
#write.csv(condition, file=paste0(savedir, "_Seurat_",mergeType,"Merge_", normalization, "_Condition1.csv"), quote = FALSE, row.names = FALSE)
#write.csv(ttype, file=paste0(savedir, "_Seurat_",mergeType,"Merge_", normalization, "_TumorType1.csv"), quote = FALSE, row.names = FALSE)


for(n in 1:(dim(toProcess)[2]-1)){
  idents <- data.frame(barcode = names(seurat.merged@active.ident), LibraryID = seurat.merged@active.ident)
  anno <- idents
  names(anno) <- c("barcode", names(toProcess)[n])
  anno[,2] <- as.character(anno[,2])
  
  for(i in 1:dim(toProcess)[1]){
    anno[anno[,2] == toProcess[i,1], 2] <- toProcess[i,n]
  }
  
  seurat.merged@meta.data[n] <- anno[,2,drop=FALSE]
  
  
  ## Save Sample Annotations
  write.csv(anno, file=paste0(savedir, "_Seurat_",mergeType,"Merge_", normalization, "_", names(toProcess)[n], ".csv"), quote = FALSE, row.names = FALSE)
  
  
}





print(Sys.time() - mid_time)

print("Running PCA, UMAP, tSNE...")
mid_time <- Sys.time()
### Create Visualizations
if(mergeType == "integration"){
  DefaultAssay(seurat.merged) <- "integrated"
} else {
  DefaultAssay(seurat.merged) <- "RNA"
}


#if(normalization == "LogNormalize"){
#  seurat.merged <- ScaleData(seurat.merged, verbose = FALSE)
#}

## Use features for PCA if provided.
if(features != ""){
  print("Subsetting PCA, UMAP, tSNE, and clustering on input features.")
  
  ## debugging issues with using a small gene list.
  #rna_assay <- seurat.merged@assays$RNA
  #rna_assay <- rna_assay[row.names(rna_assay) %in% my_feats,]
  
  #which(colSums(rna_assay) == 0)
  #rna_assay <- rna_assay[,-12]
  #rna_assay <- rna_assay[,-44]
  #Heatmap(cor(rna_assay))
  #std <- colSds(rna_assay)
  #rna_assay2 <- rna_assay[,which(std >= quantile(std)[3])]
  #colnames(rna_assay2)

  my_feats <- read.delim(features,header=FALSE)[,1]
  seurat.merged <- ScaleData(seurat.merged, verbose = FALSE, features = my_feats)
  seurat.merged <- RunPCA(seurat.merged, verbose = FALSE, features = my_feats, npcs = 50, approx=FALSE)
} else {
  seurat.merged <- RunPCA(seurat.merged, verbose = FALSE)
}

seurat.merged <- RunUMAP(seurat.merged, dims = 1:min(30,length(seurat.merged@reductions$pca)))
seurat.merged <- RunTSNE(seurat.merged, dims = 1:min(30,length(seurat.merged@reductions$pca)), check_duplicates = FALSE)

png(filename = paste0(savedir, "_pca.png"), res=150, width = 1100, height = 800)
DimPlot(seurat.merged, reduction = "pca", group.by="orig.ident")
dev.off()

png(filename = paste0(savedir, "_umap.png"), res=150, width = 1100, height = 800)
DimPlot(seurat.merged, reduction = "umap", group.by="orig.ident")
dev.off()

png(filename = paste0(savedir, "_tsne.png"), res=150, width = 1100, height = 800)
DimPlot(seurat.merged, reduction = "tsne", group.by="orig.ident")
dev.off()

if(features != ""){
  write.csv(seurat.merged@reductions$umap@cell.embeddings, file = paste0(savedir, "_UMAPCoordinates_30PCs_",mergeType,"Merge_",normalization,"_wFeatureSubset.csv"), quote = FALSE)
  write.csv(seurat.merged@reductions$tsne@cell.embeddings, file = paste0(savedir, "_tSNECoordinates_30PCs_",mergeType,"Merge_",normalization,"_wFeatureSubset.csv"), quote = FALSE)
  print(Sys.time() - mid_time)
} else {
  write.csv(seurat.merged@reductions$umap@cell.embeddings, file = paste0(savedir, "_UMAPCoordinates_30PCs_",mergeType,"Merge_",normalization,".csv"), quote = FALSE)
  write.csv(seurat.merged@reductions$tsne@cell.embeddings, file = paste0(savedir, "_tSNECoordinates_30PCs_",mergeType,"Merge_",normalization,".csv"), quote = FALSE)
  print(Sys.time() - mid_time)
}

print("SNN Clustering...")
mid_time <- Sys.time()
### Cluster the Data
seurat.merged <- FindNeighbors(seurat.merged, reduction = "pca", dims = 1:min(30,length(seurat.merged@reductions$pca)), graph.name = "merged_snn")
seurat.merged <- FindClusters(seurat.merged, resolution = 0.4, graph.name = "merged_snn")
seurat.merged <- FindClusters(seurat.merged, resolution = 0.6, graph.name = "merged_snn")

DimPlot(seurat.merged, reduction = "umap", group.by="merged_snn_res.0.6")
DimPlot(seurat.merged, reduction = "umap", group.by="merged_snn_res.0.4")


clusters1 <- data.frame("barcode" = names(seurat.merged$merged_snn_res.0.4), "SNN_res0.4_Clusters" = seurat.merged$merged_snn_res.0.4)
clusters2 <- data.frame("barcode" = names(seurat.merged$merged_snn_res.0.6), "SNN_res0.6_Clusters" = seurat.merged$merged_snn_res.0.6)

if(features != ""){
  write.csv(clusters1, file = paste0(savedir, "_SNN_Clusters_res0.4_",mergeType,"MergedData_wFeatureSubset.csv"), quote = FALSE, row.names = FALSE)
  write.csv(clusters2, file = paste0(savedir, "_SNN_Clusters_res0.6_",mergeType,"MergedData_wFeatureSubset.csv"), quote = FALSE, row.names = FALSE)
} else {
  write.csv(clusters1, file = paste0(savedir, "_SNN_Clusters_res0.4_",mergeType,"MergedData.csv"), quote = FALSE, row.names = FALSE)
  write.csv(clusters2, file = paste0(savedir, "_SNN_Clusters_res0.6_",mergeType,"MergedData.csv"), quote = FALSE, row.names = FALSE)
}

print(Sys.time() - mid_time)

print("Saving Loupe file...")
create_loupe_from_seurat(seurat.merged, output_dir = outDir, output_name = paste0(runID,"_Seurat_",mergeType,"Merge_",normalization,"_Annotated.cloupe"))


print("Saving Annotated Seurat File...")
mid_time <- Sys.time()

## Save Seurat Object for future use
if(saveH5){
  library(SeuratDisk)
  if(features != ""){
    SaveH5Seurat(seurat.merged, filename=paste0(savedir, "_Seurat_",mergeType,"Merge_",normalization,"_Annotated_wFeatureSubset.h5Seurat"), overwrite = TRUE)
  } else{
    SaveH5Seurat(seurat.merged, filename=paste0(savedir, "_Seurat_",mergeType,"Merge_",normalization,"_Annotated.h5Seurat"), overwrite = TRUE)
  }
}else{
  if(features != ""){
    save(seurat.merged, file = paste0(savedir, "_Seurat_",mergeType,"Merge_",normalization,"_Annotated_wFeatureSubset.RData"), compress = TRUE)
    saveRDS(seurat.merged, file = paste0(savedir, "_Seurat_",mergeType,"Merge_",normalization,"_Annotated_wFeatureSubset.rds"), compress = TRUE)
  } else {
    save(seurat.merged, file = paste0(savedir, "_Seurat_",mergeType,"Merge_",normalization,"_Annotated.RData"), compress = TRUE)
    saveRDS(seurat.merged, file = paste0(savedir, "_Seurat_",mergeType,"Merge_",normalization,"_Annotated.rds"), compress = TRUE)
  }
}




print(Sys.time() - mid_time)

print("Seurat merging completed in: ")
end_time <- Sys.time()
print(end_time - start_time)

