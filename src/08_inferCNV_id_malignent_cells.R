## Amy Olex
## 12.20.22
## inferCNV cluster script
## https://github.com/broadinstitute/inferCNV/wiki/Running-InferCNV
## 

library(infercnv)
library(limma)
library(Biobase)
library(optparse)
library(Seurat)

option_list = list(
  make_option(c("-r", "--runid"), type="character", default=NULL, 
              help="Required. A unique name for this analysis.", metavar="character"),
  
  make_option(c("-o", "--outdir"), type="character", 
              help="output directory for results report (must already exist). Default is current directory.",
              default = "./", metavar="character"),
  
  make_option(c("-i", "--seuratfile"), type="character", default=NULL,
              help="Seurat object saved as an RData file in a variable named seurat.merged.", metavar="character"),
  
  make_option(c("-k", "--subsetK"), type="integer", default=-1,
              help="Number of samples to include per group. Default is -1, which uses all samples.", metavar="integer"),
  
  make_option(c("-f", "--phenofile"), type="character", default=NULL,
              help="A CSV file with the sample annotations must be provided.", metavar="character"),
  
  make_option(c("-g", "--generef"), type="character", default="./genes.gtf",
              help="Location of the genes.gtf reference genome file.", metavar="character")
  
); 

opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);

# test for NULL required arguments
if (is.null(opt$runid)){
  print_help(opt_parser)
  stop("A unique analysis name must be provided.", call.=FALSE)
}

if (is.null(opt$seuratfile)){
  print_help(opt_parser)
  stop("A Seurat object saved as a RData file must be provided.", call.=FALSE)
}


if (is.null(opt$phenofile)){
  print_help(opt_parser)
  stop("A CSV formatted sample phenotype file must be provided.", call.=FALSE)
}


runID <- opt$runid  ## "inferCNV_TEST_Visvader-MiniSubset_012023"
wd <- opt$outdir  ## "/Users/alolex/Desktop/CCTR_Git_Repos/WCCTR_RNASeq_Pipeline/SingleCell/debug_files"
seuratfile <- opt$seuratfile ## "100422_Visvader-MiniSubset_30_SimpleMerge_GRCh38_Seurat_simpleMerge_LogNormalize_Annotated.RData"
phenofile <- opt$phenofile  ## "inferCNV_TEST_Visvader-MiniSubset_012023_metaData.csv"
generef <- opt$generef
subsetK <- opt$subsetK ## 3

debug=FALSE
if(debug){
runID <-  "inferCNV_TEST_Visvader-MiniSubset_012023"
wd <-  "/Users/alolex/Desktop/CCTR_Git_Repos/WCCTR_RNASeq_Pipeline/SingleCell/debug_files"
seuratfile <-  "100422_Visvader-MiniSubset_30_SimpleMerge_GRCh38_Seurat_simpleMerge_LogNormalize_Annotated.RData"
phenofile <-  "inferCNV_TEST_Visvader-MiniSubset_012023_metaData.csv"
generef <- "genes.gtf"
subsetK <- 3
}

print("Summary of input options:\n")
print(paste("Run Name: ", runID))
print(paste("Output Dir:", wd))
print(paste("Seurat File: ", seuratfile))
print(paste("Phenotype File: ", phenofile))
print(paste("Gene GTF: ", generef))
print(paste("Subset K: ", subsetK))
setwd(wd)


##### Initiate custom functions ##########

#' @title make.eset
#' @description For the expression data are transformed to a file with extension .eSet
#' @usage make.eset(expr, pdata=NULL, fdata=NULL, verbose=TRUE)
#' @param expr expression data
#' @param pdata phenotype data
#' @param fdata feature data
#' @param verbose default TRUE, Class to writing verbose messages to a connection or file.
#' @details store expression data in ExpressionSet format for convenient analysis
#' @return expression set
#' @export
make.eset <- function(expr, pdata=NULL, fdata=NULL, verbose=TRUE){
  Sys.setlocale(category = "LC_ALL", locale = "us")
  
  if (inherits(expr, "data.frame")) expr=as.matrix(expr)
  if (inherits(expr, "ExpressionSet")) expr=as.matrix(exprs(expr))
  
  #for null data
  if(is.null(pdata))    pdata=data.frame(phenoID=colnames(expr), row.names=colnames(expr))
  if(is.null(rownames(pdata))) rownames(pdata)=colnames(expr)
  
  if(is.null(fdata))  fdata=data.frame(featureID=rownames(expr), row.names=rownames(expr))
  if(is.null(rownames(fdata))) rownames(fdata)=rownames(expr)
  
  #match
  if(ncol(pdata)>1){
    p.st=match(colnames(expr), rownames(pdata))
    p.col=colnames(pdata)
    class(pdata)
    pdata=data.frame(pdata[p.st,p.col],row.names=rownames(pdata)[p.st])
    colnames(pdata)=p.col
  }
  
  
  f.st=match(rownames(expr), rownames(fdata))
  f.col=colnames(fdata)
  fdata=data.frame(fdata[f.st,], row.names=rownames(fdata)[f.st])
  colnames(fdata)=f.col
  
  #create eset
  metadata <- data.frame(labelDescription = colnames(pdata), row.names=colnames(pdata))
  phenoData<-new("AnnotatedDataFrame", data=as.data.frame(pdata), varMetadata=metadata)
  fmetadata <- data.frame(labelDescription = colnames(fdata), row.names=colnames(fdata))
  featureData<-new("AnnotatedDataFrame", data=as.data.frame(fdata), varMetadata=fmetadata)
  
  eset<-new("ExpressionSet", exprs=expr, phenoData=phenoData, featureData=featureData)
  if(verbose) print(eset)
  return(eset)
}


######### Currently not using this function.
#' @title malignant.cellTyper
#' @description A function to malignant cell typing
#' @param seurat Seurat object
#' @param rda.dir rData directory
#' @param malignant.cell.type Cell type to assign malignant cell
#' @param feature.to.test features to test as reference
#' @param cells.test_reference cells to test as reference
#' @details classification of malignant and non malignant seurat object.
#' @return Seurat object
#' @export
malignant.cellTyper <- function(seurat,
                                rda.dir = "./data",
                                malignant.cell.type="Epithelial",
                                feature.to.test = c("cell.type","tissue.type"),
                                cells.test_reference="immune"){
  
  
  message("[[",Sys.time(),"]] Run malignant.cellTyper --------")
  cell.type=as.character(seurat$cell.type)
  mal.fil.st = cell.type==malignant.cell.type
  
  if(feature.to.test=="tissue.type"){
    cnv.cut=quantile(seurat$cnv.score[seurat$cell.group %in% cells.test_reference], probs=0.90, na.rm = TRUE)
  }else if(feature.to.test=="cell.type"){
    cnv.cut=quantile(seurat$cnv.score[seurat$cell.group %in% c(cells.test_reference, "Unresolved_cell")], probs=0.90, na.rm = TRUE)
  }
  
  seurat$cnv.st=seurat$cnv.score>cnv.cut
  cnv.fil.st=seurat$cnv.st
  
  fil.st= (mal.fil.st | cnv.fil.st)
  table(fil.st)
  
  seurat$malignant.st=fil.st
  seurat$nonmalignant.st = !fil.st
  
  cell.type[fil.st]="Malignant_cell"
  seurat$cell.type=as.factor(cell.type)
  
  save(seurat, file = file.path(rda.dir, "seurat.rda"))
  message("[[",Sys.time(),"]] Finish malignant.cellTyper --------")
  return(seurat)
}

############ End custom functions ##############









## Load in data
load(seuratfile)
samples <- read.csv(phenofile)
assay = "RNA"
print("data loaded...")
print(paste("Using assay:", assay))

## create new annotation for each cell
seurat.merged@meta.data$cell.group <- samples[match(seurat.merged@meta.data$orig.ident, samples[,"Sample_ID"]),"inferCNV.normal"]

############## Need to automate: 
## Subset to chosen samples and SELECTED_Normal samples
reference_list <- unique(seurat.merged$orig.ident[seurat.merged$cell.group %in% c("REFERENCE")])
sample_list <- unique(seurat.merged$orig.ident[!(seurat.merged$cell.group %in% c("REFERENCE"))])
if(subsetK < 0){
  sample_groups <- list(sample_list)
  }else{
    sample_groups <- split(sample_list, ceiling(seq_along(sample_list)/subsetK))
  }

barcode_order <- names(seurat.merged$orig.ident)

Idents(seurat.merged) <- "orig.ident"

## create an empty DF here to hold the pData information.
pData_combined <- data.frame()
refData <- ""
iter <- 0

## Save the sample groups
write.csv(t(as.data.frame(sample_groups)), paste0(runID, "_sample_groups.csv"))

for(sg in sample_groups){
  iter <- iter+1
  print("Processing sample group: " )
  sample_subset <- c(sg, reference_list)
  print(sample_subset)
  
  ## subset Seurat object
  seurat.subset <- subset(seurat.merged, idents = sample_subset)
  
  ## Format gene expression input for inferCNV
  ## 1) create raw counts matrix
  raw.mat=as.matrix(seurat.subset[[assay]]@counts)
  
  ## 2) format annotation data
  annot=data.frame(cell.names=colnames(seurat.subset), group.st=seurat.subset@meta.data$cell.group)
  
  ## 3) get gene feature data (code from scTyper)
  ##making fdata object
  gene.ref.gtf=generef
  gtf=read.delim(gene.ref.gtf, comment.char = "#", header = F)
  gene.gtf=gtf[gtf$V3=="gene",]
  gene.annot=strsplit2(gene.gtf$V9, split = ";| ")
  gene.annot=data.frame("gene_id"=gene.annot[,2], "gene_name"=gene.annot[,8], "gene_source"=gene.annot[,11], "gene_biotype"=gene.annot[,14])
  head(gene.annot)
  mat.st=match(rownames(seurat.subset), gene.annot$gene_name)
  mat.st[is.na(mat.st)]=match(sub("\\.1", "", rownames(seurat.subset)[is.na(mat.st)]), gene.annot$gene_name)
  
  fdata=data.frame(gene.gtf[mat.st,c(1,4,5,7)], gene.annot[mat.st,])
  colnames(fdata)[1:4]=c("chr", "str", "end", "strand")
  rownames(fdata)=rownames(seurat.subset)
  
  ## adding gene annotations to Seurat obj
  seurat.subset[[assay]]@meta.features <- data.frame(seurat.subset[[assay]]@meta.features, fdata)
  names(seurat.subset[[assay]]@meta.features)
  gene.od=data.frame(gene=rownames(seurat.subset[[assay]]@meta.features), seurat.subset[[assay]]@meta.features[,c("chr", "str", "end")])
  fil.st=is.element(gene.od$chr, c(1:22, "X"))
  
  ## filtering to only primary chromosomes.
  raw.mat=raw.mat[fil.st,]
  gene.od=gene.od[fil.st,]
  gene.od$chr=factor(as.character(gene.od$chr), levels=c(1:22,"X"))
  od.st=order(gene.od$chr, gene.od$str)
  match(c(1:22, "X"), gene.od[od.st,]$chr)
  
  # create and write out infercnv formatted files
  output.dir = wd
  raw.mat.fn=file.path(output.dir, paste0(runID, "_raw_counts_matrix.txt"))
  annot.fn=file.path(output.dir, paste0(runID, "_annotations.txt"))
  gene.od.fn=file.path(output.dir, paste0(runID, "_gene_order.txt"))
  write.table(raw.mat[od.st,], raw.mat.fn, quote = FALSE, sep = "\t", row.names = TRUE, col.names = TRUE)
  write.table(annot, annot.fn, quote = FALSE, sep = "\t", row.names = FALSE, col.names = FALSE)
  write.table(gene.od[od.st,], gene.od.fn, quote = FALSE, sep = "\t", row.names = FALSE, col.names = FALSE)
  
  
  ## ###########
  ## Run inferCNV!
  ###########
  local_out_dir <- paste0("local_output_infercnv_",iter)
  
  print("creating inferCNV object...")
  cells.test_reference = c("REFERENCE")
  infercnv_obj = infercnv::CreateInfercnvObject(raw_counts_matrix=raw.mat[od.st,],
                                                annotations_file=annot.fn,
                                                gene_order_file=gene.od.fn,
                                                ref_group_names=cells.test_reference)
  ## This took 3.5 hours for about 10k cells!
  print("running inferCNV command...")
  infercnv_obj = infercnv::run(infercnv_obj,
                               cutoff=0.1,  # use 1 for smart-seq, 0.1 for 10x-genomics
                               out_dir=local_out_dir,  # dir is auto-created for storing outputs
                               cluster_by_groups=T,   # cluster
                               denoise=T,
                               HMM=T,
                               num_threads=30)
  
  save(infercnv_obj, file=file.path(paste0(local_out_dir,"/"), 'infercnv_obj.o.rda'))
  
  print("completed inferCNV...")
  
  ### Filter out by standard deviations in reference groups
  ref.sd=apply(infercnv_obj@expr.data[,infercnv_obj@reference_grouped_cell_indices[[1]]], 1, sd)
  sd.cut = 0.3
  sd.fil.st <- ref.sd < sd.cut
  infercnv_obj@expr.data=infercnv_obj@expr.data[sd.fil.st,]
  infercnv_obj@count.data=infercnv_obj@count.data[sd.fil.st,]
  infercnv_obj@gene_order=infercnv_obj@gene_order[sd.fil.st,]
  
  save(infercnv_obj, file=file.path(paste0(local_out_dir,"/"), 'infercnv_obj.sd_filtered.rda'))
  save(infercnv_obj, file=file.path(paste0(local_out_dir,"/"), 'infercnv_obj.rda'))
  
  #make filtered cnvset
  
  fil.cset=make.eset(expr = infercnv_obj@expr.data, pdata = seurat.merged@meta.data[colnames(infercnv_obj@expr.data),], fdata = infercnv_obj@gene_order)
  
  pData(fil.cset)$cnv.score=apply(exprs(fil.cset)-1, 2, function(a) sum(abs(a)))
  
  #save(fil.cset, file=file.path(paste0(local_out_dir,"/"), 'fil.cset.rda'))
  
  ######
  ## Create Z-score version of data minus the normal samples
  ##z_scores <- (data-mean(data))/sd(data)
  print("Creating ZScores..." )
  ref_sd <- sd(pData(fil.cset)[!(pData(fil.cset)$orig.ident %in% sg),"cnv.score"])
  pData(fil.cset)$cnv.Zscore <- (pData(fil.cset)$cnv.score - mean(pData(fil.cset)$cnv.score))/ref_sd
  pData(fil.cset)$cnv.ZscoreCat <- round(pData(fil.cset)$cnv.Zscore, 0)
  pData_noref <- pData(fil.cset)[pData(fil.cset)$orig.ident %in% sg,]
  
  ##merge sample DFs together
  pData_combined <- rbind(pData_combined, pData_noref)
  
  save(fil.cset, file=file.path(paste0(local_out_dir,"/"), 'fil.cset.rda'))
  
  ##save most recent ref data values (maybe average these one day?)
  refData <- pData(fil.cset)[!(pData(fil.cset)$orig.ident %in% sg),]

  
}

### Combine the pData and refData and then re-order to match the original order so we can add it as an annotation in the Seruat obj
print("Saving final files..." )

multiData <- rbind(refData, pData_combined)
multiData <- multiData[barcode_order,]

if(all(row.names(multiData) == names(seurat.merged$orig.ident))){
  ### Add to Seurat obj
  seurat.merged$cnv.score=multiData$cnv.score
  seurat.merged$cnv.Zscore=multiData$cnv.Zscore
  seurat.merged$cnv.ZscoreCat=multiData$cnv.ZscoreCat

  ##write out the annotation data and files
  write.csv(data.frame(barcode = names(seurat.merged@active.ident), inferCNV.cnvscore = seurat.merged@meta.data$cnv.score), file = paste0(runID, "_inferCNVscores.csv"), quote=FALSE, row.names=FALSE)
  write.csv(data.frame(barcode = names(seurat.merged@active.ident), inferCNV.cnvZscore = seurat.merged@meta.data$cnv.Zscore), file = paste0(runID, "_inferCNVZscores.csv"), quote=FALSE, row.names=FALSE)
  write.csv(data.frame(barcode = names(seurat.merged@active.ident), inferCNV.cnvZscoreCat = seurat.merged@meta.data$cnv.ZscoreCat), file = paste0(runID, "_inferCNVZscoreCat.csv"), quote=FALSE, row.names=FALSE)
}else{
  print("ERROR: barcode order did not match!")
}

  ###save Seurat obj with cnv scores
  print("saving Seurat...")
  saveRDS(seurat.merged, paste0(runID, ".inferCNVscores.rda"))

  print("Saving images...")
  png(filename = paste0(runID, ".inferCNVscores.png"))
    FeaturePlot(seurat.merged, features = "cnv.score")
  dev.off()

  png(filename = paste0(runID, ".inferCNVZscores.png"))
  FeaturePlot(seurat.merged, features = "cnv.Zscore")
  dev.off()

  png(filename = paste0(runID, ".inferCNVZscoreCat.png"))
  FeaturePlot(seurat.merged, features = "cnv.ZscoreCat")
  dev.off()

print("COMPLETED!")
