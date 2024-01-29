### required packages and functions ###
require(AnnotationHub)
library(dplyr)
library(GenomicFeatures)
require(diffloop)
require(DESeq2)
library(BSgenome.Hsapiens.UCSC.hg38)

getDEGsEdgeR <- function(ddata, interestGroup, baselineGroup, min_num_cpm = 3){ 
  library(edgeR)
  ddata <- ddata[, c(interestGroup, baselineGroup)]
  groups <- rep("Interest", ncol(ddata))
  names(groups) <- colnames(ddata)
  groups[names(groups) %in% baselineGroup] <- "Baseline"
  groups <- as.factor(groups)
  baseLine <- "Baseline"
  
  design <- model.matrix(~ 0 + groups)
  colnames(design) <- levels(groups)
  
  DEA <- DGEList(ddata, group = groups)
  DEA <- calcNormFactors(DEA)
  #keep genes that achieve at least une count per million (cpm) in at least ten samples
  keep <- rowSums(cpm(DEA) > 1) >= min_num_cpm
  DEA <- DEA[keep,]
  DEA$samples$lib.size <- colSums(DEA$counts)
  message("Estimating GLM Common Dispersion")
  DEA <- estimateGLMCommonDisp(DEA, design = design)
  message("Estimating GLM Tagwise Dispersion")
  DEA <- estimateGLMTagwiseDisp(DEA, design = design)
  message("Fitting the GLM Model")
  fit <- glmFit(DEA, design = design)
  contrast <- rep(1, 2)
  blId <- which(colnames(design) == baseLine)
  contrast[blId] <- -1
  message("Maximum Likelihood Estimate")
  test <- glmLRT(fit, contrast = contrast)
  ans <- topTags(test, n = nrow(DEA$counts))
  ans <- ans$table
  return(ans)
}

### Get genomic regions for selected TE ###
hub <- AnnotationHub()
repeatGr <- hub[[names(AnnotationHub::query(hub,c("rmsk", "Homo sapiens", "hg38")))]]
sel <- c("Alu","MIR","ERV1","ERVK","ERVL","ERVL-MaLR","L1","L2")
Repeats_class = repeatGr[repeatGr$repFamily %in% sel]

ah = AnnotationHub()
AnnotationHub::query(ah, c("gtf", "Homo_sapiens", "hg38"))
GRCh38.gtf<- ah[['AH28606']]
GRCh38.txdb <- makeTxDbFromGRanges(GRCh38.gtf)

ex_range = exons(GRCh38.txdb) %>% GenomicRanges::reduce() 
ex_range = addchr(ex_range)
genome(ex_range) <- "hg38"

L1 = Repeats_class
L1 = L1[!grepl("_",seqnames(L1))]
seqlevels(L1) <- seqlevelsInUse(L1)
tmp = countOverlaps(L1, ex_range)
L1 = L1[tmp == 0]

Repeats_range_selected = Repeats_class[Repeats_class$repName %in% L1$repName,]
namess = paste0(Repeats_range_selected$repName,"_",seqnames(Repeats_range_selected), ":",start(Repeats_range_selected),"-",end(Repeats_range_selected), strand(Repeats_range_selected))
names(Repeats_range_selected) = namess
file_range = rtracklayer::export(Repeats_range_selected,"TE_ranges.gtf")

### Get TE-Family count for each sample using featureCounts ###
## example:
# for file in `ls path_of_bam_files/*.bam`
#   do
#   fname=${file%-*}
#   fname=$(basename $fname)
#   featureCounts -p -t 'sequence_feature' -g 'ID' --ignoreDup -M -T 40 --fraction -a TE_ranges.gtf -o output_featurecounts/${fname}_featureCounts.txt $file 
# done

output = list.files("output_featurecounts_per_sample/", pattern = ".txt", full.names = T)
names(output) = gsub(".txt","",basename(output))
featureCounts = vector("list", length(output))
names(featureCounts) = names(output)

for(i in names(featureCounts)){
  file = output[i]
  if(length(file) == 0) next
  tmp = read.delim(file, row.names=1, comment.char="#")
  colnames(tmp) = i
  featureCounts[[i]] = tmp
}

### Get normalized count matrix per f ###

featureCounts_matrix = do.call(cbind, featureCounts)
dds <- DESeqDataSetFromMatrix(countData = round(featureCounts_matrix), 
                              colData = pdata,
                              design = ~ Sample_Type ) #Sample_Type column contains the info of Tumor/Normal samples
dds <- estimateSizeFactors(dds)
normalized_counts <- counts(dds, normalized=TRUE)


source("custom_functions.R")
#interestGroup vector of samples in interest condition (tumor)
#baselineGroup vector of samples in control condition (normal)

DEGs = getDEGsEdgeR(normalized_counts, interestGroup, baselineGroup, min_num_cpm = 3)
select_families = DEGs[DEGs$logFC >= log2(1.5)  & DEGs$FDR <= 0.05,]
select_families = rownames(select_families)

#### Get TE-name count for each sample using featureCounts ###

Repeats_range_selected = Repeats_class[Repeats_class$repName %in% select_families,]
namess = paste0(Repeats_range_selected$repName,"_",seqnames(Repeats_range_selected), ":",start(Repeats_range_selected),"-",end(Repeats_range_selected), strand(Repeats_range_selected))
names(Repeats_range_selected) = namess
file_range = rtracklayer::export(Repeats_range_selected,"sel_repName.gtf")

output = list.files("output_featurecounts/", full.names = T)
names(output) = gsub(".txt","",basename(output))
featureCounts = vector("list", length(output))
names(featureCounts) = names(output)

for(i in names(featureCounts)){
  file = output[i]
  if(length(file) == 0) next
  tmp = try(read.delim(file, row.names=1, comment.char="#"),silent = T)
  if(class(tmp) == "try-error") next
  colnames(tmp) = i
  featureCounts[[i]] = tmp
}

featureCounts_matrix = do.call(cbind, featureCounts)
dds <- DESeqDataSetFromMatrix(countData = round(featureCounts_matrix), 
                              colData = pdata,
                              design = ~ Sample_Type ) #Sample_Type column contains the info of Tumor/Normal samples
dds <- estimateSizeFactors(dds)
normalized_counts <- counts(dds, normalized=TRUE)
normalized_counts = normalized_counts[rowSums(normalized_counts) != 0,]
Repeats_range_selected = Repeats_range_selected[names(Repeats_range_selected) %in% rownames(normalized_counts)]

## Get sequence from reference fasta ###

human_hg38 = BSgenome.Hsapiens.UCSC.hg38
Repeats_seq_selected = getSeq(human_hg38, Repeats_range_selected)
names(Repeats_seq_selected) = Repeats_range_selected$repName
