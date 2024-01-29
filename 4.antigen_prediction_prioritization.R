ffiles = list.files("netmhc_predictions/",pattern="fc058.txt", full.names = T)

require(doParallel)
cl <- makeCluster(10)
registerDoParallel(cl)

res <- foreach(ii = 1:length(ffiles), .combine='c', .multicombine=TRUE) %dopar% {
  require(Biostrings)
  netMHCIIpan_predictions_SBWB = read.table(ffiles[ii], quote="\"", comment.char="")
  netMHCIIpan_predictions_SBWB = netMHCIIpan_predictions_SBWB[,-1]
  netMHCIIpan_predictions_SBWB = netMHCIIpan_predictions_SBWB[,-(ncol(netMHCIIpan_predictions_SBWB)-1)]
  netMHCIIpan_predictions_SBWB = netMHCIIpan_predictions_SBWB[,c(1:3,10,13)]
  colnames(netMHCIIpan_predictions_SBWB) = c("HLA","Neoantigen","Core","Seq","Class_Affinity")
  
  netMHCIIpan_predictions_SBWB = data.frame(netMHCIIpan_predictions_SBWB,
                                            df_predicted_peptides_filtered[netMHCIIpan_predictions_SBWB$Seq,-ncol(df_predicted_peptides_filtered)], 
                                            stringsAsFactors = F)
  dna_pos_neo_all = NULL
  for(p in 1:nrow(netMHCIIpan_predictions_SBWB)) {
    
    pep = AAString(netMHCIIpan_predictions_SBWB$peptides[p])
    neo = AAString(netMHCIIpan_predictions_SBWB$Neoantigen[p])
    dnaseq = DNAString(netMHCIIpan_predictions_SBWB$DNA_seq[p])
    dna_pos= netMHCIIpan_predictions_SBWB$genomic_localization_ORF[p]
    pos = matchPattern(pattern = neo, subject = pep, max.mismatch = 0)
    
    if(grepl("+",dna_pos)) {
      sstart = as.numeric(strsplit(strsplit(dna_pos,":")[[1]][2],"-")[[1]][1])
      a = sstart + (start(pos)*3 -3)
      b = sstart + (end(pos)*3 + 2)
      dna_pos_neo = paste0(strsplit(dna_pos,":")[[1]][1],":",a,":", b,"+")  
    } else {
      eend = as.numeric(strsplit(strsplit(dna_pos,":")[[1]][2],"-")[[1]][2])
      a = eend - (end(pos)*3 )
      b = eend - (start(pos)*3 - 3)
      dna_pos_neo = paste0(strsplit(dna_pos,":")[[1]][1],":",a,":", b,"-")  
    }
    dna_pos_neo_all = c(dna_pos_neo_all,dna_pos_neo)
  }
  netMHCIIpan_predictions_SBWB$genomic_localization_Neoantigen = dna_pos_neo_all
  df = as.list(netMHCIIpan_predictions_SBWB)
  list(df)
}
stopCluster(cl)

netMHC_list_HLA = lapply(res, function(x) as.data.frame(x))
names(netMHC_list_HLA) = unname(sapply(ffiles, function(i) strsplit(i,"_")[[1]][4]))


netMHC_list_HLA_SB = lapply(netMHC_list_HLA, function(x) x[grepl("SB",x$Class_Affinity),])
netMHC_list_HLA_SB = lapply(netMHC_list_HLA_SB, function(x) x[!grepl("_",x$genomic_localization_RE),])
all = do.call("rbind",netMHC_list_HLA_SB)
all$Id_RE = paste0(all$RE_name,"_", all$genomic_localization_RE)
tmp <- AAStringSet(all$Neoantigen)
names(tmp) = paste0(rownames(all),"__",all$Neoantigen,"|","MUT","|", 1:length(tmp))
Biostrings::writeXStringSet(tmp, "antigen_sb.fasta")

# blastp -query antigen_sb.fasta -db iedb.fasta -outfmt 5 -evalue 100000000 -matrix BLOSUM62 -gapopen 11 -gapextend 1 > antigen_iedb.xml

# Source code to compute neoantigen fitness cost from https://doi.org/10.1038/nature24473 --> output dissimilarity.csv

dissimilarity <- read.csv("dissimilarity.csv")
# hist(dissimilarity$value)
dissimilarity = dissimilarity[which(dissimilarity$value ==1),]

final = all[dissimilarity$id,]
