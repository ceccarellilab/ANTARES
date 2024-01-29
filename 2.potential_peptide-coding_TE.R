require(dplyr)
require(doParallel)
cl <- makeCluster(10)
registerDoParallel(cl)

length_ORF <- 10
length_ORF_id <- paste0("ORF", length_ORF)

res <- foreach(ii = 1:length(Repeats_seq_selected), .combine='c', .multicombine=TRUE) %dopar% {
  library(Biostrings)
  library(stringr)
  library(dplyr)
  seqs <- Repeats_seq_selected[ii]
  
  seqs_orf1 <- DNAStringSet((as.character(seqs )))
  seqs_orf2 <- DNAStringSet(str_sub(as.character(seqs), start = 2))
  seqs_orf3 <- DNAStringSet(str_sub(as.character(seqs), start = 3))
  seqs_orf4 <- reverseComplement(DNAStringSet(as.character(seqs)))
  seqs_orf5 <- reverseComplement(DNAStringSet(str_sub(as.character(seqs), end = -2)))
  seqs_orf6 <- reverseComplement(DNAStringSet(str_sub(as.character(seqs), end = -3)))
  
  list_sequences_DNA <- DNAStringSet(c(seqs_orf1, seqs_orf2, seqs_orf3, seqs_orf4, seqs_orf5, seqs_orf6))
  names(list_sequences_DNA) <- c(paste0("ORF", 1:3), paste0("ORF-", 1:3))
  
  list_sequences_AA <- AAStringSet(unlist(sapply(list_sequences_DNA, function(x) translate(x, if.fuzzy.codon = "solve"))))
  
  stop_set <- AAStringSet(x = c("*"))
  list_stop_results <- list()
  for (i in 1:6) {
    stop_results <- Biostrings::matchPDict(pdict = stop_set, subject =  list_sequences_AA[[i]], max.mismatch  = 0)
    list_stop_results[[i]] <- unlist(start(stop_results))
  }
  names(list_stop_results) <- names(list_sequences_AA)
  
  list_orf_results <- list()
  
  # canonical start codon : "ATG"
  # canonical + non canonical start codons : "ATG|CTG|GTG|TTG"
  library(ORFik)
  for (id_orf_loop in 1:6) {
    orf_df <- findORFs(list_sequences_DNA[id_orf_loop], startCodon = "ATG", minimumLength = 5, longestORF = TRUE) 
    
    if (length(orf_df) >= 1 ) {
      names(orf_df) <- names(list_sequences_DNA)[id_orf_loop]
      list_orf_results[[names(list_sequences_DNA)[id_orf_loop]]] <- orf_df[[1]]
    }
  }
  
  orf_df_test <- IRangesList(list_orf_results)
  
  orf_df  <- as.data.frame(orf_df_test) %>%
    dplyr::rename(ORF = group_name) %>%
    dplyr::filter(width >= length_ORF * 3) %>%
    dplyr::select(-group, -width) %>%
    dplyr::mutate(type = paste0("ORF", length_ORF), id = row_number()) %>%
    dplyr::ungroup()
  
  if(nrow(orf_df) == 0 ) return(NULL)
  
  DNA_seq <- subseq(list_sequences_DNA[orf_df$ORF], start = orf_df$start, end = orf_df$end)
  peptides <- translate(DNA_seq)
  
  
  five = orf_df[!grepl("-",orf_df$ORF),]
  three = orf_df[grepl("-",orf_df$ORF),]
  if(as.character(strand(Repeats_range_selected[ii])) == "+") {
    if(nrow(five) != 0) {
      tmp_five =   paste0(seqnames(Repeats_range_selected[ii]), ":",
                          (five$start - (-1*start(Repeats_range_selected[ii])  + 1)), "-",
                          (five$end - (-1*start(Repeats_range_selected[ii])  +1)), strand(Repeats_range_selected[ii]))
      names(tmp_five) = five$id
    } else {
      tmp_five = NULL
    }
    
    if(nrow(three) != 0) {
      tmp_three =   paste0(seqnames(Repeats_range_selected[ii]), ":",
                           (-1*three$end + (end(Repeats_range_selected[ii])  + 1)), "-",
                           (-1*three$start + (end(Repeats_range_selected[ii])  +1)), strand(Repeats_range_selected[ii]))
      names(tmp_three) = three$id
    } else {
      tmp_three = NULL
    }
  } else {
    if(nrow(three) != 0) {
      tmp_three = paste0(seqnames(Repeats_range_selected[ii]), ":",
                         (three$start - (-1*start(Repeats_range_selected[ii])  + 1)), "-",
                         (three$end - (-1*start(Repeats_range_selected[ii])  +1)), strand(Repeats_range_selected[ii]))
      names(tmp_three) = three$id
    } else {
      tmp_three = NULL
    }
    if(nrow(five) != 0) {
      tmp_five =   paste0(seqnames(Repeats_range_selected[ii]), ":",
                          (-1*five$end + (end(Repeats_range_selected[ii])  + 1)), "-",
                          (-1*five$start + (end(Repeats_range_selected[ii])  +1)), strand(Repeats_range_selected[ii]))
      names(tmp_five) = five$id
    } else {
      tmp_five = NULL
    }
  }
  
  tmp = c(tmp_five,tmp_three)
  
  
  df <- data.frame(RE_name = names(Repeats_seq_selected[ii]),
                   genomic_localization_RE = paste0(seqnames(Repeats_range_selected[ii]), ":",
                                                    start(Repeats_range_selected[ii]),"-",end(Repeats_range_selected[ii]), strand(Repeats_range_selected[ii])),
                   genomic_localization_ORF = tmp[as.character(orf_df$id)] ,
                   DNA_seq = as.character(DNA_seq),
                   peptides = as.character(peptides))
  df = as.list(df)
  list(df)
}
stopCluster(cl)

df_predicted_peptides = bind_rows(res)
df_predicted_peptides = df_predicted_peptides %>%
  group_by(genomic_localization_RE)%>%
  filter(!duplicated(peptides) )