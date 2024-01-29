df_predicted_peptides$Id = paste0(df_predicted_peptides$RE_name,"_",seq(1:nrow(df_predicted_peptides)))
df_predicted_peptides = as.data.frame(df_predicted_peptides)
rownames(df_predicted_peptides) = df_predicted_peptides$Id

LTR = Repeats_range_selected[Repeats_range_selected$repClass == "LTR",]
query = df_predicted_peptides[df_predicted_peptides$RE_name %in% unique(LTR$repName),]
query = query[,c("Id","peptides")]
colnames(query)[2] = "Seq"
query$Seq = gsub("\\*","",query$Seq)
tmp <- AAStringSet(query$Seq) 
names(tmp) = query$Id
Biostrings::writeXStringSet(tmp, "LTR.fasta")

LINE = Repeats_range_selected[Repeats_range_selected$repClass == "LINE",]
query = df_predicted_peptides[df_predicted_peptides$RE_name %in% unique(LTR$repName),]
query = query[,c("Id","peptides")]
colnames(query)[2] = "Seq"
query$Seq = gsub("\\*","",query$Seq)

tmp <- AAStringSet(query$Seq) 
names(tmp) = query$Id
Biostrings::writeXStringSet(tmp, "LINE.fasta")

# conda activate blast
# blastp -query LTR.fasta -db db_LTR -outfmt "6 qseqid qlen sseqid slen length pident evalue bitscore" -word_size 3 -num_threads 100 -seg no -max_hsps 1 -max_target_seqs 1 | sort -nrk 8 > blast-results2.tmp
# cat <(printf "qseqid\tqlen\tsseqid\tslen\tlength\tpident\tevalue\tbitscore\n") blast-results2.tmp > blast-LTR.tsv
# rm blast-results2.tmp

# blastp -query LINE.fasta -db db_LINE -outfmt "6 qseqid qlen sseqid slen length pident evalue bitscore" -word_size 3 -num_threads 100 -seg no -max_hsps 1 -max_target_seqs 1 | sort -nrk 8 > blast-results2.tmp
# cat <(printf "qseqid\tqlen\tsseqid\tslen\tlength\tpident\tevalue\tbitscore\n") blast-results2.tmp > blast-LINE.tsv
# rm blast-results2.tmp

blast_LTR <- read.delim("blast-LTR.tsv")
blast_LTR = blast_LTR[blast_LTR$pident >= 90 & blast_LTR$evalue <= 0.01,]
blast_LINE <- read.delim("blast-LINE.tsv")
blast_LINE = blast_LINE[blast_LINE$pident >= 90 & blast_LINE$evalue <= 0.01,]

sel = c(blast_LTR$qseqid,blast_LINE$qseqid)
df_predicted_peptides_filtered = df_predicted_peptides[sel,]

tmp <- AAStringSet(df_predicted_peptides_filtered$peptides) 
names(tmp) = df_predicted_peptides_filtered$Id
Biostrings::writeXStringSet(tmp, "predicted_peptides_filtered.fasta")