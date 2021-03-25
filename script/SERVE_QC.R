args <- commandArgs(T)
prefix=args[1]
count=args[2]

QC_dir="3_qc/"
Assem_dir="2_assem/"


overlap_file=paste0(QC_dir, prefix, "_overlap.txt")
expr_file=paste0(Assem_dir, "RSEM.isoforms.results")

overlap=read.table(overlap_file, stringsAsFactors=F)
expr=read.table(expr_file, header=T, stringsAsFactors=F)

overlap <- overlap[,c(1:5,19,7:13)]
overlap$direction <- paste0(overlap[,4], "_", overlap[,6])
overlap <- overlap[order(overlap$direction),]
overlap <- overlap[!duplicated(overlap$direction),]
overlap[duplicated(overlap[,4]),6] <- "."
overlap <- overlap[!duplicated(overlap[,4], fromLast=T),]

overlap$transcript_id <- t(as.data.frame(strsplit(overlap[,4], split=".", fixed=T)))[,1]

final=merge(overlap, expr, all.x=T, by='transcript_id', sort=F)

keep <- final$expected_count>count
final[,14] <- paste0(final[,14],";TPM=",final$TPM)
final=final[keep,c(2:14)]

out_bed=paste0(QC_dir, prefix, "_ERV.bed")

write.table(final, file=out_bed, sep='\t', row.names=F, col.names=F, quote=F)
