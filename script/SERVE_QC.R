args <- commandArgs(T)
prefix=args[1]
count=as.numeric(args[2])

QC_dir="3_qc/"
Assem_dir="2_assem/"


overlap_file=paste0(QC_dir, prefix, "_overlap.txt")
expr_file=paste0(Assem_dir, "RSEM.isoforms.results")
bed_file=paste0(QC_dir, prefix, "_ERV.bed")

overlap=read.table(overlap_file, stringsAsFactors=F)
expr=read.table(expr_file, header=T, stringsAsFactors=F)
bed=read.table(bed_file, stringsAsFactors=F)

overlap <- overlap[,c(4,12)]
overlap$direction <- paste0(overlap[,1], "_", overlap[,2])
overlap <- overlap[order(overlap$direction),]
overlap <- overlap[!duplicated(overlap$direction),]
overlap[duplicated(overlap[,1]),2] <- "."
overlap <- overlap[!duplicated(overlap[,1], fromLast=T),]
names(bed)[4]=names(overlap)[1]="ID"

bed <- merge(bed, overlap[,1:2], by="ID", sort=F)
bed[,6] <- bed[,ncol(bed)]

bed$transcript_id <- t(as.data.frame(strsplit(bed$ID, split=".", fixed=T)))[,1]

final=merge(bed, expr, all.x=T, by='transcript_id', sort=F)

keep <- final$expected_count>count
final[,14] <- paste0(final[,14],";TPM=",final$TPM)
final=final[keep,c(3:5,2,6:14)]

write.table(final, file=bed_file, sep='\t', row.names=F, col.names=F, quote=F)
