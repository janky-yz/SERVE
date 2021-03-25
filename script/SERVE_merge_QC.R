args <- commandArgs(T)
mode=as.numeric(args[1])
prefix=args[2]
ratio=as.numeric(args[3])
length=as.numeric(args[4])

compare_file=paste0(prefix, "_cuff.loci")
gtf_file=paste0(prefix, "_ERV_merge.gtf")
gff3_file=paste0(prefix, "_gmap.gff3")

if(mode==2){
	gtf=read.table(gtf_file, stringsAsFactors=F)
	gff3=read.table(gff3_file, stringsAsFactors=F, sep='\t')

	gff3_gene <- gff3[gff3[,3]=='gene',]
	gff3_gene_remove <- gff3_gene[grep(pattern='path2', x=gff3_gene[,9]),9]
	gff3_gene_remove_list <- t(as.data.frame(strsplit(gff3_gene_remove, '.', fixed=T)))[,1]
	gff3_gene_remove_list <- t(as.data.frame(strsplit(gff3_gene_remove_list, '=', fixed=T)))[,2]
	gff3_gene_remove_list <- unique(sort(gff3_gene_remove_list))
	gene_remove_list <- unique(sort(gtf[gtf[,10] %in% gff3_gene_remove_list,13]))

	gtf_filt <- gtf[!gtf[,13] %in% gene_remove_list,]
	gtf_filt[,9] <- paste0(gtf_filt[,9], ' "', gtf_filt[,10], '"; ', gtf_filt[,12], ' "', gtf_filt[,13], '";')
	write.table(gtf_filt[,1:9], file=gtf_file, sep='\t', col.names=F, row.names=F, quote=F)

	files_remove <- dir(pattern="gmap")
	unlink(files_remove, force=TRUE)
}else{
	compare=read.table(compare_file, stringsAsFactors=F)
	gtf=read.table(gtf_file, stringsAsFactors=F)

	nsample=ncol(compare)-3
	compare <- compare[compare[,3]!='-',]
	keep <- rowSums(compare[,-c(2:3)]=="-")<=ratio*nsample
	ERV_list <- unlist(strsplit(compare[keep,3], ',', fixed=T))
	ERV_list <- t(as.data.frame(strsplit(ERV_list, '|', fixed=T)))[,1]
	ERV_list <- unique(sort(ERV_list))

	gtf_filt <- gtf[gtf[,13] %in% ERV_list,]
	gtf_filt[,9] <- paste0(gtf_filt[,9], ' "', gtf_filt[,10], '"; ', gtf_filt[,12], ' "', gtf_filt[,13], '";')
	write.table(gtf_filt[,1:9], file=gtf_file, sep='\t', col.names=F, row.names=F, quote=F)

	gtf_exon <- gtf_filt[gtf_filt[,3]=='exon',]
	gtf_exon$length <- gtf_exon[,5]-gtf_exon[,4]

	transcript_length <- aggregate(gtf_exon$length, list(gtf_exon[,10]), sum)
	colnames(transcript_length) <- c("transcript_ID", "length")

	names(gtf_filt)[c(10,13)] <- names(gtf_exon)[c(10,13)] <- c("transcript_ID", "gene_ID")
	transcript_length <- merge(transcript_length, gtf_exon[!duplicated(gtf_exon[,10]),c(10,13)], by='transcript_ID')
	transcript_length <- transcript_length[order(transcript_length$gene_ID, transcript_length$length, decreasing=T),]
	transcript_length <- transcript_length[!duplicated(transcript_length$gene_ID),]
	transcript_length <- transcript_length[transcript_length$length>=length,]

	gtf_longest_transcript <- gtf_filt[gtf_filt$transcript_ID %in% transcript_length$transcript_ID,]

	gmap_file=paste0(prefix, "_gmap.gtf")
	write.table(gtf_longest_transcript[,1:9], file=gmap_file, sep='\t', col.names=F, row.names=F, quote=F)

	files_remove <- dir(pattern="cuff")
	unlink(files_remove, force=TRUE)
}
