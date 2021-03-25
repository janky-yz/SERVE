args <- commandArgs(T)

file_list=args[1]
prefix=args[2]
count_cutoff=as.numeric(args[3])
TPM_cutoff=as.numeric(args[4])
ratio_cutoff=as.numeric(args[5])

TPM_output=paste0(prefix, "_ERV_TPM.txt")
count_output=paste0(prefix, "_ERV_count.txt")

files <- read.table(file_list, stringsAsFactors=F)
files <- files[,1]

for(fid in 1:length(files))
{
	quan=read.table(files[fid], stringsAsFactors=F, header=T)
	if(fid==1){
		merge_TPM <- quan[,c(1,6)]; merge_count <- quan[,c(1,5)]
	}else{
		merge_TPM <- data.frame(merge_TPM, quan[,6])
		merge_count <- data.frame(merge_count, quan[,5])
	}
}

rownames(merge_TPM)=rownames(merge_count)=merge_TPM[,1]
merge_TPM <- merge_TPM[,-1]
merge_count <- merge_count[,-1]
colnames(merge_TPM)=colnames(merge_count)=files

keep <- rowSums(merge_TPM>=TPM_cutoff)>=length(files)*ratio_cutoff & rowSums(merge_count>count_cutoff)>=length(files)*ratio_cutoff

filt_TPM <- merge_TPM[keep,]
filt_count <- merge_count[keep,]

write.table(filt_TPM, file=TPM_output, sep='\t', col.names=T, row.names=T, quote=F)
write.table(filt_count, file=count_output, sep='\t', col.names=T, row.names=T, quote=F)

