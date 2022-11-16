#annotate SNP list
#get ensembl TSS

library(biomaRt)
ensembl = useEnsembl(biomart="ensembl", dataset="hsapiens_gene_ensembl", GRCh=37)

#positions <- getBM(attributes = c("transcription_start_site", "chromosome_name", "transcript_start", "transcript_end","strand",  "ensembl_gene_id","ensembl_transcript_id", "external_gene_name"), mart = ensembl)
positions <- getBM(attributes = c("chromosome_name", "start_position" , "end_position","strand",  "ensembl_gene_id","external_gene_name"), mart = ensembl)
tssCoords=positions

tssCoords[,1]=paste0("chr",tssCoords[,1])
ind=tssCoords[,"strand"]==1
tssCoords[ind,3]=tssCoords[ind,2]+1
ind=tssCoords[,"strand"]==-1
tssCoords[ind,2]=tssCoords[ind,3]-1
key=paste(tssCoords[,1],tssCoords[,2],sep="_")
tssCoords=tssCoords[!duplicated(key),]
d=sapply(tssCoords[,1],nchar)
tssCoords=tssCoords[d<=5,]
names(tssCoords)[c(1:4)]=c("chrom","chromstart","chromend","name")
tssCoords[,4]=seq(1,nrow(tssCoords))
write.table(positions,"refData/hg19.ENSEMBL_genes_biomart.txt",sep="\t",quote=F,row.names=F)
write.table(tssCoords,"refData/hg19.ENSEMBL_geneTSS_biomart_correct.txt",sep="\t",quote=F,row.names=F)
