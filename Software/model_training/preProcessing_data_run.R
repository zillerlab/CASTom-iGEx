#### code written by Lucia Trastulla, e-mail: lucia_trastulla@psych.mpg.de ####

options(stringsAsFactors=F)
options(max.print=1000)
suppressPackageStartupMessages(library(Matrix))
suppressPackageStartupMessages(library(argparse))
suppressPackageStartupMessages(library(biomaRt))


parser <- ArgumentParser(description="preProcessing of input data")

parser$add_argument("--geneExp_file", type = "character", help = "gene expression, complete path")
parser$add_argument("--geneList_file", type = "character", default = NULL, help = "file with name genes heritable, complete path")
parser$add_argument("--VarInfo_file", type = "character", help = "snp annotation info, complete path")
parser$add_argument("--cis_thres", type = "integer", default = 200000, help = "window (in bp) to compute distance from gene and snps")
parser$add_argument("--biomartGenePos_file", type = "character", default = NULL, help = "gene position info file, complete path, if NA recomputed from biomart")
parser$add_argument("--biomartTSS_file", type = "character", default = NULL, help = "TSS info file, complete path, if NA recomputed from biomart")
parser$add_argument("--outFold_geneExp", type = "character", help = "output folder (only for gene expression)")
parser$add_argument("--outFold_snps", type = "character", default = NULL, help = "output folder for SNP info, if NA set as outFold")
parser$add_argument("--outFold", type = "character", help = "output folder")

args <- parser$parse_args()
geneExp_file <- args$geneExp_file
geneList_file <- args$geneList_file
VarInfo_file <- args$VarInfo_file
cis_thres <- args$cis_thres
biomartGenePos_file <- args$biomartGenePos_file
biomartTSS_file <- args$biomartTSS_file
outFold <- args$outFold
outFold_snps <- args$outFold_snps
outFold_geneExp <- args$outFold_geneExp

##############################################################################################################
# geneExp_file <- '/ziller/Michael/CommonMind/SCZ/RNA-Seq_normalized/Gene/EXCLUDE\ ANCESTRY\ +\ SVA/CMC_MSSM-Penn-Pitt_DLPFC_mRNA_IlluminaHiSeq2500_gene-adjustedSVA-dataNormalization-noAncestry-adjustedLogCPM.tsv'
# geneList_file <- '/ziller/lucia/eQTL_PROJECT_CMC/INPUT_DATA_SCRIPTS_v1/list_heritableGenes.txt'
# VarInfo_file <- '/ziller/lucia/eQTL_PROJECT_CMC/INPUT_DATA/Genotyping_data/Genotype_VariantsInfo_CMC-PGC_'
# cis_thres <- 200000
# biomartGenePos_file <- '/ziller/lucia/refData/hg19.ENSEMBL_genes_biomart.txt'
# biomartTSS_file <- '/ziller/lucia/refData/hg19.ENSEMBL_geneTSS_biomart_correct.txt'
# outFold <- '/ziller/lucia/eQTL_PROJECT_CMC/OUTPUT_CMC_SCRIPTS_v1/'
# outFold_geneExp <- '/ziller/lucia/eQTL_PROJECT_CMC/INPUT_DATA_SCRIPTS_v1/RNAseq_data/'
##############################################################################################################

if(is.null(outFold_snps)){outFold_snps <- outFold}
print(outFold_snps)
print(outFold)


# compute TSS for all genes
if(is.null(biomartGenePos_file)){
  
  ensembl = useEnsembl(biomart="ensembl", dataset="hsapiens_gene_ensembl", GRCh=37)
  
  #positions <- getBM(attributes = c("transcription_start_site", "chromosome_name", "transcript_start", "transcript_end","strand",  "ensembl_gene_id","ensembl_transcript_id", "external_gene_name"), mart = ensembl)
  biomart_ann <- getBM(attributes = c("chromosome_name", "start_position" , "end_position","strand",  "ensembl_gene_id","external_gene_name"), mart = ensembl)
  biomart_annTSS=biomart_ann
  
  biomart_annTSS[,1]=paste0("chr",biomart_annTSS[,1])
  ind=biomart_annTSS[,"strand"]==1
  biomart_annTSS[ind,3]=biomart_annTSS[ind,2]+1
  ind=biomart_annTSS[,"strand"]==-1
  biomart_annTSS[ind,2]=biomart_annTSS[ind,3]-1
  key=paste(biomart_annTSS[,1],biomart_annTSS[,2],sep="_")
  biomart_annTSS=biomart_annTSS[!duplicated(key),]
  d=sapply(biomart_annTSS[,1],nchar)
  biomart_annTSS=biomart_annTSS[d<=5,]
  names(biomart_annTSS)[c(1:4)]=c("chrom","chromstart","chromend","name")
  biomart_annTSS[,4]=seq(1,nrow(biomart_annTSS))
  
  write.table(biomart_ann,sprintf("%shg19.ENSEMBL_genes_biomart.txt", outFold),sep="\t",quote=F,row.names=F)
  write.table(biomart_annTSS,sprintf("%shg19.ENSEMBL_geneTSS_biomart_correct.txt", outFold),sep="\t",quote=F,row.names=F)
  
}else{
  
  biomart_ann <- read.table(biomartGenePos_file, header = T, stringsAsFactors = F, sep = '\t')
  biomart_annTSS <- read.table(biomartTSS_file, header = T, stringsAsFactors = F, sep = '\t')
  
}


biomart_ann <- biomart_ann[biomart_ann$ensembl_gene_id %in% biomart_annTSS$ensembl_gene_id,]
if(!identical(biomart_ann$ensembl_gene_id, biomart_annTSS$ensembl_gene_id)){
  print('put biomart annotation in the same order')
  id <- sapply(biomart_annTSS$ensembl_gene_id, function(x) which(biomart_ann$ensembl_gene_id == x))
  biomart_ann <- biomart_ann[id, ]
}

# combine in a unique file
biomart_ann <- cbind(biomart_annTSS[,1:4], biomart_ann[,2:3], biomart_annTSS[,5:6])
colnames(biomart_ann)[2:3] <- c('TSS_start', 'TSS_end')

#### load gene expression file
# load rnaseq data
expData <- read.table(geneExp_file,header = T, stringsAsFactors = F, sep = '\t', check.names=FALSE)
if(any(biomart_ann$ensembl_gene_id %in% expData[, 1])){
  expInfo_tot <- biomart_ann[biomart_ann$ensembl_gene_id %in% expData[, 1], ]  
}else{
  expInfo_tot <- biomart_ann[biomart_ann$external_gene_name %in% expData[, 1], ]
}

# order based on chr and position 
tmp <- list()
for(i in paste0('chr', 1:22)){
  
  tmp[[i]] <- expInfo_tot[expInfo_tot$chrom == i,] 
  tmp[[i]] <-  tmp[[i]][order(tmp[[i]]$start_position),]
  
}

expInfo_tot <- do.call(rbind,tmp)
if(any(expInfo_tot$ensembl_gene_id %in% expData[, 1])){
  id <- sapply(expInfo_tot$ensembl_gene_id, function(x) which(x == expData[, 1]))
}else{
  id <- sapply(expInfo_tot$external_gene_name, function(x) which(x == expData[, 1]))
}

expData_filt <- expData[id, ]
expData_filt <- cbind(expInfo_tot, expData_filt[,-1])

#### load heritable genes list
# match with expData_filt and add a column indicating if belongs to the list
expData_filt <- cbind(data.frame(type = rep('not_heritable', nrow(expData_filt))), expData_filt)

if(!is.null(geneList_file)){
  geneList <- read.table(geneList_file, stringsAsFactors = F, header = T, sep = '\t')
  if(any(expData_filt$ensembl_gene_id %in% geneList[, 1])){
    expData_filt$type[expData_filt$ensembl_gene_id %in% geneList[, 1]] <- 'heritable'
  }else{
    expData_filt$type[expData_filt$external_gene_name %in% geneList[, 1]] <- 'heritable'
  }
}

# save results
write.table(expData_filt, file = sprintf('%sRNAseq_filt.txt', outFold_geneExp), col.names = T, row.names = F, quote = F, sep = '\t')


############################################################################
#### annotate for each chr separately, create gene-SNP matrix distance #####
############################################################################

for(i in 1:22){
  
  curChrom = paste0('chr', i)
  print(curChrom)
  
  snpTab <- read.table(sprintf('%s%s.txt', VarInfo_file, curChrom), header = T, stringsAsFactors = F)
  id <- which(sapply(colnames(snpTab), function(x) strsplit(x, split = 'ID')[[1]][1] == ''))
  tmp <- snpTab[,which(colnames(snpTab) %in% c('CHR', 'POS'))]
  tmp <- cbind(tmp, snpTab[, id])
  colnames(tmp) <- c("chrom","position",names(id))
  tmp$chrom <- sprintf('chr%s',tmp$chrom)
  curSnps <- tmp
  
  curProm <- expData_filt[expData_filt$chrom==curChrom,]
  
  dVals <- t(sapply(1:nrow(curProm),function(X){
    
    v <- which(abs(curSnps$position-curProm$TSS_start[X])<=cis_thres)
    
    if (length(v)==0){
      return(c(NA,NA,NA))
    }else{
      
      r <- cbind(v,X,abs(curSnps$position-curProm$TSS_start[X])[v])
      # correct if the promoter start and the snp start are the same
      r[which(r[,3]==0), 3] <- 1
      return(r)
    }
  }
  ))
  
  dVals <- do.call(rbind,dVals)
  # rm NAs
  dVals <- dVals[!is.na(dVals[,1]),]
  resMat <- sparseMatrix(dVals[,1],dVals[,2],x=dVals[,3],dims=c(nrow(curSnps),nrow(curProm)))
  
  # save
  curProm <- curProm[, colnames(curProm) %in% c('type',	'chrom', 'TSS_start','TSS_end','name','start_position','end_position','ensembl_gene_id','external_gene_name')]
  write.table(curProm,paste0(outFold, "hg19_ENSEMBL_TSS_",curChrom,"_matched.txt"),sep="\t", quote=F, row.names=F, col.names = T)
  write.table(curSnps,paste0(outFold_snps, "hg19_SNPs_", curChrom, "_matched.txt"),sep="\t",quote=F,row.names=F, col.names = T)
  writeMM(resMat,paste0(outFold, "ENSEMBL_gene_SNP_", cis_thres, "_", curChrom,"_matrix.mtx"))
  
}



