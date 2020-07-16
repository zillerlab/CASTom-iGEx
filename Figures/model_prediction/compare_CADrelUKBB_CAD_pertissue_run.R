# find genes and pathway that are significant both for CAD-UKBB and UKBB CAD related phenotypes
options(stringsAsFactors=F)
options(max.print=1000)
suppressPackageStartupMessages(library(argparse))
suppressPackageStartupMessages(library(ggplot2))
suppressPackageStartupMessages(library(ggsci))
suppressPackageStartupMessages(library(pheatmap))
suppressPackageStartupMessages(library(RColorBrewer))
suppressPackageStartupMessages(library(circlize))
suppressPackageStartupMessages(library(ggpubr))
suppressPackageStartupMessages(library(apcluster))
suppressPackageStartupMessages(library(qvalue))
suppressPackageStartupMessages(library(cowplot))
suppressPackageStartupMessages(library(ggsignif))
suppressPackageStartupMessages(library(ggpubr))
options(bitmapType = 'cairo', device = 'png')

parser <- ArgumentParser(description="find enrichment of CAD-UKBB associated element with UKBB gene-pheno or path-pheno phenotypes")
parser$add_argument("--phenoFold", type = "character", help = "Folder with results phenotype annotation UKBB")
parser$add_argument("--inputFold_CADrel", type = "character", help = "Folder with results from CAD-UKBB")
parser$add_argument("--inputFile_CADUKBB", type = "character", help = "File with results from UKBB cad related pheno")
parser$add_argument("--pval_FDR_CADrel", type = "double", default = 0.05, help = "pval threshold to filter the genes and pathways (after BH correction) in CAD related phenotypes")
parser$add_argument("--pval_FDR_CAD", type = "double", default = 0.05,  help = "pval threshold to filter the genes and pathways (after BH correction) in CAD phenotypes")
parser$add_argument("--tissue_name", type = "character", help = "tissue considered")
parser$add_argument("--pval_id", type = "integer", default = 1, help = "1: CAD_HARD, 2:CAD_SOFT") # modify code to make dependent on pval_id
parser$add_argument("--outFold", type="character", help = "Output file [basename only]")

args <- parser$parse_args()
phenoFold <- args$phenoFold
tissue_name <- args$tissue_name
inputFold_CADrel <- args$inputFold_CADrel
inputFile_CADUKBB <- args$inputFile_CADUKBB
pval_FDR_CADrel <- args$pval_FDR_CADrel
pval_FDR_CAD <- args$pval_FDR_CAD
pval_id <- args$pval_id
outFold <- args$outFold

# #########################################################################################################################
# phenoFold <- '/psycl/g/mpsziller/lucia/CAD_UKBB/eQTL_PROJECT/INPUT_DATA_GTEx/CAD/Covariates/UKBB/'
# tissue_name <- 'Liver'
# inputFold_CADrel <- '/psycl/g/mpsziller/lucia/CAD_UKBB/eQTL_PROJECT/OUTPUT_GTEx/predict_CAD/Adipose_Subcutaneous/200kb/CAD_GWAS_bin5e-2/UKBB/devgeno0.01_testdevgeno0/'
# inputFile_CADUKBB <- '/psycl/g/mpsziller/lucia/CAD_UKBB/eQTL_PROJECT/OUTPUT_GTEx/predict_CAD/Adipose_Subcutaneous/200kb/CAD_GWAS_bin5e-2/UKBB/devgeno0.01_testdevgeno0/pval_CAD_pheno_covCorr.RData'
# outFold <- '/psycl/g/mpsziller/lucia/CAD_UKBB/eQTL_PROJECT/OUTPUT_GTEx/predict_CAD/Adipose_Subcutaneous/200kb/CAD_GWAS_bin5e-2/UKBB/devgeno0.01_testdevgeno0/'
# pval_FDR_CADrel <- 0.05
# pval_FDR_CAD <- 0.05
# pval_id <- 1
# ##########################################################################################################################

### load PGC res and filter significant ###
res_pval <- get(load(inputFile_CADUKBB))
tscore_CAD <- res_pval$tscore[[pval_id]]
pathR_CAD <- res_pval$pathScore_reactome[[pval_id]]
pathGO_CAD <- res_pval$pathScore_GO[[pval_id]]

# filter out pathways with only 1 gene and recompute adjusted pvalue
pathR_CAD <- pathR_CAD[pathR_CAD$ngenes_tscore>1, ]
pathR_CAD[, 14] <- qvalue(pathR_CAD[,13])$qvalue
pathR_CAD[, 15] <- p.adjust(p = pathR_CAD[,13], method = 'BH')

pathGO_CAD <- pathGO_CAD[pathGO_CAD$ngenes_tscore>1, ]
pathGO_CAD[, 16] <- qvalue(pathGO_CAD[,15])$qvalue
pathGO_CAD[, 17] <- p.adjust(p = pathGO_CAD[,15], method = 'BH')

# filter based on BHcorr
tscore_CAD_red <- tscore_CAD[tscore_CAD[, 10] <= pval_FDR_CAD,]
pathR_CAD_red <- pathR_CAD[pathR_CAD[, 15] <= pval_FDR_CAD,]
pathGO_CAD_red <- pathGO_CAD[pathGO_CAD[, 17] <= pval_FDR_CAD,]

### load UKBB res and intersect with PGC ###
pheno_name <- c(read.table(sprintf('%smatch_cov_pheno_CADrel_filter.txt', phenoFold), h=F, stringsAsFactors = F)$V1, 'ICD10_Anaemia', 'ICD10_Circulatory_system', 'ICD10_Endocrine', 'ICD10_Respiratory_system')
pheno_name <- pheno_name[!pheno_name %in% c('Medication', 'Medical_conditions')]
pheno_input <- read.delim(sprintf('%sphenotypeDescription_PHESANTproc_CADrelatedpheno_annotated.txt', phenoFold), h=T, stringsAsFactors = F, sep = '\t')

pheno_info <- list()
df_zscore_tscore <- df_zscore_pathGO <- df_zscore_pathR <- list()
df_pval_tscore <- df_pval_pathGO <- df_pval_pathR <- list()
df_pvalcorr_tscore <- df_pvalcorr_pathGO <- df_pvalcorr_pathR <- list()

# compute hypergeometric test to see if there is an enrichment
test_tscore <- test_pathR <- test_pathGO <- list()

for(j in 1:length(pheno_name)){
  
  print(pheno_name[j])
  
  tmp <- get(load(sprintf('%s/pval_%s_pheno_covCorr.RData', inputFold_CADrel, pheno_name[j])))
  rm(final)
  
  id_keep <- 1:nrow(tmp$pheno)
 
  pheno_info[[j]] <- data.frame(pheno = tmp$pheno$pheno_id[id_keep], pheno_type = pheno_name[j], Field = tmp$pheno$Field[id_keep], meaning = tmp$pheno$Coding_meaning[id_keep])
  test_tscore[[j]] <- test_pathR[[j]] <- test_pathGO[[j]] <- data.frame(N = rep(0, length(id_keep)),
                                                                        K =rep(0, length(id_keep)), n =rep(0, length(id_keep)), k= rep(0, length(id_keep)), 
                                                                        fisher_pval = rep(NA, length(id_keep)), fisher_OR = rep(NA, length(id_keep)))
  
  # remove pathways with only 1 gene
  tmp$pathScore_reactome <- lapply(tmp$pathScore_reactome, function(x) x[x$ngenes_tscore > 1,])
  for(id in 1:length(tmp$pathScore_reactome)){
    tmp$pathScore_reactome[[id]][,14] <- qvalue(tmp$pathScore_reactome[[id]][,13])$qvalue
    tmp$pathScore_reactome[[id]][,15] <- p.adjust(tmp$pathScore_reactome[[id]][,13], method = 'BH')
  }
  
  tmp$pathScore_GO <- lapply(tmp$pathScore_GO, function(x) x[x$ngenes_tscore > 1,])
  for(id in 1:length(tmp$pathScore_GO)){
    tmp$pathScore_GO[[id]][,16] <- qvalue(tmp$pathScore_GO[[id]][,15])$qvalue
    tmp$pathScore_GO[[id]][,17] <- p.adjust(tmp$pathScore_GO[[id]][,15], method = 'BH')
  }
  
  # intersect with CAD
  tmp$tscore <- lapply(tmp$tscore[id_keep], function(x) x[x$ensembl_gene_id %in% tscore_CAD$ensembl_gene_id,])
  tmp$pathScore_reactome <- lapply(tmp$pathScore_reactome[id_keep], function(x) x[x$path %in% pathR_CAD$path,])
  tmp$pathScore_GO <- lapply(tmp$pathScore_GO[id_keep], function(x) x[x$path %in% pathGO_CAD$path,])
  
  test_tscore[[j]]$N <- sapply(tmp$tscore, nrow)
  test_pathR[[j]]$N <- sapply(tmp$pathScore_reactome, nrow)
  test_pathGO[[j]]$N <- sapply(tmp$pathScore_GO, nrow)
  
  test_tscore[[j]]$K <- sapply(tmp$tscore, function(x) nrow(tscore_CAD_red[tscore_CAD_red$ensembl_gene_id %in% x$ensembl_gene_id,]))
  test_pathR[[j]]$K <- sapply(tmp$pathScore_reactome, function(x) nrow(pathR_CAD_red[pathR_CAD_red$path %in% x$path,]))
  test_pathGO[[j]]$K <- sapply(tmp$pathScore_GO, function(x) nrow(pathGO_CAD_red[pathGO_CAD_red$path %in% x$path,]))
  
  test_tscore[[j]]$n <- sapply(tmp$tscore, function(x) length(which(x[,10]<=pval_FDR_CADrel)))
  test_pathR[[j]]$n <- sapply(tmp$pathScore_reactome, function(x) length(which(x[,15]<=pval_FDR_CADrel)))
  test_pathGO[[j]]$n <- sapply(tmp$pathScore_GO, function(x) length(which(x[,17]<=pval_FDR_CADrel)))
  
  test_tscore[[j]]$k <- sapply(tmp$tscore, function(x) nrow(x[x$ensembl_gene_id %in% tscore_CAD_red$ensembl_gene_id & x[, 10] <= pval_FDR_CADrel,]))
  test_pathR[[j]]$k <- sapply(tmp$pathScore_reactome, function(x) nrow(x[x$path %in% pathR_CAD_red$path & x[, 15] <= pval_FDR_CADrel,]))
  test_pathGO[[j]]$k <- sapply(tmp$pathScore_GO, function(x) nrow(x[x$path %in% pathGO_CAD_red$path & x[, 17] <= pval_FDR_CADrel,]))
  
  vect_CAD <- rep(0,nrow(tmp$tscore[[1]]))
  vect_CAD[which(tmp$tscore[[1]]$ensembl_gene_id %in% tscore_CAD_red$ensembl_gene_id)] <-1
  vect_UKBB <- lapply(tmp$tscore, function(x) rep(0,nrow(x)))
  for(k in 1:length(tmp$tscore)){
    vect_UKBB[[k]][which(tmp$tscore[[k]][,10]<=pval_FDR_CADrel)] <-1 
  }
  # compute fisher test
  id_notnull <- which(sapply(vect_UKBB, function(x) any(x == 1)))
  if(length(id_notnull)>0 & any(vect_CAD == 1)){
    test_tscore[[j]]$fisher_pval[id_notnull] <- sapply(vect_UKBB[id_notnull], function(x) fisher.test(x = vect_CAD, y =x, alternative = 'greater')$p.value)
    test_tscore[[j]]$fisher_OR[id_notnull] <- sapply(vect_UKBB[id_notnull], function(x) fisher.test(x = vect_CAD, y =x, alternative = 'greater')$estimate)
  }
  
  vect_CAD <- rep(0,nrow(tmp$pathScore_reactome[[1]]))
  vect_CAD[which(tmp$pathScore_reactome[[1]]$path %in% pathR_CAD_red$path)] <-1
  vect_UKBB <- lapply(tmp$pathScore_reactome, function(x) rep(0,nrow(x)))
  for(k in 1:length(tmp$pathScore_reactome)){
    vect_UKBB[[k]][which(tmp$pathScore_reactome[[k]][,15]<=pval_FDR_CADrel)] <-1 
  }
  # compute fisher test
  id_notnull <- which(sapply(vect_UKBB, function(x) any(x == 1)))
  if(length(id_notnull) > 0 & any(vect_CAD == 1)){
    test_pathR[[j]]$fisher_pval[id_notnull]<- sapply(vect_UKBB[id_notnull], function(x) fisher.test(x = vect_CAD, y =x, alternative = 'greater')$p.value)
    test_pathR[[j]]$fisher_OR[id_notnull]<- sapply(vect_UKBB[id_notnull], function(x) fisher.test(x = vect_CAD, y =x, alternative = 'greater')$estimate)
  }
  
  vect_CAD <- rep(0,nrow(tmp$pathScore_GO[[1]]))
  vect_CAD[which(tmp$pathScore_GO[[1]]$path %in% pathGO_CAD_red$path)] <-1
  vect_UKBB <- lapply(tmp$pathScore_GO, function(x) rep(0,nrow(x)))
  for(k in 1:length(tmp$pathScore_GO)){
    vect_UKBB[[k]][which(tmp$pathScore_GO[[k]][,17]<=pval_FDR_CADrel)] <-1 
  }
  # compute fisher test
  id_notnull <- which(sapply(vect_UKBB, function(x) any(x == 1)))
  if(length(id_notnull) > 0 & any(vect_CAD == 1)){
    test_pathGO[[j]]$fisher_pval[id_notnull]<- sapply(vect_UKBB[id_notnull], function(x) fisher.test(x = vect_CAD, y =x, alternative = 'greater')$p.value)
    test_pathGO[[j]]$fisher_OR[id_notnull]<- sapply(vect_UKBB[id_notnull], function(x) fisher.test(x = vect_CAD, y =x, alternative = 'greater')$estimate)
  }
  
  
  # keep only significant elements in CAD
  tmp$tscore <- lapply(tmp$tscore, function(x) x[x$ensembl_gene_id %in% tscore_CAD_red$ensembl_gene_id,])
  tmp$pathScore_reactome <- lapply(tmp$pathScore_reactome, function(x) x[x$path %in% pathR_CAD_red$path,])
  tmp$pathScore_GO <- lapply(tmp$pathScore_GO, function(x) x[x$path %in% pathGO_CAD_red$path,])
  
  # tscore
  df_zscore_tscore[[j]] <- cbind(data.frame(ensembl_gene_id =tmp$tscore[[1]]$ensembl_gene_id, 
                                            external_gene_name = tmp$tscore[[1]]$external_gene_name, stringsAsFactors = F),
                                 sapply(tmp$tscore, function(x) x[,7]))
  colnames(df_zscore_tscore[[j]])[-c(1:2)] <- pheno_info[[j]]$pheno
  
  df_pval_tscore[[j]] <- cbind(data.frame(ensembl_gene_id =tmp$tscore[[1]]$ensembl_gene_id, 
                                          external_gene_name = tmp$tscore[[1]]$external_gene_name, stringsAsFactors = F),
                               sapply(tmp$tscore, function(x) x[,8]))
  colnames(df_pval_tscore[[j]])[-c(1:2)] <- pheno_info[[j]]$pheno
  
  df_pvalcorr_tscore[[j]] <- cbind(data.frame(ensembl_gene_id =tmp$tscore[[1]]$ensembl_gene_id, 
                                              external_gene_name = tmp$tscore[[1]]$external_gene_name, stringsAsFactors = F),
                                   sapply(tmp$tscore, function(x) x[,10]))
  colnames(df_pvalcorr_tscore[[j]])[-c(1:2)] <- pheno_info[[j]]$pheno
  
  # pathR
  df_zscore_pathR[[j]] <- cbind(data.frame(path =tmp$pathScore_reactome[[1]]$path, stringsAsFactors = F),
                                sapply(tmp$pathScore_reactome, function(x) x[,12]))
  colnames(df_zscore_pathR[[j]])[-1] <- pheno_info[[j]]$pheno
  
  df_pval_pathR[[j]] <- cbind(data.frame(path =tmp$pathScore_reactome[[1]]$path,stringsAsFactors = F),
                              sapply(tmp$pathScore_reactome, function(x) x[,13]))
  colnames(df_pval_pathR[[j]])[-1] <- pheno_info[[j]]$pheno
  
  df_pvalcorr_pathR[[j]] <- cbind(data.frame(path =tmp$pathScore_reactome[[1]]$path,stringsAsFactors = F),
                                  sapply(tmp$pathScore_reactome, function(x) x[,15]))
  colnames(df_pvalcorr_pathR[[j]])[-1] <- pheno_info[[j]]$pheno
  
  # pathGO
  df_zscore_pathGO[[j]] <- cbind(data.frame(path =tmp$pathScore_GO[[1]]$path, stringsAsFactors = F),
                                 sapply(tmp$pathScore_GO, function(x) x[,14]))
  colnames(df_zscore_pathGO[[j]])[-1] <- pheno_info[[j]]$pheno
  
  df_pval_pathGO[[j]] <- cbind(data.frame(path =tmp$pathScore_GO[[1]]$path,stringsAsFactors = F),
                               sapply(tmp$pathScore_GO, function(x) x[,15]))
  colnames(df_pval_pathGO[[j]])[-1] <- pheno_info[[j]]$pheno
  
  df_pvalcorr_pathGO[[j]] <- cbind(data.frame(path =tmp$pathScore_GO[[1]]$path,stringsAsFactors = F),
                                   sapply(tmp$pathScore_GO, function(x) x[,17]))
  colnames(df_pvalcorr_pathGO[[j]])[-1] <- pheno_info[[j]]$pheno
  
}

pheno_info <- do.call(rbind, pheno_info)
pheno_info$names_field <- NA
pheno_info$names_field[is.na(pheno_info$meaning)] <- paste(pheno_info$Field, paste0('(',pheno_info$pheno,')'))[is.na(pheno_info$meaning)]
pheno_info$names_field[!is.na(pheno_info$meaning)] <- paste(pheno_info$Field, ':', pheno_info$meaning, paste0('(',pheno_info$pheno,')'))[!is.na(pheno_info$meaning)]

test_tscore <- do.call(rbind, test_tscore)
test_pathR <- do.call(rbind, test_pathR)
test_pathGO <- do.call(rbind, test_pathGO)
# correct pvalues
if(any(!is.na(test_tscore$fisher_pval))){
  test_tscore$fisher_qval <- qvalue(test_tscore$fisher_pval)$qvalue
  test_tscore$fisher_pval_BHcorr <- p.adjust(test_tscore$fisher_pval, method = 'BH')
}else{
  test_tscore$fisher_qval <- NA
  test_tscore$fisher_pval_BHcorr <- NA
}

if(any(!is.na(test_pathR$fisher_pval))){
  test_pathR$fisher_qval <- qvalue(test_pathR$fisher_pval)$qvalue
  test_pathR$fisher_pval_BHcorr <- p.adjust(test_pathR$fisher_pval, method = 'BH')
}else{
  test_pathR$fisher_qval <- NA
  test_pathR$fisher_pval_BHcorr <- NA
}

if(any(!is.na(test_pathGO$fisher_pval))){
  test_pathGO$fisher_qval <- qvalue(test_pathGO$fisher_pval)$qvalue
  test_pathGO$fisher_pval_BHcorr <- p.adjust(test_pathGO$fisher_pval, method = 'BH')
}else{
  test_pathGO$fisher_qval <- NA
  test_pathGO$fisher_pval_BHcorr <- NA
}

# add pheno name
test_tscore$names_field <- pheno_info$names_field
test_pathR$names_field <- pheno_info$names_field
test_pathGO$names_field <- pheno_info$names_field
test_tscore$pheno <- pheno_info$pheno
test_pathR$pheno <- pheno_info$pheno
test_pathGO$pheno <- pheno_info$pheno
test_tscore$pheno_type <- pheno_info$pheno_type
test_pathR$pheno_type <- pheno_info$pheno_type
test_pathGO$pheno_type <- pheno_info$pheno_type

# save test results (to be used together for all the tissues)
res_test_enrichment <- list(tscore = test_tscore, pathScore_reactome = test_pathR, pathScore_GO = test_pathGO)
save(res_test_enrichment, file = sprintf('%s/%s_enrichment_significant_pvalFDR-CAD%.2f_pvalFDR-CADrel%.2f.RData', outFold, tissue_name, pval_FDR_CAD, pval_FDR_CADrel))

# create unique matrices for all the phenotypes
# tscore
genes_int <- df_zscore_tscore[[1]][,c(1:2)] # genes are always the same!
df_zscore_tscore_all <- data.frame(genes_int)
df_pval_tscore_all <- data.frame(genes_int)
df_pvalcorr_tscore_all <- data.frame(genes_int)

if(length(genes_int)>0){
for(j in 1:length(df_zscore_tscore)){
    df_zscore_tscore_all <- cbind(df_zscore_tscore_all, df_zscore_tscore[[j]][, -c(1:2)])
    df_pval_tscore_all <- cbind(df_pval_tscore_all, df_pval_tscore[[j]][, -c(1:2)])
    df_pvalcorr_tscore_all <- cbind(df_pvalcorr_tscore_all, df_pvalcorr_tscore[[j]][, -c(1:2)])
  }
  colnames(df_zscore_tscore_all)[-c(1:2)] <- pheno_info$pheno
  colnames(df_pval_tscore_all)[-c(1:2)] <- pheno_info$pheno
  colnames(df_pvalcorr_tscore_all)[-c(1:2)] <- pheno_info$pheno
  ##### filter out based on pvalcorr #####
  id_notsign <- df_pvalcorr_tscore_all[,-c(1:2)] > pval_FDR_CADrel
  df_zscore_tscore_all[,-c(1:2)][id_notsign] <- 0
  id_keep <- list(row = rowSums(df_zscore_tscore_all[, -c(1:2)], na.rm = T) != 0, col = colSums(df_zscore_tscore_all[, -c(1:2)],  na.rm = T) != 0)
  df_zscore_tscore_all <- df_zscore_tscore_all[id_keep$row,c(T,T, id_keep$col)]
}
# pathR
path_int <- df_zscore_pathR[[1]][,1] # path are always the same!
df_zscore_pathR_all <- data.frame(path  = path_int)
df_pval_pathR_all <- data.frame(path  = path_int)
df_pvalcorr_pathR_all <- data.frame(path  = path_int)
if(length(path_int)>0){
  for(j in 1:length(df_zscore_pathR)){
    df_zscore_pathR_all <- cbind(df_zscore_pathR_all, df_zscore_pathR[[j]][, -1])
    df_pval_pathR_all <- cbind(df_pval_pathR_all, df_pval_pathR[[j]][, -1])
    df_pvalcorr_pathR_all <- cbind(df_pvalcorr_pathR_all, df_pvalcorr_pathR[[j]][, -1])
  }
  colnames(df_zscore_pathR_all)[-1] <- pheno_info$pheno
  colnames(df_pval_pathR_all)[-1] <- pheno_info$pheno
  colnames(df_pvalcorr_pathR_all)[-1] <- pheno_info$pheno
  ##### filter out based on pvalcorr #####
  id_notsign <- df_pvalcorr_pathR_all[,-c(1)] > pval_FDR_CADrel
  df_zscore_pathR_all[,-c(1)][id_notsign] <- 0
  id_keep <- list(row = rowSums(df_zscore_pathR_all[, -c(1)], na.rm = T) != 0, col = colSums(df_zscore_pathR_all[, -c(1)],  na.rm = T) != 0)
  df_zscore_pathR_all <- df_zscore_pathR_all[id_keep$row,c(T, id_keep$col)]
}


# pathGO
path_int <- df_zscore_pathGO[[1]][,1] # path are always the same!
df_zscore_pathGO_all <- data.frame(path  = path_int)
df_pval_pathGO_all <- data.frame(path  = path_int)
df_pvalcorr_pathGO_all <- data.frame(path  = path_int)
if(length(path_int)>0){
  for(j in 1:length(df_zscore_pathGO)){
    df_zscore_pathGO_all <- cbind(df_zscore_pathGO_all, df_zscore_pathGO[[j]][, -1])
    df_pval_pathGO_all <- cbind(df_pval_pathGO_all, df_pval_pathGO[[j]][, -1])
    df_pvalcorr_pathGO_all <- cbind(df_pvalcorr_pathGO_all, df_pvalcorr_pathGO[[j]][, -1])
  }
  colnames(df_zscore_pathGO_all)[-1] <- pheno_info$pheno
  colnames(df_pval_pathGO_all)[-1] <- pheno_info$pheno
  colnames(df_pvalcorr_pathGO_all)[-1] <- pheno_info$pheno
  ##### filter out based on pvalcorr #####
  id_notsign <- df_pvalcorr_pathGO_all[,-c(1)] > pval_FDR_CADrel
  df_zscore_pathGO_all[,-c(1)][id_notsign] <- 0
  id_keep <- list(row = rowSums(df_zscore_pathGO_all[, -c(1)], na.rm = T) != 0, col = colSums(df_zscore_pathGO_all[, -c(1)], na.rm = T) != 0)
  df_zscore_pathGO_all <- df_zscore_pathGO_all[id_keep$row,c(T, id_keep$col)]
}
res_zscore <- list(tscore = df_zscore_tscore_all, path_Reactome = df_zscore_pathR_all, path_GO = df_zscore_pathGO_all, pheno = pheno_info)
save(res_zscore, file = sprintf("%s/%s_Zstatistic_intersection_pvalFDR-CAD%.2f_pvalFDR-CADrel%.2f.RData", outFold, tissue_name, pval_FDR_CAD, pval_FDR_CADrel))

### make heatmap of zscore for significant elements
save_pheatmap_png <- function(x, filename, width=10, height=10, res = 500) {
  png(filename, width = width, height = height, res = res, units = 'in')
  grid::grid.newpage()
  grid::grid.draw(x$gtable)
  dev.off()
}
save_pheatmap_pdf <- function(x, filename, width=10, height=10) {
  pdf(filename, width = width, height = height)
  grid::grid.newpage()
  grid::grid.draw(x$gtable)
  dev.off()
}
# used to make heatmaps 
draw_colnames_90 <- function (coln, gaps, ...) {
  coord <- pheatmap:::find_coordinates(length(coln), gaps)
  x     <- coord$coord - 0.5 * coord$size
  res   <- grid::textGrob(
    coln, x = x, y = unit(1, "npc") - unit(3,"bigpts"),
    vjust = 0.75, hjust = 1, rot = 90, gp = grid::gpar(...)
  )
  return(res)
}
assignInNamespace(
  x = "draw_colnames",
  value = "draw_colnames_90",
  ns = asNamespace("pheatmap")
)

color <- grDevices::colors()[grep('gr(a|e)y', grDevices::colors(), invert = T)]
set.seed(19)
color_pheno_cat <- sample(color, length(pheno_name))
df_color_pheno <- data.frame(color = color_pheno_cat, pheno_type = pheno_name)
if(tissue_name == 'Adipose_Subcutaneous'){
  write.table(x = df_color_pheno, file = sprintf('%s/color_pheno_type_UKBB.txt', phenoFold), col.names = T, row.names = F, sep = '\t', quote = F)
}

plot_heatmap_split <- function(type_mat, df_zscore, CAD_res, id_to_rm, id_col, id_col_match, pheno_info, color_cat, 
                                    cap_val = 9, height_pl = 16, width_pl=16){
  
  # split blood pheno
  id_blood <- colnames(df_zscore) %in% pheno_info$pheno[pheno_info$pheno_type == 'Blood_count'] 
  df_zscore_blood <- df_zscore[, c(id_to_rm,which(id_blood))]
  df_zscore_blood <- df_zscore_blood[rowSums(df_zscore_blood[,-id_to_rm])!=0, ]
  # split ICD10
  id_icd10 <- colnames(df_zscore) %in% pheno_info$pheno[grepl('ICD10', pheno_info$pheno_type)] 
  df_zscore_icd10 <- df_zscore[, c(id_to_rm,which(id_icd10))]
  df_zscore_icd10 <- df_zscore_icd10[rowSums(df_zscore_icd10[,-id_to_rm])!=0, ]
  
  df_zscore_rest <- df_zscore[, which(!id_blood & !id_icd10)]
  df_zscore_rest <- df_zscore_rest[rowSums(df_zscore_rest[,-id_to_rm])!=0,]
  
  ####################
  #### blood part ####
  ####################
  CAD_res_int <- CAD_res[match(df_zscore_blood[,1], CAD_res[,id_col_match]), ]
  # sort according to zscore and change the order in df_zscore
  CAD_res_int <- CAD_res_int[order(CAD_res_int$CAD_HARD_z_t, decreasing = T),]
  df_zscore_blood <- df_zscore_blood[match(CAD_res_int[,id_col_match],df_zscore_blood[,1]), ]
  
  # remove duplicate names:
  dup_names <- names(which(table(df_zscore_blood[, id_col])>1))
  if(length(dup_names)>0){
    df_zscore_blood <- df_zscore_blood[!df_zscore_blood[, id_col] %in% dup_names,]
    CAD_res_int <- CAD_res_int[!CAD_res_int[, id_col] %in% dup_names,]
  }
  
  coul <- rev(colorRampPalette(brewer.pal(11, "RdBu"))(25))
  tmp_mat <- as.matrix(df_zscore_blood[, -id_to_rm])
  rownames(tmp_mat) <- df_zscore_blood[, id_col]
  val <- min(cap_val,round(max(abs(min(tmp_mat)), max(tmp_mat))))
  mat_breaks <- seq(-val, val, length.out = 25)
  pheno_tmp <- pheno_info[match(colnames(tmp_mat), pheno_info$pheno),] 
  
  # Data frame with column annotations.
  mat_col <- data.frame(pheno = pheno_tmp$pheno_type)
  rownames(mat_col) <- pheno_tmp$names_field
  
  # List with colors for each annotation.
  if(is.null(color_cat)){color_cat <- colorRampPalette(brewer.pal(8, "Set2"))(length(unique(pheno_info$pheno_type)))}
  mat_colors <- list(pheno = color_cat)
  names(mat_colors$pheno) <- unique(pheno_info$pheno_type)
  
  tmp_mat_capped <- tmp_mat
  tmp_mat_capped[tmp_mat>=cap_val] <- cap_val
  tmp_mat_capped[tmp_mat<=-cap_val] <- -cap_val
  colnames(tmp_mat_capped) <- pheno_tmp$names_field
  #id_notcol <- which(colnames(tmp_mat_capped) %in% enrichment_res$names_field[enrichment_res$fisher_pval_BHcorr>0.05])
  ## keep as pheno only the enriched ones
  #colnames(tmp_mat_capped)[id_notcol] <- rep("", length(which(enrichment_res$fisher_pval_BHcorr>0.05)))
  
  # for genes create 2 groups of z>0 and z<0
  cl <- data.frame(id = rownames(tmp_mat_capped), reg = 'up')
  cl$reg[which(CAD_res_int$CAD_HARD_z_t<0)] = 'low'
  
  mat_row <- data.frame(reg_CAD = cl$reg)
  rownames(mat_row) <- rownames(tmp_mat_capped)
  
  mat_colors_gr <- list(reg_CAD = c(coul[length(coul)-1], coul[2]))
  names(mat_colors_gr$reg_CAD) <- c('up', 'low')
  
  new_color = list(pheno = mat_colors[[1]], reg_CAD = mat_colors_gr[[1]])
  
  hm_pl <- pheatmap(mat=tmp_mat_capped, color=coul, show_colnames = T, show_rownames = T, annotation_col = mat_col,  annotation_names_col = F, annotation_names_row = F, 
                    cluster_rows=F, cluster_cols=F, border_color = NA, 
                    annotation_row = mat_row, annotation_colors = new_color, drop_levels = TRUE, fontsize_row = 8, fontsize_col = 8, fontsize = 12, 
                    main = sprintf('significant CAD-HARD and significant pheno\n in %s %s', tissue_name, type_mat),
                    breaks = mat_breaks, treeheight_row = 0, treeheight_col = 0)
  
  
  save_pheatmap_png(hm_pl, sprintf("%s/%s_%s_intersection_Blood_count_pvalFDR-CAD%.2f_pvalFDR-CADrel%.2f.png", outFold, tissue_name, type_mat, pval_FDR_CAD, pval_FDR_CADrel), height = height_pl , width =width_pl-4)
  save_pheatmap_pdf(hm_pl, sprintf("%s/%s_%s_intersection_Blood_count_pvalFDR-CAD%.2f_pvalFDR-CADrel%.2f.pdf", outFold, tissue_name, type_mat, pval_FDR_CAD, pval_FDR_CADrel), height = height_pl , width =width_pl-4)
  
  
  ###################
  #### IC10 part ####
  ###################
  
  CAD_res_int <- CAD_res[match(df_zscore_icd10[,1], CAD_res[,id_col_match]), ]
  # sort according to zscore and change the order in df_zscore
  CAD_res_int <- CAD_res_int[order(CAD_res_int$CAD_HARD_z_t, decreasing = T),]
  df_zscore_icd10 <- df_zscore_icd10[match(CAD_res_int[,id_col_match],df_zscore_icd10[,1]), ]
  
  # remove duplicate names:
  dup_names <- names(which(table(df_zscore_icd10[, id_col])>1))
  if(length(dup_names)>0){
    df_zscore_icd10 <- df_zscore_icd10[!df_zscore_icd10[, id_col] %in% dup_names,]
    CAD_res_int <- CAD_res_int[!CAD_res_int[, id_col] %in% dup_names,]
  }
  
  coul <- rev(colorRampPalette(brewer.pal(11, "RdBu"))(25))
  tmp_mat <- as.matrix(df_zscore_icd10[, -id_to_rm])
  rownames(tmp_mat) <- df_zscore_icd10[, id_col]
  val <- min(cap_val,round(max(abs(min(tmp_mat)), max(tmp_mat))))
  mat_breaks <- seq(-val, val, length.out = 25)
  pheno_tmp <- pheno_info[match(colnames(tmp_mat), pheno_info$pheno),] 
  
  # Data frame with column annotations.
  mat_col <- data.frame(pheno = pheno_tmp$pheno_type)
  rownames(mat_col) <- pheno_tmp$names_field
  
  # List with colors for each annotation.
  if(is.null(color_cat)){color_cat <- colorRampPalette(brewer.pal(8, "Set2"))(length(unique(pheno_info$pheno_type)))}
  mat_colors <- list(pheno = color_cat)
  names(mat_colors$pheno) <- unique(pheno_info$pheno_type)
  
  tmp_mat_capped <- tmp_mat
  tmp_mat_capped[tmp_mat>=cap_val] <- cap_val
  tmp_mat_capped[tmp_mat<=-cap_val] <- -cap_val
  colnames(tmp_mat_capped) <- pheno_tmp$names_field
  #id_notcol <- which(colnames(tmp_mat_capped) %in% enrichment_res$names_field[enrichment_res$fisher_pval_BHcorr>0.05])
  ## keep as pheno only the enriched ones
  #colnames(tmp_mat_capped)[id_notcol] <- rep("", length(which(enrichment_res$fisher_pval_BHcorr>0.05)))
  
  # for genes create 2 groups of z>0 and z<0
  cl <- data.frame(id = rownames(tmp_mat_capped), reg = 'up')
  cl$reg[which(CAD_res_int$CAD_HARD_z_t<0)] = 'low'
  
  mat_row <- data.frame(reg_CAD = cl$reg)
  rownames(mat_row) <- rownames(tmp_mat_capped)
  
  mat_colors_gr <- list(reg_CAD = c(coul[length(coul)-1], coul[2]))
  names(mat_colors_gr$reg_CAD) <- c('up', 'low')
  
  new_color = list(pheno = mat_colors[[1]], reg_CAD = mat_colors_gr[[1]])
  
  hm_pl <- pheatmap(mat=tmp_mat_capped, color=coul, show_colnames = T, show_rownames = T, annotation_col = mat_col,  annotation_names_col = F, annotation_names_row = F, 
                    cluster_rows=F, cluster_cols=F, border_color = NA, 
                    annotation_row = mat_row, annotation_colors = new_color, drop_levels = TRUE, fontsize_row = 8, fontsize_col = 8, fontsize = 12, 
                    main = sprintf('significant CAD-HARD and significant pheno\n in %s %s', tissue_name, type_mat),
                    breaks = mat_breaks, treeheight_row = 0, treeheight_col = 0)
  
  
  save_pheatmap_png(hm_pl, sprintf("%s/%s_%s_intersection_ICD10_pvalFDR-CAD%.2f_pvalFDR-CADrel%.2f.png", outFold, tissue_name, type_mat, pval_FDR_CAD, pval_FDR_CADrel), height = height_pl , width =width_pl-4)
  save_pheatmap_pdf(hm_pl, sprintf("%s/%s_%s_intersection_ICD10_pvalFDR-CAD%.2f_pvalFDR-CADrel%.2f.pdf", outFold, tissue_name, type_mat, pval_FDR_CAD, pval_FDR_CADrel), height = height_pl , width =width_pl-4)
  
  ###################
  #### remaining ####
  ###################
  
  CAD_res_int <- CAD_res[match(df_zscore_rest[,1], CAD_res[,id_col_match]), ]
  # sort according to zscore and change the order in df_zscore
  CAD_res_int <- CAD_res_int[order(CAD_res_int$CAD_HARD_z_t, decreasing = T),]
  df_zscore_rest <- df_zscore_rest[match(CAD_res_int[,id_col_match],df_zscore_rest[,1]), ]
  
  # remove duplicate names:
  dup_names <- names(which(table(df_zscore_rest[, id_col])>1))
  if(length(dup_names)>0){
    df_zscore_rest <- df_zscore_rest[!df_zscore_rest[, id_col] %in% dup_names,]
    CAD_res_int <- CAD_res_int[!CAD_res_int[, id_col] %in% dup_names,]
  }
  
  coul <- rev(colorRampPalette(brewer.pal(11, "RdBu"))(25))
  tmp_mat <- as.matrix(df_zscore_rest[, -id_to_rm])
  rownames(tmp_mat) <- df_zscore_rest[, id_col]
  val <- min(cap_val,round(max(abs(min(tmp_mat)), max(tmp_mat))))
  mat_breaks <- seq(-val, val, length.out = 25)
  pheno_tmp <- pheno_info[match(colnames(tmp_mat), pheno_info$pheno),] 
  
  # Data frame with column annotations.
  mat_col <- data.frame(pheno = pheno_tmp$pheno_type)
  rownames(mat_col) <- pheno_tmp$names_field
  
  # List with colors for each annotation.
  if(is.null(color_cat)){color_cat <- colorRampPalette(brewer.pal(8, "Set2"))(length(unique(pheno_info$pheno_type)))}
  mat_colors <- list(pheno = color_cat)
  names(mat_colors$pheno) <- unique(pheno_info$pheno_type)
  
  tmp_mat_capped <- tmp_mat
  tmp_mat_capped[tmp_mat>=cap_val] <- cap_val
  tmp_mat_capped[tmp_mat<=-cap_val] <- -cap_val
  colnames(tmp_mat_capped) <- pheno_tmp$names_field
  #id_notcol <- which(colnames(tmp_mat_capped) %in% enrichment_res$names_field[enrichment_res$fisher_pval_BHcorr>0.05])
  ## keep as pheno only the enriched ones
  #colnames(tmp_mat_capped)[id_notcol] <- rep("", length(which(enrichment_res$fisher_pval_BHcorr>0.05)))
  
  # for genes create 2 groups of z>0 and z<0
  cl <- data.frame(id = rownames(tmp_mat_capped), reg = 'up')
  cl$reg[which(CAD_res_int$CAD_HARD_z_t<0)] = 'low'
  
  mat_row <- data.frame(reg_CAD = cl$reg)
  rownames(mat_row) <- rownames(tmp_mat_capped)
  
  mat_colors_gr <- list(reg_CAD = c(coul[length(coul)-1], coul[2]))
  names(mat_colors_gr$reg_CAD) <- c('up', 'low')
  
  new_color = list(pheno = mat_colors[[1]], reg_CAD = mat_colors_gr[[1]])
  
  hm_pl <- pheatmap(mat=tmp_mat_capped, color=coul, show_colnames = T, show_rownames = T, annotation_col = mat_col,  annotation_names_col = F, annotation_names_row = F, 
                    cluster_rows=F, cluster_cols=F, border_color = NA, 
                    annotation_row = mat_row, annotation_colors = new_color, drop_levels = TRUE, fontsize_row = 8, fontsize_col = 8, fontsize = 12, 
                    main = sprintf('significant CAD-HARD and significant pheno\n in %s %s', tissue_name, type_mat),
                    breaks = mat_breaks, treeheight_row = 0, treeheight_col = 0)
  
  save_pheatmap_png(hm_pl, sprintf("%s/%s_%s_intersection_rest_pvalFDR-CAD%.2f_pvalFDR-CADrel%.2f.png", outFold, tissue_name, type_mat, pval_FDR_CAD, pval_FDR_CADrel), height = height_pl , width =width_pl-3)
  save_pheatmap_pdf(hm_pl, sprintf("%s/%s_%s_intersection_rest_pvalFDR-CAD%.2f_pvalFDR-CADrel%.2f.pdf", outFold, tissue_name, type_mat, pval_FDR_CAD, pval_FDR_CADrel), height = height_pl , width =width_pl-3)
  
}


# plot heatmap
if(nrow(df_zscore_tscore_all)>0){
  plot_heatmap_split(type_mat = 'tscore', id_col_match = 1, df_zscore = df_zscore_tscore_all, CAD_res = tscore_CAD_red, 
                     id_to_rm = c(1:2), id_col = 2, pheno_info = pheno_info, color_cat = color_pheno_cat)
}
if(nrow(df_zscore_pathR_all)>0){
  plot_heatmap_split(type_mat = 'pathScore_Reactome', id_col_match = 1, df_zscore = df_zscore_pathR_all, CAD_res = pathR_CAD_red, 
                     id_to_rm = c(1), id_col = 1, pheno_info = pheno_info, color_cat = color_pheno_cat, width_pl = 12, height_pl = 7)
}
if(nrow(df_zscore_pathGO_all)>0){
  plot_heatmap_split(type_mat = 'pathScore_GO', id_col_match = 2, df_zscore = df_zscore_pathGO_all, CAD_res = pathGO_CAD_red, 
                     id_to_rm = c(1), id_col = 1, pheno_info = pheno_info, color_cat = color_pheno_cat)
  
}


