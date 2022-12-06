#!/usr/bin/env Rscript

# save also beta, se, zcores
# consider results from big matrix, possibility to load more than 1 cov and pheno dat, add pheno annotation info
# specific for pathway scores

options(stringsAsFactors=F)
options(max.print=1000)
suppressPackageStartupMessages(library(argparse))
suppressPackageStartupMessages(library(bigmemory))
suppressPackageStartupMessages(library(limma))
suppressPackageStartupMessages(library(biomaRt))
suppressPackageStartupMessages(library(nnet))
suppressPackageStartupMessages(library(lattice))
suppressPackageStartupMessages(library(data.table))
suppressPackageStartupMessages(library(MASS))
suppressPackageStartupMessages(library(lmtest))
suppressPackageStartupMessages(library(doParallel))
suppressPackageStartupMessages(library(PGSEA))

parser <- ArgumentParser(description="Pathway association analysis")
parser$add_argument("--inputFile", type = "character", help = "RData to be loaded")
parser$add_argument("--inputInfoFile", type = "character",  help = "RData to be loaded, contains info on genes or path")
parser$add_argument("--sampleAnn_file", type = "character", help = "file with sample info for new genotype data, must contain Dx column (0 controls 1 cases)")
parser$add_argument("--covDat_file", type = "character", nargs = '*', help = "file with sample info for new genotype data, must contain Dx column (0 controls 1 cases), associated to phenoDat_file")
parser$add_argument("--phenoDat_file", type = "character", nargs = '*', help = "file(s) with individual_ID to match and phenotype to test association, associated to covDat_file")
parser$add_argument("--names_file", type = "character", nargs = '*', help = "for each couple of covDat/phenoDat file, associated name")
parser$add_argument("--phenoAnn_file", type = "character", help = "file with phenotype annotation (used to determine the type of regression)")
parser$add_argument("--cov_corr", type = "logical",default = F,  help = "if T column in covDat_file are use to correct association (excluded Dx and IDs)")
parser$add_argument("--functR", type = "character", help = "Rscript with functions to be used")
parser$add_argument("--path_type", type = "character", help = "reactome or GO or custom")
parser$add_argument("--ncores", type = "integer", help = "cores parallelization")
parser$add_argument("--outFile", type="character", help = "Output file [basename only]")

args <- parser$parse_args()
inputFile <- args$inputFile
inputInfoFile <- args$inputInfoFile
covDat_file <- args$covDat_file
sampleAnn_file <- args$sampleAnn_file
phenoDat_file <- args$phenoDat_file
phenoAnn_file <- args$phenoAnn_file
cov_corr <- args$cov_corr
names_file <- args$names_file
functR <- args$functR
path_type <- args$path_type
ncores <- args$ncores
outFile <- args$outFile

# #########################################################################################
# inputFile <- paste0('OUTPUT_GTEx/predict_UKBB/Pancreas/200kb/noGWAS/devgeno0.01_testdevgeno0/T1D_pheno/Pathway_Reactome_scores_splitPath1.RData')
# inputInfoFile <- 'OUTPUT_GTEx/predict_UKBB/Pancreas/200kb/noGWAS/devgeno0.01_testdevgeno0/T1D_pheno/pathScore_Reactome_info.RData'
# covDat_file <- 'INPUT_DATA/Covariates/covariateMatrix_T1D.txt'
# phenoDat_file <- 'INPUT_DATA/Covariates/phenotypeMatrix_T1D.txt'
# names_file <- 'T1D'
# sampleAnn_file <- 'INPUT_DATA/Covariates/covariateMatrix_T1D.txt'
# cov_corr <- T
# phenoAnn_file <- 'INPUT_DATA/Covariates/phenotypeDescription_T1D.txt'
# outFile <- 'OUTPUT_GTEx/predict_UKBB/Pancreas/200kb/noGWAS/devgeno0.01_testdevgeno0/T1D_pheno/'
# functR <- 'RSCRIPTS/SCRIPTS_v2/AssociationAnalysis_functions_run.R'
# path_type <- 'Reactome'
#########################################################################################


source(functR)

# load sample annotation
sampleAnn <- read.table(sampleAnn_file, h=T, stringsAsFactors=F)

sampleAnn$Temp_ID <- sampleAnn$Individual_ID
id_h <- sapply(sampleAnn$Individual_ID, function(x) length(strsplit(x, split = '-')[[1]]) == 2)
id_n <- sapply(sampleAnn$Individual_ID, function(x) !is.na(as.numeric(strsplit(x, split = '')[[1]][1])))
if(any(id_h)){
  sampleAnn$Temp_ID[id_h] <- sapply(sampleAnn$Individual_ID[id_h], function(x) paste0(strsplit(x, split = '-')[[1]], collapse = '_'))
}
if(any(id_n)){
  sampleAnn$Temp_ID[id_n] <- sapply(sampleAnn$Individual_ID[id_n], function(x) paste0('X', x))   
}

df_info <- get(load(inputInfoFile))
score_mat <- get(load(inputFile))

# filter info based on score mat
if(path_type == 'GO'){
  df_info <- df_info[df_info$path_id %in% score_mat$path, ]
  identical(df_info$path_id, score_mat$path)
}else{
  df_info <- df_info[df_info$path %in% score_mat$path, ]
  identical(df_info$path, score_mat$path)
}

pathScoreMat <- matrix(t(score_mat[, -1]), nrow = ncol(score_mat)-1, ncol = nrow(score_mat))
rownames(pathScoreMat) <- colnames(score_mat[, -1])
# match
if(!identical(rownames(pathScoreMat), sampleAnn$Temp_ID)){
  common_s <- intersect(rownames(pathScoreMat), sampleAnn$Temp_ID)
  sampleAnn <- sampleAnn[match(common_s, sampleAnn$Temp_ID),]
  pathScoreMat <- pathScoreMat[match(common_s,rownames(pathScoreMat)), ,drop=FALSE]
}

# remove sample that have NAs
id_s <- rowSums(matrix(is.na(pathScoreMat), ncol = ncol(pathScoreMat), nrow = nrow(pathScoreMat))) == 0
sampleAnn <- sampleAnn[id_s,]
samplesID_new <- sampleAnn$Temp_ID
if(!all(id_s)){pathScoreMat <- pathScoreMat[id_s, ,drop=FALSE]}

print('pathScore mat loaded')

#### load phenotype annotation ####
phenoAnn <- fread(file = phenoAnn_file, header = T, stringsAsFactors = F, data.table = F)

#######################################################################################################################
# for each phenoDat/covDat file, compute association
for(n in 1:length(phenoDat_file)){
  
  phenoDat <- fread(file = phenoDat_file[n], header = T, stringsAsFactors = F, data.table = F, check.names=F)
  covDat <- fread(file = covDat_file[n], header = T, stringsAsFactors = F, data.table = F, check.names=F)
  common_s <- intersect(samplesID_new, intersect(covDat$Individual_ID,phenoDat$Individual_ID)) # keep only intersection across all datasets
  phenoDat <- phenoDat[phenoDat$Individual_ID %in% common_s, ]
  covDat <- covDat[covDat$Individual_ID %in% common_s, ]
  
  print(paste0('############ phenoFile/covFile ', names_file[n], ' ############'))
  pheno_names <- colnames(phenoDat)[!colnames(phenoDat) %in% 'Individual_ID']
  print(pheno_names)
  if(cov_corr){
    tot_var <- merge(x = covDat[, ! colnames(covDat) %in% c('genoSample_ID', 'RNASample_ID', 'Dx')], y = phenoDat, by.x = 'Individual_ID', by.y = 'Individual_ID', sort = F)
    cov_names <- colnames(covDat[, ! colnames(covDat) %in% c('Individual_ID', 'genoSample_ID', 'RNASample_ID', 'Dx')])
  }else{
    tot_var <- phenoDat
  }
  
  # filter and maintain correct order samples
  tot_var <- tot_var[match(common_s, tot_var$Individual_ID),]
  samples_tmp <- samplesID_new[match(common_s, samplesID_new)]
  print(identical(samples_tmp, tot_var$Individual_ID))
  
  id_samples <- match(common_s, rownames(pathScoreMat))
  print(identical(rownames(pathScoreMat)[id_samples], tot_var$Individual_ID))
  
  tot_var <- tot_var[, ! colnames(tot_var) %in% 'Individual_ID']
  colnames(tot_var)[colnames(tot_var) %in% pheno_names] <- paste0('p', pheno_names)
  if(cov_corr){colnames(tot_var)[colnames(tot_var) %in% cov_names] <- paste0('c', cov_names)}
  print(colnames(tot_var))  
  
  # prepare phenotype
  phenoAnn_tmp <- phenoAnn[phenoAnn$pheno_id %in% pheno_names, c('pheno_id', 'FieldID', 'Field', 'transformed_type')]
  print(str(phenoAnn_tmp))
  names_df <- t(sapply(phenoAnn_tmp$pheno_id, function(x) paste0(x, c('_beta', '_se_beta','_z_t','_pval'))))
  # names_df_qval <- sapply(phenoAnn_tmp$pheno_id, function(x) paste0(x, c('_qval')))
  
  ###############################
  #### pathScore assocaition ####
  ###############################
  pathScoreMat_association <- cbind(matrix(pathScoreMat[id_samples, ], nrow = length(id_samples), ncol = ncol(pathScoreMat)), tot_var)
  colnames(pathScoreMat_association)[1:nrow(df_info)] <- paste0('X',1:nrow(df_info))
  
  df_corr_path <- vector(mode = 'list', length = nrow(phenoAnn_tmp))
  
  for(j in 1:nrow(phenoAnn_tmp)){
    
    print(paste0('## phenotype ', phenoAnn_tmp$pheno_id[j], ' ##'))
    
    registerDoParallel(cores = min(c(ncores, nrow(df_info))))
    output <- foreach(x=1:nrow(df_info), .combine = rbind)%dopar%{
      # print(x)
      compute_reg_pheno(id_gene = x, type_pheno = phenoAnn_tmp$transformed_type[j], id_pheno = phenoAnn_tmp$pheno_id[j], total_mat = pathScoreMat_association)
    }
    
    if(nrow(df_info) == 1){output <- matrix(output, nrow = 1, ncol = length(output))}

    colnames(output) <- names_df[j,]
    rownames(output) <- NULL
    df_corr_path[[j]] <- cbind(df_info, output)
    
    
  }
  
  rm(pathScoreMat_association) # free memory
  print('pathScore completed')
  
  
  pathScore <- df_corr_path
  
  # save results
  filename <- ifelse(cov_corr, sprintf('%s%s_covCorr.RData', outFile , names_file[n]), sprintf('%s%s.RData', outFile, names_file[n]))
  save(pathScore, file = filename)
  
  print('final result saved')
  
}


