#### code written by Lucia Trastulla, e-mail: lucia_trastulla@psych.mpg.de ####
# save also beta, se, zcores
# consider results from big matrix, possibility to load more than 1 cov and pheno dat, add pheno annotation info
# specific for tscores

options(stringsAsFactors=F)
options(max.print=1000)
suppressPackageStartupMessages(library(argparse))
suppressPackageStartupMessages(library(bigmemory))
suppressPackageStartupMessages(library(nnet))
suppressPackageStartupMessages(library(lattice))
suppressPackageStartupMessages(library(data.table))
suppressPackageStartupMessages(library(MASS))
suppressPackageStartupMessages(library(lmtest))
suppressPackageStartupMessages(library(doParallel))
suppressPackageStartupMessages(library(PGSEA))

parser <- ArgumentParser(description="T-scores association analysis")
parser$add_argument("--inputFile", type = "character", help = "RData to be loaded")
parser$add_argument("--inputInfoFile", type = "character",  help = "RData to be loaded, contains info on genes or path")
parser$add_argument("--sampleAnn_file", type = "character", help = "file with sample info for new genotype data, must contain Dx column (0 controls 1 cases)")
parser$add_argument("--covDat_file", type = "character", nargs = '*', help = "file with sample info for new genotype data, must contain Dx column (0 controls 1 cases), associated to phenoDat_file")
parser$add_argument("--phenoDat_file", type = "character", nargs = '*', help = "file(s) with individual_ID to match and phenotype to test association, associated to covDat_file")
parser$add_argument("--names_file", type = "character", nargs = '*', help = "for each couple of covDat/phenoDat file, associated name")
parser$add_argument("--phenoAnn_file", type = "character", help = "file with phenotype annotation (used to determine the type of regression)")
parser$add_argument("--cov_corr", type = "logical",default = F,  help = "if T column in covDat_file are use to correct association (excluded Dx and IDs)")
parser$add_argument("--functR", type = "character", help = "Rscript with functions to be used")
parser$add_argument("--ncores", type = "integer", help = "cores parallelization")
parser$add_argument("--split_gene_id", type="integer", help = "split partiton for genes")
parser$add_argument("--split_tot", type="integer", default = 100, help = "total number of partitions")
parser$add_argument("--outFile", type="character", help = "Output file [basename only]")

args <- parser$parse_args()
inputFile <- args$inputFile
inputInfoFile <- args$inputInfoFile
split_tot <- args$split_tot
split_gene_id <- args$split_gene_id
covDat_file <- args$covDat_file
sampleAnn_file <- args$sampleAnn_file
phenoDat_file <- args$phenoDat_file
phenoAnn_file <- args$phenoAnn_file
cov_corr <- args$cov_corr
names_file <- args$names_file
functR <- args$functR
ncores <- args$ncores
outFile <- args$outFile

# #########################################################################################
# inputFile <- paste0('OUTPUT_CMC/predict_UKBB/200kb/devgeno0.01_testdevgeno0/predictedTscores_splitGenes', 1, '.RData')
# inputInfoFile <- 'OUTPUT_CMC/predict_UKBB/200kb/devgeno0.01_testdevgeno0/tscore_info.RData'
# covDat_file <- read.table('INPUT_DATA/Covariates/match_cov_pheno_SchunkertApp.txt', h=F, stringsAsFactors = F)$V3
# phenoDat_file <- read.table('INPUT_DATA/Covariates/match_cov_pheno_SchunkertApp.txt', h=F, stringsAsFactors = F)$V2
# names_file <- read.table('INPUT_DATA/Covariates/match_cov_pheno.txt', h=F, stringsAsFactors = F)$V1
# sampleAnn_file <- 'INPUT_DATA/Covariates/covariatesMatrix.txt'
# cov_corr <- T
# phenoAnn_file <-'/psycl/g/mpsziller/lucia/CAD_UKBB/eQTL_PROJECT/INPUT_DATA_GTEx/CAD/Covariates/UKBB/phenotypeDescription_PHESANTproc_CADrelatedpheno.txt'
# outFold <- 'OUTPUT_GTEx/predict_UKBB/Brain_Cerebellum/200kb/noGWAS/devgeno0.01_testdevgeno0/Association_tscore_res/'
# functR <- 'RSCRIPTS/SCRIPTS_v2/pheno_association_functions.R'
# split_tot <- 100
# split_gene_id <- 1
# #########################################################################################

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

all_genesID <- df_info$external_gene_name
step_genes <- round(length(all_genesID)/split_tot)
split_id = as.vector(sapply(1:split_tot, function(x) rep(x, step_genes)))
if(length(split_id)<length(all_genesID)){
  split_id <- c(split_id, rep(split_tot, length(all_genesID) - length(split_id)))  
}else{
  if(length(split_id)>length(all_genesID)){
    split_id <- split_id[-((length(all_genesID)+1):length(split_id))]
  }
}
split_df <- data.frame(id = 1:length(all_genesID), split = split_id)
df_info <- df_info[split_df$id[split_df$split == split_gene_id], ]

print(paste('same gene annotation:', identical(df_info$external_gene_name, score_mat$geneId)))

# # filter geneAnn
# if(!identical(score_mat$geneId, df_info$external_gene_name)){
#    print('adjust genes annotation')
#    id <- sapply(score_mat$geneId, function(x) which(x == df_info$external_gene_name))
#    df_info <- df_info[id,]
#    print(identical(score_mat$geneId, df_info$external_gene_name))
# }

## filter info based on score mat
#df_info <- df_info[df_info$external_gene_name %in% score_mat$geneId, ]

#if(!identical(df_info$external_gene_name, score_mat$geneId)){
#  rep_gene <- names(which(table(df_info$external_gene_name)>1))
#  for(i in 1:length(rep_gene)){
#    tmp <- df_info[df_info$external_gene_name %in% rep_gene[i], ]
#    new_gene = data.frame(ensembl_gene_id = paste0(tmp$ensembl_gene_id, collapse = '_'), external_gene_name = rep_gene[i], dev_geno = mean(tmp$dev_geno), test_dev_geno = mean(tmp$test_dev_geno))
#    df_info <- df_info[! df_info$external_gene_name %in%  rep_gene[i], ]
#    df_info <- rbind(df_info, new_gene)
#  }
#  df_info <- df_info[match(score_mat$geneId, df_info$external_gene_name),]
#}

tscoreMat <- matrix(t(score_mat[, -1]), ncol = nrow(score_mat), nrow = ncol(score_mat)-1)
rownames(tscoreMat) <- colnames(score_mat)[-1]
tscoreMat <- tscoreMat[rownames(tscoreMat) %in% sampleAnn$Temp_ID,]

# match
if(!identical(rownames(tscoreMat), sampleAnn$Temp_ID)){
  sampleAnn <- sampleAnn[match(rownames(tscoreMat), sampleAnn$Temp_ID),]
}

# remove sample that have NAs
id_s <- rowSums(is.na(tscoreMat)) == 0
sampleAnn <- sampleAnn[id_s,]
samplesID_new <- sampleAnn$Temp_ID
if(!all(id_s)){
  tmp <- rownames(tscoreMat)
  tscoreMat <- matrix(tscoreMat[id_s, ], ncol = ncol(tscoreMat), nrow = length(which(id_s)))
  rownames(tscoreMat) <- tmp[id_s]
}

print(str(tscoreMat))
#print(dim(df_info))

print('Tscore mat loaded')

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
  
  id_samples <- match(common_s, rownames(tscoreMat))
  print(identical(rownames(tscoreMat)[id_samples], tot_var$Individual_ID))

  tot_var <- tot_var[, ! colnames(tot_var) %in% 'Individual_ID']
  colnames(tot_var)[colnames(tot_var) %in% pheno_names] <- paste0('p', pheno_names)
  if(cov_corr){colnames(tot_var)[colnames(tot_var) %in% cov_names] <- paste0('c', cov_names)}
  print(colnames(tot_var))  
  
  # prepare phenotype
  phenoAnn_tmp <- phenoAnn[phenoAnn$pheno_id %in% pheno_names, c('pheno_id', 'FieldID', 'Field', 'transformed_type')]
  print(str(phenoAnn_tmp))
  names_df <- t(sapply(phenoAnn_tmp$pheno_id, function(x) paste0(x, c('_beta', '_se_beta','_z_t','_pval'))))

  ############################
  #### tscore assocaition ####
  ############################
  tscoreMat_association <- cbind(matrix(tscoreMat[id_samples, ], nrow = length(id_samples), ncol = ncol(tscoreMat)), tot_var)
  colnames(tscoreMat_association)[1:nrow(df_info)] <- paste0('X',1:nrow(df_info))
  print(str(tscoreMat_association))
  
  df_corr_tscore <- vector(mode = 'list', length = nrow(phenoAnn_tmp))

  for(j in 1:nrow(phenoAnn_tmp)){
    
    print(paste0('## phenotype ', phenoAnn_tmp$pheno_id[j], ' ##'))
  
    registerDoParallel(cores = min(c(ncores, nrow(df_info))))
    output <- foreach(x=1:nrow(df_info), .combine = rbind)%dopar%{
      # print(x)
       compute_reg_pheno(id_gene = x, type_pheno = phenoAnn_tmp$transformed_type[j], id_pheno = phenoAnn_tmp$pheno_id[j], total_mat = tscoreMat_association)
    }
    
    if(nrow(df_info) == 1){output <- matrix(output, nrow = 1, ncol = length(output))}
    # if(any(sapply(output, function(x) x == 'wrong pheno type annotation'))){ stop("wrong pheno type annotation") }
    
    colnames(output) <- names_df[j,]
    rownames(output) <- NULL
    df_corr_tscore[[j]] <- cbind(df_info, output)
    
    
  }
  
  rm(tscoreMat_association) # free memory
  print('tscore completed')
  
 
  tscore <- df_corr_tscore
  
  # save results
  filename <- ifelse(cov_corr, sprintf('%s%s_covCorr.RData', outFile , names_file[n]), sprintf('%s%s.RData', outFile, names_file[n]))
  save(tscore, file = filename)
  
  print('final result saved')
  
}

