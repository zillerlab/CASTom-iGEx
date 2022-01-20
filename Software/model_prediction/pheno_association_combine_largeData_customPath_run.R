# combine splitted results association

options(stringsAsFactors=F)
options(max.print=1000)
suppressPackageStartupMessages(library(argparse))
suppressPackageStartupMessages(library(data.table))
suppressPackageStartupMessages(library(doParallel))
suppressPackageStartupMessages(library(PGSEA))
suppressPackageStartupMessages(library(qvalue))

parser <- ArgumentParser(description="combine results association analysis")
parser$add_argument("--names_file", type = "character", nargs = '*', help = "for each couple of covDat/phenoDat file, associated name")
parser$add_argument("--tscoreFold", type = "character", help = "RData to be loaded")
parser$add_argument("--pathwayCustom_name", type = "character", help = "name pathway custom")
parser$add_argument("--pathwayCustom_file", type = "character", default = NULL, help = "custom pathway structure (.RData)")
parser$add_argument("--pathScoreFold", type = "character", help = "RData to be loaded")
parser$add_argument("--n_split", type = "integer", default = 100, help = "total number of splitted parts")
parser$add_argument("--lambda_pi0", type = "double", default = 0.5,  help = "fixed lambda for pi1 computation")
parser$add_argument("--phenoAnn_file", type = "character", help = "file with phenotype annotation (used to determine the type of regression)")
parser$add_argument("--cov_corr", type = "logical",default = F,  help = "if T column in covDat_file are use to correct association (excluded Dx and IDs)")
parser$add_argument("--outFold", type="character", help = "Output file [basename only]")

args <- parser$parse_args()
names_file <- args$names_file
tscoreFold <- args$tscoreFold
pathwayCustom_name  <- args$pathwayCustom_name 
pathwayCustom_file <- args$pathwayCustom_file
lambda_pi0 <- args$lambda_pi0
pathScoreFold <- args$pathScoreFold
n_split <- args$n_split
phenoAnn_file <- args$phenoAnn_file
cov_corr <- args$cov_corr
outFold <- args$outFold

#########################################################################################
# names_file <- 'CAD'
# phenoAnn_file <- 'INPUT_DATA_GTEx/CAD/Covariates/UKBB/phenotypeDescription_CAD.txt'
# tscoreFold <- 'OUTPUT_GTEx/predict_CAD/Whole_Blood/200kb/CAD_GWAS_bin5e-2/UKBB/devgeno0.01_testdevgeno0/Association_tscore_res/'
# pathScoreFold <- 'OUTPUT_GTEx/predict_CAD/Whole_Blood/200kb/CAD_GWAS_bin5e-2/UKBB/devgeno0.01_testdevgeno0/Association_geneSets_sameLocus_res/'
# n_split <- 100
# cov_corr <- T
# lambda_pi0 <- 0.5
# outFold <- 'OUTPUT_GTEx/predict_CAD/Whole_Blood/200kb/CAD_GWAS_bin5e-2/UKBB/devgeno0.01_testdevgeno0/'
# pathwayCustom_file <- 'OUTPUT_GTEx/predict_CAD/Whole_Blood/200kb/CAD_GWAS_bin5e-2/UKBB/devgeno0.01_testdevgeno0//geneSets_sameLocus_sameSign.RData'
# pathwayCustom_name <- 'geneSets_sameLocus'
#########################################################################################

# load pheno info
phenoAnn <- fread(phenoAnn_file, data.table = F, h=T)
end_name <- ifelse(cov_corr, '_covCorr.RData', '.RData')
end_name_save <- ifelse(cov_corr, sprintf('_covCorr_customPath_%s.RData', pathwayCustom_name), 
                        sprintf('_customPath_%s.RData', pathwayCustom_name))

path_custom <- get(load(pathwayCustom_file))
path_name <- sapply(path_custom, function(x) x$name)


# load results
for(j in 1:length(names_file)){
  
  print(paste('Block', names_file[j]))
  
  tscore_res <- get(load(sprintf('%spval_tscore_splitGene1_%s%s', tscoreFold, names_file[j], end_name)))
  path_res <- get(load(sprintf('%spval_pathScore_%s_splitPath1_%s%s',
                               pathScoreFold, pathwayCustom_name, names_file[j], end_name)))
  
  # find phenoid
  pheno_id <- sapply(tscore_res, function(x) strsplit(colnames(x)[5], split = '_beta')[[1]][1])
  pheno_final <- phenoAnn[phenoAnn$pheno_id %in% pheno_id, ]
  df_pi <- data.frame(pheno_id = pheno_id, tscore = 0, pathScore = 0)
  
  for(i in 2:n_split){
    
    if(file.exists(sprintf('%spval_tscore_splitGene%i_%s%s', tscoreFold, i, names_file[j], end_name))){
      tmp <- get(load(sprintf('%spval_tscore_splitGene%i_%s%s', tscoreFold, i, names_file[j], end_name)))
      tscore_res <- mapply(function(x,y) rbind(x,y) , x = tscore_res, y = tmp, SIMPLIFY = F)  
    }else{
      print(sprintf('tscore file %i not existing', i))
    }
    
    if(file.exists(sprintf('%spval_pathScore_%s_splitPath%i_%s%s', 
                           pathScoreFold, pathwayCustom_name, i, names_file[j], end_name))){
      tmp <- get(load(sprintf('%spval_pathScore_%s_splitPath%i_%s%s', 
                              pathScoreFold, pathwayCustom_name, i, names_file[j], end_name)))
      path_res <- mapply(function(x,y) rbind(x,y) , x = path_res, y = tmp, SIMPLIFY = F)
    }else{
      print(sprintf('path %s file %i not existing', pathwayCustom_name, i))
    }
    
  }
  
  pheno_id_path <- sapply(path_res, function(x) strsplit(colnames(x)[10], split = '_beta')[[1]][1])
  path_res <- path_res[match(pheno_id, pheno_id_path)]
  
  # correct pvalues
  qval_t <- try(lapply(1:length(tscore_res), function(x) qvalue(tscore_res[[x]][, paste0(pheno_id[x], '_pval')])$qvalues))
  df_pi$tscore <- try(sapply(1:length(tscore_res), function(x) 1 - pi0est(tscore_res[[x]][, paste0(pheno_id[x], '_pval')], lambda = 0.5)$pi0))
  BH_t <- lapply(1:length(tscore_res), function(x) p.adjust(tscore_res[[x]][, paste0(pheno_id[x], '_pval')], method = 'BH'))
  
  qval_p <- try(lapply(1:length(path_res), function(x) qvalue(path_res[[x]][, paste0(pheno_id[x], '_pval')])$qvalues))
  df_pi$pathScore <- try(sapply(1:length(path_res), function(x) 1 - pi0est(path_res[[x]][, paste0(pheno_id[x], '_pval')], lambda = 0.5)$pi0))
  BH_p <- lapply(1:length(path_res), function(x) p.adjust(path_res[[x]][, paste0(pheno_id[x], '_pval')], method = 'BH'))
  
  if(class(qval_p) == 'try-error'){
    qval_p <- lapply(1:length(path_res), function(x) rep(NA, nrow(path_res[[x]])))
  }
  
  info_path <- list()
  
  for(k in 1:length(pheno_id)){
    
    print(k)
    
    tscore_res[[k]] <- cbind(tscore_res[[k]], qval_t[[k]], BH_t[[k]])
    colnames(tscore_res[[k]])[(ncol(tscore_res[[k]])-1):ncol(tscore_res[[k]])] <- paste0(pheno_id[k], c('_qval', '_BHcorr'))
    
    path_res[[k]] <- cbind(path_res[[k]], qval_p[[k]], BH_p[[k]])
    colnames(path_res[[k]])[(ncol(path_res[[k]])-1):ncol(path_res[[k]])] <- paste0(pheno_id[k], c('_qval', '_BHcorr'))
    
    # create list object for each pathway with tscore association results
    tmp <- lapply(path_res[[k]]$path, function(x) path_custom[which(path_name == x)][[1]])
    info_path[[k]] <- lapply(1:nrow(path_res[[k]]), function(x) 
      list(path = path_res[[k]][x,], genes_path = tmp[[x]]$geneIds, tscore = tscore_res[[k]][tscore_res[[k]]$external_gene_name %in% tmp[[x]]$geneIds,]))
  }
  
  # save final results
  final <- list(pheno = pheno_final, tscore = tscore_res, pathScore = path_res, pi1_lambdafixed = df_pi, 
                info_pathScore = info_path)
  save(final, file = sprintf('%spval_%s_pheno%s', outFold, names_file[j], end_name_save))
  
}

