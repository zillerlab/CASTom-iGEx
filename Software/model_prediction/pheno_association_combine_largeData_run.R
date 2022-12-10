#!/usr/bin/env Rscript
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
parser$add_argument("--reactome_file", type = "character", help = "reactome pathway anntation (.gmt)")
parser$add_argument("--GOterms_file", type = "character", help = "GO pathway anntation (.RData)")
parser$add_argument("--pathScoreFold_Reactome", type = "character", help = "RData to be loaded")
parser$add_argument("--pathScoreFold_GO", type = "character", help = "RData to be loaded")
parser$add_argument("--n_split", type = "integer", default = 100, help = "total number of splitted parts")
parser$add_argument("--lambda_pi0", type = "double", default = 0.5,  help = "fixed lambda for pi1 computation")
parser$add_argument("--phenoAnn_file", type = "character", help = "file with phenotype annotation (used to determine the type of regression)")
parser$add_argument("--cov_corr", type = "logical",default = T,  help = "if T column in covDat_file are use to correct association (excluded Dx and IDs)")
parser$add_argument("--outFold", type="character", help = "Output file [basename only]")

args <- parser$parse_args()
names_file <- args$names_file
tscoreFold <- args$tscoreFold
GOterms_file <- args$GOterms_file
reactome_file <- args$reactome_file
lambda_pi0 <- args$lambda_pi0
pathScoreFold_Reactome <- args$pathScoreFold_Reactome
pathScoreFold_GO <- args$pathScoreFold_GO
n_split <- args$n_split
phenoAnn_file <- args$phenoAnn_file
cov_corr <- args$cov_corr
outFold <- args$outFold

# #########################################################################################
# names_file <- 'T1D'
# phenoAnn_file <- 'INPUT_DATA/Covariates/phenotypeDescription_T1D.txt'
# tscoreFold <- 'OUTPUT_GTEx/predict_UKBB/Pancreas/200kb/noGWAS/devgeno0.01_testdevgeno0/T1D_pheno/Association_tscore_res/'
# pathScoreFold_Reactome <- 'OUTPUT_GTEx/predict_UKBB/Pancreas/200kb/noGWAS/devgeno0.01_testdevgeno0/T1D_pheno/Association_reactome_res/'
# pathScoreFold_GO <- 'OUTPUT_GTEx/predict_UKBB/Pancreas/200kb/noGWAS/devgeno0.01_testdevgeno0/T1D_pheno/Association_GO_res/'
# n_split <- 100
# cov_corr <- T
# lambda_pi0 <- 0.5
# outFold <- 'OUTPUT_GTEx/predict_UKBB/Pancreas/200kb/noGWAS/devgeno0.01_testdevgeno0/T1D_pheno/'
# GOterms_file <- '/psycl/g/mpsziller/lucia/refData/GOterm_geneAnnotation_allOntologies.RData'
# reactome_file <- '/psycl/g/mpsziller/lucia/refData/ReactomePathways.gmt'
# #########################################################################################

# load pheno info
phenoAnn <- fread(phenoAnn_file, data.table = F, h=T)
end_name <- ifelse(cov_corr, '_covCorr.RData', '.RData')

# load pathways objects terms
go <- get(load(GOterms_file))
go_name <- sapply(go, function(x) x$GOID)

gs <- readGmt(reactome_file)
gs=lapply(gs,function(X){
  X@ids=gsub(",1.0","",X@ids)
  return(X)
})
gs_name <- sapply(gs, function(x) x@reference)

# load results
for(j in 1:length(names_file)){
  
  print(paste('Block', names_file[j]))
  
  tscore_res <- get(load(sprintf('%spval_tscore_splitGene1_%s%s', tscoreFold, names_file[j], end_name)))
  pathR_res <- get(load(sprintf('%spval_pathScore_Reactome_splitPath1_%s%s', pathScoreFold_Reactome, names_file[j], end_name)))
  pathGO_res <- get(load(sprintf('%spval_pathScore_GO_splitPath1_%s%s', pathScoreFold_GO,  names_file[j], end_name)))
  
  # find phenoid
  pheno_id <- sapply(tscore_res, function(x) strsplit(colnames(x)[5], split = '_beta')[[1]][1])
  pheno_final <- phenoAnn[phenoAnn$pheno_id %in% pheno_id, ]
  df_pi <- data.frame(pheno_id = pheno_id, tscore = 0, pathScore_reactome = 0, pathScore_GO = 0)
  
  for(i in 2:n_split){
    
    if(file.exists(sprintf('%spval_tscore_splitGene%i_%s%s', tscoreFold, i, names_file[j], end_name))){
      tmp <- get(load(sprintf('%spval_tscore_splitGene%i_%s%s', tscoreFold, i, names_file[j], end_name)))
      tscore_res <- mapply(function(x,y) rbind(x,y) , x = tscore_res, y = tmp, SIMPLIFY = F)  
    }else{
      print(sprintf('tscore file %i not existing', i))
    }
    
    if(file.exists(sprintf('%spval_pathScore_Reactome_splitPath%i_%s%s', pathScoreFold_Reactome, i, names_file[j], end_name))){
      tmp <- get(load(sprintf('%spval_pathScore_Reactome_splitPath%i_%s%s', pathScoreFold_Reactome, i, names_file[j], end_name)))
      pathR_res <- mapply(function(x,y) rbind(x,y) , x = pathR_res, y = tmp, SIMPLIFY = F)
    }else{
      print(sprintf('path Reactome file %i not existing', i))
    }
    
    if(file.exists(sprintf('%spval_pathScore_GO_splitPath%i_%s%s', pathScoreFold_GO, i, names_file[j], end_name))){
      tmp <- get(load(sprintf('%spval_pathScore_GO_splitPath%i_%s%s', pathScoreFold_GO, i, names_file[j], end_name)))
      pathGO_res <- mapply(function(x,y) rbind(x,y) , x = pathGO_res, y = tmp, SIMPLIFY = F)
    }else{
      print(sprintf('path GO file %i not existing', i))
    }
  }
  
  # correct pvalues
  qval_t <- lapply(1:length(tscore_res), function(x) qvalue(tscore_res[[x]][, paste0(pheno_id[x], '_pval')])$qvalues)
  df_pi$tscore <- sapply(1:length(tscore_res), function(x) 1 - pi0est(tscore_res[[x]][, paste0(pheno_id[x], '_pval')], lambda = 0.5)$pi0)
  BH_t <- lapply(1:length(tscore_res), function(x) p.adjust(tscore_res[[x]][, paste0(pheno_id[x], '_pval')], method = 'BH'))
  
  qval_pR <- lapply(1:length(pathR_res), function(x) qvalue(pathR_res[[x]][, paste0(pheno_id[x], '_pval')])$qvalues)
  df_pi$pathScore_reactome <- sapply(1:length(pathR_res), function(x) 1 - pi0est(pathR_res[[x]][, paste0(pheno_id[x], '_pval')], lambda = 0.5)$pi0)
  BH_pR <- lapply(1:length(pathR_res), function(x) p.adjust(pathR_res[[x]][, paste0(pheno_id[x], '_pval')], method = 'BH'))
  
  qval_pG <- lapply(1:length(pathGO_res), function(x) qvalue(pathGO_res[[x]][, paste0(pheno_id[x], '_pval')])$qvalues)
  df_pi$pathScore_GO <- sapply(1:length(pathGO_res), function(x) 1 - pi0est(pathGO_res[[x]][, paste0(pheno_id[x], '_pval')], lambda = 0.5)$pi0)
  BH_pG <- lapply(1:length(pathGO_res), function(x) p.adjust(pathGO_res[[x]][, paste0(pheno_id[x], '_pval')], method = 'BH'))
  
  info_pathR <- list()
  info_pathGO <- list()
  
  for(k in 1:length(pheno_id)){
    
    print(k)
    
    tscore_res[[k]] <- cbind(tscore_res[[k]], qval_t[[k]], BH_t[[k]])
    colnames(tscore_res[[k]])[(ncol(tscore_res[[k]])-1):ncol(tscore_res[[k]])] <- paste0(pheno_id[k], c('_qval', '_BHcorr'))
    
    pathR_res[[k]] <- cbind(pathR_res[[k]], qval_pR[[k]], BH_pR[[k]])
    colnames(pathR_res[[k]])[(ncol(pathR_res[[k]])-1):ncol(pathR_res[[k]])] <- paste0(pheno_id[k], c('_qval', '_BHcorr'))
    
    pathGO_res[[k]] <- cbind(pathGO_res[[k]], qval_pG[[k]], BH_pG[[k]])
    colnames(pathGO_res[[k]])[(ncol(pathGO_res[[k]])-1):ncol(pathGO_res[[k]])] <- paste0(pheno_id[k], c('_qval', '_BHcorr'))
    
    # create list object for each pathway with tscore association results
    tmp <- lapply(pathR_res[[k]]$path, function(x) gs[which(gs_name == x)][[1]])
    info_pathR[[k]] <- lapply(1:nrow(pathR_res[[k]]), function(x) 
      list(path = pathR_res[[k]][x,], genes_path = tmp[[x]]@ids, tscore = tscore_res[[k]][tscore_res[[k]]$external_gene_name %in% tmp[[x]]@ids,]))
    
    tmp <- lapply(pathGO_res[[k]]$path_id, function(x) go[which(go_name == x)][[1]])
    info_pathGO[[k]] <- lapply(1:nrow(pathGO_res[[k]]), function(x) 
      list(path = pathGO_res[[k]][x,], genes_path = tmp[[x]]$geneIds, tscore = tscore_res[[k]][tscore_res[[k]]$external_gene_name %in% tmp[[x]]$geneIds,]))
    
  }
  
  # save final results
  final <- list(pheno = pheno_final, tscore = tscore_res, pathScore_reactome = pathR_res, pathScore_GO = pathGO_res, pi1_lambdafixed = df_pi, info_pathScore_reactome = info_pathR, info_pathScore_GO = info_pathGO)
  save(final, file = sprintf('%spval_%s_pheno%s', outFold, names_file[j], end_name))
  
}

