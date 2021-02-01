# for each phenotype associate R2
# if continuous, common R2 otherwise Lee libility scale R2?

options(stringsAsFactors=F)
options(max.print=1000)
suppressPackageStartupMessages(library(argparse))
suppressPackageStartupMessages(library(e1071))
suppressPackageStartupMessages(library(data.table))
suppressPackageStartupMessages(library(doParallel))
suppressPackageStartupMessages(library(qvalue))
suppressPackageStartupMessages(library(pROC))
suppressPackageStartupMessages(library(pryr))
suppressPackageStartupMessages(library(umap))
suppressPackageStartupMessages(library(igraph))
suppressPackageStartupMessages(library(Matrix))
suppressPackageStartupMessages(library(SparseM))
suppressPackageStartupMessages(library(ggplot2))
suppressPackageStartupMessages(library(RColorBrewer))
suppressPackageStartupMessages(library(pheatmap))
suppressPackageStartupMessages(library(ggplot2))
suppressPackageStartupMessages(library(data.table))
options(bitmapType = 'cairo', device = 'png')

parser <- ArgumentParser(description="R2 from phenotype - risk score")
parser$add_argument("--riskScore_file", type = "character", help = "")
parser$add_argument("--sampleAnn_file", type = "character", help = "")
parser$add_argument("--pheno_file", nargs = '*',type = "character", help = "")
parser$add_argument("--phenoAnn_file", type = "character", help = "")
parser$add_argument("--names_pheno", nargs = '*', type = "character", help = "")
parser$add_argument("--outFold", type="character", help = "Output file [basename only]")

args <- parser$parse_args()
riskScore_file <- args$riskScore_file
sampleAnn_file <- args$sampleAnn_file
pheno_file <- args$pheno_file
phenoAnn_file <- args$phenoAnn_file
names_pheno <- args$names_pheno
outFold <- args$outFold

###################################################################################################################
# sampleAnn_file <- '/psycl/g/mpsziller/lucia/CAD_UKBB/eQTL_PROJECT/INPUT_DATA_GTEx/CAD/Covariates/UKBB//covariateMatrix_latestW.txt'
# outFold <- '/psycl/g/mpsziller/lucia/CAD_UKBB/eQTL_PROJECT/OUTPUT_GTEx/predict_CAD/Liver/200kb/CAD_GWAS_bin5e-2/UKBB/devgeno0.01_testdevgeno0/tscore_corrThr0.5_'
# riskScore_file <- '/psycl/g/mpsziller/lucia/CAD_UKBB/eQTL_PROJECT/OUTPUT_GTEx/predict_CAD/Liver/200kb/CAD_GWAS_bin5e-2/UKBB/devgeno0.01_testdevgeno0/CAD_HARD_clustering/tscore_corrThr0.5_risk_score_relatedPhenotypes.txt'
# names_pheno <- c('Blood_biochemistry', 'Blood_count', 'ICD10_Endocrine', 'ICD10_Circulatory_system', 'Blood_count_ratio')
# pheno_file <- paste0('/psycl/g/mpsziller/lucia/CAD_UKBB/eQTL_PROJECT/INPUT_DATA_GTEx/CAD/Covariates/UKBB/phenotypeMatrix_',names_pheno ,'.txt')
# phenoAnn_file <- '/psycl/g/mpsziller/lucia/CAD_UKBB/eQTL_PROJECT/INPUT_DATA_GTEx/CAD/Covariates/UKBB/CAD_HARD_clustering/phenotypeDescription.txt'
##################################################################################################################

##############################################
# function to evaluate rs-phenotype
R2_risk_score <- function(fmla, fmla_cov, type_pheno, mat){
  
  R2_tot <- NA
  R2_cov <- NA
  R2_diff <- NA
  Fstat_diff <- NA
  N <- sum(!is.na(mat$pheno))
  
  res <- tryCatch(lm(fmla, data = mat),warning=function(...) NA, error=function(...) NA)
  if(is.list(res)){
    s_res <- summary(res)
    # R2_tot <- 1-s_res$deviance/s_res$null.deviance
    R2_tot <- s_res$r.squared  
  }
  
  res_cov <- tryCatch(lm(fmla_cov, data = mat),warning=function(...) NA, error=function(...) NA)
  if(is.list(res_cov)){
    s_res_cov <- summary(res_cov)
    # R2_cov <- 1-s_res$deviance/s_res$null.deviance
    R2_cov <- s_res_cov$r.squared  
  }
  
  if(is.list(res_cov) & is.list(res)){
    R2_diff <- R2_tot - R2_cov
    # Fstat_diff <- (N-nrow(s_res$coefficients))*(sum(s_res_cov$residuals^2) - sum(s_res$residuals^2))/sum(s_res$residuals^2)
    comp_mod <- anova(res_cov, res)
    Fstat_diff <-  comp_mod$F[2]

  }
  
  return(c(R2_diff, Fstat_diff, R2_tot, R2_cov))
  
}
##############################################

sampleAnn <- read.table(sampleAnn_file, h=T, stringsAsFactors = F)

phenoAnn <- fread(phenoAnn_file, h=T, stringsAsFactors = F, data.table = F)

# load risk_score_file
riskScore <- fread(riskScore_file, h=T, stringsAsFactors = F, data.table = F)

common_s <- intersect(riskScore$Individual_ID, sampleAnn$Individual_ID)
riskScore <- riskScore[match(common_s, riskScore$Individual_ID), ]
sampleAnn <- sampleAnn[match(common_s, sampleAnn$Individual_ID), ]
covDat <- sampleAnn[, !colnames(sampleAnn) %in% c('Individual_ID', 'Dx', 'genoSample_ID')]

id_pheno <- intersect(colnames(riskScore), phenoAnn$pheno_id)
phenoAnn <- phenoAnn[match(id_pheno, phenoAnn$pheno_id), ]

summary_res <- list()
fmla <- as.formula(paste('pheno~risk_score+', paste0(colnames(covDat), collapse = '+')))
fmla_cov <- as.formula(paste('pheno~', paste0(colnames(covDat), collapse = '+')))

for(i in 1:length(names_pheno)){
  
  print(names_pheno[i])
  
  pheno_tmp <-  fread(pheno_file[i], h=T, stringsAsFactors = F, data.table = F)
  pheno_tmp <- pheno_tmp[match(common_s, pheno_tmp$Individual_ID), ]
  pheno_tmp <- pheno_tmp[, colnames(pheno_tmp) %in% phenoAnn$pheno_id, drop = F]
  print(dim(pheno_tmp))
  summary_res[[i]] <- data.frame(pheno_id = colnames(pheno_tmp), R2_risk = NA, Fstat_risk = NA, 
                                 R2_tot = NA, R2_cov = NA, nsamples = NA, nsamples_T = NA)
  
  for(l in 1:ncol(pheno_tmp)){
    
    print(colnames(pheno_tmp)[l])
    
    mat <- cbind(data.frame(pheno = pheno_tmp[,l], risk_score = riskScore[, colnames(pheno_tmp)[l]]), covDat)
    # mat$pheno <- scale(mat$pheno)[,1]
    # attr(mat$pheno, 'scaled:center') <- NULL
    # attr(mat$pheno, 'scaled:scale') <- NULL
    type_pheno <- phenoAnn$transformed_type[phenoAnn$pheno_id== colnames(pheno_tmp)[l]]
    tmp <- R2_risk_score(fmla = fmla, fmla_cov = fmla_cov, type_pheno = type_pheno, mat = mat)
    summary_res[[i]][l, 2:5] <- tmp
    summary_res[[i]]$nsamples[l] <- sum(!is.na(mat$pheno))
    if(!type_pheno %in% c('CONTINUOUS', 'CAT_ORD')){
      summary_res[[i]]$nsamples_T[l] <- sum(pheno_tmp[,l] == 1 & !is.na(mat$pheno))
    }
  }
}

# save results
summary_res <- do.call(rbind, summary_res)
summary_res <- cbind(summary_res, phenoAnn[match(summary_res$pheno_id, phenoAnn$pheno_id), 
                                           c('Field', 'Coding_meaning', 'transformed_type')])

write.table(file = sprintf('%sR2_risk_score_phenotype.txt', outFold), 
            x = summary_res, quote = F, sep = '\t', col.names = T, row.names = F)


