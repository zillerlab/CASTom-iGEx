#!/usr/bin/env Rscript
# meta - analysis across different cohorts

options(stringsAsFactors=F)
options(max.print=1000)
suppressPackageStartupMessages(library(argparse))
suppressPackageStartupMessages(library(qvalue))
suppressPackageStartupMessages(library(biomaRt))


parser <- ArgumentParser(description="Meta analysis tscore and pathway scores")
parser$add_argument("--res_cohorts", type="character", nargs = '*', help = "RData file with the phenotype association results")
parser$add_argument("--phenoDatFile_cohorts", type="character", nargs = '*', help = "file with the phenotype info for each cohort")
parser$add_argument("--name_cohort", type="character", nargs = '*', help = "names of the cohorts")
parser$add_argument("--cov_corr", type = "logical", default = T,  help = "if T column in covDat_file are use to correct association (excluded Dx and IDs) [default %(default)s]")
parser$add_argument("--phenoName", type="character",  help = "phenotype group name considered")
parser$add_argument("--thr_het", type = "double", default = 0.001, help = "threshold for random effect model [default %(default)s]")
parser$add_argument("--outFold", type="character", help = "Output file [basename only]")

args <- parser$parse_args()
res_cohorts <- args$res_cohorts
name_cohort <- args$name_cohort
phenoDatFile_cohorts <- args$phenoDatFile_cohorts
cov_corr <- args$cov_corr
phenoName <- args$phenoName
thr_het <- args$thr_het
outFold <- args$outFold

###########################################################################################
# res_cohorts <- sapply(c('CG','LURIC', 'WTCCC', 'MG', paste0('German', 1:5)) , function(x) sprintf('OUTPUT_GTEx/predict_CAD/Adipose_Subcutaneous/200kb/CAD_GWAS_bin5e-2/%s/devgeno0.01_testdevgeno0/pval_CAD_pheno_covCorr.RData', x))
# name_cohort <- c('CG','LURIC', 'WTCCC', 'MG', paste0('German', 1:5))
# phenoDatFile_cohorts <- sapply(c('CG','LURIC', 'WTCCC', 'MG',paste0('German', 1:5)) , function(x) sprintf('INPUT_DATA_GTEx/CAD/Covariates/%s/phenoMatrix.txt', x))
# cov_corr = T
# thr_het <- 10^-3
# phenoName <- 'Dx'
# lambda_pi1 <- 0.5
# GOterms_file <- '/psycl/g/mpsziller/lucia/refData/GOterm_geneAnnotation_allOntologies.RData'
# reactome_file <- '/psycl/g/mpsziller/lucia/refData/ReactomePathways.gmt'
# outFold <- 'OUTPUT_GTEx/predict_CAD/Adipose_Subcutaneous/200kb/CAD_GWAS_bin5e-2/Meta_Analysis_CAD/'
################################################################################################################################


# make a list of all the possible phenotype to consider + summary samples
pheno_tab <- list()
df_len_pheno <- list()
for(i in 1:length(phenoDatFile_cohorts)){
  
  # note: some phenotype were previously removed for missingness
  pheno_tab[[i]] <- read.table(phenoDatFile_cohorts[i], h=T, stringsAsFactors = F, sep = '\t')
 
  len_pheno_single <- rep(1, ncol(pheno_tab[[i]])-1)
  
  type_val <- apply(as.matrix(pheno_tab[[i]][, -1], ncol = len_pheno_single), 2, function(x) length(unique(x[!is.na(x)])))
 
  type <- rep('binomial', length(type_val))
  type[type_val>10] <- 'gaussian'
  type[type_val<=10 & type_val>2] <- 'gaussian'
 
  
  df_len_pheno[[i]] <- data.frame(pheno = colnames(pheno_tab[[i]])[-1], len_pheno = len_pheno_single, type = type, cohort = name_cohort[i])
  pheno_tab[[i]]$cohort <- name_cohort[i]
  
}

# find phenotype that are in common and have the same amount of info
df_len_pheno <- do.call(rbind, df_len_pheno)

# summary pheno with number of cohorts
df_pheno_tot <- data.frame(pheno = unique(df_len_pheno$pheno), n_studies = NA, cohort = NA, n_samples = NA, n_cases = NA, n_controls = NA, type = NA, binom_0 = NA, binom_1 = NA)
for(i in 1:nrow(df_pheno_tot)){
  
  tmp <- df_len_pheno[df_len_pheno$pheno ==  df_pheno_tot$pheno[i],] 
  df_pheno_tot$n_studies[i] <- nrow(tmp)
  
  if(length(unique(tmp$len_pheno))>1 | length(unique(tmp$type))>1){
    df_pheno_tot$cohort[i] <- 'remove_inconsistent'
  }else{
    df_pheno_tot$type[i] <- unique(tmp$type)
    df_pheno_tot$cohort[i] <- paste(tmp$cohort, collapse = '_')
    id <- lapply(pheno_tab, function(x) !is.na(x[, colnames(x) %in% df_pheno_tot$pheno[i]]))
    df_pheno_tot$n_samples[i] <- sum(sapply(pheno_tab, function(x) length(na.omit(x[, colnames(x) %in% df_pheno_tot$pheno[i]]))))
    df_pheno_tot$n_cases[i] <- sum(mapply(function(x,y) length(which(x$Dx[y] == 1)), x = pheno_tab, y = id))
    df_pheno_tot$n_controls[i] <- sum(mapply(function(x,y) length(which(x$Dx[y] == 0)), x = pheno_tab, y = id))
    
    df_pheno_tot$binom_1[i] <- sum(mapply(function(x,y) length(which(x[y, df_pheno_tot$pheno[i]] == 1)), x = pheno_tab, y = id))
    df_pheno_tot$binom_0[i] <- sum(mapply(function(x,y) length(which(x[y, df_pheno_tot$pheno[i]] == 0)), x = pheno_tab, y = id))
    
  }
}

# save table
file_name <- sprintf('%sphenoInfo_%s_cohorts.txt', outFold, phenoName)
write.table(df_pheno_tot, file = file_name, col.names = T, row.names = F, sep = '\t', quote = F)

#####################################################################################################
## extraxt pval for each phenotype that have more tha 1 cohort or is not insonsistent
id_pheno <- which(df_pheno_tot$n_studies > 1 & df_pheno_tot$cohort != 'remove_inconsistent')
pheno_var <- lapply(id_pheno, function(x) list(pheno = 0, cohorts = 0))
for(i in 1:length(id_pheno)){
  
  pheno_var[[i]]$pheno <- df_pheno_tot$pheno[id_pheno[i]]
  pheno_var[[i]]$cohorts <- name_cohort[sapply(name_cohort, function(x) grepl(x,df_pheno_tot$cohort[id_pheno[i]]))]
  
}

tmp_tscore <- list()

for(i in 1:length(res_cohorts)){
  res <- get(load(res_cohorts[[i]]))
  tmp_tscore[[i]] <- list()
  
  for(j in 1:length(pheno_var)){
    tmp_tscore[[i]][[j]] <- res$tscore[[j]]
    tmp_tscore[[i]][[j]]$cohort <- pheno_var[[j]]$cohorts[i]
  }
}

df_tscore_cohorts <- list()

for(j in 1:length(pheno_var)){
  df_tscore_cohorts[[j]] <- do.call(rbind, lapply(tmp_tscore, function(x) x[[j]]))
}

rm(tmp_tscore)


df_tscore_all <- list()
df_pi1 <-  data.frame(pheno_id =sapply(pheno_var, function(x) x$pheno),  tscore=0, pathScore_reactome=0,  pathScore_GO=0, stringsAsFactors = F)

# for each phenotype considered
for(j in 1:length(pheno_var)){
  
  # find only genes/pathway in common (deletion could have happend during pvalue computation)
  #### tscores ####
  
  print(paste('##########', pheno_var[[j]]$pheno, '##########'))
  
  common_genes <- names(which(table(df_tscore_cohorts[[j]]$ensembl_gene_id) == length(pheno_var[[j]]$cohorts)))
  df_tscore_cohorts[[j]] <- df_tscore_cohorts[[j]][df_tscore_cohorts[[j]]$ensembl_gene_id %in% common_genes, ]
  info_tscore <- df_tscore_cohorts[[j]][!duplicated(df_tscore_cohorts[[j]]$ensembl_gene_id), c('ensembl_gene_id', 'external_gene_name','dev_geno', 'test_dev_geno')]
  df_tscore_all[[j]] <- cbind(info_tscore, matrix(0, nrow = nrow(info_tscore), ncol = 9))
  colnames(df_tscore_all[[j]])[-(1:4)] <- paste0(pheno_var[[j]]$pheno, c('_beta', '_se_beta', '_z', '_pval', '_qval', '_pval_BHcorr','_Cochran_stat', '_Cochran_pval', '_model'))
  
  for(i in 1:nrow(info_tscore)){
    
    if(i%%100 == 0){print(paste('gene', i))}
    
    model <- 'fixed'
    tmp <- df_tscore_cohorts[[j]][df_tscore_cohorts[[j]]$ensembl_gene_id %in% info_tscore$ensembl_gene_id[i],]
    w <- 1/(tmp[, sprintf('%s_se_beta', pheno_var[[j]]$pheno)])^2
    beta_all <- sum(tmp[, sprintf('%s_beta', pheno_var[[j]]$pheno)]*w)/sum(w)
    se_all <- sqrt(1/sum(w))
    z_all <- beta_all/se_all
    # p_all <- (1 - pnorm(abs(z_all), 0, 1)) * 2 # rounding approximation! use negative values
    p_all <- 2*pnorm(-abs(z_all), 0,1)
    Q <- sum(w*((beta_all - tmp[, sprintf('%s_beta', pheno_var[[j]]$pheno)])^2))
    Q_pval <- pchisq(Q, df=nrow(tmp)-1, lower.tail=FALSE) # null hypothesis consistency
    
    if(Q_pval<=thr_het){
      
      tau2 <- max(0, (Q-nrow(tmp)+1)/(sum(w) - (sum(w^2)/sum(w))))
      w_new <- 1/(tau2 + (tmp[, sprintf('%s_se_beta', pheno_var[[j]]$pheno)])^2)
      se_all <- sqrt(1/sum(w_new))
      beta_all <- sum(tmp[, sprintf('%s_beta', pheno_var[[j]]$pheno)]*w_new)/sum(w_new)
      z_all <- beta_all/se_all
      # p_all <- (1 - pnorm(abs(z_all), 0, 1)) * 2 # rounding approximation! use negative values
      p_all <- 2*pnorm(-abs(z_all), 0,1)
      model <- 'random'
    }
    
    df_tscore_all[[j]][i, paste0(pheno_var[[j]]$pheno, c('_beta', '_se_beta', '_z', '_pval', '_Cochran_stat', '_Cochran_pval'))] <- c(beta_all, se_all, z_all, p_all, Q, Q_pval)
    df_tscore_all[[j]][i, paste0(pheno_var[[j]]$pheno, '_model')] <- model
    
  }
  
  df_tscore_all[[j]][, paste0(pheno_var[[j]]$pheno, '_qval')] <- qvalue(df_tscore_all[[j]][,paste0(pheno_var[[j]]$pheno, '_pval')])$qvalue
  df_tscore_all[[j]][, paste0(pheno_var[[j]]$pheno, '_pval_BHcorr')] <- p.adjust(df_tscore_all[[j]][,paste0(pheno_var[[j]]$pheno, '_pval')], method = 'BH')
  df_pi1$tscore[j] <- 1- pi0est(df_tscore_all[[j]][,paste0(pheno_var[[j]]$pheno, '_pval')], lambda = 0.5)$pi0
  
}  

final <- list(pheno = res$pheno, tscore = df_tscore_all,,pi1_lambdafixed = df_pi1)

file_name <- ifelse(cov_corr,sprintf('%spval_%s_pheno_covCorr.RData', outFold, phenoName), sprintf('%spval_%s_pheno.RData', outFold, phenoName))
save(final, file = file_name)


