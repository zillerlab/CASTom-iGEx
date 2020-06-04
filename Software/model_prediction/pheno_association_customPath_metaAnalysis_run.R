#### code written by Lucia Trastulla, e-mail: lucia_trastulla@psych.mpg.de ####
# meta - analysis across different cohorts

options(stringsAsFactors=F)
options(max.print=1000)
suppressPackageStartupMessages(library(argparse))
suppressPackageStartupMessages(library(qvalue))


parser <- ArgumentParser(description="Meta analysis tscore and pathway scores (custom)")
parser$add_argument("--pathwayStructure_file", type = "character", help = "custom pathway anntation")
parser$add_argument("--res_cohorts", type="character", nargs = '*', help = "RData file with the phenotype association results")
parser$add_argument("--phenoDatFile_cohorts", type="character", nargs = '*', help = "file with the phenotype info for each cohort")
parser$add_argument("--name_cohort", type="character", nargs = '*', help = "names of the cohorts")
parser$add_argument("--cov_corr", type = "logical",default = F,  help = "if T column in covDat_file are use to correct association (excluded Dx and IDs)")
parser$add_argument("--phenoName", type="character",  help = "phenotype group name considered")
parser$add_argument("--thr_het", type = "double", default = 0.001, help = "threshold for random effect model")
parser$add_argument("--lambda_pi1", type = "double", default = 0.5, help = "lambda parameter to compute pi1")
parser$add_argument("--geneSetName", type = "character", help = "name pathway custom")
parser$add_argument("--abs_tscore", type = "logical", default = F, help = "if true also the pathway using absolute values of tscore is computed")
parser$add_argument("--outFold", type="character", help = "Output file [basename only]")

args <- parser$parse_args()
res_cohorts <- args$res_cohorts
name_cohort <- args$name_cohort
phenoDatFile_cohorts <- args$phenoDatFile_cohorts
cov_corr <- args$cov_corr
phenoName <- args$phenoName
thr_het <- args$thr_het
pathwayStructure_file <- args$pathwayStructure_file
geneSetName <- args$geneSetName
lambda_pi1 <- args$lambda_pi1
abs_tscore <- args$abs_tscore
outFold <- args$outFold

# ##########################################################################################
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
# ###############################################################################################################################

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
print('saved info pheno')


#####################################################################################################
## extraxt pval for each phenotype that have more tha 1 cohort or is not insonsistent
id_pheno <- which(df_pheno_tot$n_studies > 1 & df_pheno_tot$cohort != 'remove_inconsistent')
pheno_var <- lapply(id_pheno, function(x) list(pheno = 0, cohorts = 0))
for(i in 1:length(id_pheno)){
  
  pheno_var[[i]]$pheno <- df_pheno_tot$pheno[id_pheno[i]]
  pheno_var[[i]]$cohorts <- name_cohort[sapply(name_cohort, function(x) grepl(x,df_pheno_tot$cohort[id_pheno[i]]))]
  
}
print(pheno_var)

tmp_tscore <- list()
tmp_path <- list()
tmp_info_path <- list()
tmp_path_abstscore <- list()
tmp_info_path_abstscore <- list()

for(i in 1:length(res_cohorts)){
  res <- get(load(res_cohorts[[i]]))
  print(length(res))
  tmp_tscore[[i]] <- list()
  tmp_path[[i]] <- list()
  tmp_info_path[[i]] <- list()
  tmp_path_abstscore[[i]] <- list()
  tmp_info_path_abstscore[[i]] <- list()
  
  for(j in 1:length(pheno_var)){
    tmp_tscore[[i]][[j]] <- res$tscore[[j]]
    tmp_tscore[[i]][[j]]$cohort <- pheno_var[[j]]$cohorts[i]
    
    tmp_path[[i]][[j]] <- res$pathScore[[j]]
    tmp_path[[i]][[j]]$cohort <- pheno_var[[j]]$cohorts[i]
    
    tmp_info_path[[i]][[j]] <- res$info_pathScore[[j]]
    
    if(abs_tscore){
      tmp_path_abstscore[[i]][[j]] <- res$pathScore_abstscore[[j]]
      tmp_path_abstscore[[i]][[j]]$cohort <- pheno_var[[j]]$cohorts[i]
      tmp_info_path_abstscore[[i]][[j]] <- res$info_pathScore_abstscore[[j]]
    }
  }
  
}


df_tscore_cohorts <- list()
df_path_cohorts <- list()
info_path_tscore <- tmp_info_path[[1]][[1]]
if(abs_tscore){
  df_path_abstscore_cohorts <- list()
  info_path_abstscore_tscore <- tmp_info_path_abstscore[[1]][[1]]
}

for(j in 1:length(pheno_var)){
  df_tscore_cohorts[[j]] <- do.call(rbind, lapply(tmp_tscore, function(x) x[[j]]))  
  df_path_cohorts[[j]] <- do.call(rbind, lapply(tmp_path, function(x) x[[j]]))
  if(abs_tscore){
    df_path_abstscore_cohorts[[j]] <- do.call(rbind, lapply(tmp_path_abstscore, function(x) x[[j]]))
  }  
}
rm(tmp_path)
rm(tmp_tscore)
rm(tmp_info_path)
rm(tmp_path_abstscore)
rm(tmp_info_path_abstscore)
print('data loaded')


# load pathways annotation 
custom_pathway <- get(load(pathwayStructure_file))
path_name <- sapply(custom_pathway, function(x) x$name)

df_tscore_all <- list()
df_path_all <- list()
df_path_abstscore_all <- list()
df_pi1 <-  data.frame(pheno_id = sapply(pheno_var, function(x) x$pheno),  tscore=0, pathScore=0,   stringsAsFactors = F)
if(abs_tscore){df_pi1$pathScore_abstscore <- 0}
info_pathScore <- list()
info_pathScore_abstscore <- list()


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
  
  #### pathway ####
  common_path <- names(which(table(df_path_cohorts[[j]]$path) == length(pheno_var[[j]]$cohorts)))
  df_path_cohorts[[j]] <- df_path_cohorts[[j]][df_path_cohorts[[j]]$path %in% common_path, ]
  info_path <- df_path_cohorts[[j]][!duplicated(df_path_cohorts[[j]]$path), colnames(df_path_cohorts[[j]])[1:7]]
  # take the mean across the cohorts
  info_path$mean_gene_corr <- sapply(info_path$path, function(x) mean(df_path_cohorts[[j]]$mean_gene_corr[df_path_cohorts[[j]]$path %in% x]))
  info_path$sd_gene_corr <- sapply(info_path$path, function(x) mean(df_path_cohorts[[j]]$sd_gene_corr[df_path_cohorts[[j]]$path %in% x]))
  
  df_path_all[[j]] <- data.frame(cbind(info_path, matrix(0, nrow = nrow(info_path), ncol = 9)), stringsAsFactors = F)
  colnames(df_path_all[[j]])[-(1:9)] <- paste0(pheno_var[[j]]$pheno, c('_beta', '_se_beta', '_z', '_pval', '_qval', '_pval_BHcorr', '_Cochran_stat', '_Cochran_pval', '_model'))
  
  for(i in 1:nrow(df_path_all[[j]])){
    
    if(i%%100 == 0){print(paste('path', i))}
    
    model <- 'fixed'
    tmp <- df_path_cohorts[[j]][df_path_cohorts[[j]]$path %in% df_path_all[[j]]$path[i],]
    w <- 1/(tmp[, sprintf('%s_se_beta', pheno_var[[j]]$pheno)])^2
    beta_all <- sum(tmp[, sprintf('%s_beta', pheno_var[[j]]$pheno)]*w)/sum(w)
    se_all <- sqrt(1/sum(w))
    z_all <- beta_all/se_all
    #p_all <- (1 - pnorm(abs(z_all), 0, 1)) * 2
    p_all <- 2*pnorm(-abs(z_all), 0,1)
    Q <- sum(w*((beta_all - tmp[, sprintf('%s_beta', pheno_var[[j]]$pheno)])^2))
    Q_pval <- pchisq(Q, df=nrow(tmp)-1, lower.tail=FALSE) # null hypothesis consistency
    
    if(Q_pval<=thr_het){
      
      tau2 <- max(0, (Q-nrow(tmp)+1)/(sum(w) - (sum(w^2)/sum(w))))
      w_new <- 1/(tau2 + (tmp[, sprintf('%s_se_beta', pheno_var[[j]]$pheno)])^2)
      se_all <- sqrt(1/sum(w_new))
      beta_all <- sum(tmp[, sprintf('%s_beta', pheno_var[[j]]$pheno)]*w_new)/sum(w_new)
      z_all <- beta_all/se_all
      # p_all <- (1 - pnorm(abs(z_all), 0, 1)) * 2
      p_all <- 2*pnorm(-abs(z_all), 0,1)
      model <- 'random'
    }
    
    df_path_all[[j]][i, paste0(pheno_var[[j]]$pheno, c('_beta', '_se_beta', '_z', '_pval', '_Cochran_stat', '_Cochran_pval'))] <- c(beta_all, se_all, z_all, p_all,  Q, Q_pval)
    df_path_all[[j]][i, paste0(pheno_var[[j]]$pheno, '_model')] <- model
    
  }
  
  df_path_all[[j]][, paste0(pheno_var[[j]]$pheno, '_qval')] <- qvalue(df_path_all[[j]][,paste0(pheno_var[[j]]$pheno, '_pval')])$qvalue
  df_path_all[[j]][, paste0(pheno_var[[j]]$pheno, '_pval_BHcorr')] <- p.adjust(df_path_all[[j]][,paste0(pheno_var[[j]]$pheno, '_pval')], method = 'BH')
  df_pi1$pathScore[j] <- 1- pi0est(df_path_all[[j]][,paste0(pheno_var[[j]]$pheno, '_pval')], lambda = 0.5)$pi0
  
  # annotate
  info_pathScore[[j]] <- vector(mode = 'list', length = nrow(df_path_all[[j]]))
  
  for(i in 1:nrow(df_path_all[[j]])){
    
    if(i%%100 == 0){print(paste('path', i))}
    info_pathScore[[j]][[i]]$path <- df_path_all[[j]][df_path_all[[j]]$path %in% info_path_tscore[[i]]$path$path,]
    info_pathScore[[j]][[i]]$genes_path <- info_path_tscore[[i]]$genes_path
    info_pathScore[[j]][[i]]$tscore <- df_tscore_all[[j]][df_tscore_all[[j]]$ensembl_gene_id %in% info_path_tscore[[i]]$tscore$ensembl_gene_id,]
    
  }
  
  if(abs_tscore){
    
    #### pathway abstscore ####
    common_path <- names(which(table(df_path_abstscore_cohorts[[j]]$path) == length(pheno_var[[j]]$cohorts)))
    df_path_abstscore_cohorts[[j]] <- df_path_abstscore_cohorts[[j]][df_path_abstscore_cohorts[[j]]$path %in% common_path, ]
    info_path_abstscore <- df_path_abstscore_cohorts[[j]][!duplicated(df_path_abstscore_cohorts[[j]]$path), colnames(df_path_abstscore_cohorts[[j]])[1:7]]
    # take the mean across the cohorts
    info_path_abstscore$mean_gene_corr <- sapply(info_path_abstscore$path, function(x) mean(df_path_abstscore_cohorts[[j]]$mean_gene_corr[df_path_abstscore_cohorts[[j]]$path %in% x]))
    info_path_abstscore$sd_gene_corr <- sapply(info_path_abstscore$path, function(x) mean(df_path_abstscore_cohorts[[j]]$sd_gene_corr[df_path_abstscore_cohorts[[j]]$path %in% x]))
    
    df_path_abstscore_all[[j]] <- data.frame(cbind(info_path_abstscore, matrix(0, nrow = nrow(info_path_abstscore), ncol = 9)), stringsAsFactors = F)
    colnames(df_path_abstscore_all[[j]])[-(1:9)] <- paste0(pheno_var[[j]]$pheno, c('_beta', '_se_beta', '_z', '_pval', '_qval', '_pval_BHcorr', '_Cochran_stat', '_Cochran_pval', '_model'))
    
    for(i in 1:nrow(df_path_abstscore_all[[j]])){
      
      if(i%%100 == 0){print(paste('path', i))}
      
      model <- 'fixed'
      tmp <- df_path_abstscore_cohorts[[j]][df_path_abstscore_cohorts[[j]]$path %in% df_path_abstscore_all[[j]]$path[i],]
      w <- 1/(tmp[, sprintf('%s_se_beta', pheno_var[[j]]$pheno)])^2
      beta_all <- sum(tmp[, sprintf('%s_beta', pheno_var[[j]]$pheno)]*w)/sum(w)
      se_all <- sqrt(1/sum(w))
      z_all <- beta_all/se_all
      #p_all <- (1 - pnorm(abs(z_all), 0, 1)) * 2
      p_all <- 2*pnorm(-abs(z_all), 0,1)
      Q <- sum(w*((beta_all - tmp[, sprintf('%s_beta', pheno_var[[j]]$pheno)])^2))
      Q_pval <- pchisq(Q, df=nrow(tmp)-1, lower.tail=FALSE) # null hypothesis consistency
      
      if(Q_pval<=thr_het){
        
        tau2 <- max(0, (Q-nrow(tmp)+1)/(sum(w) - (sum(w^2)/sum(w))))
        w_new <- 1/(tau2 + (tmp[, sprintf('%s_se_beta', pheno_var[[j]]$pheno)])^2)
        se_all <- sqrt(1/sum(w_new))
        beta_all <- sum(tmp[, sprintf('%s_beta', pheno_var[[j]]$pheno)]*w_new)/sum(w_new)
        z_all <- beta_all/se_all
        # p_all <- (1 - pnorm(abs(z_all), 0, 1)) * 2
        p_all <- 2*pnorm(-abs(z_all), 0,1)
        model <- 'random'
      }
      
      df_path_abstscore_all[[j]][i, paste0(pheno_var[[j]]$pheno, c('_beta', '_se_beta', '_z', '_pval', '_Cochran_stat', '_Cochran_pval'))] <- c(beta_all, se_all, z_all, p_all,  Q, Q_pval)
      df_path_abstscore_all[[j]][i, paste0(pheno_var[[j]]$pheno, '_model')] <- model
      
    }
    
    df_path_abstscore_all[[j]][, paste0(pheno_var[[j]]$pheno, '_qval')] <- qvalue(df_path_abstscore_all[[j]][,paste0(pheno_var[[j]]$pheno, '_pval')])$qvalue
    df_path_abstscore_all[[j]][, paste0(pheno_var[[j]]$pheno, '_pval_BHcorr')] <- p.adjust(df_path_abstscore_all[[j]][,paste0(pheno_var[[j]]$pheno, '_pval')], method = 'BH')
    df_pi1$pathScore_abstscore[j] <- 1- pi0est(df_path_abstscore_all[[j]][,paste0(pheno_var[[j]]$pheno, '_pval')], lambda = 0.5)$pi0
    
    # annotate
    info_pathScore_abstscore[[j]] <- vector(mode = 'list', length = nrow(df_path_abstscore_all[[j]]))
    
    for(i in 1:nrow(df_path_abstscore_all[[j]])){
      
      if(i%%100 == 0){print(paste('path', i))}
      info_pathScore_abstscore[[j]][[i]]$path <- df_path_abstscore_all[[j]][df_path_abstscore_all[[j]]$path %in% info_path_abstscore_tscore[[i]]$path$path,]
      info_pathScore_abstscore[[j]][[i]]$genes_path <- info_path_abstscore_tscore[[i]]$genes_path
      info_pathScore_abstscore[[j]][[i]]$tscore <- df_tscore_all[[j]][df_tscore_all[[j]]$ensembl_gene_id %in% info_path_abstscore_tscore[[i]]$tscore$ensembl_gene_id,]
      
    }
  }
  
}  

final <- list(pheno = res$pheno, tscore = df_tscore_all, pathScore = df_path_all, 
              pi1_lambdafixed = df_pi1, info_pathScore = info_pathScore)
if(abs_tscore){
  final$pathScore_abstscore <- df_path_abstscore_all
  final$info_pathScore_abstscore <- info_pathScore_abstscore
}

file_name <- ifelse(cov_corr, sprintf('%spval_%s_pheno_covCorr_customPath_%s.RData', outFold, phenoName, geneSetName), sprintf('%spval_%s_pheno_customPath_%s.RData', outFold, phenoName, geneSetName))
save(final, file = file_name)


