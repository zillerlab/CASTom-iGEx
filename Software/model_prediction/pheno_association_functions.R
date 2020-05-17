# function used in AssociationAnalysis with bigmemory matrix
# parallelization over genes (general for all matrix type)
compute_reg_pheno <- function(id_gene, type_pheno, id_pheno, total_mat){
  
  var <- paste0('X',id_gene)
  fmla <- paste0('p', id_pheno, " ~ ", var)
  if(cov_corr){fmla <- paste(fmla, '+' , paste(paste0('c', cov_names), collapse = '+'))}
  fmla <- as.formula(fmla)
  
  if(type_pheno == 'CONTINUOUS'){
    
    res <- glm(fmla, data = total_mat, family = 'gaussian')
    output <- coef(summary(res))[rownames(coef(summary(res))) == var,1:4]
  }
  
  if(type_pheno %in% c('CAT_SINGLE_UNORDERED', 'CAT_SINGLE_BINARY', 'CAT_MUL_BINARY_VAR')){

      if(!all(unique(na.omit(total_mat[, paste0('p', id_pheno)])) %in% c(0,1))){
      
        min_id <- which(total_mat[, paste0('p', id_pheno)] == min(total_mat[, paste0('p', id_pheno)], na.rm = T))
        total_mat[min_id, paste0('p', id_pheno)] <- 0
        max_id <- which(total_mat[, paste0('p', id_pheno)] == max(total_mat[, paste0('p', id_pheno)], na.rm = T))
        total_mat[max_id, paste0('p', id_pheno)] <- 1
      
      }
    
      res <- glm(fmla, data = total_mat, family = 'binomial')
      output <- coef(summary(res))[rownames(coef(summary(res))) == var,1:4]
      
    }
  
  if(type_pheno == 'CAT_ORD'){
    
    total_mat[,paste0('p', id_pheno)] <- as.factor(total_mat[,paste0('p', id_pheno)])
    
    res <- polr(fmla, data = total_mat, Hess=TRUE)
    if(!any(is.na(res$Hess))){
      ct <- coeftest(res)  
      output <- ct[rownames(ct) == var,1:4]
    }else{
      output <- rep(NA, 4)
    }
    
  }
  
  if(! type_pheno %in% c('CAT_ORD', 'CAT_SINGLE_UNORDERED', 'CAT_SINGLE_BINARY', 'CAT_MUL_BINARY_VAR', 'CONTINUOUS')){
    output <- 'wrong pheno type annotation'
  }
  
  return(output)
  
}






# parallelization over phenotypes

compute_reg_pheno_tscore <- function(j){
  
  df <- cbind(tmp_tscore, matrix(-1, ncol = 4, nrow = nrow(tmp_tscore)))
  colnames(df)[-(1:ncol(tmp_tscore))] <-  names_df[j,]
  
  if(phenoAnn_tmp$transformed_type[j] == 'CONTINUOUS'){
    
    for(i in 1:nrow(df)){
      
      if(i%%100 == 0){print(paste('gene', i))}
      # print(i)
      var <- paste0('X',i)
      fmla <- paste0('p', phenoAnn_tmp$pheno_id[j], " ~ ", var)
      if(cov_corr){fmla <- paste(fmla, '+' , paste(cov_names, collapse = '+'))}
      fmla <- as.formula(fmla)
      res <- glm(fmla, data = tscoreMat_association, family = 'gaussian')
      df[i,-(1:ncol(tmp_tscore))] <- coef(summary(res))[rownames(coef(summary(res))) == var,1:4]
      
    }
    
  }
  
  if(phenoAnn_tmp$transformed_type[j] %in% c('CAT_SINGLE_UNORDERED', 'CAT_SINGLE_BINARY', 'CAT_MUL_BINARY_VAR')){
    
    if(!all(unique(na.omit(tscoreMat_association[, paste0('p', phenoAnn_tmp$pheno_id[j])])) %in% c(0,1))){
      
      min_id <- which(tscoreMat_association[, paste0('p', phenoAnn_tmp$pheno_id[j])] == min(tscoreMat_association[, paste0('p', phenoAnn_tmp$pheno_id[j])], na.rm = T))
      tscoreMat_association[min_id, paste0('p', phenoAnn_tmp$pheno_id[j])] <- 0
      max_id <- which(tscoreMat_association[, paste0('p', phenoAnn_tmp$pheno_id[j])] == max(tscoreMat_association[, paste0('p', phenoAnn_tmp$pheno_id[j])], na.rm = T))
      tscoreMat_association[max_id, paste0('p', phenoAnn_tmp$pheno_id[j])] <- 1
      
    }
    
    for(i in 1:nrow(df)){
      
      if(i%%100 == 0){print(paste('gene', i))}
      # print(i)
      var <- paste0('X',i)
      fmla <- paste0('p', phenoAnn_tmp$pheno_id[j], " ~ ", var)
      if(cov_corr){fmla <- paste(fmla, '+' , paste(cov_names, collapse = '+'))}
      fmla <- as.formula(fmla)
      res <- glm(fmla, data = tscoreMat_association, family = 'binomial')
      df[i,-(1:ncol(tmp_tscore))] <- coef(summary(res))[rownames(coef(summary(res))) == var,1:4]
      
    }
    
  }
  
  if(phenoAnn_tmp$transformed_type[j] == 'CAT_ORD'){
    tscoreMat_association[,paste0('p', phenoAnn_tmp$pheno_id[j])] <- as.factor(tscoreMat_association[,paste0('p', phenoAnn_tmp$pheno_id[j])])
    
    for(i in 1:nrow(df)){
      
      if(i%%100 == 0){print(paste('gene', i))}
      # print(i)
      var <- paste0('X',i)
      fmla <- paste0('p', phenoAnn_tmp$pheno_id[j], " ~ ", var)
      if(cov_corr){fmla <- paste(fmla, '+' , paste(cov_names, collapse = '+'))}
      fmla <- as.formula(fmla)
      res <- polr(fmla, data = tscoreMat_association, Hess=TRUE)
      ct <- coeftest(res)
      df[i,-(1:ncol(tmp_tscore))] <- ct[rownames(ct) == var,1:4]
      
    }
    
  }
  
  qval <- qvalue(as.vector(df[, names_df[j,4]]))
  df <- cbind(df, qval$qvalue)
  colnames(df)[ncol(df)] <-  names_df_qval[j]

  pi1 <- 1 - qval$pi0
  
  
  return(list(df = df, pi1 = pi1))
  
  
}


compute_reg_pheno_pathR <- function(j){
  
  df <- cbind(tmp_pathScore_R, matrix(-1, ncol = 4, nrow = nrow(tmp_pathScore_R)))
  colnames(df)[-(1:ncol(tmp_pathScore_R))] <-  names_df[j,]

  if(phenoAnn_tmp$transformed_type[j] == 'CONTINUOUS'){
    
    for(i in 1:nrow(df)){
      
      if(i%%100 == 0){print(paste('gene', i))}
      # print(i)
      var <- paste0('X',i)
      fmla <- paste0('p', phenoAnn_tmp$pheno_id[j], " ~ ", var)
      if(cov_corr){fmla <- paste(fmla, '+' , paste(cov_names, collapse = '+'))}
      fmla <- as.formula(fmla)
      res <- glm(fmla, data = pathScore_r_association, family = 'gaussian')
      df[i,-(1:ncol(tmp_pathScore_R))] <- coef(summary(res))[rownames(coef(summary(res))) == var,1:4]
      
    }
    
  }
  
  if(phenoAnn_tmp$transformed_type[j] %in% c('CAT_SINGLE_UNORDERED', 'CAT_SINGLE_BINARY', 'CAT_MUL_BINARY_VAR')){
    
    if(!all(unique(na.omit(pathScore_r_association[, paste0('p', phenoAnn_tmp$pheno_id[j])])) %in% c(0,1))){
      
      min_id <- which(pathScore_r_association[, paste0('p', phenoAnn_tmp$pheno_id[j])] == min(pathScore_r_association[, paste0('p', phenoAnn_tmp$pheno_id[j])], na.rm = T))
      pathScore_r_association[min_id, paste0('p', phenoAnn_tmp$pheno_id[j])] <- 0
      max_id <- which(pathScore_r_association[, paste0('p', phenoAnn_tmp$pheno_id[j])] == max(pathScore_r_association[, paste0('p', phenoAnn_tmp$pheno_id[j])], na.rm = T))
      pathScore_r_association[max_id, paste0('p', phenoAnn_tmp$pheno_id[j])] <- 1
      
    }
    
    for(i in 1:nrow(df)){
      
      if(i%%100 == 0){print(paste('gene', i))}
      # print(i)
      var <- paste0('X',i)
      fmla <- paste0('p', phenoAnn_tmp$pheno_id[j], " ~ ", var)
      if(cov_corr){fmla <- paste(fmla, '+' , paste(cov_names, collapse = '+'))}
      fmla <- as.formula(fmla)
      res <- glm(fmla, data = pathScore_r_association, family = 'binomial')
      df[i,-(1:ncol(tmp_pathScore_R))] <- coef(summary(res))[rownames(coef(summary(res))) == var,1:4]
      
    }
    
  }
  
  if(phenoAnn_tmp$transformed_type[j] == 'CAT_ORD'){
    pathScore_r_association[,paste0('p', phenoAnn_tmp$pheno_id[j])] <- as.factor(pathScore_r_association[,paste0('p', phenoAnn_tmp$pheno_id[j])])
    
    for(i in 1:nrow(df)){
      
      if(i%%100 == 0){print(paste('gene', i))}
      # print(i)
      var <- paste0('X',i)
      fmla <- paste0('p', phenoAnn_tmp$pheno_id[j], " ~ ", var)
      if(cov_corr){fmla <- paste(fmla, '+' , paste(cov_names, collapse = '+'))}
      fmla <- as.formula(fmla)
      res <- polr(fmla, data = pathScore_r_association, Hess=TRUE)
      ct <- coeftest(res)
      df[i,-(1:ncol(tmp_pathScore_R))] <- ct[rownames(ct) == var,1:4]
      
    }
    
  }
  
  qval <- qvalue(as.vector(df[, names_df[j,4]]))
  df <- cbind(df, qval$qvalue)
  colnames(df)[ncol(df)] <-  names_df_qval[j]

  pi1 <- 1 - qval$pi0
  # pi1=1
  
  return(list(df = df, pi1 = pi1))
  
}
