# functions called in Priler_ (part1 to part4 + plots) scripts
### NOTE:
### prior baseline (1 instead of 0.5)
### regression uses intercept, y and x are not centered

# combine function for parallelization
comb <- function(x, ...) {
  lapply(seq_along(x),
         function(i) c(x[[i]], lapply(list(...), function(y) y[[i]])))
}

# obtains prior from weights and prior matrix data
getPrior <- function(priorMat,betas){
  r=rowSums(t(t(priorMat)*betas))
  f=1/(1+exp(-r))
  return(2*(1-f))
}

# generate CV runs (stratified and repeated options)
generateCVRuns <- function (labels, ntimes = 10, nfold = 10, leaveOneOut = FALSE, stratified = FALSE){
  
  if (leaveOneOut) 
    ntimes <- 1
  numSamples <- length(labels)
  res <- lapply(1:ntimes, function(run) {
    if (leaveOneOut) {
      indices <- as.list(1:numSamples)
    }
    else {
      if (stratified) {
        classes <- unique(labels)
        sing.perm <- lapply(classes, function(cl) {
          index <- which(labels == cl)
          sample(index, length(index))
        })
        permut <- unlist(sing.perm)
        indices <- lapply(1:nfold, function(i) {
          c()
        })
        for (i in 1:numSamples) {
          k = i%%nfold
          if (k == 0) 
            k = nfold
          indices[[k]] <- c(indices[[k]], permut[i])
        }
      }
      else {
        permut <- sample(1:numSamples, numSamples, replace = FALSE)
        indices <- lapply(1:nfold, function(i) {
          permut[seq(i, numSamples, nfold)]
        })
      }
    }
    names(indices) <- paste("Fold ", 1:nfold)
    return(indices)
  })
  names(res) <- paste("Run ", 1:ntimes)
  return(res)
}


#### used in PriLer_part1_run.R ####
# find optimal lambda and alpha for each gene via CV, use all the samples, no prior weights,
# NOTE: assigning penalty.factor the same values is the same as not using it
expPrediction_cv_noPrior_fold <- function(X, alpha, fold, seed, nfolds){
  
  ind_SNPs <- geneSnpDist[,X]!=0 
  nSnp <- sum(ind_SNPs)
  
  if (nSnp>2){
    
    d <- which(ind_SNPs)[1]-1
    genotype <- as.matrix(cbind(genDat[fold,ind_SNPs],covDat[fold,]))
    expressionValue <- as.numeric(expDat[fold,X])
    
    # compute deviance for the model using only covariates
    cov_mod <- cbind(expressionValue, covDat[fold,])
    fit_cov <- lm(cov_mod[,1]~as.matrix(cov_mod[,-1]))
    dev_lmcov <- summary(fit_cov)$r.squared
    
    # compute lambda ranges based on alpha 
    lambda_set <- vector(mode  = 'list', length = length(alpha))
    for(i in 1:length(alpha)){
      # print(i)
      lambda_set[[i]] <- glmnet(x = genotype, y = expressionValue,  alpha = alpha[i], nlambda = 100, intercept = T, standardize = F, 
                                penalty.factor = c(rep(1, length(which(ind_SNPs))), rep(0, ncol(covDat))))$lambda
      
    }
    
    cv_res <- vector(mode = 'list', length = length(alpha)) 
    # create inner folder
    set.seed(seed)
    folds_in <- generateCVRuns(nfold = nfolds, labels = 1:nrow(genotype), ntimes = 1)[[1]]
    foldid_vect <- unlist(mapply(function(x, y) rep(x,length(y)), x =1:nfolds, y = folds_in, SIMPLIFY = F))
    foldid_vect <-  foldid_vect[order(unlist(folds_in))]
    
    for(i in 1:length(alpha)){
      
      set.seed(seed)
      print(paste('alpha val', alpha[i]))
      # some values could be constant for a fold, remove gene
      cv_res[[i]] <- tryCatch(cv.glmnet(x = genotype, y = expressionValue, alpha = alpha[i], lambda = lambda_set[[i]], foldid = foldid_vect, intercept = T, standardize = F, 
                                        penalty.factor = c(rep(1, length(which(ind_SNPs))), rep(0, ncol(covDat)))), 
                              error=function(...) list(lambda.min=NA)) 
      
    }
    if(is.na(cv_res[[1]]$lambda.min)){
      return(c(X,NA,NA,NA,NA,NA,NA,NA,NA,NA, NA, NA))
    }else{
      
      # keep only the best (lambda,alpha) pair in terms of cvm, 
      id_lambda <- sapply(cv_res, function(x) which(x$lambda == x$lambda.min))
      cvm_min <- sapply(1:length(alpha), function(x) cv_res[[x]]$cvm[id_lambda[x]])
      
      id_alpha <- which.min(cvm_min)
      alphamin <- alpha[id_alpha]
      
      # reassign keeping only the best lambda value
      res <- cv_res[[id_alpha]]
      
      ind <- which.min(res$cvm)
      selLambda <- res$lambda.min
      selAlpha <- alpha[[id_alpha]]
      betas <- res$glmnet.fit$beta[,ind]
      selIntercept <- res$glmnet.fit$a0[ind]
      r <- as.numeric(which(betas!=0))
      selBeta <- as.numeric(betas)[r]
      sInd <- r>nSnp
      v <- r[sInd]-nSnp
      r[r<=nSnp] <- r[r<=nSnp]+d
      r[sInd] <- covIndex[v]
      r <- c(r, nrow(geneSnpDist)+ncol(covDat)+1)
      r <- r[!is.na(r)]
      
      P <- length(which(ind_SNPs))
      pred_geno <- genDat[fold,ind_SNPs] %*% betas[1:P]
      pred_cov <- as.matrix(covDat[fold, ]) %*% betas[(P+1):length(betas)] + selIntercept
      expressionValue_geno <- expressionValue - selIntercept - as.matrix(covDat[fold, ]) %*% betas[(P+1):length(betas)]
      # dev_geno_old <- 1-sum((expressionValue_geno-pred_geno)^2)/sum((expressionValue_geno - mean(expressionValue_geno))^2)
      dev_lmgeno <- res$glmnet.fit$dev.ratio[ind] - dev_lmcov
      
      #pred <- genotype %*% betas + selIntercept
      # dev <- sum((pred-mean(expressionValue))^2)/sum((expressionValue - mean(expressionValue))^2) + 
      #        2*sum((expressionValue-pred)*(pred - mean(expressionValue)))/sum((expressionValue - mean(expressionValue))^2)
      
      dev_geno <- sum((pred_geno-mean(expressionValue_geno))^2 + 2*(expressionValue_geno - pred_geno)*(pred_geno - mean(expressionValue_geno)))/sum((expressionValue - mean(expressionValue))^2)
      dev_cov <- sum((pred_cov-mean(pred_cov))^2)/sum((expressionValue - mean(expressionValue))^2)
      dev_geno_cov <- sum(2*(expressionValue_geno - mean(expressionValue_geno))*(pred_cov - mean(pred_cov)))/sum((expressionValue - mean(expressionValue))^2)
      
      # correlation true vs predicted
      # consider the dataset without covariate
      cor_est <- cor.test(expressionValue_geno,as.numeric(pred_geno))$estimate
      cor_pval <- cor.test(expressionValue_geno,as.numeric(pred_geno))$p.value
      # consider the original dataset
      cor_est_noadj <- cor.test(expressionValue,as.numeric(pred_geno))$estimate
      cor_pval_noadj <- cor.test(expressionValue,as.numeric(pred_geno))$p.value
      
      return(cbind(X,r,c(selBeta,as.numeric(selIntercept)),selLambda,selAlpha, res$glmnet.fit$dev.ratio[ind], dev_geno, dev_cov, dev_geno_cov, dev_lmgeno, cor_est, cor_pval, cor_est_noadj, cor_pval_noadj))
    }
  }else{
    return(c(X,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA))
  }
}


#### used in PriLer_part1_run.R and PriLer_part2_run.R and PriLer_part4_run.R####
# function to compute the deviance on the test set
deviance_funct <- function(genDat, expDat, covDat, modelMatrix, FoldK){
  
  pExp_geno <- as.matrix(t(modelMatrix[1:ncol(genDat),]) %*% t(genDat[FoldK, ]))  
  pExp_cov <- as.matrix(t(modelMatrix[c((ncol(genDat)+1):nrow(modelMatrix)),]) %*% t(cbind(covDat[FoldK,], rep(1,length(FoldK))))) # include intercept
  # pExp_cov_mFold <- as.matrix(t(modelMatrix[c((ncol(genDat)+1):nrow(modelMatrix)),]) %*% t(cbind(ancMat[-FoldK,], rep(1,nrow(ancMat[-FoldK,])))))
  # as.matrix(t(modelMatrix[c((ncol(genDat)+1):(nrow(modelMatrix)-1)),]) %*% t(ancMat[FoldK,]))
  pExp <- pExp_geno + pExp_cov
  
  # compute total deviance on test set
  dev <- rowSums((pExp-t(expDat[FoldK,]))^2)
  dev_null <- rowSums((t(expDat[FoldK, ]) - colMeans(expDat[FoldK,]))^2) # correct dimensions
  dev.ratio <- 1-dev/dev_null
  names(dev.ratio) <- NULL
  
  # compute deviance explained by geno
  expDat_adj <- t(expDat[FoldK,]) - pExp_cov
  
  dev_geno <- rowSums((pExp_geno - rowMeans(expDat_adj))^2 + 2*(expDat_adj - pExp_geno)*(pExp_geno - rowMeans(expDat_adj)))/dev_null
  dev_cov <- rowSums((pExp_cov-rowMeans(pExp_cov))^2)/dev_null
  dev_geno_cov <- rowSums(2*(expDat_adj - rowMeans(expDat_adj))*(pExp_cov - rowMeans(pExp_cov)))/dev_null
  
  # expDat_adj_mFold <- t(expDat[-FoldK,]) - pExp_cov_mFold
  # dev_geno <- rowSums((pExp_geno-expDat_adj)^2)
  # dev_null_geno <- rowSums((expDat_adj - rowMeans(expDat_adj_mFold))^2) 
  # dev.ratio_geno_old <- 1-dev_geno/dev_null_geno
  # names(dev.ratio_geno_old) <- NULL
  
  # # compute deviance explained by geno, new version: covariance regression alone
  # dev_cov <- c()
  # for(i in 1:ncol(expDat)){
  #   cov_mod <- cbind(data.frame(y = expDat[-FoldK,i]), ancMat[-FoldK,])
  #   fmla <- "y ~ "
  #   for(l in 2:(ncol(cov_mod)-1)){fmla <- paste0(fmla, colnames(cov_mod)[l], '+')}
  #   fmla <- as.formula(paste0(fmla, colnames(cov_mod)[ncol(cov_mod)]))
  #   fit_cov <- lm(fmla, data = cov_mod)
  #   pred_cov <- predict.lm(object = fit_cov, newdata = ancMat[FoldK,])
  #   dev_cov[i] <-  1 - sum((pred_cov-expDat[FoldK,i])^2)/sum((expDat[FoldK,i] - mean(expDat[-FoldK,i]))^2) 
  # }
  
  # dev.ratio_geno <- ifelse(dev_cov <= 0, dev.ratio,dev.ratio - dev_cov) 
  
  cor_est <- c()
  cor_pval <- c()
  for(g in 1:ncol(expDat)){
    cor_est[g] <- cor.test(expDat_adj[g,],as.numeric(pExp_geno[g,]))$estimate
    cor_pval[g] <- cor.test(expDat_adj[g,],as.numeric(pExp_geno[g,]))$p.value
  }
  
  cor_est_noadj <- c()
  cor_pval_noadj <- c()
  for(g in 1:ncol(expDat)){
    cor_est_noadj[g] <- cor.test(expDat[FoldK,g],as.numeric(pExp_geno[g,]))$estimate
    cor_pval_noadj[g] <- cor.test(expDat[FoldK,g],as.numeric(pExp_geno[g,]))$p.value
  }
  
  return(cbind(dev.ratio, dev_geno, dev_cov, dev_geno_cov, cor_est, cor_pval, cor_est_noadj, cor_pval_noadj))
  
}

#### used in PriLer_part2_run.R ####
# compute regression for fixed alpha and lambda parameter for each chromosome
expPrediction_chr <- function(expDat_k, genDat_k, covDat_k, prior, d, lambda, alpha, chr, X){
  
  nSnp <- sum(ncol(genDat_k))
  
  if (nSnp>2){
    
    # # ######
    # cons <- setdiff(which(ind_SNPs)[1]:which(ind_SNPs)[length(which(ind_SNPs))], which(ind_SNPs))
    # if(length(cons)>0){
    #   print(paste0('index problem:', cons))
    # }
    # # use to check the SNPs id are consecutive
    # # ######
    
    # d=which(ind_SNPs)[1]-1
    genotype <- as.matrix(cbind(genDat_k, covDat_k))
    expressionValue <- as.numeric(expDat_k)
    
    # compute deviance for the model using only covariates
    cov_mod <- cbind(expressionValue, covDat_k)
    fit_cov <- lm(cov_mod[,1]~as.matrix(cov_mod[,-1]))
    dev_lmcov <- summary(fit_cov)$r.squared
    
    new_prior <- c(prior, rep(1,ncol(covDat_k)))
    # insert weights on the X matrix, so that they are not automatically normalized by glmnet
    new_genotype <- sweep(genotype, 2, new_prior, '/')
    
    # lambdaMax <- max(abs(expressionValue %*% new_genotype)/(nrow(new_genotype)*alphaVec[X]))
    res <- glmnet(new_genotype, expressionValue, alpha=alpha, lambda=lambda, standardize=FALSE, intercept = TRUE, 
                  penalty.factor = c(rep(1, nSnp), rep(0, ncol(covDat_k)))) 
    
    beta <- res$beta/new_prior
    selLambda <- lambda
    selAlpha <- alpha
    selIntercept <- res$a0
    r <- as.numeric(which(beta!=0))
    selBeta <- as.numeric(beta)[r]
    sInd <- r>nSnp
    v <- r[sInd]-nSnp
    r[r<=nSnp] <- r[r<=nSnp]+d
    r[sInd] <- covIndex_chr[,chr][v]
    r <- c(r, nrow(geneSnpDist[[chr]])+ncol(covDat_k)+1)
    r <- r[!is.na(r)]
    
    pred <- genotype %*% beta + selIntercept
    SquErr <- mean((expressionValue-pred)^2)
    
    P <- nSnp
    pred_geno <- genDat_k %*% beta[1:P]
    pred_cov <- as.matrix(covDat_k) %*% beta[(P+1):length(beta)] + selIntercept
    expressionValue_geno <- expressionValue - pred_cov
    # dev_geno_old <- 1-sum((expressionValue_geno-pred_geno)^2)/sum((expressionValue_geno - mean(expressionValue_geno))^2)
    
    dev_geno <- sum((pred_geno-mean(expressionValue_geno))^2 + 2*(expressionValue_geno - pred_geno)*(pred_geno - mean(expressionValue_geno)))/sum((expressionValue - mean(expressionValue))^2)
    dev_cov <- sum((pred_cov-mean(pred_cov))^2)/sum((expressionValue - mean(expressionValue))^2)
    dev_geno_cov <- sum(2*(expressionValue_geno - mean(expressionValue_geno))*(pred_cov - mean(pred_cov)))/sum((expressionValue - mean(expressionValue))^2)
    
    dev_lmgeno <- res$dev.ratio - dev_lmcov
    
    cor_est <- cor.test(expressionValue_geno,as.numeric(pred_geno))$estimate
    cor_pval <- cor.test(expressionValue_geno,as.numeric(pred_geno))$p.value
    # consider the original dataset
    cor_est_noadj <- cor.test(expressionValue,as.numeric(pred_geno))$estimate
    cor_pval_noadj <- cor.test(expressionValue,as.numeric(pred_geno))$p.value
    
    # return(cbind(X,r,c(selBeta,as.numeric(selIntercept)),selLambda,selAlpha,res$dev.ratio))
    return(cbind(X,r,c(selBeta,as.numeric(selIntercept)),selLambda,selAlpha, SquErr, res$dev.ratio, dev_geno, dev_cov, dev_geno_cov, dev_lmgeno, cor_est, cor_pval, cor_est_noadj, cor_pval_noadj))
    
  }else{
    
    
    # return(c(X,NA,NA,NA,NA,NA))
    return(c(X,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA))
  }
  
}

#### used in PriLer_part2_run.R and Priler_part3_run.R ####
# objective function for the weight optimization
getStepFunction_chr <- function(pWeight, modelMatrix, priorMat, lambdas, alphas, E_h){
  
  pen <- c()
  nchr <- length(modelMatrix)
  for(j in 1:nchr){
    
    # print(j)
    
    beta_mat <- modelMatrix[[j]][1:nrow(priorMat[[j]]),] # dim PxN
    
    p <- getPrior(priorMat[[j]], pWeight)
    
    # penalization
    temp1 <- colSums((beta_mat^2)*p) 
    temp2 <-  colSums(abs(beta_mat)*p)
    pen[j] <- sum(lambdas[[j]]*(0.5*(1-alphas[[j]])*temp1 + alphas[[j]] * temp2))  
    
  }
  
  obj <- sum(pen) + E_h*sum(pWeight^2)
  return(obj)
}

#### used in PriLer_part2_run.R and PriLer_part3_run.R ####
# gradient function for the weight optimization
getGradient_chr <- function(pWeight, modelMatrix, priorMat, lambdas, alphas, E_h){
  
  nchr <- length(modelMatrix)
  grad <- matrix(0, ncol =length(pWeight), nrow = nchr)
  
  for(j in 1:nchr){
    
    # print(j)
    
    beta_mat <- modelMatrix[[j]][1:nrow(priorMat[[j]]),] # dim PxN
    p <- getPrior(priorMat[[j]],pWeight)
    
    for (k in 1:length(pWeight)){
      
      pOrig <- priorMat[[j]][,k]
      
      temp1 <- colSums((beta_mat^2) * p * pOrig * (0.5*p-1)) #/2 added because the prior is multiplied by 2 in get prior formula
      temp2 <-  colSums(abs(beta_mat) * p * pOrig * (0.5*p-1))
      
      
      grad[j,k] <- sum(lambdas[[j]] * (0.5*(1-alphas[[j]])*temp1 + alphas[[j]] * temp2))
      
    }
    
  }
  
  grad <- colSums(grad) + 2*E_h*pWeight
  
  return(grad)
}


#### used in PriLer_part2_run.R and PriLer_part3_run.R ####
# obtain the Error componenets, can be used on both train and test set
getErrorComponents_chr <- function(lambdas, genDat, expDat, covDat, modelMatrix, pWeight, priorMat, alphas, E_h, fold, all_Chroms, gene_ann){
  
  nchr <- length(all_Chroms)
  obj <- c()
  pen <- c()
  
  for(j in 1:nchr){
    
    beta_mat <- modelMatrix[[j]][1:nrow(priorMat[[j]]),] # dim PxN
    beta_cov_mat <- modelMatrix[[j]][c((nrow(priorMat[[j]])+1):nrow(modelMatrix[[j]])),] # dim ncol(anc)XN + 1xN intercept
    
    # prediction of gene Exp
    pExp <- as.matrix(t(beta_mat) %*% t(genDat[[j]][fold,])) + as.matrix(t(beta_cov_mat) %*% t(cbind(covDat[fold,], rep(1,length(fold)))))  
    p <- getPrior(priorMat[[j]],pWeight)
    
    # regression
    obj[j] <- sum(rowMeans((pExp-t(expDat[fold, gene_ann$chrom == all_Chroms[j]]))^2))
    
    # penalization
    # exclude covariates matrix, no penalization applied to them 
    temp1 <- colSums((beta_mat^2)*p) 
    temp2 <-  colSums(abs(beta_mat)*p)
    pen[j] <- sum(lambdas[[j]]*(0.5*(1-alphas[[j]])*temp1 + alphas[[j]] * temp2))  
    
    
  }
  
  obj <- sum(obj)
  pen <- sum(pen)
  pen_par <- E_h*sum(pWeight^2)
  
  return(c(obj, pen, pen_par))
  
}

#### used in PriLer_part3_run.R ####
# find optimal lambda and alpha for each gene via CV, use all the samples, no prior weights
expPrediction_cv_noPrior_chr <- function(X, alpha, chr, genDat, seed, nfolds){
  
  ind_SNPs <- geneSnpDist[[chr]][,X]!=0 
  nSnp <- sum(ind_SNPs)
  
  if(nSnp>2){
    
    d <- which(ind_SNPs)[1]-1
    genotype <- as.matrix(cbind(genDat[,ind_SNPs],covDat))
    expressionValue <- as.numeric(expDat[,which(gene_ann$chrom == all_Chroms[chr])][,X])
    
    # compute deviance for the model using only covariates
    cov_mod <- cbind(expressionValue, covDat)
    fit_cov <- lm(cov_mod[,1]~as.matrix(cov_mod[,-1]))
    dev_lmcov <- summary(fit_cov)$r.squared
    
    # compute lambda ranges based on alpha 
    lambda_set <- vector(mode  = 'list', length = length(alpha))
    for(i in 1:length(alpha)){
      # print(i)
      lambda_set[[i]] <- glmnet(x = genotype, y = expressionValue,  alpha = alpha[i], nlambda = 100, intercept = T, standardize = F, 
                                penalty.factor = c(rep(1, length(which(ind_SNPs))), rep(0, ncol(covDat))))$lambda
      
    }
    
    cv_res <- vector(mode = 'list', length = length(alpha))  
    
    set.seed(seed)
    folds <- generateCVRuns(nfold = nfolds, labels = 1:nrow(genotype), ntimes = 1)[[1]]
    foldid_vect <- unlist(mapply(function(x, y) rep(x,length(y)), x =1:nfolds, y = folds, SIMPLIFY = F))
    foldid_vect <-  foldid_vect[order(unlist(folds))]
    
    for(i in 1:length(alpha)){
      
      cv_res[[i]] <- cv.glmnet(x = genotype, y = expressionValue, alpha = alpha[i], lambda = lambda_set[[i]], foldid = foldid_vect, 
                               penalty.factor = c(rep(1, length(which(ind_SNPs))), rep(0, ncol(covDat))), intercept = T, standardize = F) 
      
    }
    
    # keep only the best (lambda,alpha) pair in terms of cvm, 
    id_lambda <- sapply(cv_res, function(x) which(x$lambda == x$lambda.min))
    cvm_min <- sapply(1:length(alpha), function(x) cv_res[[x]]$cvm[id_lambda[x]])
    
    id_alpha <- which.min(cvm_min)
    alphamin <- alpha[id_alpha]
    
    # reassign keeping only the best lambda value
    res <- cv_res[[id_alpha]]
    
    ind <- which.min(res$cvm)
    selLambda <- res$lambda.min
    selAlpha <- alpha[[id_alpha]]
    betas <- res$glmnet.fit$beta[,ind]
    selIntercept <- res$glmnet.fit$a0[ind]
    r <- as.numeric(which(betas!=0))
    selBeta <- as.numeric(betas)[r]
    sInd <- r>nSnp
    v <- r[sInd]-nSnp
    r[r<=nSnp] <- r[r<=nSnp]+d
    r[sInd] <- covIndex_chr[,chr][v]
    r <- c(r, nrow(geneSnpDist[[chr]])+ncol(covDat)+1)
    r <- r[!is.na(r)]
    
    pred <- genotype %*% betas + selIntercept
    # SquErr <- mean((expressionValue-pred)^2)
    
    # compute dev explained by only genotype info
    P <- length(which(ind_SNPs))
    pred_geno <- genDat[, ind_SNPs] %*% betas[1:P]
    pred_cov <- selIntercept + as.matrix(covDat) %*% betas[(P+1):length(betas)]
    expressionValue_geno <- expressionValue - pred_cov
    # dev_geno_old <- 1-sum((expressionValue_geno-pred_geno)^2)/sum((expressionValue_geno - mean(expressionValue_geno))^2)
    dev_lmgeno <- res$glmnet.fit$dev.ratio[ind] - dev_lmcov
    
    dev_geno <- sum((pred_geno-mean(expressionValue_geno))^2 + 2*(expressionValue_geno - pred_geno)*(pred_geno - mean(expressionValue_geno)))/sum((expressionValue - mean(expressionValue))^2)
    dev_cov <- sum((pred_cov-mean(pred_cov))^2)/sum((expressionValue - mean(expressionValue))^2)
    dev_geno_cov <- sum(2*(expressionValue_geno - mean(expressionValue_geno))*(pred_cov - mean(pred_cov)))/sum((expressionValue - mean(expressionValue))^2)
    
    cor_est <- cor.test(expressionValue_geno,as.numeric(pred_geno))$estimate
    cor_pval <- cor.test(expressionValue_geno,as.numeric(pred_geno))$p.value
    
    # consider the original dataset
    cor_est_noadj <- cor.test(expressionValue,as.numeric(pred_geno))$estimate
    cor_pval_noadj <- cor.test(expressionValue,as.numeric(pred_geno))$p.value
    
    return(cbind(X,r,c(selBeta,as.numeric(selIntercept)),selLambda,selAlpha, res$glmnet.fit$dev.ratio[ind], dev_geno, dev_cov, dev_geno_cov, dev_lmgeno, cor_est, cor_pval, cor_est_noadj, cor_pval_noadj))
    
  }else{
    
    return(c(X,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA))
    
  }
}


#### used in PriLer_part3_run.R ####
# compute the regression using all the samples, use optimal alpha and lambda previously found
# all chr together
expPrediction_fin_chr <- function(X, prior, id_chr, genDat){
  
  gNam <- colnames(expDat[,gene_ann$chrom == all_Chroms[id_chr]])[X]
  ind_SNPs <- geneSnpDist[[id_chr]][,X]!=0 #&pDat[,2]!=0
  nSnp <- sum(ind_SNPs)
  
  # store results
  # out <- list()
  if(nSnp>2){
    
    d <- which(ind_SNPs)[1]-1
    
    # genotype=as.matrix(cbind(genDat[[id_chr]][,ind_SNPs], covDat)) #  not parallelized version
    genotype <- as.matrix(cbind(genDat[,ind_SNPs], covDat))
    expressionValue <- as.numeric(expDat[, gene_ann$chrom == all_Chroms[id_chr]][,X])
    
    # compute deviance for the model using only covariates
    cov_mod <- cbind(expressionValue, covDat)
    fit_cov <- lm(cov_mod[,1]~as.matrix(cov_mod[,-1]))
    dev_lmcov <- summary(fit_cov)$r.squared
    
    new_prior <- c(prior[ind_SNPs], rep(1,ncol(covDat)))
    new_genotype <- sweep(genotype, 2, new_prior, '/')
    
    res <- glmnet(new_genotype, expressionValue, alpha=alphaVec[[id_chr]][X], lambda=lambdaVec[[id_chr]][X], standardize=FALSE, intercept = TRUE, 
                  penalty.factor = c(rep(1,nSnp), rep(0, ncol(covDat)))) 
    beta <- res$beta/new_prior
    selLambda <- lambdaVec[[id_chr]][X]
    selAlpha <- alphaVec[[id_chr]][X]
    selIntercept <- res$a0
    r <- as.numeric(which(beta!=0))
    selBeta <- as.numeric(beta)[r]
    sInd <- r>nSnp
    v <- r[sInd]-nSnp
    r[r<=nSnp] <- r[r<=nSnp]+d
    r[sInd] <- covIndex_chr[,id_chr][v]
    r <- c(r, nrow(geneSnpDist[[id_chr]])+ncol(covDat)+1)
    r <- r[!is.na(r)]
    
    pred <- genotype %*% beta + selIntercept
    SquErr <- mean((expressionValue-pred)^2)
    # compute dev explained by only genotype info
    P <- length(which(ind_SNPs))
    # pred_geno <- genDat[[id_chr]][, ind_SNPs] %*% beta[1:P] # not parallelized version
    pred_geno <- genDat[, ind_SNPs] %*% beta[1:P]
    pred_cov <- selIntercept + as.matrix(covDat) %*% beta[(P+1):length(beta)]
    expressionValue_geno <- expressionValue - pred_cov
    # dev_geno_old <- 1-sum((expressionValue_geno-pred_geno)^2)/sum((expressionValue_geno - mean(expressionValue_geno))^2)
    dev_lmgeno <- res$dev.ratio - dev_lmcov
    
    dev_geno <- sum((pred_geno-mean(expressionValue_geno))^2 + 2*(expressionValue_geno - pred_geno)*(pred_geno - mean(expressionValue_geno)))/sum((expressionValue - mean(expressionValue))^2)
    dev_cov <- sum((pred_cov-mean(pred_cov))^2)/sum((expressionValue - mean(expressionValue))^2)
    dev_geno_cov <- sum(2*(expressionValue_geno - mean(expressionValue_geno))*(pred_cov - mean(pred_cov)))/sum((expressionValue - mean(expressionValue))^2)
    
    cor_est <- cor.test(expressionValue_geno,as.numeric(pred_geno))$estimate
    cor_pval <- cor.test(expressionValue_geno,as.numeric(pred_geno))$p.value
    # consider the original dataset
    cor_est_noadj <- cor.test(expressionValue,as.numeric(pred_geno))$estimate
    cor_pval_noadj <- cor.test(expressionValue,as.numeric(pred_geno))$p.value
    
    
    return(cbind(X,r,c(selBeta,as.numeric(selIntercept)), selLambda, selAlpha, SquErr, res$dev.ratio, dev_geno, dev_cov, dev_geno_cov, dev_lmgeno, cor_est, cor_pval, cor_est_noadj,cor_pval_noadj))
    
  }else{
    
    return(c(X,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA))
    
  }
  
}

###################################
######## REVIEW FROM HERE #########
###################################

#### used in ElNet_withPrior_part4_run.R ####
## predict on new genes
# compute regression for fixed alpha and lambda parameter for each chromosome
expPrediction_Prior_fold <- function(X, prior, lambda, alpha, fold){
  
  ind_SNPs=geneSnpDist[,X]!=0 #&pDat[,2]!=0
  nSnp=sum(ind_SNPs)
  
  if (nSnp>2 & !is.na(lambda)){
    
    # # ######
    # cons <- setdiff(which(ind_SNPs)[1]:which(ind_SNPs)[length(which(ind_SNPs))], which(ind_SNPs))
    # if(length(cons)>0){
    #   print(paste0('index problem:', cons))
    # }
    # # use to check the SNPs id are consecutive
    # # ######
    
    d=which(ind_SNPs)[1]-1
    genotype=as.matrix(cbind(genDat[fold,ind_SNPs],covDat[fold,]))
    expressionValue=as.numeric(expDat[fold,X])
    
    # compute deviance for the model using only covariates
    cov_mod <- cbind(expressionValue, covDat[fold,])
    fit_cov <- lm(cov_mod[,1]~as.matrix(cov_mod[,-1]))
    dev_lmcov <- summary(fit_cov)$r.squared
    
    new_prior <- c(prior[ind_SNPs], rep(1,ncol(covDat[fold,])))
    # insert weights on the X matrix, so that they are not automatically normalized by glmnet
    new_genotype <- sweep(genotype, 2, new_prior, '/')
    
    # lambdaMax <- max(abs(expressionValue %*% new_genotype)/(nrow(new_genotype)*alpha[X]))
    
    res <- glmnet(new_genotype, expressionValue, alpha=alpha[X], lambda=lambda[X], standardize=FALSE, intercept = TRUE, 
                  penalty.factor = c(rep(1, nSnp), rep(0, ncol(covDat[fold,])))) 
    
    
    # pred1 <- predict.glmnet(object = res, newx = genotype)
    beta <- res$beta/new_prior
    
    selLambda <- lambda[X]
    selAlpha <- alpha[X]
    # beta <- res$beta
    selIntercept <- res$a0
    r <- as.numeric(which(beta!=0))
    selBeta <- as.numeric(beta)[r]
    sInd=r>nSnp
    v=r[sInd]-nSnp
    r[r<=nSnp]=r[r<=nSnp]+d
    r[sInd]=covIndex[v]
    r <- c(r, nrow(geneSnpDist)+ncol(covDat[fold,])+1)
    r <- r[!is.na(r)]
    
    pred <- genotype %*% beta + selIntercept
    SquErr <- mean((expressionValue-pred)^2)
    
    P <- nSnp
    pred_geno <- genDat[fold,ind_SNPs] %*% beta[1:P]
    pred_cov <- as.matrix(covDat[fold,]) %*% beta[(P+1):length(beta)] + selIntercept
    expressionValue_geno <- expressionValue - pred_cov
    # dev_geno_old <- 1-sum((expressionValue_geno-pred_geno)^2)/sum((expressionValue_geno - mean(expressionValue_geno))^2)
    
    dev_geno <- sum((pred_geno-mean(expressionValue_geno))^2 + 2*(expressionValue_geno - pred_geno)*(pred_geno - mean(expressionValue_geno)))/sum((expressionValue - mean(expressionValue))^2)
    dev_cov <- sum((pred_cov-mean(pred_cov))^2)/sum((expressionValue - mean(expressionValue))^2)
    dev_geno_cov <- sum(2*(expressionValue_geno - mean(expressionValue_geno))*(pred_cov - mean(pred_cov)))/sum((expressionValue - mean(expressionValue))^2)
    
    dev_lmgeno <- res$dev.ratio - dev_lmcov
    
    cor_est <- cor.test(expressionValue_geno,as.numeric(pred_geno))$estimate
    cor_pval <- cor.test(expressionValue_geno,as.numeric(pred_geno))$p.value
    
    
    # return(cbind(X,r,c(selBeta,as.numeric(selIntercept)),selLambda,selAlpha,res$dev.ratio))
    return(cbind(X,r,c(selBeta,as.numeric(selIntercept)), selLambda, selAlpha, SquErr, res$dev.ratio, dev_geno, dev_cov, dev_geno_cov, dev_lmgeno,  cor_est, cor_pval))
    
  }else{
    
    
    return(c(X,NA,NA,NA,NA,NA,NA,NA,NA,NA, NA, NA, NA))
  }
  
}

