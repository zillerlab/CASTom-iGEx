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

# combine function for parallelization
comb <- function(x, ...) {
  lapply(seq_along(x),
         function(i) c(x[[i]], lapply(list(...), function(y) y[[i]])))
}

############
# normalize data
normalize_feat <- function(data, par = NULL){
  
  if(is.null(par)){
    
    new <- scale(data)
    par <- data.frame(id = colnames(data), sd = attr(new, "scaled:scale"), mean = attr(new, "scaled:center"))
    attr(new, "scaled:scale") <- NULL
    attr(new, "scaled:center") <- NULL
    
  }else{
    
    if(!identical(colnames(data), par$id)){stop('features train-test must be in the same order')}
    new <- scale(data, center = par$mean, scale = par$sd)
    attr(new, "scaled:scale") <- NULL
    attr(new, "scaled:center") <- NULL
    
  }
  
  return(list(data = new, par = par))
  
}

# SNF_new <- function (Wall, K = 20, t = 20){
#   check_wall_names <- function(Wall) {
#     name_match <- function(names_A, names_B) {
#       return(identical(dimnames(names_A), dimnames(names_B)))
#     }
#     return(all(unlist(lapply(Wall, FUN = name_match, Wall[[1]]))))
#   }
#   wall.name.check <- check_wall_names(Wall)
#   wall.names <- dimnames(Wall[[1]])
#   if (!wall.name.check) {
#     warning("Dim names not consistent across all matrices in Wall.\n            Returned matrix will have no dim names.")
#   }
#   LW <- length(Wall)
#   normalize <- function(X){
#   row.sum.mdiag <- (rowSums(X) - diag(X))*nrow(X)
#   row.sum.mdiag[row.sum.mdiag == 0] <- 1
#   X <- X/(2 * (row.sum.mdiag))
#   diag(X) <- 0.5
#   return(X)
#   }
#   newW <- vector("list", LW)
#   nextW <- vector("list", LW)
#   for (i in 1:LW) {
#     Wall[[i]] <- normalize(Wall[[i]])
#     Wall[[i]] <- (Wall[[i]] + t(Wall[[i]]))/2
#   }
#   for (i in 1:LW) {
#     newW[[i]] <- (.dominateset(Wall[[i]], K))
#   }
#   for (i in 1:t) {
#     for (j in 1:LW) {
#       sumWJ <- matrix(0, dim(Wall[[j]])[1], dim(Wall[[j]])[2])
#       for (k in 1:LW) {
#         if (k != j) {
#           sumWJ <- sumWJ + Wall[[k]]
#         }
#       }
#       nextW[[j]] <- newW[[j]] %*% (sumWJ/(LW - 1)) %*% 
#         t(newW[[j]])
#     }
#     for (j in 1:LW) {
#       Wall[[j]] <- normalize(nextW[[j]])
#       Wall[[j]] <- (Wall[[j]] + t(Wall[[j]]))/2
#     }
#   }
#   W <- matrix(0, nrow(Wall[[1]]), ncol(Wall[[1]]))
#   for (i in 1:LW) {
#     W <- W + Wall[[i]]
#   }
#   W <- W/LW
#   W <- normalize(W)
#   W <- (W + t(W))/2
#   if (wall.name.check) {
#     dimnames(W) <- wall.names
#   }
#   return(W)
# }


# SNF_matrix <- function(data, kNN = 10, alpha = 0.5, tpar = 20, test_data = NULL){
#   
#   if(length(data)<2){stop('number of matrices must be more than 1')}
#   if(is.null(test_data)){
#     tot_data <- data
#   }else{
#     tot_data <- mapply(function(x,y) rbind(x,y), x = data, y = test_data, SIMPLIFY = F) 
#   }
#   
#   dist_mat <-  lapply(tot_data, function(x) (dist2(as.matrix(x),as.matrix(x)))^(1/2))
#   W_single <- lapply(dist_mat, function(x) affinityMatrix(x, kNN, alpha))
#   # overall similarity matrix
#   W <- SNF(W_single, kNN, tpar)
#   
#   if(!is.null(test_data)){
#     W <- W[match(rownames(test_data[[1]]), rownames(W)), match(rownames(data[[1]]), rownames(W))]
#   }
#   return(W)
#   
# }

SNF_matrix <- function(data, kNN = 10, alpha = 0.5, tpar = 20){

  dist_mat <-  lapply(tot_data, function(x) (dist2(as.matrix(x),as.matrix(x)))^(1/2))
  W_single <- lapply(dist_mat, function(x) affinityMatrix(x, kNN, alpha))
  # overall similarity matrix
  W <- SNF(W_single, kNN, tpar)

  return(W)

}

UMAP_model <- function(mat, model = NULL, n_comp = 2, n_neigh = 30, min_dist = 0.01, seed_umap = 67){

  if(is.null(model)){

    set.seed(seed_umap)
    model_umap <- uwot::umap(mat, n_neighbors = n_neigh, n_components = n_comp, min_dist = min_dist, ret_model = T)
    # model_umap <- umap::umap(mat, custom.setting)
    # df <- as.data.frame(model_umap$layout)
    df <- as.data.frame(model_umap$embedding)
    colnames(df) <- paste0('component_', 1:n_comp)
    rownames(df) <- rownames(mat)

  }else{

    set.seed(seed_umap)
    umap_res <- umap_transform(mat, model)
    #umap_res <- predict(model,mat)
    df <- as.data.frame(umap_res)
    colnames(df) <- paste0('component_', 1:n_comp)
    rownames(df) <- rownames(mat)
    model_umap <- NULL
  }
  return(list(data = df, model = model_umap))

}

UMAP_sim <- function(mat, n_comp = 2, n_neigh = 30, min_dist = 0.01, seed_umap = 67){

  set.seed(seed_umap)
  
  dist_mat <- 0.5 - mat # specific for SNF, diagonal all 0.5
  dist_mat <- as.dist(dist_mat)
  
  model_umap <- uwot::umap(dist_mat, n_neighbors = n_neigh, n_components = n_comp, min_dist = min_dist)
  df <- as.data.frame(model_umap)
  colnames(df) <- paste0('component_', 1:n_comp)
  rownames(df) <- rownames(mat)
    
  return(df)

}


# #### error for SVM classification #####
balanced_acc_error <- function(true, pred){
  
  val_t <- sort(unique(as.numeric(as.character(true))))
  val_p <- sort(unique(as.numeric(as.character(pred))))
  min_t <- val_t[1]
  min_p <- val_p[1]
  max_t <- val_t[2]
  max_p <- val_p[2]
  
  true <- as.numeric(as.character(true))
  pred <- as.numeric(as.character(pred))
  
  TN <- length(which(true == min_t & pred == min_p))
  N <- length(which(true == min_t))
  
  TP <- length(which(true == max_p & pred == max_p))
  P <- length(which(true == max_t))
  
  sens <- TP/P
  spec <- TN/N
  
  # tab <- table(true, pred)
  # sens <- tab[2,2]/sum(tab[2,])
  # spec <- tab[1,1]/sum(tab[1,])
  return((sens+spec)/2)
  
}

#### error for SVM classification #####
MCC_err <- function(true, pred){
  
  val_t <- sort(unique(as.numeric(as.character(true))))
  val_p <- sort(unique(as.numeric(as.character(pred))))
  min_t <- val_t[1]
  min_p <- val_p[1]
  max_t <- val_t[2]
  max_p <- val_p[2]
  
  true <- as.numeric(as.character(true))
  pred <- as.numeric(as.character(pred))
  
  TN <- length(which(true == min_t & pred == min_p))
  FN <- length(which(true == max_t & pred == min_p))
  N <- length(which(true == min_t))
  
  TP <- length(which(true == max_t & pred == max_p))
  FP <- length(which(true == min_t & pred == max_p))
  P <- length(which(true == max_t))
  
  res <- (TP*TN-FP*FN)/(sqrt(TP+FP)*sqrt(TP+FN)*sqrt(TN+FP)*sqrt(TN+FN))
  
  return(res)
  
}




SVM_model <- function(data, info_sample, n_cv = 5, ncores = 5){
  
  if(!identical(rownames(data), info_sample$Individual_ID)){stop('sample name does not match')}
  data$y <- factor(info_sample$Dx)
  
  ## split data ##
  set.seed(20)
  fold_cv_out <- generateCVRuns(labels = info_sample$Dx, ntimes = 1, nfold = n_cv, stratified = T)
  val_id <- train_id <- test_cv_id <- output_nCV <- vector(mode = 'list', length = n_cv)
  gamma_set <- C_set <- 10^(-3:3)
  
  # nested 5-fold CV split
  for(n in 1:n_cv){
    
    print(paste0('outer fold', n))
    
    val_id[[n]] <- info_sample$Individual_ID[sort(fold_cv_out[[1]][[n]])]
    train_id[[n]] <- info_sample$Individual_ID[-sort(fold_cv_out[[1]][[n]])]
    
    set.seed(782)
    fold_cv <- generateCVRuns(labels = info_sample$Dx[match(train_id[[n]], info_sample$Individual_ID)], ntimes = 1, nfold = n_cv, stratified = T)
    fold_cv <- fold_cv[[1]]
    test_cv_id[[n]] <- lapply(fold_cv, function(x) train_id[[n]][x])
    
    ### 5-fold CV ###
    registerDoParallel(min(c(ncores, n_cv)))
    
    res <- foreach(id_cv=1:n_cv)%dopar%{
      #for(id_cv in 1:n_cv){
      
      print(id_cv)
      
      train_cv_id <- setdiff(train_id[[n]], test_cv_id[[n]][[id_cv]])
      tmp_train <- data[match(train_cv_id,rownames(data)),]
      tmp_test <- data[match(test_cv_id[[n]][[id_cv]], rownames(data)),]
      mat_error <- data.frame(gamma = unlist(lapply(gamma_set, function(x) rep(x, length(C_set)))),
                              C = rep(C_set, length(gamma_set)), balanced_acc = NA, mcc = NA)
      for(g in 1:length(gamma_set)){
        for(j in 1:length(C_set)){
          svm_model <- svm(y~., data = tmp_train, kernel="radial", cost=C_set[j], gamma=gamma_set[g], scale = F)
          pred <- predict(svm_model, tmp_test)
          mat_error[mat_error$gamma == gamma_set[g] & mat_error$C == C_set[j], 3:4] <- 
            c(balanced_acc_error(tmp_test$y, pred), MCC_err(tmp_test$y, pred))
        }
      }
      
      mat_error$fold <- id_cv
      mat_error
    }
    
    # best parameter combination
    balanced_acc <- sapply(res, function(x) x$balanced_acc)
    m_cv_bacc <- max(rowMeans(balanced_acc))
    sd_cv_bacc <- sd(balanced_acc[which.max(rowMeans(balanced_acc)),])
    id_max <- which.max(rowMeans(balanced_acc))
    best_par <- data.frame(gamma = res[[1]]$gamma[id_max], C = res[[1]]$C[id_max])
    
    # create final model
    tmp_train <- data[match(train_id[[n]],rownames(data)),]
    tmp_test <- data[match(val_id[[n]], rownames(data)),]
    
    # svm model
    svm_model <- svm(y~., data = tmp_train, kernel="radial", cost=best_par$C, gamma= best_par$gamma, scale = F, probability=T)
    pred <- predict(svm_model, tmp_test, probability=T)
    prob_val <- attr(pred, "probabilities")[,colnames(attr(pred, "probabilities")) == 1]
    attr(pred, "probabilities") <- NULL
    output_nCV[[n]] <- list(sample_config = list(train = train_id[[n]], val = val_id[[n]], test_cv_id = test_cv_id[[n]]), cv_performance = res, best_par_cv = best_par, 
                            train_svm = svm_model, pred_val = data.frame(y = tmp_test$y, pred = pred, prob = prob_val))
    output_nCV[[n]]$perf <- cbind(best_par, data.frame(bacc_cv_mean = m_cv_bacc, bacc_cv_sd = sd_cv_bacc, bacc_train = balanced_acc_error(tmp_train$y,svm_model$fitted), 
                                                       bacc_val = balanced_acc_error(tmp_test$y, pred), auc_val = as.numeric(roc(tmp_test$y,prob_val)$auc)))
  }
  
  tot_perf <- as.data.frame(do.call(rbind,lapply(output_nCV, function(x) x$perf)))
  tot_perf$cv <- 1:n_cv

  output <- list(nested_cv = output_nCV, performance_nestedCV = tot_perf)
  
  # create final model (5-fold CV to find best parameters)
  # 5-fold CV split
  test_cv_id <- lapply(fold_cv_out[[1]], function(x) info_sample$Individual_ID[x])
  
  registerDoParallel(min(c(ncores, n_cv)))
  res <- foreach(id_cv=1:n_cv)%dopar%{
    #for(id_cv in 1:n_cv){
    
    print(id_cv)
    
    train_cv_id <- setdiff(info_sample$Individual_ID, test_cv_id[[id_cv]])
    tmp_train <- data[match(train_cv_id,rownames(data)),]
    tmp_test <- data[match(test_cv_id[[id_cv]], rownames(data)),]
    mat_error <- data.frame(gamma = unlist(lapply(gamma_set, function(x) rep(x, length(C_set)))),
                            C = rep(C_set, length(gamma_set)), balanced_acc = NA, mcc = NA)
    for(g in 1:length(gamma_set)){
      for(j in 1:length(C_set)){
        svm_model <- svm(y~., data = tmp_train, kernel="radial", cost=C_set[j], gamma=gamma_set[g], scale = F)
        pred <- predict(svm_model, tmp_test)
        mat_error[mat_error$gamma == gamma_set[g] & mat_error$C == C_set[j], 3:4] <- 
          c(balanced_acc_error(tmp_test$y, pred), MCC_err(tmp_test$y, pred))
      }
    }
    
    mat_error$fold <- id_cv
    mat_error
  }
  
  # best parameter combination
  balanced_acc <- sapply(res, function(x) x$balanced_acc)
  m_cv_bacc <- max(rowMeans(balanced_acc))
  sd_cv_bacc <- sd(balanced_acc[which.max(rowMeans(balanced_acc)),])
  id_max <- which.max(rowMeans(balanced_acc))
  best_par <- data.frame(gamma = res[[1]]$gamma[id_max], C = res[[1]]$C[id_max])
  
  # svm model
  svm_model <- svm(y~., data = data, kernel="radial", cost=best_par$C, gamma= best_par$gamma, scale = F, probability=T)
  output$svm_model <- svm_model
  output$performance_CV <- cbind(best_par, data.frame(bacc_cv_mean = m_cv_bacc, bacc_cv_sd = sd_cv_bacc, bacc_final = balanced_acc_error(data$y, svm_model$fitted)))
  
  return(output)

}

SVM_model_test <- function(model, data, info_sample){
  
  if(!identical(rownames(data), info_sample$Individual_ID)){stop('sample name does not match')}
  pred <- predict(model, data,  probability=T)
  prob <- attr(pred,"probabilities")
  prob <- prob[,colnames(prob) == 1]
  
  perf <- data.frame(bacc_test = balanced_acc_error(info_sample$Dx, pred),  auc_test = as.numeric(roc(info_sample$Dx,prob)$auc))
  
  return(list(performance_test = perf, prediction = data.frame(pred = pred, prob = prob))) 
  
}





