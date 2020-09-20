## cluster using PG method ##
clust_PGmethod_HKsim <- function(kNN, score, type_Dx, multiple_cohorts = F, sample_info, euclDist){
  
  print(kNN)
  # 2) find kNN based on eucl_dist
  sigma_kNN <- apply(euclDist, 1, function(x) mean(x[order(x)[(2:kNN+1)]]))
  # 3) compute for each pair the customized sigma and HK similarity
  tmp <- sapply(1:(nrow(score)-1), function(x) c(rep(0,x), (sigma_kNN[x] + sigma_kNN[(x+1):nrow(score)] + euclDist[x, (x+1):nrow(score)])/3))
  tmp <- cbind(tmp, rep(0, nrow(score)))
  eps <- tmp+t(tmp)
  rm(tmp)
  W <- exp(-euclDist/(0.5*eps))
  rm(eps)
  print(mem_used())
  print('matrix W built')
  
  # 4') find shared neigbours based on W
  list_kNN_W <- apply(W, 1, function(x) order(x, decreasing = T)[1:kNN])
  rm(W)
  # 5') build weigths based on shared neighbours
  tmp <- vector(mode = 'list', length = ncol(list_kNN_W))
  id_tmp <- vector(mode = 'list', length = ncol(list_kNN_W))
  for(i in 1:(nrow(score)-1)){
    # print(i)
    tmp[[i]] <- sapply((i+1):ncol(list_kNN_W), function(x) length(intersect(list_kNN_W[,i],list_kNN_W[,x])))
    tmp[[i]] <- tmp[[i]]/(kNN*2 - tmp[[i]])
    id_tmp[[i]] <- which(tmp[[i]]!=0)+i
    tmp[[i]] <- tmp[[i]][tmp[[i]]!=0]
  }
  W_sNN <- sparseMatrix(i = unlist(lapply(1:nrow(score), function(x) rep(x, length(id_tmp[[x]])))), 
                        j = unlist(id_tmp), x = unlist(tmp), symmetric = T)
  rm(tmp)
  rm(id_tmp)
  print(mem_used())
  print('matrix W_sNN built')
  
  # 6') Louvain method
  graph_W_sNN <- graph_from_adjacency_matrix(W_sNN, weighted=TRUE, mode = 'undirected')
  louv_W_sNN <- cluster_louvain(graph_W_sNN)
  
  tmp_info <- data.frame(kNN = kNN, mod = max(louv_W_sNN$modularity), n_gr = length(unique(louv_W_sNN$membership)))
  cl <- louv_W_sNN$membership
  
  clos_gr <- lapply(sort(unique(cl)), function(x) closeness(normalized = T, induced.subgraph(graph_W_sNN, vids = which(cl == x)), mode = 'all'))
  central_node_gr <- mapply(function(x,y)  which(cl == x)[y], x = sort(unique(cl)), y = sapply(clos_gr, which.max))
  clos_central_gr <- sapply(clos_gr, max)
  print(clos_central_gr)
  invshort_path_centralnode <- 1/distances(graph_W_sNN, v = central_node_gr, to = central_node_gr, mode = 'all')
  diag(invshort_path_centralnode) <- 0
  DBindex <- sapply(1:length(unique(cl)), function(x) min((clos_central_gr[x] + clos_central_gr[-x])/invshort_path_centralnode[x, -x]))
  print(DBindex)
  
  tmp_info$DB_mean <- mean(DBindex)
  tmp_info$DB_sd <- sd(DBindex)
  if(multiple_cohorts){
    tmp_info$cohort_nmi <- compare(louv_W_sNN$membership, sample_info$cohort_id, method = 'nmi')  
  }
  
  if(type_Dx == 'All'){tmp_info$Dx_nmi <-  compare(louv_W_sNN$membership, sample_info$Dx, method = 'nmi')}
  return(list(cl = louv_W_sNN, info = tmp_info, sim_mat = W_sNN, sigma_kNN = sigma_kNN, kNN_W = t(list_kNN_W)))
  
}


project_clust_PGmethod_HKsim <- function(kNN, score, data_mod, sample_info, euclDist, cl_mod){
  
  # build graph based on HK
  # 2) find kNN based on eucl_dist
  sigma_kNN <- apply(euclDist, 1, function(x) mean(x[order(x)[(2:kNN+1)]]))
  # 3) compute for each pair the customized sigma and HK similarity
  tmp <- sapply(1:(nrow(score)-1), function(x) c(rep(0,x), (sigma_kNN[x] + sigma_kNN[(x+1):nrow(data_tot)] + ed_dist[x, (x+1):nrow(data_tot)])/3))
  tmp <- cbind(tmp, rep(0, nrow(score)))
  eps <- tmp+t(tmp)
  rm(tmp)
  W <- exp(-euclDist/(0.5*eps))
  rm(eps)
  print(mem_used())
  print('matrix W built')
  # 4) find shared neigbours based on W
  list_kNN_W <- apply(W, 1, function(x) order(x, decreasing = T)[1:kNN])
  rm(W)
  # 5) build weigths based on shared neighbours
  tmp <- vector(mode = 'list', length = ncol(list_kNN_W))
  id_tmp <- vector(mode = 'list', length = ncol(list_kNN_W))
  for(i in 1:(nrow(score)-1)){
    # print(i)
    tmp[[i]] <- sapply((i+1):ncol(list_kNN_W), function(x) length(intersect(list_kNN_W[,i],list_kNN_W[,x])))
    tmp[[i]] <- tmp[[i]]/(kNN*2 - tmp[[i]])
    id_tmp[[i]] <- which(tmp[[i]]!=0)+i
    tmp[[i]] <- tmp[[i]][tmp[[i]]!=0]
  }
  W_sNN <- sparseMatrix(i = unlist(lapply(1:nrow(score), function(x) rep(x, length(id_tmp[[x]])))), 
                        j = unlist(id_tmp), x = unlist(tmp), symmetric = T)
  rm(tmp)
  rm(id_tmp)
  print(mem_used())
  print('matrix W_sNN built')
  # 6) build final graph
  # graph_W_sNN <- graph_from_adjacency_matrix(W_sNN, weighted=TRUE, mode = 'undirected')
  # 7) build laplacian of the graph
  L <- Matrix(diag(rowSums(W_sNN)), sparse = T) - W_sNN 
  # 8) decompoase laplacian in blocks + create label Matrix
  L_u <- L[(nrow(data_mod)+1):nrow(L), (nrow(data_mod)+1):nrow(L)]
  Bt <- L[(nrow(data_mod)+1):nrow(L), 1:nrow(data_mod)]
  Q <- model.matrix(~0+factor(cl_mod$gr, levels = sort(unique(cl_mod$gr))))
  Q <- Matrix(Q)
  RHS <- -Bt %*% Q
  # 9) solve linear equation system
  P <- SparseM::solve(a = L_u, b = RHS)
  prob_mat <- as.data.frame(as.matrix(P))
  colnames(prob_mat) <- paste0('gr_', sort(unique(cl_mod$gr)))
  prob_mat <- cbind(data.frame(Individiual_ID = sample_info$Individual_ID, Dx = sample_info$Dx),prob_mat)
  
  return(list(probability = prob_mat, tot_W_sNN = W_sNN))
}


clust_PGmethod_EDdist <- function(kNN, score, type_Dx, multiple_cohorts = F, sample_info, euclDist){
  
  print(kNN)
  
  # 4') find shared neigbours based on W
  list_kNN_W <- apply(euclDist, 1, function(x) order(x)[1:kNN])
  # rm(W)
  # 5') build weigths based on shared neighbours
  tmp <- vector(mode = 'list', length = ncol(list_kNN_W))
  id_tmp <- vector(mode = 'list', length = ncol(list_kNN_W))
  for(i in 1:(nrow(score)-1)){
    # print(i)
    tmp[[i]] <- sapply((i+1):ncol(list_kNN_W), function(x) length(intersect(list_kNN_W[,i],list_kNN_W[,x])))
    tmp[[i]] <- tmp[[i]]/(kNN*2 - tmp[[i]])
    id_tmp[[i]] <- which(tmp[[i]]!=0)+i
    tmp[[i]] <- tmp[[i]][tmp[[i]]!=0]
  }
  W_sNN <- sparseMatrix(i = unlist(lapply(1:nrow(score), function(x) rep(x, length(id_tmp[[x]])))), 
                        j = unlist(id_tmp), x = unlist(tmp), symmetric = T)
  rm(tmp)
  rm(id_tmp)
  print(mem_used())
  print('matrix W_sNN built')
  
  # 6') Louvain method
  graph_W_sNN <- graph_from_adjacency_matrix(W_sNN, weighted=TRUE, mode = 'undirected')
  louv_W_sNN <- cluster_louvain(graph_W_sNN)
  
  tmp_info <- data.frame(kNN = kNN, mod = max(louv_W_sNN$modularity), n_gr = length(unique(louv_W_sNN$membership)))
  cl <- louv_W_sNN$membership
  
  clos_gr <- lapply(sort(unique(cl)), function(x) closeness(normalized = T, induced.subgraph(graph_W_sNN, vids = which(cl == x)),  mode = 'all'))
  central_node_gr <- mapply(function(x,y)  which(cl == x)[y], x = sort(unique(cl)), y = sapply(clos_gr, which.max))
  clos_central_gr <- sapply(clos_gr, max)
  invshort_path_centralnode <- 1/distances(graph_W_sNN, v = central_node_gr, to = central_node_gr, mode = 'all')
  diag(invshort_path_centralnode) <- 0
  DBindex <- sapply(1:length(unique(cl)), function(x) min((clos_central_gr[x] + clos_central_gr[-x])/invshort_path_centralnode[x, -x]))
  print(DBindex)
  
  tmp_info$DB_mean <- mean(DBindex)
  tmp_info$DB_sd <- sd(DBindex)
  if(multiple_cohorts){
    tmp_info$cohort_nmi <- compare(louv_W_sNN$membership, sample_info$cohort_id, method = 'nmi')  
  }
  
  if(type_Dx == 'All'){tmp_info$Dx_nmi <-  compare(louv_W_sNN$membership, sample_info$Dx, method = 'nmi')}
  sigma_kNN = NA
  return(list(cl = louv_W_sNN, info = tmp_info, sim_mat = W_sNN, kNN_W = t(list_kNN_W)))
  
}

# clustering using a precomputed similairty matrix
clust_PGmethod_sim <- function(kNN, sim, type_Dx, multiple_cohorts = F, sample_info){
  print(kNN)
  # 4') find shared neigbours based on W
  list_kNN_W <- apply(sim, 1, function(x) order(x, decreasing = T)[1:kNN])
  # rm(W)
  # 5') build weigths based on shared neighbours
  tmp <- vector(mode = 'list', length = ncol(list_kNN_W))
  id_tmp <- vector(mode = 'list', length = ncol(list_kNN_W))
  for(i in 1:(nrow(sim)-1)){
    # print(i)
    tmp[[i]] <- sapply((i+1):ncol(list_kNN_W), function(x) length(intersect(list_kNN_W[,i],list_kNN_W[,x])))
    tmp[[i]] <- tmp[[i]]/(kNN*2 - tmp[[i]])
    id_tmp[[i]] <- which(tmp[[i]]!=0)+i
    tmp[[i]] <- tmp[[i]][tmp[[i]]!=0]
  }
  W_sNN <- sparseMatrix(i = unlist(lapply(1:nrow(sim), function(x) rep(x, length(id_tmp[[x]])))), 
                        j = unlist(id_tmp), x = unlist(tmp), symmetric = T)
  rm(tmp)
  rm(id_tmp)
  print(mem_used())
  print('matrix W_sNN built')
  
  # 6') Louvain method
  graph_W_sNN <- graph_from_adjacency_matrix(W_sNN, weighted=TRUE, mode = 'undirected')
  louv_W_sNN <- cluster_louvain(graph_W_sNN)
  
  tmp_info <- data.frame(kNN = kNN, mod = max(louv_W_sNN$modularity), n_gr = length(unique(louv_W_sNN$membership)))
  cl <- louv_W_sNN$membership
  
  clos_gr <- lapply(sort(unique(cl)), function(x) closeness(normalized = T, induced.subgraph(graph_W_sNN, vids = which(cl == x)),  mode = 'all'))
  central_node_gr <- mapply(function(x,y)  which(cl == x)[y], x = sort(unique(cl)), y = sapply(clos_gr, which.max))
  clos_central_gr <- sapply(clos_gr, max)
  invshort_path_centralnode <- 1/distances(graph_W_sNN, v = central_node_gr, to = central_node_gr, mode = 'all')
  diag(invshort_path_centralnode) <- 0
  DBindex <- sapply(1:length(unique(cl)), function(x) min((clos_central_gr[x] + clos_central_gr[-x])/invshort_path_centralnode[x, -x]))
  print(DBindex)
  
  tmp_info$DB_mean <- mean(DBindex)
  tmp_info$DB_sd <- sd(DBindex)
  if(multiple_cohorts){
    tmp_info$cohort_nmi <- compare(louv_W_sNN$membership, sample_info$cohort_id, method = 'nmi')  
  }
  
  if(type_Dx == 'All'){tmp_info$Dx_nmi <-  compare(louv_W_sNN$membership, sample_info$Dx, method = 'nmi')}
  sigma_kNN = NA
  return(list(cl = louv_W_sNN, info = tmp_info, sim_mat = W_sNN, kNN_W = t(list_kNN_W)))
  
}

### plot
### heatmap ###
save_pheatmap_png <- function(x, filename, width=10, height=10, res = 150) {
  png(filename, width = width, height = height, res = res, units = 'in')
  grid::grid.newpage()
  grid::grid.draw(x$gtable)
  dev.off()
}
save_pheatmap_pdf <- function(x, filename, width=10, height=10, res = 150) {
  pdf(filename, width = width, height = height)
  grid::grid.newpage()
  grid::grid.draw(x$gtable)
  dev.off()
}

pheat_pl <- function(mat, cl, type_mat, height_pl = 10, width_pl = 7, outFile ){
  
  coul <- rev(colorRampPalette(brewer.pal(11, "RdBu"))(100))
  # match mat and cl
  cl <- cl[match(rownames(mat), cl$id),]
  
  tmp_mat <- as.matrix(mat)
  val <- max(abs(tmp_mat))
  mat_breaks <- seq(-val, val, length.out = 100)
  
  tmp_mat[tmp_mat>=val] <- val
  tmp_mat[tmp_mat<=-val] <- -val
  
  id <- order(cl$gr)
  P <- length(unique(cl$gr))
  tmp_mat <- tmp_mat[id,]
  
  mat_row <- data.frame(group = paste0('cl', cl$gr[id]))
  rownames(mat_row) <- rownames(tmp_mat)
  
  mat_colors_gr <- list(group = rep(c('#444444', '#C1C1C1'),P)[1:P])
  names(mat_colors_gr$group) <- paste0('cl', 1:P)
  
  
  hm_pl <- pheatmap(mat=tmp_mat, color=coul, show_colnames = F, show_rownames = F, annotation_names_col = F, annotation_names_row = F, 
                    cluster_rows=F, cluster_cols=T, border_color = NA, annotation_colors = mat_colors_gr,
                    annotation_row = mat_row, drop_levels = TRUE, fontsize_row = 6, fontsize_col = 8, fontsize = 12, 
                    main =  sprintf("%s", type_mat),
                    breaks = mat_breaks, treeheight_row = 0, treeheight_col = 0)
  
  save_pheatmap_png(hm_pl, paste0(outFile, '.png'), height =height_pl , width =width_pl)
  save_pheatmap_pdf(hm_pl, paste0(outFile, '.pdf'), height =height_pl , width =width_pl)
  
}


pheat_pl_gr <- function(mat, type_mat, height_pl = 10, width_pl = 7, color_df, outFile){
  
  coul <- rev(colorRampPalette(brewer.pal(11, "RdBu"))(100))

  tmp_mat <- as.matrix(mat[, !colnames(mat) %in% c('id', 'tissue')])
  
  val <- max(abs(tmp_mat))
  mat_breaks <- seq(-val, val, length.out = 100)
  
  tmp_mat[tmp_mat>=val] <- val
  tmp_mat[tmp_mat<=-val] <- -val
  rownames(tmp_mat) <- mat$id
  # id <- order(cl$gr)
  # P <- length(unique(cl$gr))
  # tmp_mat <- tmp_mat[id,]
  
  mat_row <- data.frame(tissue = mat$tissue)
  rownames(mat_row) <- rownames(tmp_mat)
  
  mat_colors_gr <- list(tissue = color_df$color)
  names(mat_colors_gr$tissue) <- unique(mat$tissue)
  
  
  hm_pl <- pheatmap(mat=tmp_mat, color=coul, show_colnames = T, show_rownames = T, annotation_names_col = F, annotation_names_row = F, 
                    cluster_rows=T, cluster_cols=F, border_color = NA,
                    annotation_colors = mat_colors_gr,
                    annotation_row = mat_row, drop_levels = TRUE, fontsize_row = 8, fontsize_col = 10, fontsize = 10, 
                    main =  sprintf("%s", type_mat),
                    cellwidth = 15, 
                    breaks = mat_breaks, treeheight_row = 0, treeheight_col = 0)
  
  save_pheatmap_png(hm_pl, paste0(outFile, '.png'), height =height_pl , width =width_pl)
  save_pheatmap_pdf(hm_pl, paste0(outFile, '.pdf'), height =height_pl , width =width_pl)
  
}


compute_reg_endopheno <- function(fmla, type_pheno, mat){
  
  
  if(type_pheno == 'CONTINUOUS'){
    
    res <- tryCatch(glm(fmla, data = mat, family = 'gaussian'),warning=function(...) NA, error=function(...) NA)
    if(is.list(res)){
       output <- coef(summary(res))[rownames(coef(summary(res))) == 'gr_id1',1:4]
       output[5] <- output[1]
       output[6] <-  output[1] + qnorm(0.025)*output[2]
       output[7] <-  output[1] + qnorm(0.975)*output[2]
    }else{
       output <- rep(NA, 7)
    }
  }else{
  
  if((type_pheno %in% c('CAT_SINGLE_UNORDERED', 'CAT_SINGLE_BINARY', 'CAT_MUL_BINARY_VAR')) | (type_pheno == 'CAT_ORD' & length(unique(na.omit(mat[, 'pheno']))) == 2)){
    
    if(!all(unique(na.omit(mat[, 'pheno']) %in% c(0,1)))){
      
      min_id <- which(mat[, 'pheno'] == min(mat[,'pheno'], na.rm = T))
      mat[min_id,  'pheno'] <- 0
      max_id <- which(mat[, 'pheno'] == max(mat[,  'pheno'], na.rm = T))
      mat[max_id,  'pheno'] <- 1
      
    }
    
    res <- tryCatch(glm(fmla, data = mat, family = 'binomial'),warning=function(...) NA, error=function(...) NA)
    if(is.list(res)){
    	output <- coef(summary(res))[rownames(coef(summary(res))) == 'gr_id1',1:4]
    	output[5] <- exp(output[1])
    	output[6] <- exp(output[1] + qnorm(0.025)*output[2])
    	output[7] <- exp(output[1] + qnorm(0.975)*output[2])
    }else{
	output <- rep(NA, 7)
    }
  }else{
  
  if(type_pheno == 'CAT_ORD' & length(unique(na.omit(mat[, 'pheno']))) > 2){
    
    mat$pheno <- factor(mat$pheno)
    output <- rep(NA, 7)

    res <- tryCatch(polr(fmla, data = mat, Hess=TRUE),warning=function(...) NA, error=function(...) NA)
    if(is.list(res)){
      if(!any(is.na(res$Hess))){
        ct <- coeftest(res)  
        output <- ct[rownames(ct) == 'gr_id1',1:4]
        output[5] <- exp(output[1])
        output[6] <- exp(output[1] + qnorm(0.025)*output[2])
        output[7] <- exp(output[1] + qnorm(0.975)*output[2])
      }
    }

    
  }else{
  output <- rep(NA, 7)
  }
  }
  }

  return(output)
  
}

# compute_reg_features <- function(fmla, mat){
#   
#   res <- tryCatch(glm(fmla, data = mat, family = 'gaussian'),warning=function(...) NA, error=function(...) NA)
#   if(is.list(res)){
#     output <- coef(summary(res))[rownames(coef(summary(res))) == 'gr_id1',1:4]
#   }else{
#     output <- rep(NA, 4)
#   }
#   
#   # res <- tryCatch(glm(fmla, data = mat, family = 'binomial'),warning=function(...) NA, error=function(...) NA)
#   # if(is.list(res)){
#   #   output <- coef(summary(res))[rownames(coef(summary(res))) == 'gene',1:4]
#   # }else{
#   #   output <- rep(NA, 4)
#   # }
#   
#   return(output)
#   
# }
