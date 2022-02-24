## load input data ##
load_input_matrix <- function(inputFile, sampleAnn, res_pval, split_tot, id_info){
  
  if(split_tot == 0){
    
    ### load score ###
    if(substr(inputFile, nchar(inputFile)-3, nchar(inputFile)) == '.txt'){
      tmp <- read.delim(inputFile, h=T, stringsAsFactors = F, check.names = F)
      if(type_data == 'tscore'){
        # correct sample names
        sampleID <- unname(sapply(colnames(tmp)[-1], function(x) strsplit(x, split = '.vs')[[1]][1]))
        colnames(tmp)[-1] <- sampleID
      }
      # filter out elements that are repeated twice:
      id_dup <- names(which(table(tmp[,1] ) > 1)) 
      tmp <- tmp[!tmp[,1] %in% id_dup, ]
      id_el <- intersect(tmp[,1], res_pval[, id_info])
      tmp <- tmp[match(id_el, tmp[,1]),]
      res_pval <- res_pval[match(id_el, res_pval[, id_info]),]
      elementID <- tmp[,1]
      # consider only samples in common
      common_samples <- intersect(sampleAnn$Individual_ID, colnames(tmp))
      sampleAnn <- sampleAnn[match(common_samples, sampleAnn$Individual_ID),]
      tmp <- tmp[, match(common_samples, colnames(tmp))]
      scoreMat <- as.matrix(t(tmp))
      colnames(scoreMat) <- elementID
    }else{  
      
      scoreMat <- get(load(inputFile))
      # filter out based on samples and ids
      id_el <- intersect(scoreMat[,1], res_pval[, id_info])
      scoreMat <- scoreMat[match(id_el,scoreMat[,1]), ]
      
      common_samples <- intersect(sampleAnn$Individual_ID, colnames(scoreMat))
      sampleAnn <- sampleAnn[match(common_samples, sampleAnn$Individual_ID),]
      scoreMat <- t(scoreMat[,match(common_samples,colnames(scoreMat))])
      
      rownames(scoreMat) <- common_samples
      colnames(scoreMat) <- id_el
      res_pval <- res_pval[match(id_el, res_pval[, id_info]),]
    } 
    
  }else{
    
    scoreMat_list <- vector(mode = 'list', length = split_tot)
    samplesID <- vector(mode = 'list', length = split_tot)
    elementID <- NULL
    
    for(i in 1:split_tot){
      
      print(i)
      if(file.exists(sprintf('%s%i.RData', inputFile, i))){
        tmp <- get(load(sprintf('%s%i.RData', inputFile, i)))
        elementID <- c(elementID,tmp[,1])
        samplesID[[i]] <- intersect(sampleAnn$Individual_ID, colnames(tmp))
        scoreMat_list[[i]] <- t(tmp[,match(samplesID[[i]],colnames(tmp))])
      }else{
        print(sprintf('split %i does not exist', i))
        split_tot <- split_tot - 1
      }
    }
    
    print(split_tot)
    # check samplesID always the same
    if(!all(table(unlist(samplesID)) == split_tot)){print('ERROR: wrong name annotations')}
    
    scoreMat <- do.call(cbind, scoreMat_list)
    colnames(scoreMat) <- elementID
    rm(scoreMat_list)
    
    # filter out elements that are repeated twice:
    id_dup <- names(which(table(colnames(scoreMat)) > 1)) 
    scoreMat <- scoreMat[, !colnames(scoreMat) %in% id_dup]
    
    id_el <- intersect(colnames(scoreMat),  res_pval[, id_info])
    scoreMat <- scoreMat[, match(id_el, colnames(scoreMat))]
    
    rownames(scoreMat) <- samplesID[[1]]
    # remove sample that have NAs
    id_s <- rowSums(is.na(scoreMat)) == 0
    if(!all(id_s)){scoreMat <- scoreMat[id_s, ]}
    
    common_samples <- intersect(sampleAnn$Individual_ID, rownames(scoreMat))
    sampleAnn <- sampleAnn[match(common_samples, sampleAnn$Individual_ID),]
    scoreMat <- scoreMat[match(common_samples,rownames(scoreMat)),]
    res_pval <- res_pval[match(id_el, res_pval[, id_info]),]
    
  }
  
  return(list(res_pval = res_pval, sampleAnn = sampleAnn, scoreMat = scoreMat))
  
}

## clumping ##
clumping_features <- function(res_pval, id_info, corr_feat, id_pval, corr_thr){

  feat_info <- res_pval[,c(id_info, id_pval)]
  feat_info <- feat_info[order(feat_info[,2], decreasing = F), ]
  corr_feat <- corr_feat[match(feat_info[,1], rownames(corr_feat)), 
                         match(feat_info[,1], colnames(corr_feat))]
  element_rm <- c()
  stop_cond <- F
  list_feat <- feat_info[,1]
  
  while(!stop_cond){
    
    id <- which(abs(corr_feat[,list_feat[1]])>corr_thr)
    if(length(id)>1){
      element_rm <- c(element_rm, names(id)[-which.max(feat_info[match(names(id),feat_info[,1]), 2])])
    }
    list_feat <- list_feat[!list_feat %in% names(id)]
    # print(length(list_feat))
    stop_cond <- length(list_feat) == 0
    
  }
  
  element_rm <- unique(element_rm)
  return(element_rm)
}

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
  tmp <- sapply(1:(nrow(score)-1), function(x) c(rep(0,x), (sigma_kNN[x] + sigma_kNN[(x+1):nrow(data_tot)] + euclDist[x, (x+1):nrow(data_tot)])/3))
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

pheat_pl_tscore <- function(mat_tscore, cl, info_feat_tscore = NA, test_feat_tscore = NA, pval_thr_est = 0.05, 
                            height_pl = 10, width_pl = 7, outFile, cap = NA, res_pl = 200){
  
  
  coul <- rev(colorRampPalette(brewer.pal(11, "RdBu"))(100))
  # cap mat if necessary
  tmp_mat <- as.matrix(mat_tscore)
  if(!is.na(cap)){
    val <- abs(cap)
  }else{
    val <- max(abs(tmp_mat))  
  }
  mat_breaks <- seq(-val, val, length.out = 100)
  
  tmp_mat[tmp_mat>=val] <- val
  tmp_mat[tmp_mat<=-val] <- -val
  
  # sort cl, match mat according cl
  id <- order(cl$gr)
  P <- length(unique(cl$gr))
  cl <- cl[id, ]
  tmp_mat <- tmp_mat[match(cl$id, rownames(tmp_mat)),]
  
  # order gene according location
  info_feat_tscore <- info_feat_tscore[order(info_feat_tscore$start_position), ]
  info_feat_tscore <- info_feat_tscore[order(as.numeric(sapply(info_feat_tscore$chrom, function(x) strsplit(x, split = 'chr')[[1]][2]))), ]
  keep_gene <- info_feat_tscore$external_gene_name
  tmp_mat <- tmp_mat[, match(keep_gene, colnames(tmp_mat)), drop = F]
  tmp_mat <- t(tmp_mat)
  chr_fact <- factor(info_feat_tscore$chrom, levels = unique(info_feat_tscore$chrom))
  test_feat_tscore <- test_feat_tscore[test_feat_tscore$feat %in% keep_gene, ,  drop = F]
  
  mat_colors_gr <- list(cluster = pal_d3(palette = 'category20')(P))
  names(mat_colors_gr$cluster) <- paste0('gr', 1:P)
  
  mat_colors_chr <- list(chrom = rep(c('#7C7C7C', '#C1C1C1'),length(unique(info_feat_tscore$chrom)))[1:length(unique(info_feat_tscore$chrom))])
  names(mat_colors_chr$chrom) <- unique(info_feat_tscore$chrom)
  
  zstat_col_fun = colorRamp2(c(min(c(info_feat_tscore$Zstat), na.rm = T), 0, max(c(info_feat_tscore$Zstat), na.rm = T)), 
                             c("blue","#F0F0F0", "red"))
  # add pvalue info for each group
  estimate_col_fun = colorRamp2(c(min(c(test_feat_tscore$estimates)), 0, max(c(test_feat_tscore$estimates))), 
                                c("#00677B", "#F0F0F0", "#BF443B"))
  lgd_est = Legend(title = "wilcoxon estimates", col = estimate_col_fun, 
                   at = round(c(seq(min(c(test_feat_tscore$estimates)), 0, length.out = 4), 
                                seq(0, max(c(test_feat_tscore$estimates)), length.out = 4)[-1]), digits = 2),
                   labels = as.character(round(c(seq(min(c(test_feat_tscore$estimates)), 0, length.out = 4), 
                                                 seq(0, max(c(test_feat_tscore$estimates)), length.out = 4)[-1]), digits = 2)))
  # and one for the significant p-values
  lgd_sig = Legend(pch = "*", type = "points", labels = sprintf("FDR pvalue < %s", as.character(pval_thr_est)))
  
  
  column_ha <- HeatmapAnnotation(cluster = anno_block(gp = gpar(fill = mat_colors_gr$cluster),
                                                      labels = names(mat_colors_gr$cluster),
                                                      labels_gp = gpar(col = "white", fontsize = 12,  fontface = "bold")))
  
  row_ha <- rowAnnotation(chrom = anno_block(gp = gpar(fill = mat_colors_chr$chrom),
                                             labels = factor(names(mat_colors_chr$chrom), levels = names(mat_colors_chr$chrom)),
                                             labels_gp = gpar(col = "black", fontsize = 10), 
                                             labels_rot = 0), 
                          zstat = info_feat_tscore$Zstat,
                          col = list(zstat = zstat_col_fun), 
                          annotation_label = list(zstat = sprintf('z-statistic %s', pheno_name)), 
                          annotation_name_gp = gpar(col = 'white'))
  
  df_pch <- list()
  df_est <- list()
  for(i in 1:P){
    tmp_gr <- test_feat_tscore[test_feat_tscore$comp == sprintf('gr%i_vs_all', i), ]
    df_est[[i]] <- tmp_gr$estimates[match(keep_gene, tmp_gr$feat)]
    is_sign <- tmp_gr$pval_corr[match(keep_gene, tmp_gr$feat)] < pval_thr_est
    df_pch[[i]] <- rep("*", length(is_sign))
    df_pch[[i]][!is_sign] = NA
  }
  df_est <- do.call(cbind, df_est)
  colnames(df_est) <- paste0('gr', 1:P)
  df_pch <- do.call(cbind, df_pch)
  colnames(df_pch) <- paste0('gr', 1:P)
  df_font <- matrix(rep("bold", P), nrow = 1)
  df_font <- as.data.frame(df_font)
  colnames(df_font) <- paste0('gr', 1:P)
  
  row_ha_gr <- rowAnnotation(gr = anno_simple(df_est, col = estimate_col_fun, pch = df_pch, border = T, pt_gp = gpar(fontface = df_font)), 
                             annotation_label = '', annotation_name_side = 'top', annotation_name_rot = 0, simple_anno_size_adjust = T)
  
  hm_pl <- Heatmap(tmp_mat, name = "scaled\nT-scores", col = coul, cluster_rows = FALSE, cluster_columns = FALSE,  show_column_names = F, 
                   top_annotation = column_ha, column_split = cl$gr, column_title = NULL,
                   row_names_side = "left", row_names_gp = gpar(fontsize = 10),
                   left_annotation = row_ha,  row_split  = factor(info_feat_tscore$chrom, levels = unique(info_feat_tscore$chrom)), row_title = NULL, row_gap = unit(0, "mm"),
                   right_annotation = row_ha_gr, 
                   border = TRUE, use_raster = T, show_heatmap_legend = T)
  tot_pl <- hm_pl
  
  
  # cohort info
  if('cohort' %in% colnames(cl)){
    cl$cohort <- factor(cl$cohort)
    color <- grDevices::colors()[grep('gr(a|e)y', grDevices::colors(), invert = T)]
    set.seed(24)
    color_cohort <- sample(color, length(unique(cl$cohort)))
    names(color_cohort) <- levels(cl$cohort)
    ha <- HeatmapAnnotation(cohort = factor(cl$cohort), annotation_label = list(cohort = 'Cohort'), col = list(cohort = color_cohort), border = T)  
  }
  
  if('cohort' %in% colnames(cl)){
    ht_list <- tot_pl %v% ha 
  }else{
    ht_list <- tot_pl
  }
  side_par <- round(max(sapply(rownames(tmp_mat), nchar))*0.9)
  
  png(file=paste0(outFile, '.png'), res = res_pl, units = 'in', width = width_pl, height = height_pl)
  draw(ht_list , annotation_legend_list = list(lgd_est, lgd_sig), merge_legend = T, auto_adjust = T, padding = unit(c(2, 2 + side_par, 2, 2), "mm"))
  dev.off()
  
  pdf(file=paste0(outFile, '.pdf'), width = width_pl, height = height_pl)
  draw(ht_list , annotation_legend_list = list(lgd_est, lgd_sig), merge_legend = T, auto_adjust = T, padding = unit(c(2, 2 + side_par, 2, 2), "mm"))
  dev.off()
  
}


pheat_pl_path <- function(mat, cl, info_feat = NA, test_feat = NA, height_pl = 10, width_pl = 7, outFile, cap = NA, res_pl = 200){
  
  coul <- rev(colorRampPalette(brewer.pal(11, "RdBu"))(100))
  # cap mat if necessary
  tmp_mat <- as.matrix(mat)
  if(!is.na(cap)){
    val <- abs(cap)
  }else{
    val <- max(abs(tmp_mat))  
  }
  mat_breaks <- seq(-val, val, length.out = 100)
  
  tmp_mat[tmp_mat>=val] <- val
  tmp_mat[tmp_mat<=-val] <- -val
  
  # sort cl, match mat according cl
  id <- order(cl$gr)
  P <- length(unique(cl$gr))
  cl <- cl[id, ]
  tmp_mat <- tmp_mat[match(cl$id, rownames(tmp_mat)),]
  
  # order gene according location
  keep_path <- info_feat[,1]
  tmp_mat <- tmp_mat[, match(keep_path, colnames(tmp_mat))]
  tmp_mat <- t(tmp_mat)
  test_feat <- test_feat[test_feat$feat %in% keep_path, ]
  
  # mat_col <- data.frame(group = paste0('gr', cl$gr))
  # rownames(mat_col) <- colnames(tmp_mat)
  
  mat_colors_gr <- list(cluster = pal_d3(palette = 'category20')(P))
  names(mat_colors_gr$cluster) <- paste0('gr', 1:P)
  
  column_ha <- HeatmapAnnotation(cluster = anno_block(gp = gpar(fill = mat_colors_gr$cluster),
                                                      labels = names(mat_colors_gr$cluster),
                                                      labels_gp = gpar(col = "white", fontsize = 12,  fontface = "bold")))
  
  ngenes_col_fun = colorRamp2(c(min(info_feat$ngenes_tscore),max(info_feat$ngenes_tscore)), c("white", "#7a2048"))
  perc_col_fun = colorRamp2(c(min(info_feat$ngenes_tscore/info_feat$ngenes_path), max(info_feat$ngenes_tscore/info_feat$ngenes_path)), c("white", "#316879"))
  row_ha <- rowAnnotation(n_genes = info_feat$ngenes_tscore, perc = info_feat$ngenes_tscore/info_feat$ngenes_path, 
                          annotation_label = list(n_genes = 'n. genes', perc = '% tot genes'), 
                          col = list(n_genes = ngenes_col_fun, perc = perc_col_fun))
  
  # add pvalue info for each group
  estimate_col_fun = colorRamp2(c(min(test_feat$estimates), 0, max(test_feat$estimates)), c("darkgreen", "#F0F0F0", "darkorange"))
  df_pch <- list()
  df_est <- list()
  for(i in 1:P){
    tmp_gr <- test_feat[test_feat$comp == sprintf('gr%i_vs_all', i), ]
    df_est[[i]] <- tmp_gr$estimates[match(keep_path, tmp_gr$feat)]
    is_sign <- tmp_gr$pval_corr[match(keep_path, tmp_gr$feat)] < 0.05
    df_pch[[i]] <- rep("*", length(is_sign))
    df_pch[[i]][!is_sign] = NA
  }
  df_est <- do.call(cbind, df_est)
  colnames(df_est) <- paste0('gr', 1:P)
  df_pch <- do.call(cbind, df_pch)
  colnames(df_pch) <- paste0('gr', 1:P)
  df_font <- matrix(rep("bold", P), nrow = 1)
  df_font <- as.data.frame(df_font)
  colnames(df_font) <- paste0('gr', 1:P)
  
  lgd_est = Legend(title = "wilcoxon\nestimates", col = estimate_col_fun, at = round(c(seq(min(test_feat$estimates), 0, length.out = 4), seq(0, max(test_feat$estimates), length.out = 4)[-1]), digits = 2),
                   labels = as.character(round(c(seq(min(test_feat$estimates), 0, length.out = 4), seq(0, max(test_feat$estimates), length.out = 4)[-1]), digits = 2)))
  # and one for the significant p-values
  lgd_sig = Legend(pch = "*", type = "points", labels = "FDR pvalue < 0.05")
  
  hm_pl <- Heatmap(tmp_mat, name = "scaled\nPathway-scores", col = coul, cluster_rows = F, cluster_columns = FALSE,  show_column_names = F, 
                   top_annotation = column_ha, column_split = cl$gr, column_title = NULL, 
                   row_names_side = "left", row_names_gp = gpar(fontsize = 10),
                   left_annotation = row_ha,
                   border = TRUE, use_raster = T)
  tot_pl <- hm_pl
  for(i in 1:P){
    gr_row_ha <- rowAnnotation(gr = anno_simple(df_est[,i], col = estimate_col_fun, pch = df_pch[,i], border = T, pt_gp = gpar(fontface = df_font[,i])), 
                               annotation_label = paste0('gr', i), annotation_name_side = 'top', annotation_name_rot = 0, simple_anno_size_adjust = T)
    tot_pl <- tot_pl + gr_row_ha
  }
  side_par <- round(max(sapply(rownames(tmp_mat), nchar))*0.9)
  
  png(file=paste0(outFile, '.png'), res = res_pl, units = 'in', width = width_pl, height = height_pl)
  draw(tot_pl , annotation_legend_list = list(lgd_est, lgd_sig), merge_legend = T, auto_adjust = T, padding = unit(c(2, 2 + side_par, 2, 2), "mm"))
  dev.off()
  
  pdf(file=paste0(outFile, '.pdf'), width = width_pl, height = height_pl)
  draw(tot_pl , annotation_legend_list = list(lgd_est, lgd_sig), merge_legend = T, auto_adjust = F)
  dev.off()
  
  # note: not possible to change label color of the annotation (use illustrator)
  
}

###################################################################################################

pheat_pl_tot <- function(pheno_name, mat_tscore, info_feat_tscore, test_feat_tscore, pval_thr_est = 0.05, 
                         mat_path, info_feat_path, test_feat_path,
                         cl, height_pl = 10, width_pl = 7, outFile, cap = NA, res_pl = 200){
  
  ################
  #### Tscore ####
  ################
  
  coul <- rev(colorRampPalette(brewer.pal(11, "RdBu"))(100))
  # cap mat if necessary
  tmp_mat <- as.matrix(mat_tscore)
  if(!is.na(cap)){
    val <- abs(cap)
  }else{
    val <- max(abs(tmp_mat))  
  }
  mat_breaks <- seq(-val, val, length.out = 100)
  
  tmp_mat[tmp_mat>=val] <- val
  tmp_mat[tmp_mat<=-val] <- -val
  
  # sort cl, match mat according cl
  id <- order(cl$gr)
  P <- length(unique(cl$gr))
  cl <- cl[id, ]
  tmp_mat <- tmp_mat[match(cl$id, rownames(tmp_mat)),]
  
  # order gene according location
  info_feat_tscore <- info_feat_tscore[order(info_feat_tscore$start_position), ]
  info_feat_tscore <- info_feat_tscore[order(as.numeric(sapply(info_feat_tscore$chrom, function(x) strsplit(x, split = 'chr')[[1]][2]))), ]
  keep_gene <- info_feat_tscore$external_gene_name
  tmp_mat <- tmp_mat[, match(keep_gene, colnames(tmp_mat)), drop = F]
  tmp_mat <- t(tmp_mat)
  chr_fact <- factor(info_feat_tscore$chrom, levels = unique(info_feat_tscore$chrom))
  test_feat_tscore <- test_feat_tscore[test_feat_tscore$feat %in% keep_gene, ,  drop = F]
  
  mat_colors_gr <- list(cluster = pal_d3(palette = 'category20')(P))
  names(mat_colors_gr$cluster) <- paste0('gr', 1:P)
  
  mat_colors_chr <- list(chrom = rep(c('#7C7C7C', '#C1C1C1'),length(unique(info_feat_tscore$chrom)))[1:length(unique(info_feat_tscore$chrom))])
  names(mat_colors_chr$chrom) <- unique(info_feat_tscore$chrom)
  
  zstat_col_fun = colorRamp2(c(min(c(info_feat_tscore$Zstat, info_feat_path$Zstat), na.rm = T), 0, max(c(info_feat_tscore$Zstat, info_feat_path$Zstat), na.rm = T)), 
                             c("blue","#F0F0F0", "red"))
  # add pvalue info for each group
  estimate_col_fun = colorRamp2(c(min(c(test_feat_tscore$estimates, test_feat_path$estimates)), 0, max(c(test_feat_tscore$estimates, test_feat_path$estimates))), 
                                c("#00677B", "#F0F0F0", "#BF443B"))
  lgd_est = Legend(title = "wilcoxon estimates", col = estimate_col_fun, 
                   at = round(c(seq(min(c(test_feat_tscore$estimates, test_feat_path$estimates)), 0, length.out = 4), 
                                seq(0, max(c(test_feat_tscore$estimates, test_feat_path$estimates)), length.out = 4)[-1]), digits = 2),
                   labels = as.character(round(c(seq(min(c(test_feat_tscore$estimates, test_feat_path$estimates)), 0, length.out = 4), 
                                                 seq(0, max(c(test_feat_tscore$estimates, test_feat_path$estimates)), length.out = 4)[-1]), digits = 2)))
  # and one for the significant p-values
  lgd_sig = Legend(pch = "*", type = "points", labels = sprintf("FDR pvalue < %s", as.character(pval_thr_est)))
  
  
  column_ha <- HeatmapAnnotation(cluster = anno_block(gp = gpar(fill = mat_colors_gr$cluster),
                                                      labels = names(mat_colors_gr$cluster),
                                                      labels_gp = gpar(col = "white", fontsize = 12,  fontface = "bold")))
  
  row_ha <- rowAnnotation(chrom = anno_block(gp = gpar(fill = mat_colors_chr$chrom),
                                             labels = factor(names(mat_colors_chr$chrom), levels = names(mat_colors_chr$chrom)),
                                             labels_gp = gpar(col = "black", fontsize = 10), 
                                             labels_rot = 0), 
                          zstat = info_feat_tscore$Zstat,
                          col = list(zstat = zstat_col_fun), 
                          annotation_label = list(zstat = sprintf('z-statistic %s', pheno_name)), 
                          annotation_name_gp = gpar(col = 'white'))
  
  df_pch <- list()
  df_est <- list()
  for(i in 1:P){
    tmp_gr <- test_feat_tscore[test_feat_tscore$comp == sprintf('gr%i_vs_all', i), ]
    df_est[[i]] <- tmp_gr$estimates[match(keep_gene, tmp_gr$feat)]
    is_sign <- tmp_gr$pval_corr[match(keep_gene, tmp_gr$feat)] < pval_thr_est
    df_pch[[i]] <- rep("*", length(is_sign))
    df_pch[[i]][!is_sign] = NA
  }
  df_est <- do.call(cbind, df_est)
  colnames(df_est) <- paste0('gr', 1:P)
  df_pch <- do.call(cbind, df_pch)
  colnames(df_pch) <- paste0('gr', 1:P)
  df_font <- matrix(rep("bold", P), nrow = 1)
  df_font <- as.data.frame(df_font)
  colnames(df_font) <- paste0('gr', 1:P)
  
  row_ha_gr <- rowAnnotation(gr = anno_simple(df_est, col = estimate_col_fun, pch = df_pch, border = T, pt_gp = gpar(fontface = df_font)), 
                             annotation_label = '', annotation_name_side = 'top', annotation_name_rot = 0, simple_anno_size_adjust = T)
  
  hm_pl <- Heatmap(tmp_mat, name = "scaled\nT-scores", col = coul, cluster_rows = FALSE, cluster_columns = FALSE,  show_column_names = F, 
                   top_annotation = column_ha, column_split = cl$gr, column_title = NULL,
                   row_names_side = "left", row_names_gp = gpar(fontsize = 10),
                   left_annotation = row_ha,  row_split  = factor(info_feat_tscore$chrom, levels = unique(info_feat_tscore$chrom)), row_title = NULL, row_gap = unit(0, "mm"),
                   right_annotation = row_ha_gr, 
                   border = TRUE, use_raster = T, show_heatmap_legend = F)
  tot_pl <- hm_pl
  
  
  # cohort info
  if('cohort' %in% colnames(cl)){
    cl$cohort <- factor(cl$cohort)
    color <- grDevices::colors()[grep('gr(a|e)y', grDevices::colors(), invert = T)]
    set.seed(24)
    color_cohort <- sample(color, length(unique(cl$cohort)))
    names(color_cohort) <- levels(cl$cohort)
    ha <- HeatmapAnnotation(cohort = factor(cl$cohort), annotation_label = list(cohort = 'Cohort'), col = list(cohort = color_cohort), border = T)  
  }
  
  ############################################
  ################# pathway ##################
  ############################################
  
  # cap mat if necessary
  tmp_mat <- as.matrix(mat_path)
  if(!is.na(cap)){
    val <- abs(cap)
  }else{
    val <- max(abs(tmp_mat))  
  }
  mat_breaks <- seq(-val, val, length.out = 100)
  
  tmp_mat[tmp_mat>=val] <- val
  tmp_mat[tmp_mat<=-val] <- -val
  tmp_mat <- tmp_mat[match(cl$id, rownames(tmp_mat)),]
  
  # order gene according location
  keep_path <- info_feat_path[,1]
  tmp_mat <- tmp_mat[, match(keep_path, colnames(tmp_mat)), drop = F]
  tmp_mat <- t(tmp_mat)
  test_feat_path <- test_feat_path[test_feat_path$feat %in% keep_path, , drop = F]
  
  # mat_col <- data.frame(group = paste0('gr', cl$gr))
  # rownames(mat_col) <- colnames(tmp_mat)
  
  mat_colors_gr <- list(cluster = pal_d3(palette = 'category20')(P))
  names(mat_colors_gr$cluster) <- paste0('gr', 1:P)
  
  column_ha <- HeatmapAnnotation(cluster = anno_block(gp = gpar(fill = mat_colors_gr$cluster),
                                                      labels = names(mat_colors_gr$cluster),
                                                      labels_gp = gpar(col = "white", fontsize = 12,  fontface = "bold")))
  
  if(length(unique(info_feat_path$ngenes_tscore)) == 1){
    ngenes_col_fun = colorRamp2(c(0,max(info_feat_path$ngenes_tscore)), c("white", "#035F1D"))
    perc_col_fun = colorRamp2(c(0, max(info_feat_path$ngenes_tscore/info_feat_path$ngenes_path)), c("white", "#316879"))
  }else{
    ngenes_col_fun = colorRamp2(c(min(info_feat_path$ngenes_tscore),max(info_feat_path$ngenes_tscore)), c("white", "#035F1D"))
    perc_col_fun = colorRamp2(c(min(info_feat_path$ngenes_tscore/info_feat_path$ngenes_path), max(info_feat_path$ngenes_tscore/info_feat_path$ngenes_path)), c("white", "#316879"))
  }
  
  row_ha <- rowAnnotation(n_genes = info_feat_path$ngenes_tscore, perc = info_feat_path$ngenes_tscore/info_feat_path$ngenes_path, 
                          zstat = info_feat_path$Zstat,
                          annotation_label = list(n_genes = 'n. genes', perc = '% tot genes', zstat = sprintf('z-statistic %s', pheno_name)), 
                          col = list(n_genes = ngenes_col_fun, perc = perc_col_fun, zstat = zstat_col_fun))
  
  # add pvalue info for each group
  df_pch <- list()
  df_est <- list()
  for(i in 1:P){
    tmp_gr <- test_feat_path[test_feat_path$comp == sprintf('gr%i_vs_all', i), ]
    df_est[[i]] <- tmp_gr$estimates[match(keep_path, tmp_gr$feat)]
    is_sign <- tmp_gr$pval_corr[match(keep_path, tmp_gr$feat)] < pval_thr_est
    df_pch[[i]] <- rep("*", length(is_sign))
    df_pch[[i]][!is_sign] = NA
  }
  df_est <- do.call(cbind, df_est)
  colnames(df_est) <- paste0('gr', 1:P)
  df_pch <- do.call(cbind, df_pch)
  colnames(df_pch) <- paste0('gr', 1:P)
  df_font <- matrix(rep("bold", P), nrow = 1)
  df_font <- as.data.frame(df_font)
  colnames(df_font) <- paste0('gr', 1:P)
  
  row_ha_gr <- rowAnnotation(gr = anno_simple(df_est, col = estimate_col_fun, pch = df_pch, border = T, pt_gp = gpar(fontface = df_font)), 
                             annotation_label = 'wilcoxon\nestimates', annotation_name_side = 'bottom', annotation_name_rot = 0, simple_anno_size_adjust = T)
  font_path <- rep('plain', nrow(tmp_mat))
  font_path[info_feat_path$impr] <- 'bold'
  
  hm_pl <- Heatmap(tmp_mat, name = "scaled scores", col = coul, cluster_rows = T, cluster_columns = FALSE,  show_row_dend = F, show_column_names = F, 
                   bottom_annotation = column_ha, column_split = cl$gr, column_title = NULL, 
                   row_names_side = "left", row_names_gp = gpar(fontsize = 10, col = 'black', fontface = font_path),
                   left_annotation = row_ha, right_annotation = row_ha_gr, 
                   border = TRUE, use_raster = T)
  tot_pl_p <- hm_pl
  side_par <- round(max(sapply(rownames(tmp_mat), nchar))*0.9)
  
  ########### combine #############
  if('cohort' %in% colnames(cl)){
    ht_list <- tot_pl %v% ha %v% tot_pl_p
  }else{
    ht_list <- tot_pl %v% tot_pl_p
  }
  
  png(file=paste0(outFile, '.png'), res = res_pl, units = 'in', width = width_pl, height = height_pl)
  draw(ht_list , annotation_legend_list = list(lgd_est, lgd_sig), merge_legend = T, auto_adjust = T, padding = unit(c(2, 2 + side_par, 2, 2), "mm"))
  dev.off()
  
  pdf(file=paste0(outFile, '.pdf'), width = width_pl, height = height_pl)
  draw(ht_list , annotation_legend_list = list(lgd_est, lgd_sig), merge_legend = T, auto_adjust = T, padding = unit(c(2, 2 + side_par, 2, 2), "mm"))
  dev.off()
  
}


###################################################################################################


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

compute_reg_endopheno_lmm <- function(fmla, type_pheno, mat){
  
  res <- tryCatch(lmer(fmla, data = mat),warning=function(...) NULL, error=function(...) NULL)
  if(!is.null(res)){
    ci <- confint.merMod(res, level = 0.95)
    ct <- coef(summary(res))  
    output <- data.frame(beta = ct[rownames(ct) == 'gr_id1','Estimate'], se_beta = ct[rownames(ct) == 'gr_id1','Std. Error'], z = ct[rownames(ct) == 'gr_id1','t value'])
    output$CI_low <- ci['gr_id1',1]
    output$CI_up <- ci['gr_id1',2]
  }else{
    output <- data.frame(beta = c(), se_beta = c(), z = c(), CI_low = c(), CI_up = c())
  }
  return(output)
}


compute_reg_endopheno_multi <- function(fmla, type_pheno, mat, cov_int){
  
  
  if(type_pheno == 'CONTINUOUS'){
    
    res <- tryCatch(glm(fmla, data = mat, family = 'gaussian'),warning=function(...) NA, error=function(...) NA)
    if(is.list(res)){
      output <- coef(summary(res))[match(cov_int, rownames(coef(summary(res)))) ,1:4, drop = F]
      output <- cbind(output, output[,1],  output[,1] + qnorm(0.025)*output[,2], output[,1] + qnorm(0.975)*output[,2])
      colnames(output) <- c('beta', 'se_beta', 'z', 'pvalue', 'OR_or_beta', 'CI_low', 'CI_up')
      output <- as.data.frame(output)
      rownames(output) <- cov_int
    }else{
      output <- NULL
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
        output <- coef(summary(res))[match(cov_int, rownames(coef(summary(res)))),1:4,drop = F]
        output <- cbind(output, exp(output[,1]),  exp(output[,1] + qnorm(0.025)*output[,2]),  exp(output[,1] + qnorm(0.975)*output[,2]))
        colnames(output) <- c('beta', 'se_beta', 'z', 'pvalue', 'OR_or_beta', 'CI_low', 'CI_up')
        output <- as.data.frame(output)
        rownames(output) <- cov_int
      }else{
        output <- NULL
      }
    }else{
      
      if(type_pheno == 'CAT_ORD' & length(unique(na.omit(mat[, 'pheno']))) > 2){
        
        mat$pheno <- factor(mat$pheno)
        output <- NULL
        
        res <- tryCatch(polr(fmla, data = mat, Hess=TRUE),warning=function(...) NA, error=function(...) NA)
        if(is.list(res)){
          if(!any(is.na(res$Hess))){
            ct <- coeftest(res)  
            output <- ct[match(cov_int, rownames(ct)), 1:4, drop = F]
            output <- cbind(output, exp(output[,1]),  exp(output[,1] + qnorm(0.025)*output[,2]),  exp(output[,1] + qnorm(0.975)*output[,2]))
            colnames(output) <- c('beta', 'se_beta', 'z', 'pvalue', 'OR_or_beta', 'CI_low', 'CI_up')
            output <- as.data.frame(output)
            rownames(output) <- cov_int
          }
        }
        
      }else{
        output <- NULL
      }
    }
  }
  
  return(output)
  
}


#### meta analysis endophenotypes ####
meta_analysis_res <- function(beta, se_beta, thr_het = 0.001, type_pheno = NULL){
  
  if(all(beta == 0)){
    df <- data.frame(beta = NA, se_beta = NA, z = NA, pvalue = NA, Cochran_stat = NA, Cochran_pval = NA, model = NA, OR_or_Beta = NA, CI_low = NA, CI_up = NA)
  }else{
    
    model <- 'fixed'
    w <- 1/(se_beta)^2
    beta_all <- sum(beta*w, na.rm = T)/sum(w, na.rm = T)
    se_all <- sqrt(1/sum(w, na.rm = T))
    z_all <- beta_all/se_all
    p_all <- 2*pnorm(-abs(z_all), 0,1)
    Q <- sum(w*((beta_all - beta)^2), na.rm = T)
    Q_pval <- pchisq(Q, df=length(se_beta[!is.na(se_beta)])-1, lower.tail=FALSE) # null hypothesis consistency
    
    if(Q_pval<=thr_het){
      
      tau2 <- max(0, (Q-length(beta)+1)/(sum(w, na.rm = T) - (sum(w^2, na.rm = T)/sum(w, na.rm = T))))
      w_new <- 1/(tau2 + (se_beta)^2)
      se_all <- sqrt(1/sum(w_new, na.rm = T))
      beta_all <- sum(beta*w_new, na.rm = T)/sum(w_new, na.rm = T)
      z_all <- beta_all/se_all
      p_all <- 2*pnorm(-abs(z_all), 0,1)
      model <- 'random'
    }
    
    df <- data.frame(beta = beta_all, se_beta = se_all, z = z_all, pvalue = p_all, Cochran_stat = Q, Cochran_pval = Q_pval, model = model)
    if(!is.null(type_pheno)){
      
      if(type_pheno == 'CONTINUOUS'){
        df$OR_or_Beta <- df$beta
        df$CI_low <-  df$beta + qnorm(0.025)*df$se_beta
        df$CI_up <-   df$beta + qnorm(0.975)*df$se_beta
      }else{
        df$OR_or_Beta <- exp(df$beta)
        df$CI_low <- exp(df$beta + qnorm(0.025)*df$se_beta)
        df$CI_up <- exp(df$beta + qnorm(0.975)*df$se_beta)
      }
    }
  }
  return(df)
  
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


# function merge loci
merge_loci_genes <- function(gene_table, cis_size, bp_loci, tissue = 'combined'){
  
  tmp <- gene_table
  tmp_loci <- data.frame(chrom = c(), start = c(), end = c(), ngenes = c(), gene = c(), tissue = c())
  
  # divide per chr
  chr_id <- unique(tmp$chrom)
  tmp_chr <- lapply(chr_id, function(x) tmp[tmp$chrom == x,])
  
  for(j in 1:length(chr_id)){
    print(j)
    if(nrow(tmp_chr[[j]]) == 1){
      
      tmp_loci <- rbind(tmp_loci, data.frame(chrom = chr_id[j], start = tmp_chr[[j]]$TSS_start - cis_size, end = tmp_chr[[j]]$TSS_start + cis_size, 
                                             ngenes = 1, gene = tmp_chr[[j]]$external_gene_name, ensembl_gene = tmp_chr[[j]]$ensembl_gene_id,
                                             tissue = tissue))  
    }else{
      
      tmp_chr[[j]] <- tmp_chr[[j]][order(tmp_chr[[j]]$TSS_start), ]
      reg_gene <- data.frame(start = tmp_chr[[j]]$TSS_start - cis_size,  end = tmp_chr[[j]]$TSS_start + cis_size)
      merg_cond <- sapply(reg_gene$end, function(x) abs(x-reg_gene$start) < bp_loci) # the end of the second genes is close to the start of the first gene 1Mb
      
      merge_pos <- lapply(1:nrow(merg_cond), function(x) which(merg_cond[x,]))
      merge_pos_vect <- sapply(merge_pos, function(x) paste0(x, collapse = ','))
      merge_pos_vect <- merge_pos_vect[!duplicated(merge_pos_vect)]
      
      merge_pos <- lapply(merge_pos_vect, function(x) as.numeric(strsplit(x, split = ',')[[1]]))
      new_merge_pos <- list()
      all_merg <- F
      it <- 0
      
      if(length(merge_pos)>1){
        while(!all_merg){
          
          it <- it+1
          # print(it)
          
          for(l in 1:(length(merge_pos)-1)){
            
            if(!all(is.na(merge_pos[[l]]))){
              
              if(all(!merge_pos[[l]] %in% merge_pos[[l+1]])){
                new_merge_pos <- list.append(new_merge_pos, merge_pos[[l]])
              }else{
                if(!(all(merge_pos[[l]] %in% merge_pos[[l+1]]) | all(merge_pos[[l+1]] %in% merge_pos[[l]]))){
                  new_merge_pos <- list.append(new_merge_pos, unique(c(merge_pos[[l]], merge_pos[[l+1]])))
                }else{
                  if(all(merge_pos[[l+1]] %in% merge_pos[[l]])){
                    merge_pos[[l+1]] <- NA
                    new_merge_pos <- list.append(new_merge_pos, merge_pos[[l]])
                  }
                }
              }
              
            }
            
          }
          
          new_merge_pos <- list.append(new_merge_pos, merge_pos[[length(merge_pos)]])
          
          all_merg <- all(!duplicated(unlist(new_merge_pos)))
          merge_pos <- new_merge_pos
          new_merge_pos <- list() 
          
        }
        
        # remove NA
        merge_pos <- merge_pos[!sapply(merge_pos, function(x) all(is.na(x)))]
      }
      tmp_res <-  lapply(merge_pos, function(x) data.frame(chrom = chr_id[j], start = min(tmp_chr[[j]]$TSS_start[x] - cis_size), 
                                                           end = max(tmp_chr[[j]]$TSS_start[x] + cis_size), 
                                                           ngenes = length(x),
                                                           gene = paste0(unique(tmp_chr[[j]]$external_gene_name[x]), collapse = ','), 
                                                           ensembl_gene = paste0(unique(tmp_chr[[j]]$ensembl_gene_id[x]), collapse = ','), 
                                                           tissue = tissue))
      tmp_loci <-  rbind(tmp_loci,do.call(rbind, tmp_res))
      
    }
    
  }
  tmp_loci$start[tmp_loci$start < 0] <- 0
  tmp_loci$loci_id <- paste0(tmp_loci$chrom,':',round(tmp_loci$start/1000000, digits = 1), '-', round(tmp_loci$end/1000000, digits = 1), 'Mb')
  tmp_loci$loci_complete <- paste0(tmp_loci$chrom,':',tmp_loci$start,'-',tmp_loci$end)
  
  return(tmp_loci)
}

# function to check with category ordinal to remove (used in GLMassociation)
remove_pheno_ordinal <- function(pheno_df, group, thr){
  
  name_pheno <- colnames(pheno_df)
  pheno_rm <- c()
  
  for(i in 1:ncol(pheno_df)){
    
    min_p <- min(pheno_df[,i], na.rm = T)
    n_base_gr0 <- sum(pheno_df[group == 0,i] == min_p, na.rm = T)
    n_notbase_gr0 <- sum(pheno_df[group == 0,i] > min_p, na.rm = T)
    n_base_gr1 <- sum(pheno_df[group == 1, i] == min_p, na.rm = T)
    n_notbase_gr1 <- sum(pheno_df[group == 1, i] > min_p, na.rm = T)
    
    if(any(c(n_base_gr0, n_base_gr1, n_notbase_gr0, n_notbase_gr1) < thr)){
      pheno_rm <- c(pheno_rm, name_pheno[i])
    }
  }
  
  return(pheno_rm)
  
}

