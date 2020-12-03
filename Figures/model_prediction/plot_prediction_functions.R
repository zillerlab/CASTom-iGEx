# function for prediction plot (to be used across SCZ, CAD and Asthma)

## compute correlation matrix of zscore ##
create_cor <- function(tissues_name, res, id_z){
  
  cor_res <- diag(rep(1, length(tissues_name)), ncol = length(tissues_name))
  corpval_res <- diag(rep(1, length(tissues_name)),ncol = length(tissues_name))
  percint_res <- diag(rep(1, length(tissues_name)),ncol = length(tissues_name))
  
  for(i in 1:(length(tissues_name)-1)){
    
    # print(i)
    
    ti <- res[res$tissue  %in%  tissues_name[i],]
    trest <- lapply((i+1):length(tissues_name), function(x) res[res$tissue  %in%  tissues_name[x],])
    gene_int <- lapply(trest, function(x) intersect(ti[,1], x[,1]))
    gene_union <- lapply(trest, function(x) union(ti[,1], x[,1]))
    
    # match
    ti_match <- lapply(gene_int, function(x) ti[match(x, ti[,1]), ])
    trest_match <- mapply(function(x,y) x[match(y, x[,1]), ], x = trest, y = gene_int, SIMPLIFY = F)
    
    percint_res[i, (i+1):length(tissues_name)] <- mapply(function(x,y) length(x)/length(y), x = gene_int, y = gene_union)
    
    tmp <- mapply(function(x,y) cor.test(x[,id_z],y[,id_z], method = 'spearman'), x = ti_match, y = trest_match, SIMPLIFY = F)
    
    cor_res[i, (i+1):length(tissues_name)] <- sapply(tmp, function(x) x$estimate)
    corpval_res[i, (i+1):length(tissues_name)] <- sapply(tmp, function(x) x$p.value)
    
  }
  cor_res <- cor_res + t(cor_res) - diag(diag(cor_res))
  corpval_res <- corpval_res + t(corpval_res) - diag(diag(corpval_res))
  percint_res <- percint_res + t(percint_res) - diag(diag(percint_res))
  
  rownames(cor_res) <- colnames(cor_res) <- rownames(corpval_res) <- colnames(corpval_res) <- tissues_name
  rownames(percint_res) <- colnames(percint_res) <- tissues_name
  
  return(list(cor = cor_res, pval = corpval_res, perc = percint_res))
  
}

# plots
pl_corr <- function(res_cor, type_mat, type_dat, tissues_name, df_color, outFold, width_pl = 10,height_pl = 7 ){
  
  # correlation
  diag(res_cor$cor) <- NA
  col <- colorRampPalette(brewer.pal(9, 'Oranges'))(100)
  
  ord <- corrMatOrder(res_cor$cor, order="hclust", hclust.method = 'ward.D')
  newcolours <- df_color$color[match(tissues_name, df_color$tissues)][ord]
  title_pl <- sprintf('%s %s Spearman correlation z stat', type_mat, type_dat)
  
  pdf(file = sprintf('%s/corr_zscore_%s_%s.pdf', outFold, type_mat, type_dat), width = width_pl, height = height_pl, compress = F)
  corrplot(res_cor$cor, type="upper", order = 'hclust', hclust.method = 'ward.D',
           tl.col = newcolours, title = title_pl,
           col = c(col), method = 'color', tl.srt=45, cl.align.text='c',
           addCoef.col = "black",na.label = 'square', na.label.col = 'darkgrey', is.corr = F, number.cex=0.8, mar = c(0,0,5,0))
  dev.off()
  
  png(file = sprintf('%s/corr_zscore_%s_%s.png', outFold,  type_mat, type_dat), units = 'in', width = width_pl, height = height_pl, res = 300)
  corrplot(res_cor$cor, type="upper", order = 'hclust', hclust.method = 'ward.D',
           tl.col = newcolours, title = title_pl,
           col = c(col), method = 'color', tl.srt=45, cl.align.text='c',
           addCoef.col = "black",na.label = 'square', na.label.col = 'darkgrey', is.corr = F, number.cex=0.8, mar = c(0,0,5,0))
  dev.off()
  
  # intersection
  diag(res_cor$perc) <- NA
  col <- colorRampPalette(brewer.pal(9, 'Greens'))(100)
  
  title_pl <- sprintf('%s %s percentage common', type_mat, type_dat)
  
  pdf(file = sprintf('%s/perc_zscore_%s_%s.pdf', outFold, type_mat, type_dat), width = width_pl, height = height_pl, compress = F)
  corrplot(res_cor$perc[ord,ord], type="lower", 
           tl.col = newcolours, title = title_pl,
           col = c(col), method = 'color', tl.srt=45, cl.align.text='c',
           addCoef.col = "black",na.label = 'square', na.label.col = 'darkgrey', is.corr = F, number.cex=0.8, mar = c(0,0,5,0))
  dev.off()
  
  png(file = sprintf('%s/perc_zscore_%s_%s.png', outFold,  type_mat, type_dat), units = 'in', width = width_pl, height = height_pl, res = 300)
  corrplot(res_cor$perc[ord,ord], type="lower", 
           tl.col = newcolours, title = title_pl,
           col = c(col), method = 'color', tl.srt=45, cl.align.text='c',
           addCoef.col = "black",na.label = 'square', na.label.col = 'darkgrey', is.corr = F, number.cex=0.8, mar = c(0,0,5,0))
  dev.off()
  
}

#### number of associated element ####

creat_dfnsign <- function(tissues_name, res, id_pval_corr, pval_FDR, df_color, id_pval_corr_tot){
  
  df_number <- data.frame(tissue = tissues_name, color = df_color$color[match(tissues_name, df_color$tissue)], 
                          type = df_color$type[match(tissues_name, df_color$tissue)], 
                          ntrain = df_color$nsamples_train[match(tissues_name, df_color$tissue)], stringsAsFactors = F)
  df_number$ntot <- sapply(tissues_name, function(x) nrow(res[res$tissue %in% x,]))
  df_number$nsign <- sapply(tissues_name, function(x) nrow(res[res$tissue %in% x & res[,id_pval_corr] <=pval_FDR,]))
  tmp <- lapply(tissues_name, function(x) res[res$tissue %in% x & res[,id_pval_corr] <=pval_FDR,1])
  tmp_unique <- lapply(1:length(tmp), function(x) setdiff(tmp[[x]], unique(unlist(tmp[-x]))))
  
  df_number$nsign_unique <- sapply(tmp_unique, length)
  df_number$nsign_unique_perc <- df_number$nsign_unique/df_number$nsign
  
  
  df_tot <- data.frame(tissue = tissues_name, n_tot = df_number$ntot, n_sign = df_number$nsign, n_unique = df_number$ntot, n_unique_sign = df_number$nsign)
  df_tot <- rbind(df_tot, data.frame(tissue = 'All', n_tot = nrow(res), n_sign = length(which(res[, id_pval_corr_tot] <= pval_FDR)), 
                                     n_unique = length(unique(res[,1])), n_unique_sign = length(unique(res[res[, id_pval_corr_tot] <= pval_FDR,1]))))
  
  return(list(plot = df_number, table = df_tot))
}

# plot
pl_number_function <- function(df, type_mat, outFold, type_dat){
  
  el <- ifelse(type_mat == 'tscore', 'genes', 'pathways')
  file_name <- sprintf('%s/nsign_el_%s_%s', outFold,  type_mat, type_dat)
  
  df$tissue <- factor(df$tissue, levels = df$tissue)
  df$type <- factor(df$type, levels = unique(df$type))
  
  pl <- ggplot(data = df, aes(x = nsign, y = nsign_unique_perc, size = ntrain, label = tissue, color = type))+
    geom_point(alpha = 0.8)+
    geom_text_repel(box.padding   = 0.35, point.padding = 0.5, segment.color = 'grey50', size = 3, force = 10) +
    xlab(sprintf('n. significant %s', el))+ylab(sprintf('fraction of significant %s\n tissue specific',el))+ 
    theme_bw()+ ggtitle(sprintf('%s %s',type_mat, type_dat))+theme(legend.position = 'right', plot.title = element_text(hjust = 0.5))+
    scale_color_manual(values = unique(df$color))+guides(colour=FALSE, size=guide_legend(title="n. samples\ntraining model"))
  ggsave(filename = paste0(file_name, '.png'), plot = pl, width = 6, height = 5, dpi = 500)
  ggsave(filename = paste0(file_name, '.pdf'), plot = pl, width = 6, height = 5, dpi = 500, compress = F)
  
  pl <- ggplot(data = df, aes(x = nsign, y = ntot, size = ntrain, label = tissue, color = type))+
    geom_point(alpha = 0.8)+
    geom_text_repel(box.padding   = 0.35, point.padding = 0.5, segment.color = 'grey50', size = 3, force = 10) +
    xlab(sprintf('n. significant %s', el))+ylab(sprintf('n. %s',el))+ 
    theme_bw()+ ggtitle(sprintf('%s %s',type_mat, type_dat))+theme(legend.position = 'right', plot.title = element_text(hjust = 0.5))+
    scale_color_manual(values = unique(df$color))+guides(colour=FALSE, size=guide_legend(title="n. samples\ntraining model"))
  file_name <- sprintf('%s/nsign_ngenes_el_%s_%s', outFold,  type_mat, type_dat)
  ggsave(filename = paste0(file_name, '.png'), plot = pl, width = 6, height = 5, dpi = 500)
  ggsave(filename = paste0(file_name, '.pdf'), plot = pl, width = 6, height = 5, dpi = 500, compress = F)
  
}

#### number of associated element per tissue ####

creat_dfnsign_tissueSpec <- function(tissues_name, res, id_pval_corr, pval_FDR, df_color){
  
  res_red <- res[res[,id_pval_corr]<=pval_FDR,]
  unique_el <- unique(res_red[,1])
  ntissues <- c()
  for(i in 1:length(unique_el)){
    ntissues <- c(ntissues, length(res_red$tissue[res_red[, 1] == unique_el[i]]))
  }
  
  df_tissue <- data.frame(n_tissue = 1:length(tissues_name), n_genes = sapply(1:length(tissues_name), function(x) sum(ntissues == x)))
  
  return(df_tissue)
}

# plot
pl_numberSpec_function <- function(df, type_mat, outFold, type_dat){
  
  el <- ifelse(type_mat == 'tscore', 'genes', 'pathways')
  file_name <- sprintf('%s/nsign_pertissue_el_%s_%s', outFold,  type_mat, type_dat)
  
  df$n_tissue <- factor(df$n_tissue, levels = unique(df$n_tissue))
  
  pl <- ggplot(data = df, aes(x = n_tissue, y = n_genes))+
    geom_bar(alpha = 0.7, color = 'black', stat = 'identity', width = 0.7)+
    xlab(sprintf('n. tissues'))+ylab(sprintf('n. significant %s',el))+ 
    theme_bw()+ ggtitle(sprintf('%s %s',type_mat, type_dat))+theme(legend.position = 'right', plot.title = element_text(hjust = 0.5))
  ggsave(filename = paste0(file_name, '.png'), plot = pl, width = 3.5, height = 5, dpi = 500)
  ggsave(filename = paste0(file_name, '.pdf'), plot = pl, width = 3.5, height = 5, dpi = 500, compress = F)
  
  
}

#### create dataframe for manhattan plot ####

create_df_manhattan_plot <- function(tissues_name, res, pval_FDR, df_color, id_pval, id_name, n_sign = NULL, gene = F){
  
  df <- list()
  thr_pval <- vector()
  HLA_reg <- c(26000000/1000000, 36000000/1000000)
  
  if(gene){
    tmp <- lapply(tissues_name, function(x) res[res$tissue == x, c('start_position', 'chrom')])
    df_chr <- data.frame(chr = paste0('chr', 1:22), start = 0, end = 0, pos_start = 0, stringsAsFactors = F)
    df_chr$end <- apply(sapply(tmp, function(x) sapply(1:22, function(y) max(x$start_position[x$chrom == paste0('chr', y)]))), 1, max)/1000000
    df_chr$start <- apply(sapply(tmp, function(x) sapply(1:22, function(y) min(x$start_position[x$chrom == paste0('chr', y)]))), 1, min)/1000000
    a <- cumsum(df_chr$start)
    b <- cumsum(df_chr$end)
    df_chr$pos_start <- c(0, b[-length(b)])
    # df_chr$pos_start <- c(df_chr$start)+c(0, b[-length(b)]), not correct for plot
  }else{
    df_chr <- NULL
  }
  
  for(i in 1:length(tissues_name)){
    
    tmp <- res[res$tissue %in% tissues_name[i], ]
    df[[i]] <- data.frame(tissue = tmp$tissue, pval_tr = -log10(tmp[,id_pval]), zstat = tmp[,id_pval-1], name = tmp[,id_name], stringsAsFactors = F)
    df[[i]]$type <- df_color$type[df_color$tissues == tissues_name[i]]
    
    if(gene){
      df[[i]]$id <- tmp$start_position/1000000
      df[[i]]$chr <- tmp$chr
      tmp_chr <- sapply(1:22, function(x) df[[i]][df[[i]]$chr == paste0('chr',x), 'id'])
      tmp_chr <- mapply(function(x,y) x + y, x = tmp_chr, y = df_chr$pos_start, SIMPLIFY = F)
      df[[i]]$id_pos <- unlist(tmp_chr)
    }
    
    # annotate most significant element
    # df[[i]]$name[-order(df[[i]]$pval_tr, decreasing = T)[1:n_sign]] <- ''
    df[[i]]$sign_name <- 'yes'
    if(!is.null(n_sign)){
      df[[i]]$sign_name[-order(df[[i]]$pval_tr, decreasing = T)[1:n_sign]] <- 'no'  
    }
    # remove more than 1 from HLA region
    id_tmp <- which(df[[i]]$sign_name == 'yes' & df[[i]]$chr == 'chr6' &  df[[i]]$id <= HLA_reg[2] & df[[i]]$id >= HLA_reg[1] )
    if(length(id_tmp)>0){
      df[[i]]$sign_name[id_tmp][-which.max(df[[i]]$pval_tr[id_tmp])] <- 'no'
    }
    
    df[[i]]$sign <- 'no'
    df[[i]]$sign[tmp[,id_pval+2]<=pval_FDR] <- 'yes'
    
    tmp <- tmp[order(tmp[,id_pval], decreasing = T),]
    thr_pval[i] <- tmp[which(tmp[,id_pval+2]<=pval_FDR)[1],id_pval]
    
  }
  
  color_info <- df_color[match(tissues_name, df_color$tissues), ]
  color_info$thr_pval <- thr_pval
  df <- do.call(rbind, df)
  
  if(!gene){
    all_names <- unique(df$name)
    tmp <- data.frame(all_names, 1:length(all_names))
    df$id_pos <- sapply(df$name, function(x) tmp[tmp[,1] == x,2])
  }
  
  
  df$color <- unlist(lapply(1:length(tissues_name), function(x) rep(color_info$color[color_info$tissue == tissues_name[x]],length(which(df$tissue == tissues_name[x])))))
  # name_rep <- names(which(table(df$name)>1))
  # df$color[df$name %in% name_rep] <- '#C0C0C0'
  df$color[df$sign == 'no'] <- '#C0C0C0'
  df$name[df$sign_name == 'no'] <- '' 
  df$name[df$sign_name == 'yes'] <- paste0(df$name[df$sign_name == 'yes'], '\n', df$tissue[df$sign_name == 'yes'])
  
  return(list(df = df, color = color_info, df_chr = df_chr))
  
}

# plot
pl_manhattan_function <- function(data_input, type_mat, outFold, type_dat){
  
  gene <-  type_mat == 'tscore'
  file_name <- sprintf('%smanhattan_%s_%s', outFold,  type_mat, type_dat)
  file_name_z <- sprintf('%smanhattan_zstat_%s_%s', outFold,  type_mat, type_dat)
  
  df <- data_input$df
  info_df <- data_input$color
  
  df$tissue <- factor(df$tissue, levels = data_input$color$tissues)
  df$type <- factor(df$type, levels = unique(data_input$color$type))
  df$color <- factor(df$color, levels = c('#C0C0C0', unique(data_input$color$color)))
  info_df$tissue <-  factor(info_df$tissue, levels = data_input$color$tissues)
  info_df$type <- factor(info_df$type, levels = unique(data_input$color$type))
  info_df$thr_pval_tr = -log10(info_df$thr_pval)
  
  pl <- ggplot(data = df, aes(x = id_pos, y = pval_tr, color = color, label = name))+
    geom_point(alpha = 0.8, size = 0.1)+
    geom_text_repel(segment.color = 'grey50', size = 2, min.segment.length = unit(0, 'lines'),
                    segment.alpha = 0.6,  force = 8) +
    ylab('-log10(pvalue)')+ggtitle(type_dat)+
    theme_bw()+theme(legend.position = 'none', plot.title = element_text(hjust = 0.5))+
    scale_color_manual(values = c('#C0C0C0', unique(data_input$color$color)))
  
  if(gene){
    pl <- pl+scale_x_continuous(breaks=data_input$df_chr$pos_start, labels=data_input$df_chr$chr)+xlab('chromosome')+
      theme(axis.text.x = element_text(angle = 45, vjust = 0.5))
    #geom_hline(yintercept = 8, linetype = 2, size = 0.3)
  }else{
    pl <- pl + xlab('pathways')+ theme(axis.text.x = element_blank(), axis.ticks.x = element_blank())
    #geom_hline(yintercept = 6, linetype = 2, size = 0.3)
  }
  
  ggsave(filename = paste0(file_name, '.pdf'), plot = pl, width = 11, height = 4, dpi = 500)
  ggsave(filename = paste0(file_name, '.png'), plot = pl, width = 11, height = 4, dpi = 500)
  
  # plot z-stat
  pl_z <- ggplot(data = df, aes(x = id_pos, y = zstat, color = color, label = name))+
    geom_point(alpha = 0.8, size = 0.1)+
    geom_text_repel(segment.color = 'grey50', size = 2, min.segment.length = unit(0, 'lines'),
                    segment.alpha = 0.6,  force = 8) +
    ylab('z statistic')+ggtitle(type_dat)+
    theme_bw()+theme(legend.position = 'none', plot.title = element_text(hjust = 0.5))+
    scale_color_manual(values = c('#C0C0C0', unique(data_input$color$color)))
  
  if(gene){
    pl_z <- pl_z+scale_x_continuous(breaks=data_input$df_chr$pos_start, labels=data_input$df_chr$chr)+xlab('chromosome')+
      theme(axis.text.x = element_text(angle = 45, vjust = 0.5))+
      geom_hline(yintercept = 0, linetype = 2, size = 0.3) + 
      ylim(min(df$zstat) - 3, max(df$zstat) + 3)
  }else{
    pl_z <- pl_z + xlab('pathways')+ theme(axis.text.x = element_blank(), axis.ticks.x = element_blank())
    #geom_hline(yintercept = 6, linetype = 2, size = 0.3)
  }
  ggsave(filename = paste0(file_name_z, '.pdf'), plot = pl_z, width = 11, height = 5, dpi = 500)
  ggsave(filename = paste0(file_name_z, '.png'), plot = pl_z, width = 11, height = 5, dpi = 500)
  
}

#### pathway plot with ngenes info ####
create_df_manhattan_plot_path <- function(tissues_name, res, pval_thr, pval_FDR, thr_genes = 0, df_color, id_pval, id_name, gene = F){
  
  df <- list()
  
  for(i in 1:length(tissues_name)){
    
    tmp <- res[res$tissue %in% tissues_name[i], ]
    df[[i]] <- data.frame(tissue = tmp$tissue, pval_tr = -log10(tmp[,id_pval]), name = tmp[,id_name], ngenes = tmp$ngenes_tscore, perc_genes = tmp$ngenes_tscore/tmp$ngenes_path, stringsAsFactors = F)
    df[[i]]$type <- df_color$type[df_color$tissues == tissues_name[i]]
    
    # annotate most significant element
    # df[[i]]$name[-order(df[[i]]$pval_tr, decreasing = T)[1:n_sign]] <- ''
    
    df[[i]]$sign <- 'no'
    df[[i]]$sign[tmp[,id_pval]<=pval_thr] <- 'yes'
    df[[i]]$sign_corr <- 'no'
    df[[i]]$sign_corr[tmp[,id_pval+2]<=pval_FDR] <- 'yes'
  }
  
  color_info <- df_color[match(tissues_name, df_color$tissues), ]
  
  df <- do.call(rbind, df)
  if(thr_genes==0){
    thr_genes <- quantile(df$perc_genes, probs = 0.80)
  }
  print(thr_genes)
  df$both <- 'no'
  df$both[df$perc_genes >=thr_genes & df$sign == 'yes'] <- 'yes'
  # df$name[df$sign_name == 'no'] <- '' 
  df$name[df$both == 'yes'] <- paste0(df$name[df$both == 'yes'], '\n', df$tissue[df$both == 'yes'])
  df$name[df$both == 'no'] <- ''
  
  df$color <- unlist(lapply(1:length(tissues_name), function(x) rep(color_info$color[color_info$tissue == tissues_name[x]],length(which(df$tissue == tissues_name[x])))))
  # name_rep <- names(which(table(df$name)>1))
  # df$color[df$name %in% name_rep] <- '#C0C0C0'
  df$color[df$sign_corr == 'no'] <- '#C0C0C0'
  
  return(list(df = df, color = color_info))
  
}

# plot
pl_manhattan_function_path <- function(data_input, type_mat, outFold, type_dat, pval_thr){
  
  file_name <- sprintf('%s/manhattan_ngenes_vs_pval_%s_%s', outFold,  type_mat, type_dat)
  
  df <- data_input$df
  df$ngenes <- log10(df$ngenes)
  df$perc_genes <- df$perc_genes*100
  info_df <- data_input$color
  val <- min(df$perc[df$name !='']) - 0.01
  
  df$tissue <- factor(df$tissue, levels = data_input$color$tissues)
  df$type <- factor(df$type, levels = unique(data_input$color$type))
  df$color <- factor(df$color, levels = c('#C0C0C0', unique(data_input$color$color)))
  info_df$tissue <-  factor(info_df$tissue, levels = data_input$color$tissues)
  info_df$type <- factor(info_df$type, levels = unique(data_input$color$type))
  
  pl <- ggplot(data = df, aes(x = perc_genes, y = pval_tr, color = color, label = name))+
    geom_point(alpha = 0.9, size = 0.2)+
    geom_hline(yintercept = -log10(pval_thr), linetype = 2, size = 0.3)+
    geom_vline(xintercept = val, linetype = 2, size = 0.3)+
    geom_text_repel(segment.color = 'grey', size = 2, 
                    segment.alpha = 0.5,  force = 15) +
    ylab('-log10(pvalue)')+ggtitle(type_dat)+xlab('percentage of patwhay coverage')+
    theme_bw()+theme(legend.position = 'none', plot.title = element_text(hjust = 0.5))+
    scale_color_manual(values = c('#C0C0C0', unique(data_input$color$color)))
  
  
  ggsave(filename = paste0(file_name, '.pdf'), plot = pl, width = 11, height = 5, dpi = 500, compress = F)
  ggsave(filename = paste0(file_name, '.png'), plot = pl, width = 11, height = 5, dpi = 500)
  
}

# obtain first 5 best results
# best_res_fun <- function(df, tissues, n_top = 5, id_pval){
#   
#   tmp <- lapply(tissues, function(x) df[df$tissue == x,])
#   tmp <- lapply(tmp, function(x) x[order(x[,id_pval])[1:n_top],])
#   res <- do.call(rbind, tmp)
#   res$logpval <- -log10(res[, id_pval])  
#   res$zstat <- res[, id_pval -1]
#   return(res)
# }

best_res_fun <- function(df, tissues, n_top = 3, id_pval, min_ngenes = 3, min_cov = 0.1){
  
  tmp <- lapply(tissues, function(x) df[df$tissue == x,])
  tmp <- lapply(tmp, function(x) x[x$ngenes_tscore>= min_ngenes & x$ngenes_tscore/x$ngenes_path>=min_cov,])
  tmp <- lapply(tmp, function(x) x[order(x[,id_pval])[1:n_top],])
  res <- do.call(rbind, tmp)
  res$logpval <- -log10(res[, id_pval])  
  res$zstat <- res[, id_pval -1]
  return(res)
  
}

plot_best_path <- function(best_res, tissues, color_tissues, title_plot, type_mat, type_dat, width_plot = 5.5, height_plot = 5, outFold, id_name=1, id_pval){
  
  df_tissues <- data.frame(tissue = tissues, id = NA, stringsAsFactors = F)
  tmp <- lapply(df_tissues$tissue, function(x) strsplit(x, split = '[_]')[[1]])
  tmp <- sapply(tmp, function(x) paste0(sapply(x, function(y) substr(y, start = 1, stop = 1)), collapse = ''))
  df_tissues$id <- tmp
  
  if('Brain_Hippocampus' %in% df_tissues$tissue){df_tissues$id[df_tissues$tissue == 'Brain_Hippocampus'] <- 'BHi'}
  if('Brain_Hypothalamus' %in% df_tissues$tissue){df_tissues$id[df_tissues$tissue == 'Brain_Hypothalamus'] <- 'BHy'}
  if('Brain_Cerebellum' %in% df_tissues$tissue){df_tissues$id[df_tissues$tissue == 'Brain_Cerebellum'] <- 'BCe'}
  if('Brain_Cerebellar_Hemisphere' %in% df_tissues$tissue){df_tissues$id[df_tissues$tissue == 'Brain_Cerebellar_Hemisphere'] <- 'BCeH'}
   
 
  best_res <- best_res[best_res$tissue %in% tissues, ]
  color_tmp <- color_tissues[match(tissues, color_tissues$tissues),]
  best_res$tissue <- factor(best_res$tissue, levels = tissues)
  best_res$name <- paste0(best_res$path, ' (', df_tissues$id[match(best_res$tissue, df_tissues$tissue)],')') 
  
  best_res <- best_res[order(best_res[, id_pval]),]
  # id <- which(nchar(best_res$name)>70)
  # if(length(id)>0){
  #   for(i in 1:length(id)){
  #     new <- strsplit(best_res[id[i],id_name], split = '[ ]')[[1]]
  #     split_i <- round(length(new)/2)
  #     new <- paste0(paste0(new[1:split_i], collapse = ' '), '\n', paste0(new[(split_i+1):length(new)], collapse = ' '))
  #     best_res[id[i],id_name] <- new
  #   }
  # }
  best_res$name <- factor(best_res$name, levels = rev(unique(best_res$name)))
  best_res$impr <- 'plain'
  best_res$impr[best_res$improvement_sign] <- 'bold'
  best_res$add_info <- mapply(function(x, y, z) paste0(x, '/', y, ' (', z,')'), x = best_res$ngenes_tscore, y = best_res$ngenes_path, z = round(best_res$ngenes_tscore/best_res$ngenes_path, digits = 2))
  # plot significnace
  path_pval_pl <- ggplot(best_res, aes(x=name, y=logpval, fill = tissue))+
    geom_bar(stat = "identity", color = 'black', alpha = 0.6)+theme_bw()+ylab('-log10(pvalue)')+ xlab("")+
    ggtitle(title_plot)+
    scale_fill_manual(values = color_tmp$color)+
    geom_text(aes(label = add_info), position = position_stack(0.5), color = "black", size = 3)+
    # geom_text_repel(aes(pathway, sign + 0.1, label = pathway), size = 2)+ 
    theme(plot.title = element_text(hjust = 1), strip.text = element_text(size = 8),text = element_text(size = 11), 
          axis.text.y=element_text(size = 7.5, angle = 0, hjust = 1, face = rev(best_res$impr)), legend.position = 'right') + coord_flip()
  
  file_name <- sprintf('%s/barplot_%s_best_%s', outFold,  type_mat, type_dat)
  ggsave(filename = paste0(file_name, '.pdf'), plot = path_pval_pl, width = width_plot, height = height_plot, dpi = 500, compress = F)
  ggsave(filename = paste0(file_name, '.png'), plot = path_pval_pl, width = width_plot, height = height_plot, dpi = 500)
  
  # plot significnace
  path_zstat_pl <- ggplot(best_res, aes(x=name, y=zstat, fill = tissue))+
    geom_bar(stat = "identity", color = 'black', alpha = 0.6)+theme_bw()+ylab('z-statistic')+ xlab("")+
    ggtitle(title_plot)+
    scale_fill_manual(values = color_tmp$color)+
    geom_text(aes(label = add_info), position = position_stack(0.5), color = "black", size = 3)+
    # geom_text_repel(aes(pathway, sign + 0.1, label = pathway), size = 2)+ 
    theme(plot.title = element_text(hjust = 1), strip.text = element_text(size = 8),text = element_text(size = 11), 
          axis.text.y=element_text(size = 7.5, angle = 0, hjust = 1, face = rev(best_res$impr)), legend.position = 'right') + coord_flip()
  
  file_name <- sprintf('%s/barplot_zstat_%s_best_%s', outFold,  type_mat, type_dat)
  ggsave(filename = paste0(file_name, '.pdf'), plot = path_zstat_pl, width = width_plot, height = height_plot, dpi = 500, compress = F)
  ggsave(filename = paste0(file_name, '.png'), plot = path_zstat_pl, width = width_plot, height = height_plot, dpi = 500)
  
}

# plot_best_path <- function(best_res, tissues, color_tissues, title_plot, type_mat, type_dat, width_plot = 5.5, height_plot = 5, outFold, id_name=1){
#   best_res <- best_res[best_res$tissue %in% tissues, ]
#   color_tmp <- color_tissues[match(tissues, color_tissues$tissues),]
#   best_res$tissue <- factor(best_res$tissue, levels = tissues)
#   
#   id <- which(nchar(best_res[,id_name])>70)
#   if(length(id)>0){
#     for(i in 1:length(id)){
#       new <- strsplit(best_res[id[i],id_name], split = '[ ]')[[1]]
#       split_i <- round(length(new)/2)
#       new <- paste0(paste0(new[1:split_i], collapse = ' '), '\n', paste0(new[(split_i+1):length(new)], collapse = ' '))
#       best_res[id[i],id_name] <- new
#     }
#   }
#   best_res[,id_name] <- factor(best_res[,id_name], levels = rev(unique(best_res[,id_name])))
#   
#   # plot significnace
#   path_pval_pl <- ggplot(best_res, aes(x=path, y=logpval, fill = tissue))+
#     geom_bar(stat = "identity", color = 'black', alpha = 0.8)+theme_bw()+ylab('-log10(pvalue)')+ xlab("")+
#     ggtitle(title_plot)+
#     facet_wrap(.~tissue, ncol = 1, scales = 'free_y')+
#     scale_fill_manual(values = color_tmp$color)+
#     # geom_text_repel(aes(pathway, sign + 0.1, label = pathway), size = 2)+ 
#     theme(plot.title = element_text(hjust = 1), strip.text = element_text(size = 8),text = element_text(size = 11), 
#           axis.text.y=element_text(size = 7.5, angle = 0, hjust = 1), legend.position = 'none') + coord_flip()
#   
#   file_name <- sprintf('%s/barplot_%s_best_%s', outFold,  type_mat, type_dat)
#   ggsave(filename = paste0(file_name, '.pdf'), plot = path_pval_pl, width = width_plot, height = height_plot, dpi = 500, compress = F)
#   ggsave(filename = paste0(file_name, '.png'), plot = path_pval_pl, width = width_plot, height = height_plot, dpi = 500)
#   
#   # plot significnace
#   path_zstat_pl <- ggplot(best_res, aes(x=path, y=zstat, fill = tissue))+
#     geom_bar(stat = "identity", color = 'black', alpha = 0.8)+theme_bw()+ylab('-log10(pvalue)')+ xlab("")+
#     ggtitle(title_plot)+
#     facet_wrap(.~tissue, ncol = 1, scales = 'free_y')+
#     scale_fill_manual(values = color_tmp$color)+
#     # geom_text_repel(aes(pathway, sign + 0.1, label = pathway), size = 2)+ 
#     theme(plot.title = element_text(hjust = 1), strip.text = element_text(size = 8),text = element_text(size = 11), 
#           axis.text.y=element_text(size = 7.5, angle = 0, hjust = 1), legend.position = 'none') + coord_flip()
#   
#   file_name <- sprintf('%s/barplot_zstat_%s_best_%s', outFold,  type_mat, type_dat)
#   ggsave(filename = paste0(file_name, '.pdf'), plot = path_zstat_pl, width = width_plot, height = height_plot, dpi = 500, compress = F)
#   ggsave(filename = paste0(file_name, '.png'), plot = path_zstat_pl, width = width_plot, height = height_plot, dpi = 500)
# }
# 

# Venn diagram
venn_plot <- function(gwas_known_file, tscore, pval_FDR, type_dat, type_mat){
  
  file_name <-  sprintf('%sVennDiag_list_genes_%s_%s', fold, type_mat, type_dat)
  
  gwas_summary <- read.table(gwas_known_file, h=T, stringsAsFactors = F, sep = '\t')
  genes_tot <- unique(unlist(lapply(gwas_summary$gene_int, function(x) strsplit(x, split = ',')[[1]])))
  
  genes_new <- unique(tscore$external_gene_name[tscore[, 10] <= pval_FDR])
  
  x <- list('Genes in GWAS loci' = genes_tot , 'PriLer Genes' = genes_new)
  png(sprintf('%s.png',file_name), units = 'in',width = 4, height = 4, res = 500)
  v0 <- venn.diagram( x, filename=NULL,fill = c("grey60", "cornflowerblue"),
                      alpha = c(0.5, 0.5), cat.cex = 1.3, cex=1.3, cat.pos = 0)
  overlaps <- calculate.overlap(x)
  # extract indexes of overlaps from list names
  indx <- as.numeric(substr(names(overlaps),2,2))
  # v0[[7]]$label  <- paste(overlaps[[3]], collapse="\n")  
  grid.newpage()
  grid.draw(v0)
  dev.off()
  
  pdf(sprintf('%s.pdf',file_name),width = 4, height = 4)
  v0 <- venn.diagram( x, filename=NULL,fill = c("grey60", "cornflowerblue"),
                      alpha = c(0.5, 0.5), cat.cex = 1.3, cex=1.3, cat.pos = 0)
  overlaps <- calculate.overlap(x)
  # extract indexes of overlaps from list names
  indx <- as.numeric(substr(names(overlaps),2,2))
  # v0[[7]]$label  <- paste(overlaps[[3]], collapse="\n")  
  grid.newpage()
  grid.draw(v0)
  dev.off()
  
}

# # # results from latest GWAS
# latest_res <- read.xls('/psycl/g/mpsziller/lucia/refData/41588_2017_BFng3913_MOESM2_ESM.xlsx', h=T, sheet = 4, skip=1, stringsAsFactors = F)
# latest_res <- latest_res[-nrow(latest_res), ]
# genes_list <- latest_res[!duplicated(latest_res[,8]),]
# genes_list$chr <- sapply(genes_list[,3], function(x) strsplit(x, split = ':')[[1]][1])
# genes_list$pos <- sapply(genes_list[,3], function(x) strsplit(x, split = ':')[[1]][2])
# 
# tmp <- tscore[tscore[, 10] <= pval_FDR,]
# tmp <- tmp[!duplicated(tmp[,1]),]
# tmp$chr_new <- sapply(tmp$chr, function(x) as.numeric(strsplit(x, split = 'chr')[[1]][2]))
# tmp <- tmp[order(tmp$start_position),]
# tmp <- tmp[order(tmp$chr_new),]
# tmp$new_loci <- F
# # find new loci
# for(i in 1:nrow(tmp)){
#   new <- latest_res[tmp$chr_new[i] - as.numeric(latest_res$chr) == 0 & abs(tmp$TSS_start[i] - as.numeric(latest_res$pos))<500000 ,]
#   if(nrow(new) == 0){tmp$new_loci[i] <- T}
# }

