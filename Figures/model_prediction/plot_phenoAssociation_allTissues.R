library(ggplot2)
library(ggpubr)
library(qvalue)
library(corrplot)
library(RColorBrewer)
library(ggrepel)
library(argparse)
library(Matrix)


Sys.setlocale("LC_NUMERIC", "C")
options(bitmapType = 'cairo', device = 'png')

parser <- ArgumentParser(description="plot pheno association")
parser$add_argument("--fold", type = "character", help = "fold output")
parser$add_argument("--fold_tissue", type = "character", help = "fold input case example enrichment")
parser$add_argument("--fold_geno_input", type = "character", default = 'NA',help = "fold genotype input info (for case pathway)")
parser$add_argument("--train_fold", type = "character", help = "train model fold")
parser$add_argument("--color_file", type = "character", help = "file with tissues color code")
parser$add_argument("--pheno", type = "character", help = "name phenotype")
parser$add_argument("--type_dat", type = "character", help = "name to append to plot file")
parser$add_argument("--pval_FDR", type = "double", default = 0.05, help = "pvalue threshold")


args <- parser$parse_args()
fold <- args$fold
fold_tissue <- args$fold_tissue
train_fold <- args$train_fold
color_file <- args$color_file
pheno <- args$pheno
type_dat <- args$type_dat
pval_FDR <- args$pval_FDR
fold_geno_input <- args$fold_geno_input

# fold <- '/psycl/g/mpsziller/lucia/CAD_UKBB/eQTL_PROJECT/OUTPUT_GTEx/predict_CAD/AllTissues/200kb/CAD_GWAS_bin5e-2/UKBB/'
# fold_tissue <- '/psycl/g/mpsziller/lucia/CAD_UKBB/eQTL_PROJECT/OUTPUT_GTEx/predict_CAD/Artery_Aorta/200kb/CAD_GWAS_bin5e-2/UKBB/devgeno0.01_testdevgeno0/'
# color_tissues <- read.table('/psycl/g/mpsziller/lucia/color_tissues.txt', h=T, stringsAsFactors = F)
# pval_FDR <- 0.05
# train_fold <- '/psycl/g/mpsziller/lucia/CAD_UKBB/eQTL_PROJECT/OUTPUT_GTEx/train_GTEx/'
# pheno <- 'CAD_HARD'
# type_dat <- 'CAD_HARD-UKBB'
# fold_geno_input <- '/psycl/g/mpsziller/lucia/CAD_UKBB/eQTL_PROJECT/INPUT_DATA_GTEx/GTEX_v6/Genotyping_data/Genotype_VariantsInfo_GTEx-PGCgwas-CADgwas-CADall-UKBB_'

color_tissues <- read.table(color_file, h=T, stringsAsFactors = F)

# load results
tscore <- read.delim(sprintf('%stscore_pval_%s_covCorr.txt', fold, pheno), h=T, stringsAsFactors = F, sep = '\t')
pathR <- read.delim(sprintf('%spath_Reactome_pval_%s_covCorr_filt.txt', fold, pheno), h=T, stringsAsFactors = F, sep = '\t')
pathGO <- read.delim(sprintf('%spath_GO_pval_%s_covCorr_filt.txt', fold, pheno), h=T, stringsAsFactors = F, sep = '\t')

tissues <- unique(tscore$tissue)
if(grepl('CAD', pheno)){
  train_fold <- paste0(train_fold, tissues, '/200kb/CAD_GWAS_bin5e-2/')  
}else{
  train_fold <- paste0(train_fold, tissues, '/')
}


# gene location
tscore$start_position <- NA
tscore$chrom <- NA

for(i in 1:length(train_fold)){
  
  tmp <- read.table(sprintf('%s/resPrior_regEval_allchr.txt', train_fold[i]), h=T,stringsAsFactors = F)
  tmp <- tmp[match(tscore$ensembl_gene_id[tscore$tissue == tissues[i]], tmp$ensembl_gene_id),]
  tscore$start_position[tscore$tissue == tissues[i]] <- tmp$start_position
  tscore$chrom[tscore$tissue == tissues[i]] <- tmp$chrom
  
}

#####################################################################################################
### functions ###

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
pl_corr <- function(res_cor, type_mat, type_dat, tissues_name, df_color, outFold){
  
  # correlation
  diag(res_cor$cor) <- NA
  col <- colorRampPalette(brewer.pal(9, 'Oranges'))(100)
  
  ord <- corrMatOrder(res_cor$cor, order="hclust", hclust.method = 'ward.D')
  newcolours <- df_color$color[match(tissues_name, df_color$tissues)][ord]
  title_pl <- sprintf('%s %s Spearman correlation z stat', type_mat, type_dat)
  
  pdf(file = sprintf('%s/corr_zscore_%s_%s.pdf', outFold, type_mat, type_dat), width = 11, height = 8, compress = F)
  corrplot(res_cor$cor, type="upper", order = 'hclust', hclust.method = 'ward.D',
           tl.col = newcolours, title = title_pl,
           col = c(col), method = 'color', tl.srt=45, cl.align.text='c',
           addCoef.col = "black",na.label = 'square', na.label.col = 'darkgrey', is.corr = F, number.cex=0.7, mar = c(0,0,5,0))
  dev.off()
  
  png(file = sprintf('%s/corr_zscore_%s_%s.png', outFold,  type_mat, type_dat), units = 'in', width = 11, height = 8, res = 300)
  corrplot(res_cor$cor, type="upper", order = 'hclust', hclust.method = 'ward.D',
           tl.col = newcolours, title = title_pl,
           col = c(col), method = 'color', tl.srt=45, cl.align.text='c',
           addCoef.col = "black",na.label = 'square', na.label.col = 'darkgrey', is.corr = F, number.cex=0.7, mar = c(0,0,5,0))
  dev.off()
  
  # intersection
  diag(res_cor$perc) <- NA
  col <- colorRampPalette(brewer.pal(9, 'Greens'))(100)
  
  title_pl <- sprintf('%s %s percentage common', type_mat, type_dat)
  
  pdf(file = sprintf('%s/perc_zscore_%s_%s.pdf', outFold, type_mat, type_dat), width = 11, height = 8, compress = F)
  corrplot(res_cor$perc[ord,ord], type="lower", 
           tl.col = newcolours, title = title_pl,
           col = c(col), method = 'color', tl.srt=45, cl.align.text='c',
           addCoef.col = "black",na.label = 'square', na.label.col = 'darkgrey', is.corr = F, number.cex=0.7, mar = c(0,0,5,0))
  dev.off()
  
  png(file = sprintf('%s/perc_zscore_%s_%s.png', outFold,  type_mat, type_dat), units = 'in', width = 11, height = 8, res = 300)
  corrplot(res_cor$perc[ord,ord], type="lower", 
           tl.col = newcolours, title = title_pl,
           col = c(col), method = 'color', tl.srt=45, cl.align.text='c',
           addCoef.col = "black",na.label = 'square', na.label.col = 'darkgrey', is.corr = F, number.cex=0.7, mar = c(0,0,5,0))
  dev.off()
  
}

#### number of associated element ####

creat_dfnsign <- function(tissues_name, res, id_pval_corr, pval_FDR, df_color){
  
  df_number <- data.frame(tissue = tissues_name, color = df_color$color[match(tissues_name, df_color$tissue)], 
                          type = df_color$type[match(tissues_name, df_color$tissue)], 
                          ntrain = df_color$nsamples_train[match(tissues_name, df_color$tissue)], stringsAsFactors = F)
  df_number$ntot <- sapply(tissues_name, function(x) nrow(res[res$tissue %in% x,]))
  df_number$nsign <- sapply(tissues_name, function(x) nrow(res[res$tissue %in% x & res[,id_pval_corr] <=pval_FDR,]))
  tmp <- lapply(tissues_name, function(x) res[res$tissue %in% x & res[,id_pval_corr] <=pval_FDR,1])
  tmp_unique <- lapply(1:length(tmp), function(x) setdiff(tmp[[x]], unique(unlist(tmp[-x]))))
  
  df_number$nsign_unique <- sapply(tmp_unique, length)
  df_number$nsign_unique_perc <- df_number$nsign_unique/df_number$nsign
  
  return(df_number)
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
  
}

#### create dataframe for manhattan plot ####

create_df_manhattan_plot <- function(tissues_name, res, pval_FDR, df_color, id_pval, id_name, n_sign, gene = F){
  
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
    df_chr$pos_start <- c(df_chr$start)+c(0, b[-length(b)])
  }else{
    df_chr <- NULL
  }
  
  for(i in 1:length(tissues_name)){
    
    tmp <- res[res$tissue %in% tissues_name[i], ]
    df[[i]] <- data.frame(tissue = tmp$tissue, pval_tr = -log10(tmp[,id_pval]), name = tmp[,id_name], stringsAsFactors = F)
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
    df[[i]]$sign_name[-order(df[[i]]$pval_tr, decreasing = T)[1:n_sign]] <- 'no'
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
  file_name <- sprintf('%s/manhattan_%s_%s', outFold,  type_mat, type_dat)
  
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
    geom_text_repel(segment.color = 'grey50', size = 2, 
                    segment.alpha = 0.6,  force = 15) +
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
  
  ggsave(filename = paste0(file_name, '.pdf'), plot = pl, width = 11, height = 4, dpi = 500, compress = F)
  ggsave(filename = paste0(file_name, '.png'), plot = pl, width = 11, height = 4, dpi = 500)
  
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
best_res_fun <- function(df, tissues, n_top = 5, id_pval){
  
  tmp <- lapply(tissues, function(x) df[df$tissue == x,])
  tmp <- lapply(tmp, function(x) x[order(x[,id_pval])[1:n_top],])
  res <- do.call(rbind, tmp)
  res$logpval <- -log10(res[, id_pval])  
  return(res)
  
}

plot_best_path <- function(best_res, tissues, color_tissues, title_plot, type_mat, type_dat, width_plot = 5.5, height_plot = 5, outFold){
  best_res <- best_res[best_res$tissue %in% tissues, ]
  color_tmp <- color_tissues[match(tissues, color_tissues$tissues),]
  best_res$tissue <- factor(best_res$tissue, levels = tissues)
  
  id <- which(nchar(best_res[,1])>70)
  if(length(id)>0){
    for(i in 1:length(id)){
      new <- strsplit(best_res[id[i],1], split = '[ ]')[[1]]
      split_i <- round(length(new)/2)
      new <- paste0(paste0(new[1:split_i], collapse = ' '), '\n', paste0(new[(split_i+1):length(new)], collapse = ' '))
      best_res[id[i],1] <- new
    }
  }
  best_res[,1] <- factor(best_res[,1], levels = rev(unique(best_res[,1])))
  
  # plot significnace
  path_pval_pl <- ggplot(best_res, aes(x=path, y=logpval, fill = tissue))+
    geom_bar(stat = "identity", color = 'black', alpha = 0.8)+theme_bw()+ylab('-log10(pvalue)')+ xlab("")+
    ggtitle(title_plot)+
    facet_wrap(.~tissue, ncol = 1, scales = 'free_y')+
    scale_fill_manual(values = color_tmp$color)+
    # geom_text_repel(aes(pathway, sign + 0.1, label = pathway), size = 2)+ 
    theme(plot.title = element_text(hjust = 1), strip.text = element_text(size = 8),text = element_text(size = 11), 
          axis.text.y=element_text(size = 7.5, angle = 0, hjust = 1), legend.position = 'none') + coord_flip()
  
  file_name <- sprintf('%s/barplot_%s_best_%s', outFold,  type_mat, type_dat)
  ggsave(filename = paste0(file_name, '.pdf'), plot = path_pval_pl, width = width_plot, height = height_plot, dpi = 500, compress = F)
  ggsave(filename = paste0(file_name, '.png'), plot = path_pval_pl, width = width_plot, height = height_plot, dpi = 500)
  
  
  
}


###################################################################################################################################################

### correlation tissues ###
tscore_cor <- create_cor(tissues_name = tissues, res = tscore, id_z = 7)
pathR_cor <- create_cor(tissues_name = tissues, res = pathR, id_z = 12)
pathGO_cor <- create_cor(tissues_name = tissues, res = pathGO, id_z = 14)

pl_corr(tscore_cor, type_mat = 'tscore', type_dat = type_dat, tissues_name = tissues, df_color = color_tissues, outFold = fold)
pl_corr(pathR_cor, type_mat = 'path_Reactome', type_dat = type_dat, tissues_name = tissues, df_color = color_tissues, outFold = fold)
pl_corr(pathGO_cor, type_mat = 'path_GO', type_dat = type_dat, tissues_name = tissues, df_color = color_tissues, outFold = fold)

### number of associated elements ###
tscore_nsgin <- creat_dfnsign(tissues_name = tissues, res = tscore, id_pval_corr = 10, pval_FDR = pval_FDR, df_color = color_tissues)
pathR_nsgin <- creat_dfnsign(tissues_name = tissues, res = pathR, id_pval_corr = 15, pval_FDR = pval_FDR, df_color = color_tissues)
pathGO_nsgin <- creat_dfnsign(tissues_name = tissues, res = pathGO, id_pval_corr = 17, pval_FDR = pval_FDR, df_color = color_tissues)

pl_number_function(df = tscore_nsgin, type_mat = 'tscore', outFold = fold, type_dat = type_dat)
pl_number_function(df = pathR_nsgin, type_mat = 'path_Reactome', outFold = fold, type_dat = type_dat)
pl_number_function(df = pathGO_nsgin, type_mat = 'path_GO', outFold = fold, type_dat = type_dat)

### manhattan plot ###
if(grepl('SCZ', pheno)){n_sign=10}
if(grepl('CAD', pheno)){n_sign=3}

tscore_df <- create_df_manhattan_plot(tissues_name = tissues, res = tscore, id_pval = 8, pval_FDR = pval_FDR, df_color = color_tissues, id_name = 2, n_sign = n_sign, gene = T)
pathR_df <- create_df_manhattan_plot(tissues_name = tissues, res = pathR, id_pval = 13, pval_FDR = pval_FDR, df_color = color_tissues, id_name = 1, n_sign = 2)
pathGO_df <- create_df_manhattan_plot(tissues_name = tissues, res = pathGO, id_pval = 15, pval_FDR = pval_FDR, df_color = color_tissues, id_name = 2, n_sign = 2)

pl_manhattan_function(data_input = tscore_df, type_mat = 'tscore', outFold = fold, type_dat = type_dat)
pl_manhattan_function(data_input = pathR_df, type_mat = 'path_Reactome', outFold = fold, type_dat = type_dat)
pl_manhattan_function(data_input = pathGO_df, type_mat = 'path_GO', outFold = fold, type_dat = type_dat)

### pathway plot ngenes info ###
pathR_df <- create_df_manhattan_plot_path(tissues_name = tissues,  res = pathR, id_pval = 13, thr_genes = 0.1, pval_thr = 10^-4, pval_FDR = pval_FDR, df_color = color_tissues, id_name = 1)
pathGO_df <- create_df_manhattan_plot_path(tissues_name = tissues, res = pathGO, id_pval = 15, pval_thr = 10^-4, pval_FDR = pval_FDR,df_color = color_tissues, id_name = 2)

pl_manhattan_function_path(data_input = pathR_df, type_mat = 'path_Reactome', outFold = fold, type_dat = type_dat, pval_thr = 10^-4)
pl_manhattan_function_path(data_input = pathGO_df, type_mat = 'path_GO', outFold = fold, type_dat = type_dat, pval_thr = 10^-4)

###################################################################################################################################################

# specific to phenotype
if(grepl('CAD', pheno)){
  
  best_pathR <- best_res_fun(pathR, tissues, id_pval = 13)
  best_pathGO <- best_res_fun(pathGO, tissues, id_pval = 15)
  
  plot_best_path(best_res = best_pathR, color_tissues = color_tissues, title_plot = sprintf('Reactome pathways %s', pheno), type_mat = 'path_Reactome', outFold = fold, type_dat = type_dat, 
                 tissues = c('Adipose_Visceral_Omentum', 'Artery_Aorta', 'Colon_Sigmoid', 'Heart_Left_Ventricle', 'Liver'), height_plot = 6.5)
  plot_best_path(best_res = best_pathR, color_tissues = color_tissues, title_plot = sprintf('Reactome pathways %s', pheno), type_mat = 'path_Reactome', outFold = fold, type_dat = paste0(type_dat, 'v2'),
                 tissues = c('Adipose_Subcutaneous', 'Artery_Coronary', 'Colon_Transverse', 'Heart_Atrial_Appendage', 'Adrenal_Gland'), height_plot = 6.5)

  
  plot_best_path(best_res = best_pathGO, color_tissues = color_tissues, title_plot = sprintf('GO pathways %s', pheno),  type_mat = 'path_GO', outFold = fold, type_dat = type_dat, 
                 tissues = c('Adipose_Visceral_Omentum', 'Artery_Aorta', 'Colon_Sigmoid', 'Heart_Left_Ventricle', 'Liver'), height_plot = 6.5)
  plot_best_path(best_res = best_pathGO, color_tissues = color_tissues, title_plot = sprintf('GO pathways %s', pheno),  type_mat = 'path_GO', outFold = fold, type_dat =  paste0(type_dat, 'v2'),
                 tissues = c('Adipose_Subcutaneous', 'Artery_Coronary', 'Colon_Transverse', 'Heart_Atrial_Appendage', 'Adrenal_Gland'), height_plot = 6.5)

  
}

if(grepl('SCZ', pheno)){

  best_pathR <- best_res_fun(pathR, tissues, id_pval = 13)
  best_pathGO <- best_res_fun(pathGO, tissues, id_pval = 15)
  plot_best_path(best_res = best_pathR, color_tissues = color_tissues, title_plot = sprintf('Reactome pathways %s', pheno), type_mat = 'path_Reactome',outFold = fold, type_dat = type_dat,
  tissues = c('DLPC_CMC', 'Brain_Cerebellum', 'Brain_Hypothalamus', 'Cells_EBV-transformed_lymphocytes', 'Colon_Transverse'), height_plot = 6.5)
  plot_best_path(best_res = best_pathGO, color_tissues = color_tissues, title_plot = sprintf('GO pathways %s', pheno),  type_mat = 'path_GO', outFold = fold, type_dat = type_dat,
  tissues = c('DLPC_CMC', 'Brain_Cerebellum', 'Brain_Hypothalamus', 'Cells_EBV-transformed_lymphocytes', 'Colon_Transverse'), height_plot = 6.5)

}

### example pathway enrichment
if(grepl('CAD', pheno)){
  
  if(pheno == 'CAD_HARD'){id_pval <- 1}
  if(pheno == 'CAD_SOFT'){id_pval <- 2}
  
  tissue <- 'Artery_Aorta'
  pathway <- 'Death Receptor Signalling'
  new_path <- paste0(strsplit(pathway, split = '[ ]')[[1]], collapse = '_')
  color_tmp <- color_tissues$color[color_tissues$tissues == tissue]
  
  # load
  info_res <- get(load(sprintf('%spval_CAD_pheno_covCorr.RData', fold_tissue)))
  id <- which(info_res$pathScore_reactome[[id_pval]]$path == pathway)
  genes_path <- info_res$info_pathScore_reactome[[id_pval]][[id]]
  gene_res <- info_res$tscore[[id_pval]]
  gene_info <- read.table(sprintf('%s/resPrior_regEval_allchr.txt', train_fold[grepl(tissue, train_fold)]), h=T,stringsAsFactors = F)

  id <- sapply(gene_res$external_gene_name, function(x) which(gene_info$external_gene_name == x))
  if(any(sapply(id, length)>1)){
    rm_id <- names(which(sapply(id, length)>1))
    gene_res <- gene_res[! gene_res$external_gene_name %in% rm_id, ]
    id <- id[-which(sapply(id, length)>1)]
    id <- unlist(id)
  }
  gene_info <- gene_info[id, ] 
  identical(gene_info$external_gene_name, gene_res$external_gene_name)
  
  id <- which(colnames(gene_info) == 'train_dev')-1
  gene_res <- cbind(gene_res, gene_info[,1:id])
  
  gene_tmp <- lapply(paste0('chr', 1:22), function(x) gene_res[gene_res$chrom==x,])
  start_pos <- sapply(gene_tmp, function(x) min(x$start_position))
  end_pos <- sapply(gene_tmp, function(x) max(x$end_position))
  df_add <- data.frame(start = start_pos, end = end_pos)
  df_add$start_plot <- cumsum(c(0,df_add$end[-nrow(df_add)]))+df_add$start
  df_add$add_plot <- cumsum(c(0,df_add$end[-nrow(df_add)]))
  new_pos <- mapply(function(x, y) x + y$start_position, x = df_add$add_plot, y = gene_tmp, SIMPLIFY = T)
  new_pos <- unlist(new_pos)
  
  gene_tmp <- do.call(rbind,gene_tmp)
  gene_tmp$new_pos <- new_pos
  
  new_df <- data.frame(pos = new_pos, val = -log10(gene_tmp[,8]), name = gene_tmp$external_gene_name, path = 0, stringsAsFactors = F, 
                       chr = sapply(gene_tmp$chrom, function(x) strsplit(x, split = 'chr')[[1]][2]))
  new_df$chr <- as.numeric(new_df$chr)
  new_df$path[new_df$name %in% genes_path$tscore$external_gene_name] <- 1 
  new_df$path[!(new_df$name %in% genes_path$tscore$external_gene_name) & (new_df$chr %% 2 ==0)] <- 2
  new_df$path <- factor(new_df$path, levels = c(1,0,2))
  new_df$type <- 0
  new_df$type[new_df$name%in% genes_path$tscore$external_gene_name] <- 1 
  new_df$type <- factor(new_df$type, levels = c(1,0))
  
  int_val = -log10(genes_path$path[,13])
  pl_manh <- ggplot(new_df, aes(x = pos, y = val, color = path, alpha = type)) + 
    geom_abline(intercept = int_val, slope = 0, color=color_tmp, linetype="dashed")+
    geom_point() + ggtitle('Tscore association with CAD')+
    scale_color_manual(values = c(color_tmp, "grey90", 'grey70')) +
    ylim(0, max(c(int_val, min(c(max(new_df$val), 8)))))+
    geom_text(x=new_df$pos[nrow(new_df)-1000], y=int_val+0.5, label=genes_path$path$path, size = 5, color = color_tmp)+
    geom_text_repel(
      data = subset(new_df, path == 1), aes(label = new_df$name[new_df$path==1]), size = 3.9, color = 'black', box.padding = unit(0.35, "lines"), point.padding = unit(0.3, "lines"))+
    xlab('chromosome') +  ylab('-log10(pvalue)')+ theme_bw() +
    theme(plot.title = element_text(hjust = 0.5), text = element_text(size = 15),
          axis.text.x=element_text(size = 11, angle = 45),
          axis.text.y=element_text(size = 11), legend.position = 'none')+
    scale_size_discrete(range = c(1.2, 0.3, 0.3))+
    scale_alpha_discrete(range = c(1, 0.5))+
    scale_x_continuous(breaks = df_add$start_plot, labels = c(1:22))
  ggsave(filename = sprintf('%smanhattanPlot_tscore_genes_%s_%s_%s.png', fold, pheno, tissue, new_path), plot = pl_manh, width = 10, height = 5, dpi=500, device = 'png')
  ggsave(filename = sprintf('%smanhattanPlot_tscore_genes_%s_%s_%s.pdf', fold, pheno, tissue, new_path), plot = pl_manh, width = 10, height = 5, device = 'pdf')
  
  ### plot SNPs gwas info for genes in the pathway ###
  df_genePath <- gene_tmp[gene_tmp$external_gene_name %in%  new_df$name[new_df$path == 1],]
  resBeta <- get(load(sprintf('%sresPrior_regCoeffSnps_allchr.RData',  train_fold[grepl(tissue, train_fold)])))
  # load only certain chr
  beta_res <- NULL
  id_chr <- unique(df_genePath$chrom)
  for(i in id_chr){
    
    print(i)
    ind_chr <- as.numeric(strsplit(i, split = 'chr')[[1]][2])
    gene_pos <- read.table(sprintf('%shg19_ENSEMBL_TSS_%s_matched.txt', train_fold[grepl(tissue, train_fold)], i), header = T, stringsAsFactors = F)
    snp_pos <- read.table(sprintf('%s%s.txt', fold_geno_input, i), header = T, stringsAsFactors = F)
    id <- which(gene_pos$ensembl_gene_id %in% df_genePath$ensembl_gene_id)
    tmp <- resBeta[[ind_chr]][,id]
    if(length(id)==1){
      beta_res <- rbind(beta_res, cbind(data.frame(stringsAsFactors = F, val = tmp[tmp!=0], id = which(tmp!=0), 
                                                   gene = gene_pos$external_gene_name[id]), snp_pos[which(tmp!=0), ]))  
    }else{
      for(j in id){
        tmp <- resBeta[[ind_chr]][,j]
        beta_res <- rbind(beta_res, cbind(data.frame(stringsAsFactors = F,val = tmp[tmp!=0], id = which(tmp!=0), 
                                                     gene = gene_pos$external_gene_name[j]), snp_pos[which(tmp!=0), ]))    
      }
    }
    
  }
  
  # find double of beta_res and merge
  dup_snp <- names(which(table(beta_res$ID)>1))
  if(length(dup_snp)>1){
    for(i in 1:length(dup_snp)){
      id <- which(beta_res$ID == dup_snp[i])
      beta_res <- rbind(beta_res, beta_res[id[1], ])
      beta_res$val[nrow(beta_res)] <- mean(beta_res$val[id])
      beta_res$gene[nrow(beta_res)] <- paste(beta_res$gene[id], sep = '', collapse = '_')
      beta_res <- beta_res[-id, ]
    }
  }
  
  # plot only gwas results: SNPs that influece have not so much relevance
  df_start_end <- data.frame(chr = 1:22, start = 0, end = 0)
  for(i in 1:22){
    print(i)
    snp_pos <- read.table(sprintf('%shg19_SNPs_chr%i_matched.txt', train_fold[grepl(tissue, train_fold)], i), header = T, stringsAsFactors = F)
    df_start_end[i,2:3] <- c(snp_pos$position[1], snp_pos$position[nrow(snp_pos)]) 
  }
  
  beta_res <- lapply(1:22, function(x) beta_res[beta_res$CHR == x, ])
  df_start_end$add_plot <- cumsum(c(0,df_start_end$end[-nrow(df_start_end)]))
  new_pos <- mapply(function(x, y) x + y$POS, x = df_start_end$add_plot, y = beta_res, SIMPLIFY = T)
  new_pos <- unlist(new_pos)
  beta_res <- do.call(rbind, beta_res)
  beta_res$new_pos <- new_pos
  beta_res$transf_pvalue <- -log10(beta_res$CAD_p_dgc)
  df_start_end$start_plot <- cumsum(c(0,df_start_end$end[-nrow(df_start_end)]))+df_start_end$start
  
  pl_manh_snps <- ggplot(beta_res, aes(x = new_pos, y = transf_pvalue, color = val)) + 
    geom_point(size = 1) + ggtitle('GWAS pvalues for CAD')+
    xlab('chromosome') +  ylab('-log10(pvalue)')+ theme_bw() +
    theme(plot.title = element_text(hjust = 0.5), text = element_text(size = 11),
          axis.text.x=element_text(size = 11, angle = 45),
          axis.text.y=element_text(size = 11), legend.position = 'bottom')+
    scale_x_continuous(breaks = df_start_end$start_plot, labels = c(1:22), limits = c(df_start_end$start_plot[1], sum(df_start_end$end)))+
    scale_color_gradient2(midpoint=0, low="blue", mid="grey",
                          high="red", space ="Lab")+
    labs(color = "reg. coefficient")
  
  ggsave(filename = sprintf('%smanhattanPlot_GWASsnps_genes_%s_%s_%s.png', fold, pheno, tissue, new_path), plot = pl_manh_snps, width = 10, height = 3, dpi=500, device = 'png')
  ggsave(filename = sprintf('%smanhattanPlot_GWASsnps_genes_%s_%s_%s.pdf', fold, pheno, tissue, new_path), plot = pl_manh_snps, width = 10, height = 3, device = 'pdf')
}
