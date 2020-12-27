# analysis conservation score and DNAase

inputFold <- '/psycl/g/mpsziller/lucia/PriLer_GTEx_model/AllTissues/200kb/noGWAS/'
DNAaseFold <- '/psycl/g/mpsziller/Analysis/SysMed/SNPannotation/'

# function to create dataframe for conservation score
create_df <- function(tissue_ann_file, dnase_file, name_type){
  
  snpsAnn <- read.table(gzfile(tissue_ann_file), h=T, stringsAsFactors = F)
  tissues_name <- c('AllTissues', colnames(snpsAnn[, !colnames(snpsAnn) %in% c('ID', 'ID_GTEx', 'ID_CAD', 'ID_PGC', 'chrom', 'position')]))
  
  dnaseAnn <- read.table(dnase_file, h=F, stringsAsFactors = F)
  colnames(dnaseAnn) <- c("chrom","chromstart","chromend","ID","chrom", "chromstart","chromend","identifier","mean_signal","numsamples","nOrder","Overlap")
  id_snpAnn <- paste(snpsAnn[,1], snpsAnn[,2], snpsAnn[,3], sep = '_')
  id_dnaaseAnn <- paste(dnaseAnn[,1], dnaseAnn[,2], dnaseAnn[,4], sep = '_')
  print(identical(id_snpAnn, id_dnaaseAnn))
  
  if(!identical(id_snpAnn, id_dnaaseAnn)){
    dnaseAnn$chrom_id <- as.numeric(sapply(dnaseAnn[,1], function(x) strsplit(x , split = 'chr')[[1]][2]))
    dnaseAnn <- dnaseAnn[order(dnaseAnn[,2]), ]
    dnaseAnn <- dnaseAnn[order(dnaseAnn$chrom_id), ]
    snpsAnn$chrom_id <-  as.numeric(sapply(snpsAnn$chrom, function(x) strsplit(x , split = 'chr')[[1]][2]))
    snpsAnn <- snpsAnn[order(snpsAnn$position), ]
    snpsAnn <- snpsAnn[order(snpsAnn$chrom_id), ]
    id_snpAnn <- paste(snpsAnn[,1], snpsAnn[,2], snpsAnn[,3], sep = '_')
    id_dnaaseAnn <- paste(dnaseAnn[,1], dnaseAnn[,2], dnaseAnn[,4], sep = '_')
    print(identical(id_snpAnn, id_dnaaseAnn))
  }
  
  dnaseAnn$numsamples_new <- dnaseAnn$numsamples
  dnaseAnn$numsamples_new[is.na(dnaseAnn$numsamples_new)] <- 0
  
  tmp <- lapply(tissues_name[-1], function(x) dnaseAnn$numsamples_new[snpsAnn[, x] == 1])
  tmp_tot <- dnaseAnn$numsamples_new[rowSums(snpsAnn[, tissues_name[-1]])>0]
  df <- data.frame(number_types = c(tmp_tot, unlist(tmp)), tissue = c(rep(tissues_name[1], length(tmp_tot)), 
                                                                      unlist(mapply(function(x,y) rep(x, length(y)), x = tissues_name[-1], y = tmp, SIMPLIFY = F))), stringsAsFactors = F)
  df$type <- name_type
  
  return(df)
  
}

compute_perc <- function(df, n_id){
  
  tissues_name <- unique(df$tissue)
  perc_df <- matrix(nrow = length(n_id), ncol = length(tissues_name))
  for(i in 1:length(n_id)){
    perc_df[i,] <- sapply(tissues_name, function(x) sum(df$tissue == x & df$number_types>n_id[i])/sum(df$tissue == x))
  }
  colnames(perc_df) <- tissues_name
  rownames(perc_df) <- n_id
  
  return(perc_df)
}


# Priler
df_Priler <- create_df(tissue_ann_file = sprintf('%sresPrior_regSNPs_annotation.txt.gz', inputFold), 
                       dnase_file = sprintf('%sresPrior_regSNPs_annotationregion_dnaseAnnotated.bed', DNAaseFold), name_type = 'PriLer')

# PriLer (reliableGenes)
df_Priler_relGenes <- create_df(tissue_ann_file = sprintf('%sresPrior_regSNPs_reliableGenes_annotation.txt.gz', inputFold), 
                                dnase_file = sprintf('%sresPrior_regSNPs_reliableGenes_annotationregion_dnaseAnnotated.bed', DNAaseFold), name_type = 'PriLer (reliable Genes)')
# PriLer (heritableGenes)
df_Priler_herGenes <- create_df(tissue_ann_file = sprintf('%sresPrior_regSNPs_heritableGenes_annotation.txt.gz', inputFold), 
                                dnase_file = sprintf('%sresPrior_regSNPs_heritableGenes_annotationregion_dnaseAnnotated.bed', DNAaseFold), name_type = 'PriLer (heritable Genes)')
# el-net
df_elnet <- create_df(tissue_ann_file = sprintf('%sresNoPrior_regSNPs_annotation.txt.gz', inputFold), 
                      dnase_file = sprintf('%sresNoPrior_regSNPs_annotationregion_dnaseAnnotated.bed', DNAaseFold), name_type = 'el-net')
# el-net (reliableGenes)
df_elnet_relGenes <- create_df(tissue_ann_file = sprintf('%sresNoPrior_regSNPs_reliableGenes_annotation.txt.gz', inputFold), 
                               dnase_file = sprintf('%sresNoPrior_regSNPs_reliableGenes_annotationregion_dnaseAnnotated.bed', DNAaseFold), name_type = 'el-net (reliable Genes)')
# el-net (heritableGenes)
df_elnet_herGenes <- create_df(tissue_ann_file = sprintf('%sresNoPrior_regSNPs_heritableGenes_annotation.txt.gz', inputFold), 
                               dnase_file = sprintf('%sresNoPrior_regSNPs_heritableGenes_annotationregion_dnaseAnnotated.bed', DNAaseFold), name_type = 'el-net (heritable Genes)')
# TWAS
df_TWAS <- create_df(tissue_ann_file = sprintf('%sTWAS_regSNPs_annotation.txt.gz', inputFold), 
                     dnase_file = sprintf('%sTWAS_regSNPs_annotationregion_dnaseAnnotated.bed', DNAaseFold), name_type = 'TWAS')
# prediXcan
df_prediXcan <- create_df(tissue_ann_file = sprintf('%sprediXcan_regSNPs_annotation.txt.gz', inputFold), 
                          dnase_file = sprintf('%sprediXcan_regSNPs_annotationregion_dnaseAnnotated.bed', DNAaseFold), name_type = 'prediXcan')

n_id <- c(0, 1, 2, 3, 4, 5, seq(10, 100, 10), seq(150, 700, 50))

# PriLer
tmp <- compute_perc(df_Priler, n_id)
df_Priler_perc <- data.frame(n_samples_thr = rep(n_id, length(ncol(tmp))), tissue = as.vector(sapply(colnames(tmp), function(x) rep(x,length(n_id)))), 
                             perc_higher = as.vector(tmp))
df_Priler_perc$type <- 'PriLer'

# PriLer (reliableGenes)
tmp <- compute_perc(df_Priler_relGenes, n_id)
df_Priler_relGenes_perc <- data.frame(n_samples_thr = rep(n_id, length(ncol(tmp))), tissue = as.vector(sapply(colnames(tmp), function(x) rep(x,length(n_id)))), 
                                      perc_higher = as.vector(tmp))
df_Priler_relGenes_perc$type <- 'PriLer (reliable Genes)'

# PriLer (heritableGenes)
tmp <- compute_perc(df_Priler_herGenes, n_id)
df_Priler_herGenes_perc <- data.frame(n_samples_thr = rep(n_id, length(ncol(tmp))), tissue = as.vector(sapply(colnames(tmp), function(x) rep(x,length(n_id)))), 
                                      perc_higher = as.vector(tmp))
df_Priler_herGenes_perc$type <- 'PriLer (heritable Genes)'

# el-net
tmp <- compute_perc(df_elnet, n_id)
df_elnet_perc <- data.frame(n_samples_thr = rep(n_id, length(ncol(tmp))), tissue = as.vector(sapply(colnames(tmp), function(x) rep(x,length(n_id)))), 
                            perc_higher = as.vector(tmp))
df_elnet_perc$type <- 'el-net'

# el-net (reliableGenes)
tmp <- compute_perc(df_elnet_relGenes, n_id)
df_elnet_relGenes_perc <- data.frame(n_samples_thr = rep(n_id, length(ncol(tmp))), tissue = as.vector(sapply(colnames(tmp), function(x) rep(x,length(n_id)))), 
                                     perc_higher = as.vector(tmp))
df_elnet_relGenes_perc$type <- 'el-net (reliable Genes)'

# el-net (heritableGenes)
tmp <- compute_perc(df_elnet_herGenes, n_id)
df_elnet_herGenes_perc <- data.frame(n_samples_thr = rep(n_id, length(ncol(tmp))), tissue = as.vector(sapply(colnames(tmp), function(x) rep(x,length(n_id)))), 
                                     perc_higher = as.vector(tmp))
df_elnet_herGenes_perc$type <- 'el-net (heritable Genes)'

# TWAS
tmp <- compute_perc(df_TWAS, n_id)
df_TWAS_perc <- data.frame(n_samples_thr = rep(n_id, length(ncol(tmp))), tissue = as.vector(sapply(colnames(tmp), function(x) rep(x,length(n_id)))), 
                           perc_higher = as.vector(tmp))
df_TWAS_perc$type <- 'TWAS'

# prediXcan
tmp <- compute_perc(df_prediXcan, n_id)
df_prediXcan_perc <- data.frame(n_samples_thr = rep(n_id, length(ncol(tmp))), tissue = as.vector(sapply(colnames(tmp), function(x) rep(x,length(n_id)))), 
                                perc_higher = as.vector(tmp))
df_prediXcan_perc$type <- 'prediXcan'

#### create general tables (at least 1 sample) ####
df_tot <- rbind(df_Priler_perc[df_Priler_perc$n_samples_thr == 0, -1], 
                df_Priler_relGenes_perc[df_Priler_relGenes_perc$n_samples_thr == 0, -1],
                df_Priler_herGenes_perc[df_Priler_herGenes_perc$n_samples_thr == 0, -1],
                df_elnet_perc[df_elnet_perc$n_samples_thr == 0, -1], 
                df_elnet_relGenes_perc[df_elnet_relGenes_perc$n_samples_thr == 0, -1], 
                df_elnet_herGenes_perc[df_elnet_herGenes_perc$n_samples_thr == 0, -1], 
                df_prediXcan_perc[df_prediXcan_perc$n_samples_thr == 0, -1], 
                df_TWAS_perc[df_TWAS_perc$n_samples_thr == 0, -1])

tissues <- as.character(unique(df_tot$tissue))[-1]
type_analysis <- as.character(unique(df_tot$type))
df_tot_mean <- data.frame(type = type_analysis, mean = sapply(type_analysis, function(x) mean(df_tot$perc_higher[df_tot$type == x & df_tot$tissue %in% tissues])),
                          sd = sapply(type_analysis, function(x) sd(df_tot$perc_higher[df_tot$type == x & df_tot$tissue %in% tissues])))

# save
write.table(sprintf('%spercentage_DNAase_intersection_summary.txt', inputFold), x = df_tot_mean, quote = F, col.names = T, row.names = F, sep = '\t')
write.table(sprintf('%spercentage_DNAase_intersection_tissueSpec.txt', inputFold), x = df_tot, quote = F, col.names = T, row.names = F, sep = '\t')

# plot distribution
library(ggplot2)
library(ggpubr)
options(bitmapType = 'cairo', device = 'png')


df_all <- rbind(df_Priler_relGenes, df_elnet_relGenes, df_TWAS, df_prediXcan)

quantile_sep <- seq(0, 1, by = 0.001)
df_plot <- data.frame(n_samples = c(), quant =c(), tissues = c(), type = c(), stringsAsFactors = F)
type_analysis <- c('PriLer (reliable Genes)','el-net (reliable Genes)','prediXcan', 'TWAS')
df_test <- data.frame(tissues = c(), type = c(), wilcox_pvalue = c(), mean = c(), stringsAsFactors = F)

for(i in 1:length(type_analysis)){
  print(i)
  for(j in 1:length(tissues)){
    
    tmp <- quantile(df_all$number_types[df_all$tissue == tissues[j] & df_all$type == type_analysis[i]], probs = quantile_sep)
    df_plot <- rbind(df_plot, data.frame(stringsAsFactors = F, n_samples = tmp, quant = quantile_sep*100, tissues = rep(tissues[j], length(tmp)), 
                                         type = rep(type_analysis[i], length(tmp))))
    
    # test difference with PriLer method
    df_test <- rbind(df_test, data.frame(stringsAsFactors = F,tissues = tissues[j], type = type_analysis[i], 
                                         wilcox_pvalue = wilcox.test(df_all$number_types[df_all$tissue == tissues[j] & df_all$type == type_analysis[i]], 
                                                                     df_all$number_types[df_all$tissue == tissues[j] & df_all$type == type_analysis[1]])$p.value, 
                                         mean = mean(df_all$number_types[df_all$tissue == tissues[j] & df_all$type == type_analysis[i]])))
    
  }
  
  tmp <- quantile(df_all$number_types[df_all$tissue == 'AllTissues' & df_all$type == type_analysis[i]], probs = quantile_sep)
  df_plot <- rbind(df_plot, data.frame(n_samples = tmp, quant = quantile_sep*100, tissues = rep( 'AllTissues', length(tmp)), 
                                       type = rep(type_analysis[i], length(tmp))))
  # test difference with PriLer method
  df_test <- rbind(df_test, data.frame(stringsAsFactors = F,tissues = 'AllTissues', type = type_analysis[i], 
                                       wilcox_pvalue = wilcox.test(df_all$number_types[df_all$tissue == 'AllTissues' & df_all$type == type_analysis[i]], 
                                                                   df_all$number_types[df_all$tissue == 'AllTissues' & df_all$type == type_analysis[1]])$p.value, 
                                       mean = mean(df_all$number_types[df_all$tissue == 'AllTissues' & df_all$type == type_analysis[i]])))
  
}

df_plot$tissues <- factor(df_plot$tissues, levels = c('AllTissues', tissues))
df_plot$type[df_plot$type == type_analysis[1]] <- 'PriLer'
df_plot$type[df_plot$type == type_analysis[2]] <- 'el-net'
df_plot$type <- factor(df_plot$type, levels = c('PriLer', 'el-net', 'prediXcan', 'TWAS'))
df_plot$log_nsamples <- log2(df_plot$n_samples + 1)
df_test_ks <- data.frame(tissues = c(), type = c(), thr_quant = c(), ks_pvalue = c(), stringsAsFactors = F)
type_analysis <- c('PriLer', 'el-net', 'prediXcan', 'TWAS')
for(i in 2:length(type_analysis)){
  
  print(i)
  for(j in 1:length(tissues)){
    tmp <-  df_plot[df_plot$type %in% c(type_analysis[1],type_analysis[i]) & df_plot$tissue == tissues[j] & df_plot$n_samples>0, ]
    thr <- min(tmp$quant)
    df_test_ks <- rbind(df_test_ks, data.frame(tissues = tissues[j], type = type_analysis[i], thr_quant = thr, 
                                               ks_pvalue = ks.test(x = df_plot$n_samples[df_plot$quant >= thr & df_plot$type == type_analysis[1] & df_plot$tissue == tissues[j]],
            y = df_plot$n_samples[df_plot$quant >= thr & df_plot$type == type_analysis[i] & df_plot$tissue == tissues[j]], alternative = 'less')$p.value))
  }
  
  tmp <-  df_plot[df_plot$type %in% c(type_analysis[1],type_analysis[i]) & df_plot$tissue == 'AllTissues' & df_plot$n_samples>0, ]
  thr <- min(tmp$quant)
  df_test_ks <- rbind(df_test_ks, data.frame(tissues = 'AllTissues', 
                                             type = type_analysis[i],  thr_quant = thr, 
                                             ks_pvalue = ks.test(x = df_plot$n_samples[df_plot$quant >= thr & df_plot$type == type_analysis[1] & df_plot$tissue == 'AllTissues'],
                                                                 y = df_plot$n_samples[df_plot$quant >= thr & df_plot$type == type_analysis[i] & df_plot$tissue == 'AllTissues'], alternative = 'less')$p.value))
  
}
write.table(x = df_test_ks, file = sprintf('%sKS_test_distribution_notzero.txt', inputFold), quote = F, col.names = T, row.names = F, sep = '\t')

df_test_ks$value <- paste(df_test_ks$type,  formatC(df_test_ks$ks_pvalue, format = "e", digits = 1))
df_test_ks$position <- 74
df_test_ks$position[df_test_ks$type == 'prediXcan'] <- 68
df_test_ks$position[df_test_ks$type == 'TWAS'] <- 62
df_test_ks$tissue <- factor(df_test_ks$tissues, levels = c('AllTissues', tissues))
df_test_ks$x <- 6.5
df_test_ks$type <- factor(df_test_ks$type, level = type_analysis)

# color_type <- c('#1F77B4FF', '#339733', '#FF7F0EFF', '#D62728FF')
color_type <- c('#1F77B4FF', 'darkgrey', '#FF7F0EFF', '#D62728FF')

pl <- ggplot(data = subset(df_plot, tissues == 'AllTissues'), aes(y = quant, x = log_nsamples, color = type))+
  geom_line()+
  theme_bw()+ ylim(58, 100)+
  ylab('percentile')+xlab('log2(n. biosamples + 1)')+ 
  theme(legend.position = 'right', plot.title = element_text(size=9))+
  annotate("text", x = 7.5, y = df_test_ks$position[df_test_ks$tissue == 'AllTissues'], label = df_test_ks$value[df_test_ks$tissue == 'AllTissues'], color = c('darkgrey', '#FF7F0EFF', '#D62728FF'))+
  scale_color_manual(values=color_type)
ggsave(filename = sprintf('%spercentile_nsamples_DNAase_AllTissues.png', inputFold), width = 6, height = 5, plot = pl, device = 'png')
ggsave(filename = sprintf('%spercentile_nsamples_DNAase_AllTissues.pdf', inputFold), width = 6, height = 5, plot = pl, device = 'pdf')

pl <- ggplot(data = subset(df_plot, tissues != 'AllTissues'), aes(y = quant, x = log_nsamples, color = type))+
  geom_line()+
  theme_bw()+ ylim(58, 100)+
  facet_wrap(.~tissues, ncol = 6)+
  ylab('percentile')+xlab('log2(n. biosamples + 1)')+
  geom_text(data = subset(df_test_ks, tissue != 'AllTissues'), mapping = aes(x = x, y = position, label = value, color = type), size = 3)+
  theme(legend.position = 'bottom', strip.text = element_text(size = 8))+
  scale_color_manual(values=color_type)
ggsave(filename = sprintf('%spercentile_nsamples_DNAase_TissueSpec.png', inputFold), width = 12, height = 10, plot = pl, device = 'png')
ggsave(filename = sprintf('%spercentile_nsamples_DNAase_TissueSpec.pdf', inputFold), width = 12, height = 10, plot = pl, device = 'pdf')


