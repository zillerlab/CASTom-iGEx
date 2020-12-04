# random prior plots (only in GTEx Artery_Coronary and Brain_Cortex)

library(ggplot2)
library(ggrepel)
library(gridExtra)
library(latex2exp)

setwd('/mnt/lucia/PriLer_PROJECT_GTEx/OUTPUT_SCRIPTS_v2')
df_info <- data.frame(tissue = c('Artery_Coronary', 'Brain_Cortex'), fold_epi = c('200kb/CAD_GWAS_bin5e-2_randomEpi/', '200kb/PGC_GWAS_bin1e-2_randomEpi/'), 
                      fold_var = c('200kb/CAD_GWAS_bin5e-2_randomVar/', '200kb/PGC_GWAS_bin1e-2_randomVar/'),
                      fold_gwas = c('200kb/CAD_GWAS_bin5e-2_randomGWAS/', '200kb/PGC_GWAS_bin1e-2_randomGWAS/'), 
                      file_prior = c('priorName_CADgwas_withIndex.txt', 'priorName_PGCgwas_withIndex.txt'), stringsAsFactors = F)

random_prior <- 'priorMatrix_random_Ctrl_150_allPeaks_allRanger_heart_left_ventricle_GWAS_withRep_'
nrep_rV <- 50
nrep_rE <- 50
nrep_rG <- 10
conv_par <- 0.25
outFold <- '/mnt/lucia/PriLer_TRAIN_PLOTS/randomPrior_res/'

randomPriorMat <- list()
for(chr in 1:22){
  print(chr)
  randomPriorMat[[chr]] <- read.table(gzfile(sprintf('%schr%i.txt.gz', random_prior, chr)), h=T, stringsAsFactors = F)
}

fixed_prior <- list()
for(i in 1:nrow(df_info)){
  
  tmp_file <- read.table(sprintf('%s/%s/rep1/%s', df_info$tissue[i], df_info$fold_var[i], df_info$file_prior[i]), stringsAsFactors = F, h=F)
  id <- sapply(tmp_file$V2, function(x) length(strsplit(x, split = 'random')[[1]]) == 1)
  fixed_prior[[i]] <- data.frame(prior = tmp_file$V2[id], stringsAsFactors = F)
  
}


# count number of snps with fixed prior
nchr <- list()
for(chr in 1:22){
  
  tmp <- lapply(fixed_prior, function(x) randomPriorMat[[chr]][, x$prior])
  nchr[[chr]] <- lapply(tmp, function(x) colSums(x!=0))
  
}
nprior_fixed <- list()
for(i in 1:nrow(df_info)){
  
  nprior_fixed[[i]] <- rowSums(sapply(nchr, function(x) x[[i]]))
  
}


####################
#### random Var ####
####################

rV_prior <- list()
for(i in 1:nrow(df_info)){
  
  tmp_file <- read.table(sprintf('%s/%s/rep%i/%s', df_info$tissue[i], df_info$fold_var[i], 1, df_info$file_prior[i]), stringsAsFactors = F, h=F)
  id <- sapply(tmp_file$V2, function(x) length(strsplit(x, split = 'random')[[1]]) > 1)
  rV_prior[[i]] <- matrix(tmp_file$V2[id], ncol=1)
  
  for(j in 2:nrep_rV){
    tmp_file <- read.table(sprintf('%s/%s/rep%i/%s', df_info$tissue[i], df_info$fold_var[i], j, df_info$file_prior[i]), stringsAsFactors = F, h=F)
    id <- sapply(tmp_file$V2, function(x) length(strsplit(x, split = 'random')[[1]]) > 1)
    rV_prior[[i]] <- cbind(rV_prior[[i]], tmp_file$V2[id])
  }
}

# count number of snps with randomVar prior
nchr_rV <- list()
for(chr in 1:22){
  
  nchr_rV[[chr]] <- lapply(rV_prior, function(x) apply(x, 2,  function(y) colSums(randomPriorMat[[chr]][,y]!=0)))
  
}

nprior_rV <- list()
for(i in 1:nrow(df_info)){
  
  nprior_rV[[i]] <- matrix(ncol = nrep_rV, nrow = nrow(rV_prior[[i]]))
  tmp <- lapply(nchr_rV, function(x) x[[i]])
  for(j in 1:nrep_rV){
    nprior_rV[[i]][,j] <- rowSums(sapply(tmp, function(x) x[,j]))    
  }
}

# count the number of snps that are shared across the epi for the random comfiguration
df_int <- list()
perc_int <- list()
for(i in 1:nrow(df_info)){
  
  df_int[[i]] <- matrix(ncol = nrep_rV, nrow = nrow(rV_prior[[i]]))
  
  tmp_fixed <- do.call(rbind, lapply(randomPriorMat, function(x) x[,fixed_prior[[i]]$prior]!=0))
  tmp_fixed_all <- apply(tmp_fixed, 1, any)
  
  for(j in 1:nrep_rV){
    
    tmp_rV <-  do.call(rbind, lapply(randomPriorMat, function(x) x[,rV_prior[[i]][,j]]!=0))
    df_int[[i]][,j] <- apply(tmp_rV, 2, function(x) length(which(x & tmp_fixed_all)))
    
  }
  
  perc_int[[i]] <- df_int[[i]]/nprior_rV[[i]]
  
}

# for each repetiton load final results
df_res <- list()
for(i in 1:nrow(df_info)){
  
  tmp_file <- read.table(sprintf('%s/%s/rep1/%s', df_info$tissue[i], df_info$fold_var[i], df_info$file_prior[i]), stringsAsFactors = F, h=F)
  names_p <- as.vector(sapply(tmp_file$V2, function(x) strsplit(x, split = '_r1')[[1]][1]))
  tmp <- matrix(ncol = nrep_rV, nrow = length(names_p))
  
  for(j in 1:nrep_rV){
    print(j)
    res_w <- get(load(sprintf('%s/%s/rep%i/resPrior_EOpt_HeritableGenes_allchr.RData', df_info$tissue[i], df_info$fold_var[i], j)))
    
    tmp[,j] <- res_w$weights
  }
  
  df_res[[i]] <- data.frame(prior = names_p, stringsAsFactors = F)
  df_res[[i]]$mean_w <- rowMeans(tmp)*conv_par
  df_res[[i]]$sd_w <- apply(tmp, 1, sd)
  
}

# weights iteration
df_it <- list()
for(i in 1:nrow(df_info)){
  
  tmp_file <- read.table(sprintf('%s/%s/rep1/%s', df_info$tissue[i], df_info$fold_var[i], df_info$file_prior[i]), stringsAsFactors = F, h=F)
  names_p <- as.vector(sapply(tmp_file$V2, function(x) strsplit(x, split = '_r1')[[1]][1]))
  tmp <- array(0,dim=c(20,length(names_p),nrep_rV))
  
  
  for(j in 1:nrep_rV){
    print(j)
    res_it <- get(load(sprintf('%s/%s/rep%i/resPrior_EOpt_Iteration_HeritableGenes_allchr.RData', df_info$tissue[i], df_info$fold_var[i], j)))
    tmp[,,j] <- as.matrix(res_it$pWeight, ncol = length(names_p))
  }
  
  df_it[[i]] <- data.frame(mean_w = as.vector(apply(tmp, c(1,2), mean))*conv_par, sd_w = as.vector(apply(tmp, c(1,2), sd)), stringsAsFactors = F)
  df_it[[i]]$it <- rep(1:20, length(names_p))
  df_it[[i]]$prior <- as.vector(sapply(names_p, function(x) rep(x, 20)))
  
}

# build data frame
# f <- function(pal) brewer.pal(brewer.pal.info[pal, "maxcolors"], pal)
for(i in 1:nrow(df_info)){
  
  df_rV_n <- data.frame(prior = c(fixed_prior[[i]]$prior, as.vector(sapply(rV_prior[[i]][,1], function(x) strsplit(x, split = '_r1')[[1]][1]))), stringsAsFactors = F)
  df_rV_n$nsnps <- c(nprior_fixed[[i]], rowMeans(nprior_rV[[i]]))
  df_rV_n$sd <- rep(0, length(nprior_fixed[[i]])+2)
  df_rV_n$prior <- factor(df_rV_n$prior, levels = rev(df_rV_n$prior))
  
  df_perc <- data.frame(prior = rep(as.vector(sapply(rV_prior[[i]][,1], function(x) strsplit(x, split = '_r1')[[1]][1])), 2), stringsAsFactors = F)
  df_perc$perc <- c(rowMeans(perc_int[[i]]), rowMeans(1- perc_int[[i]]))
  df_perc$sd <- c(apply(perc_int[[i]], 1, sd), apply(1 -perc_int[[i]], 1, sd))
  df_perc$type <- c(rep('shared', nrow(perc_int[[i]])), rep('unique', nrow(perc_int[[i]])))  
  df_perc$prior <- factor( df_perc$prior, levels =  rev(unique(df_perc$prior)))
  
  if(df_info$tissue[i] == 'Artery_Coronary'){
    col <- rev(c(rep('grey', 2), '#FA8072', rep('grey', 4), rep('#FA8072', 2)))
    col_per <- c('#fac9c3', '#FA8072')
    # shape_it <- c(rep(16, 7), 15,17)
    height_pl <- 3
  }
  if(df_info$tissue[i] == 'Brain_Cortex'){
    col <- rev(c(rep('grey', 13), '#4169E1', 'grey', rep('#4169E1', 2)))
    col_per <- c('#a3b3e3', '#4169E1')
    height_pl <- 4.5
    # shape_it <- c(rep(16, 15), 15,17)
  }
  
  
  pl_nsnps <- ggplot(df_rV_n, aes(x = prior, y=nsnps, fill = prior)) + 
    geom_bar(stat="identity", color="black") +ggtitle('Random Variants')+
    ylab('n. SNPs with prior')+ xlab('')+ theme_classic()+
    geom_errorbar(aes(ymin=nsnps-sd, ymax=nsnps+sd), width=.2, position=position_dodge(.9))+
    theme(plot.title = element_text(hjust = 0.5), text = element_text(size = 10),
          axis.text.x=element_text(size = 10, angle = 0, hjust = 1),
          axis.text.y=element_text(size = 10), legend.position = 'none')+
    scale_fill_manual(values = col)+coord_flip()
  ggsave(filename = sprintf('%snSNPswithPrior_rV_%s.png', outFold, df_info$tissue[i]), plot = pl_nsnps, width = 6, height = height_pl, dpi = 500)
  ggsave(filename = sprintf('%snSNPswithPrior_rV_%s.pdf', outFold, df_info$tissue[i]), plot = pl_nsnps, width = 6, height = height_pl)
  
  pl_perc <- ggplot(df_perc, aes(x = prior, y=perc, fill = type)) + 
    geom_bar(stat="identity", color="black") +ggtitle('Random Variants')+
    ylab('mean percentage')+ xlab('')+ theme_classic()+ylim(0,1)+
    # geom_errorbar(aes(ymin=perc-sd, ymax=perc+sd), width=.2)+
    theme(plot.title = element_text(hjust = 0.5), text = element_text(size = 10),
          axis.text.x=element_text(size = 10, angle = 0, hjust = 1),
          axis.text.y=element_text(size = 10), legend.position = 'bottom')+
    scale_fill_manual(values = col_per)+
    labs(fill = " ") + coord_flip()
  ggsave(filename = sprintf('%spercShared_SNPswithPrior_rV_%s.png', outFold,  df_info$tissue[i]), plot = pl_perc, width = 6, height = 2, dpi = 500)
  ggsave(filename = sprintf('%spercShared_SNPswithPrior_rV_%s.pdf', outFold,  df_info$tissue[i]), plot = pl_perc, width = 6, height = 2)
  
  df_res[[i]]$prior <- factor(df_res[[i]]$prior, levels = rev(df_rV_n$prior)) 
  pl_w <- ggplot(df_res[[i]], aes(x = prior, y=mean_w, fill = prior)) + 
    geom_bar(stat="identity", color="black") +ggtitle('Random Variants')+
    ylab('Maximum contribution to prior coefficient')+ xlab('')+ theme_classic()+
    geom_errorbar(aes(ymin=mean_w-sd_w, ymax=mean_w+sd_w), width=.2, position=position_dodge(.9))+
    theme(plot.title = element_text(hjust = 0.5), text = element_text(size = 10),
          axis.text.x=element_text(size = 10, angle = 0, hjust = 1),
          axis.text.y=element_text(size = 10), legend.position = 'none')+
    scale_fill_manual(values = col)+coord_flip()
  ggsave(filename = sprintf('%smaxContrPrior_rV_%s.png', outFold, df_info$tissue[i]), plot = pl_w, width = 6, height = height_pl, dpi = 500)
  ggsave(filename = sprintf('%smaxContrPrior_rV_%s.pdf', outFold, df_info$tissue[i]), plot = pl_w, width = 6, height = height_pl)
  
  df_it[[i]]$prior <- factor(df_it[[i]]$prior, levels = df_rV_n$prior)
  df_it[[i]]$text <- ""
  id <- df_it[[i]]$prior %in% df_rV_n$prior[rev(col)!='grey'] & df_it[[i]]$it == 1
  df_it[[i]]$text[id] <- as.character(df_it[[i]]$prior)[id]
  
  pl_it <- ggplot(df_it[[i]], aes(x = it, y = mean_w, color = prior, label = text)) + 
    geom_point(size=0.7, alpha=0.8) + geom_line(alpha = 0.8)+ ggtitle('Random Variants')+
    geom_errorbar(aes(ymin=mean_w-sd_w, ymax=mean_w+sd_w), width=.2, position=position_dodge(0.05))+
    geom_text_repel(size = 2.5, segment.size = 0.2, segment.color = 'black', force = 10, nudge_x = 10, aes(fontface = 2))+
    theme_classic() + xlab('iteration')+  
    ylab('Maximum contribution to prior coefficient') + theme(legend.position = 'none', plot.title = element_text(hjust = 0.5))+
    scale_color_manual(values = rev(col))+labs(color = "")
  ggsave(filename = sprintf('%smaxContrPrior_iteration_rV_%s.png', outFold, df_info$tissue[i]), plot = pl_it, width = 4, height = 4, dpi=500)
  ggsave(filename = sprintf('%smaxContrPrior_iteration_rV_%s.pdf', outFold, df_info$tissue[i]), plot = pl_it, width = 4, height = 4)
  
  res <- list(weights = df_res[[i]], nsnps = df_rV_n, it = df_it[[i]], perc_shared = df_perc)
  save(res, file = sprintf('%stotRes_rV_%s.RData', outFold, df_info$tissue[i]))
  
}


####################
#### random Epi ####
####################

rE_prior <- list()
for(i in 1:nrow(df_info)){
  
  add_name <- ifelse(df_info$tissue[i] == 'Artery_Coronary', 'Ctrl_150_allPeaks_cellRanger', 'heart_left_ventricle')
  tmp_file <- read.table(sprintf('%s/%s/rep%i/%s', df_info$tissue[i], df_info$fold_epi[i], 1, df_info$file_prior[i]), stringsAsFactors = F, h=F)
  id <- sapply(tmp_file$V2, function(x) length(strsplit(x, split = 'random')[[1]]) > 1) | tmp_file$V2 %in% add_name
  rE_prior[[i]] <- matrix(tmp_file$V2[id], ncol=1)
  
  for(j in 2:nrep_rE){
    tmp_file <- read.table(sprintf('%s/%s/rep%i/%s', df_info$tissue[i], df_info$fold_epi[i], j, df_info$file_prior[i]), stringsAsFactors = F, h=F)
    id <- sapply(tmp_file$V2, function(x) length(strsplit(x, split = 'random')[[1]]) > 1) | tmp_file$V2 %in% add_name
    rE_prior[[i]] <-cbind(rE_prior[[i]], tmp_file$V2[id])
  }
}

# count number of snps with randomVar prior
nchr_rE <- list()
for(chr in 1:22){
  
  nchr_rE[[chr]] <- lapply(rE_prior, function(x) apply(x, 2,  function(y) colSums(randomPriorMat[[chr]][,y]!=0)))
  
}

nprior_rE <- list()
for(i in 1:nrow(df_info)){
  
  nprior_rE[[i]] <- matrix(ncol = nrep_rE, nrow = nrow(rE_prior[[i]]))
  tmp <- lapply(nchr_rE, function(x) x[[i]])
  for(j in 1:nrep_rE){
    nprior_rE[[i]][,j] <- rowSums(sapply(tmp, function(x) x[,j]))    
  }
}

# count the number of snps that are shared across the epi for the random comfiguration
df_int <- list()
perc_int <- list()
for(i in 1:nrow(df_info)){
  
  df_int[[i]] <- matrix(ncol = nrep_rE, nrow = nrow(rE_prior[[i]]))
  
  tmp_fixed <- do.call(rbind, lapply(randomPriorMat, function(x) x[,fixed_prior[[i]]$prior]!=0))
  tmp_fixed_all <- apply(tmp_fixed, 1, any)
  
  for(j in 1:nrep_rE){
    
    tmp_rE <-  do.call(rbind, lapply(randomPriorMat, function(x) x[,rE_prior[[i]][,j]]!=0))
    df_int[[i]][,j] <- apply(tmp_rE, 2, function(x) length(which(x & tmp_fixed_all)))
    
  }
  
  perc_int[[i]] <- df_int[[i]]/nprior_rE[[i]]
  
}

# for each repetiton load final results
df_res <- list()
for(i in 1:nrow(df_info)){
  
  tmp_file <- read.table(sprintf('%s/%s/rep1/%s', df_info$tissue[i], df_info$fold_epi[i], df_info$file_prior[i]), stringsAsFactors = F, h=F)
  names_p <- as.vector(sapply(tmp_file$V2, function(x) strsplit(x, split = '_r1')[[1]][1]))
  tmp <- matrix(ncol = nrep_rV, nrow = length(names_p))
  
  for(j in 1:nrep_rV){
    print(j)
    res_w <- get(load(sprintf('%s/%s/rep%i/resPrior_EOpt_HeritableGenes_allchr.RData', df_info$tissue[i], df_info$fold_epi[i], j)))
    
    tmp[,j] <- res_w$weights
  }
  
  df_res[[i]] <- data.frame(prior = names_p, stringsAsFactors = F)
  df_res[[i]]$mean_w <- rowMeans(tmp)*conv_par
  df_res[[i]]$sd_w <- apply(tmp, 1, sd)
  
}

# weights iteration
df_it <- list()
for(i in 1:nrow(df_info)){
  
  tmp_file <- read.table(sprintf('%s/%s/rep1/%s', df_info$tissue[i], df_info$fold_epi[i], df_info$file_prior[i]), stringsAsFactors = F, h=F)
  names_p <- as.vector(sapply(tmp_file$V2, function(x) strsplit(x, split = '_r1')[[1]][1]))
  tmp <- array(0,dim=c(20,length(names_p),nrep_rV))
  
  
  for(j in 1:nrep_rV){
    print(j)
    res_it <- get(load(sprintf('%s/%s/rep%i/resPrior_EOpt_Iteration_HeritableGenes_allchr.RData', df_info$tissue[i], df_info$fold_epi[i], j)))
    tmp[,,j] <- as.matrix(res_it$pWeight, ncol = length(names_p))
  }
  
  df_it[[i]] <- data.frame(mean_w = as.vector(apply(tmp, c(1,2), mean))*conv_par, sd_w = as.vector(apply(tmp, c(1,2), sd)), stringsAsFactors = F)
  df_it[[i]]$it <- rep(1:20, length(names_p))
  df_it[[i]]$prior <- as.vector(sapply(names_p, function(x) rep(x, 20)))
  
}

# build data frame
# f <- function(pal) brewer.pal(brewer.pal.info[pal, "maxcolors"], pal)
for(i in 1:nrow(df_info)){
  
  df_rE_n <- data.frame(prior = c(fixed_prior[[i]]$prior, as.vector(sapply(rE_prior[[i]][,1], function(x) strsplit(x, split = '_r1')[[1]][1]))), stringsAsFactors = F)
  df_rE_n$nsnps <- c(nprior_fixed[[i]], rowMeans(nprior_rE[[i]]))
  df_rE_n$sd <- c(rep(0, length(nprior_fixed[[i]])+1), apply(nprior_rE[[i]][-1,], 1, sd))
  df_rE_n$prior <- factor(df_rE_n$prior, levels = rev(df_rE_n$prior))
  
  df_perc <- data.frame(prior = rep(as.vector(sapply(rE_prior[[i]][,1], function(x) strsplit(x, split = '_r1')[[1]][1])), 2), stringsAsFactors = F)
  df_perc$perc <- c(rowMeans(perc_int[[i]]), rowMeans(1-perc_int[[i]]))
  df_perc$sd <- c(apply(perc_int[[i]], 1, sd), apply(1-perc_int[[i]], 1, sd))
  df_perc$type <- c('shared (other)',rep('shared', nrow(perc_int[[i]])-1), 'unique (other)',rep('unique', nrow(perc_int[[i]])-1))  
  df_perc$type <- factor(df_perc$type, levels = c('shared (other)', 'shared', 'unique (other)', 'unique'))
  df_perc$prior <- factor( df_perc$prior, levels =  rev(unique(df_perc$prior)))
  
  if(df_info$tissue[i] == 'Artery_Coronary'){
    col <- rev(c(rep('grey', 2), '#FA8072', rep('grey', 4), '#4169E1', rep('#FA8072', 3)))
    col_per <- c('#a3b3e3', '#fac9c3',  '#4169E1', '#FA8072')
    # shape_it <- c(rep(16, 7), 15,17)
    height_pl <- 3
  }
  if(df_info$tissue[i] == 'Brain_Cortex'){
    col <- rev(c(rep('grey', 13), '#4169E1', 'grey', '#FA8072', rep('#4169E1', 3)))
    col_per <- c('#fac9c3', '#a3b3e3', '#FA8072', '#4169E1')
    height_pl <- 4.5
    # shape_it <- c(rep(16, 15), 15,17)
  }
  
  pl_nsnps <- ggplot(df_rE_n, aes(x = prior, y=nsnps, fill = prior)) + 
    geom_bar(stat="identity", color="black") +ggtitle('Random GREs')+
    ylab('n. SNPs with prior')+ xlab('')+ theme_classic()+
    geom_errorbar(aes(ymin=nsnps-sd, ymax=nsnps+sd), width=.2, position=position_dodge(.9))+
    theme(plot.title = element_text(hjust = 0.5), text = element_text(size = 10),
          axis.text.x=element_text(size = 10, angle = 0, hjust = 1),
          axis.text.y=element_text(size = 10), legend.position = 'none')+
    scale_fill_manual(values = col)+coord_flip()
  ggsave(filename = sprintf('%snSNPswithPrior_rE_%s.png',outFold, df_info$tissue[i]), plot = pl_nsnps, width = 6, height = height_pl, dpi = 500)
  ggsave(filename = sprintf('%snSNPswithPrior_rE_%s.pdf',outFold, df_info$tissue[i]), plot = pl_nsnps, width = 6, height = height_pl)
  
  pl_perc <- ggplot(df_perc, aes(x = prior, y=perc, fill = type)) + 
    geom_bar(stat="identity", color="black") +ggtitle('Random GREs')+
    ylab('mean percentage')+ xlab('')+ theme_classic()+ylim(0,1)+
    # geom_errorbar(aes(ymin=perc-sd, ymax=perc+sd), width=.2)+
    theme(plot.title = element_text(hjust = 0.5), text = element_text(size = 10),
          axis.text.x=element_text(size = 10, angle = 0, hjust = 1),
          axis.text.y=element_text(size = 10), legend.position = 'bottom')+
    scale_fill_manual(values = col_per)+
    labs(fill = " ") + coord_flip()+guides(fill=guide_legend(nrow=2,byrow=TRUE))
  ggsave(filename = sprintf('%spercShared_SNPswithPrior_rE_%s.png', outFold, df_info$tissue[i]), plot = pl_perc, width = 6, height = 2.6, dpi = 500)
  ggsave(filename = sprintf('%spercShared_SNPswithPrior_rE_%s.pdf', outFold, df_info$tissue[i]), plot = pl_perc, width = 6, height = 2.6)
  
  df_res[[i]]$prior <- factor(df_res[[i]]$prior, levels = rev(df_rE_n$prior)) 
  pl_w <- ggplot(df_res[[i]], aes(x = prior, y=mean_w, fill = prior)) + 
    geom_bar(stat="identity", color="black") +ggtitle('Random GREs')+
    ylab('Maximum contribution to prior coefficient')+ xlab('')+ theme_classic()+
    geom_errorbar(aes(ymin=mean_w-sd_w, ymax=mean_w+sd_w), width=.2, position=position_dodge(.9))+
    theme(plot.title = element_text(hjust = 0.5), text = element_text(size = 10),
          axis.text.x=element_text(size = 10, angle = 0, hjust = 1),
          axis.text.y=element_text(size = 10), legend.position = 'none')+
    scale_fill_manual(values = col)+coord_flip()
  ggsave(filename = sprintf('%smaxContrPrior_rE_%s.png', outFold, df_info$tissue[i]), plot = pl_w, width = 6, height = height_pl, dpi = 500)
  ggsave(filename = sprintf('%smaxContrPrior_rE_%s.pdf', outFold, df_info$tissue[i]), plot = pl_w, width = 6, height = height_pl)
  
  df_it[[i]]$prior <- factor(df_it[[i]]$prior, levels = df_rE_n$prior)
  df_it[[i]]$text <- ""
  id <- df_it[[i]]$prior %in% df_rE_n$prior[rev(col)!='grey'] & df_it[[i]]$it == 1
  df_it[[i]]$text[id] <- as.character(df_it[[i]]$prior)[id]
  
  pl_it <- ggplot(df_it[[i]], aes(x = it, y = mean_w, color = prior, label = text)) + 
    geom_point(size=0.7, alpha=0.8) + geom_line(alpha = 0.8)+ ggtitle('Random GREs')+
    geom_errorbar(aes(ymin=mean_w-sd_w, ymax=mean_w+sd_w), width=.2, position=position_dodge(0.05))+
    geom_text_repel(size = 2.5, segment.size = 0.2, segment.color = 'black', force = 10, nudge_x = 10, aes(fontface = 2))+
    theme_classic() + xlab('iteration')+  
    ylab('Maximum contribution to prior coefficient') + theme(legend.position = 'none', plot.title = element_text(hjust = 0.5))+
    scale_color_manual(values = rev(col))+labs(color = "")
  ggsave(filename = sprintf('%smaxContrPrior_iteration_rE_%s.png', outFold, df_info$tissue[i]), plot = pl_it, width = 4, height = 4, dpi=500)
  ggsave(filename = sprintf('%smaxContrPrior_iteration_rE_%s.pdf', outFold, df_info$tissue[i]), plot = pl_it, width = 4, height = 4)
  
  res <- list(weights = df_res[[i]], nsnps = df_rE_n, it = df_it[[i]], perc_shared = df_perc)
  save(res, file = sprintf('%stotRes_rE_%s.RData', outFold, df_info$tissue[i]))
  
}

#####################
#### random GWAS ####
#####################

rG_prior <- list()
for(i in 1:nrow(df_info)){
  
  tmp_file <- read.table(sprintf('%s/%s/rep%i/%s', df_info$tissue[i], df_info$fold_gwas[i], 1, df_info$file_prior[i]), stringsAsFactors = F, h=F)
  id <- sapply(tmp_file$V2, function(x) length(strsplit(x, split = 'random')[[1]]) > 1)
  rG_prior[[i]] <- matrix(tmp_file$V2[id], ncol=1)
  
  for(j in 2:nrep_rG){
    tmp_file <- read.table(sprintf('%s/%s/rep%i/%s', df_info$tissue[i], df_info$fold_gwas[i], j, df_info$file_prior[i]), stringsAsFactors = F, h=F)
    id <- sapply(tmp_file$V2, function(x) length(strsplit(x, split = 'random')[[1]]) > 1) 
    rG_prior[[i]] <-cbind(rG_prior[[i]], tmp_file$V2[id])
  }
}

# count number of snps with randomVar prior
nchr_rG <- list()
for(chr in 1:22){
  
  nchr_rG[[chr]] <- lapply(rG_prior, function(x) apply(x, 2,  function(y) colSums(randomPriorMat[[chr]][,y]!=0)))
  
}

nprior_rG <- list()
for(i in 1:nrow(df_info)){
  
  nprior_rG[[i]] <- matrix(ncol = nrep_rG, nrow = nrow(rG_prior[[i]]))
  tmp <- lapply(nchr_rG, function(x) x[[i]])
  for(j in 1:nrep_rG){
    nprior_rG[[i]][,j] <- rowSums(sapply(tmp, function(x) x[,j]))    
  }
}

# count the number of snps that are shared across the epi for the random comfiguration
df_int <- list()
df_int_g <- list()
perc_int <- list()
perc_int_g <- list()
for(i in 1:nrow(df_info)){
  
  df_int[[i]] <- matrix(ncol = nrep_rG, nrow = nrow(rG_prior[[i]]))
  df_int_g[[i]] <- matrix(ncol = nrep_rG, nrow = nrow(rG_prior[[i]]))
  
  gwas_name <- fixed_prior[[i]]$prior[sapply(fixed_prior[[i]]$prior, function(x) length(strsplit(x, split = 'gwas')[[1]])>1)]
  
  tmp_fixed <- do.call(rbind, lapply(randomPriorMat, function(x) x[,fixed_prior[[i]]$prior]!=0))
  tmp_fixed_all <- apply(tmp_fixed, 1, any)
  tmp_fixed_gwas <- tmp_fixed[,gwas_name]
  
  for(j in 1:nrep_rG){
    print(j)
    tmp_rG <-  do.call(rbind, lapply(randomPriorMat, function(x) x[,rG_prior[[i]][,j]]!=0))
    df_int[[i]][,j] <- apply(tmp_rG, 2, function(x) length(which(x & tmp_fixed_all)))
    df_int_g[[i]][,j] <- apply(tmp_rG, 2, function(x) length(which(x & tmp_fixed_gwas)))
    
  }
  
  perc_int[[i]] <- df_int[[i]]/nprior_rG[[i]]
  perc_int_g[[i]] <- df_int_g[[i]]/nprior_rG[[i]]
  
}

# for each repetiton load final results
df_res <- list()
for(i in 1:nrow(df_info)){
  
  tmp_file <- read.table(sprintf('%s/%s/rep1/%s', df_info$tissue[i], df_info$fold_gwas[i], df_info$file_prior[i]), stringsAsFactors = F, h=F)
  names_p <- as.vector(sapply(tmp_file$V2, function(x) strsplit(x, split = '_r1')[[1]][1]))
  tmp <- matrix(ncol = nrep_rG, nrow = length(names_p))
  
  for(j in 1:nrep_rG){
    print(j)
    res_w <- get(load(sprintf('%s/%s/rep%i/resPrior_EOpt_HeritableGenes_allchr.RData', df_info$tissue[i], df_info$fold_gwas[i], j)))
    
    tmp[,j] <- res_w$weights
  }
  
  df_res[[i]] <- data.frame(prior = names_p, stringsAsFactors = F)
  df_res[[i]]$mean_w <- rowMeans(tmp)*conv_par
  df_res[[i]]$sd_w <- apply(tmp, 1, sd)
  
}

# weights iteration
df_it <- list()
for(i in 1:nrow(df_info)){
  
  tmp_file <- read.table(sprintf('%s/%s/rep1/%s', df_info$tissue[i], df_info$fold_gwas[i], df_info$file_prior[i]), stringsAsFactors = F, h=F)
  names_p <- as.vector(sapply(tmp_file$V2, function(x) strsplit(x, split = '_r1')[[1]][1]))
  tmp <- array(0,dim=c(20,length(names_p),nrep_rG))
  
  
  for(j in 1:nrep_rG){
    print(j)
    res_it <- get(load(sprintf('%s/%s/rep%i/resPrior_EOpt_Iteration_HeritableGenes_allchr.RData', df_info$tissue[i], df_info$fold_gwas[i], j)))
    tmp[,,j] <- as.matrix(res_it$pWeight, ncol = length(names_p))
  }
  
  df_it[[i]] <- data.frame(mean_w = as.vector(apply(tmp, c(1,2), mean))*conv_par, sd_w = as.vector(apply(tmp, c(1,2), sd)), stringsAsFactors = F)
  df_it[[i]]$it <- rep(1:20, length(names_p))
  df_it[[i]]$prior <- as.vector(sapply(names_p, function(x) rep(x, 20)))
  
}

# build data frame
# f <- function(pal) brewer.pal(brewer.pal.info[pal, "maxcolors"], pal)
for(i in 1:nrow(df_info)){
  
  df_rG_n <- data.frame(prior = c(fixed_prior[[i]]$prior, as.vector(sapply(rG_prior[[i]][,1], function(x) strsplit(x, split = '_r1')[[1]][1]))), stringsAsFactors = F)
  df_rG_n$nsnps <- c(nprior_fixed[[i]], rowMeans(nprior_rG[[i]]))
  df_rG_n$sd <- c(rep(0, length(nprior_fixed[[i]])), apply(nprior_rG[[i]], 1, sd))
  df_rG_n$prior <- factor(df_rG_n$prior, levels = rev(df_rG_n$prior))
  
  df_perc <- data.frame(prior = rep(as.vector(sapply(rG_prior[[i]][,1], function(x) strsplit(x, split = '_r1')[[1]][1])), 2), stringsAsFactors = F)
  df_perc$perc <- c(rowMeans(perc_int[[i]]), rowMeans(1-perc_int[[i]]))
  df_perc$sd <- c(apply(perc_int[[i]], 1, sd), apply(1-perc_int[[i]], 1, sd))
  df_perc$type <- c(rep('shared', nrow(perc_int[[i]])),rep('unique', nrow(perc_int[[i]])))  
  df_perc$type <- factor(df_perc$type, levels = c( 'shared',  'unique'))
  df_perc$prior <- factor( df_perc$prior, levels =  rev(unique(df_perc$prior)))
  
  if(df_info$tissue[i] == 'Artery_Coronary'){
    col <- rev(c(rep('grey', 6), rep('#FFA500', 4)))
    col_per <- c('#ffda96', '#FFA500')
    # shape_it <- c(rep(16, 7), 15,17)
    height_pl <- 3
  }
  if(df_info$tissue[i] == 'Brain_Cortex'){
    col <- rev(c(rep('grey', 14), rep('#8A2BE2', 4)))
    col_per <- c('#b791db', '#8A2BE2')
    height_pl <- 4.5
    # shape_it <- c(rep(16, 15), 15,17)
  }
  
  
  pl_nsnps <- ggplot(df_rG_n, aes(x = prior, y=nsnps, fill = prior)) + 
    geom_bar(stat="identity", color="black") +ggtitle('Random GWAS')+
    ylab('n. SNPs with prior')+ xlab('')+ theme_classic()+
    geom_errorbar(aes(ymin=nsnps-sd, ymax=nsnps+sd), width=.2, position=position_dodge(.9))+
    theme(plot.title = element_text(hjust = 0.5), text = element_text(size = 10),
          axis.text.x=element_text(size = 10, angle = 0, hjust = 1),
          axis.text.y=element_text(size = 10), legend.position = 'none')+
    scale_fill_manual(values = col)+coord_flip()
  ggsave(filename = sprintf('%snSNPswithPrior_rG_%s.png',outFold, df_info$tissue[i]), plot = pl_nsnps, width = 6, height = height_pl, dpi = 500)
  ggsave(filename = sprintf('%snSNPswithPrior_rG_%s.pdf',outFold, df_info$tissue[i]), plot = pl_nsnps, width = 6, height = height_pl)
  
  pl_perc <- ggplot(df_perc, aes(x = prior, y=perc, fill = type)) + 
    geom_bar(stat="identity", color="black") +ggtitle('Random GWAS')+
    ylab('mean percentage')+ xlab('')+ theme_classic()+ylim(0,1)+
    # geom_errorbar(aes(ymin=perc-sd, ymax=perc+sd), width=.2)+
    theme(plot.title = element_text(hjust = 0.5), text = element_text(size = 10),
          axis.text.x=element_text(size = 10, angle = 0, hjust = 1),
          axis.text.y=element_text(size = 10), legend.position = 'bottom')+
    scale_fill_manual(values = col_per)+
    labs(fill = " ") + coord_flip()
  ggsave(filename = sprintf('%spercShared_SNPswithPrior_rG_%s.png', outFold, df_info$tissue[i]), plot = pl_perc, width = 6, height = 2.1, dpi = 500)
  ggsave(filename = sprintf('%spercShared_SNPswithPrior_rG_%s.pdf', outFold, df_info$tissue[i]), plot = pl_perc, width = 6, height = 2.1)
  
  df_res[[i]]$prior <- factor(df_res[[i]]$prior, levels = rev(df_rG_n$prior)) 
  pl_w <- ggplot(df_res[[i]], aes(x = prior, y=mean_w, fill = prior)) + 
    geom_bar(stat="identity", color="black") +ggtitle('Random GWAS')+
    ylab('Maximum contribution to prior coefficient')+ xlab('')+ theme_classic()+
    geom_errorbar(aes(ymin=mean_w-sd_w, ymax=mean_w+sd_w), width=.2, position=position_dodge(.9))+
    theme(plot.title = element_text(hjust = 0.5), text = element_text(size = 10),
          axis.text.x=element_text(size = 10, angle = 0, hjust = 1),
          axis.text.y=element_text(size = 10), legend.position = 'none')+
    scale_fill_manual(values = col)+coord_flip()
  ggsave(filename = sprintf('%smaxContrPrior_rG_%s.png', outFold, df_info$tissue[i]), plot = pl_w, width = 6, height = height_pl, dpi = 500)
  ggsave(filename = sprintf('%smaxContrPrior_rG_%s.pdf', outFold, df_info$tissue[i]), plot = pl_w, width = 6, height = height_pl)
  
  df_it[[i]]$prior <- factor(df_it[[i]]$prior, levels = df_rG_n$prior)
  df_it[[i]]$text <- ""
  id <- df_it[[i]]$prior %in% df_rG_n$prior[rev(col)!='grey'] & df_it[[i]]$it == 1
  df_it[[i]]$text[id] <- as.character(df_it[[i]]$prior)[id]
  
  pl_it <- ggplot(df_it[[i]], aes(x = it, y = mean_w, color = prior, label = text)) + 
    geom_point(size=0.7, alpha=0.8) + geom_line(alpha = 0.8)+ ggtitle('Random GWAS')+
    geom_errorbar(aes(ymin=mean_w-sd_w, ymax=mean_w+sd_w), width=.2, position=position_dodge(0.05))+
    geom_text_repel(size = 2.5, segment.size = 0.2, segment.color = 'black', force = 10, nudge_x = 10, aes(fontface = 2))+
    theme_classic() + xlab('iteration')+  
    ylab('Maximum contribution to prior coefficient') + theme(legend.position = 'none', plot.title = element_text(hjust = 0.5))+
    scale_color_manual(values = rev(col))+labs(color = "")
  ggsave(filename = sprintf('%smaxContrPrior_iteration_rG_%s.png', outFold, df_info$tissue[i]), plot = pl_it, width = 4, height = 4, dpi=500)
  ggsave(filename = sprintf('%smaxContrPrior_iteration_rG_%s.pdf', outFold, df_info$tissue[i]), plot = pl_it, width = 4, height = 4)
  
  res <- list(weights = df_res[[i]], nsnps = df_rG_n, it = df_it[[i]], perc_shared = df_perc)
  save(res, file = sprintf('%stotRes_rG_%s.RData', outFold, df_info$tissue[i]))
  
}


