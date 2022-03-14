# plot prediction cluster

options(stringsAsFactors=F)
options(max.print=1000)
suppressPackageStartupMessages(library(argparse))
suppressPackageStartupMessages(library(data.table))
suppressPackageStartupMessages(library(doParallel))
suppressPackageStartupMessages(library(qvalue))
suppressPackageStartupMessages(library(pROC))
suppressPackageStartupMessages(library(pryr))
suppressPackageStartupMessages(library(umap))
suppressPackageStartupMessages(library(igraph))
suppressPackageStartupMessages(library(Matrix))
suppressPackageStartupMessages(library(SparseM))
suppressPackageStartupMessages(library(RColorBrewer))
suppressPackageStartupMessages(library(pheatmap))
suppressPackageStartupMessages(library(ggplot2))
suppressPackageStartupMessages(library(MASS))
suppressPackageStartupMessages(library(lmtest))
suppressPackageStartupMessages(library(ggsci))
options(bitmapType = 'cairo', device = 'png')

parser <- ArgumentParser(description="predict cluster probability for new samples")
parser$add_argument("--cohort_name", nargs = '*', type = "character", help = "")
parser$add_argument("--model_name", type = "character", help = "")
parser$add_argument("--phenoNew_file", nargs = '*', type = "character", help = "")
parser$add_argument("--type_cluster", type = "character",default = 'All', help = "All, Cases, Controls")
parser$add_argument("--type_data", type = "character", help = "pathway or tscore")
parser$add_argument("--clustFile", type = "character", help = "file cluster results")
parser$add_argument("--featRel_model", type = "character", default = NULL, help = "file association features wilcoxon test")
parser$add_argument("--clustFile_new", type = "character", nargs = '*', help = "file cluster results")
parser$add_argument("--featRel_predict", type = "character", default = NULL, nargs = '*', help = "file association features wilcoxon test")
parser$add_argument("--functR", type = "character", help = "functions to be used")
parser$add_argument("--type_input", type = "character", default = 'original', help = "original or zscaled")
parser$add_argument("--tissues_name", type = "character", help = "name tissue")
parser$add_argument("--geneLoci_summ", type = "character", default = NULL, help = "file with group summary divided per loci")
parser$add_argument("--outFold", type="character", help = "Output file [basename only]")

args <- parser$parse_args()
cohort_name <- args$cohort_name
model_name <- args$model_name
clustFile <- args$clustFile
clustFile_new <- args$clustFile_new
functR <- args$functR
type_data <- args$type_data
type_input <- args$type_input
type_cluster <- args$type_cluster
tissues_name <- args$tissues_name
phenoNew_file <- args$phenoNew_file
featRel_predict <- args$featRel_predict
featRel_model <- args$featRel_model
geneLoci_summ <- args$geneLoci_summ
outFold <- args$outFold

###################################################################################################################
# functR <- '/psycl/g/mpsziller/lucia/castom-igex/Software/model_clustering/clustering_functions.R'
# cohort_name = c(paste0('German', 1:5), c('MG', 'WTCCC', 'LURIC', 'CG'))
# type_data <- 'tscore_corrPCs'
# type_input <- 'zscaled'
# type_cluster <- 'Cases'
# clustFile_new <- sprintf('OUTPUT_GTEx/predict_CAD/Liver/200kb/CAD_GWAS_bin5e-2/%s/devgeno0.01_testdevgeno0/CAD_HARD_clustering/update_corrPCs/tscore_corrPCs_zscaled_predictClusterCases_PGmethod_HKmetric.RData', cohort_name)
# clustFile <- 'OUTPUT_GTEx/predict_CAD/Liver/200kb/CAD_GWAS_bin5e-2/UKBB/devgeno0.01_testdevgeno0/CAD_HARD_clustering/update_corrPCs/tscore_corrPCs_zscaled_clusterCases_PGmethod_HKmetric.RData'
# tissues_name <- 'Liver'
# phenoNew_file <-  ''
# outFold <- ''
# featRel_model <- 'OUTPUT_GTEx/predict_CAD/Liver/200kb/CAD_GWAS_bin5e-2/UKBB/devgeno0.01_testdevgeno0/CAD_HARD_clustering/update_corrPCs/tscoreOriginal_corrPCs_tscoreClusterCases_featAssociation.RData'
# featRel_predict <- sprintf('OUTPUT_GTEx/predict_CAD/Liver/200kb/CAD_GWAS_bin5e-2/%s/devgeno0.01_testdevgeno0/CAD_HARD_clustering/update_corrPCs/tscoreOriginal_corrPCs_tscoreClusterCases_featAssociation.RData', cohort_name)
# geneLoci_summ <- 'OUTPUT_GTEx/predict_CAD/Liver/200kb/CAD_GWAS_bin5e-2/UKBB/devgeno0.01_testdevgeno0/CAD_HARD_clustering/update_corrPCs/tscore_corrPCs_zscaled_clusterCases_summary_geneLoci_allTissues.txt'
# #################################################################################################################

source(functR)

tmp <- get(load(clustFile)) 
clust <- tmp$cl_best
P <- length(unique(clust$gr))
cl_name <- sort(unique(clust$gr))

df <- data.frame(dataset = rep(model_name, P), type = rep('model', P), gr = paste0('gr_', sort(unique(clust$gr))))
df$n <- sapply(sort(unique(clust$gr)), function(x) length(which(clust$gr == x)))
df$percentage <- sapply(sort(unique(clust$gr)), function(x) length(which(clust$gr == x))/nrow(clust))
data <- tmp$input_data
mean_gr <- tmp$gr_input$mean

# load associated features
if(!is.null(featRel_model)){
  tmp <- get(load(featRel_model))
  featRel <- tmp$test_feat
  featRel <- do.call(rbind, featRel)
  featRel <- featRel[featRel$pval_corr <= 0.05, ]
  featRel$new_id <- paste(featRel$feat, featRel$comp, featRel$tissue, sep = '_')
}

if(!is.null(geneLoci_summ)){
  geneLoci <- read.delim(geneLoci_summ, h=T, stringsAsFactors = F, sep = '\t') 
}

sampleAnn_new <- list()
phenoDat_new <- list()
clust_new <- list()
data_new <- list()
df_new <- list()
mean_gr_new <- list()
df_corr <- list()
df_corr_rel <- list()
df_perc_loci <- list()

for(i in 1:length(cohort_name)){
  print(cohort_name[i])
  tmp <- get(load(clustFile_new[i])) 
  
  sampleAnn_new[[i]] <- tmp$sampleAnn

  clust_new[[i]] <- tmp$cl_new
  data_new[[i]] <- tmp$data_new
  mean_gr_new[[i]] <- tmp$gr_input$mean
  
  df_new[[i]] <- data.frame(dataset = rep(cohort_name[i], P), type = rep('predict', P), gr = df$gr)
  df_new[[i]]$n <- sapply(sort(unique(clust$gr)), function(x) length(which(clust_new[[i]]$gr == x)))
  df_new[[i]]$percentage <- sapply(sort(unique(clust$gr)), function(x) length(which(clust_new[[i]]$gr == x))/nrow(clust_new[[i]]))
  
  df_corr[[i]] <- data.frame(dataset = rep(cohort_name[i], P), 
                             gr = df$gr, corr = rep(NA, P), 
                             pvalue = rep(NA, P), CI_low= rep(NA, P),  CI_up= rep(NA, P))
  df_corr[[i]] <- df_corr[[i]][df_corr[[i]]$gr %in% colnames(mean_gr_new[[i]]), ]
  df_corr[[i]]$corr <- sapply(df$gr[df$gr %in% colnames(mean_gr_new[[i]])], 
                              function(x) cor.test(mean_gr_new[[i]][,x], mean_gr[,x])$estimate)
  df_corr[[i]]$pvalue <- sapply(df$gr[df$gr %in% colnames(mean_gr_new[[i]])], 
                                function(x) cor.test(mean_gr_new[[i]][,x], mean_gr[,x])$p.value)
  df_corr[[i]]$CI_low <- sapply(df$gr[df$gr %in% colnames(mean_gr_new[[i]])], 
                                function(x) cor.test(mean_gr_new[[i]][,x], mean_gr[,x])$conf.int[1])
  df_corr[[i]]$CI_up <- sapply(df$gr[df$gr %in% colnames(mean_gr_new[[i]])], 
                               function(x) cor.test(mean_gr_new[[i]][,x], mean_gr[,x])$conf.int[2])
  
  if(!is.null(featRel_predict)){
    
    comp <- sort(unique(featRel$comp))
    tmp <- get(load(featRel_predict[i]))
    featRel_new <- tmp$test_feat
    featRel_new <- do.call(rbind, featRel_new)
    featRel_new$new_id <- paste(featRel_new$feat, featRel_new$comp, featRel_new$tissue, sep = '_')
    common_f <- intersect(featRel_new$new_id, featRel$new_id)
    featRel_new <- featRel_new[match(common_f, featRel_new$new_id), ]
    tmp <- featRel[match(common_f, featRel$new_id),]
    comp_new <- sort(unique(featRel_new$comp))
    # compute spearman correlation
    df_corr_rel[[i]] <- data.frame(dataset = rep(cohort_name[i], P), 
                                   gr = df$gr, corr = rep(NA, P), pvalue = rep(NA, P))
    df_corr_rel[[i]] <- df_corr_rel[[i]][df_corr_rel[[i]]$gr %in% colnames(mean_gr_new[[i]]), ]
    df_corr_rel[[i]]$corr <- sapply(comp[comp %in% comp_new], 
                                    function(x) cor.test(featRel_new$estimates[featRel_new$comp == x], tmp$estimates[tmp$comp == x], method = 'spearman')$estimate)
    df_corr_rel[[i]]$pvalue <- sapply(comp[comp %in% comp_new], 
                                      function(x) cor.test(featRel_new$estimates[featRel_new$comp == x], tmp$estimates[tmp$comp == x], method = 'spearman')$p.value)
  }
  
  if(!is.null(geneLoci_summ)){
    
    id <- sapply(df$gr, function(x) paste0(strsplit(x, split = '_')[[1]], collapse = ''))
    df_perc_loci[[i]] <- data.frame(dataset = rep(cohort_name[i], P), gr = df$gr)
    df_perc_loci[[i]]$nloci <- sapply(id, function(x) sum(grepl(x, geneLoci$comp_sign)))
    
    for(j in 1:nrow(df)){
      tmp <- geneLoci[grepl(id[j], geneLoci$comp_sign), ]
      genes <- lapply(tmp$gene, function(x) strsplit(x, split = ',')[[1]])
      id_keep <- sapply(genes, function(x) which.min(featRel$pval[featRel$comp == comp[j] & 
                                                                    featRel$feat %in% x]))
      genes <- mapply(function(x, y) featRel$new_id[featRel$comp == comp[j] & featRel$feat %in% x][y], 
                      x = genes, y = id_keep)
      tmp_mod <- featRel[match(genes, featRel$new_id), ]
      tmp_pred <- featRel_new[match(genes, featRel_new$new_id), ]
      df_perc_loci[[i]]$nloci_rep_sign[j] <- sum(sign(tmp_pred$estimates) == sign(tmp_mod$estimates))
      df_perc_loci[[i]]$nloci_rep[j] <- sum(sign(tmp_pred$estimates) == sign(tmp_mod$estimates) 
                                            & tmp_pred$pval <= 0.05)
      
    }
    df_perc_loci[[i]]$nperc_loci_rep_sign <- df_perc_loci[[i]]$nloci_rep_sign/df_perc_loci[[i]]$nloci
  }
  
}

df_tot <- rbind(df, do.call(rbind, df_new))
df_corr_tot <- do.call(rbind, df_corr)
if(!is.null(geneLoci_summ)){
  df_perc_loci <- do.call(rbind, df_perc_loci)
}
if(!is.null(featRel_model)){
  df_corr_rel <- do.call(rbind,  df_corr_rel)
}


# save and plot
write.table(df_tot, file = sprintf('%s%s_%s_cluster%s_percentageGropus_prediction_model%s.txt', outFold, type_data, type_input, type_cluster, model_name),quote = F, 
            col.names = T, row.names = T, sep = '\t')

write.table(df_corr_tot, file = sprintf('%s%s_%s_cluster%s_correlationMeanGroups_prediction_model%s.txt', outFold, type_data, type_input, type_cluster, model_name),quote = F, 
            col.names = T, row.names = T, sep = '\t')

if(!is.null(featRel_model)){
  write.table(df_corr_rel, file = sprintf('%s%s_%s_cluster%s_correlationSpear_WMWestSign_Groups_prediction_model%s.txt', outFold, type_data, type_input, type_cluster, model_name),quote = F, 
              col.names = T, row.names = T, sep = '\t')
}
if(!is.null(geneLoci_summ)){
  write.table(df_perc_loci, file = sprintf('%s%s_%s_cluster%s_numberLociRep_prediction_model%s.txt', outFold, type_data, type_input, type_cluster, model_name),quote = F, 
              col.names = T, row.names = T, sep = '\t')
}

###
df_tot$new_id <- paste0(df_tot$dataset, '\n(', df_tot$type, ')')
df_tot$new_id <- factor(df_tot$new_id, levels = c(paste0(model_name, '\n(model)'),  paste0(cohort_name, '\n(predict)')))
df_tot$gr <- factor(df_tot$gr, levels = paste0('gr_', sort(unique(clust$gr))))
df$gr <- factor(df$gr, levels = paste0('gr_', sort(unique(clust$gr))))

gr_color <- pal_d3(palette = 'category20')(P)

pl <- ggplot(df_tot, aes(x = new_id, y = percentage, color = gr, group = gr))+
  geom_point(size = 2, ,position = position_dodge(width = 0.3))+
  theme_bw()+ 
  geom_segment(data = df, aes(x = 2, y = percentage, xend = length(cohort_name)+1, yend = percentage, group = gr), linetype = 2, alpha = 0.6)+
  ylab('Fraction of Cases')+ 
  theme(legend.position = 'right', axis.title.x = element_blank(), 
        axis.text.x = element_text(angle = 45, hjust = 1))+
  scale_color_manual(values = gr_color)
# scale_shape_manual(values=c(1, 19))+
w=ifelse(length(cohort_name)==1, 3, 1+length(cohort_name)*0.6)

ggsave(filename = sprintf('%s%s_%s_cluster%s_percentageGropus_prediction_model%s.png', outFold, type_data, type_input, type_cluster, model_name), width = w, height = 3.5, plot = pl, device = 'png')
ggsave(filename = sprintf('%s%s_%s_cluster%s_percentageGropus_prediction_model%s.pdf', outFold, type_data, type_input, type_cluster, model_name), width = w, height = 3.5, plot = pl, device = 'pdf')

###
df_corr_tot$dataset <-factor(df_corr_tot$dataset, levels = cohort_name)
df_corr_tot$gr <- factor(df_corr_tot$gr, levels = paste0('gr_', sort(unique(clust$gr))))

pl <- ggplot(df_corr_tot, aes(x = dataset, y = corr, fill = gr, group = gr))+
  geom_bar(stat = 'identity',width = 0.7, color = 'black', alpha = 0.7, position = position_dodge())+
  geom_errorbar(aes(ymin=CI_low, ymax=CI_up), width=.2, position=position_dodge(.75))+
  theme_bw()+ 
  coord_cartesian(ylim = c(ifelse(min(df_corr_tot$corr)<0.6, 0, 0.6),1)) +
  ylab(sprintf('correlation mean scores\nwith %s (model)', model_name))+ 
  theme(legend.position = 'right', axis.title.x = element_blank(), 
  axis.text.x = element_text(angle = 45, hjust = 1))+
  scale_fill_manual(values = gr_color)
# scale_shape_manual(values=c(1, 19))+
ggsave(filename = sprintf('%s%s_%s_cluster%s_correlationMeanGroups_prediction_model%s.png', outFold, type_data, type_input, type_cluster, model_name), width = w, height = 3.5, plot = pl, device = 'png')
ggsave(filename = sprintf('%s%s_%s_cluster%s_correlationMeanGroups_prediction_model%s.pdf', outFold, type_data, type_input, type_cluster, model_name), width = w, height = 3.5, plot = pl, device = 'pdf')

###
if(!is.null(featRel_model)){
  
  df_corr_rel$dataset <- factor(df_corr_rel$dataset, levels = cohort_name)
  df_corr_rel$gr <- factor(df_corr_rel$gr, levels = paste0('gr_', sort(unique(clust$gr))))
  
  pl <- ggplot(df_corr_rel, aes(x = dataset, y = corr, fill = gr, group = gr))+
    geom_bar(stat = 'identity',width = 0.7, color = 'black', alpha = 0.7, position = position_dodge())+
    # geom_errorbar(aes(ymin=CI_low, ymax=CI_up), width=.2, position=position_dodge(.75))+
    theme_bw()+ 
    coord_cartesian(ylim = c(ifelse(min(df_corr_rel$corr)<0.6, 0, 0.6),1)) +
    ylab(sprintf('Spearman corr. from WMW estimates\nwith %s (model)', model_name))+ 
    theme(legend.position = 'right', axis.title.x = element_blank(),
          axis.text.x = element_text(angle = 45, hjust = 1))+
    scale_fill_manual(values = gr_color)
  # scale_shape_manual(values=c(1, 19))+
  ggsave(filename = sprintf('%s%s_%s_cluster%s_correlationSpear_WMWestSign_Groups_prediction_model%s.png', outFold, type_data, type_input, type_cluster, model_name), width = w, height = 3.5, plot = pl, device = 'png')
  ggsave(filename = sprintf('%s%s_%s_cluster%s_correlationSpear_WMWestSign_Groups_prediction_model%s.pdf', outFold, type_data, type_input, type_cluster, model_name), width = w, height = 3.5, plot = pl, device = 'pdf')
  
}

###
if(!is.null(geneLoci_summ)){
  
  df_perc_loci$dataset <- factor(df_perc_loci$dataset, levels = cohort_name)
  df_perc_loci$gr <- factor(df_perc_loci$gr, levels = paste0('gr_', sort(unique(clust$gr))))
  
  if(length(cohort_name)>1){  
    pl <- ggplot(df_perc_loci,aes(x = nloci, y = nloci_rep, color = gr))+
      geom_point(alpha = 0.7, size = 2)+
      facet_wrap(.~dataset, nrow = 1)+
      theme_bw()+ 
      ylab('n. of loci reproduced')+ xlab('n. of loci')+
      theme(legend.position = 'right')+
      scale_color_manual(values = gr_color)
    # scale_shape_manual(values=c(1, 19))+
    width_plot=1.5*length(cohort_name)
  }else{
    pl <- ggplot(df_perc_loci,aes(x = nloci, y = nloci_rep, color = gr))+
      geom_point(alpha = 0.7, size = 2)+
      theme_bw()+
      ylab('n. of loci reproduced')+ xlab('n. of loci')+
      theme(legend.position = 'right')+
      scale_color_manual(values = gr_color)
    width_plot=3
  }
  
  ggsave(filename = sprintf('%s%s_%s_cluster%s_numberLociRep_Groups_prediction_model%s.png', outFold, type_data, type_input, type_cluster, model_name), width = width_plot, height = 2, plot = pl, device = 'png')
  ggsave(filename = sprintf('%s%s_%s_cluster%s_numberLociRep_Groups_prediction_model%s.pdf', outFold, type_data, type_input, type_cluster, model_name), width = width_plot, height = 2, plot = pl, device = 'pdf')
  
}


#### endophenotype association ####
################
## gri vs grj ##
################

if(any(sapply(phenoNew_file, file.exists))){
  suppressPackageStartupMessages(library(RNOmni))
}

phenoInfo_new <- list()
tot_bin_reg <- list()

for(i in 1:length(cohort_name)){
  
  if(file.exists(phenoNew_file[i])){        
    
    phenoDat_new[[i]] <- read.table(phenoNew_file[i], h=T, stringsAsFactors = F)
    
    print(all(phenoDat_new[[i]]$Individual_ID %in% sampleAnn_new[[i]]$Individual_ID))
    phenoDat_new[[i]] <- phenoDat_new[[i]][match(sampleAnn_new[[i]]$Individual_ID, phenoDat_new[[i]]$Individual_ID),]
    
    sampleAnn_new[[i]] <- sampleAnn_new[[i]][, colnames(sampleAnn_new[[i]]) %in% c('Individual_ID', paste0('C', 1:10))]
    sampleAnn_new[[i]]$Age <- phenoDat_new[[i]]$Age
    sampleAnn_new[[i]]$Gender <- phenoDat_new[[i]]$Gender
    
    P <- length(unique(clust_new[[i]]$gr))
    gr_names <- sort(unique(clust_new[[i]]$gr))
    cl <- clust_new[[i]]$gr
    
    covDat <- sampleAnn_new[[i]][, !colnames(sampleAnn_new[[i]]) %in% c('Individual_ID', 'genoSample_ID', 'Dx')]
    fmla  <- as.formula(paste('pheno~gr_id+', paste0(colnames(covDat), collapse = '+')))
    phenoDat <- phenoDat_new[[i]][, !colnames(phenoDat_new[[i]]) %in%  c('Individual_ID', 'Dx', 'Age', 'Gender')]
    
    if(any(table(cl)<=10)){
      rm_id <- names(which(table(cl)<=10))
      P <- P-length(rm_id)
      gr_names <- gr_names[!gr_names %in% rm_id]
      covDat <- covDat[!cl %in% rm_id,]
      phenoDat <- phenoDat[!cl %in% rm_id,]
      cl <- cl[!cl %in% rm_id]
    }
    
    phenoInfo_new[[i]] <- data.frame(pheno_id = colnames(phenoDat))
    phenoInfo_new[[i]]$type_pheno <- 'CONTINUOUS'
    phenoInfo_new[[i]]$type_pheno[sapply(1:ncol(phenoDat), function(x) is.integer(phenoDat[,x]) & length(unique(na.omit(phenoDat[,x]))) == 2)] <- 'CAT_SINGLE_BINARY'
    phenoInfo_new[[i]]$type_pheno[sapply(1:ncol(phenoDat), function(x) is.integer(phenoDat[,x]) & length(unique(na.omit(phenoDat[,x]))) > 2)] <- 'CAT_ORD'
    for(j in 1:ncol(phenoDat)){
      if(phenoInfo_new[[i]]$type_pheno[j] == 'CONTINUOUS'){
        tmp <- phenoDat[!is.na(phenoDat[,j]),j]
        phenoDat[!is.na(phenoDat[,j]),j] <- rankNorm(tmp)
      }
    }
    
    bin_reg <- vector(mode = 'list', length = length(gr_names)-1)
    for(k in 1:(length(gr_names)-1)){
      
      print(paste0('group', gr_names[k], '_vs_groupj'))
      
      # j vs all
      pheno_case_tmp <- lapply(gr_names[k:length(gr_names)], function(x) phenoDat[cl == x,])
      covDat_tmp <- lapply(gr_names[k:length(gr_names)], function(x) covDat[cl == x,])
      bin_reg[[k]] <-  vector(mode = 'list', length = length(pheno_case_tmp)-1)
      
      for(j in 2:length(pheno_case_tmp)){
        
        print(j)
        
        new <- rbind(pheno_case_tmp[[1]], pheno_case_tmp[[j]])
        gr_id <- factor(c(rep(0, nrow(pheno_case_tmp[[1]])), rep(1, nrow(pheno_case_tmp[[j]]))))
        
        # remove pheno with constant values
        p_rm <- apply(new, 2, function(x) length(unique(x)) == 1)
        new <- new[,!p_rm]
        
        new_cov <- rbind(covDat_tmp[[1]], covDat_tmp[[j]])
        res_glm <- matrix(nrow = ncol(new), ncol = 7)
        for(l in 1:ncol(new)){
          type_pheno <- phenoInfo_new[[i]]$type_pheno[phenoInfo_new[[i]]$pheno_id == colnames(new)[l]]
          tmp_dat <- cbind(data.frame(pheno = new[, l], gr_id = gr_id), new_cov)
          res_glm[l,] <- compute_reg_endopheno(mat = tmp_dat, fmla = fmla, type_pheno = type_pheno)
        }
        colnames(res_glm) <- c('beta', 'se_beta', 'z', 'pvalue', 'OR_or_beta', 'CI_low', 'CI_up')
        res_glm <- as.data.frame(res_glm)
        
        phenoInfo_tmp <- phenoInfo_new[[i]][match(colnames(new), phenoInfo_new[[i]]$pheno_id),]
        
        bin_reg[[k]][[j-1]] <- cbind(data.frame(pheno_id = phenoInfo_tmp$pheno_id, type_pheno = phenoInfo_tmp$type_pheno), res_glm)
        bin_reg[[k]][[j-1]]$pval_corr <- p.adjust(bin_reg[[k]][[j-1]]$pvalue, method = 'BH')
        bin_reg[[k]][[j-1]]$comp <- sprintf('gr%i_vs_gr%i',  gr_names[k:length(gr_names)][j], gr_names[k])
        
      }
      
    }
    tot_bin_reg[[i]] <- do.call(rbind, do.call(c,bin_reg))
    tot_bin_reg[[i]]$pval_corr_overall <-  p.adjust(tot_bin_reg[[i]]$pvalue, method = 'BY')
  }  
}

if(any(sapply(phenoNew_file, file.exists))){
  # save results
  output <- list(bin_reg = tot_bin_reg, cl = clust_new, phenoDat = phenoDat_new, phenoInfo = phenoInfo_new)
  save(output, file = sprintf('%s%s_%s_cluster%s_phenoAssociationGLMpairwise_prediction_model%s.RData', outFold, type_data, type_input, type_cluster, model_name))
}

################
## gri vs all ##
################

phenoInfo_new <- list()
tot_bin_reg <- list()

for(i in 1:length(cohort_name)){
  
  if(file.exists(phenoNew_file[i])){        
    
    phenoDat_new[[i]] <- read.table(phenoNew_file[i], h=T, stringsAsFactors = F)
    
    print(all(phenoDat_new[[i]]$Individual_ID %in% sampleAnn_new[[i]]$Individual_ID))
    phenoDat_new[[i]] <- phenoDat_new[[i]][match(sampleAnn_new[[i]]$Individual_ID, phenoDat_new[[i]]$Individual_ID),]
    
    sampleAnn_new[[i]] <- sampleAnn_new[[i]][, colnames(sampleAnn_new[[i]]) %in% c('Individual_ID', paste0('C', 1:10))]
    sampleAnn_new[[i]]$Age <- phenoDat_new[[i]]$Age
    sampleAnn_new[[i]]$Gender <- phenoDat_new[[i]]$Gender
    
    P <- length(unique(clust_new[[i]]$gr))
    gr_names <- sort(unique(clust_new[[i]]$gr))
    cl <- clust_new[[i]]$gr
    
    covDat <- sampleAnn_new[[i]][, !colnames(sampleAnn_new[[i]]) %in% c('Individual_ID', 'genoSample_ID', 'Dx')]
    fmla  <- as.formula(paste('pheno~gr_id+', paste0(colnames(covDat), collapse = '+')))
    phenoDat <- phenoDat_new[[i]][, !colnames(phenoDat_new[[i]]) %in%  c('Individual_ID', 'Dx', 'Age', 'Gender')]
    
    if(any(table(cl)<=10)){
      rm_id <- names(which(table(cl)<=10))
      P <- P-length(rm_id)
      gr_names <- gr_names[!gr_names %in% rm_id]
      covDat <- covDat[!cl %in% rm_id,]
      phenoDat <- phenoDat[!cl %in% rm_id,]
      cl <- cl[!cl %in% rm_id]
    }
    
    phenoInfo_new[[i]] <- data.frame(pheno_id = colnames(phenoDat))
    phenoInfo_new[[i]]$type_pheno <- 'CONTINUOUS'
    phenoInfo_new[[i]]$type_pheno[sapply(1:ncol(phenoDat), function(x) is.integer(phenoDat[,x]) & length(unique(na.omit(phenoDat[,x]))) == 2)] <- 'CAT_SINGLE_BINARY'
    phenoInfo_new[[i]]$type_pheno[sapply(1:ncol(phenoDat), function(x) is.integer(phenoDat[,x]) & length(unique(na.omit(phenoDat[,x]))) > 2)] <- 'CAT_ORD'
    for(j in 1:ncol(phenoDat)){
      if(phenoInfo_new[[i]]$type_pheno[j] == 'CONTINUOUS'){
        tmp <- phenoDat[!is.na(phenoDat[,j]),j]
        phenoDat[!is.na(phenoDat[,j]),j] <- rankNorm(tmp)
      }
    }
    
    
    bin_reg <- vector(mode = 'list', length = length(gr_names))
    for(k in 1:length(gr_names)){
      
      print(paste0('group', gr_names[k], '_vs_all'))
      
      # j vs all
      pheno_case_tmp <- list(phenoDat[cl == gr_names[k],], phenoDat[cl != gr_names[k],])
      covDat_tmp <- list(covDat[cl == gr_names[k],], covDat[cl != gr_names[k],]) 
      
      new <- do.call(rbind, pheno_case_tmp)
      gr_id <- factor(c(rep(1, nrow(pheno_case_tmp[[1]])), rep(0, nrow(pheno_case_tmp[[2]]))))
      
      # remove pheno with constant values
      p_rm <- apply(new, 2, function(x) length(unique(x)) == 1)
      new <- new[,!p_rm]
      
      new_cov <- do.call(rbind, covDat_tmp)
      res_glm <- matrix(nrow = ncol(new), ncol = 7)
      for(l in 1:ncol(new)){
        # print(l)  
        type_pheno <- phenoInfo_new[[i]]$type_pheno[phenoInfo_new[[i]]$pheno_id == colnames(new)[l]]
        tmp_dat <- cbind(data.frame(pheno = new[, l], gr_id = gr_id), new_cov)
        res_glm[l,] <- compute_reg_endopheno(mat = tmp_dat, fmla = fmla, type_pheno = type_pheno)
      }
      
      colnames(res_glm) <- c('beta', 'se_beta', 'z', 'pvalue', 'OR_or_Beta', 'CI_low', 'CI_up')
      res_glm <- as.data.frame(res_glm)
      
      phenoInfo_tmp <- phenoInfo_new[[i]][match(colnames(new), phenoInfo_new[[i]]$pheno_id),]
      bin_reg[[k]] <- cbind(data.frame(pheno_id = phenoInfo_tmp$pheno_id, type_pheno = phenoInfo_tmp$type_pheno), res_glm)
      bin_reg[[k]]$pval_corr <- p.adjust(bin_reg[[k]]$pvalue, method = 'BH')
      bin_reg[[k]]$comp <- sprintf('gr%i_vs_all', gr_names[k])
      
    }
    
    tot_bin_reg[[i]] <- do.call(rbind, bin_reg)
    tot_bin_reg[[i]]$pval_corr_overall <-  p.adjust(tot_bin_reg[[i]]$pvalue, method = 'BY')
  }
}

if(any(sapply(phenoNew_file, file.exists))){
  output <- list(bin_reg = tot_bin_reg, cl = clust_new, phenoDat = phenoDat_new, phenoInfo = phenoInfo_new)
  save(output, file = sprintf('%s%s_%s_cluster%s_phenoAssociationGLM_prediction_model%s.RData', outFold, type_data, type_input, type_cluster, model_name))
}


