# plot specific endophenotype associaition (original data)
# add nominal pvalues

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
suppressPackageStartupMessages(library(ggpubr))
suppressPackageStartupMessages(library(ggplot2))
suppressPackageStartupMessages(library(MASS))
suppressPackageStartupMessages(library(lmtest))
suppressPackageStartupMessages(library(ggsci))
options(bitmapType = 'cairo', device = 'png')

parser <- ArgumentParser(description="plot nominal association")
parser$add_argument("--pheno_plot", nargs = '*', type = "character", help = "")
parser$add_argument("--type_cluster", type = "character",default = 'All', help = "All, Cases, Controls")
parser$add_argument("--type_data", type = "character", help = "pathway or tscore")
parser$add_argument("--phenoRegFile", type = "character", help = "")
parser$add_argument("--functR", type = "character", help = "functions to be used")
parser$add_argument("--type_input", type = "character", default = 'original', help = "original or zscaled")
parser$add_argument("--tissues_name", type = "character", help = "name tissue")
parser$add_argument("--outFold", type="character", help = "Output file [basename only]")

args <- parser$parse_args()
cohort_name_id <- args$cohort_name_id
phenoRegFile <- args$phenoRegFile
functR <- args$functR
type_data <- args$type_data
type_input <- args$type_input
type_cluster <- args$type_cluster
tissues_name <- args$tissues_name
pheno_plot <- args$pheno_plot
outFold <- args$outFold

###################################################################################################################
# functR <- '/psycl/g/mpsziller/lucia/priler_project/Software/model_clustering/clustering_functions.R'
# type_data <- 'tscore'
# type_input <- 'zscaled'
# type_cluster <- 'Cases'
# phenoRegFile <- paste0('/psycl/g/mpsziller/lucia/CAD_UKBB/eQTL_PROJECT/OUTPUT_GTEx/predict_CAD/Liver/200kb/CAD_GWAS_bin5e-2/UKBB/devgeno0.01_testdevgeno0/CAD_HARD_clustering/nominalAnalysis_tscore_zscaled_clusterCases_PGmethod_HKmetric_phenoAssociation_GLM.RData')
# # clustFile <- paste0('/psycl/g/mpsziller/lucia/CAD_UKBB/eQTL_PROJECT/OUTPUT_GTEx/predict_CAD/Liver/200kb/CAD_GWAS_bin5e-2/',cohort_name,'/devgeno0.01_testdevgeno0/CAD_HARD_clustering/tscore_zscaled_predictClusterCases_PGmethod_HKmetric.RData')
# # tissues_name <- 'Adipose_Visceral_Omentum'
# # pheno_plot <- c('History_bleeding', 'Smoking', 'Age_stroke', 'Hypertension', 'UAP')
# # tissues_name <- 'Colon_Sigmoid'
# # pheno_plot <- c('History_bleeding', 'Atherosclerotic_heart_disease', 'Coronary_artery_bypass_graft', 'Age_stroke)
# tissues_name <- 'Liver'
# pheno_plot <- c('Hyperlipidemia', 'Coronary_artery_bypass_graft', 'Age_stroke', 'Chronic_obstructive_pulmonary_disease', 'Peripheral_vascular_disease', 'Atherosclerotic_heart_disease')
# # phenoNew_file <-  paste0('/psycl/g/mpsziller/lucia/CAD_UKBB/eQTL_PROJECT/INPUT_DATA_GTEx/CAD/Covariates/',cohort_name,'/phenotypeMatrix_CADrel_Cases.txt')
# outFold <- '/psycl/g/mpsziller/lucia/CAD_UKBB/eQTL_PROJECT/OUTPUT_GTEx/predict_CAD/Liver/200kb/CAD_GWAS_bin5e-2/UKBB/devgeno0.01_testdevgeno0/CAD_HARD_clustering/nominalAnalysis_'
# #################################################################################################################

source(functR)

# load results
tmp <- get(load(phenoRegFile))

reg_res <- tmp$bin_reg
phenoDat <- tmp$phenoDat
cl <- tmp$cl
P <- length(unique(cl$gr))
gr_id <- sort(unique(cl$gr))
gr_color <- pal_d3(palette = 'category20')(P)

# phenoDat <- phenoDat[match(cl$id,phenoDat$Individual_ID),]
phenoDat <- phenoDat[, colnames(phenoDat) %in% pheno_plot]

phenoDat <- cbind(data.frame(gr = paste0('gr', cl$gr), stringsAsFactors = F), phenoDat)

list_pl <- list()

for(i in 1:length(pheno_plot)){
  
  tmp_reg <- reg_res[reg_res$pheno_id == pheno_plot[i],]
  tmp_reg <- tmp_reg[match(paste0('gr',gr_id, '_vs_all'),tmp_reg$comp),]
  tmp_reg$star_ann <- 'ns'
  tmp_reg$star_ann[ tmp_reg$pvalue <= 0.05 & tmp_reg$pvalue > 0.01] <- '*'
  tmp_reg$star_ann[ tmp_reg$pvalue <= 0.01 & tmp_reg$pvalue > 0.001] <- '**'
  tmp_reg$star_ann[ tmp_reg$pvalue <= 0.001 & tmp_reg$pvalue > 0.0001] <- '***'
  tmp_reg$star_ann[ tmp_reg$pvalue <= 0.00001] <- '****'
  
  new <- phenoDat[!is.na(phenoDat[, pheno_plot[i]]), ]

  if(is.integer(new[, pheno_plot[i]]) & length(unique(new[, pheno_plot[i]])) == 2){
    
    print(i)
    
    val <- max(unique(new[, pheno_plot[i]]))
    new_df <- data.frame(perc = as.vector(sapply(gr_id,  function(y) sum(new[, pheno_plot[i]] == val & new$gr == paste0('gr',y))/sum(new$gr == paste0('gr',y)))))
    new_df$gr <- paste0('gr',gr_id)
    new_df$gr <- factor(new_df$gr, levels = paste0('gr', gr_id))
    
    list_pl[[i]] <- ggplot(new_df, aes(x = gr, y = perc, fill = gr))+
      geom_bar(size = 1, alpha = 0.8, stat = 'identity', color = 'black', width = 0.7)+
      xlab('')+ ylab(paste0('Percentage'))+ 
      scale_fill_manual(values = gr_color)+
      guides(fill = FALSE)+
      ggtitle(pheno_plot[i])+
      # annotate("text", x = 1:P, y = rep(max(new_df$perc)+0.02, P), label =formatC(tmp_reg$pvalue,format="e", digits = 2), size = 3)+
      annotate("text", x = 1:P, y = rep(max(new_df$perc)+ max(new_df$perc)*0.03, P), label =tmp_reg$star_ann, size = 4)+
      theme_bw()+theme(legend.position = 'bottom', plot.title=element_text(hjust = 0.5, size = 10))
    
  }else{
    if(is.integer(new[, pheno_plot[i]]) & length(unique(new[, pheno_plot[i]])) <= 10){
    
    val <- sort(unique(new[, pheno_plot[i]]))
    new_df <- data.frame(perc = as.vector(sapply(gr_id,  function(y) sapply(val, function(x) sum(new[, pheno_plot[i]] == x & new$gr == paste0('gr',y))/sum(new$gr == paste0('gr',y))))))
    new_df$gr <- unlist(lapply(paste0('gr',gr_id), function(x) rep(x, length(val))))    
    new_df$value <- rep(val, P) 
    new_df$gr <- factor(new_df$gr, levels = paste0('gr', gr_id))
    new_df$value <- factor(new_df$value, levels = val)
    val_max <- max(sapply(gr_id, function(x) sum(new_df$perc[new_df$gr == paste0('gr', x)])))
    
    
    list_pl[[i]] <- ggplot(new_df, aes(x = gr, y = perc, fill = value, color = gr))+
      geom_bar(size = 1, alpha = 0.8, stat = 'identity', width = 0.7)+
      xlab('')+ ylab('Percentage')+
      scale_color_manual(values = gr_color)+
      scale_fill_grey(start=0.9, end=0.1)+
      guides(color = FALSE)+
      labs(fill = pheno_plot[i])+
      # annotate("text", x = 1:P, y = rep(1.05, P), label =formatC(tmp_reg$pvalue,format="e", digits = 2), size = 3)+
      annotate("text", x = 1:P, y = rep(val_max+ val_max*0.03, P), label =tmp_reg$star_ann, size = 4)+
      theme_bw()+theme(legend.position = 'bottom')
      
    }else{
      if(is.numeric(phenoDat[, pheno_plot[i]]) | (is.integer(phenoDat[, pheno_plot[i]]) & length(unique(phenoDat[, pheno_plot[i]])) >= 10)){
        
        new <- phenoDat[!is.na(phenoDat[, pheno_plot[i]]), ]
        new_df <- data.frame(value = new[, pheno_plot[i]], gr = new$gr)
        new_df$gr <- factor(new_df$gr, levels = paste0('gr', gr_id))
      
        list_pl[[i]] <- ggplot(new_df, aes(x = gr, y = value, fill = gr))+
          geom_violin(alpha = 0.8)+
          geom_boxplot(width=0.2, fill="white")+
          xlab('')+ ylab('')+
          scale_fill_manual(values = gr_color, drop = F)+
          ggtitle(pheno_plot[i])+
          # annotate("text", x = 1:P, y = rep(max(new_df$value) + 3, P), label =formatC(tmp_reg$pvalue,format="e", digits = 2), size = 3)+
          annotate("text", x = 1:P, y = rep(max(new_df$value) + max(new_df$value)*0.03, P), label = tmp_reg$star_ann, size = 4)+
          theme_bw()+theme(legend.position = 'none', plot.title=element_text(hjust = 0.5, size = 10))
      }
    }
  }
  
  ggsave(filename = sprintf('%s%s_%s_cluster%s_PGmethod_HKmetric_phenoSelected%s.png', outFold, type_data, type_input, type_cluster, pheno_plot[i]), width = 3.5, height = 3.7, plot = list_pl[[i]], device = 'png')
  ggsave(filename = sprintf('%s%s_%s_cluster%s_PGmethod_HKmetric_phenoSelected%s.pdf', outFold, type_data, type_input, type_cluster,  pheno_plot[i]), width = 3.5, height = 3.7, plot = list_pl[[i]], device = 'pdf')
  
}

# # tot_pl <- ggarrange(plotlist = list_pl, ncol = 1, align='v')
# ggsave(filename = sprintf('%s%s_%s_cluster%s_PGmethod_HKmetric_phenoSelected.png', outFold, type_data, type_input, type_cluster, type_sim), width = 3.5, height = 3.7, plot = tot_pl, device = 'png')
# ggsave(filename = sprintf('%s%s_%s_cluster%s_PGmethod_HKmetric_phenoSelected.pdf', outFold, type_data, type_input, type_cluster, type_sim), width = 3.5, height = 3.7, plot = tot_pl, device = 'pdf')
# 




