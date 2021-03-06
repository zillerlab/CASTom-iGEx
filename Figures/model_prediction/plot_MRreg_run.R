# plot specific dotplot for MR analysis

options(stringsAsFactors=F)
options(max.print=1000)
suppressPackageStartupMessages(library(argparse))
suppressPackageStartupMessages(library(data.table))
suppressPackageStartupMessages(library(doParallel))
suppressPackageStartupMessages(library(ggplot2))
suppressPackageStartupMessages(library(RColorBrewer))
suppressPackageStartupMessages(library(ggpubr))
suppressPackageStartupMessages(library(ggsci))
suppressPackageStartupMessages(library(qvalue))
suppressPackageStartupMessages(library(ggrepel))
options(bitmapType = 'cairo', device = 'png')

parser <- ArgumentParser(description="plot MR single association")
parser$add_argument("--MRRes_Egg_file", type = "character", help = "")
parser$add_argument("--MRRes_IVW_file", type = "character", help = "")
parser$add_argument("--corrRes_file", type = "character", help = "")
parser$add_argument("--tissue", type = "character", help = "")
parser$add_argument("--name_exposure", type = "character", help = "")
parser$add_argument("--name_outcome", type = "character", help = "")
parser$add_argument("--exposure_file", type = "character", help = "")
parser$add_argument("--outcome_file", type = "character", help = "")
parser$add_argument("--type_data", type = "character", help = "")
parser$add_argument("--pval_FDR_rel", type = "double",default = 0.05, help = "")
parser$add_argument("--outFold", type="character", help = "Output file [basename only]")

args <- parser$parse_args()
MRRes_Egg_file <- args$MRRes_Egg_file
MRRes_IVW_file <- args$MRRes_IVW_file
corrRes_file <- args$corrRes_file
tissue <- args$tissue
name_exposure <- args$name_exposure
name_outcome <- args$name_outcome
exposure_file <- args$exposure_file
outcome_file <- args$outcome_file
type_data <- args$type_data
pval_FDR_rel <- args$pval_FDR_rel
outFold <- args$outFold

###################################################################################################################
# tissue <- 'Whole_Blood'
# MRRes_file <- paste0('OUTPUT_GTEx/predict_CAD/',tissue,'/200kb/CAD_GWAS_bin5e-2/UKBB/devgeno0.01_testdevgeno0/enrichment_CADHARD_res/Mendelian_randomization_tscore_pvalFDRrel0.05.txt')
# corrRes_file <- paste0('OUTPUT_GTEx/predict_CAD/',tissue,'/200kb/CAD_GWAS_bin5e-2/UKBB/devgeno0.01_testdevgeno0/enrichment_CADHARD_res/perc0.3_correlation_enrich_CAD_HARD_relatedPheno.RData')
# name_exposure <- '30710'
# name_outcome <- 'CAD_HARD'
# exposure_file <- paste0('OUTPUT_GTEx/predict_CAD/',tissue,'/200kb/CAD_GWAS_bin5e-2/UKBB/devgeno0.01_testdevgeno0/pval_Blood_biochemistry_withMed_pheno_covCorr.RData')
# outcome_file <- paste0('OUTPUT_GTEx/predict_CAD/',tissue,'/200kb/CAD_GWAS_bin5e-2/UKBB/devgeno0.01_testdevgeno0/pval_CAD_pheno_covCorr.RData')
# type_data <- 'tscore'
# pval_FDR_rel <- 0.05
# outFold <- paste0('OUTPUT_GTEx/predict_CAD/',tissue,'/200kb/CAD_GWAS_bin5e-2/UKBB/devgeno0.01_testdevgeno0/enrichment_CADHARD_res/perc0.3_')
################################################################################################################

load_data <- function(data_file, pval_FDR, tissue_name, pval_id){
  
  tscore <- list()
  pathR <- list()
  pathGO <- list()
  
  for(i in 1:length(tissue_name)){
    
    # print(tissue_name[i])
    res_pval <- get(load(data_file[i]))
    
    tscore[[i]] <- list()
    pathR[[i]] <- list()
    pathGO[[i]] <- list()
    
    for(j in 1:length(pval_id)){
      
      tscore[[i]][[j]] <- res_pval$tscore[[pval_id[j]]]
      tscore[[i]][[j]] <- res_pval$tscore[[pval_id[j]]]
      pathR[[i]][[j]] <- res_pval$pathScore_reactome[[pval_id[j]]]
      pathGO[[i]][[j]] <- res_pval$pathScore_GO[[pval_id[j]]]
      tscore[[i]][[j]]$tissue <- pathR[[i]][[j]]$tissue <- pathGO[[i]][[j]]$tissue <- tissue_name[i]
      # filter out pathways with only 1 gene and recompute adjusted pvalue
      pathR[[i]][[j]] <- pathR[[i]][[j]][pathR[[i]][[j]]$ngenes_tscore>1, ]
      pathR[[i]][[j]][, 14] <- qvalue(pathR[[i]][[j]][,13])$qvalue
      pathR[[i]][[j]][, 15] <- p.adjust(p = pathR[[i]][[j]][,13], method = 'BH')
      
      pathGO[[i]][[j]] <- pathGO[[i]][[j]][pathGO[[i]][[j]]$ngenes_tscore>1, ]
      pathGO[[i]][[j]][, 16] <- qvalue(pathGO[[i]][[j]][,15])$qvalue
      pathGO[[i]][[j]][, 17] <- p.adjust(p = pathGO[[i]][[j]][,15], method = 'BH')
      
    }
  }
  
  pheno_tscore <- pheno_pathR <- pheno_pathGO <- list()
  tscore_red <- pathR_red <- pathGO_red <- list()
  for(j in 1:length(pval_id)){
    
    pheno_tscore[[j]] <- do.call(rbind, lapply(tscore, function(x) x[[j]]))
    pheno_pathR[[j]] <- do.call(rbind, lapply(pathR, function(x) x[[j]]))
    pheno_pathGO[[j]] <- do.call(rbind, lapply(pathGO, function(x) x[[j]]))
    pheno_tscore[[j]]$new_id <- paste0(pheno_tscore[[j]][,1],'_tissue_',pheno_tscore[[j]]$tissue)
    pheno_pathR[[j]]$new_id <- paste0(pheno_pathR[[j]][,1],'_tissue_',pheno_pathR[[j]]$tissue)
    pheno_pathGO[[j]]$new_id <- paste0(pheno_pathGO[[j]][,2],'_tissue_',pheno_pathGO[[j]]$tissue)
    
    tscore_red[[j]] <- pheno_tscore[[j]][pheno_tscore[[j]][, 10] <= pval_FDR,]
    pathR_red[[j]] <- pheno_pathR[[j]][pheno_pathR[[j]][, 15] <= pval_FDR,]
    pathGO_red[[j]] <- pheno_pathGO[[j]][pheno_pathGO[[j]][, 17] <= pval_FDR,]
    
  }
  return(list(tscore = pheno_tscore, tscore_red = tscore_red, pathR = pheno_pathR, pathR_red = pathR_red, 
              pathGO = pheno_pathGO, pathGO_red = pathGO_red))
  
}

MR_Egg_res <- read.delim(MRRes_Egg_file, h=T, stringsAsFactors = F, sep = '\t')
MR_IVW_res <- read.delim(MRRes_IVW_file, h=T, stringsAsFactors = F, sep = '\t')
corrRes <- get(load(corrRes_file))

# load association results
exposure_res <- get(load(exposure_file))
outcome_res <- get(load(outcome_file))

id <- match(name_exposure, exposure_res$pheno$pheno_id)

exposure_spec <- load_data(exposure_file, pval_FDR = pval_FDR_rel, tissue_name = tissue, pval_id = id)

if(type_data == 'tot_path'){
  
  pathR_pheno <- exposure_spec$pathR_red
  pathGO_pheno <- exposure_spec$pathGO_red
  tot_path_pheno <- lapply(1:length(id), function(x) rbind(cbind(pathR_pheno[[x]], data.frame(type = rep('Reactome', nrow(pathR_pheno[[x]])))), 
                                                           cbind(pathGO_pheno[[x]][, !colnames(pathGO_pheno[[x]]) %in% c('path_id', 'path_ont')], data.frame(type = rep('GO', nrow(pathGO_pheno[[x]]))))))
  for(i in 1:length(tot_path_pheno)){
    tot_path_pheno[[i]]$new_id <- paste0(tot_path_pheno[[i]]$new_id, '_type_', tot_path_pheno[[i]]$type)
    # match with path_ann
    common_p <- intersect(corrRes$path_ann$new_id, tot_path_pheno[[i]]$new_id)
    tot_path_pheno[[i]] <- tot_path_pheno[[i]][match(common_p,tot_path_pheno[[i]]$new_id),]
  }
  res_pheno <- tot_path_pheno
}else{
  
  tscore_pheno <- exposure_spec$tscore_red
  common_g <- lapply(tscore_pheno, function(x) intersect(corrRes$gene_ann$new_id, x$new_id))
  tscore_pheno <- mapply(function(x, y) x[match(y,x$new_id),], x = tscore_pheno, y = common_g, SIMPLIFY = F)
  res_pheno <- tscore_pheno
}


id <- which(outcome_res$pheno$pheno_id == name_outcome)
outcome_spec <- load_data(outcome_file, pval_FDR = pval_FDR_rel, tissue_name = tissue, pval_id = id)
if(type_data == 'tot_path'){
  
  pathR_out <- outcome_spec$pathR[[1]]
  pathGO_out <- outcome_spec$pathGO[[1]]
  tot_path_out <- rbind(cbind(pathR_out, data.frame(type = rep('Reactome', nrow(pathR_out)))), 
                        cbind(pathGO_out[, !colnames(pathGO_out) %in% c('path_id', 'path_ont')], 
                              data.frame(type = rep('GO', nrow(pathGO_out)))))
  
  tot_path_out$new_id <- paste0(tot_path_out$new_id, '_type_', tot_path_out$type)
  common_p <- lapply(tot_path_pheno, function(x) intersect(x$new_id, tot_path_out$new_id))
  tot_path_out_red <- lapply(common_p, function(x) 
    tot_path_out[match(x,tot_path_out$new_id),])
  res_out_red <- tot_path_out_red
  id_beta <- 10 
  name_id <- 'path'
  
}else{
  
  tscore_out <- outcome_spec$tscore[[1]]
  common_p <- lapply(tscore_pheno, function(x) intersect(x$new_id,tscore_out$new_id))
  tscore_out_red <- lapply(common_p, function(x) tscore_out[match(x,tscore_out$new_id),])
  res_out_red <- tscore_out_red
  id_beta <- 5 
  name_id <- 'external_gene_name'
  
}

MREgg_res_red <- MR_Egg_res[match(name_exposure, MR_Egg_res$pheno),]
MRIVW_res_red <- MR_IVW_res[match(name_exposure, MR_IVW_res$pheno),]

for(i in 1:nrow(MREgg_res_red)){
  
  # plot betas
  new_name <-  MREgg_res_red$names_field[i]
  df <- data.frame(exp_beta = res_pheno[[i]][, id_beta], exp_se = res_pheno[[i]][, id_beta+1], 
                   out_beta = res_out_red[[i]][, id_beta], out_se = res_out_red[[i]][, id_beta+1], 
                   name = res_pheno[[i]][, name_id], stringsAsFactors = F)
  df$exp_min <- df$exp_beta - df$exp_se
  df$exp_max <- df$exp_beta + df$exp_se
  df$out_min <- df$out_beta - df$out_se
  df$out_max <- df$out_beta + df$out_se
  df$sign <- sign(df$exp_beta*df$out_beta) 
  if(sign(MRIVW_res_red$MRIVW_est[i]) ==1){
    df$sign <- factor(df$sign, levels = c(-1, 1))
  }else{
    df$sign <- factor(df$sign, levels = c(1, -1))
  }  
  df$name_plot <- df$name
  # df$name_plot[res_out_red[[i]][, id_beta] > thr_plot+3] <- ''
  df$name_plot[order(abs(df$out_beta))[1:(nrow(df)-25)]] <- ''
  if(sign(MRIVW_res_red$MRIVW_est[i]) ==1){
    # remove the not concordant sign
    df$name_plot[sign(df$exp_beta*df$out_beta) == -1] <- ''
  }else{
    df$name_plot[sign(df$exp_beta*df$out_beta) == 1] <- ''
  }
  
  pl <-  ggplot(df, aes(x = exp_beta, y = out_beta, label = name_plot, color = sign))+
    theme_bw()+ 
    ylab(sprintf('Association with outcome\n%s', name_outcome))+
    xlab(sprintf('Association with exposure\n%s', new_name))+
    geom_hline(yintercept = 0, color = 'blue', size = 0.5, alpha = 0.7)+
    geom_vline(xintercept = 0, color = 'blue', size = 0.5,  alpha = 0.7)+
    geom_abline(slope = MRIVW_res_red$MRIVW_est[i], intercept = 0,
                color = 'red')+
    geom_abline(slope = MRIVW_res_red$MRIVW_est[i] - MRIVW_res_red$MRIVW_est_se[i], intercept = 0,
                color = 'red', linetype = 'dashed')+
    geom_abline(slope = MRIVW_res_red$MRIVW_est[i] + MRIVW_res_red$MRIVW_est_se[i], intercept = 0,
                color = 'red', linetype = 'dashed')+
    geom_abline(slope = MREgg_res_red$MREgg_est[i], intercept = MREgg_res_red$MREgg_int[i],
                color = 'orange')+
    geom_abline(slope = MREgg_res_red$MREgg_est[i] - MREgg_res_red$MREgg_est_se[i], intercept = MREgg_res_red$MREgg_int[i] - MREgg_res_red$MREgg_int_se[i],
                color = 'orange', linetype = 'dashed')+
    geom_abline(slope = MREgg_res_red$MREgg_est[i] + MREgg_res_red$MREgg_est_se[i], intercept = MREgg_res_red$MREgg_int[i] + MREgg_res_red$MREgg_int_se[i],
                color = 'orange', linetype = 'dashed')+
    geom_point(size = 0.5, alpha = 0.7)+
    geom_errorbar(aes(ymin=out_min, ymax=out_max), size = 0.5,  alpha = 0.7)+
    geom_errorbar(aes(xmin=exp_min, xmax=exp_max), size = 0.5,  alpha = 0.7)+
    geom_text_repel(force = 20, box.padding = 0.5, size = 2.5, min.segment.length = 0, 
                    arrow = arrow(length = unit(0.15, 'cm'), type = 'open')) +
    theme(legend.position = 'none', plot.title = element_text(size=9), axis.title.y = element_text(size=10), axis.title.x = element_text(size=10),
          axis.text.x = element_text(size = 9), axis.text.y = element_text(size = 9))+
    scale_color_manual(values = c('grey80', 'grey30'))
  
  ggsave(filename = sprintf('%sMRplot_%s_out%s_exp%s.png', outFold, type_data, name_outcome, name_exposure[i]), width = 5, height = 5, plot = pl, device = 'png')
  ggsave(filename = sprintf('%sMRplot_%s_out%s_exp%s.pdf', outFold, type_data, name_outcome, name_exposure[i]), width = 5, height = 5, plot = pl, device = 'pdf')
  
}

