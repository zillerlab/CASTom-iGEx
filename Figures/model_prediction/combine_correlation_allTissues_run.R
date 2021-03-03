# combine results overall tissue
# correlation, MR IVW and MREgg

options(stringsAsFactors=F)
options(max.print=1000)
suppressPackageStartupMessages(library(argparse))
options(bitmapType = 'cairo', device = 'png')

parser <- ArgumentParser(description="combine results MR abd correlations across all tissues")
parser$add_argument("--corrRes_tissue_file", type = "character", nargs = '*', help = "results for correlation analysis tissue spec")
parser$add_argument("--thr_pval", type = "double", default = 0.05, help = "")
parser$add_argument("--tissue_name", type = "character",nargs = '*', help = "tissues")
parser$add_argument("--type_data", type = "character", help = "")
parser$add_argument("--outFold", type="character", help = "Output file [basename only]")

args <- parser$parse_args()
corrRes_tissue_file <- args$corrRes_tissue_file
tissue_name <- args$tissue_name
type_data <- args$type_data
thr_pval <- args$thr_pval
outFold <- args$outFold


########################################################################################################################
# tissue_name <- c('Adipose_Subcutaneous','Adipose_Visceral_Omentum', 'Adrenal_Gland','Artery_Coronary',
#                  'Artery_Aorta', 'Colon_Sigmoid', 'Colon_Transverse', 'Heart_Atrial_Appendage',
#                  'Heart_Left_Ventricle', 'Liver', 'Whole_Blood')
# corrRes_tissue_file <- paste0('/psycl/g/mpsziller/lucia/CAD_UKBB/eQTL_PROJECT/OUTPUT_GTEx/predict_CAD/', tissue_name,'/200kb/CAD_GWAS_bin5e-2/UKBB/devgeno0.01_testdevgeno0/enrichment_CADHARD_res/CADCardioG_perc0.3_dist200000b_correlation_enrich_CADCardioG_relatedPheno.RData')
# mrRes_tissue_file <- paste0('/psycl/g/mpsziller/lucia/CAD_UKBB/eQTL_PROJECT/OUTPUT_GTEx/predict_CAD/', tissue_name,'/200kb/CAD_GWAS_bin5e-2/UKBB/devgeno0.01_testdevgeno0/enrichment_CADHARD_res/CADCardioG_dist200000b_Mendelian_randomization_Egger_tscore_pvalFDRrel0.05.txt')
# mr_type <- 'Egger'
# type_data <- 'tscore'
######################################################################################################################


feat_corr <- feat_corr_pval <- list()
for(i in 1:length(tissue_name)){
  
  t <- tissue_name[i]
  tmp <- get(load(corrRes_tissue_file[i]))
  if(i==1){
    pheno_info <- tmp$pheno
  }
  
  if(type_data == 'tscore'){
    feat_corr[[i]] <- tmp$tscore$cor_spearman
    feat_corr_pval[[i]] <- tmp$tscore$cor_pval
  }
  
  if(type_data == 'tot_path'){
    feat_corr[[i]] <- tmp$pathScore$cor_spearman
    feat_corr_pval[[i]] <- tmp$pathScore$cor_pval
  }
  
}

feat_corr <- do.call(cbind, feat_corr)
feat_corr_pval <- do.call(cbind, feat_corr_pval)
colnames(feat_corr) <- colnames(feat_corr_pval) <- tissue_name
# matrix with only significant corr
feat_corr_sign <- feat_corr
for(i in 1:length(tissue_name)){
  feat_corr_sign[feat_corr_pval[,i]> 0.05, i] <- NA
}

feat_corr <- cbind(feat_corr, pheno_info)
feat_corr_pval <- cbind(feat_corr_pval, pheno_info)
feat_corr_sign <- cbind(feat_corr_sign, pheno_info)
feat_corr_sign <- feat_corr_sign[rowSums(!is.na(feat_corr_sign[, tissue_name])) > 0,]

# save
write.table(x = feat_corr, file = sprintf('%s%s_correlation_allTissues.txt', outFold, type_data), 
            quote = F, sep = '\t', col.names = T, row.names = F)

write.table(x = feat_corr_pval, file = sprintf('%s%s_correlation_pvalue_allTissues.txt', outFold, type_data), 
            quote = F, sep = '\t', col.names = T, row.names = F)

write.table(x = feat_corr_sign, file = sprintf('%s%s_correlation_sign_pval%s_allTissues.txt', outFold, type_data, as.character(thr_pval)), 
            quote = F, sep = '\t', col.names = T, row.names = F)

