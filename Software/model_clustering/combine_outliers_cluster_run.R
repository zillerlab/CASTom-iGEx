# cluster (multiple cohort combined)

options(stringsAsFactors=F)
options(max.print=1000)
suppressPackageStartupMessages(library(argparse))
options(bitmapType = 'cairo', device = 'png')


parser <- ArgumentParser(description="create a unique file with samples to remove (across tissues and analysis total/noMHC)")
parser$add_argument("--sampleFiles", type = "character", default = 'NA', nargs = '*', help = "")
parser$add_argument("--type_cluster", type = "character", default = 'All', help = "All, Cases, Controls")
parser$add_argument("--type_data", type = "character", help = "tscore, path_Reactome or path_GO")
parser$add_argument("--type_input", type = "character", default = 'original', help = "original or zscaled")
parser$add_argument("--outFold", type="character", help = "Output file [basename only]")

args <- parser$parse_args()
sampleFiles <- args$sampleFiles
type_cluster <- args$type_cluster
type_data <- args$type_data
type_input <- args$type_input
outFold <- args$outFold

####################################################################################################################
# tissues <- read.table('/home/luciat/eQTL_PROJECT/OUTPUT_GTEx/Tissue_PGCgwas', h=F, stringsAsFactors = F)$V1
# folds = paste0('/home/luciat/eQTL_PROJECT/OUTPUT_GTEx/predict_PGC/', tissues, '/200kb/PGC_GWAS_bin1e-2/Meta_Analysis_SCZ/devgeno0.01_testdevgeno0/')
# sampleFiles <- c(paste0(folds, 'tscore_zscaled_clusterCases_PGmethod_umap_oultiers.txt'), paste0(folds, 'excludeMHC_tscore_zscaled_clusterCases_PGmethod_umap_oultiers.txt'))
# sampleFiles <- c('/home/luciat/eQTL_PROJECT/OUTPUT_CMC/predict_PGC/200kb/Meta_Analysis_SCZ/devgeno0.01_testdevgeno0/tscore_zscaled_clusterCases_PGmethod_umap_oultiers.txt',
#                  '/home/luciat/eQTL_PROJECT/OUTPUT_CMC/predict_PGC/200kb/Meta_Analysis_SCZ/devgeno0.01_testdevgeno0/excludeMHC_tscore_zscaled_clusterCases_PGmethod_umap_oultiers.txt', 
#                  sampleFiles)
#####################################################################################################################

sampleAnn <- vector(mode = 'list', length = length(sampleFiles))
for(i in 1:length(sampleFiles)){
  if(file.exists(sampleFiles[i])){
    print(i)
    tmp <- read.table(sampleFiles[i], h=T, stringsAsFactors = F, sep= '\t')
    sampleAnn[[i]] <- tmp
  }
}
sampleAnn <- do.call(rbind, sampleAnn)
sampleAnn <- sampleAnn[!duplicated(sampleAnn$Temp_ID), ]

write.table(file = sprintf('%ssamples_to_remove_outliersUMAP_%s_%s_cluster%s.txt', outFold, type_data, type_input, type_cluster), x = sampleAnn, 
            col.names = T, row.names = F, sep = '\t', quote = F)










