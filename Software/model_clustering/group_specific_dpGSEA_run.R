# perform drug enrichment analysis

options(stringsAsFactors=F)
options(max.print=1000)
suppressPackageStartupMessages(library(argparse))
suppressPackageStartupMessages(library(data.table))
suppressPackageStartupMessages(library(doParallel))
suppressPackageStartupMessages(library(tidyverse))
options(bitmapType = 'cairo', device = 'png')

parser <- ArgumentParser(description="adjust features table")
parser$add_argument("--featRelFile", type = "character", default = NULL, help = "file association features wilcoxon test")
parser$add_argument("--tissue", type="character", help = "tissue to focus on for features")
parser$add_argument("--type_input", type="character", help = "")
parser$add_argument("--python_fold", type="character", help = "")
parser$add_argument("--env_name", type="character", help = "")
parser$add_argument("--miniconda_call", type="character", help = "for clustering")
parser$add_argument("--outFold", type="character", help = "Output file [basename only]")

args <- parser$parse_args()
featRelFile <- args$featRelFile
outFold <- args$outFold
type_input <- args$type_input
tissue <- args$tissue
python_fold <- args$python_fold
miniconda_call <- args$miniconda_call
env_name <- args$env_name

##################################################################################
# featRelFile <- 'OUTPUT_GTEx/predict_CAD/Liver/200kb/CAD_GWAS_bin5e-2/UKBB/devgeno0.01_testdevgeno0/CAD_HARD_clustering/update_corrPCs/tscoreOriginal_corrPCs_tscoreClusterCases_featAssociation.RData'
# outFold <- 'OUTPUT_GTEx/predict_CAD/Liver/200kb/CAD_GWAS_bin5e-2/UKBB/devgeno0.01_testdevgeno0/CAD_HARD_clustering/update_corrPCs/'
# tissue <- 'Liver'
# type_input <- 'tscore'
# python_fold <- '/psycl/g/mpsziller/lucia/drug_targeting/'
# miniconda_call <- '/u/luciat/miniconda3/etc/profile.d/conda.sh'
# env_name <- 'dpGSEA_py'
##################################################################################

featRel <- get(load(featRelFile))
id_t <- which(featRel$tissues == tissue)
featRel_t <- featRel$test_feat[[id_t]]

featRel_t <- featRel_t %>% 
  filter(pval_corr <= 0.05) %>%
  mutate(t = estimates/(CI_up - CI_low)*3.92) %>%
  select(comp, feat, pval, estimates, pval_corr, t) %>% 
  rename(SYMBOLS = feat, P.value = pval, logFC = estimates, adj.P.Val = pval_corr)

system(paste('source', miniconda_call))
comp_name <- unique(featRel_t$comp) 

drug_enrich <- list()
for(i in 1:length(comp_name)){
  
  tmp <- featRel_t %>% filter(comp == comp_name[i]) %>% arrange(P.value) %>%
     select(-comp)
    
  file_gr <- sprintf('%s%s_tscoreOriginal_corrPCs_%sClusterCases_featAssociation_%s.csv', 
                     outFold, tissue, type_input, comp[i])
  
  write.table(as.data.frame(tmp), file = file_gr, 
              col.names = T, row.names = F, sep = ',', quote = F)
  
  system(sprintf('conda run -n %s python %sdpGSEA.py -tt %s -dr %spms/L1K_P50.csv -i 1000 -sd 10 -o %sdpGSEA_%s.tsv',
                env_name, python_fold, file_gr, python_fold, outFold, comp[i]))

  drug_enrich[[i]] <- read.delim(sprintf('%sdpGSEA_%s.tsv', outFold, comp[i]), h=T, stringsAsFactors = F, sep = '\t')
  drug_enrich[[i]] <- drug_enrich[[i]] %>% arrange(ES_p)
  
}

WHO_ATC <- read.csv(sprintf('%sWHO ATC-DDD 2021-12-03.csv', python_fold), h=T, stringsAsFactors = F)
atc_class <- WHO_ATC[nchar(WHO_ATC$atc_code) <= 5,]

drug_atc <- list()  
for(i in 1:length(comp_name)){
  
  drug_name <- drug_enrich[[i]] %>% filter(NES_95 == 1) %>% select(drug)
  drug_name <- unique(sapply(str_split(drug_name$drug, "_"), function(x) x[1]))
  
  atc_id <- WHO_ATC[WHO_ATC$atc_name %in% drug_name,]
  atc_id$class1 <- sapply(atc_id$atc_code, function(x) atc_class$atc_name[atc_class$atc_code == substr(x, 1, 1)]) 
  atc_id$class3 <- sapply(atc_id$atc_code, function(x) atc_class$atc_name[atc_class$atc_code == substr(x, 1, 3)]) 
  atc_id$class4 <- sapply(atc_id$atc_code, function(x) atc_class$atc_name[atc_class$atc_code == substr(x, 1, 4)]) 
  atc_id$class5 <- sapply(atc_id$atc_code, function(x) atc_class$atc_name[atc_class$atc_code == substr(x, 1, 5)]) 
  
  drug_atc[[i]] <- atc_id  %>% distinct(atc_code, .keep_all = T) 
  drug_atc[[i]]$comp <- comp_name[i]
  #%>% 
   # filter(class1 %in% c('CARDIOVASCULAR SYSTEM', 'BLOOD AND BLOOD FORMING ORGANS', 'ALIMENTARY TRACT AND METABOLISM'))
}

drug_atc  <- do.call(rbind, drug_atc)
file_save <- sprintf('%sdpGSEA_ATCcodes_allgroups.tsv', outFold)
write.table(as.data.frame(drug_atc), file = file_save, 
            col.names = T, row.names = F, sep = '\t', quote = F)




