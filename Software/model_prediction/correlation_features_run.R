# compute correlation matrix for Tscores/pathways for a specific tissue
# use a random samples setting 

options(stringsAsFactors=F)
options(max.print=1000)
suppressPackageStartupMessages(library(argparse))
suppressPackageStartupMessages(library(Matrix))
options(bitmapType = 'cairo', device = 'png')

parser <- ArgumentParser(description="correlation among features, estimates")
parser$add_argument("--inputFile", type = "character", nargs = '*', help = "file to be loaded (predicted tscore or pathScore)")
parser$add_argument("--sampleAnnFile", type = "character", help = "file with samples to be used")
parser$add_argument("--tissues_name", type = "character", help = "name tissue")
parser$add_argument("--split_tot", type = "integer", default = 0, help = "if 0 then inpuntFile load alone, otherwise splitted version")
parser$add_argument("--type_data", type = "character", nargs = '*', help = "tscore, path_Reactome, path_GO or tot_path")
parser$add_argument("--outFold", type="character", help = "Output file [basename only]")

args <- parser$parse_args()
inputFile <- args$inputFile
sampleAnnFile <- args$sampleAnnFile
tissue_name <- args$tissue_name
split_tot <- args$split_tot
type_data <- args$type_data
outFold <- args$outFold

#########################################################################################################################
# inputFile <- c('/psycl/g/mpsziller/lucia/CAD_UKBB/eQTL_PROJECT/OUTPUT_GTEx/predict_CAD/Liver/200kb/CAD_GWAS_bin5e-2/UKBB/devgeno0.01_testdevgeno0/Pathway_Reactome_scores.RData',
#               '/psycl/g/mpsziller/lucia/CAD_UKBB/eQTL_PROJECT/OUTPUT_GTEx/predict_CAD/Liver/200kb/CAD_GWAS_bin5e-2/UKBB/devgeno0.01_testdevgeno0/Pathway_GO_scores.RData')
# sampleAnnFile <- '/psycl/g/mpsziller/lucia/CAD_UKBB/eQTL_PROJECT/INPUT_DATA_GTEx/CAD/Covariates/UKBB/covariateMatrix_forCorrelation.txt'
# split_tot <- 0
# type_data <- c('Reactome', 'GO')
# outFold <- '/psycl/g/mpsziller/lucia/CAD_UKBB/eQTL_PROJECT/OUTPUT_GTEx/predict_CAD/Liver/200kb/CAD_GWAS_bin5e-2/UKBB/devgeno0.01_testdevgeno0/'
# tissues_name <- 'Liver'
##########################################################################################################################

sampleAnn <- read.table(sampleAnnFile, h=T, stringsAsFactors = F, check.names = F)

# load input matrix 
if(split_tot == 0){
  
  if(length(inputFile) == 1){
    scoreMat <- get(load(inputFile))  
    # filter out based on samples and ids
    id_el <- scoreMat[,1]
    scoreMat <- scoreMat[match(id_el,scoreMat[,1]), ]
    
    common_samples <- intersect(sampleAnn$Individual_ID, colnames(scoreMat))
    sampleAnn <- sampleAnn[match(common_samples, sampleAnn$Individual_ID),]
    scoreMat <- t(scoreMat[,match(common_samples,colnames(scoreMat))])
    
    rownames(scoreMat) <- common_samples
    colnames(scoreMat) <- id_el
  }else{
    
    tmp <- list()
    for(i in 1:length(inputFile)){
      tmp[[i]] <- get(load(inputFile[i]))  
      tmp[[i]][, 1] <- paste0(tmp[[i]][, 1], '_type_', type_data[i])
      id_el <- tmp[[i]][,1]
      tmp[[i]] <- tmp[[i]][match(id_el,tmp[[i]][,1]), ]
      
      common_samples <- intersect(sampleAnn$Individual_ID, colnames(tmp[[i]]))
      sampleAnn <- sampleAnn[match(common_samples, sampleAnn$Individual_ID),]
      tmp[[i]] <- t(tmp[[i]][,match(common_samples,colnames(tmp[[i]]))])
      
      rownames(tmp[[i]]) <- common_samples
      colnames(tmp[[i]]) <- id_el
      
    }
    scoreMat <- do.call(cbind, tmp)
    
  }

}else{
  
  ###### load score Mat #######
  scoreMat_list <- vector(mode = 'list', length = split_tot)
  samplesID <- vector(mode = 'list', length = split_tot)
  elementID <- NULL
  
  for(i in 1:split_tot){
    
    print(i)
    if(file.exists(sprintf('%s%i.RData', inputFile, i))){
      tmp <- get(load(sprintf('%s%i.RData', inputFile, i)))
      elementID <- c(elementID,tmp[,1])
      samplesID[[i]] <- intersect(sampleAnn$Individual_ID, colnames(tmp))
      scoreMat_list[[i]] <- t(tmp[,match(samplesID[[i]],colnames(tmp))])
    }else{
      print(sprintf('split %i does not exist', i))
      split_tot <- split_tot - 1
    }
  }
  
  print(split_tot)
  # check samplesID always the same
  if(!all(table(unlist(samplesID)) == split_tot)){print('ERROR: wrong name annotations')}
  
  scoreMat <- do.call(cbind, scoreMat_list)
  colnames(scoreMat) <- elementID
  rm(scoreMat_list)
  
  # filter out elements that are repeated twice:
  id_dup <- names(which(table(colnames(scoreMat)) > 1)) 
  scoreMat <- scoreMat[, !colnames(scoreMat) %in% id_dup]
  
  rownames(scoreMat) <- samplesID[[1]]
  # remove sample that have NAs
  id_s <- rowSums(is.na(scoreMat)) == 0
  if(!all(id_s)){scoreMat <- scoreMat[id_s, ]}
  
  common_samples <- intersect(sampleAnn$Individual_ID, rownames(scoreMat))
  sampleAnn <- sampleAnn[match(common_samples, sampleAnn$Individual_ID),]
  scoreMat <- scoreMat[match(common_samples,rownames(scoreMat)),]

}

# remove higly correlated features, keep highest association
cor_score <- cor(scoreMat)

# save results in RData object
cor_res <- list(cor = cor_score, sample = sampleAnn)
if(length(type_data) == 2){
  type_data <- 'tot_path'
}
save(cor_res, file = sprintf('%scorrelation_estimate_%s.RData', outFold, type_data))






