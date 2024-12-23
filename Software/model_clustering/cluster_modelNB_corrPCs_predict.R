#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(argparse)
  library(dplyr)
  library(e1071)
  library(Matrix)
  library(pROC)
})


parser <- argparse::ArgumentParser()

parser$add_argument(
  "--inputFile",
  type = "character",
  help = "Full path to the predicted T-scores for the new cohort"
)

parser$add_argument(
  "--modelNB",
  type = "character",
  help = "Full path to the Naive Bayes model from a reference dataset"
)

parser$add_argument(
  "--sampleAnnNew_file",
  type = "character",
  help = "Path to the sample annotation file for the new cohort"
)

parser$add_argument(
  "--functR",
  type = "character",
  help = "Path to additional functions"
)

parser$add_argument(
  "--type_cluster",
  type = "character",
  help = "Type of clustering - Cases, Control or All (default)",
  default = "All"
)

parser$add_argument(
  "--outFold",
  type = "character",
  help = "Folder to save the results",
  default = NULL
)


args <- parser$parse_args()

source(args$functR)

model <- get(load(args$modelNB))


sampleAnn_new <- read.table(args$sampleAnnNew_file, h=T, stringsAsFactors = F)
name_cov <- setdiff(colnames(sampleAnn_new),c('Individual_ID', 'genoSample_ID', 'RNASample_ID', 'Dx', 'Sex', 'Age', 'Array'))


if(args$type_cluster == 'Cases'){
  if('Dx' %in% colnames(sampleAnn_new)){
    sampleAnn_new <- sampleAnn_new[sampleAnn_new$Dx == 1,]
  }
} else {
  if(args$type_cluster == 'Controls'){
    if('Dx' %in% colnames(sampleAnn_new)){
      sampleAnn_new <- sampleAnn_new[sampleAnn_new$Dx == 0,]
    }
    
  }else{
    if(args$type_cluster != 'All')
      stop('type_cluster must be either Cases or Controls or All')
  }
}


type_data <- 'tscore'

load_output <- load_input_matrix(inputFile = args$inputFile, 
                                 sampleAnn = sampleAnn_new, 
                                 res_pval = model$res_pval, 
                                 split_tot = 0, 
                                 id_info = 1)

scoreMat <- load_output$scoreMat
sampleAnn_new <- load_output$sampleAnn

input_data_notcorr <- scale(scoreMat)
attr(input_data_notcorr, "scaled:scale") <- NULL
attr(input_data_notcorr, "scaled:center") <- NULL

input_data <- matrix(ncol = ncol(input_data_notcorr), nrow = nrow(input_data_notcorr))
rownames(input_data) <- rownames(input_data_notcorr)
colnames(input_data) <- colnames(input_data_notcorr)

fmla <- as.formula(paste('g ~', paste0(name_cov, collapse = '+')))

for(i in seq_len(ncol(input_data_notcorr))){
  # print(i)
  tmp <- data.frame(g = input_data_notcorr[,i], sampleAnn_new[, name_cov])
  if(any(is.na(tmp$g))){
    # only needed when imputed gene expression is computed in non-harmonized data
    # some genes might have variance zero
    input_data[,i] <- 0
  }else{
    reg <- lm(fmla, data = tmp)
    input_data[,i] <- reg$residuals
  }
}

input_data <- sapply(seq_len(ncol(input_data)), 
                     function(x) input_data[, x]*model$res_pval[model$res_pval[, 1] == colnames(input_data)[x], 2])
colnames(input_data) <- colnames(scoreMat)


pred <- predict(model$model, input_data, type = "raw")

colnames(pred) <- paste0("gr_", colnames(pred))

res <- data.frame(pred = max.col(pred), prob = apply(pred, 1, max))
res$confident <- with(res, prob >= model$thresh)

n_confident <- nrow(res[res$confident, ])


clust_new <- list(
  probability = cbind(sampleAnn_new[, c("Individual_ID", "Dx")], pred),
  sampleAnn = sampleAnn_new,
  data_new = input_data,
  cl_new = data.frame(id = sampleAnn_new[, "Individual_ID"], gr = res$pred, confident = res$confident)
)


save(
  clust_new,
  file = paste0(
    args$outFold, ifelse(is.null(args$outFold), "", "/"),
    sprintf("tscore_corrPCs_zscaled_predictCluster%s_PGmethod_HKmetric.RData", args$type_cluster)
  )
)

message(paste0(
  "Prediction for the new cohort finished. ",
  round(n_confident / nrow(res) * 100, 2), "% samples are confidently assigned"
))
