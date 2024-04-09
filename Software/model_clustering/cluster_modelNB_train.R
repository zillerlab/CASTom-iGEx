#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(argparse)
  library(e1071)
  library(Matrix)
  library(pROC)
})


parser <- argparse::ArgumentParser()

parser$add_argument(
  "--clustFile",
  type = "character",
  help = "Full path to the clustering object"
)

parser$add_argument(
  "--nFolds",
  type = "numeric",
  help = "Number of folds for cross-validation",
  default = 5
)

parser$add_argument(
  "--outFold",
  type = "character",
  help = "Folder to save the results",
  default = ""
)


args <- parser$parse_args()


base <- sub(".RData", "", basename(args$clustFile))


clust <- get(load(args$clustFile))

clust_input <- as.data.frame(clust$input_data)
clust_gr <- clust$cl_best$gr


fold_size <- round(nrow(clust_input) / args$nFolds)

labels <- seq_len(nrow(clust_input))
labels_fold <- list()

for (fold in 1:(args$nFolds - 1)) {
  lf <- sample(labels, size = fold_size)

  labels_fold[[fold]] <- lf
  labels <- labels[!(labels %in% lf)]
}

labels_fold[[args$nFolds]] <- labels


model <- list()
res <- list()

for (fold in 1:args$nFolds) {
  fold_test_i <- labels_fold[[fold]]

  input_train <- clust_input[-fold_test_i, ]
  input_train["gr"] <- clust_gr[-fold_test_i]

  input_test <- clust_input[fold_test_i, ]
  gr_test <- clust_gr[fold_test_i]

  model[[fold]] <- e1071::naiveBayes(gr ~ ., input_train)

  test_res <- predict(model[[fold]], input_test, type = "raw")

  res[[fold]] <- data.frame(
    fold = fold,
    pred = max.col(test_res),
    prob = apply(test_res, 1, max),
    test = gr_test,
    correct = max.col(test_res) == gr_test
  )
}


auc <- sapply(
  res,
  function(x) pROC::auc(pROC::roc(x, response = correct, predictor = prob, levels = c(FALSE, TRUE), direction = "<"))
)


model <- model[[which.max(auc)]]

thresh <- pROC::roc(
  res[[which.max(auc)]],
  response = correct,
  predictor = prob,
  levels = c(FALSE, TRUE),
  direction = "<"
)

thresh <- pROC::coords(thresh, ret = c("fdr", "threshold"))
thresh <- subset(thresh, fdr <= 0.05)[1, 2]


model_best <- list(
  model = model,
  res_pval = clust$res_pval[, c(2, 7)],
  thresh = thresh
)


save(model_best, file = paste0(args$outFold, "/", base, "_modelNB.Rdata"))


message(paste("NB model training finished. Best AUC =", round(max(auc), 2)))
