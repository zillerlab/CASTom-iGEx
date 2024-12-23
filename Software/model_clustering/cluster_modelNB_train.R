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
  "--seed",
  type = "numeric",
  help = "Random seed for test set partition",
  default = 42
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


set.seed(args$seed)
test_ids <- sample(1:nrow(clust_input), size = round(nrow(clust_input) * 0.2))

input_train <- clust_input[-test_ids, ]
input_train["gr"] <- clust_gr[-test_ids]

input_test <- clust_input[test_ids, ]
gr_test <- clust_gr[test_ids]

model <- e1071::naiveBayes(gr ~ ., input_train)

test_res <- stats::predict(model, input_test, type = "raw")

res <- data.frame(
  pred = max.col(test_res),
  prob = apply(test_res, 1, max),
  test = gr_test,
  correct = max.col(test_res) == gr_test
)

auc <- pROC::auc(pROC::roc(
  res,
  response = correct,
  predictor = prob,
  levels = c(FALSE, TRUE),
  direction = "<"
))

thresh <- pROC::roc(
  res,
  response = correct,
  predictor = prob,
  levels = c(FALSE, TRUE),
  direction = "<"
)

thresh <- pROC::coords(thresh, ret = c("fdr", "threshold"))
# Can't pick the exact coordinate at which FDR <= 0.05, so pick the closest
thresh <- subset(thresh, fdr <= 0.05)[1, 2]


model <- list(
  model = model,
  res_pval = clust$res_pval[, c(2, 7)],
  thresh = thresh
)

save(model, file = paste0(args$outFold, "/", base, "_modelNB.Rdata"))

message(paste("NB model training finished. Test AUC =", round(max(auc), 2)))
