requirements_cran <- c(
  "argparse",
  "bigmemory",
  "BiocManager",
  "circlize",
  "clustAnalytics",
  "coin",
  "cowplot",
  "data.table",
  "doParallel",
  "ggExtra",
  "ggpubr",
  "ggrepel",
  "ggsci",
  "ggsignif",
  "glmnet",
  "gridExtra",
  "igraph",
  "lattice",
  "lme4",
  "lmtest",
  "MASS",
  "Matrix",
  "matrixStats",
  "nloptr",
  "nnet",
  "pheatmap",
  "pROC",
  "pryr",
  "RColorBrewer",
  "rlist",
  "RNOmni",
  "rstatix",
  "SparseM",
  "tidyverse"
)

requirements_bioc <- c(
  "biomaRt",
  "GO.db",
  "limma",
  "qvalue",
  "sva",
  "umap"
)


requirements_cran <- requirements_cran[!requirements_cran %in% installed.packages()]

install.packages(requirements_cran, repos = "https://cloud.r-project.org")

BiocManager::install(requirements_bioc)
