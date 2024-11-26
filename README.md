# CASTom-iGEx pipeline: Cases stratification from imputed gene expression
CASTom-iGEx is a multiple step pipeline aiming at converting individual level genotype data into meaningful biological entities such as genes and molecular pathways in order to aggregate single variant effects and detect distinct pathomechanisms in specific patient subgroups. 

## About the project
The framework is divided in 3 separate modules:
1. **Gene-expression model estimation**: based on reference panels with matching RNAseq and genotype data, gene expression activity is estimated based on cis-genetic effects integrating additional biological knowledge on variants. 
2. **Imputation of genes and pathways individual levels on genotype-only cohorts leading to TWAS (transcriptome-wide association) and PALAS (pathway activity level association) studies**: based on the previous module, gene expression is imputed on large-scale cohorts and additionally converted into individual pathway activity levels. Genes and pathways are then tested for association with the disease of interest (TWAS and PALAS). 
3. **Patients stratification** based on imputed gene expression, it quantifies differences in genetic liability distribution
across directly interpretable biological process and pathways as well as clinical and endophenotypic parameters.

![](./overview.png)

## Installation
Requires R >= 4.0. Install using Git:
```bash
git clone https://github.com/zillerlab/CASTom-iGEx.git
```

Required R packages for the complete pipeline can be installed with a simple R script from the CASTom-iGEx directory:
```bash
Rscript Software/install_requirements.R
```

<details>
<summary><b>List of required R packages for the complete pipeline</b></summary>

- argparse
- bigmemory
- biomaRt
- circlize
- coin
- cowplot
- data.table
- doParallel
- gep2pep
- ggExtra
- ggpubr
- ggrepel
- ggsci
- ggsignif
- glmnet
- GO.db
- gridExtra
- igraph
- lattice
- limma
- lme4
- lmtest
- MASS
- Matrix
- matrixStats
- nloptr
- nnet
- pROC
- pheatmap
- pryr
- qvalue
- RColorBrewer
- rlist
- RNOmni
- rstatix
- SparseM
- sva
- tidyverse
- umap

</details>

## Usage
Guide on how to preprocess custom genetic data to work with the pipeline can be found here:
https://github.com/zillerlab/CASTom-iGEx/wiki/Processing-genetic-data-to-work-with-CASTom%E2%80%90iGEx.

For details of each module of the pipeline please refer to 
* [Module 1](https://github.com/zillerlab/CASTom-iGEx/tree/master/Software/model_training)
* [Module 2](https://github.com/zillerlab/CASTom-iGEx/tree/master/Software/model_prediction)
* [Module 3](https://github.com/zillerlab/CASTom-iGEx/tree/master/Software/model_clustering)

## Reference PriLer models and example workflow
Pretrained PriLer models are available for a large number of tissues and can be accessed here: 
https://figshare.com/account/projects/163249/articles/22347574

Example workflow based on simulated genetic data can be found here: 
https://figshare.com/account/projects/163249/articles/22347574

## References
The pipeline and its application is described in details in:
[Trastulla, L., Dolgalev, G., Moser, S. et al. Distinct genetic liability profiles define clinically relevant patient strata across common diseases. Nat Commun 15, 5534 (2024)](https://doi.org/10.1038/s41467-024-49338-2)