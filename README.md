# CASTom-iGEx pipeline: Cases stratification from imputed gene expression
CASTom-iGEx is a multiple step pipeline aiming at convert individualing level genotype data into meaningful biological entities such as genes and molecular pathways in order to 1) aggregate single variant effects 2) investigate intermediate mechanisms leading to the disease etiology through shared molecular pathways and 3) detect distinct pathomechanisms in specific patient subgroups. 

## About the project
The framework is divided in 3 separate modules:
1. **Gene-expression model estimation**: based on reference panels with matching RNAseq and genotype data, gene expression activity is estimated based on cis-genetic effects integrating additional biological knowledge on variants (prior). 
2. **Imputation of genes and pathways individual levels on genotype-only cohorts leading TWA (transcriptome-wide association) and PALA (pathway activity level association) studies**: based on previous Module, gene expression is imputed on large-scale cohorts and additionally converted into individual pathway activity levels. Genes and pathways are then tested for association with the disease of interest (TWAS and PALAS). If multiple phenotypes are available, this module also estimates the causality of endophenotypes for a certain disease etiology based on genes and pathway association via Mendelian Randomization.
3. **Patients stratification** based on imputed gene expression quantifying differences in genetic liability distribution
across directly interpretable biological process and pathways as well as clinical and endophenotypic parameters.

![](./overview.png)


