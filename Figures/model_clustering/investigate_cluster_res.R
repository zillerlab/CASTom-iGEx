# load total result and check most relevant info
gene_feat <- read.delim('OUTPUT_GTEx/predict_CAD/AllTissues/200kb/CAD_GWAS_bin5e-2/UKBB/CAD_HARD_clustering/clLiver_tscoreOriginal_tscoreClusterCases_featAssociation.txt', h=T, stringsAsFactors = F, sep = '\t')
gene_info <- read.delim('OUTPUT_GTEx/predict_CAD/AllTissues/200kb/CAD_GWAS_bin5e-2/UKBB/CAD_HARD_clustering/clLiver_tscoreOriginal_tscoreClusterCases_infoGenes.txt', h=T, stringsAsFactors = F, sep = '\t')
path_feat <- read.delim('OUTPUT_GTEx/predict_CAD/AllTissues/200kb/CAD_GWAS_bin5e-2/UKBB/CAD_HARD_clustering/clLiver_path_ReactomeOriginal_tscoreClusterCases_featAssociation.txt', h=T, stringsAsFactors = F, sep = '\t')
path_info <- read.delim('OUTPUT_GTEx/predict_CAD/AllTissues/200kb/CAD_GWAS_bin5e-2/UKBB/CAD_HARD_clustering/clLiver_path_ReactomeOriginal_tscoreClusterCases_infoGenes.txt', h=T, stringsAsFactors = F, sep = '\t')
path_go_feat <- read.delim('OUTPUT_GTEx/predict_CAD/AllTissues/200kb/CAD_GWAS_bin5e-2/UKBB/CAD_HARD_clustering/clLiver_path_GOOriginal_tscoreClusterCases_featAssociation.txt', h=T, stringsAsFactors = F, sep = '\t')
path_go_info <- read.delim('OUTPUT_GTEx/predict_CAD/AllTissues/200kb/CAD_GWAS_bin5e-2/UKBB/CAD_HARD_clustering/clLiver_path_GOOriginal_tscoreClusterCases_infoGenes.txt', h=T, stringsAsFactors = F, sep = '\t')


P <- length(unique(path_feat$comp))
tmp_impr <- matrix(F, nrow = nrow(path_info), ncol = P)
for(i in 1:nrow(path_info)){
  print(i)
  id <- strsplit(path_info$genes_id[i], split = '[,]')[[1]]
  tmp <- gene_feat[gene_feat$tissue == path_info$tissue[i] & gene_feat$feat %in% id, ]
  gr_path <- path_feat$pval[path_feat$feat == path_info$path[i] & path_feat$tissue == path_info$tissue[i]]
  tmp_impr[i,] <- sapply(1:P, function(p) all(gr_path[p] < tmp$pval[tmp$comp == unique(path_feat$comp)[p]]))
}
colnames(tmp_impr) <- unique(path_feat$comp)
tmp_impr <- as.data.frame(tmp_impr)

tmp_impr_go <- matrix(F, nrow = nrow(path_go_info), ncol = P)
for(i in 1:nrow(path_go_info)){
  print(i)
  id <- strsplit(path_go_info$genes_id[i], split = '[,]')[[1]]
  tmp <- gene_feat[gene_feat$tissue == path_go_info$tissue[i] & gene_feat$feat %in% id, ]
  gr_path <- path_go_feat$pval[path_go_feat$feat == path_go_info$path[i] & path_go_feat$tissue == path_go_info$tissue[i]]
  tmp_impr_go[i,] <- sapply(1:P, function(p) all(gr_path[p] < tmp$pval[tmp$comp == unique(path_go_feat$comp)[p]]))
}
colnames(tmp_impr_go) <- unique(path_go_feat$comp)
tmp_impr_go <- as.data.frame(tmp_impr_go)


### look at the results:
check_res <- function(name_path, comp_name, path_feat, path_info, gene_feat, gene_info, tissue = NA){
  
  if(is.na(tissue)){
    new_path_feat <- path_feat[path_feat$comp %in% comp_name & path_feat$feat %in% name_path,]
    new_path_info <- path_info[path_info$path %in% name_path,]
  }else{
    new_path_feat <- path_feat[path_feat$comp %in% comp_name & path_feat$feat %in% name_path & path_feat$tissue %in% tissue,]
    new_path_info <- path_info[path_info$path %in% name_path & path_info$tissue %in% tissue,]
  }
  
  # check genes
  gene_names <- unique(unlist(lapply(new_path_info$genes_id, function(x) strsplit(x, split = '[,]')[[1]])))
  if(is.na(tissue)){
    new_gene_feat <- gene_feat[gene_feat$feat %in% gene_names & gene_feat$comp %in% comp_name, ]
    new_gene_info <- gene_info[gene_info$external_gene_name %in% gene_names, ]
  }else{
    new_gene_feat <- gene_feat[gene_feat$feat %in% gene_names & gene_feat$comp %in% comp_name & gene_feat$tissue %in% tissue, ]
    new_gene_info <- gene_info[gene_info$external_gene_name %in% gene_names & gene_info$tissue %in% tissue, ]
  }
  return(list(path = new_path_feat, path_info = new_path_info, gene = new_gene_feat[order(new_gene_feat$pval),], gene_info = new_gene_info))
}

increas_path <- function(impr_table, path_info, path_feat, comp_name, tissue = NA, pval_thr=0.001){
  
  if(is.na(tissue)){
    new <- path_feat[path_feat$comp %in% comp_name,]  
  }else{
    new <- path_feat[path_feat$comp %in% comp_name & path_feat$tissue %in% tissue,]  
  }
  
  new <- new[match(paste0(path_info$path, path_info$tissue), paste0(new$feat, new$tissue)),]
  new_impr <- new[which(impr_table[,comp_name] & new$pval<= pval_thr),]
  new_impr <- new_impr[order(new_impr$pval),]
  return(new_impr)
  
}

check_res(name_path = name_path, comp_name = comp_names[i], path_feat = path_feat, path_info = path_info, gene_feat = gene_feat, gene_info =  gene_info, tissue = 'Whole_Blood')
check_res(name_path = name_path, comp_name = comp_names[i], path_feat = path_go_feat, path_info = path_go_info, gene_feat = gene_feat, gene_info =  gene_info, tissue = 'Whole_Blood')

### check results ###
comp_names <-  unique(path_feat$comp)
for(i in 1:P){
  
  # reactome
  increas_path(impr_table = tmp_impr, path_info = path_info, path_feat = path_feat, comp_name = comp_names[i],  tissue = 'Whole_Blood')
  # go
  increas_path(impr_table = tmp_impr_go, path_info = path_go_info, path_feat = path_go_feat, comp_name = comp_names[i], tissue = 'Whole_Blood')
}





