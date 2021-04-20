suppressPackageStartupMessages(library(argparse))
suppressPackageStartupMessages(library(biomaRt))
suppressPackageStartupMessages(library(GO.db))

parser <- ArgumentParser(description="Annotate gene ontology")
parser$add_argument("--folder", type = "character", help = "folder to save results")

folder <- args$folder

####################################################################
# create annotation for go terms
ensembl <- useEnsembl(biomart="ensembl", dataset="hsapiens_gene_ensembl",host="ensembl.org", GRCh=37)
# ens2Name <- getBM(attributes=c('ensembl_gene_id',"hgnc_symbol"),mart=ensembl)
# ens2Name <- ens2Name[!duplicated(ens2Name[,1]),]

GO_list <- as.list(GOTERM)

GO_list_new <- vector(mode = 'list', length = length(GO_list))
id <- 1

# for each GO term, add genes (both annotation)
for(i in id:length(GO_list)){

  print(i)
  GO_list_new[[i]]$GOID <- GO_list[[i]]@GOID
  GO_list_new[[i]]$Term <- GO_list[[i]]@Term
  GO_list_new[[i]]$Ontology <-GO_list[[i]]@Ontology
  GO_list_new[[i]]$Synonym <-GO_list[[i]]@Synonym
  GO_list_new[[i]]$Secondary <- GO_list[[i]]@Secondary

  # note: in this way for each term the genes from child terms are not included
  tmp <- getBM(attributes = c("hgnc_symbol", "ensembl_gene_id"), filters = 'go', values = GO_list[[i]]@GOID, mart = ensembl)
  GO_list_new[[i]]$geneIds <- tmp$hgnc_symbol
  GO_list_new[[i]]$geneIds2 <- tmp$ensembl_gene_id

}

# save result (RData format)
# remove GO terms with no genes
id_rem <- which(sapply(GO_list_new, function(x) length(x$geneIds2)==0))
GO_list_new <- GO_list_new[-id_rem]

GO_term_annotated <- GO_list_new
save(GO_term_annotated,file = sprintf('%sGOterm_geneAnnotation_allOntologies.RData', folder))

#################
# create df
df <- data.frame(go_id = sapply(GO_term_annotated, function(x) x$GOID), 
                 go_name = sapply(GO_term_annotated, function(x) x$Term), go_type = sapply(GO_term_annotated, function(x) x$Ontology), stringsAsFactors = F)
write.table(x = df, file = sprintf('%sGOterm_names.txt', folder), sep = '\t', col.names = T, row.names = F, quote = F)





