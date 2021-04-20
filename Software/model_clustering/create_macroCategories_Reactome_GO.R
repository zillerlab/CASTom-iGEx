library(igraph)

setwd('/psycl/g/mpsziller/lucia/refData/')

df_R <- read.table('ReactomePathwaysRelation_Human.txt', h=F, stringsAsFactors = F, sep = '\t')
colnames(df_R) <- c('Parent', 'Child')
tot_names <- unique(df_R$Parent, df_R$Child)
class_names <- tot_names[tot_names %in% df_R$Parent & !tot_names %in% df_R$Child]
list_path <- read.delim('ReactomePathways_Human.txt', h=F, stringsAsFactors = F, sep = '\t')
class_names <- list_path[list_path$V1 %in% class_names, ]

gr_R = graph_from_data_frame(df_R, directed = TRUE, vertices = NULL)
distMatrix <- shortest.paths(gr_R, v=class_names$V1, to=V(gr_R))
group_class <- apply(distMatrix, 1, function(x) names(which(x!=Inf)))
Reactome_macro <- list()
for(i in 1:length(group_class)){
  Reactome_macro[[i]] <- list(macro = class_names[i, 1:2], 
                              subclasses = list_path[list_path$V1 %in% group_class[[i]], 1:2])
  Reactome_macro[[i]]$subclasses$shortestPath <- distMatrix[i,match(Reactome_macro[[i]]$subclasses$V1,colnames(distMatrix))]
}

save(Reactome_macro, file = 'ReactomePathways_macro_2021.RData')


