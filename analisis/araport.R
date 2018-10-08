araport <- function(genes){
  genes_araport <- data.frame(matrix(ncol=5, nrow=0))
  colnames(genes_araport) <- c("briefDescription", "curatorSummary", "symbol", "annotation", "synonym")
  for(gen in genes){
    print(gen)
    ara <- thalemine(gen)
    ara$annotation[1] <- paste(unique(ara$annotation[grep("GO:", ara$identifier)]), collapse=", ")
    ara$synonym[1] <- paste(unique(ara$synonym), collapse=", ")
    ara <- ara[1, c("briefDescription", "curatorSummary", "symbol", "annotation", "synonym")]
    genes_araport <- rbind(genes_araport, ara)
  }
  return(genes_araport)
}

library(reticulate)
library(igraph)
setwd("~/doctorado/programacion/redes_mixtas/")

(load("pipeline_archivos/refT27/3_transcriptoma.Rdata"))
g_genes_t27 <- g_genes
(load("pipeline_archivos/3_transcriptoma.Rdata"))
g_genes_t22 <- g_genes

use_python("~/miniconda3/bin/python", required = T)
source_python("~/doctorado/programacion/api_araport/query_araport.py")
genes_diferentes <- names(V(g_genes_t22))[!(names(V(g_genes_t22)) %in% names(V(g_genes_t27)))]
a<-araport(genes_diferentes)
genes_t22<-araport(names(V(g_genes_t22)))
genes_t27<-araport(names(V(g_genes_t27)))
genes_ego<-araport(strsplit2(names(V(g_ego)),":")[, 1])
