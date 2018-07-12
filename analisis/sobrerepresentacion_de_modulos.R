##############################################################################################
#ANALISIS DE FUNCIONALIDAD DE MODULOS 
#Autor: Andrés Rabinovich
#Creación: 06/07/2018
#Última modificación: XX/XX/XXXX (por XXX)
##############################################################################################

#Librerías que necesita
library(org.At.tair.db)
library(GO.db)
library(igraph)
library(limma) #Para incluir strsplit2
setwd("/home/arabinov/doctorado/programacion/redes_mixtas/")

#Levantamos las redes
(load("pipeline_archivos/6_red_mixta.Rdata"))

#Palabras clave en análisis de funcionalidad
keywords <- c("cold", "splic", "circadian", "rRNA", " RNA", "negative regulation")

#Encuentra los enlaces inter e intra redes
enlaces           <- ends(g, E(g))
enlaces_interred  <- setdiff(grep(":", enlaces[, 2]), grep(":", enlaces[, 1]))
genes             <- unique(enlaces[enlaces_interred, 1])
enlaces_interred  <- enlaces[enlaces_interred, ]

#Buscamos la función de cada uno de nuestros genes
evidencia   <- c("EXP", "IDA", "IPI", "IMP", "IGI", "IEP")
goterms     <- unlist(Term(GOTERM))
go          <- as.data.frame(org.At.tairGO)
go          <- go[go$gene_id %in% genes, ]
go          <- go[go$Ontology == "BP", ]
go          <- go[go$Evidence %in% evidencia, ]
go          <- unique(go[, c("gene_id", "go_id")])
go$funcion  <- goterms[go$go_id]
go          <- aggregate(funcion ~ gene_id, go, paste, collapse=", ")
funciones   <- sapply(keywords, function(x){return(grep(x, go[, 2]))})
funciones

go_completo     <-as.data.frame(org.At.tairGO2ALLTAIRS)
genes_de_bines  <- strsplit2(enlaces_interred[, 2], ":")[, 1]
todos_los_genes <- unique(c(genes, genes_de_bines))
#todos_los_genes <- unique(c(names(V(g))[-grep(":", names(V(g)))], genes_de_bines))
go_completo     <- go_completo[go_completo$gene_id %in% todos_los_genes, ]
go_completo     <- go_completo[go_completo$Ontology == "BP", ]
go_completo     <- go_completo[go_completo$Evidence %in% c("RCA", "ISS", evidencia), ]
go_completo     <- unique(go_completo[, c("gene_id", "go_id")])


n_universo          <- length(todos_los_genes)
sobrerepresentacion <- setNames(lapply(genes, function(gen){
  blancas                   <- c(gen, genes_de_bines[enlaces_interred[, 1] == gen])
  n_blancas                 <- length(blancas)
  n_negras                  <- n_universo - n_blancas
  terminos_a_analizar       <- go_completo$go_id[go_completo$gene_id %in% blancas]
  n_blancas_en_cada_termino <- table(terminos_a_analizar)
  terminos_a_analizar       <- unique(terminos_a_analizar)
  pvalue                    <- setNames(rep(1, length(terminos_a_analizar)), terminos_a_analizar)
  for(j in seq_along(terminos_a_analizar)){
    n_en_el_termino           <- sum(go_completo$go_id==terminos_a_analizar[j])
    pvalue[j]                 <- phyper(q=n_blancas_en_cada_termino[terminos_a_analizar[j]]-1, 
                                        m=n_blancas, 
                                        n=n_negras,
                                        k = n_en_el_termino, 
                                        lower.tail = FALSE)
  }
  return(pvalue)
}), genes)

sobrerepresentacion_ajustada   <- lapply(sobrerepresentacion, p.adjust, method="fdr")
sobrerepresentacion_de_modulos <- Filter(length, lapply(sobrerepresentacion_ajustada, function(x){goterms[names(which(x < 0.05))]}))
