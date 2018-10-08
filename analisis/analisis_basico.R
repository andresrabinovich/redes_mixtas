#Librerías que necesita
library(igraph)
library(limma)
setwd("~/doctorado/programacion/redes_mixtas/")

#Funciones varias para grafos
source("pipeline/funciones_grafos.R")

#Cargamos la lista de reguladores generada a partir de varios repositorios de datos.
(load("pipeline_archivos/reguladores.Rdata"))
(load("pipeline_archivos/3_transcriptoma.Rdata"))

#Nos quedamos con los reguladores de expresión génica y cofactores de transcripcion
reguladores_de_expresion <- unique(c(poi[["TF"]], 
                                     poi[["TCF"]],
                                     reguladores[reguladores[,"tipo_de_regulador"]=="TF","gene_id"])
)

#Levantamos los genes relacionados con splicing en el transcriptoma
proteinas_relacionadas_con_splicing <- unique(c(poi$SP, reguladores$gene_id[reguladores$tipo_de_regulador=="RBP" | reguladores$tipo_de_regulador=="NO CLASIFICADOS"]))

#Los genes que son tanto factores de transcripcion como proteinas de splicing
proteinas_mixtas <- intersect(proteinas_relacionadas_con_splicing, reguladores_de_expresion)

#Levantamos los simbolos de los genes
at_simbolos        <- at_to_symbol <- read.table("pipeline_archivos/at_to_symbol_map", sep = "\t", header = F, stringsAsFactors = F)
at_simbolos        <- at_to_symbol$V2
names(at_simbolos) <- at_to_symbol$V1

V(g_genes)$type                                                  <- "target"
V(g_genes)$type[names(V(g_genes)) %in% reguladores_de_expresion] <- "reg"
V(g_genes)$type[names(V(g_genes)) %in% proteinas_mixtas]         <- "mixta"
table(V(g_genes)$type)
mean(degree(g_genes))

plot(setGraphPlotOptions(g_genes))
cg_genes<-infomap.community(componente_gigante(g_genes))
sort(table(cg_genes$membership), decreasing = T)
cluster <- paste(cg_genes$names[cg_genes$membership == 1], collapse = ", ")

(load("pipeline_archivos/5_spliceoma.Rdata"))
V(g_bines)$type                                                                       <- "bin"
V(g_bines)$type[strsplit2(names(V(g_bines)), ":")[, 1] %in% proteinas_relacionadas_con_splicing]         <- "mixta"
V(g_bines)$name[V(g_bines)$type == "bin"]   <- ""
V(g_bines)$name[V(g_bines)$type == "mixta"] <- at_simbolos[strsplit2(V(g_bines)$name[V(g_bines)$type == "mixta"], ":")[, 1]]
plot(setGraphPlotOptions(g_bines, nombres=T))
length(grep(":I", names(V(g_bines))))
length(grep(":E", names(V(g_bines))))
table(V(g_bines)$type)
V(g_bines)$name[V(g_bines)$type == "bin"] <- ""
names(V(g_bines))[V(g_bines)$type == "mixta"]
mean(degree(g_bines))
cg_bines<-infomap.community(componente_gigante(g_bines))
sort(table(cg_bines$membership), decreasing = T)
cluster <- paste(strsplit2(unique(cg_bines$names[cg_bines$membership == 1]), ":")[, 1], collapse = ", ")

(load("pipeline_archivos/6_red_mixta.Rdata"))
srp_bin <- ends(g, E(g))[grepl(":", ends(g, E(g))[, 2]) & !grepl(":", ends(g, E(g))[, 1]), ]
unique(srp_bin[, 1])
unique(strsplit2(srp_bin[, 2], ":")[, 1])
sort(table(srp_bin[, 1]), decreasing = T)
