##############################################################################################
#ANALISIS DE RED MIXTA
#Autor: Andrés Rabinovich
#Creación: 18/06/2018
#Última modificación: XX/XX/XXXX (por XXX)
##############################################################################################

#Librerías que necesita
library(limma)
library(igraph)
setwd("/home/arabinov/doctorado/programacion/redes_mixtas/")

#Funciones varias para grafos
source("pipeline/funciones_grafos.R")

#Levantamos las redes
(load("pipeline_archivos/2_genes_prefiltrados.Rdata"))
(load("pipeline_archivos/3_transcriptoma.Rdata"))
(load("pipeline_archivos/4_bines_prefiltrados.Rdata"))
(load("pipeline_archivos/5_spliceoma.Rdata"))
(load("pipeline_archivos/6_red_mixta.Rdata"))

#Perfiles 
perfiles <- rbind(perfiles_genes, perfiles_bines)

#Estadísticas de redes
vcount(g_genes)
ecount(g_genes)

vcount(g_bines)
ecount(g_bines)

vcount(g)
ecount(g)

#Encuentra los enlaces intrared (bin contra gen)
enlaces                <- ends(g, E(g))
enlaces_intrared       <- setdiff(grep(":", enlaces[, 2]), grep(":", enlaces[, 1]))
enlaces_interred_bines <- intersect(grep(":", enlaces[, 1]), grep(":", enlaces[, 2]))
enlaces_interred_genes <- setdiff(1:ecount(g), c(enlaces_intrared, enlaces_interred_bines))
enlaces_intrared       <- enlaces[enlaces_intrared, ]
enlaces_interred_bines <- enlaces[enlaces_interred_bines, ]
enlaces_interred_genes <- enlaces[enlaces_interred_genes, ]

#Moduladores
l_enlaces_intrared <- setNames(lapply(unique(enlaces_intrared[, 1]), function(x){
  return(enlaces_intrared[which(enlaces_intrared[, 1] == x), 2])
}), unique(enlaces_intrared[, 1]))
layout(matrix(1:16, ncol=4, nrow=4))
for(m in names(l_enlaces_intrared)){
  perfiles   <- rbind(perfiles_genes[m, ], perfiles_bines[l_enlaces_intrared[[m]], ])
  perfiles_e <- t(apply(perfiles, 1, estandarizar))
  matplot(t(perfiles_e), type="l", xlab="Tiempo", ylab="Expresión", main=m)
  matplot(apply(perfiles_e, 2, mean), type = "l", lwd = 3, add = TRUE)
}
layout(1)

#Redes ego a k=1
g_interes <- setNames(make_ego_graph(g, 1, V(g)[which(names(V(g)) %in% enlaces_intrared[, 1])]), names(V(g))[(names(V(g)) %in% enlaces_intrared[, 1])])
layout(matrix(1:16, ncol=4, nrow=4))
for(m in names(g_interes)){
  perfiles_e <- t(apply(perfiles[names(V(g_interes[[m]])), ], 1, estandarizar))
  matplot(t(perfiles_e), type="l", xlab="Tiempo", ylab="Expresión", main=m)
  matplot(apply(perfiles_e, 2, mean), type = "l", lwd = 3, add = TRUE)
}
layout(1)

#Redes ego a k=2
g_interes <- setNames(make_ego_graph(g, 2, V(g)[which(names(V(g)) %in% enlaces_intrared[, 1])]), names(V(g))[(names(V(g)) %in% enlaces_intrared[, 1])])
layout(matrix(1:16, ncol=4, nrow=4))
for(m in names(g_interes)){
  perfiles_e <- t(apply(perfiles[names(V(g_interes[[m]])), ], 1, estandarizar))
  matplot(t(perfiles_e), type="l", xlab="Tiempo", ylab="Expresión", main=m)
  matplot(apply(perfiles_e, 2, mean), type = "l", lwd = 5, add = TRUE)  
}
layout(1)

#Grafico solo los modulos con perfiles prototipicos
perfiles_prototipicos <- data.frame(
                                    modulos = c("AT1G02780", "AT1G09140", "AT4G19660", 
                                                "AT5G03740", "AT1G49590", "AT1G17745", "AT1G71860",
                                                "AT3G09650", "AT2G03510", "AT3G13460", "AT2G07698", 
                                                "AT3G17390", "AT3G54500", "AT5G16180"
                                    ), 
                                    tipos =  c("Crecimiento modulado por reloj", "Circadiano", "Crecimiento modulado por reloj",
                                               "Crecimiento modulado por reloj", "Crecimiento modulado por reloj", "Circadiano",
                                               "Circadiano", "Disminución modulada por reloj", "Crecimiento modulado por reloj",
                                               "Circadiano", "Circadiano", "Circadiano", "Circadiano", "Crecimiento modulado por reloj"
                                    ), stringsAsFactors = F)
perfiles_prototipicos <- perfiles_prototipicos[order(perfiles_prototipicos[, 2]), ]

layout(matrix(1:16, ncol=4, nrow=4, byrow = T))
for(p in 1:nrow(perfiles_prototipicos)){
  m <- perfiles_prototipicos[p, 1]
  perfiles   <- rbind(perfiles_genes[m, ], perfiles_bines[l_enlaces_intrared[[m]], ])
  perfiles_e <- t(apply(perfiles, 1, estandarizar))
  matplot(t(perfiles_e), type="l", xlab="Tiempo", ylab="Expresión", main=perfiles_prototipicos[p, 2])
  matplot(apply(perfiles_e, 2, mean), type = "l", lwd = 3, add = TRUE)
}
layout(1)

