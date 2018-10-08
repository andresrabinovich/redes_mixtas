#!/usr/bin/env Rscript
#args = commandArgs(trailingOnly=TRUE)
#if (length(args)!=2) stop("Se deben pasar dos argumentos.", call.=FALSE)

#Elijo el directorio de trabajo
#setwd("~/doctorado/programacion/rna_seq/viejo/")

#Cargo los datos
#(load("glasso_ariel.Rdata"))

##############################################################################################
#ANALISIS DE REDES MIXTAS EGO
#Autor: Andrés Rabinovich
#Creación: 08/06/2018
#Última modificación: 11/06/2018 (por Andrés Rabinovich)
##############################################################################################

#Librerías que necesita
library(limma)
library(igraph)
library(glasso)
setwd("/home/arabinov/doctorado/programacion/redes_mixtas/")

#Funciones varias para grafos
source("pipeline/funciones_grafos.R")

#Levantamos las redes
(load("pipeline_archivos/reguladores.Rdata"))
(load("pipeline_archivos/2_genes_prefiltrados.Rdata"))
(load("pipeline_archivos/3_transcriptoma.Rdata"))
(load("pipeline_archivos/4_bines_prefiltrados.Rdata"))
(load("pipeline_archivos/5_spliceoma.Rdata"))
(load("pipeline_archivos/6_red_mixta.Rdata"))

#Levantamos los genes relacionados con splicing en el transcriptoma
proteinas_relacionadas_con_splicing <- unique(c(poi$SP, reguladores$gene_id[reguladores$tipo_de_regulador=="RBP" | reguladores$tipo_de_regulador=="NO CLASIFICADOS"]))
#ego <- "AT2G21660"
#g_ego <- make_ego_graph(g, 1, V(g)[names(V(g)) == ego])[[1]]
g_ego <- g
#Calculamos la correlación entre bines y proteinas relacionadas con splicing
perfiles_bines  <- perfiles_bines[rownames(perfiles_bines) %in% names(V(g_ego)), ]
perfiles_genes  <- perfiles_genes[c(sample(rownames(perfiles_genes), 0), intersect(rownames(perfiles_genes), names(V(g_ego)))), ]

proteinas_relacionadas_con_splicing <- intersect(rownames(perfiles_genes), proteinas_relacionadas_con_splicing)
ft                                  <- intersect(rownames(perfiles_genes), reguladores_de_expresion)

sd   <- apply(perfiles_genes, 1, sd)
m    <- apply(perfiles_genes, 1, mean)
z    <- (perfiles_genes - m)/sd
sd   <- apply(perfiles_bines, 1, sd)
m    <- apply(perfiles_bines, 1, mean)
zb   <- (perfiles_bines - m)/sd
#ccov <- cov(t(rbind(z, zb)))
#ccov <- cov(t(z))
#ccov <- cov(t(zb))

inestabilidades <- g_glasso <- penalizaciones <- ccov <- setNames(vector("list", 3), c("gg", "bb", "gb"))

ccov[["gg"]]           <- cov(t(z))
ccov[["bb"]]           <- cov(t(zb))
ccov[["gb"]]           <- cov(t(rbind(z[proteinas_relacionadas_con_splicing, ], zb)))

penalizaciones[["gg"]]                               <- matrix(ncol=ncol(ccov[["gg"]] ), nrow=nrow(ccov[["gg"]]), 1)
colnames(penalizaciones[["gg"]])                     <- colnames(ccov[["gg"]])
rownames(penalizaciones[["gg"]])                     <- rownames(ccov[["gg"]])
penalizaciones[["gg"]][ft, rownames(perfiles_genes)] <- "gg"
penalizaciones[["gg"]][rownames(perfiles_genes), ft] <- "gg"

penalizaciones[["bb"]]                                                     <- matrix(ncol=ncol(ccov[["bb"]] ), nrow=nrow(ccov[["bb"]]), 1)
colnames(penalizaciones[["bb"]])                                           <- colnames(ccov[["bb"]])
rownames(penalizaciones[["bb"]])                                           <- rownames(ccov[["bb"]])
penalizaciones[["bb"]][rownames(perfiles_bines), rownames(perfiles_bines)] <- "bb"

penalizaciones[["gb"]]                                                                <- matrix(ncol=ncol(ccov[["gb"]] ), nrow=nrow(ccov[["gb"]]), 1)
colnames(penalizaciones[["gb"]])                                                      <- colnames(ccov[["gb"]])
rownames(penalizaciones[["gb"]])                                                      <- rownames(ccov[["gb"]])
penalizaciones[["gb"]][proteinas_relacionadas_con_splicing, rownames(perfiles_bines)] <- "gb"
penalizaciones[["gb"]][rownames(perfiles_bines), proteinas_relacionadas_con_splicing] <- "gb"

stars <- function(ccov, penalizaciones, pasos = 1000, subsample_ratio = NULL, beta = 0.1, numero_de_repeticiones = 20){
  n <- nrow(ccov)  
  if(is.null(subsample_ratio)){
    if(n > 144){
      subsample_ratio = 10*sqrt(n)/n
    }else{
      subsample_ratio = 0.8
    }
  }
  b = round(subsample_ratio*n)
  paths <- vector("list", length=numero_de_repeticiones)
  for(i in 1:numero_de_repeticiones){
    print(paste0("Subsample ", i, " de ", numero_de_repeticiones))
    subsample <- sample(1:n, replace = FALSE, size = b)
    paths[[i]] <- barrido_glasso(ccov[subsample, subsample], penalizaciones[subsample, subsample], pasos, verbose = T)
    paths[[i]] <- lapply(paths[[i]], function(p){
      completa <- matrix(0, ncol=n, nrow=n)
      completa[subsample, subsample] <- p
      return(completa)
    })
  }
  inestabilidades <- vector("list", length=pasos)
  for(i in 1:pasos){
    tita <- matrix(0, ncol=n, nrow=n)
    for(j in 1:numero_de_repeticiones){
      tita <- tita + paths[[j]][[i]]
    }
    tita <- tita/numero_de_repeticiones
    seda <- 2*tita*(1-tita)
    inestabilidades[[i]] <- sum(seda[upper.tri(seda)])/choose(b, 2)
  }
  inestabilidades <- unlist(inestabilidades)
  inestabilidades_monotonizado <- inestabilidades
  for(i in length(inestabilidades_monotonizado):1){
    inestabilidades_monotonizado[i] <- max(inestabilidades_monotonizado[i:length(inestabilidades_monotonizado)])
  }
  return(inestabilidades_monotonizado)
}

#paths <- barrido_glasso(ccov, penalizaciones, pasos=20, verbose = T)
#densidades <- 0.5*unlist(lapply(paths, sum))/choose(nrow(ccov), 2)
#plot(hexbin(densidades, inestabilidades))

# barrido_glasso_triple <- function(ccov, penalizaciones, pasos = 10, verbose = TRUE){
#   rho                        <- matrix(ncol=ncol(ccov), nrow=nrow(ccov), 1)
#   rho[penalizaciones == "1"] <- 1
#   path <- vector("list", length=pasos^3)
#   i = 1
#   for(penalizacion_gg in seq(0.1, 1, length.out = pasos)){
#     if(verbose == TRUE){
#       print(paste0("Barrido ", i, " de ", pasos^3))        
#     }
#     for(penalizacion_bb in seq(0.1, 1, length.out = pasos)){
#       start       <- "cold"
#       for(penalizacion_gb in seq(0.1, 1, length.out = pasos)){
#         rho[penalizaciones == "gg"] <- penalizacion_gg
#         rho[penalizaciones == "bb"] <- penalizacion_bb
#         rho[penalizaciones == "gb"] <- penalizacion_gb
#         me        <- glasso(ccov, rho=rho, thr = 1e-4, start = start, w.init = me$w, wi.init = me$wi)
#         path[[i]] <- ifelse(me$wi!=0 & row(me$wi)!=col(me$wi), 1, 0)
#         start     <- "warm"
#         i = i + 1
#       }
#     }
#     
#   }
#   return(path)
# }

barrido_glasso <- function(ccov, penalizaciones, pasos = 1000, verbose = TRUE){
  rho                        <- matrix(ncol=ncol(ccov), nrow=nrow(ccov), 1)
  rho[penalizaciones == "1"] <- 1
  path <- vector("list", length=pasos)
  i = 1
  start       <- "cold"  
  for(penalizacion in seq(0.1, 1, length.out = length(path))){
    if(verbose == TRUE & i %% 1000 == 0){
      print(paste0("Barrido ", i, " de ", pasos))        
    }
    rho[penalizaciones != "1"] <- penalizacion
    me        <- glasso(ccov, rho=rho, thr = 1e-4, start = start, w.init = me$w, wi.init = me$wi)
    path[[i]] <- ifelse(me$wi!=0 & row(me$wi)!=col(me$wi), 1, 0)
    start     <- "warm"
    i = i + 1
  }
  return(path)
}

g_glasso <- setNames(vector("list", 3), c("gg", "bb", "gb"))
for(tipo in c("gg", "gb", "bb")){
  print(paste("Barriendo", tipo))
  #inestabilidades[[tipo]] <- stars(ccov[[tipo]], penalizaciones[[tipo]], pasos = 1000, numero_de_repeticiones = 10)
  lambdas <- seq(0.1, 1, length.out = length(inestabilidades[[tipo]]))
  plot(lambdas, inestabilidades[[tipo]])
  lambda <- which.max(inestabilidades[[tipo]] < 0.01)
  penalizaciones[[tipo]][penalizaciones[[tipo]] != "1"] <- lambdas[lambda]
  
  me <- glasso(ccov[[tipo]], penalizaciones[[tipo]])
  path <- ifelse(me$wi!=0 & row(me$wi)!=col(me$wi), 1, 0)
  colnames(path) <- rownames(path) <- rownames(ccov[[tipo]])
  g_glasso[[tipo]] <- graph_from_adjacency_matrix(path, mode="undirected")
  plot(setGraphPlotOptions(g_glasso[[tipo]], nombres = T))
}
g_glasso <- componente_gigante(Reduce(union, g_glasso))
g_glasso 
intersection(g_glasso, g_ego)
layout_g                 <- layout_nicely(g_glasso)
rownames(layout_g)       <- names(V(g_glasso))
layout_g[intersect(rownames(layout_g), names(V(g_genes))), 1] <- 10*normalizar(layout_g[intersect(rownames(layout_g), names(V(g_genes))), 1])
layout_g[intersect(rownames(layout_g), names(V(g_genes))), 2] <- 10*normalizar(layout_g[intersect(rownames(layout_g), names(V(g_genes))), 2])
layout_g[intersect(rownames(layout_g), names(V(g_bines))), 1] <- 10*normalizar(layout_g[intersect(rownames(layout_g), names(V(g_bines))), 1])
layout_g[intersect(rownames(layout_g), names(V(g_bines))), 2] <- 10*normalizar(layout_g[intersect(rownames(layout_g), names(V(g_bines))), 2])
layout_g[intersect(rownames(layout_g), names(V(g_genes))), 2] <- layout_g[intersect(rownames(layout_g), names(V(g_genes))), 2] + 20
layout_g[intersect(rownames(layout_g), proteinas_relacionadas_con_splicing), 2]         <- 0.5*(min(layout_g[intersect(rownames(layout_g), names(V(g_genes))), 2])-max(layout_g[intersect(rownames(layout_g), names(V(g_bines))), 2])) + max(layout_g[intersect(rownames(layout_g), names(V(g_bines))), 2])
g_interes_plot <- setGraphPlotOptions(g_glasso, nombres = T)
V(g_interes_plot)$color[names(V(g_interes_plot)) %in% proteinas_relacionadas_con_splicing]                                     <- rgb(1, 0, 0, 0.5)
nombres <- strsplit2(names(V(g_interes_plot)), ":")
#V(g_interes_plot)$label <- paste0(names(V(g_interes_plot)), " (", at_simbolos[nombres[, 1]], ")")

#Encuentra los enlaces inter e intra redes y los grafica con otros colores
enlaces                <- ends(g_interes_plot, E(g_interes_plot))
enlaces_interred       <- setdiff(grep(":", enlaces[, 2]), grep(":", enlaces[, 1]))
enlaces_intrared_bines <- intersect(grep(":", enlaces[, 1]), grep(":", enlaces[, 2]))
enlaces_intrared_genes <- setdiff(1:ecount(g_interes_plot), c(enlaces_interred, enlaces_intrared_bines))
E(g_interes_plot)$color[enlaces_interred]       <- rgb(0.5, 0, 0, 0.5)
E(g_interes_plot)$color[enlaces_intrared_bines] <- rgb(0, 0.5, 0, 0.5)
E(g_interes_plot)$color[enlaces_intrared_genes] <- rgb(0, 0, 0.5, 0.5)

V(g_interes_plot)$shape[names(V(g_interes_plot)) %in% reguladores_de_expresion]  <- "square"
V(g_interes_plot)$shape[!names(V(g_interes_plot)) %in% reguladores_de_expresion] <- "circle"
V(g_interes_plot)$shape[grep(":", names(V(g_interes_plot)))]                     <- "triangle"

V(g_interes_plot)$size        <- 8
V(g_interes_plot)$frame.color <- "black"

plot(g_interes_plot, layout=layout_g)   
legend(-1.4,-1.4, legend = c("SRP", "Otras proteínas", "Bines"), col = c(rgb(0.5,0,0,0.5), rgb(0,0.5,0,0.5), rgb(0,0,0.5,0.5)),    pch = c(15,19,17), bty = "n", pt.cex = 2, cex = 1.2, text.col = "black", horiz = F, inset = c(0.1, 0.1, 0.1))   



# 
# lambda <- 4000#which(inestabilidades_monotonizado > 0.78 & inestabilidades_monotonizado < 0.8)[1]
# i = 1
# for(penalizacion_gg in seq(0.1, 1, length.out = pasos)){
#   for(penalizacion_bb in seq(0.1, 1, length.out = pasos)){
#     for(penalizacion_gb in seq(0.1, 1, length.out = pasos)){
#       if(i == lambda){
#         print(paste0(penalizacion_gg, ",", penalizacion_bb, ",", penalizacion_gb))
#         break 
#       }
#       i = i + 1
#     }
#     if(i == lambda){
#       break
#     }
#   }
#   if(i == lambda){
#     break
#   }  
# }
# 
# 
# paths <- barrido_glasso(ccov, penalizaciones, pasos=20, verbose = T)
# 
# rho                         <- matrix(ncol=ncol(ccov), nrow=nrow(ccov), 1)
# rho[penalizaciones == "1"]  <- 1
# rho[penalizaciones == "gg"] <- seq(0.1, 1, length.out = pasos^3)[lambda]
# rho[penalizaciones == "bb"] <- seq(0.1, 1, length.out = pasos^3)[lambda]
# rho[penalizaciones == "gb"] <- 0.1*seq(0.1, 1, length.out = pasos^3)[lambda]
# me   <- glasso(ccov, rho=rho, thr = 1e-4)
# 
# numero_a_letra <- function(num){
#   #letras <- c("a", "b", "c", "d", "e", "f", "g", "h", "i", "j", "k", "l", "m", "n", "o", "p", "q", "r", "s", "t", "u", "v", "w", "x", "y", "z")
#   resultado <- c()
#   #if(num == 0) return("a")
#   while (num > 0) {
#     remainder <- num %% 26
#     resultado <- c(chr(remainder + 97), resultado)
#     num <- (num - remainder) / 26
#   }
#   resultado <- paste(resultado, collapse = "")
#   return(resultado)
# }
# numero_a_letra(28)
# 
# 
# paths <- barrido_glasso(ccov, penalizaciones, pasos=10, verbose = T)
# layout_g <- NULL
# for(i in seq_along(paths)){
#   path <- ifelse(paths[[i]] !=0 & row(paths[[i]])!=col(paths[[i]]), 1, 0)
#   colnames(path) <- rownames(path) <- rownames(ccov)
#   #plot(setGraphPlotOptions(graph_from_adjacency_matrix(path, mode="undirected"), nombres = T))
#   g_interes <- graph_from_adjacency_matrix(path, mode="undirected")
#   
#   if(is.null(layout_g)){
#     layout_g                 <- layout_nicely(g_interes)
#     rownames(layout_g)       <- names(V(g_interes))
#     layout_g[intersect(rownames(layout_g), names(V(g_genes))), 1] <- 10*normalizar(layout_g[intersect(rownames(layout_g), names(V(g_genes))), 1])
#     layout_g[intersect(rownames(layout_g), names(V(g_genes))), 2] <- 10*normalizar(layout_g[intersect(rownames(layout_g), names(V(g_genes))), 2])
#     layout_g[intersect(rownames(layout_g), names(V(g_bines))), 1] <- 10*normalizar(layout_g[intersect(rownames(layout_g), names(V(g_bines))), 1])
#     layout_g[intersect(rownames(layout_g), names(V(g_bines))), 2] <- 10*normalizar(layout_g[intersect(rownames(layout_g), names(V(g_bines))), 2])
#     layout_g[intersect(rownames(layout_g), names(V(g_genes))), 2] <- layout_g[intersect(rownames(layout_g), names(V(g_genes))), 2] + 20
#     layout_g[intersect(rownames(layout_g), ego), 2]               <- 0.5*(min(layout_g[intersect(rownames(layout_g), names(V(g_genes))), 2])-max(layout_g[intersect(rownames(layout_g), names(V(g_bines))), 2])) + max(layout_g[intersect(rownames(layout_g), names(V(g_bines))), 2])
#   }
#   g_interes_plot <- setGraphPlotOptions(g_interes, nombres = T)
#   #V(g_interes_plot)$color[strsplit2(names(V(g_interes_plot)), ":")[, 1] %in% genes_de_interes_marcelo] <- rgb(0, 0, 1, 0.5)
#   V(g_interes_plot)$color[names(V(g_interes_plot)) %in% ego] <- rgb(1, 0, 0, 0.5)
#   nombres <- strsplit2(names(V(g_interes_plot)), ":")
#   V(g_interes_plot)$label <- paste0(names(V(g_interes_plot)), " (", at_simbolos[nombres[, 1]], ")")
#   #  enlaces_interred
#   
#   #Encuentra los enlaces inter e intra redes y los grafica con otros colores
#   enlaces                <- ends(g_interes_plot, E(g_interes_plot))
#   enlaces_interred       <- setdiff(grep(":", enlaces[, 2]), grep(":", enlaces[, 1]))
#   enlaces_intrared_bines <- intersect(grep(":", enlaces[, 1]), grep(":", enlaces[, 2]))
#   enlaces_intrared_genes <- setdiff(1:ecount(g_interes_plot), c(enlaces_interred, enlaces_intrared_bines))
#   E(g_interes_plot)$color[enlaces_interred]       <- rgb(0.5, 0, 0, 0.5)
#   E(g_interes_plot)$color[enlaces_intrared_bines] <- rgb(0, 0.5, 0, 0.5)
#   E(g_interes_plot)$color[enlaces_intrared_genes] <- rgb(0, 0, 0.5, 0.5)
#   
#   V(g_interes_plot)$shape[names(V(g_interes_plot)) %in% reguladores_de_expresion]  <- "square"
#   V(g_interes_plot)$shape[!names(V(g_interes_plot)) %in% reguladores_de_expresion] <- "circle"
#   V(g_interes_plot)$shape[grep(":", names(V(g_interes_plot)))]                     <- "triangle"
# 
#   V(g_interes_plot)$size        <- 8
#   V(g_interes_plot)$frame.color <- "black"
#   #nombre_de_archivo <- numero_a_letra(i)
#   tiff(filename = paste0("imagenes/", sprintf("imagen_%05d", i), ".tiff"), type = "cairo")
#   plot(g_interes_plot, layout=layout_g)   
#   legend(-1.4,-1, legend = c(paste("gg:", length(enlaces_intrared_genes)),
#                                paste("gb:", length(enlaces_interred)),
#                                paste("bb:", length(enlaces_intrared_bines)),
#                                 paste("lambda:",  seq(0.1, 1, length.out = pasos^3)[i])))
#   dev.off()
#   if(i %% 100 == 0){
#     print(paste(sprintf("imagen_%05d", i), i))
#     #print(i)
#   }
# }
# 
# system("convert-im6 ~/doctorado/programacion/redes_mixtas/imagenes/*.tiff ~/doctorado/programacion/redes_mixtas/animacion.gif")

# path <- ifelse(me$wi!=0 & row(me$wi)!=col(me$wi), 1, 0)
# colnames(path) <- rownames(path) <- rownames(ccov)
# plot(setGraphPlotOptions(graph_from_adjacency_matrix(path, mode="undirected"), nombres = T))
# g_interes <- graph_from_adjacency_matrix(path, mode="undirected")

# layout_g                 <- layout_nicely(g_glasso)
# rownames(layout_g)       <- names(V(g_glasso))
# layout_g[intersect(rownames(layout_g), names(V(g_genes))), 1] <- 10*normalizar(layout_g[intersect(rownames(layout_g), names(V(g_genes))), 1])
# layout_g[intersect(rownames(layout_g), names(V(g_genes))), 2] <- 10*normalizar(layout_g[intersect(rownames(layout_g), names(V(g_genes))), 2])
# layout_g[intersect(rownames(layout_g), names(V(g_bines))), 1] <- 10*normalizar(layout_g[intersect(rownames(layout_g), names(V(g_bines))), 1])
# layout_g[intersect(rownames(layout_g), names(V(g_bines))), 2] <- 10*normalizar(layout_g[intersect(rownames(layout_g), names(V(g_bines))), 2])
# layout_g[intersect(rownames(layout_g), names(V(g_genes))), 2] <- layout_g[intersect(rownames(layout_g), names(V(g_genes))), 2] + 20
# layout_g[intersect(rownames(layout_g), ego), 2]         <- 0.5*(min(layout_g[intersect(rownames(layout_g), names(V(g_genes))), 2])-max(layout_g[intersect(rownames(layout_g), names(V(g_bines))), 2])) + max(layout_g[intersect(rownames(layout_g), names(V(g_bines))), 2])
# g_interes_plot <- setGraphPlotOptions(g_glasso, nombres = T)
# #V(g_interes_plot)$color[strsplit2(names(V(g_interes_plot)), ":")[, 1] %in% genes_de_interes_marcelo] <- rgb(0, 0, 1, 0.5)
# V(g_interes_plot)$color[names(V(g_interes_plot)) %in% ego]                                     <- rgb(1, 0, 0, 0.5)
# nombres <- strsplit2(names(V(g_interes_plot)), ":")
# V(g_interes_plot)$label <- paste0(names(V(g_interes_plot)), " (", at_simbolos[nombres[, 1]], ")")
# #  enlaces_interred
# 
# #Encuentra los enlaces inter e intra redes y los grafica con otros colores
#  enlaces                <- ends(g_interes_plot, E(g_interes_plot))
#  enlaces_interred       <- setdiff(grep(":", enlaces[, 2]), grep(":", enlaces[, 1]))
#  enlaces_intrared_bines <- intersect(grep(":", enlaces[, 1]), grep(":", enlaces[, 2]))
#  enlaces_intrared_genes <- setdiff(1:ecount(g_interes_plot), c(enlaces_interred, enlaces_intrared_bines))
#  E(g_interes_plot)$color[enlaces_interred]       <- rgb(0.5, 0, 0, 0.5)
#  E(g_interes_plot)$color[enlaces_intrared_bines] <- rgb(0, 0.5, 0, 0.5)
#  E(g_interes_plot)$color[enlaces_intrared_genes] <- rgb(0, 0, 0.5, 0.5)
# #   
#  V(g_interes_plot)$shape[names(V(g_interes_plot)) %in% reguladores_de_expresion]  <- "square"
#  V(g_interes_plot)$shape[!names(V(g_interes_plot)) %in% reguladores_de_expresion] <- "circle"
#  V(g_interes_plot)$shape[grep(":", names(V(g_interes_plot)))]                                <- "triangle"
# #   
#    V(g_interes_plot)$size        <- 8
#    V(g_interes_plot)$frame.color <- "black"
# #   
#    plot(g_interes_plot, layout=layout_g)   
#    legend(-1.4,-1.4, legend = c("SRP", "Otras proteínas", "Bines"), col = c(rgb(0.5,0,0,0.5), rgb(0,0.5,0,0.5), rgb(0,0,0.5,0.5)),    pch = c(15,19,17), bty = "n", pt.cex = 2, cex = 1.2, text.col = "black", horiz = F, inset = c(0.1, 0.1, 0.1))   
# #   }
# 
# 
# 
# intersection(g_ego, g_glasso)


