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

g_ego <- make_ego_graph(g, 1, V(g)[names(V(g)) == "AT2G21660"])[[1]]
a_ego <- as_adjacency_matrix(g_ego, sparse = F)
a_ego <- melt(a_ego)
a_ego <- a_ego[a_ego[, 3] == 1, 1:2]
a_ego <- apply(a_ego, 1, function(s){return(paste(sort(s), collapse = ","))})

#Calculamos la correlación entre bines y proteinas relacionadas con splicing
perfiles_bines  <- perfiles_bines[rownames(perfiles_bines) %in% names(V(g_ego)), ]
perfiles_genes  <- perfiles_genes[c(sample(rownames(perfiles_genes), 0), intersect(rownames(perfiles_genes), names(V(g_ego)))), ]
sd <- apply(perfiles_genes, 1, sd)
m <- apply(perfiles_genes, 1, mean)
z <- (perfiles_genes - m)/sd

sd <- apply(perfiles_bines, 1, sd)
m <- apply(perfiles_bines, 1, mean)
zb <- (perfiles_bines - m)/sd

proteinas_relacionadas_con_splicing <- intersect(rownames(perfiles_genes), proteinas_relacionadas_con_splicing)
ft                                  <- intersect(rownames(perfiles_genes), reguladores_de_expresion)
#ccov <- cor(t(rbind(perfiles_bines, perfiles_genes)))
#ccov <- cor(t(perfiles_genes))
ccov <- cov(t(rbind(z, zb)))
#ccov[proteinas_relacionadas_con_splicing, proteinas_relacionadas_con_splicing] <- 0
#ccov <- cor(t(perfiles_genes[sample(rownames(perfiles_genes), 30), -1]))
#diag(ccor)      <- 0

#est <- huge(t(rbind(z, zb)), method = "glasso", nlambda = 100)

#est <- huge(ccov, method = "glasso", nlambda = 100)
#sel <- huge.select(est, criterion = "stars", rep.num = 10, stars.thresh = 0.05)
#plot(sel)

#sel$opt.lambda
#sel$variability
#Construyo la matriz de penalizaciones
rho           <- matrix(ncol=ncol(ccov), nrow=nrow(ccov), 1)
colnames(rho) <- colnames(ccov)
rownames(rho) <- rownames(ccov)

#rho[proteinas_relacionadas_con_splicing, proteinas_relacionadas_con_splicing] <- 1
#rho[setdiff(rownames(perfiles_genes), proteinas_relacionadas_con_splicing), intersect(rownames(ccov), rownames(perfiles_bines))] <- 1
#rho[intersect(rownames(ccov), rownames(perfiles_bines)), setdiff(rownames(perfiles_genes), proteinas_relacionadas_con_splicing)] <- 1


# ceros         <- cbind(rep(which(colnames(ccov) %in% setdiff(rownames(perfiles_genes), proteinas_relacionadas_con_splicing)), each=nrow(perfiles_bines)), 
#                        which(colnames(ccov) %in% rownames(perfiles_bines)))
# 
# ceros         <- rbind(ceros, 
#                        cbind(rep(which(colnames(ccov) %in% setdiff(rownames(perfiles_genes), proteinas_relacionadas_con_splicing)), each=nrow(perfiles_bines)), 
#                              which(colnames(ccov) %in% rownames(perfiles_bines)))
nparametros <- ncol(ccov)^2
pasos       <- 100
bic         <- matrix(0, ncol = 6, nrow = pasos^3)
#path        <- vector("list", length=nrow(bic))
#bic         <- matrix(0, ncol = 4, nrow = pasos)
i           <- 1

for(penalizacion_gg in seq(0.1, 1, length.out = pasos)){
  for(penalizacion_bb in seq(0.1, 1, length.out = pasos)){
    print(paste0("Barrido ", i, " de ", nrow(bic)))    
    start       <- "cold"
    for(penalizacion_gb in seq(0.1, 1, length.out = pasos)){
      #print(penalizacion_gb)
      #print(paste0("Barrido ", i, " de ", nrow(bic), " con gg=", penalizacion_gg, ", bb=", penalizacion_bb, " y gb=", penalizacion_gb))
      penalizacion_gg <- 0.9272727 
      penalizacion_bb <- 0.9090909
      penalizacion_gb <- 0.9454545
          
      rho[ft, rownames(perfiles_genes)]                                  <- penalizacion_gg
      rho[rownames(perfiles_genes), ft]                                  <- penalizacion_gg

      rho[rownames(perfiles_bines), rownames(perfiles_bines)]            <- penalizacion_bb
  
      rho[proteinas_relacionadas_con_splicing, rownames(perfiles_bines)] <- penalizacion_gb
      rho[rownames(perfiles_bines), proteinas_relacionadas_con_splicing] <- penalizacion_gb
      #rho <- penalizacion_gb
      #me <- glasso(ccov, rho=rho, thr = 1e-4, zero = ceros, start = start, w.init = me$w, wi.init = me$wi)
      me <- glasso(ccov, rho=rho, thr = 1e-4, start = start, w.init = me$w, wi.init = me$wi)
      #p_off_d <- sum(me$wi[upper.tri(me$wi)] != 0)
      #for(j in c(1e-8, 1e-9, 1e-10)){
      #me <- glasso(ccov, rho=rho, thr = 1e-4)
      p_off_d   <- sum(me$wi!=0 & col(ccov)<row(ccov))
      path <- ifelse(me$wi!=0 & row(me$wi)!=col(me$wi), 1, 0)
      colnames(path) <- rownames(path) <- rownames(ccov)
      m_path <- melt(path)
      m_path <- m_path[m_path[, 3] == 1, 1:2]
      m_path <- apply(m_path, 1, function(s){return(paste(sort(s), collapse = ","))})
      
      bic[i, ]  <- c(-2*(me$loglik) + p_off_d*log(nparametros), sum(a_ego %in% m_path), 
                     sum(!m_path %in% a_ego), penalizacion_gg, penalizacion_bb, penalizacion_gb)
      #bic[i, ]  <- c(rho, p_off_d*log(nparametros)-2*(me$loglik), p_off_d, me$loglik)
      i <- i + 1
      #}
      #start <- "warm"
    }
  }
}
#plot(bic[, 1], bic[, 2]/bic[, 3])
plot(bic[, 3], bic[, 2])
hist(bic[bic[, 3] == 0, 2])
bic[which(bic[, 2] == 312 & bic[, 3] == 0), ]
abline(v=20, col="red")
plot(setGraphPlotOptions(graph_from_adjacency_matrix(path, mode="undirected"), nombres = T))

0.9368421 0.8526316 0.9157895

seleccionados <- bic[, 2]/bic[, 3] > 2
#rainbow(1000)
points(bic[seleccionados, 3], bic[seleccionados, 2], col="red")

plot(bic[, 2], bic[, 1])
plot(bic[, 1], bic[, 3])
plot(bic[, 1], bic[, 4])
#  }
#}
#Grabo la matriz calculada
save(bic, file="pipeline_archivos/bic.Rdata")
library(plot3D)
pairs(bic)
plot3d(bic[, 2], bic[, 3], bic[, 1])
#Calculo glasso
A <- ifelse(me$wi!=0 & row(me$wi)!=col(me$wi),1,0)

g <- igraph::simplify(graph_from_adjacency_matrix(A, mode = "undirected", weighted = NULL), remove.multiple = T, remove.loops = T)


rownames(me$wi) <- rownames(rho)
colnames(me$wi) <- colnames(rho)

A <- ifelse(me$wi!=0 & row(me$wi)!=col(me$wi),1,0)

g <- igraph::simplify(graph_from_adjacency_matrix(A, mode = "undirected", weighted = NULL), remove.multiple = T, remove.loops = T)
plot(g)
V(g)$type              <- "target"
V(g)[names(V(g)) %in% reguladores_seleccionados]$type <- "reg"
V(red_ego_g)$type              <- "target"
V(red_ego_g)[names(V(red_ego_g)) %in% reguladores_seleccionados]$type <- "reg"

plot(setGraphPlotOptions(red_ego_g, ecolor="gray", nombres="reg"))
plot(setGraphPlotOptions(g, ecolor="gray", nombres="reg"))

#Grabo la matriz calculada
save(me, file=paste("me_", lambda, "_", factor_penalizacion_no_reguladores, ".Rdata", sep=""))

# NOT RUN {
library(rgl)
open3d()
x <- sort(rnorm(1000))
y <- rnorm(1000)
z <- rnorm(1000) + atan2(x, y)
plot3d(x, y, z, col = rainbow(1000))
# }
