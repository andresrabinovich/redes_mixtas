##############################################################################################
#Graficos 3d de redes multicapa
#Autor: Andrés Rabinovich
#Creación: 13/07/2018
#Última modificación: 20/07/2018 (por Andrés Rabinovich)
##############################################################################################

#Librerías que necesita
library(limma)
library(igraph)
library(rgl)
setwd("~/doctorado/programacion/redes_mixtas/")

#Funciones varias para grafos
source("pipeline/funciones_grafos.R")

#Levantamos las redes
(load("pipeline_archivos/reguladores.Rdata"))
(load("pipeline_archivos/2_genes_prefiltrados.Rdata"))
(load("pipeline_archivos/3_transcriptoma.Rdata"))
(load("pipeline_archivos/5_spliceoma.Rdata"))
(load("pipeline_archivos/6_red_mixta.Rdata"))
(load("pipeline_archivos/8_red_genes_regulados.Rdata"))

#Levantamos los genes relacionados con splicing en el transcriptoma
proteinas_relacionadas_con_splicing <- unique(c(poi$SP, reguladores$gene_id[reguladores$tipo_de_regulador=="RBP" | reguladores$tipo_de_regulador=="NO CLASIFICADOS"]))

#Levantamos los simbolos de los genes
at_simbolos        <- at_to_symbol <- read.table("pipeline_archivos/at_to_symbol_map", sep = "\t", header = F, stringsAsFactors = F)
at_simbolos        <- at_to_symbol$V2
names(at_simbolos) <- at_to_symbol$V1

#Función para generar enlaces curvos que mejoran la visualización
encontrar_curva <- function(x, y, n){
  if(n < 2) n <- 2
  dy <- diff(y)
  dx <- diff(x)
  
  #Hacemos curvas con radios grandes
  minrad <- (1/2)*sqrt(dx^2 + dy^2)
  R <- 5*minrad
  
  #Resolvemos la ecuación paramétrica del círculo  
  A <- (dy^2 + dx^2)
  B <- 0
  C <- (1/4)*(dx^2 + dy^2) - R^2
  k <- suppressWarnings(as.numeric(polyroot(c(C, B, A))))
  
  #Encontramos el centro del círculo
  mn <- c(mean(x), mean(y))
  cc <- rbind(mn + k[1]*c(dy, -dx),
              mn + k[2]*c(dy, -dx))
  
  colnames(cc) = c("x", "y")
  
  #Buscamos los ángulos del punto inicial y del final
  th1     <- acos((x[1]-cc[1, "x"])/R)
  if(abs(cc[1, "y"]+R*sin(th1) - y[1]) >1e-10) th1 = -th1
  th2     <- acos((x[2]-cc[1, "x"])/R)
  if(abs(cc[1, "y"]+R*sin(th2) - y[2]) >1e-10) th2 = -th2
  
  #Corregimos por ángulos extremos
  if(th1*th2 < 0 & abs(th1) > 0.5*pi & abs(th2) > 0.5*pi){
    th <- seq(max(th1, th2), pi, length.out = floor(0.5*n))
    th <- c(th, seq(-pi,min(th1, th2), length.out = ceiling(0.5*n)))
  }else{
    th <- seq(min(th1, th2), max(th1, th2), length.out = n)
  }
  
  #Generamos los puntos de la curva
  circulo <- cbind(cc[1, "x"]+R*cos(th), cc[1, "y"]+R*sin(th))
  
  #Devolvemos la lista de rectas como inicio y fin
  if(n == 2){
    inicio <- matrix(circulo[-n, ], ncol=2, byrow=T)
    fin    <- matrix(circulo[-1, ], ncol=2, byrow=T)
  }else{
    inicio <- circulo[-n, ]
    fin    <- circulo[-1, ]
  }
  return(cbind(inicio, fin))
}

#Graficamos la red por capas
igraph_plot3d_por_capas <- function(g){
  g <- g_completo
  l_objetos_geometricos <- list("proteina"                                                   = icosahedron3d(),
                                "regulador_de_expresion"                                     = cube3d(),
                                "proteina_relacionada_con_splicing"                          = tetrahedron3d(),
                                "regulador_de_expresion_y_proteina_relacionada_con_splicing" = cuboctahedron3d()
                                )
  
  colores <- setNames(c("black", "red", "green", "blue"), c("proteina", "regulador_de_expresion", "proteina_relacionada_con_splicing", "regulador_de_expresion_y_proteina_relacionada_con_splicing"))

  capas <- unique(V(g)$capa)
  l_capas <- setNames(vector("list", length(capas)), capas)
  for(capa in capas){
    l_capas[[capa]] <- names(V(g))[V(g)$capa == capa]
  }
  
  colores_intercapa <- t(combn(1:length(l_capas), 2))
  rownames(colores_intercapa) <- rainbow(nrow(colores_intercapa))
  enlaces <- ends(g, E(g))
  enlaces_inter_capas <- matrix(ncol=3, nrow=0)
  for(i in 1:(length(capas)-1)){
    for(j in (i+1):length(capas)){
      enlaces_inter_capas <- rbind(enlaces_inter_capas,
                                   cbind(enlaces[enlaces[, 1] %in% l_capas[[i]] & enlaces[, 2] %in% l_capas[[j]] |
                                   enlaces[, 2] %in% l_capas[[i]] & enlaces[, 1] %in% l_capas[[j]], ], rownames(colores_intercapa)[(colores_intercapa[, 1] == i & colores_intercapa[, 2] == j) | (colores_intercapa[, 1] == i & colores_intercapa[, 2] == j)]))
      #                             enlaces[enlaces[, 1] %in% l_capas[[i]] & enlaces[, 2] %in% unlist(l_capas[-i]) |
      #                             enlaces[, 2] %in% l_capas[[i]] & enlaces[, 1] %in% unlist(l_capas[-i]), ])      
    }
  }
  enlaces_inter_capas[, 1:2] <- t(apply(enlaces_inter_capas[, 1:2], 1, sort))
  enlaces_inter_capas <- unique(enlaces_inter_capas)
  colnames(enlaces_inter_capas) <- c("gen1", "gen2", "color")
  rgl.open() # Open a new RGL device
  aspect3d(1,1,1)
  
  rgl.bg(color = "white") # Setup the background color
  ncapa <- 0
  layout_total <- matrix(0, ncol=4, nrow = vcount(g))
  #layout_total  <- data.frame(layout_nicely(g), z=0)
  colnames(layout_total) <- c("x", "y", "z", "capa")
  rownames(layout_total) <- names(V(g))
  
  nill <- lapply(l_capas, function(capa){
    
    sg                    <- names(V(g)) %in% capa
    g_capa                <- induced_subgraph(g, V(g)[sg], impl = "create_from_scratch")
    layout_capa           <- data.frame(layout_nicely(g_capa))
    #layout_capa           <- layout_total[sg, ]
    #Centramos la capa
    layout_capa$X1        <- layout_capa$X1 - mean(layout_capa$X1)
    layout_capa$X2        <- layout_capa$X2 - mean(layout_capa$X2)
    
    #Movemos las capas impares para un costado
    if(length(ncapa) %% 2 == 0){
      #layout_capa$X1 <- layout_capa$X1 + 1*(max(layout_capa$X1) - min(layout_capa$X1))
      #layout_capa$X2 <- layout_capa$X2 + (max(layout_capa$X2) - min(layout_capa$X2))
    }
    colnames(layout_capa) <- c("x", "y")
    rownames(layout_capa) <- names(V(g_capa))
    layout_capa$z         <- ncapa[length(ncapa)]
    layout_capa$capa      <- length(ncapa)
    enlaces               <- ends(g_capa, E(g_capa))
    elementos             <- data.frame(x = layout_capa[enlaces, "x"], 
                                        y = layout_capa[enlaces, "y"], 
                                        z = layout_capa[enlaces, "z"])
    db<-as.matrix(dist(layout_capa))
    curvas <- matrix(ncol=4, nrow=0)
    if(nrow(enlaces) > 0){
      null <- apply(enlaces, 1, function(e){
        p1 <- layout_capa[e[1], c("x", "y")]
        p2 <- layout_capa[e[2], c("x", "y")]
        curvas <<- rbind(curvas, encontrar_curva(c(p1$x, p2$x), c(p1$y, p2$y), n=5))
      })
      curvas <- cbind(cbind(curvas, ncapa[length(ncapa)]), ncapa[length(ncapa)])
      colnames(curvas) <- c("x_inicio", "y_inicio", "x_fin", "y_fin", "z_inicio", "z_fin")
      rgl.lines(c(t(cbind(curvas[, "x_inicio"], curvas[, "x_fin"]))), 
                c(t(cbind(curvas[, "y_inicio"], curvas[, "y_fin"]))),
                c(t(cbind(curvas[, "z_inicio"], curvas[, "z_fin"]))),
                color = "black", alpha=0.2)      
    }
    for(tipo in names(l_objetos_geometricos)){
      keep <- V(g_capa)$tipo == tipo
      if(sum(keep) > 0) {
        shapelist3d(l_objetos_geometricos[[tipo]], layout_capa[keep, "x"], layout_capa[keep, "y"], layout_capa[keep, "z"], size=0.3*min(db[db != 0]), col=colores[tipo], alpha = 1)        
      }
    }
    ncapa <<- c(ncapa, ncapa[length(ncapa)] + max(max(ncapa), 0.5*max(db[db != 0])))
    layout_total[rownames(layout_capa), ] <<- as.matrix(layout_capa)
  })
  #colores <- setNames(c("red", "green", "blue"), c(3, 4, 5))
  #print(layout_total[t(enlaces_inter_capas), "capa"])
  #2*rowMeans(matrix(layout_total[t(enlaces_inter_capas), "capa"], ncol=2, byrow=TRUE))
  for(color in unique(enlaces_inter_capas[, "color"])){
    rgl.lines(layout_total[t(enlaces_inter_capas[enlaces_inter_capas[, "color"] == color, c("gen1", "gen2")]), "x"],
              layout_total[t(enlaces_inter_capas[enlaces_inter_capas[, "color"] == color, c("gen1", "gen2")]), "y"],
              layout_total[t(enlaces_inter_capas[enlaces_inter_capas[, "color"] == color, c("gen1", "gen2")]), "z"],
              color = enlaces_inter_capas[enlaces_inter_capas[, "color"] == color, "color"], alpha=0.1)
  }
  for(n in 1:length(l_capas)){
    rgl.quads( x = rep(c(min(layout_total[layout_total[, "capa"] == n, "x"]), max(layout_total[layout_total[, "capa"] == n, "x"])), each = 2),
             y = c(min(layout_total[, "y"]), max(layout_total[, "y"]), max(layout_total[, "y"]), min(layout_total[, "y"])), 
             z = rep(ncapa[n], times=4), color="black", alpha=0.1)
  }
  
  #while(interactive()){
  #  keep <- identify3d(layout_total[, "x"], layout_total[, "y"], layout_total[, "z"], n=1, plot=F)  
  #  if(!is.null(keep)) rgl.texts(jitter(layout_total[keep, "x"]), jitter(layout_total[keep, "y"]), layout_total[keep, "z"], rownames(layout_total)[keep], color="black" )
  #}
  #rgl.texts(jitter(layout_total[, "x"]), jitter(layout_total[, "y"]), jitter(layout_total[, "z"]), rownames(layout_total) )
  
  if(F){
    center <- c(0.5, 0.5) %*% matrix(par3d("bbox"), 2, 3) 
    prev <- c(0, 0, 0)  # initially, no translation
    while(interactive()) {   
      keep <- identify3d(layout_total[, "x"], layout_total[, "y"], layout_total[, "z"], n=1, plot=F)  
      if(!is.null(keep)){
        pt <- c(layout_total[keep, "x"], layout_total[keep, "y"], layout_total[keep, "z"])
        xlat <- center - pt  # undo automatic centering and center on the point
        userMatrix <- par3d("userMatrix")
        par3d(userMatrix = userMatrix %*% t(translationMatrix(xlat[1] - prev[1], 
                                                              xlat[2] - prev[2], 
                                                              xlat[3] - prev[3])))
        prev <- xlat
      }
    }  
  }
  return(enlaces_inter_capas)
}

# gcor         <- abs(cor(t(perfiles_genes[setdiff(names(V(g_regulacion_completa)), c(names(V(g_bines)), names(V(g_genes)))), ])))
# diag(gcor)   <- 0
# adyacencia   <- setNames(melt(gcor), c('gen1', 'gen2', 'pesos'))
# enlaces      <- adyacencia[adyacencia$pesos > 0.95, c('gen1', 'gen2')] 
# enlaces      <- unique(enlaces)
# enlaces      <- as.character(melt(t(enlaces))$value)
# 
# g_regulacion_completa <- add_edges(g_regulacion_completa, edges = enlaces)
# g_completo <- g %u% as.undirected(g_regulacion_completa)

#g_completo <- g
g_completo <- g_glasso

V(g_completo)$tipo                                                                                                     <- "proteina"
V(g_completo)[names(V(g_completo)) %in% reguladores_de_expresion]$tipo                                                 <- "regulador_de_expresion"
V(g_completo)[names(V(g_completo)) %in% proteinas_relacionadas_con_splicing]$tipo                                      <- "proteina_relacionada_con_splicing"
V(g_completo)[names(V(g_completo)) %in% intersect(reguladores_de_expresion, proteinas_relacionadas_con_splicing)]$tipo <- "regulador_de_expresion_y_proteina_relacionada_con_splicing"
V(g_completo)$capa <- "FT"
V(g_completo)[names(V(g_completo)) %in% names(V(g_bines))]$capa <- "bines"
V(g_completo)[names(V(g_completo)) %in% setdiff(names(V(g_regulacion_completa)),c(names(V(g_genes)), names(V(g_bines))))]$capa <- "genes"  
V(g_completo)[names(V(g_completo)) %in% proteinas_relacionadas_con_splicing]$capa <- "SRP"

enlaces_inter_capas <- igraph_plot3d_por_capas(g_completo)

# a <- (rbind(enlaces_inter_capas[enlaces_inter_capas[, 1] %in% names(V(g_completo))[V(g_completo)$capa == "genes"] & !enlaces_inter_capas[, 2] %in% names(V(g_completo))[V(g_completo)$capa == "genes"], c(1,3)],
# enlaces_inter_capas[enlaces_inter_capas[, 2] %in% names(V(g_completo))[V(g_completo)$capa == "genes"] & !enlaces_inter_capas[, 1] %in% names(V(g_completo))[V(g_completo)$capa == "genes"], c(2,3)]))
# a <- aggregate(color ~ gen1, a, table)
# genes_duales <- a[a$color[, "blue"] > 0 & a$color[, "red"] > 0, 1]
# 
# g_completo
