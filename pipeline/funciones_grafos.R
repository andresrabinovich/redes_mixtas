library(igraph)

#Genera un grafo a partir de una lista de reguladores y sus conexiones.
#
# Args:
#   mta: lista de reguladores y sus conexiones
#   plot: si TRUE imprime un grafo de la componente gigante
#
# Returns:
#   Un grafo
graphFromMTAs<-function(mta, plot.gigante=TRUE){
  require(igraph)
  
  #Agrega una columna más a cada elemento de la lista con el nombre del regulador al que pertenece de forma auxiliar
  #y la transforma en un data.frame
  laux<-lapply(names(mta),function(x){
    return(data.frame(reg=rep(x,nrow(mta[[x]])),mta[[x]]))
  })
  df    <- do.call("rbind",laux)
  
  #Prepara los dos nodos de cada enlace
  vs    <- union(df[,1],df[,2])
  
  #Marca cual es regulador y cual un target del regulador
  vtype <- ifelse(vs%in%regGX,"reg","target")
  df1 <- data.frame(name=vs,type=vtype)
  rownames(df1)<-vs
  
  #Genera el grafo
  g0 <- graph.data.frame(df,directed=TRUE,vertices=df1)
  g <- simplify(g0)
  
  #Plotea el grafo de la componente gigante
  if(plot.gigante){
    clg <- clusters(g)
    icl<-names(sort(table(clg$membership),decreasing=TRUE)[1])
    giant<-induced_subgraph(g,names(clg$membership)[clg$membership==icl])
    giant<-setGraphPlotOptions(giant,size=c(7,5, 0),ecolor="gray")
    plot(giant)
  }
  return(g)
}

# triangle vertex shape
mytriangle <- function(coords, v=NULL, params) {
  vertex.color <- params("vertex", "color")
  if (length(vertex.color) != 1 && !is.null(v)) {
    vertex.color <- vertex.color[v]
  }
  vertex.size <- 1/200 * params("vertex", "size")
  if (length(vertex.size) != 1 && !is.null(v)) {
    vertex.size <- vertex.size[v]
  }
  vertex.frame.color <- params("vertex", "frame.color")
  if (length(vertex.frame.color) != 1 && !is.null(v)) {
    vertex.frame.color <- vertex.frame.color[v]
  }  
  symbols(x=coords[,1], y=coords[,2], bg=vertex.color, fg=vertex.frame.color,
          stars=cbind(vertex.size, vertex.size, vertex.size),
          add=TRUE, inches=FALSE)
}

# clips as a circle
add_shape("triangle", clip=shapes("circle")$clip, plot=mytriangle)

# star vertex shape
mystar <- function(coords, v=NULL, params) {
  vertex.color <- params("vertex", "color")
  if (length(vertex.color) != 1 && !is.null(v)) {
    vertex.color <- vertex.color[v]
  }
  vertex.size <- 1/200 * params("vertex", "size")
  if (length(vertex.size) != 1 && !is.null(v)) {
    vertex.size <- vertex.size[v]
  }
  vertex.frame.color <- params("vertex", "frame.color")
  if (length(vertex.frame.color) != 1 && !is.null(v)) {
    vertex.frame.color <- vertex.frame.color[v]
  }
  symbols(x=coords[,1], y=coords[,2], bg=vertex.color, fg=vertex.frame.color,
          stars=cbind(vertex.size, vertex.size, vertex.size, vertex.size, vertex.size, vertex.size),
          add=TRUE, inches=FALSE)
}

# clips as a circle
add_shape("star", clip=shapes("circle")$clip, plot=mystar)

setGraphPlotOptions<-function(gg,size=c(reg=3,target=3,bin=4), ecolor=rgb(.7,.7,.7,.9), alfa=0.5, nombres=FALSE, lcolor="black", lsize=0.8, ego=NULL, sacar_desconectados = FALSE, graficar_grado = FALSE, clusters = NULL){
  
  if(is.null(V(gg)$label)) V(gg)$label <- names(V(gg))
  if(nombres == FALSE) {
    V(gg)$label <- NA
  }else if(nombres == "reg"){
    for(i in 1:vcount(gg)){
      if(V(gg)[i]$type!="reg") V(gg)[i]$label <- ""
    }
  }
  if(is.null(V(gg)$type)) V(gg)$type <- "target"
  V(gg)$label.color <- lcolor
  V(gg)$label.cex   <- lsize
  #V(gg)$size  <- ifelse(V(gg)$type=="reg",size[1],size[2])
  names(size) <- c("reg", "target", "bin")
  V(gg)$size <- size[V(gg)$type]
  V(gg)$frame.color <- "white"
  
  formas <- c(reg="square", target="circle", bin_reg="triangle", bin="star")
  colores <- c(reg=rgb(0,.6,0,alfa), target=rgb(.6,0,0,alfa), bin=rgb(0,.5, 0,alfa))
  #names(size) <- names(colores) <- names(formas) <- c("reg", "target", "bin")
  V(gg)$shape<-formas[V(gg)$type]
  V(gg)$color<-colores[V(gg)$type]
  E(gg)$curved      <- 0.2
  E(gg)$color       <- ecolor
  E(gg)$arrow.size  <- 0.1
  if(!is.null(E(gg)$weight)){ #Alfa de edges entre 0.2 y 1
      alfa = (E(gg)$weight-min(E(gg)$weight))/(max(E(gg)$weight)-min(E(gg)$weight)) * (1-0.2) + 0.2
      #E(gg)$color <- rgb(.7,.7,.7, alfa)
  }
  if(sacar_desconectados){
    V(gg)[degree(gg) == 0]$color <- "white"
    V(gg)[degree(gg) == 0]$label.color <- "white"
  }
  #if(is.null(pesos)){
  #E(gg)$color <- ecolor
  #}
  if(graficar_grado){
    V(gg)$size  <- ifelse(V(gg)$type=="reg",round(log2(igraph::degree(gg)+1)), round(log2(0.5*igraph::degree(gg)+1)))
    #V(gg)$size  <- round(log2(igraph::degree(gg)+1))
  }
  if(!is.null(clusters)){
    #colores <- c()
    #for(i in 1:nrow(brewer.pal.info)){
    #  colores <- c(colores, brewer.pal(brewer.pal.info[i, 1], rownames(brewer.pal.info)[i]))
    #}
    V(gg)$frame.color <- rgb(0,0,0,1)
    V(gg)$color <- cg$membership+1
  }
  if(!is.null(ego)) {
    #V(gg)[1]$shape <- "csquare"
    V(gg)[ego]$color <- rgb(0,1,0,alfa)
  }
  return(gg) 
}


red_ego<-function(g, nodos = TRUE, unica=TRUE){
  vertices <- ego(g, nodes = nodos)
  redes <- lapply(vertices, function(x){
    return(induced_subgraph(g, x))
  })
  if(unica){
    red <- redes[[1]]
    if(length(redes) > 1){
      for(i in 2:length(redes)){
        red <- red %u% redes[[i]]
      }
    }
    redes <- simplify(red)
  }
  return(redes)
}

componente_gigante<-function(g){
  componentes <- components(g)$membership
  cg <- which.max(table(componentes))
  return(induced_subgraph(g, vids = V(g)[componentes == cg]))
}

sacar_nodo<-function(g, nodo){
  if(is.numeric(nodo)){
    return(induced_subgraph(g, vids = V(g)[-nodo]))
  }else{
    nodo = which(names(V(g)) == nodo)
    return(induced_subgraph(g, vids = V(g)[-nodo]))
  }
}

#Función para normalizar  al plotear
normalizar <- function(x){
  if(max(x) == min(x)) return(x)
  return((x-min(x))/(max(x)-min(x)))
}

#Función para estandarizar
estandarizar <- function(x){
  return((x-mean(x))/sd(x))
}