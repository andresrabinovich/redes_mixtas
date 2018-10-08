##############################################################################################
#SCRIPT QUE GENERA UNA RED DE REGULACION ENTRE BINES Y GENES
#Autor: Andrés Rabinovich
#Creación: 24/07/2018
#Última modificación: 24/07/2018 (por Andrés Rabinovich)
##############################################################################################

#Librerías que necesita
library(igraph)
library(limma) #Para strsplit2
library(shiny)
library(reshape2) #Para meltear la matrix de cor
setwd("~/doctorado/programacion/redes_mixtas/")

#Funciones varias para grafos
source("pipeline/funciones_grafos.R")

#Elegimos las condiciones para la que vamos a generar la red
(load("pipeline_archivos/1_seleccion_de_condiciones.Rdata"))

#Cargamos archivos de reguladores
(load("pipeline_archivos/reguladores.Rdata"))

#Levantamos las redes de genes y bines y los perfiles de ambos para armar la red mixta
(load("pipeline_archivos/2_genes_prefiltrados.Rdata"))
(load("pipeline_archivos/3_transcriptoma.Rdata"))
(load("pipeline_archivos/4_bines_prefiltrados.Rdata"))
(load("pipeline_archivos/5_spliceoma.Rdata"))
(load("pipeline_archivos/6_red_mixta.Rdata"))


#Calculamos la correlación entre bines de la red spliceoma y proteinas en general que no estén en la red transcriptoma
perfiles_bines       <- perfiles_bines[names(V(g_bines)), ]
perfiles_tf          <- perfiles_genes[rownames(perfiles_genes) %in% intersect(names(V(g_genes)), reguladores_de_expresion), ]
perfiles_genes       <- perfiles_genes[!rownames(perfiles_genes) %in% names(V(g_genes)), ]
ccor                 <- cor(t(perfiles_bines), t(perfiles_genes))
adyacencia_bg        <- setNames(melt(ccor), c('bines', 'genes', 'pesos'))
adyacencia_bg        <- adyacencia_bg[strsplit2(adyacencia_bg$bines, ":")[, 1] == adyacencia_bg$genes, ]

ccor            <- abs(cor(t(perfiles_tf), t(perfiles_genes)))
adyacencia_ftg  <- setNames(melt(ccor), c('ft', 'genes', 'pesos'))

#interfaz de shiny
ui <- fluidPage(
  titlePanel("Red regulación de expresión por splicing"),
  
  sidebarLayout(
    sidebarPanel(
      
      sliderInput("correlacion_minima", 
                  label = "Correlación mínima:",
                  min = 0.65, max = 1, value=0.95, step=0.01),
      #sliderInput("correlacion_minima_ft", 
      #            label = "Correlación mínima FT:",
      #            min = 0.65, max = 1, value=0.81, step=0.01),
      br(),
      actionButton("siguiente", label = "Siguiente")
    ),
    
    mainPanel(
      plotOutput("cantidad_de_enlaces_gen_bin", height = "250px"),
      plotOutput("cantidad_de_enlaces_gen_ft", height = "250px"),
      plotOutput("cantidad_de_genes", height = "250px")
    )
  )
)

#Server de shiny para filtrar y actualizar los datos
server <- function(input, output) {
  
  #Inicializa las variables reactive para que sean accesibles desde otros elementos
  rv <- reactiveValues(bines_con_uso_diferencial=NULL) 
  
  #Plotea el gráfico de cantidad de enlaces gen bin en función de la correlación
  output$cantidad_de_enlaces_gen_bin <- renderPlot({
    enlaces <- c()
    correlaciones <- signif(seq(0.65, 1, 0.01), 2)
    for(correlacion_minima in correlaciones){
      enlaces <- c(enlaces, sum(adyacencia_bg$pesos > correlacion_minima))
    }
    plot(correlaciones, enlaces, xlab="Correlación", ylab="Enlaces bin - gen", cex.lab=1.5)
    grid(col = "gray", lty = "dotted")
    cantidad_de_enlaces_elegidos <- enlaces[correlaciones == input$correlacion_minima]
    text(x=0.95, y=0.9*max(enlaces), labels=paste0("(", input$correlacion_minima, ", ", cantidad_de_enlaces_elegidos, ")"), cex=0.9, pos=2)
    abline(v=input$correlacion_minima, col="red")
  })
  
  #Plotea el gráfico de cantidad de enlaces gen bin en función de la correlación
  output$cantidad_de_enlaces_gen_ft <- renderPlot({
    enlaces <- c()
    correlaciones <- signif(seq(0.65, 1, 0.01), 2)
    for(correlacion_minima in correlaciones){
      enlaces <- c(enlaces, sum(adyacencia_ftg$pesos > correlacion_minima))
    }
    plot(correlaciones, enlaces, xlab="Correlación", ylab="Enlaces FT - gen", cex.lab=1.5)
    grid(col = "gray", lty = "dotted")
    cantidad_de_enlaces_elegidos <- enlaces[correlaciones == input$correlacion_minima]
    text(x=0.95, y=0.9*max(enlaces), labels=paste0("(", input$correlacion_minima, ", ", cantidad_de_enlaces_elegidos, ")"), cex=0.9, pos=2)
    abline(v=input$correlacion_minima, col="red")
  })  
  
  #Plotea el gráfico de cantidad de enlaces gen bin en función de la correlación
  #output$cantidad_de_bines <- renderPlot({
  #  bines <- c()
  #  correlaciones <- signif(seq(0.65, 1, 0.01), 2)
  #  for(correlacion_minima in correlaciones){
  #    bines <- c(bines, length(unique(adyacencia[adyacencia$pesos > correlacion_minima, "bines"])))
  #  }
  #  plot(correlaciones, bines, xlab="Correlación", ylab="Bines", cex.lab=1.5)
  #  grid(col = "gray", lty = "dotted")
  #  cantidad_de_bines_elegidos <- bines[correlaciones == input$correlacion_minima]
  #  text(x=0.7, y=0.8*max(bines), labels=paste0("(", input$correlacion_minima, ", ", cantidad_de_bines_elegidos, ")"), cex=0.9, pos=2)
  #  abline(v=input$correlacion_minima, col="red")
  #})
  
  #Plotea el gráfico de cantidad de enlaces gen bin en función de la correlación
  output$cantidad_de_genes <- renderPlot({
    genes <- c()
    correlaciones <- signif(seq(0.65, 1, 0.01), 2)
    for(correlacion_minima in correlaciones){
      genes <- c(genes, length(unique(adyacencia_ftg[adyacencia_ftg$pesos > correlacion_minima, "genes"]))+
                        length(unique(adyacencia_bg[adyacencia_bg$pesos > correlacion_minima, "genes"])))
    }
    plot(correlaciones, genes, xlab="Correlación", ylab="Genes", cex.lab=1.5)
    grid(col = "gray", lty = "dotted")
    cantidad_de_genes_elegidos <- genes[correlaciones == input$correlacion_minima]
    text(x=0.95, y=0.9*max(genes), labels=paste0("(", input$correlacion_minima, ", ", cantidad_de_genes_elegidos, ")"), cex=0.9, pos=2)
    abline(v=input$correlacion_minima, col="red")
  })
  
  #Agrega los nuevos enlaces a la red
  observeEvent(input$siguiente, {
    enlaces_totales                          <- adyacencia_bg[adyacencia_bg$pesos > input$correlacion_minima, ]
    g_genes_regulados_por_splicing           <- graph_from_edgelist(as.matrix(enlaces_totales[, c("bines", "genes")]), directed = TRUE)
    E(g_genes_regulados_por_splicing)$weight <- enlaces_totales$pesos
    
    enlaces_totales                          <- adyacencia_ftg[adyacencia_ftg$pesos > 0.99, ]#input$correlacion_minima, ]
    g_genes_regulados_por_ft                 <- graph_from_edgelist(as.matrix(enlaces_totales[, c("ft", "genes")]), directed = TRUE)
    E(g_genes_regulados_por_ft)$weight       <- enlaces_totales$pesos
    g_regulacion_completa                    <- g_genes_regulados_por_ft %u% g_genes_regulados_por_splicing
    readme <- data.frame(dato="correlacion_minima",
                         valor=input$correlacion_minima)
    save(g_regulacion_completa, readme, file="pipeline_archivos/8_red_genes_regulados.Rdata")
    stopApp()
  })
}
shinyApp(ui = ui, server = server)