##############################################################################################
#SCRIPT QUE GENERA UNA RED DE COSPLICING PARA UNA DADA TEMPERATURA USANDO COMO REFERENCIA T22.
#Autor: Andrés Rabinovich en base al script red_de_cosplicing.R de Ariel Chernomoretz
#Creación: 08/06/2018
#Última modificación: 11/06/2018 (por Andrés Rabinovich)
##############################################################################################

#Librerías que necesita
library(limma)
library(igraph)
library(shiny)
setwd("/home/arabinov/doctorado/programacion/redes_cosplicing/pipeline_archivos/")

#Funciones varias para grafos
source("../pipeline/funciones_grafos.R")

#Levantamos las redes
(load("3_transcriptoma.Rdata"))
(load("5_spliceoma.Rdata"))
(load("6_red_mixta.Rdata"))

#Levantamos los simbolos de los genes
at_simbolos        <- at_to_symbol <- read.table("at_to_symbol_map", sep = "\t", header = F, stringsAsFactors = F)
at_simbolos        <- at_to_symbol$V2
names(at_simbolos) <- at_to_symbol$V1

#Encuentra los enlaces intrared (bin contra gen)
enlaces               <- ends(g, E(g))
enlaces_bin_gen       <- enlaces[setdiff(grep(":", enlaces[, 1]), grep(":", enlaces[, 2])), ]
genes_ego             <- enlaces_bin_gen[, 2] 
names(genes_ego)      <- paste0(enlaces_bin_gen[, 2], " (", at_simbolos[enlaces_bin_gen[, 2]], ")")

#interfaz de shiny
ui <- fluidPage(
  titlePanel("Red Ego"),
  
  sidebarLayout(
    sidebarPanel(

      #valores límite para filtrar genes
      selectInput("k_vecinos", 
                  label = "Cuántos vecinos más cercanos: ",
                  choices = 1:3,
                  selected = 1),
  
      #valores límite para filtrar genes
      selectInput("ego", 
                  label = "Centrar la red en: ",
                  choices = genes_ego,
                  selected = 1),
      
      textAreaInput("genes_de_interes_marcelo", 
                    value=paste0(c("AT3G09600",
                                   "AT5G64170", 
                                   "AT3G54500", 
                                   "AT1G09140", 
                                   "AT5G02810", 
                                   "AT5G02840", 
                                   "AT2G46830", 
                                   "AT1G01060", 
                                   "AT2G21660"), collapse=", "),
                    label = "Genes de interes",
                    width = '100%',
                    height = '150px'),
      
      uiOutput("modificar_posicion")      

    ),
    
    mainPanel(
      #htmlOutput("genes_de_interes"),
      plotOutput("red_ego", height = "700px", width = "500px"),
      textOutput("genes_en_la_red"),
      textOutput("bines_en_la_red")
    )
  )
)

#Server de shiny para filtrar y actualizar los datos
server <- function(input, output) {
  
  #Inicializa las variables reactive para que sean accesibles desde otros elementos
  rv <- reactiveValues(bines_con_uso_diferencial=NULL, g_ego=NULL, genes_de_interes_marcelo=NULL) 

  observeEvent(input$genes_de_interes_marcelo, {
    rv$genes_de_interes_marcelo <- input$genes_de_interes_marcelo
  })    

  #Plotea la red ego
  output$red_ego <- renderPlot({
    genes_de_interes_marcelo <- trimws(strsplit2(input$genes_de_interes_marcelo, ","))
    g_interes                <- make_ego_graph(g, input$k_vecinos, V(g)[which(names(V(g)) %in% input$ego)])[[1]]
    rv$g_ego                 <- names(V(g_interes))
    layout_g                 <<- layout_nicely(g_interes)
    rownames(layout_g)       <<- names(V(g_interes))
    layout_g[intersect(rownames(layout_g), names(V(g_genes))), 1] <<- 10*normalizar(layout_g[intersect(rownames(layout_g), names(V(g_genes))), 1])
    layout_g[intersect(rownames(layout_g), names(V(g_genes))), 2] <<- 10*normalizar(layout_g[intersect(rownames(layout_g), names(V(g_genes))), 2])
    layout_g[intersect(rownames(layout_g), names(V(g_bines))), 1] <<- 10*normalizar(layout_g[intersect(rownames(layout_g), names(V(g_bines))), 1])
    layout_g[intersect(rownames(layout_g), names(V(g_bines))), 2] <<- 10*normalizar(layout_g[intersect(rownames(layout_g), names(V(g_bines))), 2])
    layout_g[intersect(rownames(layout_g), names(V(g_genes))), 2] <<- layout_g[intersect(rownames(layout_g), names(V(g_genes))), 2] + 20
    g_interes_plot <<- setGraphPlotOptions(g_interes, nombres = T)
    V(g_interes_plot)$color[strsplit2(names(V(g_interes_plot)), ":")[, 1] %in% genes_de_interes_marcelo] <<- rgb(0, 0, 1, 0.5)
    V(g_interes_plot)$color[names(V(g_interes_plot)) %in% input$ego]                                     <<- rgb(1, 0, 0, 0.5)
    nombres <- strsplit2(names(V(g_interes_plot)), ":")
    V(g_interes_plot)$label <<- paste0(names(V(g_interes_plot)), " (", at_simbolos[nombres[, 1]], ")")
    plot(g_interes_plot, layout=layout_g)    
  })

  output$genes_en_la_red <- renderText({ 
    paste("Genes:", paste0(intersect(rv$g_ego, names(V(g_genes))), collapse = ", "))
  })     

  output$bines_en_la_red <- renderText({ 
    paste("Bines:", paste0(intersect(rv$g_ego, names(V(g_bines))), collapse = ", "))
  })     

}
shinyApp(ui = ui, server = server)

