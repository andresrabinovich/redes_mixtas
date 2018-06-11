##############################################################################################
#SCRIPT QUE GENERA UNA RED DE COSPLICING PARA UNA DADA TEMPERATURA USANDO COMO REFERENCIA T22.
#Autor: Andrés Rabinovich en base al script red_de_cosplicing.R de Ariel Chernomoretz
#Creación: 16/03/2018
#Última modificación: 05/06/2018 (por Andrés Rabinovich)
##############################################################################################

#Librerías que necesita
library(igraph)
library(shiny)
setwd("/home/arabinov/doctorado/programacion/redes_cosplicing/pipeline_archivos/")

#Funciones varias para grafos
source("../pipeline/funciones_grafos.R")

#Elegimos las condiciones para la que vamos a generar la red
(load("1_seleccion_de_condiciones.Rdata"))

#Levantamos las cuentas de los genes 
(load("4_bines_prefiltrados.Rdata"))

#Usamos correlación entre perfiles de bines con uso diferencial para armar la red
correlaciones_bines <- cor(t(perfiles_bines))

#interfaz de shiny
ui <- fluidPage(
  titlePanel("Cospliceoma"),
  
  sidebarLayout(
    sidebarPanel(
      helpText("Filtar bines y armar red splicing en base a: "),
      
      #valores límite para filtrar genes
      selectInput("cantidad_de_tiempos_limite", 
                  label = "Cantidad de puntos temporales con uso diferencial",
                  choices = 1:12,
                  selected = 2),
      
      sliderInput("qvalue_limite", 
                  label = "Máximo qvalue:",
                  min = 0, max = 0.9, value = 0.24, step=0.01),
      
      sliderInput("lfchange_limite", 
                  label = "Mínimo cambio a detectar (%):",
                  min = 0, max = 200, value=50, step=25),
      
      sliderInput("correlacion_minima", 
                  label = "Correlación mínima:",
                  min = 0.5, max = 1, value=0.73, step=0.01),
      
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
      
      br(),
      actionButton("siguiente", label = "Siguiente")
    ),
    
    mainPanel(
      textOutput("cantidad_de_bines"),
      htmlOutput("genes_de_interes"),
      plotOutput("componente_gigante", height = "600px")
    )
  )
)

#Server de shiny para filtrar y actualizar los datos
server <- function(input, output) {
  
  #Inicializa las variables reactive para que sean accesibles desde otros elementos
  rv <- reactiveValues(bines_con_uso_diferencial=NULL, g_bines=NULL) 
  
  #Filtra y actualiza los genes de interes que pasaron el filtro    
  filtrar_bines <- reactive({
    #Pide que el gen haya cambiado al menos en cantidad_de_tiempos_limite tiempos en un lfchange_limite
    bines_con_alto_lfchange <- rownames(lfchange_bines)[apply(lfchange_bines, 1, function(x){
      sum(abs(x) > log2(as.numeric(input$lfchange_limite)/100+1)) >= as.numeric(input$cantidad_de_tiempos_limite)
    })]

    #Pide que el gen tenga un qvalue < qvalue_limite en cantidad_de_tiempos_limite tiempos para validad el fold change
    bines_con_qvalue_bajo  <- rownames(qvalues_bines)[apply(qvalues_bines, 1, function(x){
      sum(x<as.numeric(input$qvalue_limite)) >= as.numeric(input$cantidad_de_tiempos_limite)
    })]
    
    #Me quedo con los bines que pasen ambos filtros
    rv$bines_con_uso_diferencial <- sort(intersect(bines_con_alto_lfchange, bines_con_qvalue_bajo))
  })
  
  #Filtra los bines
  output$cantidad_de_bines <- renderText({ 
    filtrar_bines()
    paste("Total de bines que pasaron el filtro:", length(rv$bines_con_uso_diferencial))
  })

  observeEvent(input$genes_de_interes_marcelo, {
    rv$genes_de_interes_marcelo <- input$genes_de_interes_marcelo 
  })    
  
  #Filtra los genes de interes
  output$genes_de_interes <- renderText({ 
    genes_de_interes_marcelo                       <- trimws(strsplit2(input$genes_de_interes_marcelo, ","))
    genes_de_bines_con_uso_diferencial             <- strsplit2(rv$bines_con_uso_diferencial, ":")[, 1]
    genes_de_interes_marcelo_que_pasaron_el_filtro <- intersect(genes_de_bines_con_uso_diferencial, genes_de_interes_marcelo)
    HTML(paste0("Genes de interés que pasaron el filtro (", length(genes_de_interes_marcelo_que_pasaron_el_filtro), "):</br>", paste0(genes_de_interes_marcelo_que_pasaron_el_filtro, collapse="</br>")))
  })   

  #Plotea el gráfico de cantidad de elementos en la componente gigante en función de la correlación
  output$componente_gigante <- renderPlot({
    genes_de_interes_marcelo           <- trimws(strsplit2(input$genes_de_interes_marcelo, ","))
    enlaces                            <- c()
    vertices                           <- c()
    cantidad_genes_de_interes          <- c()
    cantidad_bines_de_genes_de_interes <- c()
    correlaciones <- signif(seq(.5, 1, 0.01), 2)
    for(correlacion_minima in correlaciones){
      adyacencia                                                               <- correlaciones_bines[rv$bines_con_uso_diferencial, rv$bines_con_uso_diferencial]
      adyacencia[adyacencia < correlacion_minima]                              <- 0
      diag(adyacencia)                                                         <- 0

      #Generamos el grafo a partir de la matriz de adjacencia
      g_bines                            <- componente_gigante(graph_from_adjacency_matrix(adyacencia, "undirected", weighted=TRUE, diag=FALSE))
      genes_de_bines                     <- strsplit2(names(V(g_bines)), ":")
      cantidad_genes_de_interes          <- c(cantidad_genes_de_interes, length(intersect(genes_de_bines[, 1], genes_de_interes_marcelo)))
      cantidad_bines_de_genes_de_interes <- c(cantidad_bines_de_genes_de_interes, sum(genes_de_bines[, 1] %in% genes_de_interes_marcelo))
      enlaces                            <- c(enlaces, ecount(g_bines))
      vertices                           <- c(vertices, vcount(g_bines))
      if(input$correlacion_minima == correlacion_minima){
        rv$g_bines <- g_bines
      }      
    }
    layout(matrix(1:4, nrow=4, ncol=1))
    plot(correlaciones, vertices, xlab="Correlación", ylab="Vértices en la red", cex.lab=1.5)
    grid(col = "gray", lty = "dotted")    
    abline(v=input$correlacion_minima, col="red")
    y <- vertices[correlaciones == input$correlacion_minima]    
    text(x=1, y=round(max(vertices)*0.85), labels=paste0("(", input$correlacion_minima, ", ", y, ")"), cex=1.5, pos=2)
    
    plot(correlaciones, enlaces, xlab="Correlación", ylab="Enlaces en la red", cex.lab=1.5)
    grid(col = "gray", lty = "dotted")    
    abline(v=input$correlacion_minima, col="red")
    y <- enlaces[correlaciones == input$correlacion_minima]    
    text(x=1, y=round(max(enlaces)*0.85), labels=paste0("(", input$correlacion_minima, ", ", y, ")"), cex=1.5, pos=2)    
    
    plot(correlaciones, cantidad_genes_de_interes, xlab="Correlación", ylab="Genes de interés en la red", cex.lab=1.5)
    grid(col = "gray", lty = "dotted")    
    abline(v=input$correlacion_minima, col="red")
    y <- cantidad_genes_de_interes[correlaciones == input$correlacion_minima]    
    text(x=1, y=round(max(cantidad_genes_de_interes)*0.85), labels=paste0("(", input$correlacion_minima, ", ", y, ")"), cex=1.5, pos=2)
    
    plot(correlaciones, cantidad_bines_de_genes_de_interes, xlab="Correlación", ylab="Bines de interés en la red", cex.lab=1.5)
    grid(col = "gray", lty = "dotted")
    abline(v=input$correlacion_minima, col="red")
    y <- cantidad_bines_de_genes_de_interes[correlaciones == input$correlacion_minima]    
    text(x=1, y=round(max(cantidad_bines_de_genes_de_interes)*0.85), labels=paste0("(", input$correlacion_minima, ", ", y, ")"), cex=1.5, pos=2)
    
    print(paste(vcount(rv$g_bines), ecount(rv$g_bines)))
  })

  #Guarda los perfiles
  observeEvent(input$siguiente, {
    readme <- data.frame(dato=c("lfchange_limite", "qvalue_limite", "cantidad_de_tiempos_limite", "correlacion_minima"),
                         valor=c(log2(input$lfchange_limite/100+1), input$qvalue_limite, input$cantidad_de_tiempos_limite, input$correlacion_minima))
    g_bines <- rv$g_bines
    save(g_bines, readme, file="5_spliceoma.Rdata")
    stopApp()
  })
}
shinyApp(ui = ui, server = server)

