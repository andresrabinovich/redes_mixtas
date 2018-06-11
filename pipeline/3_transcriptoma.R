#Librerías que necesita
library(shiny)
library(igraph)
library(limma) #Para usar strsplit2
library(WGCNA) #Para usar bicor

#Directorio de trabajo
setwd("/home/arabinov/doctorado/programacion/redes_cosplicing/pipeline_archivos/")

#Funciones varias para grafos
source("../pipeline/funciones_grafos.R")

#Elegimos las condiciones para la que vamos a generar la red
(load("1_seleccion_de_condiciones.Rdata"))

#Levantamos las cuentas de los genes 
(load("2_genes_prefiltrados.Rdata"))

#Levantamos los perfiles unicamente de la temperatura (sin la referencia)
iTemp    <- grep(paste0("at_",temperatura,"_"), colnames(cuentas_genes))

#Filtra y genera un perfil del log de cuentas por millon
perfiles <- apply(cuentas_genes[, iTemp],2,function(x){log(x/sum(x)*1e6+0.01)})

#Promedia replicas
perfiles <- 0.5*(perfiles[,(1:ncol(perfiles))%%2==1]+perfiles[,(1:ncol(perfiles))%%2==0])

#Calcula la bicorrelacion entre los perfiles. Si el MAD es cero, usa la correlación de pearson. 
#Hay pocos genes en esa situación
bcor <- bicor(t(perfiles))

#interfaz de shiny
ui <- fluidPage(
  titlePanel("Transcriptoma"),
  
  sidebarLayout(
    sidebarPanel(
      helpText("Filtar genes y armar red transcriptoma en base a: "),
      
      #valores límite para filtrar genes
      selectInput("cantidad_de_tiempos_limite", 
                  label = "Cantidad de puntos temporales con expresión diferencial",
                  choices = 1:12,
                  selected = 1),
      
      sliderInput("qvalue_limite", 
                  label = "Máximo qvalue:",
                  min = 0, max = 1, value = 0.05, step=0.01),
      
      sliderInput("lfchange_limite", 
                  label = "Mínimo cambio a detectar (%):",
                  min = 0, max = 200, value=50, step=25),

      sliderInput("correlacion_minima", 
                  label = "Correlación mínima:",
                  min = 0.7, max = 1, value=0.85, step=0.01),
      
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
      textOutput("cantidad_de_genes"),
      textOutput("cantidad_de_reguladores"),
      htmlOutput("genes_de_interes"),
      plotOutput("componente_gigante", height = "600px")
    )
  )
)

#Server de shiny para filtrar y actualizar los datos
server <- function(input, output) {
  
  #Inicializa las variables reactive para que sean accesibles desde otros elementos
  rv <- reactiveValues(genes_con_expresion_diferencial=NULL) 
  
  #Filtra y actualiza los genes de interes que pasaron el filtro    
  filtrar_genes <- reactive({
    #Pide que el gen haya cambiado al menos en cantidad_de_tiempos_limite tiempos en un lfchange_limite
    genes_con_alto_lfchange <- rownames(lfchange_genes)[apply(lfchange_genes, 1, function(x){
      sum(abs(x) > log2(as.numeric(input$lfchange_limite)/100+1)) >= as.numeric(input$cantidad_de_tiempos_limite)
    })]
    
    #Pide que el gen tenga un qvalue < qvalue_limite en cantidad_de_tiempos_limite tiempos para validad el fold change
    genes_con_qvalue_bajo  <- rownames(qvalues_genes)[apply(qvalues_genes, 1, function(x){
      sum(x<as.numeric(input$qvalue_limite)) >= as.numeric(input$cantidad_de_tiempos_limite)
    })]
    
    #Me quedo con los genes que pasen ambos filtros
    rv$genes_con_expresion_diferencial <- sort(intersect(genes_con_alto_lfchange, genes_con_qvalue_bajo))
    
  })

  #Filtra los genes
  output$cantidad_de_genes <- renderText({ 
    filtrar_genes()
    paste("Total de genes que pasaron el filtro:", length(rv$genes_con_expresion_diferencial))
  })

  #Reguladores que pasaron el filtro  
  output$cantidad_de_reguladores <- renderText({ 
    paste("Total de reguladores que pasaron el filtro:", length(intersect(rv$genes_con_expresion_diferencial, reguladores_de_expresion)))
  })

  observeEvent(input$genes_de_interes_marcelo, {
    rv$genes_de_interes_marcelo <- input$genes_de_interes_marcelo 
  })    

  #Filtra los genes
  output$genes_de_interes <- renderText({ 
    genes_de_interes_marcelo                       <- trimws(strsplit2(input$genes_de_interes_marcelo, ","))
    genes_de_interes_marcelo_que_pasaron_el_filtro <- intersect(rv$genes_con_expresion_diferencial, genes_de_interes_marcelo)
    HTML(paste0("Genes de interés que pasaron el filtro (", length(genes_de_interes_marcelo_que_pasaron_el_filtro), "):</br>", paste0(genes_de_interes_marcelo_que_pasaron_el_filtro, collapse="</br>")))
  })   
  
  #Plotea el gráfico de cantidad de elementos en la componente gigante en función de la correlación
  output$componente_gigante <- renderPlot({ 
    genes_de_interes_marcelo  <- trimws(strsplit2(input$genes_de_interes_marcelo, ","))
    enlaces                   <- c()
    vertices                  <- c()
    cantidad_genes_de_interes <- c()
    correlaciones <- signif(seq(.7, 1, 0.01), 2)
    for(correlacion_minima in correlaciones){
      adyacencia                                                               <- bcor[rv$genes_con_expresion_diferencial, rv$genes_con_expresion_diferencial]
      adyacencia[adyacencia < correlacion_minima]                              <- 0
      diag(adyacencia)                                                         <- 0
      
      adyacencia_columna <- adyacencia_fila <- adyacencia
      adyacencia_fila[!rownames(adyacencia_fila) %in% reguladores_de_expresion, ] <- 0
      adyacencia_columna[, !colnames(adyacencia_columna) %in% reguladores_de_expresion] <- 0
      adyacencia <- .5*(adyacencia_columna+adyacencia_fila)
      adyacencia[lower.tri(adyacencia)] <- 0

      #Generamos el grafo a partir de la matriz de adjacencia
      g_genes                   <- componente_gigante(graph_from_adjacency_matrix(adyacencia, "directed", weighted=TRUE, diag=FALSE))
      cantidad_genes_de_interes <- c(cantidad_genes_de_interes, sum(names(V(g_genes)) %in% genes_de_interes_marcelo))
      enlaces                   <- c(enlaces, ecount(g_genes))
      vertices                  <- c(vertices, vcount(g_genes))
      if(input$correlacion_minima == correlacion_minima){
        rv$g_genes <- g_genes
      }
    }
    layout(matrix(1:3, nrow=3, ncol=1))

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
    
  })  

    #Guarda los perfiles
  observeEvent(input$siguiente, {
    perfiles_genes <- perfiles[rv$genes_con_expresion_diferencial, ]
    readme <- data.frame(dato=c("lfchange_limite", "qvalue_limite", "cantidad_de_tiempos_limite", "correlacion_minima"),
                         valor=c(log2(input$lfchange_limite/100+1), input$qvalue_limite, input$cantidad_de_tiempos_limite, input$correlacion_minima))
    g_genes <- rv$g_genes
    save(perfiles_genes, g_genes, readme, file="3_transcriptoma.Rdata")
    stopApp()
  })  
}
shinyApp(ui = ui, server = server)

#(load("../transcriptoma.RData"))
