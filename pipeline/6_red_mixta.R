##############################################################################################
#SCRIPT QUE GENERA UNA RED MIXTA
#Autor: Andrés Rabinovich
#Creación: 16/03/2018
#Última modificación: 05/06/2018 (por Andrés Rabinovich)
##############################################################################################

#Librerías que necesita
library(igraph)
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

#Levantamos los simbolos de los genes
at_simbolos        <- at_to_symbol <- read.table("pipeline_archivos/at_to_symbol_map", sep = "\t", header = F, stringsAsFactors = F)
at_simbolos        <- at_to_symbol$V2
names(at_simbolos) <- at_to_symbol$V1

#Levantamos los genes relacionados con splicing en el transcriptoma
proteinas_relacionadas_con_splicing <- unique(c(poi$SP, reguladores$gene_id[reguladores$tipo_de_regulador %in% c("RBP", "NO CLASIFICADOS")], reguladores_de_splicing_extra))
#proteinas_relacionadas_con_splicing <- c(proteinas_relacionadas_con_splicing, "AT4G31120")

#Unimos las dos redes en un solo grafo pero por ahora no están vinculadas entre si
g <- as.undirected(g_genes) + g_bines
E(g)$weight  <- c(E(g)$weight_2[!is.na(E(g)$weight_2)], E(g)$weight_1[!is.na(E(g)$weight_1)])
#g           <- graph_from_edgelist(rbind(ends(g_genes, E(g_genes)), ends(g_bines, E(g_bines))), directed = F)
#E(g)$weight <- c(E(g_genes)$weight, E(g_bines)$weight)
#g           <- simplify(g)

#Calculamos la correlación entre bines y proteinas relacionadas con splicing
perfiles_bines  <- perfiles_bines[names(V(g_bines)), ]
perfiles_genes  <- perfiles_genes[intersect(names(V(g_genes)), proteinas_relacionadas_con_splicing), ]
#ccor           <- abs(cor(t(perfiles_bines), t(perfiles_genes)))
ccor            <- abs(cor(t(perfiles_bines), t(perfiles_genes)))
noabsccor       <- setNames(melt(cor(t(rbind(perfiles_bines, perfiles_genes)))), c('targets', 'sf', 'cor'))
adyacencia      <- setNames(melt(ccor), c('targets', 'sf', 'pesos'))

#interfaz de shiny
ui <- fluidPage(
  titlePanel("Red Mixta"),
  
  sidebarLayout(
    sidebarPanel(

      sliderInput("correlacion_minima", 
                  label = "Correlación mínima:",
                  min = 0.65, max = 1, value=0.81, step=0.01),
      
      textAreaInput("genes_de_interes_marcelo", 
                    value=paste0(c("AT3G09600",
                                   "AT5G64170", 
                                   "AT3G54500", 
                                   "AT1G09140", 
                                   "AT5G02810", 
                                   "AT5G02840", 
                                   "AT2G46830", 
                                   "AT1G01060", 
                                   "AT2G21660",
                                   "AT2G18740"), collapse=", "),
                    label = "Genes de interes",
                    width = '100%',
                    height = '150px'),
      
      br(),
      actionButton("siguiente", label = "Siguiente")
    ),
    
    mainPanel(
      htmlOutput("genes_de_interes"),
      plotOutput("cantidad_de_enlaces_gen_bin", height = "600px")
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
      enlaces <- c(enlaces, sum(adyacencia$pesos > correlacion_minima))
    }
    plot(correlaciones, enlaces, xlab="Correlación", ylab="Enlaces inter redes", cex.lab=1.5)
    grid(col = "gray", lty = "dotted")
    cantidad_de_enlaces_elegidos <- enlaces[correlaciones == input$correlacion_minima]
    text(x=input$correlacion_minima, y=cantidad_de_enlaces_elegidos, labels=paste0("(", input$correlacion_minima, ", ", cantidad_de_enlaces_elegidos, ")"), cex=0.9, pos=2)
    abline(v=input$correlacion_minima, col="red")
  })
  
  observeEvent(input$genes_de_interes_marcelo, {
    rv$genes_de_interes_marcelo <- input$genes_de_interes_marcelo 
  })    
  
  #Filtra los genes de interes
  output$genes_de_interes <- renderText({ 
    genes_de_interes_marcelo <- trimws(strsplit2(input$genes_de_interes_marcelo, ","))
    salida <- ""
    for(i in genes_de_interes_marcelo){
      if(i %in% proteinas_relacionadas_con_splicing){
        salida <- paste0(salida, i, " (gen): ", sum(ccor[, colnames(ccor) %in% i] > input$correlacion_minima), " - ")
      }
      salida <- paste0(salida, i, " (bin): ", sum(ccor[grep(i, rownames(ccor)), ] > input$correlacion_minima), "</br>")
    }
    HTML(salida)
  })   
  
  #Agrega los nuevos enlaces a la red
  observeEvent(input$siguiente, {
    enlaces_totales     <- which(adyacencia$pesos > input$correlacion_minima)
    g                   <<- add_edges(g, c(t(adyacencia[enlaces_totales, c(2,1)])), attr=list(weight = adyacencia$pesos[enlaces_totales]))
    readme <- data.frame(dato="correlacion_minima",
                         valor=input$correlacion_minima)
    save(g, readme, file="pipeline_archivos/6_red_mixta.Rdata")
    stopApp()
  })
}
shinyApp(ui = ui, server = server)
