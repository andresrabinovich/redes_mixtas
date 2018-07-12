##############################################################################################
#ANALISIS DE RED MIXTA
#Autor: Andrés Rabinovich
#Creación: 08/06/2018
#Última modificación: XX/XX/XXXX (por XXX)
##############################################################################################

#Librerías que necesita
library(limma)
library(igraph)
library(shiny)
library(Biostrings)
library(ASpli)
library(AnnotationDbi)
library(readr)
setwd("/home/arabinov/doctorado/programacion/redes_mixtas/")

#Funciones varias para grafos
source("pipeline/funciones_grafos.R")

#Levantamos las redes
(load("pipeline_archivos/3_transcriptoma.Rdata"))
(load("pipeline_archivos/5_spliceoma.Rdata"))
(load("pipeline_archivos/6_red_mixta.Rdata"))

fasta    <- readDNAStringSet("/xdata1/tcourseT/ATtcourseT/comp/rabi/genoma/genome.fa")
txdb     <- loadDb("/xdata1/tcourseT/ATtcourseT/comp/rabi/genoma/genoma.sqlite")
features <- binGenome(txdb)

#Levantamos los simbolos de los genes
at_simbolos        <- at_to_symbol    <- read.table("pipeline_archivos/at_to_symbol_map", sep = "\t", header = F, stringsAsFactors = F)
at_simbolos        <- at_to_symbol$V2
names(at_simbolos) <- at_to_symbol$V1

#Encuentra los enlaces intrared (bin contra gen)
enlaces                <- ends(g, E(g))
enlaces_interred       <- setdiff(grep(":", enlaces[, 2]), grep(":", enlaces[, 1]))
enlaces_intrared_bines <- intersect(grep(":", enlaces[, 1]), grep(":", enlaces[, 2]))
enlaces_intrared_genes <- setdiff(1:ecount(g), c(enlaces_interred, enlaces_intrared_bines))
enlaces_interred       <- enlaces[enlaces_interred, ]
enlaces_intrared_bines <- enlaces[enlaces_intrared_bines, ]
enlaces_intrared_genes <- enlaces[enlaces_intrared_genes, ]
datos_bines            <- data.frame(features@bins[unique(enlaces_interred[, 2]), ])
rownames(datos_bines)  <- unique(enlaces_interred[, 2])

#Función para generar todas las subsecuencias de tamaño dado a partir de una secuencia más grande
subsecuencias <- function(min_len, max_len, secuencia){
  secuencias <- c()
  for(longitud in min_len:max_len){
    for(comienzo in 1:(nchar(secuencia)-longitud+1)){
      secuencias <- c(secuencias, substr(secuencia, comienzo, comienzo+longitud-1))
    }
  }
  return(unique(secuencias))
}

#Modulos
l_enlaces_interred <- setNames(lapply(unique(enlaces_interred[, 1]), function(x){
  return(enlaces_interred[which(enlaces_interred[, 1] == x), 2])
}), unique(enlaces_interred[, 1]))

#Secuencias de los modulados
l_secuencias <- setNames(lapply(unique(enlaces_interred[, 2]), function(n){
  a          <- datos_bines[n, c("seqnames", "start", "end", "strand")]
  s          <- subseq(fasta[a$seqnames[1]], start=a$start[1], end=a$end[1])
  if(a$strand == "-") s <- reverseComplement(s)
  return(as.character(s))
}), unique(enlaces_interred[, 2]))

#Subsecuencias entre 6 y 12 caracteres de los módulos
l_subsecuencias2 <- lapply(l_secuencias, function(s){
  secuencias <- subsecuencias(min_len = 6, max_len = 12, secuencia  = s)
  return(secuencias)
})

#Secuencias con máxima frecuencia de aparición en cada módulo
l_secuencias_con_maxima_frecuencia <- lapply(l_enlaces_interred, function(modulo){
  if(length(modulo) < 3) return(NULL)
  secuencias <- table(unlist(l_subsecuencias[modulo]))
  #Devuelve la/s secuencia/s con máxima frecuencia de aparición, en cuántas aparece/n y la cantidad de elementos total
  return(c(c(max(secuencias), length(modulo)), 
           names(which(secuencias == max(secuencias))))
         )
})
l_secuencias_con_maxima_frecuencia <- Filter(length, l_secuencias_con_maxima_frecuencia)


#Calculamos pvalues para subsecuencias con maxima frecuencia a partir de distribución nula. Para generar la distribución nula
#replicamos 10000 veces grupos de bines de un tamaño igual al tamaño del módulo.
replicaciones <- 10000
n_modulos     <- length(l_secuencias_con_maxima_frecuencia)
i             <- 0
#universo      <- unique(enlaces_interred[, 2]) #Unique o no unique? Qué representa mejor la distribución nula?
universo      <- enlaces_interred[, 2] #Unique o no unique? Qué representa mejor la distribución nula?
l_pval <- lapply(l_secuencias_con_maxima_frecuencia, function(x){
  i <<- i + 1
  print(paste0(i, "/", n_modulos))
  tamano_modulo <- as.numeric(x[2])
  hits          <- as.numeric(x[1])
  pval          <- c()
  for(patron in x[-c(1, 2)]){
    #Calcula un pvalue a las secuencias con máxima frecuencia de aparición haciendo bootstrap
    #lll <- sample(unique(enlaces_interred[, 2]), tamano_modulo, replace = F)
    b<-replicate(replicaciones, 
               length(grep(patron, l_secuencias[sample(universo, tamano_modulo, replace = F)]))
               )
    pval = c(pval, sum(b >= hits)/replicaciones)
  }
  names(pval) <- x[-c(1, 2)]
  return(pval)
})

l_qval <- lapply(l_pval, p.adjust, method="fdr")

#Usamos tomtom para ver si hay motif conocidos para las secuencias que encontramos sobrerepresentadas.
#Bases de datos de motifs
dbs                               <- c("ArabidopsisDAPv1.meme", "ArabidopsisPBM_20140210.meme")

#El template para el archivo de secuencias para pasarle a tomtom
template_meme                     <- read_file("~/meme/secuencias/template_meme.txt")

#Template para generar la matriz de probabilidades de secuencias
template_matriz_de_probabilidades <- data.frame(A=c(1, 0, 0, 0),
                                                C=c(0, 1, 0, 0),
                                                G=c(0, 0, 1, 0),
                                                T=c(0, 0, 0, 1))

#Buscamos cada una de las secuencias en la base de datos con tomtom
l_tomtom <- lapply(l_qval, function(secuencias){
  #Matriz de resultados vacia
  tomtom           <- matrix(ncol=7, nrow=0)
  colnames(tomtom) <- c("secuencia", "target", "optimal_offset", "pvalue", "qvalue", "target_consensus", "db")
  
  #Iteramos cada secuencia de la ista
  for(secuencia in names(secuencias)){
    
    #Generamos la matriz de probabilidades
    matriz_de_probabilidades <- t(template_matriz_de_probabilidades[, strsplit2(secuencia, "")])
    
    #La transformamos en un string con el formato del archivo de secuencias de tomtom
    matriz_de_probabilidades <- paste0(apply(matriz_de_probabilidades, 1, function(s){
      return(paste0(" ", paste(s, collapse="  "), "\n"))
    }), collapse="")
    
    #Reemplazamos en el template la secuencia
    meme <- gsub("@secuencia", secuencia, template_meme)
    
    #Reemplazamos en el template la matriz de probabilidades
    meme <- gsub("@matriz_de_probabilidades", matriz_de_probabilidades, meme)
    
    #Generamos el archivo de secuencias para tomtom
    write_file(meme, "~/meme/secuencias/meme.txt")
    
    #Iteramos en las bases de datos
    for(db in dbs){
      #Corremos tomtom para esta secuencia y la base de datos seleccionada
      system(paste0("~/meme/bin/tomtom -no-ssc -oc ~/meme/salidas -verbosity 1 -min-overlap 5 -mi 1 -dist pearson -evalue -thresh 10.0 ~/meme/secuencias/meme.txt ~/meme/motif_databases/ARABD/", db))
      
      #Cargamos la salida de tomtom
      salida <- read.table("~/meme/salidas/tomtom.tsv", header = T)
      
      #Guardamos los resultados de tomtom en la matriz de resultados
      if(nrow(salida) > 0){
        tomtom <- rbind(tomtom, data.frame(secuencia=secuencia, target=salida$Target_ID, optimal_offset=salida$Optimal_offset,
                                         pvalue=salida$p.value, qvalue=salida$q.value, target_consensus=salida$Target_consensus, 
                                         db=db))
      }
    }
  }
  #Devolvemos la matriz de resultados
  return(tomtom)
})
#Limpiamos los archivos de tomtom
limpiar <- file.remove("~/meme/salidas/tomtom.tsv", "~/meme/salidas/tomtom.xml", "~/meme/salidas/tomtom.html", "~/meme/secuencias/meme.txt")



