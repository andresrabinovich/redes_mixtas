#############################################################################################
#SCRIPT QUE GENERA UNA LISTA CON TODAS LAS CUENTAS POR BIN PARA CADA UNO DE LOS TRATAMIENTOS.                   
#Autor: Andrés Rabinovich
#Creación: 01/03/2017
#Última modificación: 22/01/2018 (por Andrés Rabinovich)
#############################################################################################

#Instalador de versión local de aspli
#install.packages("~/doctorado/programacion/rna_seq/paquetes/ASpli_1.3.12.2002.tar.gz", repos = NULL, type="source")

#Incluimos ASpli para leer las cuentas
library(AnnotationDbi)
library(ASpli)

#Guardamos el directorio de trabajo actual y levantamos el directorio de los archivos bam
wd <- getwd()
setwd("/xdata1/tcourseT/ATtcourseT/comp/04_DE/")

#Levantamos las anotaciones
anotaciones <- "genome.sqlite"
txdb        <- loadDb(anotaciones)
print("Levantando anotaciones")
features    <- binGenome(txdb)

#Inicializamos los data.frame que van a contener las cuentas
cuentas_genes    <- cuentas_bines <- NULL 


#Cargamos los bam, 12 puntos temporales, dos réplicas, cuatro temperaturas y generamos un archivo de cuentas por temp.
puntos_temporales <- 1:12
replicas          <- c("A", "B") 
temperaturas      <- c("22", "12", "17", "27")
cuentas_junturas  <- setNames(vector("list", length(temperaturas)), temperaturas)

print("Cargando los bam por temperatura")
for(temperatura in temperaturas){
  print(paste("Procesando archivo de", temperatura, "grados"))
  
  #Leemos todos los archivos bam de una misma temperatura (los 12 puntos temporales x las 2 réplicas)
  bamFiles   <- paste(rep("at_", time=length(puntos_temporales)*length(replicas)), 
                      rep(temperatura, each=length(puntos_temporales)*length(replicas)),
                      rep("_", time=length(puntos_temporales)*length(replicas)),      
                      rep(puntos_temporales, each=length(replicas)), 
                      rep("_", time=length(puntos_temporales)*length(replicas)),
                      rep(replicas, time=length(puntos_temporales)),
                      rep(".bam", time=length(puntos_temporales)*length(replicas)),
                      sep="")
  bkg        <- rep(temperatura, time=length(puntos_temporales)*length(replicas))
  time       <- rep(paste(rep("T", time=length(puntos_temporales)), puntos_temporales, sep=""), each=length(replicas))
  targets    <- data.frame(bam=bamFiles,
                        condition=paste(bkg, time, sep="."),
                        row.names=sub("\\.bam$", "", bamFiles),
                        stringsAsFactors = FALSE)
  bam        <- loadBAM(targets)
  
  #Lee las cuentas que caen sobre cada feature 
  counts     <- readCounts(features, bam, targets, readLength = 100, maxISize = 5000)

  #Armamos los data.frame de cuentas para los genes y para los bines  
  if(is.null(cuentas_genes) | is.null(cuentas_bines)){
    cuentas_genes    <- countsg(counts)
    cuentas_bines    <- countsb(counts)

  }else{
    cuentas_genes    <- cbind(cuentas_genes, ASpli:::.extractCountColumns(countsg(counts), targets))
    cuentas_bines    <- cbind(cuentas_bines, ASpli:::.extractCountColumns(countsb(counts), targets))
  }  
  cuentas_junturas[[temperatura]] <- countsj(counts)

}

#Guardamos la lista de cuentas en un Rdata
save(cuentas_genes, cuentas_bines, cuentas_junturas, file = paste0(wd, "/cuentas/cuentas.Rdata"))

#Limpiamos todo
rm(list=c("bam", "counts"))
gc()
setwd(wd)
print(paste0("Todos los archivos procesados correctamente. Archivo ", wd, "/cuentas/cuentas.Rdata creado."))
